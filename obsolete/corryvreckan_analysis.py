#!/usr/bin/python
import os
import math
import subprocess
import array
import numpy as np
import ROOT
from ROOT import *
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
from scipy.optimize import curve_fit
from skspatial.objects import Line, Sphere
from skspatial.plotting import plot_3d

import config
from config import *
import utils
from utils import *
import svd_fit
from svd_fit import *
import chi2_fit
from chi2_fit import *
import hists
from hists import *

import argparse
parser = argparse.ArgumentParser(description='corryvreckan_analysis.py...')
parser.add_argument('-f', metavar='input file', required=True,  help='full path to input file')
parser.add_argument('-r', metavar='run type', required=True,  help='source/cosmics/beam')
parser.add_argument('-s', metavar='is simulation', required=True,  help='is simulation? [1/0]')
parser.add_argument('-v', metavar='do vertex?', required=False,  help='do vertex? [1/0]')
parser.add_argument('-p', metavar='plot?', required=False,  help='plot? [1/0]')
parser.add_argument('-n', metavar='n events to process', required=False,  help='n events to process')
parser.add_argument('-c', metavar='is CVMFS?', required=False,  help='is CVMFS? [1/0]')
argus = parser.parse_args()
inrootfilename = argus.f
runtype = argus.r
nmaxtoprocess  = int(argus.n) if(argus.n is not None) else -1
doVtx  = True if(argus.v=="1") else False
doplot = True if(argus.p=="1") else False
### important - this will reset the config
mc    = True if(argus.s=="1") else False
cvmfs = True if(argus.c is not None and argus.c=="1") else False
reconfig(mc,cvmfs)

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetOptStat(0)

print("-------------------")
print("Monte Carlo run") if(isMC) else print("Data run")
print("Run type:",runtype)
print("-------------------")

### randomness?
rnd = TRandom()
rnd.SetSeed()
doRnd = False

### globals
absRes  = 0.15
absChi2 = 20
if(runtype=="source"):
    absRes  *= 100
    absChi2 *= 100


###########################################################
###########################################################
###########################################################

### get noisy pixels:
def getPixels2Mask(h1D,h2D,threshold):
    masked = {}
    for bx in range(h2D.GetNbinsX()+1):
        for by in range(h2D.GetNbinsY()+1):
            x = int(h2D.GetXaxis().GetBinCenter(bx))
            y = int(h2D.GetYaxis().GetBinCenter(by))
            hpix = h2D.GetBinContent(bx,by)
            ipix = h2D.FindBin(x,y)
            if(hpix>threshold):
                masked.update({ipix:[x,y]})
    return masked


### check each pixel if it is noise
def isNoise(pix_col,pix_row,h2D,masked):
    isnoise = False
    for i in range(pix_col.size()):
        ipix = h2D.FindBin(pix_col[i],pix_row[i])
        if(ipix in masked):
            isnoise = True
            break
    return isnoise


def fillOcc(pix_col,pix_row,h1D,h2D,masked,doMask):
    for p in range(pix_col.size()):
        ixpix = pix_col[p]
        iypix = pix_row[p]
        ipix  = h2D.FindBin(ixpix,iypix)
        if(doMask):
            if(ipix in masked):
                continue
        h1D.AddBinContent(ipix,1)
        h2D.Fill(ixpix,iypix)


def fillFitOcc(params,hname2,hname3):
    for det in detectors:
        x,y,z = line(rdetectors[det][2],params)
        histos[hname2+"_"+det].Fill(x,y)
        histos[hname3].Fill(x,y,z)


def fillClsHists(det,clsnpix,clsncol,clsnrow,postfit=False):
    pf = "postfit_" if(postfit) else ""
    histos["h_ncls_"+pf+det].Fill(clsncol.size())
    for c in range(clsncol.size()):
        histos["h_cls_size_"+pf+det].Fill(clsnpix[c])
        histos["h_cls_size_ncol_"+pf+det].Fill(clsncol[c])
        histos["h_cls_size_nrow_"+pf+det].Fill(clsnrow[c])


def SVD_candidate(clsx,clsy,clsz,clsncol,clsnrow,vtx=[],evtx=[]):
    isvtx = (len(vtx)>0 and len(evtx)>0)
    clusters = [ [xVtx, yVtx, zVtx] ]  if(isvtx) else []
    clerrors = [ [exVtx,eyVtx,ezVtx] ] if(isvtx) else []
    for det in detectors:
        clusters.append( [clsx[det]+offsets_x[det], clsy[det]+offsets_y[det], clsz[det]] )
        clerrors.append( [clsncol[det]*pix_x/2.,    clsnrow[det]*pix_y/2.,    ezCls] )
    points = np.array(clusters)
    errors = np.array(clerrors)
    return points,errors


def Chi2_candidate(clsx,clsy,clsz,clsncol,clsnrow,vtx=[],evtx=[]):
    isvtx = (len(vtx)>0 and len(evtx)>0)
    clusters_x = [xVtx]  if(isvtx) else []
    clusters_y = [yVtx]  if(isvtx) else []
    clusters_z = [zVtx]  if(isvtx) else []
    clerrors_x = [exVtx] if(isvtx) else []
    clerrors_y = [eyVtx] if(isvtx) else []
    clerrors_z = [ezVtx] if(isvtx) else []
    for det in detectors:
        clusters_x.append( clsx[det]+offsets_x[det] )
        clusters_y.append( clsy[det]+offsets_y[det] )
        clusters_z.append( clsz[det] )
        clerrors_x.append( clsncol[det]*pix_x/2. )
        clerrors_y.append( clsnrow[det]*pix_y/2. )
        clerrors_z.append( ezCls )
    points = np.array([ clusters_x,clusters_y,clusters_z ])
    errors = np.array([ clerrors_x,clerrors_y,clerrors_z ])
    return points,errors


def fill_trk2cls_residuals(points,direction,centroid,hname):
    for det in detectors:
        dx,dy = res_track2cluster(det,points,direction,centroid)
        histos[hname+"_x_"+det].Fill(dx)
        histos[hname+"_y_"+det].Fill(dy)


def fill_trk2tru_residuals(mcx,mcy,points,direction,centroid,hname):
    mcparticles = [ [mcx[det]-chipX/2.+pix_x/2., mcy[det]-chipY/2.+pix_y/2.],
                    [mcx[det]-chipX/2.+pix_x/2., mcy[det]-chipY/2.+pix_y/2.],
                    [mcx[det]-chipX/2.+pix_x/2., mcy[det]-chipY/2.+pix_y/2.] ]
    for det in detectors:
        dx,dy = res_track2truth(det,mcparticles,points,direction,centroid)
        histos[hname+"_x_"+det].Fill(dx)
        histos[hname+"_y_"+det].Fill(dy)


def fill_trk2vtx_residuals(vtx,direction,centroid,hname):
    dxv,dyv = res_track2vertex(vtx,direction,centroid)
    histos[hname+"_x"].Fill(dxv)
    histos[hname+"_y"].Fill(dyv)



### get data
def Run(inrootfilename,outrootfilename):
    tfi = TFile(inrootfilename,"READ")
    tt  = tfi.Get("tt")
    
    noise_threshold = {}
    masked          = {}
    h1D_noise       = {}
    h2D_noise       = {}
    for det in detectors:
        noise_threshold.update( {det:-1} )
        masked.update( {det:{}} )
        h1D_noise.update( {det:tfi.Get("h_pixocc1D_"+det)} )
        h2D_noise.update( {det:tfi.Get("h_pixocc2D_"+det)} )
    
    for det in detectors:    
        # noise_threshold.update( {det:-1} )
        avg,std,threshold = getNoiseThreshold(h1D_noise[det],pTrim,nSigma,zeroSupp)
        print(det,": avg,std:",avg,std,"--> threshold:",threshold,"(pTrim=",pTrim,")")
        noise_threshold[det] = threshold if(threshold>noise_threshold[det]) else noise_threshold[det]
        print("Final noise threshold for",det,"is:",noise_threshold[det])
        masked.update( {det:getPixels2Mask(h1D_noise[det],h2D_noise[det],noise_threshold[det])} )
        print("Masked pixels for threshold of",noise_threshold[det],"in",det,"is:",masked[det])
    
    sphere_points_a_x = []
    sphere_points_a_y = []
    sphere_points_a_z = []
    sphere_points_b_x = []
    sphere_points_b_y = []
    sphere_points_b_z = []
    
    
    tfo = TFile(outrootfilename,"RECREATE")
    tfo.cd()
    book_histos(absRes,absChi2)
    print("Entries=",tt.GetEntries())
    nevents = 0
    for event in tt:
        
        orig_evt = event.evt[0]
        
        ### some general hostos
        histos["h_events"].Fill(0.5)
        histos["h_cutflow"].Fill( cuts.index("All") )
        histos["h_npix"].Fill(event.ALPIDE_2_pix_col.size())
        histos["h_npix"].Fill(event.ALPIDE_1_pix_col.size())
        histos["h_npix"].Fill(event.ALPIDE_0_pix_col.size())
        
        ### fill occupancy plots for all events, not only the ones participating in the tracking WITHOUT masking
        fillOcc(event.ALPIDE_2_pix_col,event.ALPIDE_2_pix_row,histos["h_pix_occ_1D_ALPIDE_2"],histos["h_pix_occ_2D_ALPIDE_2"],masked["ALPIDE_2"],False)
        fillOcc(event.ALPIDE_1_pix_col,event.ALPIDE_1_pix_row,histos["h_pix_occ_1D_ALPIDE_1"],histos["h_pix_occ_2D_ALPIDE_1"],masked["ALPIDE_1"],False)
        fillOcc(event.ALPIDE_0_pix_col,event.ALPIDE_0_pix_row,histos["h_pix_occ_1D_ALPIDE_0"],histos["h_pix_occ_2D_ALPIDE_0"],masked["ALPIDE_0"],False)
        
        ### fill occupancy plots for all events, not only the ones participating in the tracking WITH masking
        fillOcc(event.ALPIDE_2_pix_col,event.ALPIDE_2_pix_row,histos["h_pix_occ_1D_masked_ALPIDE_2"],histos["h_pix_occ_2D_masked_ALPIDE_2"],masked["ALPIDE_2"],True)
        fillOcc(event.ALPIDE_1_pix_col,event.ALPIDE_1_pix_row,histos["h_pix_occ_1D_masked_ALPIDE_1"],histos["h_pix_occ_2D_masked_ALPIDE_1"],masked["ALPIDE_1"],True)
        fillOcc(event.ALPIDE_0_pix_col,event.ALPIDE_0_pix_row,histos["h_pix_occ_1D_masked_ALPIDE_0"],histos["h_pix_occ_2D_masked_ALPIDE_0"],masked["ALPIDE_0"],True)
        
        ### check if there's no noise
        isnoise2 = isNoise(event.ALPIDE_2_pix_col,event.ALPIDE_2_pix_row,histos["h_pix_occ_2D_ALPIDE_2"],masked["ALPIDE_2"])
        isnoise1 = isNoise(event.ALPIDE_1_pix_col,event.ALPIDE_1_pix_row,histos["h_pix_occ_2D_ALPIDE_1"],masked["ALPIDE_1"])
        isnoise0 = isNoise(event.ALPIDE_0_pix_col,event.ALPIDE_0_pix_row,histos["h_pix_occ_2D_ALPIDE_0"],masked["ALPIDE_0"])
        isnoise = (isnoise2 or isnoise1 or isnoise0)
        if(isnoise): continue ### CUT!!!
        histos["h_cutflow"].Fill( cuts.index("No-noise") )
        
        ### check number of clusters per layer
        # nclusters_per_event = (event.ALPIDE_2_cls_glo_x.size()>0)+(event.ALPIDE_1_cls_glo_x.size()>0)+(event.ALPIDE_0_cls_glo_x.size()>0) ### at least one cluster
        nclusters_per_event = (event.ALPIDE_2_cls_glo_x.size()==1)+(event.ALPIDE_1_cls_glo_x.size()==1)+(event.ALPIDE_0_cls_glo_x.size()==1) ### exactly one cluster
        if(nclusters_per_event!=3): continue ### CUT!!!
        histos["h_cutflow"].Fill( cuts.index("N_{cls/det}==1") )
        histos["h_cls_3D"].Fill( event.ALPIDE_2_cls_glo_x[0],event.ALPIDE_2_cls_glo_y[0],event.ALPIDE_2_cls_glo_z[0] )
        histos["h_cls_3D"].Fill( event.ALPIDE_1_cls_glo_x[0],event.ALPIDE_1_cls_glo_y[0],event.ALPIDE_1_cls_glo_z[0] )
        histos["h_cls_3D"].Fill( event.ALPIDE_0_cls_glo_x[0],event.ALPIDE_0_cls_glo_y[0],event.ALPIDE_0_cls_glo_z[0] )
        
        ### fill cluster histos before fits
        fillClsHists("ALPIDE_0",event.ALPIDE_0_cls_npix,event.ALPIDE_0_cls_ncol,event.ALPIDE_0_cls_nrow) ### TODO: this function should get the cluster indices in the vectors
        fillClsHists("ALPIDE_1",event.ALPIDE_1_cls_npix,event.ALPIDE_1_cls_ncol,event.ALPIDE_1_cls_nrow) ### TODO: this function should get the cluster indices in the vectors
        fillClsHists("ALPIDE_2",event.ALPIDE_2_cls_npix,event.ALPIDE_2_cls_ncol,event.ALPIDE_2_cls_nrow) ### TODO: this function should get the cluster indices in the vectors
        

        ### fit
        vtx  = [xVtx,yVtx,zVtx]    if(doVtx) else []
        evtx = [exVtx,eyVtx,ezVtx] if(doVtx) else []
        
        mcx     = {"ALPIDE_2":event.ALPIDE_2_mc_loc_end_x, "ALPIDE_1":event.ALPIDE_1_mc_loc_end_x, "ALPIDE_0":event.ALPIDE_0_mc_loc_end_x}
        mcy     = {"ALPIDE_2":event.ALPIDE_2_mc_loc_end_y, "ALPIDE_1":event.ALPIDE_1_mc_loc_end_y, "ALPIDE_0":event.ALPIDE_0_mc_loc_end_y}
        nmc = 0
        for det in detectors: nmc += (mcx[det].size()>0)
        mcOK = (nmc==len(detectors))
        
        ### loop on all cluster combinations
        best_SVD  = {}
        best_Chi2 = {}
        best_value_SVD  = +1e10
        best_value_Chi2 = +1e10
        
        for i2 in range(event.ALPIDE_2_cls_glo_x.size()):
            for i1 in range(event.ALPIDE_1_cls_glo_x.size()):
                for i0 in range(event.ALPIDE_0_cls_glo_x.size()):
                    
                    clst    = {"ALPIDE_2":event.ALPIDE_2_cls_time[i2],  "ALPIDE_1":event.ALPIDE_1_cls_time[i1],  "ALPIDE_0":event.ALPIDE_0_cls_time[i0]}
                    clsx    = {"ALPIDE_2":event.ALPIDE_2_cls_glo_x[i2], "ALPIDE_1":event.ALPIDE_1_cls_glo_x[i1], "ALPIDE_0":event.ALPIDE_0_cls_glo_x[i0]}
                    clsy    = {"ALPIDE_2":event.ALPIDE_2_cls_glo_y[i2], "ALPIDE_1":event.ALPIDE_1_cls_glo_y[i1], "ALPIDE_0":event.ALPIDE_0_cls_glo_y[i0]}
                    clsz    = {"ALPIDE_2":event.ALPIDE_2_cls_glo_z[i2], "ALPIDE_1":event.ALPIDE_1_cls_glo_z[i1], "ALPIDE_0":event.ALPIDE_0_cls_glo_z[i0]}
                    clsncol = {"ALPIDE_2":event.ALPIDE_2_cls_ncol[i2],  "ALPIDE_1":event.ALPIDE_1_cls_ncol[i1],  "ALPIDE_0":event.ALPIDE_0_cls_ncol[i0]}
                    clsnrow = {"ALPIDE_2":event.ALPIDE_2_cls_nrow[i2],  "ALPIDE_1":event.ALPIDE_1_cls_nrow[i1],  "ALPIDE_0":event.ALPIDE_0_cls_nrow[i0]}
                    
                    ### to check timing #TODO
                    dt20 = abs(clst["ALPIDE_2"]-clst["ALPIDE_0"])
                    dt21 = abs(clst["ALPIDE_2"]-clst["ALPIDE_1"])
                    dt10 = abs(clst["ALPIDE_1"]-clst["ALPIDE_0"])
                    # print("dt20=",dt20,"dt21=",dt21,"dt10=",dt10)
                    
                    ### perform the SVD fit
                    points_SVD,errors_SVD = SVD_candidate(clsx,clsy,clsz,clsncol,clsnrow,vtx,evtx)
                    chisq,ndof,direction_SVD,centroid_SVD = fit_3d_SVD(points_SVD,errors_SVD)
                    chi2ndof_SVD = chisq/ndof if(ndof>0) else 99999
                    if(chi2ndof_SVD<best_value_SVD): ### happens allways!!!
                        best_value_SVD = chi2ndof_SVD
                        best_SVD.update( {"points":points_SVD} )
                        best_SVD.update( {"errors":errors_SVD} )
                        best_SVD.update( {"direction":direction_SVD} )
                        best_SVD.update( {"centroid":centroid_SVD} )
                        best_SVD.update( {"chi2ndof":chi2ndof_SVD} )
                    
                    ### perform the chi2 w/err fit
                    points_Chi2,errors_Chi2 = Chi2_candidate(clsx,clsy,clsz,clsncol,clsnrow,vtx,evtx)
                    chisq,ndof,direction_Chi2,centroid_Chi2,params_Chi2,success_Chi2 = fit_3d_chi2err(points_Chi2,errors_Chi2)
                    chi2ndof_Chi2 = chisq/ndof if(ndof>0) else 99999
                    if(success_Chi2 and chi2ndof_Chi2<best_value_Chi2): ### happens only when success_Chi2==True
                        best_value_Chi2 = chi2ndof_Chi2
                        best_Chi2.update( {"points":points_Chi2} )
                        best_Chi2.update( {"errors":errors_Chi2} )
                        best_Chi2.update( {"direction":direction_Chi2} )
                        best_Chi2.update( {"centroid":centroid_Chi2} )
                        best_Chi2.update( {"chi2ndof":chi2ndof_Chi2} )
                        best_Chi2.update( {"params":params_Chi2} )

        
        ### get the best SVD fit
        points_SVD    = best_SVD["points"]
        errors_SVD    = best_SVD["errors"]
        direction_SVD = best_SVD["direction"]
        centroid_SVD  = best_SVD["centroid"]
        chi2ndof_SVD  = best_SVD["chi2ndof"]
        plot_3d_SVD(orig_evt,points_SVD,direction_SVD,centroid_SVD,doplot)
        ### important!!!
        points = points_SVD
        ### fill some histos
        histos["h_SVDchi2"].Fill(chi2ndof_SVD)
        ### SVD track to cluster residuals
        fill_trk2cls_residuals(points,direction_SVD,centroid_SVD,"h_SVDfit_res_trk2cls")
        ### SVD track to truth residuals
        if(isMC and mcOK): fill_trk2tru_residuals(mcx,mcy,points,direction_SVD,centroid_SVD,"h_SVDfit_res_trk2tru")
        ### SVD track to vertex residuals
        if(doVtx): fill_trk2vtx_residuals(vtx,direction_SVD,centroid_SVD,"h_SVDfit_res_trk2vtx")
        
        
        ### fit successful
        passFit = (len(best_Chi2)>0)
        if(passFit):            
            ### get the best Chi2 fit
            points_Chi2    = best_Chi2["points"]
            errors_Chi2    = best_Chi2["errors"]
            direction_Chi2 = best_Chi2["direction"]
            centroid_Chi2  = best_Chi2["centroid"]
            chi2ndof_Chi2  = best_Chi2["chi2ndof"]
            params_Chi2    = best_Chi2["params"]
            plot_3d_chi2err(orig_evt,points_Chi2,params_Chi2,doplot)
            ### fill some histos
            histos["h_3Dchi2err"].Fill(chi2ndof_Chi2)
            histos["h_3Dchi2err_zoom"].Fill(chi2ndof_Chi2)
            histos["h_cutflow"].Fill( cuts.index("Chi2 Fitted") )
            dx = direction_Chi2[0]
            dy = direction_Chi2[1]
            dz = direction_Chi2[2]
            theta = np.arctan(np.sqrt(dx*dx+dy*dy)/dz)
            phi   = np.arctan(dy/dx)
            histos["h_Chi2_phi"].Fill(phi)
            histos["h_Chi2_theta"].Fill(theta)
            if(chi2ndof_Chi2<=450):
                histos["h_cutflow"].Fill( cuts.index("Fit #chi^{2}/N_{DoF}#leq450") )
                ### fill cluster histos after fit            
                fillClsHists("ALPIDE_0",event.ALPIDE_0_cls_npix,event.ALPIDE_0_cls_ncol,event.ALPIDE_0_cls_nrow,True) ### TODO: this function should get the cluster indices in the vectors
                fillClsHists("ALPIDE_1",event.ALPIDE_1_cls_npix,event.ALPIDE_1_cls_ncol,event.ALPIDE_1_cls_nrow,True) ### TODO: this function should get the cluster indices in the vectors
                fillClsHists("ALPIDE_2",event.ALPIDE_2_cls_npix,event.ALPIDE_2_cls_ncol,event.ALPIDE_2_cls_nrow,True) ### TODO: this function should get the cluster indices in the vectors
            ### Chi2 track to cluster residuals
            fill_trk2cls_residuals(points,direction_Chi2,centroid_Chi2,"h_Chi2fit_res_trk2cls")
            ### Chi2 track to truth residuals
            if(isMC and mcOK): fill_trk2tru_residuals(mcx,mcy,points,direction_Chi2,centroid_Chi2,"h_Chi2fit_res_trk2tru")
            ### Chi2 fit points on laters
            fillFitOcc(params_Chi2,"h_fit_occ_2D", "h_fit_3D")
            ### Chi2 track to vertex residuals
            Chi2_diff2vtx = 9999
            if(doVtx):
                fill_trk2vtx_residuals(vtx,direction_Chi2,centroid_Chi2,"h_Chi2fit_res_trk2vtx")
                dxv,dyv = res_track2vertex(vtx,direction_Chi2,centroid_Chi2)
                Chi2_diff2vtx = np.sqrt(dxv*dxv + dyv*dyv)            

            # if((runtype!="source" and chi2ndof_Chi2<=450) or (runtype=="source" and Chi2_diff2vtx<3.0)):
            if((runtype!="source" and chi2ndof_Chi2<=450) or (runtype=="source" and True)):
                ### test sphere for good fits only
                ## https://scikit-spatial.readthedocs.io/en/latest/gallery/intersection/plot_sphere_line.html#sphx-glr-download-gallery-intersection-plot-sphere-line-py
                x0,y0,z0 = line(rdetectors["ALPIDE_0"][2], params_Chi2) if(not doVtx) else line(0, params_Chi2)
                x1,y1,z1 = line(rdetectors["ALPIDE_2"][2], params_Chi2)
                # sphere_center = sphere_center_point
                # sphere_radius = sphere_radius_size
                sphere = Sphere(sphere_center_point,sphere_radius_size)
                trackline = Line([x0,y0,z0],direction_Chi2)
                point_a, point_b = sphere.intersect_line(trackline)
                sphere_points_a_x.append(point_a[0])
                sphere_points_a_y.append(point_a[1])
                sphere_points_a_z.append(point_a[2])
                sphere_points_b_x.append(point_b[0])
                sphere_points_b_y.append(point_b[1])
                sphere_points_b_z.append(point_b[2])
                histos["h_3Dsphere"].Fill(point_a[0],point_a[1],point_a[2])
                histos["h_3Dsphere"].Fill(point_b[0],point_b[1],point_b[2])
                histos["h_3Dsphere_a"].Fill(point_a[0],point_a[1],point_a[2])
                histos["h_3Dsphere_b"].Fill(point_b[0],point_b[1],point_b[2])
                if(doplot):
                    L1verts = getChips()
                    fig, ax = plot_3d( sphere.plotter(alpha=0.2), point_a.plotter(c='k', s=20), point_b.plotter(c='k', s=20), )
                    ax.add_collection3d(Poly3DCollection(L1verts, facecolors='green', linewidths=1, edgecolors='g', alpha=.20))
                    ax.scatter(points_Chi2[0], points_Chi2[1], points_Chi2[2], c='r', marker='o')
                    ax.plot([x0, x1], [y0, y1], [z0, z1], c='b')
                    ax.set_xlabel("x [mm]")
                    ax.set_ylabel("y [mm]")
                    ax.set_zlabel("z [mm]")
                    ax.axes.set_aspect('equal') if(not isCVMFS) else ax.axes.set_aspect('auto')
                    plt.show()
            

        ### CVR built-in fits
        if(event.trk_chi2ndof.size()>0):
            # histos["h_cutflow"].Fill( cuts.index("CVR Fitted") )
            
            for i in range(event.trk_chi2ndof.size()): histos["h_CVRchi2"].Fill( event.trk_chi2ndof[i] )

            ### centroid and direction from CVR fit
            centroid_CVR  = [event.trk_state_x[0], event.trk_state_y[0], event.trk_state_z[0]]
            direction_CVR = [event.trk_dir_x[0], event.trk_dir_y[0], event.trk_dir_z[0]]
            
            ### CVR track to cluster residuals
            fill_trk2cls_residuals(points,direction_CVR,centroid_CVR,"h_CVRfit_res_trk2cls")
            ### CVR track to truth residuals
            if(isMC and mcOK): fill_trk2tru_residuals(mcx,mcy,points,direction_CVR,centroid_CVR,"h_CVRfit_res_trk2tru")
            ### CVR track to vertex residuals
            if(doVtx): fill_trk2vtx_residuals(vtx,direction_CVR,centroid_CVR,"h_CVRfit_res_trk2vtx")

        
        ### monitoring
        if(nevents%nprintout==0 and nevents>0): print("processed events:",nevents)
        if(nmaxtoprocess>0 and nevents>=nmaxtoprocess): break
        nevents += 1
    
    
    ### post processing before writing out the histos:
    
    ### finish
    tfi.Close()
    tfo.cd()
    tfo.Write()
    tfo.Close()
    
    ### sphere plot
    if(not isCVMFS and doplot):
        sphere_center = [rdetectors["ALPIDE_1"][0],rdetectors["ALPIDE_1"][1],rdetectors["ALPIDE_1"][2]]
        sphere_radius = (rdetectors["ALPIDE_2"][2]-rdetectors["ALPIDE_0"][2])*0.7
        sphere = Sphere(sphere_center,sphere_radius)
        L1verts = getChips()
        fig, ax = plot_3d( sphere.plotter(alpha=0.2) )
        ax.add_collection3d(Poly3DCollection(L1verts, facecolors='green', linewidths=1, edgecolors='g', alpha=.20))
        ax.scatter(sphere_points_a_x, sphere_points_a_y, sphere_points_a_z, c='k', marker='o')
        ax.scatter(sphere_points_b_x, sphere_points_b_y, sphere_points_b_z, c='b', marker='o')
        ax.axes.set_aspect('equal') if(not isCVMFS) else ax.axes.set_aspect('auto')
        plt.title("Intersection of tracks", fontdict=None, loc='center', pad=None)
        ax.set_xlabel("x [mm]")
        ax.set_ylabel("y [mm]")
        ax.set_zlabel("z [mm]")
        plt.show()


#############################################################################
#############################################################################
#############################################################################

### run
outrootfilename = inrootfilename.replace("_flattree","_analysis")
Run(inrootfilename,outrootfilename)



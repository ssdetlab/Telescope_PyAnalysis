#!/usr/bin/python
import os
import os.path
import math
import time
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

import objects
from objects import *
import pixels
from pixels import *
import clusters
from clusters import *
import truth
from truth import *
import noise
from noise import *
import candidate
from candidate import *



ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(0)
# ROOT.gStyle.SetOptStat(0)

print("-----------------------------------------------------------------------------------")
print("Need to add TelescopeEvent lib and CVR libs:")
print("export LD_LIBRARY_PATH=$HOME/telescope_event:$LD_LIBRARY_PATH")
print("export LD_LIBRARY_PATH=$HOME/corryvreckan/corryvreckan-master/lib:$LD_LIBRARY_PATH")
print("-----------------------------------------------------------------------------------")

print("---- start loading libs")
### see https://root.cern/manual/python/
gInterpreter.AddIncludePath('~/telescope_event/')
gSystem.Load('libtel_event_dict.dylib')
gInterpreter.AddIncludePath('~/corryvreckan/corryvreckan-master/src/objects/')
gSystem.Load('libCorryvreckanObjects.dylib')
print("---- finish loading libs")

###############################################################
###############################################################
###############################################################
doVtx = False
runtype = "cosmics"
pdgIdMatch = 13
nmax2process = 90
doplot = False
doDiagnostics = False
doNoiseScan = False
isCVRroot = False
mc = False
cvmfs = False
reconfig(mc,cvmfs)
# tfilenamein = "~/Downloads/data_telescope/eudaq/Apr03/source_vbb6_dv15/tree_vbb6_sr90_120_Apr03_dv15.root" ## threshold is ~120e(?)
# tfilenamein = "~/Downloads/data_telescope/eudaq/Mar26/source_vbb6/tree_vbb6_src_120_Mar24.root" ## high threshold: ~400e(?)
# tfilenamein = "~/Downloads/data_telescope/Mar26/sim_out_mar_26_newstats/mc_sim_hadd_mar26_flattree.root"
# tfilenamein = "~/Downloads/data_telescope/eudaq/Mar26/cosmics_vbb6/tree.root"
# tfilenamein = "~/Downloads/data_telescope/Mar26/sim_cosmics_Apr11_thresh400e/out_corry_TelescopeRunCosmics_telescope_cosmic_mu_0_flattree.root"
# tfilenamein = "~/Downloads/data_telescope/eudaq/Apr23/cosmics_vbb6_24h_threshold7dv/tree_vbb6_cosmic_Apr10_dv7_vresetd200_clip70_run696.root"
# tfilenamein = "~/Downloads/data_telescope/Mar26/sim_out_mar_27_cosmics/corry_out_telescope_cosmic_mu_0_flattree.root"
# tfilenamein = "~/Downloads/data_telescope/Mar26/sim_src_Apr7_thresh400e/out_corry_allrun_telescope_400eThreshold_ALPIDEEField.root"
tfilenamein = "~/Downloads/data_telescope/eudaq/Apr24/source_vbb3_dv9/tree_vbb3_sr_dv9_vresetd147_clip60_run699.root"
# tfilenamein = "~/Downloads/data_telescope/eudaq/Apr25/cosmics_sim_threshold400_cvr_root/out_structured_corry_TelescopeRunCosmics_telescope_cosmic_mu_0.root"
# tfilenamein = "~/Downloads/data_telescope/eudaq/Apr25/cosmics_sim_threshold120_cvr_root/out_structured_corry_TelescopeRunCosmics_telescope_cosmic_mu_0_120e.root"



### globals
absRes  = 0.15
absChi2 = 20
if(runtype=="source"):
    absRes  *= 100
    absChi2 *= 100
histos = {}


# def fillPixOcc(det,pixels,masked):
#     for pix in pixels:
#         i = histos["h_pix_occ_2D_"+det].FindBin(pix.x,pix.y)
#         histos["h_pix_occ_1D_"+det].AddBinContent(i,1)
#         histos["h_pix_occ_2D_"+det].Fill(pix.x,pix.y)
#         if(i not in masked):
#             histos["h_pix_occ_1D_masked_"+det].AddBinContent(i,1)
#             histos["h_pix_occ_2D_masked_"+det].Fill(pix.x,pix.y)
#
#
# def fillClsHists(det,clusters,masked):
#     histos["h_ncls_"+det].Fill(len(clusters))
#     for c in clusters:
#         noisy = False
#         for pix in c.pixels:
#             i = histos["h_pix_occ_2D_"+det].FindBin(pix.x,pix.y)
#             if(i in masked):
#                 noisy = True
#                 break
#         ### not masked
#         histos["h_cls_size_"+det].Fill(len(c.pixels))
#         histos["h_cls_size_ncol_"+det].Fill(c.dx)
#         histos["h_cls_size_nrow_"+det].Fill(c.dy)
#         histos["h_cls_occ_2D_"+det].Fill(c.xmm,c.ymm)
#         if(not noisy):
#             histos["h_cls_size_masked_"+det].Fill(len(c.pixels))
#             histos["h_cls_size_ncol_masked_"+det].Fill(c.dx)
#             histos["h_cls_size_nrow_masked_"+det].Fill(c.dy)
#             histos["h_cls_occ_2D_masked_"+det].Fill(c.xmm,c.ymm)
#
#
# def fillFitOcc(params,hname2,hname3):
#     for det in detectors:
#         x,y,z = line(rdetectors[det][2],params)
#         histos[hname2+"_"+det].Fill(x,y)
#         histos[hname3].Fill(x,y,z)
#
#
# def fill_trk2cls_residuals(points,direction,centroid,hname):
#     for det in detectors:
#         dx,dy = res_track2cluster(det,points,direction,centroid)
#         histos[hname+"_x_"+det].Fill(dx)
#         histos[hname+"_y_"+det].Fill(dy)
#
#
# def fill_trk2vtx_residuals(vtx,direction,centroid,hname):
#     dxv,dyv = res_track2vertex(vtx,direction,centroid)
#     histos[hname+"_x"].Fill(dxv)
#     histos[hname+"_y"].Fill(dyv)
#
#
# def fill_trk2tru_residuals(mcparticles,pdgIdMatch,points,direction,centroid,hname):
#     for det in detectors:
#         dx,dy = res_track2truth(det,mcparticles,pdgIdMatch,points,direction,centroid)
#         # print(dy,offsets_y[det])
#         histos[hname+"_x_"+det].Fill(dx)
#         histos[hname+"_y_"+det].Fill(dy)

    
#####################################################################################
#####################################################################################
#####################################################################################


def GetTree(tfilename):
    tfile = TFile(tfilename,"READ")
    ttree = None
    if(not isMC): ttree = tfile.Get("MyTree")
    else:
        if(isCVRroot): ttree = tfile.Get("Pixel")
        else:          ttree = tfile.Get("tt")
    print("Events in tree:",ttree.GetEntries())
    if(nmax2process>0): print("Will process only",nmax2process,"events")
    return tfile,ttree


def RunNoiseScan(tfilename,tfnoisename):
    tfilenoise = TFile(tfnoisename,"RECREATE")
    tfilenoise.cd()
    h1D_noise       = {}
    h2D_noise       = {}
    for det in detectors:
        h1D_noise.update( { det:TH1D("h_noisescan_pix_occ_1D_"+det,";Pixel;Hits",npix_x*npix_y,1,npix_x*npix_y+1) } )
        h2D_noise.update( { det:TH2D("h_noisescan_pix_occ_2D_"+det,";Pixel;Hits",npix_x+1,-0.5,npix_x+0.5, npix_y+1,-0.5,npix_y+0.5) } )

    ### get the tree
    tfile,ttree = GetTree(tfilename)
    
    nprocevents = 0
    for evt in ttree:
        if(nmax2process>0 and nprocevents>nmax2process): break
        ### get the pixels
        n_active_planes,pixels = get_all_pixles(evt,h2D_noise,isCVRroot)
        for det in detectors:
            for pix in pixels[det]:
                i = h2D_noise[det].FindBin(pix.x,pix.y)
                h1D_noise[det].AddBinContent(i,1)
                h2D_noise[det].Fill(pix.x,pix.y)
        if(nprocevents%100000==0 and nprocevents>0): print("event:",nprocevents)
        nprocevents += 1
    ### finish
    tfilenoise.Write()
    tfilenoise.Close()
    print("Noise scan histos saved in:",tfnoisename)



#####################################################################################
#####################################################################################
#####################################################################################


def Run(tfilename,tfnoisename,tfo):
    ### get the tree
    tfile,ttree = GetTree(tfilename)
    
    truth_tree = None
    if(isCVRroot):
        truth_tree = tfile.Get("MCParticle")
    
    masked = GetNoiseMask(tfnoisename)
    if(isMC):
        for det in detectors:
            masked.update( {det:{}} )
    
    hPixMatix = GetPixMatrix()
    
    largest_clster = {}
    for det in detectors:
        largest_clster.update({det:Cls([],det)})
    
    nprocevents = 0
    norigevents = -1
    ientry      = 0 ### impoortant!!
    for evt in ttree:
        ### before anything else
        if(nmax2process>0 and nprocevents>nmax2process): break
        histos["h_events"].Fill(0.5)
        histos["h_cutflow"].Fill( cuts.index("All") )
        norigevents += 1
        
        ### truth particles
        mcparticles = {}
        if(isCVRroot and truth_tree is not None):
            mcparticles = get_truth_cvr(truth_tree,ientry)
            for det in detectors:
                xtru,ytru,ztru = getTruPos(det,mcparticles,pdgIdMatch)
                histos["h_tru_3D"].Fill( xtru,ytru,ztru )
                histos["h_tru_occ_2D_"+det].Fill( xtru,ytru )
        ientry += 1 ### important!
        
        ### get the pixels
        n_active_planes, pixels = get_all_pixles(evt,hPixMatix,isCVRroot)
        for det in detectors: fillPixOcc(det,pixels[det],masked[det]) ### fill pixel occupancy
        if(n_active_planes!=len(detectors)): continue ### CUT!!!
        histos["h_cutflow"].Fill( cuts.index("N_{hits/det}>0") )
        
        ### check if there's no noise
        isnoise = False
        pixels_save = {}  ### to hold a copy of all pixels
        for det in detectors:
            goodpixels = getGoodPixels(det,pixels[det],masked[det],hPixMatix[det])
            pixels[det] = goodpixels
            pixels_save.update({det:goodpixels.copy()})

        ### run clustering
        clusters = {}
        nclusters = 0
        for det in detectors:
            det_clusters = GetAllClusters(pixels[det],det)
            clusters.update( {det:det_clusters} )
            fillClsHists(det,clusters[det],masked[det])
            if(len(det_clusters)==1): nclusters += 1
        
        ### find the largest cluster
        for det in detectors:
            for c in clusters[det]:
                if(len(c.pixels)>len(largest_clster[det].pixels)): largest_clster[det] = c
        
        ### exactly one cluster per layer
        if(nclusters!=len(detectors)): continue ### CUT!!!
        histos["h_cutflow"].Fill( cuts.index("N_{cls/det}==1") )
        for det in detectors:
            histos["h_cls_3D"].Fill( clusters[det][0].xmm,clusters[det][0].ymm,clusters[det][0].zmm )

        ### diagnostics, also with truth
        if(len(mcparticles)>0 and doDiagnostics):
            for det in detectors:
                print("-------"+det+":")
                for pr in mcparticles[det]:
                    print("["+str(mcparticles[det].index(pr))+"]:",pr)
                for px in pixels_save[det]:
                    print(px)
                for cl in clusters[det]:
                    print(cl)


        # ### TODO: trying to see what is the characteristics of events with 3 single-pixel clusters alone
        # singlepixel = True
        # for det in detectors:
        #     if(len(clusters[det][0].pixels)>1):
        #         singlepixel = False
        #         break
        # if(not singlepixel): continue
        

        ### run tracking
        vtx  = [xVtx,yVtx,zVtx]    if(doVtx) else []
        evtx = [exVtx,eyVtx,ezVtx] if(doVtx) else []
        best_Chi2 = {}
        best_value_Chi2 = +1e10
        ### loop on all cluster combinations
        for i2 in range(len(clusters["ALPIDE_2"])):
            for i1 in range(len(clusters["ALPIDE_1"])):
                for i0 in range(len(clusters["ALPIDE_0"])):

                    clsx  = {"ALPIDE_2":clusters["ALPIDE_2"][i2].xmm,  "ALPIDE_1":clusters["ALPIDE_1"][i1].xmm,  "ALPIDE_0":clusters["ALPIDE_0"][i0].xmm}
                    clsy  = {"ALPIDE_2":clusters["ALPIDE_2"][i2].ymm,  "ALPIDE_1":clusters["ALPIDE_1"][i1].ymm,  "ALPIDE_0":clusters["ALPIDE_0"][i0].ymm}
                    clsz  = {"ALPIDE_2":clusters["ALPIDE_2"][i2].zmm,  "ALPIDE_1":clusters["ALPIDE_1"][i1].zmm,  "ALPIDE_0":clusters["ALPIDE_0"][i0].zmm}
                    clsdx = {"ALPIDE_2":clusters["ALPIDE_2"][i2].dxmm, "ALPIDE_1":clusters["ALPIDE_1"][i1].dxmm, "ALPIDE_0":clusters["ALPIDE_0"][i0].dxmm}
                    clsdy = {"ALPIDE_2":clusters["ALPIDE_2"][i2].dymm, "ALPIDE_1":clusters["ALPIDE_1"][i1].dymm, "ALPIDE_0":clusters["ALPIDE_0"][i0].dymm}

                    #############################
                    ### to check timing #TODO ###
                    #############################

                    points_SVD,errors_SVD = SVD_candidate(clsx,clsy,clsz,clsdx,clsdy,vtx,evtx)
                    points_Chi2,errors_Chi2 = Chi2_candidate(clsx,clsy,clsz,clsdx,clsdy,vtx,evtx)
                    chisq,ndof,direction_Chi2,centroid_Chi2,params_Chi2,success_Chi2 = fit_3d_chi2err(points_Chi2,errors_Chi2)
                    chi2ndof_Chi2 = chisq/ndof if(ndof>0) else 99999
                    if(success_Chi2 and chi2ndof_Chi2<best_value_Chi2): ### happens only when success_Chi2==True
                        best_value_Chi2 = chi2ndof_Chi2
                        best_Chi2.update( {"svd_points":points_SVD} )
                        best_Chi2.update( {"points":points_Chi2} )
                        best_Chi2.update( {"errors":errors_Chi2} )
                        best_Chi2.update( {"direction":direction_Chi2} )
                        best_Chi2.update( {"centroid":centroid_Chi2} )
                        best_Chi2.update( {"chi2ndof":chi2ndof_Chi2} )
                        best_Chi2.update( {"params":params_Chi2} )
        
        ### fit successful
        passFit = (len(best_Chi2)>0)
        if(passFit):
            ### get the best Chi2 fit
            points_SVD     = best_Chi2["svd_points"]
            points_Chi2    = best_Chi2["points"]
            errors_Chi2    = best_Chi2["errors"]
            direction_Chi2 = best_Chi2["direction"]
            centroid_Chi2  = best_Chi2["centroid"]
            chi2ndof_Chi2  = best_Chi2["chi2ndof"]
            params_Chi2    = best_Chi2["params"]
            plot_3d_chi2err(norigevents,points_Chi2,params_Chi2,doplot)
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
            if(abs(np.sin(theta))>1e-10): histos["h_Chi2_theta_weighted"].Fill( theta,abs(1/(2*np.pi*np.sin(theta))) )
            if(chi2ndof_Chi2<=450): histos["h_cutflow"].Fill( cuts.index("Fit #chi^{2}/N_{DoF}#leq450") )
            ### Chi2 track to cluster residuals
            fill_trk2cls_residuals(points_SVD,direction_Chi2,centroid_Chi2,"h_Chi2fit_res_trk2cls",histos)
            ### Chi2 track to truth residuals
            if(isMC): fill_trk2tru_residuals(mcparticles,pdgIdMatch,points_SVD,direction_Chi2,centroid_Chi2,"h_Chi2fit_res_trk2tru",histos)
            ### Chi2 fit points on laters
            fillFitOcc(params_Chi2,"h_fit_occ_2D", "h_fit_3D",histos)
            ### Chi2 track to vertex residuals
            if(doVtx): fill_trk2vtx_residuals(vtx,direction_Chi2,centroid_Chi2,"h_Chi2fit_res_trk2vtx",histos)

            ### fill cluster size vs true position
            if(isCVRroot and truth_tree is not None):
                for det in detectors:
                    xtru,ytru,ztru = getTruPos(det,mcparticles,pdgIdMatch)
                    wgt = clusters[det][0].n
                    posx = ((xtru-pix_x/2.)%(2*pix_x))
                    posy = ((ytru-pix_y/2.)%(2*pix_y))
                    histos["h_csize_vs_trupos"].Fill(posx,posy,wgt)
                    histos["h_ntrks_vs_trupos"].Fill(posx,posy)
                    histos["h_csize_vs_trupos_"+det].Fill(posx,posy,wgt)
                    histos["h_ntrks_vs_trupos_"+det].Fill(posx,posy)
                    ### divide into smaller sizes
                    strcsize = str(wgt) if(wgt<5) else "n"
                    histos["h_csize_"+strcsize+"_vs_trupos"].Fill(posx,posy,wgt)
                    histos["h_ntrks_"+strcsize+"_vs_trupos"].Fill(posx,posy)
                    histos["h_csize_"+strcsize+"_vs_trupos_"+det].Fill(posx,posy,wgt)
                    histos["h_ntrks_"+strcsize+"_vs_trupos_"+det].Fill(posx,posy)
                
                    # if(det=="ALPIDE_0"): print("Size:",wgt,"Tru:",xtru,ytru,"Residuals:",(xtru%pix_x),(ytru%pix_y))
        
        ### event counter
        if(nprocevents%10==0 and nprocevents>0): print("processed event:",nprocevents,"out of",norigevents,"events read")
        nprocevents += 1


    #######################
    ### post processing ###
    #######################
    
    
    ### cluster mean size vs position
    tfo.cd()
    hname = "h_csize_vs_trupos"
    hnewname = hname.replace("csize","mean")
    hdenname = hname.replace("csize","ntrks")
    histos.update( {hnewname:histos[hname].Clone(hnewname)} )
    histos[hnewname].Divide(histos[hdenname])
    for det in detectors:
        tfo.cd(det)
        hname = "h_csize_vs_trupos_"+det
        hnewname = hname.replace("csize","mean")
        hdenname = hname.replace("csize","ntrks")
        histos.update( {hnewname:histos[hname].Clone(hnewname)} )
        histos[hnewname].Divide(histos[hdenname])
    for j in range(1,6):
        tfo.cd()
        strcsize = str(j) if(j<5) else "n"
        hname = "h_csize_"+strcsize+"_vs_trupos"
        hnewname = hname.replace("csize","mean")
        hdenname = hname.replace("csize","ntrks")
        histos.update( {hnewname:histos[hname].Clone(hnewname)} )
        histos[hnewname].Divide(histos[hdenname])
        for det in detectors:
            tfo.cd(det)
            hname = "h_csize_"+strcsize+"_vs_trupos_"+det
            hnewname = hname.replace("csize","mean")
            hdenname = hname.replace("csize","ntrks")
            histos.update( {hnewname:histos[hname].Clone(hnewname)} )
            histos[hnewname].Divide(histos[hdenname])
    
    ### largest clusters
    for det in detectors:    
        for pix in largest_clster[det].pixels:
            histos["h_big_cls_2D_"+det].Fill(pix.x,pix.y)
        

#############################################################################
#############################################################################
#############################################################################

# get the start time
st = time.time()


tfnoisename = tfilenamein.replace(".root","_noise.root")
isnoisefile = os.path.isfile(os.path.expanduser(tfnoisename))
print("Running on:",tfilenamein)
if(doNoiseScan):
    print("Noise run file exists?:",isnoisefile)
    if(isnoisefile):
        redonoise = input("Noise file exists - do you want to rederive it?[y/n]:")
        if(redonoise=="y" or redonoise=="Y"):
            RunNoiseScan(tfilenamein,tfnoisename)
            masked = GetNoiseMask(tfnoisename)
        else:
            print("Option not understood - please try again.")
    else:
        RunNoiseScan(tfilenamein,tfnoisename)
        masked = GetNoiseMask(tfnoisename)
    quit()
else:
    if(not isnoisefile):
        print("Noise file",tfnoisename,"not found")
        print("Generate first by setting doNoiseScan=True")
        quit()

tfilenameout = tfilenamein.replace(".root","_histograms.root")
tfo = TFile(tfilenameout,"RECREATE")
tfo.cd()
# book_histos(tfo,absRes,absChi2)
histos = book_histos(absRes,absChi2,tfo)
Run(tfilenamein,tfnoisename,tfo)
tfo.cd()
tfo.Write()
tfo.Close()

# get the end time
et = time.time()
# get the execution time
elapsed_time = et - st
print('Execution time:', elapsed_time, 'seconds')




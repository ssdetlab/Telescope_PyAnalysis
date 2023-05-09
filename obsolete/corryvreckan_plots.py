#!/usr/bin/python
import os
import math
import subprocess
import array
import numpy as np
import ROOT
from ROOT import *

import config
from config import *


import argparse
parser = argparse.ArgumentParser(description='corryvreckan_analysis.py...')
parser.add_argument('-fd', metavar='input data file', required=True,  help='full path to input data file')
parser.add_argument('-fs', metavar='input simulation file', required=True,  help='full path to input simulation file')
argus = parser.parse_args()
inrootfilename_data = argus.fd
inrootfilename_simu = argus.fs

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPadBottomMargin(0.15)
ROOT.gStyle.SetPadLeftMargin(0.13)
ROOT.gStyle.SetPadRightMargin(0.13)

### histograms
histos = {}

#############################################################################
#############################################################################
#############################################################################

hnames_overall   = ["h_events","h_cutflow","h_SVDchi2","h_3Dchi2err","h_CVRchi2","h_npix","h_3Dchi2err_zoom",
                    "h_SVDfit_res_trk2vtx_x","h_SVDfit_res_trk2vtx_y","h_Chi2fit_res_trk2vtx_x","h_Chi2fit_res_trk2vtx_y", "h_CVRfit_res_trk2vtx_x", "h_CVRfit_res_trk2vtx_y",
                    "h_Chi2_phi","h_Chi2_theta","h_3Dsphere","h_3Dsphere_a","h_3Dsphere_b","h_fit_3D" ]
                        
hnames_detectors = ["h_pix_occ_1D","h_pix_occ_1D_masked","h_pix_occ_2D","h_pix_occ_2D_masked","h_fit_occ_2D",
                    "h_cls_occ_2D","h_cls_occ_2D_masked",
                    "h_SVDfit_res_trk2cls_x","h_SVDfit_res_trk2cls_y","h_Chi2fit_res_trk2cls_x","h_Chi2fit_res_trk2cls_y","h_CVRfit_res_trk2cls_x","h_CVRfit_res_trk2cls_y",
                    "h_SVDfit_res_trk2tru_x","h_SVDfit_res_trk2tru_y","h_Chi2fit_res_trk2tru_x","h_Chi2fit_res_trk2tru_y","h_CVRfit_res_trk2tru_x","h_CVRfit_res_trk2tru_y",
                    "h_ncls","h_cls_size","h_cls_size_ncol","h_cls_size_nrow",
                    # "h_ncls_postfit","h_cls_size_postfit","h_cls_size_ncol_postfit","h_cls_size_nrow_postfit",
                    "h_pix_occ_1D","h_pix_occ_1D_masked",
                    ]


### book histos
def book_histos(fdataname,fsimuname):
    
    fdata = TFile(fdataname,"READ")
    fdata_tree = TFile(fdataname.replace("_analysis","_flattree"),"READ")
    fsimu = TFile(fsimuname,"READ")
    
    files = {"data":fdata,"simu":fsimu}
    
    for name,f in files.items():
        for hname in hnames_overall:
            histos.update( { hname+"_"+name : f.Get(hname).Clone(hname+"_"+name) } )
            histos[hname+"_"+name].SetDirectory(0)
            
    
        for hname in hnames_detectors:
            for det in detectors:
                if(hname=="h_pix_occ_1D" and "eudaq" not in fdataname):
                    hname1 = "h_pixocc1D"
                    histos.update( { hname+"_"+det+"_"+name : fdata_tree.Get(hname1+"_"+det).Clone(hname+"_"+det+"_"+name) } )
                else:
                    histos.update( { hname+"_"+det+"_"+name : f.Get(hname+"_"+det).Clone(hname+"_"+det+"_"+name) } )
                histos[hname+"_"+det+"_"+name].SetDirectory(0)


def plotdatamasking2(hname,hnamemasked,det,pdfname,logy=True):
    histos[hname+"_data"].SetLineColor(ROOT.kBlack)
    histos[hname+"_data"].SetLineWidth(1)
    histos[hnamemasked+"_data"].SetLineColor(ROOT.kBlack)
    histos[hnamemasked+"_data"].SetLineWidth(1)
    cnv = TCanvas("cnv","",1400,500)
    cnv.Divide(2,1)
    cnv.cd(1)
    gPad.SetTicks(1,1)
    if(logy): gPad.SetLogy()
    histos[hname+"_data"].Draw("hist")
    cnv.cd(2)
    gPad.SetTicks(1,1)
    # if(logy): gPad.SetLogy()
    histos[hnamemasked+"_data"].Draw("hist")
    cnv.SaveAs(pdfname.replace(".pdf","")+"_"+hname+".pdf")
    cnv.SaveAs(pdfname)


def overlay2(hname,det,pdfname,norm=True,logy=True):
    histos[hname+"_data"].SetLineColor(ROOT.kBlack)
    histos[hname+"_simu"].SetLineColor(ROOT.kRed)
    histos[hname+"_data"].SetMarkerColor(ROOT.kBlack)
    histos[hname+"_simu"].SetMarkerColor(ROOT.kRed)
    histos[hname+"_data"].SetMarkerSize(1)
    histos[hname+"_simu"].SetMarkerSize(1)
    histos[hname+"_data"].SetMarkerStyle(20)
    histos[hname+"_simu"].SetMarkerStyle(24)
    
    leg = TLegend(0.60,0.70,0.85,0.87) if(det!="") else TLegend(0.70,0.70,0.85,0.87)
    leg.SetFillStyle(4000) # will be transparent
    leg.SetFillColor(0)
    leg.SetTextFont(42)
    leg.SetBorderSize(0)
    leg.AddEntry(histos[hname+"_data"],det+" Data" if(det!="") else "Data","lp")
    leg.AddEntry(histos[hname+"_simu"],det+" Sim." if(det!="") else "Sim.","lp")
    
    if(norm):
        histos[hname+"_simu"].Scale(1./histos[hname+"_simu"].Integral())
        histos[hname+"_data"].Scale(1./histos[hname+"_data"].Integral())
    ymax_simu = histos[hname+"_simu"].GetMaximum()
    ymax_data = histos[hname+"_data"].GetMaximum()
    ymax = ymax_simu if(ymax_simu>ymax_data) else ymax_data
    histos[hname+"_simu"].SetMaximum(ymax*1.5)
    histos[hname+"_data"].SetMaximum(ymax*1.5)
    
    cnv = TCanvas("cnv","",1000,500)
    cnv.SetTicks(1,1)
    if(logy): cnv.SetLogy()
    # histos[hname+"_simu"].DrawNormalized("e1p") if(norm) else histos[hname+"_simu"].Draw("e1p")
    histos[hname+"_simu"].Draw("e1p")
    # histos[hname+"_data"].DrawNormalized("e1p same") if(norm) else histos[hname+"_data"].Draw("e1p same")
    histos[hname+"_data"].Draw("e1p same")
    leg.Draw("same")
    cnv.SaveAs(pdfname.replace(".pdf","")+"_"+hname+".pdf")
    cnv.SaveAs(pdfname)


def overlay2_xy(hname,det,pdfname,norm=True):
    name = "_"+det if(det!="") else ""
    
    histos[hname+"_x"+name+"_data"].SetLineColor(ROOT.kBlack)
    histos[hname+"_x"+name+"_simu"].SetLineColor(ROOT.kRed)
    histos[hname+"_x"+name+"_data"].SetMarkerColor(ROOT.kBlack)
    histos[hname+"_x"+name+"_simu"].SetMarkerColor(ROOT.kRed)
    histos[hname+"_x"+name+"_data"].SetMarkerSize(1)
    histos[hname+"_x"+name+"_simu"].SetMarkerSize(1)
    histos[hname+"_x"+name+"_data"].SetMarkerStyle(20)
    histos[hname+"_x"+name+"_simu"].SetMarkerStyle(24)
    
    histos[hname+"_y"+name+"_data"].SetLineColor(ROOT.kBlack)
    histos[hname+"_y"+name+"_simu"].SetLineColor(ROOT.kRed)
    histos[hname+"_y"+name+"_data"].SetMarkerColor(ROOT.kBlack)
    histos[hname+"_y"+name+"_simu"].SetMarkerColor(ROOT.kRed)
    histos[hname+"_y"+name+"_data"].SetMarkerSize(1)
    histos[hname+"_y"+name+"_simu"].SetMarkerSize(1)
    histos[hname+"_y"+name+"_data"].SetMarkerStyle(20)
    histos[hname+"_y"+name+"_simu"].SetMarkerStyle(24)
    
    leg = TLegend(0.55,0.70,0.85,0.87)
    leg.SetFillStyle(4000) # will be transparent
    leg.SetFillColor(0)
    leg.SetTextFont(42)
    leg.SetBorderSize(0)
    leg.AddEntry(histos[hname+"_x"+name+"_data"],det+" Data","lp")
    leg.AddEntry(histos[hname+"_x"+name+"_simu"],det+" Sim.","lp")
    
    if(norm):
        histos[hname+"_x"+name+"_simu"].Scale(1./histos[hname+"_x"+name+"_simu"].Integral())
        histos[hname+"_y"+name+"_simu"].Scale(1./histos[hname+"_y"+name+"_simu"].Integral())
        histos[hname+"_x"+name+"_data"].Scale(1./histos[hname+"_x"+name+"_data"].Integral())
        histos[hname+"_y"+name+"_data"].Scale(1./histos[hname+"_y"+name+"_data"].Integral())
    ymax_simu = histos[hname+"_x"+name+"_simu"].GetMaximum()
    ymax_data = histos[hname+"_x"+name+"_data"].GetMaximum()
    ymax = ymax_simu if(ymax_simu>ymax_data) else ymax_data
    histos[hname+"_x"+name+"_simu"].SetMaximum(ymax*1.5)
    histos[hname+"_x"+name+"_data"].SetMaximum(ymax*1.5)
    ymax_simu = histos[hname+"_y"+name+"_simu"].GetMaximum()
    ymax_data = histos[hname+"_y"+name+"_data"].GetMaximum()
    ymax = ymax_simu if(ymax_simu>ymax_data) else ymax_data
    histos[hname+"_y"+name+"_simu"].SetMaximum(ymax*1.5)
    histos[hname+"_y"+name+"_data"].SetMaximum(ymax*1.5)
    
    cnv = TCanvas("cnv","",1400,500)
    cnv.Divide(2,1)
    cnv.cd(1)
    gPad.SetTicks(1,1)
    # histos[hname+"_x"+name+"_simu"].DrawNormalized("e1p")      if(norm) else histos[hname+"_x"+name+"_simu"].Draw("e1p")
    histos[hname+"_x"+name+"_simu"].Draw("e1p")
    # histos[hname+"_x"+name+"_data"].DrawNormalized("e1p same") if(norm) else histos[hname+"_x"+name+"_data"].Draw("e1p same")
    histos[hname+"_x"+name+"_data"].Draw("e1p same")
    leg.Draw("same")
    cnv.cd(2)
    gPad.SetTicks(1,1)
    # histos[hname+"_y"+name+"_simu"].DrawNormalized("e1p")      if(norm) else histos[hname+"_y"+name+"_simu"].Draw("e1p")
    histos[hname+"_y"+name+"_simu"].Draw("e1p")
    # histos[hname+"_y"+name+"_data"].DrawNormalized("e1p same") if(norm) else histos[hname+"_y"+name+"_data"].Draw("e1p same")
    histos[hname+"_y"+name+"_data"].Draw("e1p same")
    leg.Draw("same")
    
    cnv.SaveAs(pdfname.replace(".pdf","")+"_"+hname+name+".pdf")
    cnv.SaveAs(pdfname)


def plot2_3D(hname,pdfname):
    cnv = TCanvas("cnv","",1200,1200)
    cnv.Divide(2,1)
    cnv.cd(1)
    gPad.SetTicks(1,1)
    histos[hname+"_data"].Draw("BOX2")
    cnv.cd(2)
    gPad.SetTicks(1,1)
    histos[hname+"_simu"].Draw("BOX2")
    cnv.SaveAs(pdfname.replace(".pdf","")+"_"+hname+".pdf")
    cnv.SaveAs(pdfname)


def plot2_2D(det,hname,pdfname,axestitle=""):
    stitle = det+" Sim."
    dtitle = det+" Data"
    if(axestitle!=""):
        stitle = stitle+" "+axestitle
        dtitle = dtitle+" "+axestitle
    cnv = TCanvas("cnv","",1400,500)
    cnv.Divide(2,1)
    cnv.cd(1)
    gPad.SetTicks(1,1)
    histos[hname+"_data"].SetTitle(dtitle)
    histos[hname+"_data"].Draw("col")
    cnv.cd(2)
    gPad.SetTicks(1,1)
    histos[hname+"_simu"].SetTitle(stitle)
    histos[hname+"_simu"].Draw("col")
    cnv.SaveAs(pdfname.replace(".pdf","")+"_"+hname+".pdf")
    cnv.SaveAs(pdfname)


### start running
book_histos(inrootfilename_data,inrootfilename_simu)
# print(histos)

dataname = inrootfilename_data.split("/")[-1].replace(".root","")
simuname = inrootfilename_simu.split("/")[-1].replace(".root","")
suffix = "data-"+dataname+"_simul-"+simuname
pdfname = "comparison_"+suffix+".pdf"

### start plotting:
cnv = TCanvas("cnv","",500,500)
cnv.SaveAs(pdfname+"(")


plot2_3D("h_fit_3D",pdfname)
for det in detectors: plotdatamasking2("h_pix_occ_1D_"+det,"h_pix_occ_1D_masked_"+det,det,pdfname)
histos["h_cutflow_data"].Scale(1./histos["h_cutflow_data"].GetBinContent(3))
histos["h_cutflow_simu"].Scale(1./histos["h_cutflow_simu"].GetBinContent(3))
histos["h_cutflow_data"].GetYaxis().SetTitle("Events normalize")
histos["h_cutflow_simu"].GetYaxis().SetTitle("Events normalize")
overlay2("h_cutflow","",pdfname,False,True)
histos["h_3Dchi2err_data"].Rebin(5)
histos["h_3Dchi2err_simu"].Rebin(5)
overlay2("h_3Dchi2err","",pdfname,True,False)
if("h_3Dchi2err_zoom" in hnames_overall):
    histos["h_3Dchi2err_zoom_data"].Rebin(2)
    histos["h_3Dchi2err_zoom_simu"].Rebin(2)
    overlay2("h_3Dchi2err_zoom","",pdfname,True,False)
# overlay2("h_CVRchi2","",pdfname,True,True)
overlay2("h_Chi2_phi","",pdfname,True,False)
overlay2("h_Chi2_theta","",pdfname,True,False)
for det in detectors: plot2_2D(det,"h_pix_occ_2D_masked_"+det,pdfname,"Hit occupancy;Pixel column index;Pixel row index;Hits")
for det in detectors: plot2_2D(det,"h_cls_occ_2D_masked_"+det,pdfname,"Cluster occupancy;x [mm];y [mm];Clusters")
for det in detectors: plot2_2D(det,"h_fit_occ_2D_"+det,pdfname,"Track occupancy;x [mm];y [mm];Tracks")
# for det in detectors: overlay2("h_ncls_"+det,det,pdfname)
for det in detectors: overlay2("h_cls_size_"+det,det,pdfname,True,True)
for det in detectors: overlay2("h_cls_size_ncol_"+det,det,pdfname,True,True)
for det in detectors: overlay2("h_cls_size_nrow_"+det,det,pdfname,True,True)
for det in detectors: overlay2_xy("h_Chi2fit_res_trk2cls",det,pdfname)
if("cosmics" not in inrootfilename_data and "cosmics" not in inrootfilename_simu): overlay2_xy("h_Chi2fit_res_trk2vtx","",pdfname)

### end plotting
cnv = TCanvas("cnv","",500,500)
cnv.SaveAs(pdfname+")")

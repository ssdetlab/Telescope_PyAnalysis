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

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetOptStat(0)

cols = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2, ]
files = {
    "run74x":"~/Downloads/data_telescope/eudaq/Jun05/vbb6_dv10_vresetd200_clip70_run74x/tree_vbb6_dv10_vresetd200_clip70_run74x_multiprocess_histograms.root", ##delay=4.0us, strobe=100ns
    # "run75x":"~/Downloads/data_telescope/eudaq/Jun12/vbb6_dv10_vresetd200_clip70_run75x/tree_vbb6_dv10_vresetd200_clip70_run75x_multiprocess_histograms.root", ##delay=4.7us, strobe=100ns
    # "run756":"~/Downloads/data_telescope/eudaq/Jun17/vbb6_dv10_vresetd200_clip70_run756/tree_vbb6_dv10_vresetd200_clip70_run756_multiprocess_histograms.root", ##delay=165ns, strobe=12us
    # "run759":"~/Downloads/data_telescope/eudaq/Jun18/vbb6_dv10_vresetd200_clip70_run759/tree_vbb6_dv10_vresetd200_clip70_run759_multiprocess_histograms.root", ##delay=165ns, strobe=100us
}
detectors = ["ALPIDE_0", "ALPIDE_1", "ALPIDE_2", "ALPIDE_3"]
histprefx = ["h_cls_size", "h_Chi2fit_res_trk2cls_x", "h_Chi2fit_res_trk2cls_y", ]
# histnames = ["h_cls_size_ALPIDE_0", "h_cls_size_ALPIDE_1", "h_cls_size_ALPIDE_2", "h_cls_size_ALPIDE_3",
#              "h_Chi2fit_res_trk2cls_x_ALPIDE_0", "h_Chi2fit_res_trk2cls_x_ALPIDE_1", "h_Chi2fit_res_trk2cls_x_ALPIDE_2", "h_Chi2fit_res_trk2cls_x_ALPIDE_3",
#              "h_Chi2fit_res_trk2cls_y_ALPIDE_0", "h_Chi2fit_res_trk2cls_y_ALPIDE_1", "h_Chi2fit_res_trk2cls_y_ALPIDE_2", "h_Chi2fit_res_trk2cls_y_ALPIDE_3",
#             ]
histos = {}
runs = []
runscol = {}
for run,fname in files.items():
    runs.append(run)
    runscol.update({run:cols[runs.index(run)]})



def book_histos(tfo):
    tfo.cd()
    for run,fname in files.items():
        for prefx in histprefx:
            for det in detectors:
                hname = prefx+"_"+det
                hist = det+"/"+hname
                name = run+"_"+hname
                tfi = TFile(fname,"READ")
                histos.update({name:tfi.Get(hist).Clone(name)})
                histos[name].SetDirectory(0)


def write_histos(tfo):
    tfo.cd()
    for hname,hist in histos.items():
        hist.Write()


def plot_2x2_histos(pdf,prefix):
    ymax = -1e10
    for run in runs:
        for det in detectors:
            hname = prefix+"_"+det
            histos[run+"_"+hname].SetLineColor(runscol[run])
            histos[run+"_"+hname].SetMarkerColor(runscol[run])
            histos[run+"_"+hname].SetMarkerSize(1)
            histos[run+"_"+hname].SetMarkerStyle(20)
            histos[run+"_"+hname].Scale(1./histos[run+"_"+hname].Integral())
            histos[run+"_"+hname].SetTitle(hname.replace("h_cls_size_",""))
            histos[run+"_"+hname].GetYaxis().SetTitle("Normalized")
            tmax = histos[run+"_"+hname].GetMaximum()
            ymax = tmax if(tmax>ymax) else ymax

    for run in runs:
        for det in detectors:
            hname = prefix+"_"+det
            histos[run+"_"+hname].SetMaximum(ymax*1.5)
    
    leg = TLegend(0.65,0.60,0.85,0.87)
    leg.SetFillStyle(4000) # will be transparent
    leg.SetFillColor(0)
    leg.SetTextFont(42)
    leg.SetBorderSize(0)
    for run in runs:
        leg.AddEntry(histos[run+"_"+prefix+"_ALPIDE_0"],run,"lp")
    
    cnv = TCanvas("cnv","",1200,1000)
    cnv.Divide(2,2)
    for count1,det in enumerate(detectors):
        p = cnv.cd(count1+1)
        p.SetTicks(1,1)
        hname = prefix+"_"+det
        for count2,run in enumerate(runs):
            if(count2==0): histos[run+"_"+hname].Draw("e1p")
            else:          histos[run+"_"+hname].Draw("e1p same")
        leg.Draw("same")
    cnv.SaveAs(pdf)
    
#####################################################################################
#####################################################################################
#####################################################################################

tfilenameout = "compare.root"
tfo = TFile(tfilenameout,"RECREATE")
book_histos(tfo)

plot_2x2_histos(tfilenameout.replace("root","pdf("),"h_cls_size")
plot_2x2_histos(tfilenameout.replace("root","pdf"),"h_Chi2fit_res_trk2cls_x")
plot_2x2_histos(tfilenameout.replace("root","pdf)"),"h_Chi2fit_res_trk2cls_y")

write_histos(tfo)
tfo.Write()
tfo.Close()
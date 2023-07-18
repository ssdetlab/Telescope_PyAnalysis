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

import argparse
parser = argparse.ArgumentParser(description='serial_analyzer.py...')
parser.add_argument('-run', metavar='run type', required=True,  help='run type [source/cosmics]')
argus = parser.parse_args()
runtype = argus.run

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetOptStat(0)

mm2um = 1000
cols = [ROOT.kBlack, ROOT.kRed, ROOT.kGreen+2, ROOT.kOrange, ROOT.kBlue]
mrks = [20,          24,        32,         25,            23]
files = {
    # "run74x":"~/Downloads/data_telescope/eudaq/Jun05/vbb6_dv10_vresetd200_clip70_run74x/tree_vbb6_dv10_vresetd200_clip70_run74x_multiprocess_histograms.root", ##delay=4.0us, strobe=100ns
    # "run75x":"~/Downloads/data_telescope/eudaq/Jun12/vbb6_dv10_vresetd200_clip70_run75x/tree_vbb6_dv10_vresetd200_clip70_run75x_multiprocess_histograms.root", ##delay=4.7us, strobe=100ns
    # "run75y":"~/Downloads/data_telescope/eudaq/Jun17/vbb6_dv10_vresetd200_clip70_run75y/tree_vbb6_dv10_vresetd200_clip70_run75y_multiprocess_histograms.root", ##delay=165ns, strobe=12us
    # "run76x":"~/Downloads/data_telescope/eudaq/Jun27/vbb6_dv10_vresetd200_clip70_run76x/tree_vbb6_dv10_vresetd200_clip70_run76x_multiprocess_histograms.root", ##delay=1.5us, strobe=10us
    # "run760":"~/Downloads/data_telescope/eudaq/Jun22/vbb6_dv10_vresetd200_clip70_run760/tree_vbb6_dv10_vresetd200_clip70_run760_multiprocess_histograms.root", ##delay=165ns, strobe=100ns
    # # "run759":"~/Downloads/data_telescope/eudaq/Jun18/vbb6_dv10_vresetd200_clip70_run759/tree_vbb6_dv10_vresetd200_clip70_run759_multiprocess_histograms.root", ##delay=165ns, strobe=100us
    
    "run77x":"~/Downloads/data_telescope/eudaq/Jul08/vbb0_dv10_vresetd147_clip60_run77x/tree_vbb0_dv10_vresetd147_clip60_run77x_multiprocess_histograms.root",
    "run77y":"~/Downloads/data_telescope/eudaq/Jul12/vbb0_dv10_vresetd147_clip60_run77y/tree_vbb0_dv10_vresetd147_clip60_run77y_multiprocess_histograms.root",
}
dely = {
    "run74x":"4.0 #mus",
    "run75x":"4.7 #mus",
    "run75y":"165 ns",
    "run76x":"1.5 #mus",
    "run760":"165 ns",
    "run77x":"150 ns",
    "run77y":"150 ns",
}
strb = {
    "run74x":"100 ns",
    "run75x":"100 ns",
    "run75y":"12 #mus",
    "run76x":"10 #mus",
    "run760":"100 ns",
    "run77x":"10 #mus",
    "run77y":"10 #mus",
}
detectors = ["ALPIDE_0", "ALPIDE_1", "ALPIDE_2"] 
if(runtype=="cosmics"): detectors.append("ALPIDE_3")
histprefx = ["h_cls_size", "h_Chi2fit_res_trk2cls_x", "h_Chi2fit_res_trk2cls_y", ]
histos = {}
runs = []
runscol = {}
runsmrk = {}
for run,fname in files.items():
    runs.append(run)
    runscol.update({run:cols[runs.index(run)]})
    runsmrk.update({run:mrks[runs.index(run)]})
run2fit = "run75y"



def book_histos(tfo):
    tfo.cd()
    for run,fname in files.items():
        for prefx in histprefx:
            for det in detectors:
                hname = prefx+"_"+det
                hist = det+"/"+hname
                name = run+"_"+hname
                tfi = TFile(fname,"READ")
                print("From file:",fname,"getting histogram named:",hist)
                histos.update({name:tfi.Get(hist).Clone(name)})
                if(det in histos[name].GetTitle()): histos[name].SetTitle( det )
                histos[name].SetDirectory(0)


def write_histos(tfo):
    tfo.cd()
    for hname,hist in histos.items():
        hist.Write()

def fit1(h,col,xmin,xmax):
    g1 = TF1("g1", "gaus", xmin,xmax)
    # f1 = TF1("f1", "gaus(0)", xmin,xmax)
    g1.SetLineColor(col)
    # f1.SetLineColor(col)
    h.Fit(g1,"EMRS")
    # f1.SetParameter(0,g1.GetParameter(0))
    # f1.SetParameter(1,g1.GetParameter(1))
    # f1.SetParameter(2,g1.GetParameter(2))
    # chi2dof = f1.GetChisquare()/f1.GetNDF() if(f1.GetNDF()>0) else -1
    chi2dof = g1.GetChisquare()/g1.GetNDF() if(g1.GetNDF()>0) else -1
    print("g1 chi2/Ndof=",chi2dof)
    return g1

def fit2(h,col):
    g1 = TF1("g1", "gaus", xmin,xmax)
    g2 = TF1("g2", "gaus", xmin,xmax)
    f1 = TF1("f2", "gaus(0)+gaus(3)", xmin,xmax)
    g1.SetLineColor(col)
    g2.SetLineColor(col)
    f2.SetLineColor(col)
    h.Fit(g1,"EMRS")
    h.Fit(g2,"EMRS")
    f2.SetParameter(0,g1.GetParameter(0))
    f2.SetParameter(1,g1.GetParameter(1))
    f2.SetParameter(2,g1.GetParameter(2))
    f2.SetParameter(3,g2.GetParameter(0))
    f2.SetParameter(4,g2.GetParameter(1))
    f2.SetParameter(5,g2.GetParameter(2))
    chi2dof = f2.GetChisquare()/f2.GetNDF() if(f2.GetNDF()>0) else -1
    print("f2 chi2/Ndof=",chi2dof)
    return f2

def fit3(h,col):
    g1 = TF1("g1", "gaus", xmin,xmax)
    g2 = TF1("g2", "gaus", xmin,xmax)
    g3 = TF1("g3", "gaus", xmin,xmax)
    f3 = TF1("f3", "gaus(0)+gaus(3)+gaus(6)", xmin,xmax)
    g1.SetLineColor(col)
    g2.SetLineColor(col)
    g3.SetLineColor(col)
    f3.SetLineColor(col)
    h.Fit(g1,"EMRS")
    h.Fit(g2,"EMRS")
    h.Fit(g3,"EMRS")
    f3.SetParameter(0,g1.GetParameter(0))
    f3.SetParameter(1,g1.GetParameter(1))
    f3.SetParameter(2,g1.GetParameter(2))
    f3.SetParameter(3,g2.GetParameter(0))
    f3.SetParameter(4,g2.GetParameter(1))
    f3.SetParameter(5,g2.GetParameter(2))
    f3.SetParameter(6,g3.GetParameter(0))
    f3.SetParameter(7,g3.GetParameter(1))
    f3.SetParameter(8,g3.GetParameter(2))
    chi2dof = f3.GetChisquare()/f3.GetNDF() if(f3.GetNDF()>0) else -1
    print("f3 chi2/Ndof=",chi2dof)
    return f3


def plot_2x2_histos(pdf,prefix):
    ymax = -1e10
    for run in runs:
        for det in detectors:
            hname = prefix+"_"+det
            histos[run+"_"+hname].SetLineColor(runscol[run])
            histos[run+"_"+hname].SetMarkerColor(runscol[run])
            histos[run+"_"+hname].SetMarkerSize(1)
            histos[run+"_"+hname].SetMarkerStyle(runsmrk[run])
            histos[run+"_"+hname].Scale(1./histos[run+"_"+hname].Integral())
            histos[run+"_"+hname].SetTitle(det)
            histos[run+"_"+hname].GetYaxis().SetTitle("Normalized")
            tmax = histos[run+"_"+hname].GetMaximum()
            ymax = tmax if(tmax>ymax) else ymax

    # factor = 2 if("cls_size" in prefix) else 1.2
    factor = 1.2
    for run in runs:
        for det in detectors:
            hname = prefix+"_"+det
            histos[run+"_"+hname].SetMaximum(ymax*factor)
    
    leg = TLegend(0.53,0.50,0.87,0.87)
    leg.SetFillStyle(4000) # will be transparent
    leg.SetFillColor(0)
    leg.SetTextFont(42)
    leg.SetBorderSize(0)
    for run in runs:
        label = run+": Trig="+dely[run]+", Strb="+strb[run]
        leg.AddEntry(histos[run+"_"+prefix+"_ALPIDE_0"],label,"lp")
    
    cnv = TCanvas("cnv","",1200,1000)
    cnv.Divide(2,2)
    for count1,det in enumerate(detectors):
        p = cnv.cd(count1+1)
        p.SetTicks(1,1)
        # if("cls_size" in prefix): p.SetLogy()
        hname = prefix+"_"+det
        for count2,run in enumerate(runs):                
            if(count2==0): histos[run+"_"+hname].Draw("e1p")
            else:          histos[run+"_"+hname].Draw("e1p same")
            ### fit 
            if(run==run2fit and "h_Chi2fit_res_trk2cls" in prefix):
                func = fit1(histos[run+"_"+hname],runscol[run],-0.01,+0.01)
                s = ROOT.TLatex()
                s.SetNDC(1);
                s.SetTextAlign(13);
                s.SetTextColor(ROOT.kBlack)
                s.SetTextFont(22)
                s.SetTextSize(0.045)
                s.DrawLatex(0.17,0.85,ROOT.Form("Mean: %.2f #mum" % (mm2um*func.GetParameter(1))))
                s.DrawLatex(0.17,0.78,ROOT.Form("Sigma: %.2f #mum" % (mm2um*func.GetParameter(2))))
                s.DrawLatex(0.2,0.71,ROOT.Form("#chi^{2}/N_{DOF}: %.2f" % (func.GetChisquare()/func.GetNDF())))
            
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
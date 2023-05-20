#!/usr/bin/python
import multiprocessing as mp
# from multiprocessing.pool import ThreadPool
import time
import os
import os.path
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
import pickle
from pathlib import Path
import ctypes

import argparse
parser = argparse.ArgumentParser(description='serial_analyzer.py...')
parser.add_argument('-conf', metavar='config file', required=True,  help='full path to config file')
parser.add_argument('-det', metavar='detector to align', required=True,  help='detector to align')
argus = parser.parse_args()
configfile = argus.conf
aligndet = argus.det

import config
from config import *
### must be called here (first) and only once!
init_config(configfile,False)


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

### defined below as global
allhistos = {}


def getfileslist(directory,pattern,suff):
    files = Path( os.path.expanduser(directory) ).glob(pattern+'*'+suff)
    ff = []
    for f in files: ff.append(f)
    return ff


def getfiles(tfilenamein):
    words = tfilenamein.split("/")
    directory = ""
    for w in range(len(words)-1):
        directory += words[w]+"/"
    strippedname = words[-1].split(".pkl")[0]
    words = strippedname.split("_")
    pattern = ""
    for w in range(len(words)):
        word = words[w].replace(".root","")
        pattern += word+"_"
    print("directory:",directory)
    print("pattern:",pattern)
    files = getfileslist(directory,pattern,".pkl")
    return files


def fitSVD(event,aligndet,dx,dy,theta=999):
    clsx = {}
    clsy = {}
    clsz = {}
    clsdx = {}
    clsdy = {}
    for det in cfg["detectors"]:
        x = event.clusters[det][0].xmm
        y = event.clusters[det][0].ymm
        z = event.clusters[det][0].zmm
        if(det==aligndet):
            if(abs(theta)<np.pi):
                x,y = rotate(theta,x,y)
            x = x+dx
            y = y+dy
        clsx.update({det:x})
        clsy.update({det:y})
        clsz.update({det:z})
        clsdx.update({det:event.clusters[det][0].dxmm})
        clsdy.update({det:event.clusters[det][0].dymm})
    vtx  = [cfg["xVtx"], cfg["yVtx"],  cfg["zVtx"]]  if(cfg["doVtx"]) else []
    evtx = [cfg["exVtx"],cfg["eyVtx"], cfg["ezVtx"]] if(cfg["doVtx"]) else []
    points_SVD,errors_SVD = SVD_candidate(clsx,clsy,clsz,clsdx,clsdy,vtx,evtx)
    chisq,ndof,direction,centroid = fit_3d_SVD(points_SVD,errors_SVD)
    chi2ndof_SVD = chisq/ndof if(ndof>0) else 99999
    return chi2ndof_SVD

    
def analyze(fpkl,aligndet,suff):
    lock = mp.Lock()
    lock.acquire()
    
    tfoname = str(fpkl).replace(".pkl","_alignment.root")
    tfo = TFile(tfoname,"RECREATE")
    tfo.cd()
    histos = book_alignment_histos(tfo)
    for name,hist in histos.items():
        hist.SetName(name+suff)
        hist.SetDirectory(0)
    
    with open(fpkl, 'rb') as handle:
        data = pickle.load(handle)
        for event in data:
            ### chi2
            origchi2dof = event.track.chi2ndof
            chi2dof = fitSVD(event,aligndet,0,0,999)
            histos["hSVDchi2dof"].Fill(chi2dof)
            histos["hChi2dof"].Fill(origchi2dof)
            if(chi2dof>cfg["maxChi2align"]): continue
            ### scan x-y-theta misalignment
            for bx in range(1,histos["hTransform"].GetNbinsX()+1):
                for by in range(1,histos["hTransform"].GetNbinsY()+1):
                    for bt in range(1,histos["hTransform"].GetNbinsZ()+1):
                        dx = histos["hTransform"].GetXaxis().GetBinCenter(bx)
                        dy = histos["hTransform"].GetYaxis().GetBinCenter(by)
                        dt = histos["hTransform"].GetZaxis().GetBinCenter(bt)
                        chi2dof = fitSVD(event,aligndet,dx,dy,dt)
                        histos["hTransform"].Fill(dx,dy,dt,chi2dof)
    print("Worker of",fpkl,"is done!")
    lock.release()
    return histos


def collect_errors(error):
    ### https://superfastpython.com/multiprocessing-pool-error-callback-functions-in-python/
    print(f'Error: {error}', flush=True)

def collect_histos(histos):
    ### https://www.machinelearningplus.com/python/parallel-processing-python/
    global allhistos ### defined above!!!
    for name,hist in allhistos.items():
        hist.Add(histos[name])


if __name__ == "__main__":
    # get the start time
    st = time.time()
    
    ### architecture depndent
    nCPUs = mp.cpu_count()
    print("nCPUs:",nCPUs)
    
    # print config once
    show_config()
    if(aligndet not in cfg["detectors"]):
        print("Unknown detector:",aligndet," --> quitting")
        quit()
    
    # Create a pool of workers
    pool = mp.Pool(nCPUs)
    
    tfilenamein = cfg["inputfile"]
    files = getfiles(tfilenamein)
    
    ### histos
    tfoname = tfilenamein.replace(".root","_alignment_"+aligndet+".root")
    tfo = TFile(tfoname,"RECREATE")
    tfo.cd()
    allhistos = book_alignment_histos(tfo)
    
    for fpkl in files:
        suff = str(fpkl).split("_")[-1].replace(".pkl","")
        print("Sending job for",fpkl)
        pool.apply_async(analyze, args=(fpkl,aligndet,suff), callback=collect_histos, error_callback=collect_errors)
        
    ### Wait for all the workers to finish
    pool.close()
    pool.join()

    hXY = allhistos["hTransform"].Project3D("yx")
    hXT = allhistos["hTransform"].Project3D("zx")
    hYT = allhistos["hTransform"].Project3D("zy")
    hX = allhistos["hTransform"].Project3D("x")
    hY = allhistos["hTransform"].Project3D("y")
    hT = allhistos["hTransform"].Project3D("z")
    
    hXT.SetName("hTransform_tx")
    hYT.SetName("hTransform_ty")
    hT.SetName("hTransform_t")

    bx = ctypes.c_int(-1)
    by = ctypes.c_int(-1)
    bt = ctypes.c_int(-1)
    allhistos["hTransform"].GetMinimumBin(bx,by,bt)
    x = allhistos["hTransform"].GetXaxis().GetBinCenter(bx.value)
    y = allhistos["hTransform"].GetYaxis().GetBinCenter(by.value)
    t = allhistos["hTransform"].GetZaxis().GetBinCenter(bt.value)
    print("3D misalignment for "+aligndet+" in x is:",x,"[mm] (or in 1D:",hX.GetBinCenter(hX.GetMinimumBin()),")")
    print("3D misalignment for "+aligndet+" in y is:",y,"[mm] (or in 1D:",hY.GetBinCenter(hY.GetMinimumBin()),")")
    print("3D misalignment for "+aligndet+" in theta is:",t," (or in 1D:",hT.GetBinCenter(hT.GetMinimumBin()),")")
    print("")
    

    tfo.Write()
    tfo.Close()
    
    # get the end time
    et = time.time()
    # get the execution time
    elapsed_time = et - st
    print('Execution time:', elapsed_time, 'seconds')
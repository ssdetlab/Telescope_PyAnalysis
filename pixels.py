#!/usr/bin/python
import os
import math
import array
import numpy as np
import ROOT
from ROOT import *

import config
from config import *
import objects
from objects import *


def get_all_pixles_eudaq(evt,hPixMatrix):
    planes = evt.event.ch_ev_buffer
    pixels = {}
    raws   = {}
    for det in cfg["detectors"]:
        pixels.update({det:[]})
        raws.update({det:[]})
    n_active_planes = 0
    for ipln in range(planes.size()):
        planeid = planes[ipln].plane_id
        detector = cfg["detectorslist"][planeid]
        nhits = planes[ipln].hits.size()
        n_active_planes += (nhits>0)
        for ipix in range(nhits):
            iy,ix = planes[ipln].hits[ipix] ##TODO: check if it is swapped
            raw = hPixMatrix[detector].FindBin(ix,iy)
            if(raw not in raws[detector]):
                raws[detector].append(raw)
                pixels[detector].append( Hit(detector,ix,iy,raw) )
    return n_active_planes,pixels

    
def get_all_pixles_mc(evt,hPixMatrix):
    pixels = {}
    raws = {}
    
    def get_det_pixles_mc(det,col,row):
        for i in range(col.size()):
            ix = col[i]
            iy = row[i]
            raw = hPixMatrix[det].FindBin(ix,iy)
            if(raw not in raws[det]):
                raws[det].append(raw)
                pixels[det].append( Hit(det,ix,iy,raw) )
    
    for det in cfg["detectors"]:
        pixels.update({det:[]})
        raws.update({det:[]})
    ndet = len(cfg["detectors"])
    n_active_planes = (evt.ALPIDE_0_pix_col.size()>0) ## at least one... 
    if(ndet==2): n_active_planes += (evt.ALPIDE_1_pix_col.size()>0)
    if(ndet==3): n_active_planes += (evt.ALPIDE_2_pix_col.size()>0)
    if(ndet==4): n_active_planes += (evt.ALPIDE_3_pix_col.size()>0)
    if(ndet==5): n_active_planes += (evt.ALPIDE_4_pix_col.size()>0)
    if(ndet==6): n_active_planes += (evt.ALPIDE_5_pix_col.size()>0)
    if(ndet==7): n_active_planes += (evt.ALPIDE_6_pix_col.size()>0)
    if(ndet==8): n_active_planes += (evt.ALPIDE_7_pix_col.size()>0)
    if(ndet>0): get_det_pixles_mc("ALPIDE_0",evt.ALPIDE_0_pix_col,evt.ALPIDE_0_pix_row)
    if(ndet>1): get_det_pixles_mc("ALPIDE_1",evt.ALPIDE_1_pix_col,evt.ALPIDE_1_pix_row)
    if(ndet>2): get_det_pixles_mc("ALPIDE_2",evt.ALPIDE_2_pix_col,evt.ALPIDE_2_pix_row)
    if(ndet>3): get_det_pixles_mc("ALPIDE_3",evt.ALPIDE_3_pix_col,evt.ALPIDE_3_pix_row)
    return n_active_planes,pixels



def get_all_pixles_cvr(evt,hPixMatrix):
    pixels = {}
    raws = {}
    
    def get_det_pixles_cvr(det,evtsensor):
        for i in range(evtsensor.size()):
            ix = evtsensor[i].column()
            iy = evtsensor[i].row()
            q  = evtsensor[i].charge()
            raw = hPixMatrix[det].FindBin(ix,iy)
            if(raw not in raws[det]):
                raws[det].append(raw)
                pixels[det].append( Hit(det,ix,iy,raw,q) )
    
    for det in cfg["detectors"]:
        pixels.update({det:[]})
        raws.update({det:[]})
    ndet = len(cfg["detectors"])
    n_active_planes = (evt.ALPIDE_0.size()>0) ## at least one...
    if(ndet==2): n_active_planes += (evt.ALPIDE_1.size()>0)
    if(ndet==3): n_active_planes += (evt.ALPIDE_2.size()>0)
    if(ndet==4): n_active_planes += (evt.ALPIDE_3.size()>0)
    if(ndet==5): n_active_planes += (evt.ALPIDE_4.size()>0)
    if(ndet==6): n_active_planes += (evt.ALPIDE_5.size()>0)
    if(ndet==7): n_active_planes += (evt.ALPIDE_6.size()>0)
    if(ndet==8): n_active_planes += (evt.ALPIDE_7.size()>0)
    if(ndet>0): get_det_pixles_cvr("ALPIDE_0",evt.ALPIDE_0)
    if(ndet>1): get_det_pixles_cvr("ALPIDE_1",evt.ALPIDE_1)
    if(ndet>2): get_det_pixles_cvr("ALPIDE_2",evt.ALPIDE_2)
    if(ndet>3): get_det_pixles_cvr("ALPIDE_3",evt.ALPIDE_3)
    if(ndet>4): get_det_pixles_cvr("ALPIDE_4",evt.ALPIDE_4)
    if(ndet>5): get_det_pixles_cvr("ALPIDE_5",evt.ALPIDE_5)
    if(ndet>6): get_det_pixles_cvr("ALPIDE_6",evt.ALPIDE_6)
    if(ndet>7): get_det_pixles_cvr("ALPIDE_7",evt.ALPIDE_7)
    return n_active_planes,pixels



def get_all_pixles(evt,hPixMatrix,isCVRroot=False):
    n_active_planes = -1
    pixels = {}
    if(not cfg["isMC"]): n_active_planes,pixels = get_all_pixles_eudaq(evt,hPixMatrix)
    else:
        if(isCVRroot):   n_active_planes,pixels = get_all_pixles_cvr(evt,hPixMatrix)
        else:            n_active_planes,pixels = get_all_pixles_mc(evt,hPixMatrix)
    return n_active_planes,pixels

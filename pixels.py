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
    for det in cfg["detectors"]:
        pixels.update({det:[]})
        raws.update({det:[]})
    n_active_planes = (evt.ALPIDE_0_pix_col.size()>0) + (evt.ALPIDE_1_pix_col.size()>0) + (evt.ALPIDE_2_pix_col.size()>0)
    for i in range(evt.ALPIDE_0_pix_col.size()):
        det = "ALPIDE_0"
        ix = evt.ALPIDE_0_pix_col[i]
        iy = evt.ALPIDE_0_pix_row[i]
        raw = hPixMatrix[det].FindBin(ix,iy)
        if(raw not in raws[det]):
            raws[det].append(raw)
            pixels[det].append( Hit(det,ix,iy,raw) )
    for i in range(evt.ALPIDE_1_pix_col.size()):
        det = "ALPIDE_1"
        ix = evt.ALPIDE_0_pix_col[i]
        iy = evt.ALPIDE_0_pix_row[i]
        raw = hPixMatrix[det].FindBin(ix,iy)
        if(raw not in raws[det]):
            raws[det].append(raw)
            pixels[det].append( Hit(det,ix,iy,raw) )
    for i in range(evt.ALPIDE_2_pix_col.size()):
        det = "ALPIDE_2"
        ix = evt.ALPIDE_0_pix_col[i]
        iy = evt.ALPIDE_0_pix_row[i]
        raw = hPixMatrix[det].FindBin(ix,iy)
        if(raw not in raws[det]):
            raws[det].append(raw)
            pixels[det].append( Hit(det,ix,iy,raw) )
    return n_active_planes,pixels


def get_all_pixles_cvr(evt,hPixMatrix):
    pixels = {}
    raws = {}
    for det in cfg["detectors"]:
        pixels.update({det:[]})
        raws.update({det:[]})
    n_active_planes = (evt.ALPIDE_0.size()>0) + (evt.ALPIDE_1.size()>0) + (evt.ALPIDE_2.size()>0)
    for i in range(evt.ALPIDE_0.size()):
        det = "ALPIDE_0"
        ix = evt.ALPIDE_0[i].column()
        iy = evt.ALPIDE_0[i].row()
        q  = evt.ALPIDE_0[i].charge()
        raw = hPixMatrix[det].FindBin(ix,iy)
        if(raw not in raws[det]):
            raws[det].append(raw)
            pixels[det].append( Hit(det,ix,iy,raw,q) )
    for i in range(evt.ALPIDE_1.size()):
        det = "ALPIDE_1"
        ix = evt.ALPIDE_1[i].column()
        iy = evt.ALPIDE_1[i].row()
        q  = evt.ALPIDE_1[i].charge()
        raw = hPixMatrix[det].FindBin(ix,iy)
        if(raw not in raws[det]):
            raws[det].append(raw)
            pixels[det].append( Hit(det,ix,iy,raw,q) )
    for i in range(evt.ALPIDE_2.size()):
        det = "ALPIDE_2"
        ix = evt.ALPIDE_2[i].column()
        iy = evt.ALPIDE_2[i].row()
        q  = evt.ALPIDE_2[i].charge()
        raw = hPixMatrix[det].FindBin(ix,iy)
        if(raw not in raws[det]):
            raws[det].append(raw)
            pixels[det].append( Hit(det,ix,iy,raw,q) )
    return n_active_planes,pixels


def get_all_pixles(evt,hPixMatrix,isCVRroot=False):
    n_active_planes = -1
    pixels = {}
    if(not cfg["isMC"]): n_active_planes,pixels = get_all_pixles_eudaq(evt,hPixMatrix)
    else:
        if(isCVRroot):   n_active_planes,pixels = get_all_pixles_cvr(evt,hPixMatrix)
        else:            n_active_planes,pixels = get_all_pixles_mc(evt,hPixMatrix)
    return n_active_planes,pixels

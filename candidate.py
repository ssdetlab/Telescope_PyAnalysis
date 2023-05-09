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


def SVD_candidate(clsx,clsy,clsz,clsdx,clsdy,vtx=[],evtx=[]):
    isvtx = (len(vtx)>0 and len(evtx)>0)
    clusters = [ [xVtx, yVtx, zVtx] ]  if(isvtx) else []
    clerrors = [ [exVtx,eyVtx,ezVtx] ] if(isvtx) else []
    for det in detectors:
        clusters.append( [clsx[det]+offsets_x[det], clsy[det]+offsets_y[det], clsz[det]] )
        clerrors.append( [clsdx[det],               clsdy[det],               ezCls] )
    points = np.array(clusters)
    errors = np.array(clerrors)
    return points,errors


def Chi2_candidate(clsx,clsy,clsz,clsdx,clsdy,vtx=[],evtx=[]):
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
        clerrors_x.append( clsdx[det] )
        clerrors_y.append( clsdy[det] )
        clerrors_z.append( ezCls )
    points = np.array([ clusters_x,clusters_y,clusters_z ])
    errors = np.array([ clerrors_x,clerrors_y,clerrors_z ])
    return points,errors
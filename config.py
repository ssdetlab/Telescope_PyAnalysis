#!/usr/bin/python
import os
import math
import array
import numpy as np
import sys

isMC    = True
isCVMFS = False

npix_x = 1024
npix_y = 512

pix_x = 0.02924
pix_y = 0.02688

chipX = npix_x*pix_x
chipY = npix_y*pix_y

xVtx = 0
yVtx = 0
zVtx = 0

exVtx = 1.0
eyVtx = 1.0
ezVtx = 0.05

## redundant
ezCls = 0.1 # mm. fixed

lineScaleUp = 70
lineScaleDn = 50

nprintout = 20

### for noize masking
pTrim    = 0.01
zeroSupp = True
nSigma   = 25

### detectors
detectors = ["ALPIDE_0", "ALPIDE_1", "ALPIDE_2"]
detectorslist = list(detectors)
rdetectorslist = detectorslist.reverse()
rdetectors = {"ALPIDE_0":[0,0,29.0], "ALPIDE_1":[0,0,54.8], "ALPIDE_2":[0,0,80.6]}
### for the fit line edges
zFirst = zVtx*0.9
zLast  = rdetectors["ALPIDE_2"][2]*1.1
### world dimensions for plots
world = {"x":[-1.2*chipX, +1.2*chipX], "y":[-1.5*chipY, +1.5*chipY], "z":[zFirst,zLast]}

### for the sphere
sphere_center_point = [ rdetectors["ALPIDE_1"][0], rdetectors["ALPIDE_1"][1], rdetectors["ALPIDE_1"][2] ]
sphere_radius_size  = (rdetectors["ALPIDE_2"][2]-rdetectors["ALPIDE_0"][2])*0.7

### offsets
offsets_x = {}
offsets_y = {}
for det in detectors:
    offsets_x.update( {det:rdetectors[det][0]} )
    offsets_y.update( {det:rdetectors[det][1]} )

def reconfig(mc,cvmfs):
    isMC = mc
    isCVMFS = cvmfs
    # print("isMC=",isMC)
    # print("isCVMFS=",isCVMFS)
    ### reset path
    if(isCVMFS and isMC): sys.path.insert(0, '/storage/agrp/rouxs/analysis_sr90')
    ### correct offset in y for simulation
    if(isMC):
        for det in detectors:
            rdetectors[det][1] = 0.61872
            ### correct the offset in y
            offsets_y[det] = rdetectors[det][1]
        ### correct the sphere in y
        sphere_center_point[1] = rdetectors["ALPIDE_1"][1]

### cuts
cuts = ["All", "N_{hits/det}>0", "N_{cls/det}==1", "Chi2 Fitted", "Fit #chi^{2}/N_{DoF}#leq450"]
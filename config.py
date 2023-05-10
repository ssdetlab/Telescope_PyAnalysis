#!/usr/bin/python
import os
import math
import array
import numpy as np
import sys
import configparser

##########################################################
##########################################################
##########################################################

### should be called once from main
def init_config(fname,show):
    ConfigCls = Config(fname,show)
    global cfg
    cfg = ConfigCls.map

def show_config():
    print("Configuration map:",cfg)
    print("")

### config file looks like that:
# [SECTION_NAME]
# key1 = value1
# key2 = value2
class Config:
    def __init__(self,fname,doprint=False):
        self.fname = fname
        self.doprint = doprint
        self.configurator = configparser.RawConfigParser()
        self.configurator.optionxform = str ### preserve case sensitivity
        self.map = {} ### the config map
        self.set(fname,doprint)
        
    def read(self,fname):
        if(self.doprint): print("Reading configuration from: ",fname)
        self.configurator.read(fname)
        
    def getF(self,section,var):
        expr = dict(self.configurator.items(section))[var]
        if(not expr.isnumeric()):
            return float(eval(expr))
        return float(dict(self.configurator.items(section))[var])
       
    def getI(self,section,var):
        expr = dict(self.configurator.items(section))[var]
        if(not expr.isnumeric()):
            return int(eval(expr))
        return int(dict(self.configurator.items(section))[var])
    
    def getB(self,section,var):
        return True if(int(dict(self.configurator.items(section))[var])==1) else False
    
    def getS(self,section,var):
        return str(dict(self.configurator.items(section))[var])

    def getArrS(self,section,var):
        s = self.getS(section,var)
        return s.split(" ")
    
    def getArrI(self,section,var):
        s = self.getS(section,var).split(",")
        i = [int(x) for x in s]
        return i
        
    def getArrF(self,section,var):
        s = self.getS(section,var).split(",")
        f = [float(x) for x in s]
        return f
    
    def getMap2ArrF(self,section,var):
        s = self.getS(section,var).split(" ")
        m = {}
        for x in s:
            x = x.split(":")
            name = x[0]
            sarr = x[1].split(",")
            farr = [float(n) for n in sarr]
            m.update({name:farr})
        return m
    
    def add(self,name,var):
        self.map.update( {name:var} )
    
    def set(self,fname,doprint=False):
        ### read
        self.read(fname)
        ### set
        self.add("isMC", self.getB('RUN','isMC'))
        self.add("isCVMFS", self.getB('RUN','isCVMFS'))
        self.add("doVtx", self.getB('RUN','doVtx'))
        self.add("runtype", self.getS('RUN','runtype'))
        self.add("pdgIdMatch", self.getI('RUN','pdgIdMatch'))
        self.add("nmax2process", self.getI('RUN','nmax2process'))
        self.add("doplot", self.getB('RUN','doplot'))
        self.add("doDiagnostics", self.getB('RUN','doDiagnostics'))
        self.add("doNoiseScan", self.getB('RUN','doNoiseScan'))
        self.add("isCVRroot", self.getB('RUN','isCVRroot'))
        self.add("nprintout", self.getI('RUN','nprintout'))
        self.add("inputfile", self.getS('RUN','inputfile'))

        self.add("npix_x", self.getI('CHIP','npix_x'))
        self.add("npix_y", self.getI('CHIP','npix_y'))
        self.add("pix_x",  self.getF('CHIP','pix_x'))
        self.add("pix_y",  self.getF('CHIP','pix_y'))
        self.add("chipX",  self.map["npix_x"]*self.map["pix_x"])
        self.add("chipY",  self.map["npix_y"]*self.map["pix_y"])

        self.add("xVtx", self.getB('VTX','xVtx'))
        self.add("yVtx", self.getB('VTX','yVtx'))
        self.add("zVtx", self.getB('VTX','zVtx'))
        self.add("exVtx", self.getF('VTX','exVtx'))
        self.add("eyVtx", self.getF('VTX','exVtx'))
        self.add("ezVtx", self.getF('VTX','exVtx'))

        self.add("ezCls", self.getF('CLUSTER','ezCls'))

        self.add("lineScaleUp", self.getF('WORLD','lineScaleUp'))
        self.add("lineScaleDn", self.getF('WORLD','lineScaleDn'))

        self.add("pTrim", self.getF('NOISE','pTrim'))
        self.add("zeroSupp", self.getB('NOISE','zeroSupp'))
        self.add("nSigma", self.getF('NOISE','nSigma'))

        self.add("detectors", self.getArrS('DETECTOR','detectors'))
        self.add("rdetectors", self.getMap2ArrF('DETECTOR','rdetectors'))
        self.add("worldmargins", self.getF('DETECTOR','worldmargins'))
        self.add("zFirst", self.map["zVtx"]*(1-self.map["worldmargins"]))
        self.add("zLast",  self.map["rdetectors"]["ALPIDE_2"][2]*(1+self.map["worldmargins"]))
        self.add("worldscales", self.getMap2ArrF('DETECTOR','worldscales'))
        self.add("worldcenter", self.getArrF('DETECTOR','worldcenter'))
        self.add("worldradius",  self.getF('DETECTOR','worldradius'))
        world = {}
        for axis,scales in self.map["worldscales"].items():
            bounds = -9999
            if(axis=="x"): bounds = [ -self.map["chipX"]*scales[0],+self.map["chipX"]*scales[1] ]
            if(axis=="y"): bounds = [ -self.map["chipY"]*scales[0],+self.map["chipY"]*scales[1] ]
            if(axis=="z"): bounds = [ self.map["zFirst"]*scales[0], self.map["zLast"]*scales[1] ]
            world.update( {axis:bounds} )
        self.add("world", world)
        
        offsets_x = {}
        offsets_y = {}
        for det in self.map["detectors"]:
            offsets_x.update( {det:self.map["rdetectors"][det][0]} )
            offsets_y.update( {det:self.map["rdetectors"][det][1]} )
        self.add("offsets_x", offsets_x)
        self.add("offsets_y", offsets_y)
        
        self.add("cuts", self.getArrS('CUTS','cuts'))
    
        if(doprint):
            print("Configuration map:",self.map)
            print("")

    def __str__(self):
        return f"Config map: {self.map}"



##########################################################
##########################################################
##########################################################


# isMC    = True
# isCVMFS = False
#
# npix_x = 1024
# npix_y = 512
#
# pix_x = 0.02924
# pix_y = 0.02688
#
# chipX = npix_x*pix_x
# chipY = npix_y*pix_y
#
# xVtx = 0
# yVtx = 0
# zVtx = 0
#
# exVtx = 1.0
# eyVtx = 1.0
# ezVtx = 0.05
#
# ## redundant
# ezCls = 0.1 # mm. fixed
#
# lineScaleUp = 70
# lineScaleDn = 50
#
# nprintout = 20
#
# ### for noize masking
# pTrim    = 0.01
# zeroSupp = True
# nSigma   = 25

### detectors
# detectors = ["ALPIDE_0", "ALPIDE_1", "ALPIDE_2"]
### detectorslist = list(detectors)
### rdetectorslist = detectorslist.reverse()
# rdetectors = {"ALPIDE_0":[0,0,29.0], "ALPIDE_1":[0,0,54.8], "ALPIDE_2":[0,0,80.6]}
### for the fit line edges
# zFirst = zVtx*0.9
# zLast  = rdetectors["ALPIDE_2"][2]*1.1
### world dimensions for plots
# world = {"x":[-1.2*chipX, +1.2*chipX], "y":[-1.5*chipY, +1.5*chipY], "z":[zFirst,zLast]}

### for the sphere
# sphere_center_point = [ rdetectors["ALPIDE_1"][0], rdetectors["ALPIDE_1"][1], rdetectors["ALPIDE_1"][2] ]
# sphere_radius_size  = (rdetectors["ALPIDE_2"][2]-rdetectors["ALPIDE_0"][2])*0.7

# ### offsets
# offsets_x = {}
# offsets_y = {}
# for det in detectors:
#     offsets_x.update( {det:rdetectors[det][0]} )
#     offsets_y.update( {det:rdetectors[det][1]} )

# def reconfig(mc,cvmfs):
#     isMC = mc
#     isCVMFS = cvmfs
#     # print("isMC=",isMC)
#     # print("isCVMFS=",isCVMFS)
#     ### reset path
#     if(isCVMFS and isMC): sys.path.insert(0, '/storage/agrp/rouxs/analysis_sr90')
#     ### correct offset in y for simulation
#     if(isMC):
#         for det in detectors:
#             rdetectors[det][1] = 0.61872
#             ### correct the offset in y
#             offsets_y[det] = rdetectors[det][1]
#         ### correct the sphere in y
#         sphere_center_point[1] = rdetectors["ALPIDE_1"][1]

### cuts
# cuts = ["All", "N_{hits/det}>0", "N_{cls/det}==1", "Chi2 Fitted", "Fit #chi^{2}/N_{DoF}#leq450"]
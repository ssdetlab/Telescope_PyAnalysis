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
parser = argparse.ArgumentParser(description='corryvreckan_ascii2root.py...')
parser.add_argument('-f', metavar='input file', required=True,  help='full path to input file')
parser.add_argument('-a', metavar='write all?', required=False, help='write all? [0/1]')
parser.add_argument('-l', metavar='nmax lines to read', required=False,  help='nmax lines to read')
parser.add_argument('-n', metavar='nmax events to read', required=False,  help='nmax events to read')
argus = parser.parse_args()
txtfilename = argus.f
lmaxtoread  = -1 if(argus.l is None) else int(argus.l)
nmaxtoread  = -1 if(argus.n is None) else int(argus.n)
writeAll    =  0 if(argus.a is None) else int(argus.a)

print("--------------------------")
print("Settings:")
print("  nmaxtoread:",nmaxtoread)
print("  writeAll:",writeAll)
print("--------------------------")

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(0);
ROOT.gStyle.SetOptStat(0);

### histograms
histos = {}

### brnches
branches = {}

#############################################################################
#############################################################################
#############################################################################

def init_event():
    event = {}
    event.update( {"evt":-1} )
    for det in detectors:
        event.update( {det:{"mc_pdgId":[], "mc_loc_start":[], "mc_loc_end":[],
                            "pix_col":[],"pix_row":[],"pix_raw":[],"pix_crg":[],"pix_time":[],
                            "cls_loc_x":[], "cls_loc_y":[], "cls_glo_r":[], "cls_crg":[], "cls_split":[], "cls_npix":[], "cls_ncol":[], "cls_nrow":[], "cls_time":[],
                           }
                      }
                    )
    event.update( {"trk_state":[], "trk_dir":[], "trk_chi2":[], "trk_ndof":[], "trk_chi2ndof":[], "trk_time":[]} )
    return event



def set_branches():
    branches.update({ "evt":ROOT.std.vector( int )() })
    branches.update({ "trk_state_x":ROOT.std.vector( float )() })
    branches.update({ "trk_state_y":ROOT.std.vector( float )() })
    branches.update({ "trk_state_z":ROOT.std.vector( float )() })
    branches.update({ "trk_dir_x":ROOT.std.vector( float )() })
    branches.update({ "trk_dir_y":ROOT.std.vector( float )() })
    branches.update({ "trk_dir_z":ROOT.std.vector( float )() })
    branches.update({ "trk_chi2":ROOT.std.vector( float )() })
    branches.update({ "trk_ndof":ROOT.std.vector( int )() })
    branches.update({ "trk_chi2ndof":ROOT.std.vector( float )() })
    branches.update({ "trk_time":ROOT.std.vector( float )() })
    for det in detectors:
        branches.update({ det+"_mc_pdgId":ROOT.std.vector( int )() })
        branches.update({ det+"_mc_loc_start_x":ROOT.std.vector( float )() })
        branches.update({ det+"_mc_loc_start_y":ROOT.std.vector( float )() })
        branches.update({ det+"_mc_loc_start_z":ROOT.std.vector( float )() })
        branches.update({ det+"_mc_loc_end_x":ROOT.std.vector( float )() })
        branches.update({ det+"_mc_loc_end_y":ROOT.std.vector( float )() })
        branches.update({ det+"_mc_loc_end_z":ROOT.std.vector( float )() })
        
        branches.update({ det+"_npix":-1 })
        branches.update({ det+"_ncls":-1 })
        
        branches.update({ det+"_pix_col":ROOT.std.vector( int )() })
        branches.update({ det+"_pix_row":ROOT.std.vector( int )() })
        branches.update({ det+"_pix_raw":ROOT.std.vector( int )() })
        branches.update({ det+"_pix_crg":ROOT.std.vector( float )() })
        branches.update({ det+"_pix_time":ROOT.std.vector( float )() })
        
        branches.update({ det+"_cls_loc_x":ROOT.std.vector( float )() })
        branches.update({ det+"_cls_loc_y":ROOT.std.vector( float )() })
        branches.update({ det+"_cls_glo_x":ROOT.std.vector( float )() })
        branches.update({ det+"_cls_glo_y":ROOT.std.vector( float )() })
        branches.update({ det+"_cls_glo_z":ROOT.std.vector( float )() })
        branches.update({ det+"_cls_crg":ROOT.std.vector( float )() })
        branches.update({ det+"_cls_split":ROOT.std.vector( int )() })
        branches.update({ det+"_cls_npix":ROOT.std.vector( int )() })
        branches.update({ det+"_cls_ncol":ROOT.std.vector( int )() })
        branches.update({ det+"_cls_nrow":ROOT.std.vector( int )() })
        branches.update({ det+"_cls_time":ROOT.std.vector( float )() })

def clean_branches():
    for branch in branches:
        if("_npix" in branch and "_cls_npix" not in branch): continue
        if("_ncls" in branch): continue
        branches[branch].clear()
        
def init_branches(tree):
    for branch in branches:
        tree.Branch(branch, branches[branch])

def fill_branches(event):
    ### global parameters
    branches["evt"].push_back( event["evt"] )
    ### global tracks
    for t in range(len(event["trk_state"])):
        branches["trk_state_x"].push_back( event["trk_state"][t].X() )
        branches["trk_state_y"].push_back( event["trk_state"][t].Y() )
        branches["trk_state_z"].push_back( event["trk_state"][t].Z() )
        branches["trk_dir_x"].push_back( event["trk_dir"][t].X() )
        branches["trk_dir_y"].push_back( event["trk_dir"][t].Y() )
        branches["trk_dir_z"].push_back( event["trk_dir"][t].Z() )
        branches["trk_chi2"].push_back( event["trk_chi2"][t] )
        branches["trk_ndof"].push_back( event["trk_ndof"][t] )
        branches["trk_chi2ndof"].push_back( event["trk_chi2ndof"][t] )
        branches["trk_time"].push_back( event["trk_time"][t] )
    ### detector specific
    for det in detectors:
        ### global
        branches[det+"_npix"] = len(event[det]["pix_col"])
        branches[det+"_ncls"] = len(event[det]["cls_loc_x"])
        ### mc particles
        for m in range(len(event[det]["mc_pdgId"])):
            branches[det+"_mc_pdgId"].push_back( event[det]["mc_pdgId"][m] )
            branches[det+"_mc_loc_start_x"].push_back( event[det]["mc_loc_start"][m].X() )
            branches[det+"_mc_loc_start_y"].push_back( event[det]["mc_loc_start"][m].Y() )
            branches[det+"_mc_loc_start_z"].push_back( event[det]["mc_loc_start"][m].Z() )            
            branches[det+"_mc_loc_end_x"].push_back( event[det]["mc_loc_end"][m].X() )
            branches[det+"_mc_loc_end_y"].push_back( event[det]["mc_loc_end"][m].Y() )
            branches[det+"_mc_loc_end_z"].push_back( event[det]["mc_loc_end"][m].Z() )
        ### pixels
        for p in range(len(event[det]["pix_col"])):
            name = "pix_col";  branches[det+"_"+name].push_back( event[det][name][p] )
            name = "pix_row";  branches[det+"_"+name].push_back( event[det][name][p] )
            name = "pix_raw";  branches[det+"_"+name].push_back( event[det][name][p] )
            name = "pix_crg";  branches[det+"_"+name].push_back( event[det][name][p] )
            name = "pix_time"; branches[det+"_"+name].push_back( event[det][name][p] )
        ### clusters
        for c in range(len(event[det]["cls_loc_x"])):
            name = "cls_loc_x"; branches[det+"_"+name].push_back( event[det][name][c] )
            name = "cls_loc_y"; branches[det+"_"+name].push_back( event[det][name][c] )
            name = "cls_crg";   branches[det+"_"+name].push_back( event[det][name][c] )
            name = "cls_split"; branches[det+"_"+name].push_back( event[det][name][c] )
            name = "cls_npix";  branches[det+"_"+name].push_back( event[det][name][c] )
            name = "cls_ncol";  branches[det+"_"+name].push_back( event[det][name][c] )
            name = "cls_nrow";  branches[det+"_"+name].push_back( event[det][name][c] )
            name = "cls_time";  branches[det+"_"+name].push_back( event[det][name][c] )
            branches[det+"_cls_glo_x"].push_back( event[det]["cls_glo_r"][c].X() )
            branches[det+"_cls_glo_y"].push_back( event[det]["cls_glo_r"][c].Y() )
            branches[det+"_cls_glo_z"].push_back( event[det]["cls_glo_r"][c].Z() )


### book histos
def book_histos():
    histos.update( { "h_events" : TH1D("h_events",";;Events",1,0,1) } )
    histos["h_events"].GetXaxis().SetBinLabel(1,"")
    
    histos.update( { "h_chi2dof" : TH1D("h_chi2dof",";#chi^{2}/N_{dof};Tracks",100,0,1) } )
    histos.update( { "h_npix"    : TH1D("h_npix",";N_{pixels} per detector;Events",20,0,20) } )
    histos.update( { "h_ncls"    : TH1D("h_ncls",";N_{clusters} per detector;Events",50,0,50) } )
    for det in detectors:
        histos.update( { "h_pixocc1D_"+det : TH1D("h_pixocc1D_"+det,";Pixel;Hits",npix_x*npix_y,1,npix_x*npix_y+1) } )
        histos.update( { "h_pixocc2D_"+det : TH2D("h_pixocc2D_"+det,";x;y;Hits",npix_x+1,-0.5,npix_x+0.5, npix_y+1,-0.5,npix_y+0.5) } )
    
    for hname,hist in histos.items():
        hist.SetLineColor(ROOT.kBlack)
        hist.Sumw2()




def get_words(line):
    line = line.replace("\t"," ")
    words = line.split(" ")
    for w in range(len(words)):
        words[w] = words[w].replace(","," ").replace("\n","").replace("(","").replace(")","")
    return words


def fill_mc(event,detector,words):
    # print(words[1],words[2],words[3],words[4],words[5])
    event[detector]["mc_pdgId"].append( int(words[1]) )
    subwords = words[2].split(" ")
    event[detector]["mc_loc_start"].append( TVector3( float(subwords[0]), float(subwords[1]), float(subwords[2]) ) )
    subwords = words[3].split(" ")
    event[detector]["mc_loc_end"].append( TVector3( float(subwords[0]), float(subwords[1]), float(subwords[2]) ) )

def fill_pixel(event,detector,words):
    # print(words[1],words[2],words[3],words[4],words[5])
    event[detector]["pix_col"].append( int(words[1]) )
    event[detector]["pix_row"].append( int(words[2]) )
    event[detector]["pix_raw"].append( int(words[3]) )
    event[detector]["pix_crg"].append( float(words[4]) )
    event[detector]["pix_time"].append( float(words[5]) )

def fill_cluster(event,detector,words):
    # print(words[1],words[2],words[3],words[4],words[5],words[6],words[7],words[8],words[9])
    event[detector]["cls_loc_x"].append( float(words[1]) )
    event[detector]["cls_loc_y"].append( float(words[2]) )
    subwords = words[3].split(" ")
    event[detector]["cls_glo_r"].append( TVector3( float(subwords[0]), float(subwords[1]), float(subwords[2]) ) )
    event[detector]["cls_crg"].append( float(words[4]) )
    event[detector]["cls_split"].append( int(words[5]) )
    event[detector]["cls_npix"].append( int(words[6]) )
    event[detector]["cls_ncol"].append( int(words[7]) )
    event[detector]["cls_nrow"].append( int(words[8]) )
    event[detector]["cls_time"].append( float(words[9]) )

def fill_track(event,words):
    # print(words[1],words[2],words[3],words[4],words[5],words[6])
    subwords = words[1].split(" ")
    event["trk_state"].append( TVector3( float(subwords[0]), float(subwords[1]), float(subwords[2]) ) )
    subwords = words[2].split(" ")
    event["trk_dir"].append( TVector3( float(subwords[0]), float(subwords[1]), float(subwords[2]) ) )
    event["trk_chi2"].append( float(words[3]) )
    event["trk_ndof"].append( int(words[4]) )
    event["trk_chi2ndof"].append( float(words[5]) )
    event["trk_time"].append( float(words[6]) )



def Parser(lines):
    nLinesInEvent = 0
    detector = ""
    event = init_event()
    for line in lines:
        ### count lines anyway
        nLinesInEvent += 1
        ### lines to ignore
        if(line==""):                continue
        if("#"             in line): continue
        if("Start: "       in line): continue
        if("End: "         in line): continue
        if("Trigger list:" in line): continue
        ### parse the line into words
        words = get_words(line)
        ### start filling data
        if(words[0]=="==="):
            evt = int(words[1])
            if(evt!=event["evt"] and event["evt"]!=-1): break ### go to next event
            event["evt"] = evt
        if(words[0]=="---"):
            if(words[1]!="<global>"): detector = words[1]
            else:                     continue
        if(words[0]=="MCParticle"):
            fill_mc(event,detector,words)
        if(words[0]=="Pixel"):
            fill_pixel(event,detector,words)
        if(words[0]=="Cluster"):
            fill_cluster(event,detector,words)
        if(words[0]=="StraightLineTrack"):
            fill_track(event,words)
    ### end of function
    return event,nLinesInEvent


### get data
def readASCII(txtfilename,rootfilename):
    tf = TFile(rootfilename,"RECREATE")
    tt = TTree("tt","tt")
    set_branches()
    init_branches(tt)
    book_histos()
    
    with open(txtfilename) as fp:
        # lines = fp.readlines()
        lines = [next(fp) for _ in range(lmaxtoread)] if(lmaxtoread>0) else fp.readlines()
        nlines = len(lines)
        nevents = 0
        nselect = 0
        print("Going to read",nlines,"lines")
        while(len(lines)>1):
            
            ### get the event and update the text file
            event,nlines = Parser(lines)
            del lines[:nlines-1]
            
            ### some on-the-fly trivial analysis
            histos["h_events"].Fill(0.5) # for bookeeping
            ncls = 0
            npix = 0
            for det in detectors:
                ncls += (len(event[det]["cls_loc_x"])>0)
                npix += (len(event[det]["pix_col"])>0)
                for p in range(len(event[det]["pix_col"])):
                    ixpix = event[det]["pix_col"][p]
                    iypix = event[det]["pix_row"][p]
                    ipix  = histos["h_pixocc2D_"+det].FindBin(ixpix,iypix)
                    histos["h_pixocc1D_"+det].AddBinContent(ipix,1)
                    histos["h_pixocc2D_"+det].Fill(ixpix,iypix)
            histos["h_npix"].Fill(npix)
            histos["h_ncls"].Fill(ncls)
            for x in event["trk_chi2ndof"]: histos["h_chi2dof"].Fill(x)

            ### selection
            at_least_3_clusters = (len(event[detectors[0]]["cls_loc_x"])>0 and len(event[detectors[1]]["cls_loc_x"])>0 and len(event[detectors[2]]["cls_loc_x"])>0)
            
            ### output tree
            if(at_least_3_clusters or writeAll):
                clean_branches() ### important!!!
                fill_branches(event)
                tt.Fill()
                nselect += 1

            ### monitoring
            nevents += 1
            nlines = len(lines)
            if(nevents%(10*nprintout)==0):
                print("Processed",nevents,"events, Selected",nselect,"events, Lines remaining:",nlines)
            
            ### stopping criteria
            if(nmaxtoread>0 and nselect>=nmaxtoread): break
            
    tf.Write()
    tf.Close()


#############################################################################
#############################################################################
#############################################################################

### run
rootfilename = txtfilename.replace(".txt","_flattree.root")
readASCII(txtfilename,rootfilename)



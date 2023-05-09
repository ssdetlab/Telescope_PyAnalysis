#!/usr/bin/python
import os
import math
import subprocess
import array
import numpy as np
import ROOT
from ROOT import *

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetOptStat(0)

print("-----------------------------------------------------------------------------------")
print("Need to add CVR libs:")
print("export LD_LIBRARY_PATH=$HOME/corryvreckan/corryvreckan-master/lib:$LD_LIBRARY_PATH")
print("-----------------------------------------------------------------------------------")

### histograms
histos = {}

### brnches
branches = {}

### detectors
detectors = ["Stave00", "Stave10", "Stave20"]

print("---- start loading libs")
### see https://root.cern/manual/python/
gInterpreter.AddIncludePath('~/corryvreckan/corryvreckan-master/src/objects/')
gSystem.Load('libCorryvreckanObjects.dylib')
print("---- finish loading libs")

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
    event.update( {"trk_state":[], "trk_intrcpt":[], "trk_dir":[], "trk_chi2":[], "trk_ndof":[], "trk_chi2ndof":[], "trk_time":[]} )
    return event



def init_branches():
    branches.update({ "event":-1 })
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

def set_branches(tree):
    for branch in branches:
        tree.Branch(branch, branches[branch])

def clean_branches():
    for branch in branches:
        if(branch=="event"): continue
        if("_npix" in branch and "_cls_npix" not in branch): continue
        if("_ncls" in branch): continue
        branches[branch].clear()

def fill_branches(event):
    ### global parameters
    branches["event"] = event["evt"]
    ### global tracks
    for t in range(len(event["trk_state"])):
        branches["trk_state_x"].push_back( event["trk_state"][t].X() )
        branches["trk_state_y"].push_back( event["trk_state"][t].Y() )
        branches["trk_state_z"].push_back( event["trk_state"][t].Z() )
        branches["trk_intrcpt_x"].push_back( event["trk_intrcpt"][t].X() )
        branches["trk_intrcpt_y"].push_back( event["trk_intrcpt"][t].Y() )
        branches["trk_intrcpt_z"].push_back( event["trk_intrcpt"][t].Z() )
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
    
    for hname,hist in histos.items():
        hist.SetLineColor(ROOT.kBlack)
        hist.Sumw2()


# def fill_mc(event,detector,words):
#     # print(words[1],words[2],words[3],words[4],words[5])
#     event[detector]["mc_pdgId"].append( int(words[1]) )
#     subwords = words[2].split(" ")
#     event[detector]["mc_loc_start"].append( TVector3( float(subwords[0]), float(subwords[1]), float(subwords[2]) ) )
#     subwords = words[3].split(" ")
#     event[detector]["mc_loc_end"].append( TVector3( float(subwords[0]), float(subwords[1]), float(subwords[2]) ) )



npix_x = 1024
npix_y = 512

pix_x = 0.02924
pix_y = 0.02688

chipX = npix_x*pix_x
chipY = npix_y*pix_y

def get_event(ievt,tCls,tPix,tTrk,tMCp):
    
    event = init_event()
    
    tCls.GetEntry(ievt)
    tPix.GetEntry(ievt)
    tTrk.GetEntry(ievt)
    tMCp.GetEntry(ievt)
    
    for det in detectors:
        
        ### get clusters
        clusters = None
        if  (det=="Stave00"): clusters = tCls.Stave00
        elif(det=="Stave10"): clusters = tCls.Stave10
        elif(det=="Stave20"): clusters = tCls.Stave20
        else:
            print("unknown detector:",det)
            quit()
        for i in range(clusters.size()):
            
            errorX = clusters[i].errorX()
            errorX = clusters[i].errorY()
            
            # event[det]["cls_glo_r"].append( clusters[i].global() ) ### ---> cannot do since "global" is a python3 keyword
            event[det]["cls_glo_r"].append( TVector3(clusters[i].local().X()-chipX/2.+pix_x/2., clusters[i].local().Y()-chipY/2.+pix_y/2., clusters[i].local().Z()) )
            event[det]["cls_loc_x"].append( clusters[i].local().X() )
            event[det]["cls_loc_y"].append( clusters[i].local().Y() )
            event[det]["cls_crg"].append( clusters[i].charge() )
            event[det]["cls_split"].append( clusters[i].isSplit() )
            event[det]["cls_npix"].append( clusters[i].size() )
            event[det]["cls_ncol"].append( clusters[i].columnWidth() )
            event[det]["cls_nrow"].append( clusters[i].rowWidth() )
            event[det]["cls_time"].append( clusters[i].timestamp() )

        ### get pixels
        pixels = None
        if  (det=="Stave00"): pixels = tPix.Stave00
        elif(det=="Stave10"): pixels = tPix.Stave10
        elif(det=="Stave20"): pixels = tPix.Stave20
        else:
            print("unknown detector:",det)
            quit()
        for i in range(pixels.size()):
            event[det]["pix_col"].append( pixels[i].column() )
            event[det]["pix_row"].append( pixels[i].row() )
            event[det]["pix_raw"].append( pixels[i].raw() )
            event[det]["pix_crg"].append( pixels[i].charge() )
            event[det]["pix_time"].append( pixels[i].timestamp() )
        
        ### get MCParticles
        mcprts = None
        if  (det=="Stave00"): mcprts = tMCp.Stave00
        elif(det=="Stave10"): mcprts = tMCp.Stave10
        elif(det=="Stave20"): mcprts = tMCp.Stave20
        else:
            print("unknown detector:",det)
            quit()
        for i in range(mcprts.size()):
            event[det]["mc_pdgId"].append( mcprts[i].getID() )
            event[det]["mc_loc_start"].append( mcprts[i].getLocalStart() )
            event[det]["mc_loc_end"].append( mcprts[i].getLocalEnd() )

        # ### get tracks
        # tracks = tTrk.global
        # for i in range(tracks.size()):
        #     event[det]["trk_state"].append( tracks[i].getState() )
        #     event[det]["trk_intrcpt"].append( tracks[i].getIntercept() )
        #     event[det]["trk_dir"].append( tracks[i].getDirection() )
        #     event[det]["trk_chi2"].append( tracks[i].getChi2() )
        #     event[det]["trk_ndof"].append( tracks[i].getNdof() )
        #     event[det]["pix_time"].append( tracks[i].timestamp() )


    return event


### get data
def readROOT(inrootfilename,outrootfilename):
    tfo = TFile(outrootfilename,"RECREATE")
    tt = TTree("tt","tt")
    init_branches()
    set_branches(tt)
    book_histos()
    
    tfi  = TFile(inrootfilename,"READ")
    tCls = tfi.Get("Cluster")
    tPix = tfi.Get("Pixel")
    tTrk = tfi.Get("Track")
    tMCp = tfi.Get("MCParticle")
    
    nevents = tCls.GetEntries()
    print("Going to process",nevents,"events")
    n3 = 0
    for ievt in range(nevents):
        clean_branches() ### important!!!
        event = get_event(ievt,tCls,tPix,tTrk,tMCp)
        fill_branches(event)
        tt.Fill()
    
        ### for bookeeping
        histos["h_events"].Fill(0.5)
        
        ### some on-the-fly trivial analysis
        ncls = (len(event["Stave00"]["cls_glo_r"])>0)+(len(event["Stave10"]["cls_glo_r"])>0)+(len(event["Stave20"]["cls_glo_r"])>0)
        histos["h_ncls"].Fill(ncls)
        if(ncls>=3): n3+=1  
        
    print("n3",n3)
    tfo.Write()
    tfo.Close()


#############################################################################
#############################################################################
#############################################################################

### run
inrootfilename = "data/corry_out_100_new_unsh.root"
outrootfilename = inrootfilename.replace(".root","_flat.root")
readROOT(inrootfilename,outrootfilename)



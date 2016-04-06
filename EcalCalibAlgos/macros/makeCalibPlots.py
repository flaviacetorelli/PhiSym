# import ROOT in batch mode
from __future__ import division
import sys
import re
import argparse
oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
import ROOT
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv

from PhiSym.EcalCalibAlgos.EcalCalibAnalysis import *
from pluginCondDBPyInterface import *
from pluginEcalPyUtils import *

# Read command line options
parser = argparse.ArgumentParser (description = 'Draw PhiSym intercalibration plots: channels maps --- Eta,Phi / X,Y profiles')
parser.add_argument('inputfile' , default="phisym_intercalibs.root", help='analyze this file')
parser.add_argument('-b', '--block' , default=-1, type=int, help='analyze this block')
parser.add_argument('--plots' , action="store_true", default=False, help='produce .png and .pdf plots')
parser.add_argument('-d', '--outdir' , default="./", help='store path for .png, .pdf and .root')

opts = parser.parse_args ()

# Load the CalibrationFile format
ROOT.gSystem.Load("libPhiSymEcalCalibDataFormats.so");
ROOT.AutoLibraryLoader.enable()

# Setup the style
ROOT.gErrorIgnoreLevel=ROOT.kInfo + 1
ROOT.gStyle.SetOptStat("e")
ROOT.gStyle.SetOptFit(111)
ROOT.gStyle.SetLabelSize(0.05, "XY");
ROOT.gStyle.SetTickLength(0.03, "XYZ");
ROOT.gStyle.SetNdivisions(510, "XYZ");
ROOT.gStyle.SetTitleColor(1, "XYZ");
ROOT.gStyle.SetTitleSize(0.07, "XYZ");
ROOT.gStyle.SetTitleXOffset(0.8);
ROOT.gStyle.SetTitleYOffset(0.7);
ROOT.gStyle.SetPadTopMargin(0.13);
ROOT.gStyle.SetPadBottomMargin(0.13);
ROOT.gStyle.SetPadLeftMargin(0.1);
ROOT.gStyle.SetPadRightMargin(0.15);

# initialize variables
nLumis=-1

# Get the crystals trees
inFile = ROOT.TFile(opts.inputfile)
bareTree = inFile.Get("eb_xstals")
ebTree = ROOT.CrystalsEBTree(bareTree)
bareTree = inFile.Get("ee_xstals")
eeTree = ROOT.CrystalsEETree(bareTree)

## Book the histos ##
histos={}
profiles={}
maps={}
## 1D histo
# EB
histos["EB_ic_ch"]=ROOT.TH1F("histo_EB_ic_ch", "IC_{ch}-2015 EB --- ABSOLUTE IC", 500, 0.5, 1.5)
histos["EB_ic_ch_corr"]=ROOT.TH1F("histo_EB_ic_ch_corr", "IC_{ch_corr}-2015 EB --- ABSOLUTE IC", 500, 0.5, 1.5)
histos["EB_ratio_ic_ch"]=ROOT.TH1F("histo_EB_ratio_ic_ch", "IC_{ch}-2015 / IC-2012D EB --- ABSOLUTE IC", 500, 0.8, 1.2)
histos["EB_ratio_ic_ch_corr"]=ROOT.TH1F("histo_EB_ratio_ic_ch_corr", "IC_{ch_corr}-2015 / IC-2012D EB --- ABSOLUTE IC", 500, 0.8, 1.2)
# EE
histos["EE_ic_ch"]=ROOT.TH1F("histo_EE_ic_ch", "IC_{ch}-2015 EE --- ABSOLUTE IC", 500, 0.3, 1.7)
histos["EE_ratio_ic_ch"]=ROOT.TH1F("histo_EE_ratio_ic_ch", "IC_{ch}-2015 / IC-2012D EE --- ABSOLUTE IC", 500, 0.5, 1.5)

## Profiles
# EB
profiles["prPhi_EB_ic_ch"]=ROOT.TH1F("prPhi_EB_ic_ch", "EB ICs with channel based k-factors --- PROFILE #phi;i#phi",
                                     360, -0.5, 359.5)
profiles["prPhi_EB_ic_ring"]=ROOT.TH1F("prPhi_EB_ic_ring", "EB ICs with ring based k-factors --- PROFILE #phi;i#phi",
                                       360, -0.5, 359.5)
profiles["prEta_EB_ic_ch"]=ROOT.TH1F("prEta_EB_ic_ch", "EB ICs with channel based k-factors --- PROFILE #eta;i#eta",
                                     171, -85.5, 85.5)
profiles["prEta_EB_ic_ring"]=ROOT.TH1F("prEta_EB_ic_ring", "EB ICs with ring based k-factors --- PROFILE #eta;i#eta",
                                       171, -85.5, 85.5)
profiles["prPhi_EB_err_ic_ch"]=ROOT.TH1F("prPhi_EB_err_ic_ch", "EB ICs statistical error --- PROFILE #phi;i#phi",
                                         360, -0.5, 359.5)
profiles["prEta_EB_err_ic_ch"]=ROOT.TH1F("prEta_EB_err_ic_ch", "EB ICs statistical error --- PROFILE #eta;i#eta",
                                         171, -85.5, 85.5)
profiles["prEta_EE_err_ic_ch"]=ROOT.TH1F("prEta_EE_err_ic_ch", "EE ICs statistical error --- PROFILE #eta;i#eta",
                                         81, -40.5, 40.5)

## maps
# EB
maps["EB_ic_ring"]=ROOT.TH2F("map_EB_ic_ring", "EB ICs with ring-based k-factors --- MAP", 360, -0.5, 359.5, 171, -85.5, 85.5)
maps["EB_ic_ch"]=ROOT.TH2F("map_EB_ic_ch", "EB ICs with channel-based k-factors --- MAP", 360, -0.5, 359.5, 171, -85.5, 85.5)
maps["EB_ic_ch_corr"]=ROOT.TH2F("map_EB_ic_ch_corr", "EB ICs with channel-based k-factors and geo corrections --- MAP",
                                360, -0.5, 359.5, 171, -85.5, 85.5)
maps["EB_ic_diff"]=ROOT.TH2F("map_EB_ic_diff", "EB ICs relative difference (ch-ring)/ring--- MAP",
                                   360, -0.5, 359.5, 171, -85.5, 85.5)
maps["EB_k_ring"]=ROOT.TH2F("map_EB_k_ring", "EB ring-based k-factors --- MAP", 360, -0.5, 359.5, 171, -85.5, 85.5)
maps["EB_k_ch"]=ROOT.TH2F("map_EB_k_ch", "EB channel-based k-factors --- MAP", 360, -0.5, 359.5, 171, -85.5, 85.5)
maps["EB_k_diff"]=ROOT.TH2F("map_EB_k_diff", "EB k-factors relative difference (ch-ring)/ring --- MAP",
                                  360, -0.5, 359.5, 171, -85.5, 85.5)
maps["EB_n_hits"]=ROOT.TH2F("map_EB_n_hits", "Number of hits/lumi EB--- MAP", 360, -0.5, 359.5, 171, -85.5, 85.5)
maps["EB_old_ic"]=ROOT.TH2F("map_EB_old_ic", "IC-2012D EB--- MAP", 360, -0.5, 359.5, 171, -85.5, 85.5)
maps["EB_ratio_ic_ring"]=ROOT.TH2F("map_EB_ratio_ic_ring", "IC_{ring}-2015 / IC-2012D EB--- MAP", 360, -0.5, 359.5, 171, -85.5, 85.5)
maps["EB_ratio_ic_ch"]=ROOT.TH2F("map_EB_ratio_ic_ch", "IC_{ch}-2015 / IC-2012D EB--- MAP", 360, -0.5, 359.5, 171, -85.5, 85.5)
maps["EB_ratio_ic_ch_corr"]=ROOT.TH2F("map_EB_ratio_ic_ch_corr", "IC_{ch_corr}-2015 / IC-2012D EB--- MAP",
                                      360, -0.5, 359.5, 171, -85.5, 85.5)
# EE plus
maps["EEp_ic_ring"]=ROOT.TH2F("map_EEp_ic_ring", "EE plus ICs with ring-based k-factors --- MAP",
                                    100, 0.5, 100.5, 100, 0.5, 100.5)
maps["EEp_ic_ch"]=ROOT.TH2F("map_EEp_ic_ch", "EE plus ICs with channel-based k-factors --- MAP",
                                  100, 0.5, 100.5, 100, 0.5, 100.5)
maps["EEp_ic_diff"]=ROOT.TH2F("map_EEp_ic_diff", "EE plus ICs relative difference (ch-ring)/ring --- MAP",
                                    100, 0.5, 100.5, 100, 0.5, 100.5)
maps["EEp_k_ring"]=ROOT.TH2F("map_EEp_k_ring", "EE plus ring-based k-factors --- MAP", 100, 0.5, 100.5, 100, 0.5, 100.5)
maps["EEp_k_ch"]=ROOT.TH2F("map_EEp_k_ch", "EE plus channel-based k-factors --- MAP", 100, 0.5, 100.5, 100, 0.5, 100.5)
maps["EEp_k_diff"]=ROOT.TH2F("map_EEp_k_diff", "EE plus k-factors relative difference (ch-ring)/ring --- MAP",
                               100, 0.5, 100.5, 100, 0.5, 100.5)
maps["EEp_n_hits"]=ROOT.TH2F("map_EEp_n_hits", "Number of hits/lumi EE+ --- MAP", 100, 0.5, 100.5, 100, 0.5, 100.5)
maps["EEp_ratio_ic_ring"]=ROOT.TH2F("map_EEp_ratio_ic_ring", "IC_{ring}-2015 / IC-2012D EE+ --- MAP",
                                    100, 0.5, 100.5, 100, 0.5, 100.5)
maps["EEp_ratio_ic_ch"]=ROOT.TH2F("map_EEp_ratio_ic_ch", "IC_{ch}-2015 / IC-2012D EE+ --- MAP",
                                    100, 0.5, 100.5, 100, 0.5, 100.5)
# EE minus
maps["EEm_ic_ring"]=ROOT.TH2F("map_EEm_ic_ring", "EE minus ICs with channel-based k-factors --- MAP",
                                    100, 0.5, 100.5, 100, 0.5, 100.5)
maps["EEm_ic_ch"]=ROOT.TH2F("map_EEm_ic_ch", "EE minus ICs with channel-based k-factors --- MAP",
                                  100, 0.5, 100.5, 100, 0.5, 100.5)
maps["EEm_ic_diff"]=ROOT.TH2F("map_EEm_ic_diff", "EE minus ICs relative difference (ch-ring)/ring --- MAP",
                                    100, 0.5, 100.5, 100, 0.5, 100.5)
maps["EEm_k_ring"]=ROOT.TH2F("map_EEm_k_ring", "EE minus ring-based k-factors --- MAP", 100, 0.5, 100.5, 100, 0.5, 100.5)
maps["EEm_k_ch"]=ROOT.TH2F("map_EEm_k_ch", "EE minus channel-based k-factors --- MAP", 100, 0.5, 100.5, 100, 0.5, 100.5)
maps["EEm_k_diff"]=ROOT.TH2F("map_EEm_k_diff", "EE minus k-factors relative difference (ch-ring)/ring --- MAP",
                               100, 0.5, 100.5, 100, 0.5, 100.5)
maps["EEm_n_hits"]=ROOT.TH2F("map_EEm_n_hits", "Number of hits/lumi EE- --- MAP", 100, 0.5, 100.5, 100, 0.5, 100.5)
maps["EEm_ratio_ic_ring"]=ROOT.TH2F("map_EEm_ratio_ic_ring", "IC_{ring}-2015 / IC-2012D EE- --- MAP",
                                    100, 0.5, 100.5, 100, 0.5, 100.5)
maps["EEm_ratio_ic_ch"]=ROOT.TH2F("map_EEm_ratio_ic_ch", "IC_{ch}-2015 / IC-2012D EE- --- MAP",
                                    100, 0.5, 100.5, 100, 0.5, 100.5)

## VARIABLES map
varReMap={}
varReMap["err_ic_ch"]="ic_ch_err*ic_ch"
varReMap["err_ic_ring"]="ic_ring_err*ic_ring"

## geometrical corr
ones = [1 for i in range (61200)]
ebUnCorrIC= [0 for i in range (61200)]
ebOldIC={}
ebAbsIC={}
for ieta in range (-85,86):
    ebOldIC[str(ieta)]=ROOT.TH1F(str(ieta)+"old", "", 100, 0.5, 1.5)
    ebAbsIC[str(ieta)]=ROOT.TH1F(str(ieta)+"abs", "", 100, 0.5, 1.5)
    
while ebTree.NextEntry():
    if opts.block != -1 and ebTree.block != opts.block:
        continue
    ebUnCorrIC[hashedIndex(ebTree.ieta, ebTree.iphi)]=(ebTree.ic_ch)    
    ebOldIC[str(ebTree.ieta)].Fill(ebTree.ic_old)
    ebAbsIC[str(ebTree.ieta)].Fill(ebTree.ic_abs*ebTree.ic_ch)

diffp,rmsp,diffm,rmsm= calculatedifferences(ebUnCorrIC, ones)
geoCorrIC = applycorrections(ebUnCorrIC, diffp, diffm)

ebICCorrOld={}
ebICCorrAbs={}
for ring in ebOldIC.keys():
    ebICCorrOld[ring]=ebOldIC[ring].GetMean()
    ebICCorrAbs[ring]=ebAbsIC[ring].GetMean()

## FILL THE HISTOS ##
## EB
# 1D & 2D
while ebTree.NextEntry():
    if opts.block != -1 and ebTree.block != opts.block:
        continue
    if ebTree.rec_hit.GetSumEt() == 0:
        continue    
    if nLumis==-1:
        nLumis = ebTree.n_lumis
    maps["EB_ic_ring"].Fill(ebTree.iphi, ebTree.ieta, ebTree.ic_ring)
    maps["EB_ic_ch"].Fill(ebTree.iphi, ebTree.ieta, ebTree.ic_ch)
    maps["EB_ic_diff"].Fill(ebTree.iphi, ebTree.ieta, (ebTree.ic_ch-ebTree.ic_ring)/ebTree.ic_ring)
    maps["EB_k_ring"].Fill(ebTree.iphi, ebTree.ieta, ebTree.k_ring)
    maps["EB_k_ch"].Fill(ebTree.iphi, ebTree.ieta, ebTree.k_ch)
    maps["EB_k_diff"].Fill(ebTree.iphi, ebTree.ieta, (ebTree.k_ch-ebTree.k_ring)/ebTree.k_ring)
    maps["EB_n_hits"].Fill(ebTree.iphi, ebTree.ieta, ebTree.rec_hit.GetNhits()/ebTree.n_lumis)
    if ebTree.ic_old != 0 and ebTree.ic_ch != 0:
        abs_ic_norm = ebTree.ic_abs*ebTree.ic_ch/ebICCorrAbs[str(ebTree.ieta)]
        old_ic_norm = ebTree.ic_old/ebICCorrOld[str(ebTree.ieta)]
        maps["EB_old_ic"].Fill(ebTree.iphi, ebTree.ieta, ebTree.ic_old)
        maps["EB_ratio_ic_ring"].Fill(ebTree.iphi, ebTree.ieta, ebTree.ic_abs*ebTree.ic_ring/ebTree.ic_old)
        maps["EB_ratio_ic_ch"].Fill(ebTree.iphi, ebTree.ieta, abs_ic_norm/old_ic_norm) 
        maps["EB_ratio_ic_ch_corr"].Fill(ebTree.iphi, ebTree.ieta, ebTree.ic_abs*geoCorrIC[hashedIndex(ebTree.ieta, ebTree.iphi)]/
                                         ebTree.ic_old)
        maps["EB_ic_ch_corr"].Fill(ebTree.iphi, ebTree.ieta, ebTree.ic_abs*geoCorrIC[hashedIndex(ebTree.ieta, ebTree.iphi)])
        histos["EB_ic_ch"].Fill(ebTree.ic_abs*ebTree.ic_ch)
        histos["EB_ic_ch_corr"].Fill(ebTree.ic_abs*geoCorrIC[hashedIndex(ebTree.ieta, ebTree.iphi)])    
        histos["EB_ratio_ic_ch"].Fill(abs_ic_norm/old_ic_norm)
        histos["EB_ratio_ic_ch_corr"].Fill(ebTree.ic_abs*geoCorrIC[hashedIndex(ebTree.ieta, ebTree.iphi)]/ebTree.ic_old)

# profiles
tmpHisto=ROOT.TH1F("tmp", "tmp", 1000, -100, 100)
for key in profiles:
    if "Phi" in key:
        varRange = range(0, 360)
        cutBase = "iphi=="
    if "Eta" in key:
        varRange = range(-85, 85)
        cutBase = "ieta=="
    if "EE" in key:
        varRange = range(-40, 40)
        cutBase = "iring=="
    for prVar in varRange:
        var = key[9:]+">>tmp"
        if key[9:] in varReMap.keys():
            var = varReMap[key[9:]]+">>tmp"
        cut = cutBase+str(prVar)
        if opts.block != -1:
            cut += " && block=="+str(opts.block)
        if "EE" in key:
            eeTree.Draw(var, cut, "goff")
        else:
            ebTree.Draw(var, cut, "goff")
        profiles[key].SetBinContent(profiles[key].FindBin(prVar), tmpHisto.GetMean())
        profiles[key].SetBinError(profiles[key].FindBin(prVar), tmpHisto.GetRMS())

# ugly workaround to get white bins for non existing channels
for xbin in range(0, 101):
    for ybin in range(0, 101):
        for key in maps.keys():
            if "EE" in key:
                if not ROOT.EEDetId.validDetId(xbin, ybin, 1):
                    maps[key].Fill(xbin, ybin, -10)
                
while eeTree.NextEntry():
    if opts.block != -1 and eeTree.block != opts.block:
        continue
    if eeTree.iring > 0:
        subdet = "EEp_"
    else:
        subdet = "EEm_"
    maps[subdet+"ic_ring"].Fill(eeTree.ix, eeTree.iy, eeTree.ic_ring)
    maps[subdet+"ic_ch"].Fill(eeTree.ix, eeTree.iy, eeTree.ic_ch)
    maps[subdet+"ic_diff"].Fill(eeTree.ix, eeTree.iy, (eeTree.ic_ch-eeTree.ic_ring)/eeTree.ic_ring)
    maps[subdet+"k_ring"].Fill(eeTree.ix, eeTree.iy, eeTree.k_ring)
    maps[subdet+"k_ch"].Fill(eeTree.ix, eeTree.iy, eeTree.k_ch)
    maps[subdet+"k_diff"].Fill(eeTree.ix, eeTree.iy, (eeTree.k_ch-eeTree.k_ring)/eeTree.k_ring)
    maps[subdet+"n_hits"].Fill(eeTree.ix, eeTree.iy, eeTree.rec_hit.GetNhits())#/eeTree.n_lumis)
    if eeTree.ic_old>0 :
        maps[subdet+"ratio_ic_ring"].Fill(eeTree.ix, eeTree.iy, eeTree.ic_abs*eeTree.ic_ring/eeTree.ic_old)
        maps[subdet+"ratio_ic_ch"].Fill(eeTree.ix, eeTree.iy, eeTree.ic_abs*eeTree.ic_ch/eeTree.ic_old)
        histos["EE_ic_ch"].Fill(eeTree.ic_abs*eeTree.ic_ch)
        histos["EE_ratio_ic_ch"].Fill(eeTree.ic_abs*eeTree.ic_ch/eeTree.ic_old)

## DRAW PLOTS ##    
# Define z-axis ranges: with regexp matching ;)
ranges={}
ranges[re.compile("^EB_ic_((?!diff).)*$")]=[0.9, 1.1]
ranges[re.compile("^EB_ic_diff$")]=[-0.001, 0.001]
ranges[re.compile(".*EB_ratio.*")]=[0.9, 1.1]
ranges[re.compile("prEta_EB_err.*")]=[0, 0.0001]
ranges[re.compile(".*EB_k_((?!diff).)*$")]=[1.5, 2.6]
ranges[re.compile(".*EB_k_diff$")]=[-0.05, 0.05]
ranges[re.compile("^EE.*_ic_((?!diff).)*$")]=[0.5, 1.5]
ranges[re.compile("^EE.*_ic_diff$")]=[-0.01, 0.01]
ranges[re.compile("prEta_EE_err.*")]=[0, 0.0005]
ranges[re.compile(".*EE.*_k_((?!diff).)*$")]=[1, 1.7]
ranges[re.compile(".*EE.*_k_diff$")]=[-0.05, 0.05]
ranges[re.compile(".*_n_hits$")]=[0, 500]
ranges[re.compile(".*_old_ic$")]=[0.9, 1.1]

# Define axis titles
axisTitles={}
axisTitles[re.compile("^EB.*")]=["i#phi", "i#eta"]
axisTitles[re.compile("^EE.*")]=["ix", "iy"]
axisTitles[re.compile("^pr.*err$")]="#sigma_{IC}/IC"

outFile=ROOT.TFile(opts.outdir+"phiSymCalibCheck.root", "RECREATE")
# draw MAPS
for key in maps.keys():
    for title in axisTitles:
        if title.match(key):
            maps[key].GetXaxis().SetTitle(axisTitles[title][0])
            maps[key].GetYaxis().SetTitle(axisTitles[title][1])
    for cat in ranges.keys():
        if cat.match(key):
            maps[key].SetAxisRange(ranges[cat][0], ranges[cat][1], "Z")
    maps[key].Write()
    # png, pdf
    if opts.plots:
        if "EB" in key:
            canvas=ROOT.TCanvas(key, key, 800, 380)
        else:
            canvas=ROOT.TCanvas(key, key, 900, 800)
        maps[key].Draw("COLZ")
        canvas.Print(opts.outdir+maps[key].GetName()+".png", "png")
        canvas.Print(opts.outdir+maps[key].GetName()+".pdf", "pdf")

# draw PROFILES
for key in profiles.keys():    
    profiles[key].GetYaxis().SetTitle(key[9:])
    for title in axisTitles:
        if title.match(key):
            profiles[key].GetYaxis().SetTitle(axisTitles[title])
    for cat in ranges.keys():
        if cat.match(key):
            profiles[key].SetAxisRange(ranges[cat][0], ranges[cat][1], "Y")
    profiles[key].Write()
    # png. pdf
    if opts.plots:
        canvas=ROOT.TCanvas(key, key, 800, 400)
        canvas.SetRightMargin(0.07)
        profiles[key].Draw()
        canvas.Print(opts.outdir+profiles[key].GetName()+".png", "png")
        canvas.Print(opts.outdir+profiles[key].GetName()+".pdf", "pdf")

# draw histos
for key in histos.keys():    
    # histos[key].GetYaxis().SetTitle(key[9:])
    # for title in axisTitles:
    #     if title.match(key):
    #         histos[key].GetYaxis().SetTitle(axisTitles[title])
    # for cat in ranges.keys():
    #     if cat.match(key):
    #         histos[key].SetAxisRange(ranges[cat][0], ranges[cat][1], "Y")
    histos[key].Write()
    # png. pdf
    if opts.plots:
        canvas=ROOT.TCanvas(key, key, 800, 400)
        canvas.SetRightMargin(0.07)
        ROOT.gStyle.SetOptStat("mr")
        ROOT.gStyle.SetOptFit(111)
        ROOT.gROOT.ForceStyle()
        histos[key].Fit("gaus")
        histos[key].Draw()
        canvas.Print(opts.outdir+histos[key].GetName()+".png", "png")
        canvas.Print(opts.outdir+histos[key].GetName()+".pdf", "pdf")
        
outFile.Close()

print "SUMMARY:"
print "--------"
print "analyzed lumis: "+str(nLumis)

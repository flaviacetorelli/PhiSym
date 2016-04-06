from __future__ import division
import lxml.etree as ET
import sys
import re
import time
import argparse
oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
import ROOT
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv

from PhiSym.EcalCalibAlgos.EcalCalibAnalysis import *

parser = argparse.ArgumentParser (description = 'Dump PhiSym IC correction from a ROOT file')
parser.add_argument('inputfile' , default="phisym_intercalibs.root", help='analyze this file')
parser.add_argument('-b', '--block' , default=-1, type=int, help='analyze this block')
parser.add_argument('-v', '--version' , default=-1, type=int, help='ic version')
parser.add_argument('-f', '--outputFile' , default="phisym_ic.xml", help='store the ICs in this file')
parser.add_argument('-c', '--correctionsFile' , default="", help='load geo&material correction from this file')
parser.add_argument('-k', '--kType' , default="ch", help='k-Factor computation type: ch / ring')
parser.add_argument('--rel' , action="store_true", default=False, help='dump phisym ICs correction only')

opts = parser.parse_args ()

# Load the CalibrationFile format
ROOT.gSystem.Load("libDataFormatsEcalDetId.so");
ROOT.gSystem.Load("libPhiSymEcalCalibDataFormats.so");
ROOT.AutoLibraryLoader.enable()

# Get the crystals trees
inFile = ROOT.TFile(opts.inputfile)
bareTree = inFile.Get("eb_xstals")
ebTree = ROOT.CrystalsEBTree(bareTree)
bareTree = inFile.Get("ee_xstals")
eeTree = ROOT.CrystalsEETree(bareTree)

ebIC=[-999 for i in range(ROOT.EBDetId.kSizeForDenseIndexing)]
eeIC=[-999 for i in range(ROOT.EEDetId.kSizeForDenseIndexing)]
ebK=[-999 for i in range(ROOT.EBDetId.kSizeForDenseIndexing)]
eeK=[-999 for i in range(ROOT.EEDetId.kSizeForDenseIndexing)]
ebN=[0 for i in range(ROOT.EBDetId.kSizeForDenseIndexing)]
eeN=[0 for i in range(ROOT.EEDetId.kSizeForDenseIndexing)]
ebICerr=[999 for i in range(ROOT.EBDetId.kSizeForDenseIndexing)]
eeICerr=[999 for i in range(ROOT.EEDetId.kSizeForDenseIndexing)]
ebICsys=[999 for i in range(ROOT.EBDetId.kSizeForDenseIndexing)]
eeICsys=[999 for i in range(ROOT.EEDetId.kSizeForDenseIndexing)]
ebCorr=[1 for i in range(ROOT.EBDetId.kSizeForDenseIndexing)]
eeCorr=[1 for i in range(ROOT.EEDetId.kSizeForDenseIndexing)]

#EB
while ebTree.NextEntry():
    if opts.block != -1 and ebTree.block != opts.block:
        continue;
    index = ROOT.EBDetId(ebTree.ieta, ebTree.iphi).hashedIndex()    
    if opts.kType == "ch":
        ic = ebTree.ic_ch
        ic_err = ebTree.ic_ch_err
    else:
        ic = ebTree.ic_ring
        ic_err = ebTree.ic_ring_err
    if not opts.rel:
        ic = ebTree.ic_abs*ic
        ic_err = ebTree.ic_abs*ic_err
    if ic <= 0:
        ic = -1
        ic_err = 999    

    ebIC[index] = ic
    ebICerr[index] = ic_err
    ebICsys[index] = ebTree.ic_err_sys
    ebK[index] = ebTree.k_ch
    ebN[index] = ebTree.rec_hit.GetNhits()

#EE
while eeTree.NextEntry():
    if opts.block != -1 and eeTree.block != opts.block:
        continue;
    if eeTree.iring > 0:
        side = 1
    else:
        side = -1
    index = ROOT.EEDetId(eeTree.ix, eeTree.iy, side).hashedIndex()    
    if opts.kType == "ch":
        ic = eeTree.ic_ch
        ic_err = eeTree.ic_ch_err
    else:
        ic = eeTree.ic_ring
        ic_err = eeTree.ic_ring_err
    if not opts.rel:
        ic = eeTree.ic_abs*ic
        ic_err = eeTree.ic_abs*ic_err
    if ic <= 0:
        ic = -1
        ic_err = 999    

    eeIC[index] = ic
    eeICerr[index] = ic_err
    eeICsys[index] = eeTree.ic_err_sys
    eeK[index] = eeTree.k_ch
    eeN[index] = eeTree.rec_hit.GetNhits()

# read correction file if specified
with open(opts.correctionsFile) as corrections:
    channels = corrections.readlines()
    for channel in channels:
        tokens = channel.split()
        if int(tokens[2]) == 0:
            ebCorr[ROOT.EBDetId(int(tokens[0]), int(tokens[1])).hashedIndex()] = float(tokens[3])
        elif ROOT.EEDetId.validDetId(int(tokens[0]), int(tokens[1]), int(tokens[2])):
            eeCorr[ROOT.EEDetId(int(tokens[0]), int(tokens[1]), int(tokens[2])).hashedIndex()] = float(tokens[3])

# write the XML file
container = ET.Element("EcalFloatCondObjectContainer")
header = ET.SubElement(container, "EcalCondHeader")
ET.SubElement(header, "method").text = "PhiSymmetry"
ET.SubElement(header, "version").text = str(opts.version)
ET.SubElement(header, "datasource").text = opts.inputfile
ET.SubElement(header, "since").text = "1"
ET.SubElement(header, "tag").text = "unknown"
ET.SubElement(header, "date").text = time.strftime("%c")

for index in range(61200):    
    if ebN[index]==0 or ebIC[index] <= 0:
        ebIC[index] = -1
        ebICerr[index] = 999
        ebICsys[index] = 999

    ieta = ROOT.EBDetId(ROOT.EBDetId.detIdFromDenseIndex(index)).ieta()
    iphi = ROOT.EBDetId(ROOT.EBDetId.detIdFromDenseIndex(index)).iphi()
    if ebICerr[index] == 999:
        err = 999
    else:
        err = sqrt(ebICerr[index]*ebICerr[index] + ebICsys[index]*ebICsys[index])

    cell = ET.SubElement(container, "cell", iEta=str(ieta), iPhi=str(iphi))
    ET.SubElement(cell, "Value").text = str(ebIC[index]*ebCorr[index])

for index in range(14648):    
    if ebN[index]==0 or eeIC[index] <= 0:
        eeIC[index] = -1
        eeICerr[index] = 999
        eeICsys[index] = 999
        
    ix = ROOT.EEDetId(ROOT.EEDetId.detIdFromDenseIndex(index)).ix()
    iy = ROOT.EEDetId(ROOT.EEDetId.detIdFromDenseIndex(index)).iy()
    zside = ROOT.EEDetId(ROOT.EEDetId.detIdFromDenseIndex(index)).zside()
    if eeICerr[index] == 999:
        err = 999
    else:
        err = sqrt(eeICerr[index]*eeICerr[index] + eeICsys[index]*eeICsys[index])

    cell = ET.SubElement(container, "cell", ix=str(ix), iy=str(iy), zside=str(zside))
    ET.SubElement(cell, "Value").text = str(eeIC[index]*eeCorr[index])
    
# write the xml tree to tho out file
out = open(opts.outputFile, "w")
out.write(ET.tostring(container, pretty_print=True))
out.close()

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

parser = argparse.ArgumentParser (description = 'Dump PhiSym IC correction from a ROOT file')
parser.add_argument('inputfile' , default="phisym_intercalibs.root", help='analyze this file')
parser.add_argument('-b', '--block' , default=-1, type=int, help='analyze this block')
parser.add_argument('-f', '--outputFile' , default=-1, type=int, help='store the ICs in this file')
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

ebIC=[-999 for i in range(61200)]
eeIC=[-999 for i in range(14648)]
ebK=[-999 for i in range(61200)]
eeK=[-999 for i in range(14648)]
ebN=[0 for i in range(61200)]
eeN=[0 for i in range(14648)]
ebICerr=[999 for i in range(61200)]
eeICerr=[999 for i in range(14648)]
ebICsys=[999 for i in range(61200)]
eeICsys=[999 for i in range(14648)]
ebCorr=[1 for i in range(61200)]
eeCorr=[1 for i in range(14648)]

# read correction file if specified
with open(opts.correctionsFile) as corrections:
    channels = corrections.readlines()
    for channel in channels:
        tokens = channel.split()
        if int(tokens[2]) == 0:
            ebCorr[ROOT.EBDetId(int(tokens[0]), int(tokens[1])).hashedIndex()] = float(tokens[3])
        elif ROOT.EEDetId.validDetId(int(tokens[0]), int(tokens[1]), int(tokens[2])):
            eeCorr[ROOT.EEDetId(int(tokens[0]), int(tokens[1]), int(tokens[2])).hashedIndex()] = float(tokens[3])

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

print "f"
print "# status = used/unused (1/0) in ring sumEt avarage computation, or bad channel (-1)"
print "# ieta(ix) -- iphi(iy) -- zside -- IC -- IC_err -- IC_err_stat -- IC_err_sys -- status -- n_hits -- k-factor"
    
for index in range(61200):    
    if ebN[index]==0 or ebIC[index] <= 0:
        ebIC[index] = -1
        status = -1
        ebICerr[index] = 999
        ebICsys[index] = 999
    else:
        status = 1

    ieta = ROOT.EBDetId(ROOT.EBDetId.detIdFromDenseIndex(index)).ieta()
    iphi = ROOT.EBDetId(ROOT.EBDetId.detIdFromDenseIndex(index)).iphi()
    if ebICerr[index] == 999:
        err = 999
    else:
        err = sqrt(ebICerr[index]*ebICerr[index] + ebICsys[index]*ebICsys[index])
    print ieta, iphi, 0, "%.5f" % (ebIC[index]*ebCorr[index]), "%.7f" % err, "%.7f" % ebICerr[index], "%.7f" % ebICsys[index],
    print status, ebN[index], ebK[index]

for index in range(14648):    
    if ebN[index]==0 or eeIC[index] <= 0:
        eeIC[index] = -1
        status = -1
        eeICerr[index] = 999
        eeICsys[index] = 999
    else:
        status = 1
        
    ix = ROOT.EEDetId(ROOT.EEDetId.detIdFromDenseIndex(index)).ix()
    iy = ROOT.EEDetId(ROOT.EEDetId.detIdFromDenseIndex(index)).iy()
    zside = ROOT.EEDetId(ROOT.EEDetId.detIdFromDenseIndex(index)).zside()
    if eeICerr[index] == 999:
        err = 999
    else:
        err = sqrt(eeICerr[index]*eeICerr[index] + eeICsys[index]*eeICsys[index])
    print ix, iy, zside, "%.5f" % (eeIC[index]*eeCorr[index]), "%.7f" % err, "%.7f" % eeICerr[index], "%.7f" % eeICsys[index],
    print status, eeN[index], eeK[index]


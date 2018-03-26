import subprocess
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing


debug = True
eosdirs = ['crab_PHISYM-CMSSW_9_4_0-multifit-94X_dataRun2_ReReco_EOY17_v2-IC2017_RunB_v9_Merged_v2/180206_104913/0000/', 'crab_PHISYM-CMSSW_9_4_0-multifit-94X_dataRun2_ReReco_EOY17_v2-IC2017_RunC_v9_Merged_v2/180206_105059/0000/', 'crab_PHISYM-CMSSW_9_4_0-multifit-94X_dataRun2_ReReco_EOY17_v2-IC2017_RunD_v9_Merged_v2/180206_105134/0000/', 'crab_PHISYM-CMSSW_9_4_0-multifit-94X_dataRun2_ReReco_EOY17_v2-IC2017_RunE_v9_Merged_v2/180206_105212/0000/']

files = []
for eosdir in eosdirs:
    if eosdir[-1] != '/':
        eosdir += '/'
    if "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/phiSymmetry/AlCaPhiSym/" not in eosdir:
        eosdir = "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/phiSymmetry/AlCaPhiSym/"+eosdir
    print('>> Creating list of files from: \n'+eosdir)
    lsCmd = subprocess.Popen(['ls '+eosdir+'*.root'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    str_files, err = lsCmd.communicate()
    files.extend(['root://eoscms/'+ifile for ifile in str_files.split("\n")])
    files.pop()
    if debug:
        for ifile in files:
            print(ifile)
    
process = cms.Process('Calibration')


process.IOVBounds = cms.PSet(
    startingIOV     = cms.int32(0),
    nIOVs           = cms.int32(-1),
    manualSplitting = cms.bool(True),
    beginRuns       = cms.vint32(0),
    endRuns         = cms.vint32(9999999 ),
    IOVMaps         = cms.vstring('')
)

if debug:
    print(process.IOVBounds)

process.ioFilesOpt = cms.PSet(    
    outputFile = cms.string('summed_'),
    
    oldConstantsFiles = cms.vstring(''),
    
    recoConstantsFile = cms.string('/afs/cern.ch/work/s/spigazzi/ECAL/CMSSW_8_0_17/src/PhiSym/EcalCalibAlgos/data/EcalIntercalibConstants_Prompt2016.dat'),
    
    inputFiles = cms.vstring(files)
)


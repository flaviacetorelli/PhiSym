import subprocess
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')
options.register('eosdirs',
                 '',
                 VarParsing.multiplicity.list,
                 VarParsing.varType.string,
                 "Merged ntuples location(s) on EOS")
options.register('iovmaps',
                 '',
                 VarParsing.multiplicity.list,
                 VarParsing.varType.string,
                 "List of IOVMap files")
options.register('iovbounds',
                 '',
                 VarParsing.multiplicity.list,
                 VarParsing.varType.int,
                 "List of IOV bounds: begin1, end1, begin2, end2, ...")
options.register('firstiov',
                 0,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 "First IOV")
options.register('niovs',
                 -1,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 "Number of IOVs")
options.register('debug',
                 True,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Print debug messages")
options.parseArguments()

files = []
for eosdir in options.eosdirs:
    if eosdir[-1] != '/':
        eosdir += '/'
    if "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/phiSymmetry/AlCaPhiSym/" not in eosdir:
        eosdir = "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/phiSymmetry/AlCaPhiSym/"+eosdir
    print('>> Creating list of files from: \n'+eosdir)
    lsCmd = subprocess.Popen(['ls '+eosdir+'*.root'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    str_files, err = lsCmd.communicate()
    files.extend(['root://eoscms/'+ifile for ifile in str_files.split("\n")])
    files.pop()
    if options.debug:
        for ifile in files:
            print(ifile)
    
process = cms.Process('Calibration')

manual_splitting = True if len(options.iovbounds) > 0 else False

process.IOVBounds = cms.PSet(
    startingIOV     = cms.int32(options.firstiov),
    nIOVs           = cms.int32(options.niovs),
    manualSplitting = cms.bool(manual_splitting),
    beginRuns       = cms.vint32( [options.iovbounds[i] for i in range(0, len(options.iovbounds), 2)] ),
    endRuns         = cms.vint32( [options.iovbounds[i] for i in range(1, len(options.iovbounds), 2)] ),
    IOVMaps         = cms.vstring(options.iovmaps)
)

if options.debug:
    print(process.IOVBounds)
    print 'ciao'
process.ioFilesOpt = cms.PSet(    
    outputFile = cms.string('summed_'),

    
    oldConstantsFiles = cms.vstring(''),
    
    recoConstantsFile = cms.string('/afs/cern.ch/work/s/spigazzi/ECAL/CMSSW_8_0_17/src/PhiSym/EcalCalibAlgos/data/EcalIntercalibConstants_Prompt2016.dat'),
    
    inputFiles = cms.vstring(files)

)

print 'ciao'


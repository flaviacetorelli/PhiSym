import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.AlCa.GlobalTag import GlobalTag

process = cms.Process("PHISYMstep2")

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100000

process.load('Configuration.Geometry.GeometryExtended2015Reco_cff')

# skip bad events
process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

# Global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag = GlobalTag(process.GlobalTag, '74X_dataRun2_Prompt_v4')

# Input source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                                #'/store/user/spigazzi/AlCaPhiSym/crab_PHISYM-CMSSW_746-weights-74X_dataRun2_Prompt_v0-Run2015B_v1/150716_200712/0000/phisym_weights_1lumis_1.root'
                                'root://xrootd-cms.infn.it:1194//store/user/spigazzi/AlCaPhiSym/crab_PHISYM-CMSSW_7415-weights-74X_dataRun2_Prompt_v4-Run2015B_v1/151109_102420/0000/phisym_weights_1lumis_219.root'
                            )
)                                

# PHISYM Calib
process.load('PhiSym.EcalCalibAlgos.PhiSymMerger_cfi')
process.PhiSymMerger.blocksToSum = 1000
process.PhiSymMerger.IOVfile = cms.untracked.string("IOVmap.root")

# Output TFile
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("phisym_merged.root"))

process.path = cms.Path(process.PhiSymMerger)

process.schedule = cms.Schedule(process.path)

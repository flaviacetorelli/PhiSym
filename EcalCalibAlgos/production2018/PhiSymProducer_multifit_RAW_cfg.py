import subprocess
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.AlCa.GlobalTag import GlobalTag

# parse commad line options
options = VarParsing('analysis')
options.maxEvents = -1
options.outputFile = 'phisym_multifit_1lumis.root'
options.register('datasets',
                 '',
                 VarParsing.multiplicity.list,
                 VarParsing.varType.string,
                 "Input dataset(s)")
options.register('debug',
                 False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Print debug messages")
options.parseArguments()

process=cms.Process("PHISYM")

rawTag    = cms.InputTag('source')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.Geometry.GeometryExtended2017Reco_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('RecoLuminosity.LumiProducer.bunchSpacingProducer_cfi')
process.load('RecoLocalCalo.EcalRecProducers.ecalMultiFitUncalibRecHit_cfi')
process.load('RecoLocalCalo.EcalRecProducers.ecalUncalibRecHit_cfi')
process.load('RecoLocalCalo.EcalRecProducers.ecalRecHit_cfi')
process.load('RecoVertex.BeamSpotProducer.BeamSpot_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("RecoLocalCalo.Configuration.hcalLocalReco_cff")
process.load("RecoLocalCalo.Configuration.ecalLocalRecoSequence_cff")
process.load("EventFilter.HcalRawToDigi.HcalRawToDigi_cfi")
process.load("EventFilter.EcalRawToDigi.EcalUnpackerData_cfi")
process.load("RecoLuminosity.LumiProducer.bunchSpacingProducer_cfi")
process.load('FWCore/MessageService/MessageLogger_cfi')

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.cerr.default = cms.untracked.PSet(
    limit = cms.untracked.int32(10000000),
    reportEvery = cms.untracked.int32(5000)
)

# import of standard configurations
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# skip bad events
process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound'),
)

# Input source
files = []
for dataset in options.datasets:
    print('>> Creating list of files from: \n'+dataset)
    query = "--query='file instance=prod/global dataset="+dataset+"'"
    if options.debug:
        print(query)
    lsCmd = subprocess.Popen(['das_client.py '+query+' --limit=0'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    str_files, err = lsCmd.communicate()
    files.extend(['root://cms-xrd-global.cern.ch/'+ifile for ifile in str_files.split("\n")])
    files.pop()
    if options.debug:
        for ifile in files:
            print(ifile)

process.source = cms.Source("PoolSource",
#                             inputCommands = cms.untracked.vstring(
#                                 'keep *',
#                                 'drop *_hltEcalDigis_*_*',
#                                 'drop *_hltTriggerSummaryAOD_*_*'
#                             ),
                            #fileNames = cms.untracked.vstring(files)
                            fileNames = cms.untracked.vstring(
"/store/mc/RunIISummer17PrePremix/Neutrino_E-10_gun/GEN-SIM-DIGI-RAW/MCv2_correctPU_94X_mc2017_realistic_v9-v1/20025/FED79EE9-CE15-E811-99DF-0CC47A7C3572.root"
#"/store/data/Run2018A/SingleMuon/AOD/PromptReco-v3/000/316/717/00000/CA9D04B2-EB65-E811-8372-FA163EE0B5EF.root"
#                                "/store/data/Run2018A/AlCaPhiSym/RAW/v1/000/315/259/00000/9443EE2A-9F49-E811-9B86-FA163E7121A5.root"
#                               "/store/data/Run2018A/SingleMuon/RAW/v1/000/316/944/00000/18C27320-8865-E811-91E4-FA163EF26B2D.root"
                            )
                        )

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.20 $'),
    annotation = cms.untracked.string('step_PHISYM nevts:'+str(options.maxEvents)),
    name = cms.untracked.string('PhiSymProducer')
)

#ecalMultiFitUncalibRecHit
process.ecalMultiFitUncalibRecHit.EBdigiCollection = cms.InputTag("hltEcalPhiSymFilter","phiSymEcalDigisEB")
process.ecalMultiFitUncalibRecHit.EEdigiCollection = cms.InputTag("hltEcalPhiSymFilter","phiSymEcalDigisEE")

#ecalRecHit (no ricovery)
process.ecalRecHit.killDeadChannels = cms.bool( False )
process.ecalRecHit.recoverEBVFE = cms.bool( False )
process.ecalRecHit.recoverEEVFE = cms.bool( False )
process.ecalRecHit.recoverEBFE = cms.bool( False )
process.ecalRecHit.recoverEEFE = cms.bool( False )
process.ecalRecHit.recoverEEIsolatedChannels = cms.bool( False )
process.ecalRecHit.recoverEBIsolatedChannels = cms.bool( False )

# PHISYM producer
process.load('PhiSym.EcalCalibAlgos.PhiSymProducer_cfi')
# process.PhiSymProducer.makeSpectraTreeEB = True
# process.PhiSymProducer.makeSpectraTreeEE = True
#process.PhiSymProducer.eThreshold_barrel = 1.1
#process.PhiSymProducer.barrelHitCollection = cms.InputTag('ecalRecHit', 'EcalRecHitsEB', 'PHISYM')
#EBDigiCollection                      "selectDigi"                "selectedEcalEBDigiCollection"   "RECO"    



# Output definition
PHISYM_output_commands = cms.untracked.vstring(
    "drop *",
    "keep *_PhiSymProducer_*_*")

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
                                         splitLevel = cms.untracked.int32(0),
                                         outputCommands = PHISYM_output_commands,
                                         fileName = cms.untracked.string(options.outputFile)
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("phisym_spectra.root"))

# GLOBAL-TAG
from CondCore.DBCommon.CondDBSetup_cfi import *
process.GlobalTag = cms.ESSource("PoolDBESSource",
                                 CondDBSetup,
                                 connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
#                                 globaltag = cms.string('101X_dataRun2_Prompt_v9'),             #RunD
#                                 globaltag = cms.string('102X_dataRun2_Sep2018Rereco_v1'), #RunABC
                                 globaltag = cms.string('80X_dataRun2_2016LegacyRepro_v4'), #2016
                                 # Get individual tags (template)
                                 toGet = cms.VPSet(
#                                     cms.PSet(record = cms.string("EcalIntercalibConstantsRcd"),
#                                              tag = cms.string("EcalIntercalibConstants_2017_2015_at_high_eta"),
#                                              connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
#                                          ),
                                     cms.PSet(record = cms.string("EcalLaserAPDPNRatiosRcd"),
                                              tag = cms.string("EcalLaserAPDPNRatios_rereco2016_v1"),
                                              connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
                                          ),
                                     cms.PSet(record = cms.string("EcalPedestalsRcd"),
                                              tag = cms.string("EcalPedestals_timestamp_UltraLegacy_2016_v1"),
                                              connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
                                          ),
                                     cms.PSet(record = cms.string("EcalPulseShapesRcd"),
                                              tag = cms.string("EcalPulseShapes_Ultimate2016_calib"),
                                              connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
                                          ),
                                 )
)

# SCHEDULE

#process.hcalDigis.silent = cms.untracked.bool(False)
process.hcalDigis.InputLabel = rawTag
process.ecalDigis = process.ecalEBunpacker.clone()
process.ecalDigis.InputLabel = rawTag
#process.hbheprerecogpu.processQIE11 = cms.bool(True)

process.finalize = cms.EndPath(process.Out)

process.bunchSpacing = cms.Path(
    process.bunchSpacingProducer
)

process.digiPath = cms.Path(
    #process.hcalDigis
    process.ecalDigis
)



process.reconstruction_step = cms.Sequence( process.bunchSpacingProducer * (process.ecalMultiFitUncalibRecHit + process.ecalRecHit) )


process.schedule = cms.Schedule(
    process.bunchSpacing,
    process.digiPath,
    process.reconstruction_step,
#    process.ecalecalLocalRecoSequence,
    #process.TreeProducer_step,
)

process.p *= process.offlineBeamSpot
process.p *= process.PhiSymProducer

process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)
process.schedule = cms.Schedule(process.p, process.RECOSIMoutput_step)








import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.AlCa.GlobalTag import GlobalTag

# parse commad line options
options = VarParsing('analysis')
options.maxEvents = -1
options.outputFile = 'phisym_weights_1lumis.root'
options.parseArguments()

process=cms.Process("PHISYM")

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.Geometry.GeometryExtended2015Reco_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('RecoLocalCalo.EcalRecProducers.ecalMultiFitUncalibRecHit_cfi')
process.load('RecoLocalCalo.EcalRecProducers.ecalUncalibRecHit_cfi')
process.load('RecoLocalCalo.EcalRecProducers.ecalRecHit_cfi')
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cff")

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
process.source = cms.Source("PoolSource",
                            inputCommands = cms.untracked.vstring(
                                'keep *',
                                'drop *_hltEcalDigis_*_*',
                                'drop *_hltTriggerSummaryAOD_*_*'
                            ),
                            fileNames = cms.untracked.vstring(
                                "/store/data/Run2015D/AlCaPhiSym/RAW/v1/000/256/587/00000/C4E4637E-FF5B-E511-B847-02163E013735.root"
                                #"root://cmsxrootd-site.fnal.gov//store/data/Run2015B/AlCaPhiSym/RAW/v1/000/251/562/00000/0014158C-7728-E511-8847-02163E0122C2.root",
))

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.20 $'),
    annotation = cms.untracked.string('step_PHISYM nevts:'+str(options.maxEvents)),
    name = cms.untracked.string('PhiSymProducer')
)

isStream=True
runMultiFit=True
isBX50ns=False

if (runMultiFit):
    if (isBX50ns):
        process.ecalMultiFitUncalibRecHit.algoPSet = cms.PSet(
            useLumiInfoRunHeader = cms.bool(False),
            activeBXs = cms.vint32(-4,-2,0,2,4)
        )
    else:
        process.ecalMultiFitUncalibRecHit.algoPSet = cms.PSet(
            useLumiInfoRunHeader = cms.bool(False),
            activeBXs = cms.vint32(-5,-4,-3,-2,-1,0,1,2,3,4)
        ) 

#ecalUncalibRecHit
if (isStream):
    process.ecalUncalibRecHit.EBdigiCollection = cms.InputTag("hltEcalPhiSymFilter","phiSymEcalDigisEB")
    process.ecalUncalibRecHit.EEdigiCollection = cms.InputTag("hltEcalPhiSymFilter","phiSymEcalDigisEE")

#ecalMultiFitUncalibRecHit
if (isStream):
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
if (not runMultiFit):
    process.ecalRecHit.EBuncalibRecHitCollection = cms.InputTag("ecalUncalibRecHit","EcalUncalibRecHitsEB")
    process.ecalRecHit.EEuncalibRecHitCollection = cms.InputTag("ecalUncalibRecHit","EcalUncalibRecHitsEE")

# PHISYM producer
process.load('PhiSym.EcalCalibAlgos.PhiSymProducer_cfi')
#process.PhiSymProducer.makeSpectraTreeEB = True
#process.PhiSymProducer.makeSpectraTreeEE = True

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
                                 globaltag = cms.string('74X_dataRun2_Prompt_v2'),
                                 toGet = cms.VPSet(
                                     cms.PSet(record = cms.string("EcalIntercalibConstantsRcd"),
                                              tag = cms.string("EcalIntercalibConstants_2012ABCD_offline"),
                                              connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_ECAL"),
                                          ),
                                     cms.PSet(record = cms.string("EcalPulseShapesRcd"),
                                              tag = cms.string("EcalPulseShapes_data"),
                                              connect = cms.untracked.string("sqlite_file:ecaltemplates_popcon_weekly_best.db"),
                                          ),
                                     cms.PSet(record = cms.string("EcalPulseCovariancesRcd"),
                                              tag = cms.string("EcalPulseCovariances_data"),
                                              connect = cms.untracked.string("sqlite_file:ecalcovariances_popcon_weekly_best.db"),
                                          ),
                                     cms.PSet(record = cms.string("EBAlignmentRcd"),
                                              tag = cms.string("EBAlignment_measured_v10_offline"),
                                              connect = cms.untracked.string("frontier://FrontierProd/CMS_CONDITIONS"),
                                          ),
                                     cms.PSet(record = cms.string("EEAlignmentRcd"),
                                              tag = cms.string("EEAlignment_measured_v10_offline"),
                                              connect = cms.untracked.string("frontier://FrontierProd/CMS_CONDITIONS"),
                                          ),
                                     cms.PSet(record = cms.string("ESAlignmentRcd"), # only Bon!
                                              tag = cms.string("ESAlignment_measured_v08_offline"),
                                              connect = cms.untracked.string("frontier://FrontierProd/CMS_CONDITIONS"),
                                          )
                                 )
                            )

# SCHEDULE
if (not runMultiFit):
    process.reconstruction_step = cms.Sequence( process.ecalUncalibRecHit + process.ecalRecHit )
else:
    process.reconstruction_step = cms.Sequence( process.ecalMultiFitUncalibRecHit + process.ecalRecHit )

if (isStream):
    process.p = cms.Path(process.reconstruction_step)
    process.p *= process.offlineBeamSpot
    process.p *= process.PhiSymProducer
else:
    process.p = cms.Path(process.RawToDigi) 
    process.p *= process.L1Reco
    process.p *= process.reconstruction_step
    process.p *= process.offlineBeamSpot
    process.p *= process.PhiSymProducer

process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)
process.schedule = cms.Schedule(process.p, process.RECOSIMoutput_step)

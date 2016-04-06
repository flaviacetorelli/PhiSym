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
    input = cms.untracked.int32(50000)
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
                                #"root://cmsxrootd-site.fnal.gov//store/data/Run2015A/AlCaPhiSym/RAW/v1/000/247/720/00000/4C0AF78B-4810-E511-8C09-02163E0143CB.root",
                                #"root://cmsxrootd-site.fnal.gov//store/data/Run2015B/AlCaPhiSym/RAW/v1/000/251/562/00000/0014158C-7728-E511-8847-02163E0122C2.root"

                                # "root://cmsxrootd-site.fnal.gov//store/data/Run2015B/AlCaPhiSym/RAW/v1/000/251/561/00000/141BB8BB-8928-E511-A6C6-02163E011B19.root",
                                # "root://cmsxrootd-site.fnal.gov//store/data/Run2015B/AlCaPhiSym/RAW/v1/000/251/561/00000/143169BD-8928-E511-9223-02163E0126E1.root",
                                # "root://cmsxrootd-site.fnal.gov//store/data/Run2015B/AlCaPhiSym/RAW/v1/000/251/561/00000/1A616841-7128-E511-837E-02163E0138B3.root",
                                # "root://cmsxrootd-site.fnal.gov//store/data/Run2015B/AlCaPhiSym/RAW/v1/000/251/561/00000/226B28BB-8928-E511-A130-02163E0133D1.root",
                                # "root://cmsxrootd-site.fnal.gov//store/data/Run2015B/AlCaPhiSym/RAW/v1/000/251/561/00000/26DE7ABA-8928-E511-9718-02163E01412F.root",
                                # "root://cmsxrootd-site.fnal.gov//store/data/Run2015B/AlCaPhiSym/RAW/v1/000/251/561/00000/28520B3B-7128-E511-B938-02163E013674.root",
                                # "root://cmsxrootd-site.fnal.gov//store/data/Run2015B/AlCaPhiSym/RAW/v1/000/251/561/00000/36EE68BB-8928-E511-A60E-02163E013901.root",
                                # "root://cmsxrootd-site.fnal.gov//store/data/Run2015B/AlCaPhiSym/RAW/v1/000/251/561/00000/3AA78CBC-8928-E511-B690-02163E014531.root",
                                # "root://cmsxrootd-site.fnal.gov//store/data/Run2015B/AlCaPhiSym/RAW/v1/000/251/561/00000/40073738-7128-E511-B7DF-02163E011A46.root",
                                # "root://cmsxrootd-site.fnal.gov//store/data/Run2015B/AlCaPhiSym/RAW/v1/000/251/561/00000/4CBE81BC-8928-E511-AB82-02163E012965.root",
                                # "root://cmsxrootd-site.fnal.gov//store/data/Run2015B/AlCaPhiSym/RAW/v1/000/251/561/00000/52044E3A-7128-E511-93B9-02163E012BD2.root",
                                # "root://cmsxrootd-site.fnal.gov//store/data/Run2015B/AlCaPhiSym/RAW/v1/000/251/561/00000/541653BD-8928-E511-B6FD-02163E012283.root",
                                # "root://cmsxrootd-site.fnal.gov//store/data/Run2015B/AlCaPhiSym/RAW/v1/000/251/561/00000/5EABF842-6528-E511-8BDF-02163E011D88.root",
                                # "root://cmsxrootd-site.fnal.gov//store/data/Run2015B/AlCaPhiSym/RAW/v1/000/251/561/00000/64DEFF48-6528-E511-8CF8-02163E011DCE.root",
                                # "root://cmsxrootd-site.fnal.gov//store/data/Run2015B/AlCaPhiSym/RAW/v1/000/251/561/00000/686420CA-8928-E511-98E4-02163E011DDE.root",
                                # "root://cmsxrootd-site.fnal.gov//store/data/Run2015B/AlCaPhiSym/RAW/v1/000/251/561/00000/68EF982E-4F28-E511-94C1-02163E0133F2.root",
                                # "root://cmsxrootd-site.fnal.gov//store/data/Run2015B/AlCaPhiSym/RAW/v1/000/251/561/00000/6A593EB5-8928-E511-875E-02163E01192D.root",
                                # "root://cmsxrootd-site.fnal.gov//store/data/Run2015B/AlCaPhiSym/RAW/v1/000/251/561/00000/76A185BA-8928-E511-A817-02163E01416E.root",
                                # "root://cmsxrootd-site.fnal.gov//store/data/Run2015B/AlCaPhiSym/RAW/v1/000/251/561/00000/7CAA69BD-8928-E511-9438-02163E013406.root",
                                # "root://cmsxrootd-site.fnal.gov//store/data/Run2015B/AlCaPhiSym/RAW/v1/000/251/561/00000/96D53EBC-8928-E511-B0B0-02163E0134FD.root",
                                # "root://cmsxrootd-site.fnal.gov//store/data/Run2015B/AlCaPhiSym/RAW/v1/000/251/561/00000/9CF4BE3E-7128-E511-98C0-02163E0127D3.root",
                                # "root://cmsxrootd-site.fnal.gov//store/data/Run2015B/AlCaPhiSym/RAW/v1/000/251/561/00000/A2D7D6BB-8928-E511-935D-02163E0141EF.root",
                                # "root://cmsxrootd-site.fnal.gov//store/data/Run2015B/AlCaPhiSym/RAW/v1/000/251/561/00000/AAC052D1-8928-E511-AA97-02163E014629.root",
                                # "root://cmsxrootd-site.fnal.gov//store/data/Run2015B/AlCaPhiSym/RAW/v1/000/251/561/00000/AC2028BD-8928-E511-9EFE-02163E012AA4.root",
                                # "root://cmsxrootd-site.fnal.gov//store/data/Run2015B/AlCaPhiSym/RAW/v1/000/251/561/00000/BE0CB3BB-8928-E511-B062-02163E012601.root",
                                # "root://cmsxrootd-site.fnal.gov//store/data/Run2015B/AlCaPhiSym/RAW/v1/000/251/561/00000/CAB0F4BA-8928-E511-9772-02163E0117FF.root",
                                # "root://cmsxrootd-site.fnal.gov//store/data/Run2015B/AlCaPhiSym/RAW/v1/000/251/561/00000/CC009E3C-7128-E511-911C-02163E011D23.root",
                                # "root://cmsxrootd-site.fnal.gov//store/data/Run2015B/AlCaPhiSym/RAW/v1/000/251/561/00000/FC9B1038-6528-E511-BB04-02163E0139A2.root",
                                #"root://cmsxrootd-site.fnal.gov//store/data/Run2015B/AlCaPhiSym/RAW/v1/000/251/561/00000/FE3699BB-8928-E511-B0B2-02163E011B42.root"
                                #2015D
                                "root://eoscms//eos/cms/store/data/Run2015D/AlCaPhiSym/RAW/v1/000/258/159/00000/0434665B-C969-E511-86CD-02163E014200.root"
                            )
)

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

process.PhiSymProducer.makeSpectraTreeEB = True
process.PhiSymProducer.makeSpectraTreeEE = True

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
#process.GlobalTag = GlobalTag(process.GlobalTag, 'GR_P_V56')
process.GlobalTag = GlobalTag(process.GlobalTag, '74X_dataRun2_Prompt_v0')
# process.GlobalTag.toGet = cms.VPSet(
#     cms.PSet(record = cms.string("EcalIntercalibConstantsRcd"),
#              tag = cms.string("EcalIntercalibConstants_2012ABCD_offline"),
#              connect = cms.string("frontier://PromptProd/CMS_COND_31X_ECAL")
#          )
# )

# L1 filter for Lone bunch studies
process.triggerSelectionLoneBunch = cms.EDFilter( "TriggerResultsFilter",
                                                  triggerConditions = cms.vstring('L1_AlwaysTrue'),
                                                  hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
                                                  l1tResults = cms.InputTag( "hltGtDigis" ),
                                                  l1tIgnoreMask = cms.bool( False ),
                                                  l1techIgnorePrescales = cms.bool( False ),
                                                  daqPartitions = cms.uint32( 1 ),
                                                  throw = cms.bool( True )
)

# SCHEDULE
if (not runMultiFit):
    process.reconstruction_step = cms.Sequence( process.ecalUncalibRecHit + process.ecalRecHit )
else:
    process.reconstruction_step = cms.Sequence( process.ecalMultiFitUncalibRecHit + process.ecalRecHit )

process.phisymSequence = cms.Sequence()
    
if (isStream):
    process.phisymSequence += process.reconstruction_step
    process.phisymSequence += process.offlineBeamSpot
    process.phisymSequence += process.PhiSymProducer
else:
    process.phisymSequence = cms.Path(process.RawToDigi) 
    process.phisymSequence *= process.L1Reco
    process.phisymSequence *= process.reconstruction_step
    process.phisymSequence *= process.offlineBeamSpot
    process.phisymSequence *= process.PhiSymProducer

process.path = cms.Path(process.triggerSelectionLoneBunch*process.phisymSequence)
process.RECOSIMoutput.SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('path'))

process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)
process.schedule = cms.Schedule(process.path, process.RECOSIMoutput_step)


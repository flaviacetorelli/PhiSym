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

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.Geometry.GeometryExtended2017Reco_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('RecoLuminosity.LumiProducer.bunchSpacingProducer_cfi')
process.load('RecoLocalCalo.EcalRecProducers.ecalMultiFitUncalibRecHit_cfi')
process.load('RecoLocalCalo.EcalRecProducers.ecalUncalibRecHit_cfi')
process.load('RecoLocalCalo.EcalRecProducers.ecalRecHit_cfi')
process.load('RecoVertex.BeamSpotProducer.BeamSpot_cff')

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
                                "root://eoscms//eos/cms/store/data/Run2016B/AlCaPhiSym/RAW/v1/000/272/798/00000/0040CAC7-8114-E611-B271-02163E0135C0.root"
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
process.PhiSymProducer.thrEEmod = 14.

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
                                 globaltag = cms.string('80X_dataRun2_2016LegacyRepro_Candidate_v2'),
                                 toGet = cms.VPSet(
                                     cms.PSet(record = cms.string("EcalPedestalsRcd"),
                                              tag = cms.string("EcalPedestals_timestamp_2016"),
                                              connect = cms.string("frontier://FrontierPrep/CMS_CONDITIONS"),
                                          ),
                                     cms.PSet(
                                         record = cms.string('EcalLaserAlphasRcd'),
                                         tag = cms.string('EcalLaserAlphas_EB_1.52Russian_1.5Chinese'),
                                         connect = cms.string('frontier://FrontierPrep/CMS_CONDITIONS')
                                     ),
                                     cms.PSet(record = cms.string("EcalLinearCorrectionsRcd"),
                                              tag = cms.string("EcalLinearCorrections_from2011_offline"),
                                              connect = cms.string("frontier://FrontierPrep/CMS_CONDITIONS"),
                                          ),
                                 #     globaltag = cms.string('92X_dataRun2_Prompt_v9'),
                                 # # Get individual tags (template)
                                 # toGet = cms.VPSet(
                                 #     cms.PSet(record = cms.string("EcalADCToGeVConstantRcd"),
                                 #              tag = cms.string("EcalADCToGeVConstant_plus_2.4prct_in_EE"),
                                 #              connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
                                 #          ),
                                 #     cms.PSet(record = cms.string("EcalIntercalibConstantsRcd"),
                                 #              tag = cms.string("EcalIntercalibConstants_2017_2015_at_high_eta"),
                                 #              connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
                                 #          ),
                                 #     # cms.PSet(record = cms.string("EcalPedestalsRcd"),
                                 #     #          tag = cms.string("EcalPedestals_Legacy2017_time_v1"),
                                 #     #          connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
                                 #     #      ),
                                 #     cms.PSet(
                                 #         record = cms.string('EcalLaserAPDPNRatiosRcd'),
                                 #         tag = cms.string('EcalLaserAPDPNRatios_offline_2017pp_v4_for_tests'),
                                 #         connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS')
                                 #     ),
                                 #     cms.PSet(record = cms.string("EcalPulseShapesRcd"),
                                 #              tag = cms.string("EcalPulseShapes_October2017_rereco_v1"),
                                 #              connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
                                 #     ),
                                 )
)

# SCHEDULE
process.reconstruction_step = cms.Sequence( process.bunchSpacingProducer * (process.ecalMultiFitUncalibRecHit + process.ecalRecHit) )

process.p = cms.Path(process.reconstruction_step)
process.p *= process.offlineBeamSpot
process.p *= process.PhiSymProducer

process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)
process.schedule = cms.Schedule(process.p, process.RECOSIMoutput_step)

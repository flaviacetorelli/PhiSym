import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.IOVBounds = cms.PSet(
    startingIOV     = cms.int32(0),
    nIOVs           = cms.int32(-1),
    manualSplitting = cms.bool(True),
    beginRuns       = cms.vint32(251022, 251023, 251024, 251025, 251026, 251027, 251028,
                                 251244, 251251, 251252, 251521, 251522, 251523, 251548,
                                 251559, 251560, 251561, 251562),
    endRuns         = cms.vint32(251022, 251023, 251024, 251025, 251026, 251027, 251028,
                                 251244, 251251, 251252, 251521, 251522, 251523, 251548,
                                 251559, 251560, 251561, 251562),
    IOVMaps         = cms.vstring(["$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/production2015/IOVmap.root"])
)

process.ioFilesOpt = cms.PSet(    
    outputFile = cms.string('summed_'),
    
    oldConstantsFiles = cms.vstring('/afs/cern.ch/work/s/spigazzi/ECAL/CMSSW_7_4_6_patch6/src/PhiSym/EcalCalibAlgos/data/EcalIntercalibConstants_2012D_newThr.dat'),
    #oldConstantsFiles = cms.vstring('/afs/cern.ch/work/s/spigazzi/ECAL/CMSSW_7_4_6_patch6/src/PhiSym/EcalCalibAlgos/data/giulia_ic.dat'),
    
    recoConstantsFile = cms.string('/afs/cern.ch/work/s/spigazzi/ECAL/CMSSW_7_4_6_patch6/src/PhiSym/EcalCalibAlgos/data/EcalIntercalibConstants_2012ABCD_25iov.dat'),
    #recoConstantsFile = cms.string('/afs/cern.ch/work/s/spigazzi/ECAL/CMSSW_7_4_12_patch4/src/PhiSym/EcalCalibAlgos/data/EcalIntercalibConstants_2015B_Abs.dat'),
    
    inputFiles = cms.vstring([        
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415_merged/phisym_merged_1.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415_merged/phisym_merged_16.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415_merged/phisym_merged_4.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415_merged/phisym_merged_10.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415_merged/phisym_merged_17.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415_merged/phisym_merged_5.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415_merged/phisym_merged_11.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415_merged/phisym_merged_18.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415_merged/phisym_merged_6.root",        
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415_merged/phisym_merged_12.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415_merged/phisym_merged_19.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415_merged/phisym_merged_7.root",        
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415_merged/phisym_merged_13.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415_merged/phisym_merged_2.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415_merged/phisym_merged_8.root",        
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415_merged/phisym_merged_14.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415_merged/phisym_merged_9.root",        
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415_merged/phisym_merged_15.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415_merged/phisym_merged_3.root"        

        "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415_m/phisym_intercalibs_1000blocks_1.root",
        "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415_m/phisym_intercalibs_1000blocks_16.root",
        "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415_m/phisym_intercalibs_1000blocks_4.root",
        "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415_m/phisym_intercalibs_1000blocks_10.root",
        "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415_m/phisym_intercalibs_1000blocks_17.root",
        "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415_m/phisym_intercalibs_1000blocks_5.root",
        "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415_m/phisym_intercalibs_1000blocks_11.root",
        "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415_m/phisym_intercalibs_1000blocks_18.root",
        "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415_m/phisym_intercalibs_1000blocks_6.root",        
        "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415_m/phisym_intercalibs_1000blocks_12.root",
        "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415_m/phisym_intercalibs_1000blocks_19.root",
        "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415_m/phisym_intercalibs_1000blocks_7.root",        
        "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415_m/phisym_intercalibs_1000blocks_13.root",
        "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415_m/phisym_intercalibs_1000blocks_2.root",
        "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415_m/phisym_intercalibs_1000blocks_8.root",        
        "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415_m/phisym_intercalibs_1000blocks_14.root",
        "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415_m/phisym_intercalibs_1000blocks_9.root",        
        "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415_m/phisym_intercalibs_1000blocks_15.root",
        "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415_m/phisym_intercalibs_1000blocks_3.root"        

        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415/phisym_intercalibs_1000blocks_1.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415/phisym_intercalibs_1000blocks_16.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415/phisym_intercalibs_1000blocks_4.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415/phisym_intercalibs_1000blocks_10.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415/phisym_intercalibs_1000blocks_17.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415/phisym_intercalibs_1000blocks_5.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415/phisym_intercalibs_1000blocks_11.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415/phisym_intercalibs_1000blocks_18.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415/phisym_intercalibs_1000blocks_6.root",        
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415/phisym_intercalibs_1000blocks_12.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415/phisym_intercalibs_1000blocks_19.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415/phisym_intercalibs_1000blocks_7.root",        
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415/phisym_intercalibs_1000blocks_13.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415/phisym_intercalibs_1000blocks_2.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415/phisym_intercalibs_1000blocks_8.root",        
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415/phisym_intercalibs_1000blocks_14.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415/phisym_intercalibs_1000blocks_9.root",        
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415/phisym_intercalibs_1000blocks_15.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415/phisym_intercalibs_1000blocks_3.root"        
    ])
)


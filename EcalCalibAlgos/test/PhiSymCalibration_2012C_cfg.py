import FWCore.ParameterSet.Config as cms

process = cms.PSet()

nIOVs=1

process.IOVBounds = cms.PSet(
    startingIOV = cms.int32(0),
    nIOVs       = cms.int32(nIOVs),    
    beginRuns   = cms.vint32(200781),
    beginLumis  = cms.vint32(-1),
    endRuns     = cms.vint32(200798),
    endLumis    = cms.vint32(-1)
)

process.ioFilesOpt = cms.PSet(    
    outputFile = cms.string('summed_'),
    
    oldConstantsFiles = cms.vstring(""),
        
    recoConstantsFile = cms.string('/afs/cern.ch/work/s/spigazzi/ECAL/CMSSW_7_4_6_patch6/src/PhiSym/EcalCalibAlgos/data/EcalIntercalibConstants_2012ABCD_17iov.dat'),
        
    inputFiles = cms.vstring([        
        "/afs/cern.ch/user/s/spigazzi/work/ECAL/CMSSW_7_4_12_patch4/src/PhiSym/EcalCalibAlgos/ntuples/2012C_v1/phisym_intercalibs_1000blocks_1.root",
        "/afs/cern.ch/user/s/spigazzi/work/ECAL/CMSSW_7_4_12_patch4/src/PhiSym/EcalCalibAlgos/ntuples/2012C_v1/phisym_intercalibs_1000blocks_10.root",
        "/afs/cern.ch/user/s/spigazzi/work/ECAL/CMSSW_7_4_12_patch4/src/PhiSym/EcalCalibAlgos/ntuples/2012C_v1/phisym_intercalibs_1000blocks_11.root",
        "/afs/cern.ch/user/s/spigazzi/work/ECAL/CMSSW_7_4_12_patch4/src/PhiSym/EcalCalibAlgos/ntuples/2012C_v1/phisym_intercalibs_1000blocks_12.root",
        "/afs/cern.ch/user/s/spigazzi/work/ECAL/CMSSW_7_4_12_patch4/src/PhiSym/EcalCalibAlgos/ntuples/2012C_v1/phisym_intercalibs_1000blocks_13.root",
        "/afs/cern.ch/user/s/spigazzi/work/ECAL/CMSSW_7_4_12_patch4/src/PhiSym/EcalCalibAlgos/ntuples/2012C_v1/phisym_intercalibs_1000blocks_14.root",
        "/afs/cern.ch/user/s/spigazzi/work/ECAL/CMSSW_7_4_12_patch4/src/PhiSym/EcalCalibAlgos/ntuples/2012C_v1/phisym_intercalibs_1000blocks_2.root",
        "/afs/cern.ch/user/s/spigazzi/work/ECAL/CMSSW_7_4_12_patch4/src/PhiSym/EcalCalibAlgos/ntuples/2012C_v1/phisym_intercalibs_1000blocks_3.root",
        "/afs/cern.ch/user/s/spigazzi/work/ECAL/CMSSW_7_4_12_patch4/src/PhiSym/EcalCalibAlgos/ntuples/2012C_v1/phisym_intercalibs_1000blocks_4.root",
        "/afs/cern.ch/user/s/spigazzi/work/ECAL/CMSSW_7_4_12_patch4/src/PhiSym/EcalCalibAlgos/ntuples/2012C_v1/phisym_intercalibs_1000blocks_5.root",
        "/afs/cern.ch/user/s/spigazzi/work/ECAL/CMSSW_7_4_12_patch4/src/PhiSym/EcalCalibAlgos/ntuples/2012C_v1/phisym_intercalibs_1000blocks_6.root",
        "/afs/cern.ch/user/s/spigazzi/work/ECAL/CMSSW_7_4_12_patch4/src/PhiSym/EcalCalibAlgos/ntuples/2012C_v1/phisym_intercalibs_1000blocks_7.root",
        "/afs/cern.ch/user/s/spigazzi/work/ECAL/CMSSW_7_4_12_patch4/src/PhiSym/EcalCalibAlgos/ntuples/2012C_v1/phisym_intercalibs_1000blocks_8.root",
        "/afs/cern.ch/user/s/spigazzi/work/ECAL/CMSSW_7_4_12_patch4/src/PhiSym/EcalCalibAlgos/ntuples/2012C_v1/phisym_intercalibs_1000blocks_9.root"
        ])
)

# for iov in range(0, 1):
#     process.ioFilesOpt.oldConstantsFiles.append('/afs/cern.ch/cms/CAF/CMSALCA/ALCA_ECALCALIB/energy-calibration-repository/phiSymmetry/2012D_LTCorrJan13/EcalIntercalibConstants_Oct13Jump_Corr_2012D_'+str(process.IOVBounds.beginRuns[iov])+'_'+str(process.IOVBounds.endRuns[iov])+'_Absolute.txt')
    #print process.ioFilesOpt.oldConstantsFiles[iov]
    #process.ioFilesOpt.oldConstantsFiles.append('/afs/cern.ch/cms/CAF/CMSALCA/ALCA_ECALCALIB/energy-calibration-repository/phiSymmetry/2012D_NewLaserTag/EcalIntercalibConstants_NLT_2012D_'+str(process.IOVBounds.beginRuns[iov])+'_'+str(process.IOVBounds.endRuns[iov])+'_NewCodeChecked_Absolute.txt')


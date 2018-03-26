#!/bin/sh

export X509_USER_PROXY=/afs/cern.ch/work/a/abeschi/PhySym/CMSSW_9_4_0/src/PhiSym/EcalCalibAlgos/ntuples/x509up_u81383 

cp /afs/cern.ch/work/a/abeschi/PhySym/CMSSW_9_4_0/src/PhiSym/EcalCalibAlgos/test/PhiSymCalibration_cfg.py . 
cd /afs/cern.ch/work/a/abeschi/PhySym/CMSSW_9_4_0 
eval `scramv1 runtime -sh` 
cd - 

PhiSymCalibration PhiSymCalibration_cfg.py eosdirs=crab_PHISYM-CMSSW_9_4_0-multifit-94X_dataRun2_ReReco_EOY17_v2-IC2017_RunB_v9_Merged_v2/180206_104913/0000/,crab_PHISYM-CMSSW_9_4_0-multifit-94X_dataRun2_ReReco_EOY17_v2-IC2017_RunC_v9_Merged_v2/180206_105059/0000/,crab_PHISYM-CMSSW_9_4_0-multifit-94X_dataRun2_ReReco_EOY17_v2-IC2017_RunD_v9_Merged_v2/180206_105134/0000/,crab_PHISYM-CMSSW_9_4_0-multifit-94X_dataRun2_ReReco_EOY17_v2-IC2017_RunE_v9_Merged_v2/180206_105212/0000/ iovbounds=1,999999 outputFile=test.root

cp *.root /afs/cern.ch/work/a/abeschi/PhySym/CMSSW_9_4_0/src/PhiSym/EcalCalibAlgos/test

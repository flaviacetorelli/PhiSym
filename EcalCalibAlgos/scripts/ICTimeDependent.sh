#!/bin/bash

filename="$1"
prefix="summed_"
suffixROOT=".root"
suffixTXT=".txt"

ICFolder="$CMSSW_BASE/src/PhiSym/TimeDependentIC/"                               #This should be passed as input parameter
correctionFolder="$CMSSW_BASE/src/PhiSym/CorrectionFiles/"                       #This should be passed as input parameter
plotFolder="/afs/cern.ch/user/a/abeschi/www/ECAL/PhiSym/IC2017UltraLegacy/"      #This should be passed as input parameter
doFullYearIC=false                                                                #This should be passed as input parameter
doTimeDependentIC=true                                                          #This should be passed as input parameter
eosDirs=/eos/cms/store/group/dpg_ecal/alca_ecalcalib/phiSymmetry/AlCaPhiSym/crab_PHISYM-CMSSW_9_4_0-multifit-94X_dataRun2_ReReco_EOY17_v6-Run2017B_UltraLegacy_v2_Merged_v1/190117_142021/0000/,/eos/cms/store/group/dpg_ecal/alca_ecalcalib/phiSymmetry/AlCaPhiSym/crab_PHISYM-CMSSW_9_4_0-multifit-94X_dataRun2_ReReco_EOY17_v6-Run2017C_UltraLegacy_v2_Merged_v1/190117_142001/0000/,/eos/cms/store/group/dpg_ecal/alca_ecalcalib/phiSymmetry/AlCaPhiSym/crab_PHISYM-CMSSW_9_4_0-multifit-94X_dataRun2_ReReco_EOY17_v6-Run2017D_UltraLegacy_v2_Merged_v1/190117_142139/0000/,/eos/cms/store/group/dpg_ecal/alca_ecalcalib/phiSymmetry/AlCaPhiSym/crab_PHISYM-CMSSW_9_4_0-multifit-94X_dataRun2_ReReco_EOY17_v6-Run2017E_UltraLegacy_v2_Merged_v1/190117_142256/0000/,/eos/cms/store/group/dpg_ecal/alca_ecalcalib/phiSymmetry/AlCaPhiSym/crab_PHISYM-CMSSW_9_4_0-multifit-94X_dataRun2_ReReco_EOY17_v6-Run2017F_UltraLegacy_v2_Merged_v1/190117_142427/0000/

if [ ! -d "$ICFolder" ]; then
  mkdir $ICFolder
fi

if [ ! -d "$correctionFolder" ]; then
  mkdir $correctionFolder
fi

if [ ! -d "$plotFolder" ]; then
  mkdir $plotFolder
  cp /afs/cern.ch/user/a/abeschi/www/index.php $plotFolder
fi

if $doFullYearIC ; then
	echo "Deriving IC for the whole year"
	PhiSymCalibration $CMSSW_BASE/src/PhiSym/EcalCalibAlgos/scripts/PhiSymCalibration_cfg.py eosdirs=$eosDirs iovbounds=1,999999 outputFile=IC2017.root
	MaterialCorrection $CMSSW_BASE/src/PhiSym/EcalCalibAlgos/test/MaterialCorrection_cfg.py
        python $CMSSW_BASE/src/PhiSym/EcalCalibAlgos/scripts/phisymROOT2TXT.py $CMSSW_BASE/src/PhiSym/EcalCalibAlgos/scripts/summed_1_999999.root -c $CMSSW_BASE/src/PhiSym/EcalCalibAlgos/scripts/corrections_1_999999.txt --rel > $CMSSW_BASE/src/PhiSym/EcalCalibAlgos/scripts/IC2017.txt
	#$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/scripts/MaterialCorrectionPlots $correctionFileROOT $plotFolderTmp
        mv $CMSSW_BASE/src/PhiSym/EcalCalibAlgos/scripts/IC2017.txt $correctionFolder
fi

if $doTimeDependentIC ; then
	echo "Time dependent IC are being calculated"
	while IFS='' read -r line || [[ -n "$line" ]]; do
		echo "Reading file: $line"
		sed -i "/.root/c \'$line\'" $CMSSW_BASE/src/PhiSym/EcalCalibAlgos/test/MaterialCorrectionIC_cfg.py
		MaterialCorrection $CMSSW_BASE/src/PhiSym/EcalCalibAlgos/test/MaterialCorrectionIC_cfg.py
		iov=${line##*/}
		iov=${iov#"$prefix"}
		iov=${iov%"$suffixROOT"}
		correctionFile="corrections_"
		correctionFileTXT="$correctionFile$iov$suffixTXT"
		correctionFileROOT="$correctionFile$iov$suffixROOT"
		ICFileName="IC_$iov$suffixTXT"
		python $CMSSW_BASE/src/PhiSym/EcalCalibAlgos/scripts/phisymROOT2TXT.py $line -c $correctionFileTXT --rel > $ICFileName
		mv $ICFileName $ICFolder
		mv $correctionFileTXT $correctionFolder
		plotFolderTmp="$plotFolder$iov/"
		if [ ! -d "$plotFolderTmp" ]; then
			mkdir $plotFolderTmp
			indexFile="$plotFolder/index.php"
			cp $indexFile $plotFolderTmp
		fi
	
	$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/scripts/MaterialCorrectionPlots $correctionFileROOT $plotFolderTmp
		mv $correctionFileROOT $correctionFolder
	done < "$1" 
fi


echo "Lufthansa, partner of Star Alliance, thanks you for choosign our company, we hope to see you again on board of our aircrafts"


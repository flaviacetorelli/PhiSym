import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.absoluteICs = cms.bool(False)
process.ioFilesOpt = cms.PSet(
    getRunsFromFileName = cms.bool(True),
    userOutputName = cms.bool(False),
    outputFiles = cms.vstring(""),
    outputFileBase = cms.string('corrections_'),
    inputFiles = cms.vstring(
        # # 2012C Boff
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012C_v1/summed_200781_200798.root",
        
        # # 2015A Boff
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015A_v2/summed_250896_250896.root"

        # 2012D Bon
        #"$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012D_newThr/summed_208538_208686.root",

        # 2015B Bon
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415/summed_251562_251562.root"

        # 2015D
        #"$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015D_7415/summed_256843_256843.root"

        # Comm2016
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/Comm2016_weights/summed_268930-1_269000-30.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/Comm2016_weights/summed_269012-1_269073-68.root"

        # Comm2017
        "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/test/ICRunF_noCorrection.root"
    )
)
        


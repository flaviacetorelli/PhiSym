import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.ioFilesOpt = cms.PSet(
    getRunsFromFileName = cms.bool(False),
    outputFileBase = cms.string('corrections_'),
    inputFiles = cms.vstring(
        # # 2012C Boff
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012C_v1/summed_200781_200798.root",
        
        # # 2015A Boff
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015A_v2/summed_250896_250896.root"

        # 2012D Bon
        "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012D_newThr/summed_208538_208686.root",

        # 2015B Bon
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415/summed_251562_251562.root"

        # 2015D
        "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015D_7415/summed_256843_256843.root"
        
    )
)
        


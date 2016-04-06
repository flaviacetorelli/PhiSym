import FWCore.ParameterSet.Config as cms

PhiSymMerger = cms.EDAnalyzer(
    "PhiSymMerger",
    infoTag = cms.untracked.InputTag("PhiSymProducer", "", "PHISYM"),
    recHitEBTag = cms.untracked.InputTag("PhiSymProducer", "EB", "PHISYM"),
    recHitEETag = cms.untracked.InputTag("PhiSymProducer", "EE", "PHISYM"),
    blocksToSum = cms.untracked.int32(100)
)

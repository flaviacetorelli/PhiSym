import FWCore.ParameterSet.Config as cms

PhiSymProducer = cms.EDProducer(
    "PhiSymProducer",
    barrelHitCollection = cms.InputTag('ecalRecHit', 'EcalRecHitsEB', 'PHISYM'),
    endcapHitCollection = cms.InputTag('ecalRecHit', 'EcalRecHitsEE', 'PHISYM'),
    beamspot = cms.InputTag('offlineBeamSpot'),
    eThreshold_barrel = cms.double(10*1.75*0.04),    
    etCut_barrel = cms.double(1), #this is actually summed to eThr in order to define the upper bound    
    etCut_endcap = cms.double(1), #this is actually summed to eThr in order to define the upper bound    
    A = cms.vdouble(4.87, 35.42), # B + A*ring 2016 thr are defined as two linear cut (one for iring<30 and one above)
    B = cms.vdouble(98.12, 787.79),
    thrEEmod = cms.double(1.),
    nMisCalib = cms.int32(10), # <= 10; even; central value does not count
    misCalibRangeEB = cms.vdouble(0.95, 1.05),
    misCalibRangeEE = cms.vdouble(0.90, 1.10),
    lumisToSum = cms.int32(1),          
    statusThreshold = cms.int32(0),
    makeSpectraTreeEB = cms.untracked.bool(False),
    makeSpectraTreeEE = cms.untracked.bool(False)
)

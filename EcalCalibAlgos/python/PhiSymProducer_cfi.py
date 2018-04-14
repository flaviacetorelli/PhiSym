import FWCore.ParameterSet.Config as cms

# Noise parametrization in EE:
#   if (iRing<31) (A*iRing)^2 + B
#   else           A*iRing    + B
# Thresholds are set in GeV
# theEEmod is a multiplier to change the trheshold by a constant factor, now is set to 1

PhiSymProducer = cms.EDProducer(
    "PhiSymProducer",
    barrelHitCollection = cms.InputTag('ecalRecHit', 'EcalRecHitsEB', 'PHISYM'),
    endcapHitCollection = cms.InputTag('ecalRecHit', 'EcalRecHitsEE', 'PHISYM'),
    beamspot = cms.InputTag('offlineBeamSpot'),
    eThresholds_barrel = cms.vdouble(
	1.29300, 1.27493, 1.26427, 1.27347, 1.26648, 1.25665, 1.25295, 1.24255, 1.25091, 1.24271,
	1.22783, 1.22645, 1.21825, 1.21550, 1.20936, 1.20175, 1.19644, 1.18338, 1.18541, 1.18404,
	1.19315, 1.18655, 1.17628, 1.17594, 1.16888, 1.16148, 1.15566, 1.14691, 1.14357, 1.14112,
	1.12880, 1.11827, 1.11261, 1.11111, 1.10938, 1.09581, 1.08796, 1.07868, 1.07913, 1.07543,
	1.06980, 1.06133, 1.05405, 1.05279, 1.04524, 1.04282, 1.03482, 1.02966, 1.02448, 1.01996,
	1.01693, 1.00512, 1.00165, 1.00270, 1.00317, 0.99470, 0.98861, 0.98488, 0.98534, 0.98458,
	0.98964, 0.98543, 0.98067, 0.98153, 0.97725, 0.97784, 0.97228, 0.96898, 0.96820, 0.96645,
	0.96555, 0.96404, 0.96314, 0.96746, 0.96637, 0.96283, 0.95876, 0.95732, 0.95394, 0.95482,
	0.94690, 0.94567, 0.94436, 0.94289, 0.94559, 0.94337, 0.94068, 0.93937, 0.94191, 0.94485,
	0.95557, 0.95403, 0.95205, 0.95691, 0.95625, 0.95814, 0.95904, 0.95931, 0.96337, 0.96601,
	0.96328, 0.96732, 0.96763, 0.97438, 0.97443, 0.97891, 0.97720, 0.97861, 0.98305, 0.99100,
	0.99281, 0.99258, 0.99011, 0.99744, 0.99849, 1.00693, 1.01120, 1.01004, 1.01722, 1.02164,
	1.02606, 1.03054, 1.03379, 1.04375, 1.04749, 1.05195, 1.05350, 1.05771, 1.06185, 1.06814,
	1.07002, 1.06818, 1.07013, 1.08236, 1.08526, 1.09418, 1.09726, 1.09779, 1.10780, 1.11819,
	1.12716, 1.13133, 1.13869, 1.14657, 1.14952, 1.16214, 1.16154, 1.16424, 1.17543, 1.18310,
	1.17199, 1.17836, 1.17830, 1.18707, 1.20023, 1.21305, 1.21402, 1.21616, 1.22748, 1.22796,
	1.23891, 1.24301, 1.24454, 1.25348, 1.25823, 1.26981, 1.26422, 1.26380, 1.26953, 1.30068
    ), 
    etCut_barrel = cms.double(1), #this is actually summed to eThr in order to define the upper bound    
    etCut_endcap = cms.double(1), #this is actually summed to eThr in order to define the upper bound    
    A = cms.vdouble(0,05, 2.85), 
    B = cms.vdouble(1.6, -87),
    thrEEmod = cms.double(1.),
    nMisCalib = cms.int32(10), # <= 10; even; central value does not count
    misCalibRangeEB = cms.vdouble(0.95, 1.05),
    misCalibRangeEE = cms.vdouble(0.90, 1.10),
    lumisToSum = cms.int32(1),          
    statusThreshold = cms.int32(0),
    makeSpectraTreeEB = cms.untracked.bool(False),
    makeSpectraTreeEE = cms.untracked.bool(False)
)

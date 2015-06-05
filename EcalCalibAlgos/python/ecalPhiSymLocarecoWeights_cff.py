import FWCore.ParameterSet.Config as cms

from RecoLocalCalo.EcalRecProducers.ecalMultiFitUncalibRecHit_cfi import *
from RecoLocalCalo.EcalRecProducers.ecalUncalibRecHit_cfi import *
from RecoLocalCalo.EcalRecProducers.ecalRecHit_cfi import *


#if (isBX50ns):
ecalMultiFitUncalibRecHit.algoPSet = cms.PSet(
    useLumiInfoRunHeader = cms.bool(False),
    activeBXs = cms.vint32(-4,-2,0,2,4)
    )

def custom_25ns(process):
    process.ecalMultiFitUncalibRecHit.algoPSet = cms.PSet(
        useLumiInfoRunHeader = cms.bool(False),
        activeBXs = cms.vint32(-5,-4,-3,-2,-1,0,1,2,3,4)
        ) 
    

ecalUncalibRecHit.EBdigiCollection = cms.InputTag("HLTEcalPhiSymFilter","phiSymEcalDigisEB")
ecalUncalibRecHit.EEdigiCollection = cms.InputTag("HLTEcalPhiSymFilter","phiSymEcalDigisEE")

ecalMultiFitUncalibRecHit.EBdigiCollection = cms.InputTag("HLTEcalPhiSymFilter","phiSymEcalDigisEB")
ecalMultiFitUncalibRecHit.EEdigiCollection = cms.InputTag("HLTEcalPhiSymFilter","phiSymEcalDigisEE")

ecalRecHit.killDeadChannels = cms.bool( False )
ecalRecHit.recoverEBVFE = cms.bool( False )
ecalRecHit.recoverEEVFE = cms.bool( False )
ecalRecHit.recoverEBFE = cms.bool( False )
ecalRecHit.recoverEEFE = cms.bool( False )
ecalRecHit.recoverEEIsolatedChannels = cms.bool( False )
ecalRecHit.recoverEBIsolatedChannels = cms.bool( False )

#if (not runMultiFit):
#    process.ecalRecHit.EBuncalibRecHitCollection = cms.InputTag("ecalUncalibRecHit","EcalUncalibRecHitsEB")
#    process.ecalRecHit.EEuncalibRecHitCollection = cms.InputTag("ecalUncalibRecHit","EcalUncalibRecHitsEE")




reconstruction_step_weights = cms.Sequence( process.ecalUncalibRecHit + process.ecalRecHit )
reconstruction_step_multiFit = cms.Sequence( process.ecalMultiFitUncalibRecHit + process.ecalRecHit )

#ifndef _PHISYM_PRODUCER_
#define _PHISYM_PRODUCER_

#include <iostream>
#include <fstream>
#include <memory>

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

#include "Calibration/Tools/interface/EcalRingCalibrationTools.h"
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbService.h"
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbRecord.h"

#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"

#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatusCode.h"

#include "PhiSym/EcalCalibDataFormats/interface/PhiSymRecHit.h"
#include "PhiSym/EcalCalibDataFormats/interface/PhiSymInfo.h"
#include "PhiSym/EcalCalibDataFormats/interface/PhiSymFile.h"

using namespace std;

//****************************************************************************************
		
class PhiSymProducer : public edm::one::EDProducer<edm::EndLuminosityBlockProducer,
						   edm::one::WatchLuminosityBlocks,
						   edm::Accumulator>
{
public:
    explicit PhiSymProducer(const edm::ParameterSet& pSet);
    ~PhiSymProducer();

private:

    //---methods
    virtual void beginJob();
    virtual void beginLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& setup);
    virtual void endLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& setup) {};
    virtual void endLuminosityBlockProduce(edm::LuminosityBlock& lumi, edm::EventSetup const& setup) override;
    virtual void accumulate(edm::Event const& event, edm::EventSetup const& setup) override;
    virtual void endJob();

    //---input 
    edm::ESHandle<EcalChannelStatus>     chStatus_;
    edm::Handle<EBRecHitCollection>      barrelRecHitsHandle_;
    edm::Handle<EERecHitCollection>      endcapRecHitsHandle_;
    edm::Handle<reco::BeamSpot>          recoBeamSpotHandle_;
    edm::EDGetTokenT<EBRecHitCollection> ebToken_;
    edm::EDGetTokenT<EBRecHitCollection> eeToken_;
    edm::EDGetTokenT<reco::BeamSpot>     beamSpotToken_;
    float          etCutEB_;
    vector<double> eThresholdsEB2016_;
    vector<double> eThresholdsEB_;
    float          etCutEE_;
    vector<double> A2016_;
    vector<double> B2016_;
    vector<double> A_;
    vector<double> B_;
    float          thrEEmod2016_;
    float          thrEEmod_;
    int            nMisCalib_;
    vector<double> misCalibRangeEB_;
    float          misCalibStepsEB_[11];    
    vector<double> misCalibRangeEE_;
    float          misCalibStepsEE_[11];
    int            lumisToSum_;
    int            statusThreshold_;
    int            nLumis_;
    //---geometry
    EcalRingCalibrationTools       calibRing_;
    static const short             kNRingsEB = EcalRingCalibrationTools::N_RING_BARREL;
    static const short             kNRingsEE = EcalRingCalibrationTools::N_RING_ENDCAP;
    static const short             ringsInOneEE = kNRingsEE/2;
    float                          etCutsEB_[kNRingsEB];
    float                          etCutsEE_[kNRingsEE];
    float                          eThresholdsEE_[kNRingsEE];

    //---output edm
    unique_ptr<PhiSymInfoCollection>   lumiInfo_;
    unique_ptr<PhiSymRecHitCollection> recHitCollEB_;
    unique_ptr<PhiSymRecHitCollection> recHitCollEE_;
    //---output bookeeping json
    ofstream outLSInfoJson_;
    //---output plain tree
    bool makeSpectraTreeEB_;
    bool makeSpectraTreeEE_;
    bool use2016Thresholds_;
    EBTree outEBTree_;
    EETree outEETree_;    
    edm::Service<TFileService> fs_;
};

//----------IMPLEMENTATION----------------------------------------------------------------

PhiSymProducer::PhiSymProducer(const edm::ParameterSet& pSet):
    ebToken_(consumes<EBRecHitCollection>(pSet.getParameter<edm::InputTag>("barrelHitCollection"))),
    eeToken_(consumes<EBRecHitCollection>(pSet.getParameter<edm::InputTag>("endcapHitCollection"))),
    beamSpotToken_(consumes<reco::BeamSpot>(pSet.getParameter<edm::InputTag>("beamspot"))),
    etCutEB_(pSet.getParameter<double>("etCut_barrel")),
    eThresholdsEB2016_(pSet.getParameter<vector<double> >("eThresholds_barrel2016")),
    eThresholdsEB_(pSet.getParameter<vector<double> >("eThresholds_barrel")),
    etCutEE_(pSet.getParameter<double>("etCut_endcap")),
    A2016_(pSet.getParameter<vector<double> >("A2016")),
    B2016_(pSet.getParameter<vector<double> >("B2016")),
    A_(pSet.getParameter<vector<double> >("A")),
    B_(pSet.getParameter<vector<double> >("B")),
    thrEEmod2016_(pSet.getParameter<double>("thrEEmod2016")),
    thrEEmod_(pSet.getParameter<double>("thrEEmod")),
    nMisCalib_(pSet.getParameter<int>("nMisCalib")/2),
    misCalibRangeEB_(pSet.getParameter<vector<double> >("misCalibRangeEB")),
    misCalibRangeEE_(pSet.getParameter<vector<double> >("misCalibRangeEE")),
    lumisToSum_(pSet.getParameter<int>("lumisToSum")),
    statusThreshold_(pSet.getParameter<int>("statusThreshold")),
    nLumis_(0),
    makeSpectraTreeEB_(pSet.getUntrackedParameter<bool>("makeSpectraTreeEB")),
    makeSpectraTreeEE_(pSet.getUntrackedParameter<bool>("makeSpectraTreeEE")),
    use2016Thresholds_(pSet.getUntrackedParameter<bool>("use2016Thresholds"))
{    
    //---register the product
    produces<PhiSymInfoCollection,   edm::Transition::EndLuminosityBlock>();
    produces<PhiSymRecHitCollection, edm::Transition::EndLuminosityBlock>("EB");
    produces<PhiSymRecHitCollection, edm::Transition::EndLuminosityBlock>("EE");

    //---create spectra output file
    if(makeSpectraTreeEB_)
    {
        outEBTree_ = EBTree("eb", "EB sprectra");
        outEBTree_.GetTTreePtr()->SetDirectory(fs_->file().GetDirectory(0));
    }
    if(makeSpectraTreeEE_)
    {
        outEETree_ = EETree("ee", "EE sprectra");
        outEETree_.GetTTreePtr()->SetDirectory(fs_->file().GetDirectory(0));        
    }
}

PhiSymProducer::~PhiSymProducer()
{}

void PhiSymProducer::beginJob()
{
    //---set E thresholds, Et cuts and miscalib steps
    //---spectrum window: E > thr && Et < cut
    //---NOTE: etCutsEE need the geometry, so it is set later in beginLumi
    for(int iRing=0; iRing<kNRingsEB; ++iRing)
        etCutsEB_[iRing] = -1;
    for(int iRing=0; iRing<ringsInOneEE; ++iRing)
    {
		if(use2016Thresholds_)
		{
		    if(iRing < 30)
		        eThresholdsEE_[iRing] = thrEEmod2016_*(B2016_[0] + A2016_[0]*iRing)/1000;
		    else
		        eThresholdsEE_[iRing] = thrEEmod2016_*(B2016_[1] + A2016_[1]*iRing)/1000;
		}

		else
		{
			if(iRing < 31)
			{	eThresholdsEE_[iRing] = thrEEmod_*(B_[0] + (A_[0]*iRing)*(A_[0]*iRing));
			}

			else
			{	eThresholdsEE_[iRing] = thrEEmod_*(B_[1] + A_[1]*iRing);
			}
		}

		eThresholdsEE_[iRing+ringsInOneEE] = eThresholdsEE_[iRing];
		etCutsEE_[iRing] = -1;
		etCutsEE_[iRing+ringsInOneEE] = -1;
    }
	
	for(unsigned int i=0; i<eThresholdsEB_.size(); ++i)
		cout << eThresholdsEB_[i] << endl;
	if(use2016Thresholds_)
	{	eThresholdsEB_.clear();
		eThresholdsEB_ = eThresholdsEB2016_;
	}
    //---misCalib value init (nMisCalib is half oj the correct value!)
    float misCalibStepEB = fabs(misCalibRangeEB_[1]-misCalibRangeEB_[0])/(nMisCalib_*2);
    float misCalibStepEE = fabs(misCalibRangeEE_[1]-misCalibRangeEE_[0])/(nMisCalib_*2);
    for(int iMis=-nMisCalib_; iMis<=nMisCalib_; ++iMis)
    {
        //--- 0 -> 0; -i -> [1...n/2]; +i -> [n/2+1...n]
        int index = iMis > 0 ? iMis+nMisCalib_ : iMis == 0 ? 0 : iMis+nMisCalib_+1;
        misCalibStepsEB_[index] = iMis*misCalibStepEB;
        misCalibStepsEE_[index] = iMis*misCalibStepEE;
    }

    //---reset spectra output file
    if(makeSpectraTreeEB_)
        outEBTree_.Reset();
    if(makeSpectraTreeEE_)            
        outEETree_.Reset();
}

void PhiSymProducer::endJob()
{
    //---close LS info json file
    outLSInfoJson_ << endl << "}" << endl;
    outLSInfoJson_.close();
}

void PhiSymProducer::beginLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& setup)
{
    //---update plain tree run and lumi info
    if(makeSpectraTreeEB_)
    {
        outEBTree_.run = lumi.luminosityBlockAuxiliary().run();
        outEBTree_.lumi = lumi.luminosityBlockAuxiliary().luminosityBlock();
    }
    if(makeSpectraTreeEE_)
    {
        outEETree_.run = lumi.luminosityBlockAuxiliary().run();
        outEETree_.lumi = lumi.luminosityBlockAuxiliary().luminosityBlock();
    }
    
    //---reset the RecHit and LumiInfo collection
    if(nLumis_ == 0)
    {
        lumiInfo_ = std::make_unique<PhiSymInfoCollection>();
        lumiInfo_->push_back(PhiSymInfo());
        lumiInfo_->back().SetStartLumi(lumi);
        recHitCollEB_ = std::make_unique<PhiSymRecHitCollection>();
        recHitCollEE_ = std::make_unique<PhiSymRecHitCollection>();
	//---get the ecal geometry
        edm::ESHandle<CaloGeometry> geoHandle;
        setup.get<CaloGeometryRecord>().get(geoHandle);
        const CaloSubdetectorGeometry* barrelGeometry =
            geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);;
        const CaloSubdetectorGeometry* endcapGeometry =
            geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalEndcap);;
        calibRing_.setCaloGeometry(&(*geoHandle));
	
	std::vector<DetId> barrelDetIds=barrelGeometry->getValidDetIds(DetId::Ecal, EcalBarrel);
	std::vector<DetId> endcapDetIds=endcapGeometry->getValidDetIds(DetId::Ecal, EcalEndcap);
	recHitCollEB_->resize(barrelDetIds.size());
	recHitCollEE_->resize(endcapDetIds.size());
	for(auto& ebDetId : barrelDetIds)
        {
	    EBDetId myId(ebDetId);
	    recHitCollEB_->at(myId.denseIndex())=PhiSymRecHit(ebDetId.rawId(), 0);
            int ring = calibRing_.getRingIndex(myId);
            //---set etCut if first lumi
            if(etCutsEB_[ring] == -1 && myId.iphi()==1)
            {
                auto cellGeometry = barrelGeometry->getGeometry(myId);
                float eta=cellGeometry->getPosition().eta();
                etCutsEB_[ring] = eThresholdsEB_[ring]/cosh(eta) + etCutEB_;
            }
        }
	for(auto& eeDetId : endcapDetIds)
        {
	    EEDetId myId(eeDetId);
            int ring = calibRing_.getRingIndex(myId) - kNRingsEB;
	    recHitCollEE_->at(myId.denseIndex())=PhiSymRecHit(eeDetId.rawId(), 0);
            //---set eCutEE if first lumi
            if(ring < ringsInOneEE && etCutsEE_[ring] == -1 && myId.ix() == 50)
            {
                auto cellGeometry = endcapGeometry->getGeometry(myId);
                etCutsEE_[ring] = eThresholdsEE_[ring]/cosh(cellGeometry->getPosition().eta()) + etCutEE_;
                etCutsEE_[ring+ringsInOneEE] = etCutsEE_[ring];
            }
        }
    }
    ++nLumis_;
}

void PhiSymProducer::endLuminosityBlockProduce(edm::LuminosityBlock& lumi, edm::EventSetup const& setup)    
{
    //---put the collection in the LuminosityBlocks tree
    if(nLumis_ == lumisToSum_)
    {
        lumiInfo_->back().SetEndLumi(lumi);

        //---dump LS information into json
        if(!outLSInfoJson_.is_open())
        {
            //---CRAB3 only transfers files ending in .root, so name this .root even though is a json
            outLSInfoJson_.open("phisym_lumi_info_json.root", ofstream::out);
            outLSInfoJson_ << "{" << endl;
        }
        else
            outLSInfoJson_ << "," << endl;

        outLSInfoJson_ << "\"" << lumi.luminosityBlockAuxiliary().beginTime().unixTime() << "\" : " 
                       << "["
                       << lumiInfo_->back().getStartLumi().run() << ","
                       << lumiInfo_->back().getStartLumi().luminosityBlock() << ","
                       << lumiInfo_->back().GetTotHitsEB() << ","
                       << lumiInfo_->back().GetNEvents() 
                       << "]";


	cout << lumiInfo_->size() << endl;
        
        lumi.put(std::move(lumiInfo_), "");
        lumi.put(std::move(recHitCollEB_), "EB");
        lumi.put(std::move(recHitCollEE_), "EE");

        nLumis_ = 0;
    }
}

void PhiSymProducer::accumulate(edm::Event const& event, edm::EventSetup const& setup)
{
    uint64_t totHitsEB=0;
    uint64_t totHitsEE=0;

    //---get recHits collections 
    event.getByToken(ebToken_, barrelRecHitsHandle_);  
    event.getByToken(eeToken_, endcapRecHitsHandle_);

    //---get the ecal geometry
    edm::ESHandle<CaloGeometry> geoHandle;
    setup.get<CaloGeometryRecord>().get(geoHandle);
    const CaloSubdetectorGeometry* barrelGeometry =
        geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);;
    const CaloSubdetectorGeometry* endcapGeometry =
        geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalEndcap);;
    setup.get<EcalChannelStatusRcd>().get(chStatus_);
    
    //---get the beamspot
    event.getByToken(beamSpotToken_, recoBeamSpotHandle_);
    
    //---get the laser corrections
    edm::Timestamp evtTimeStamp(event.time().value());
    edm::ESHandle<EcalLaserDbService> laser;
    setup.get<EcalLaserDbRecord>().get(laser);

    //---EB---
    for(auto& recHit : *barrelRecHitsHandle_.product())
    {
        float energy = recHit.energy();
        EBDetId ebHit = EBDetId(recHit.id());
        int ring = calibRing_.getRingIndex(ebHit);
        //---if recHit energy is below thr even with the highest miscalib skip this recHit
        if(energy*misCalibRangeEB_[1] < eThresholdsEB_[ring] && !makeSpectraTreeEB_)
            continue;
        float eta=barrelGeometry->getGeometry(ebHit)->getPosition().eta();
        
        //---check channel status
        if((*chStatus_)[ebHit].getStatusCode() > statusThreshold_)
            lumiInfo_->back().SetBadChannel(recHit.id(), (*chStatus_)[ebHit].getStatusCode());

        //---compute et + miscalibration
        float etValues[11];
        //---one can do this in one for loop from -nMis to +nMis but in this way the
        //---program is faster
        //---NOTE: nMisCalib is half on the value set in the cfg python
        etValues[0] = recHit.energy()/cosh(eta);
        for(int iMis=-nMisCalib_; iMis<0; ++iMis)
        {
            //--- 0 -> 0; -i -> [1...n/2]; +i -> [n/2+1...n]
            int index = iMis + nMisCalib_ + 1; 
            etValues[index] = etValues[0]*(1+misCalibStepsEB_[index]);
            //---set et to zero if out of range [e_thr, et_thr+1]
            if(etValues[index]*cosh(eta) < eThresholdsEB_[ring] || etValues[index] > etCutsEB_[ring])
                etValues[index] = 0;
        }
        for(int iMis=1; iMis<=nMisCalib_; ++iMis)
        {
            //--- 0 -> 0; -i -> [1...n/2]; +i -> [n/2+1...n]
            int index = iMis + nMisCalib_; 
            etValues[index] = etValues[0]*(1+misCalibStepsEB_[index]);
            //---set et to zero if out of range [e_thr, et_thr+1]
            if(etValues[index]*cosh(eta) < eThresholdsEB_[ring] || etValues[index] > etCutsEB_[ring])
                etValues[index] = 0;
        }
        //---set et to zero if out of range [e_thr, et_thr+1]
        if(energy < eThresholdsEB_[ring] || etValues[0] > etCutsEB_[ring])
            etValues[0] = 0;
        else
            ++totHitsEB;
        //---update the rechHit sumEt
        recHitCollEB_->at(ebHit.denseIndex()).AddHit(etValues,
        					     laser.product()->getLaserCorrection(recHit.id(), evtTimeStamp));
     
        //---fill the plain tree
        if(makeSpectraTreeEB_)
        {
            outEBTree_.ieta = ebHit.ieta();
            outEBTree_.iphi = ebHit.iphi();
            outEBTree_.et = recHit.energy()/cosh(eta);
            outEBTree_.GetTTreePtr()->Fill();
        }
    }

    //---EE---
    for(auto& recHit : *endcapRecHitsHandle_.product())
    {
        EEDetId eeHit = EEDetId(recHit.id());
        int ring = calibRing_.getRingIndex(eeHit) - kNRingsEB;
        float energy = recHit.energy();
        //---if recHit energy is below thr even with the highest miscalib skip this recHit
        if(energy*misCalibRangeEE_[1] < eThresholdsEE_[ring] && makeSpectraTreeEE_)
            continue;

        //---check channel status
        float eta=endcapGeometry->getGeometry(eeHit)->getPosition().eta();
        if((*chStatus_)[eeHit].getStatusCode() > statusThreshold_)
            lumiInfo_->back().SetBadChannel(recHit.id(), (*chStatus_)[eeHit].getStatusCode());

        //---compute et + miscalibration
        float etValues[11];
        //---one can do this in one for loop from -nMis to +nMis but in this way the
        //---program is faster
        //---NOTE: nMisCalib is half on the value set in the cfg python
        etValues[0] = recHit.energy()/cosh(eta);
        for(int iMis=-nMisCalib_; iMis<0; ++iMis)
        {
            //--- 0 -> 0; -i -> [1...n/2]; +i -> [n/2+1...n]
            int index = iMis + nMisCalib_ + 1; 
            etValues[index] = etValues[0]*(1+misCalibStepsEE_[index]);
            //---set et to zero if out of range [e_thr, et_thr+1]
            if(etValues[index]*cosh(eta) < eThresholdsEE_[ring] || etValues[index] > etCutsEE_[ring])
                etValues[index] = 0;
        }
        for(int iMis=1; iMis<=nMisCalib_; ++iMis)
        {
            //--- 0 -> 0; -i -> [1...n/2]; +i -> [n/2+1...n]
            int index = iMis + nMisCalib_;
            etValues[index] = etValues[0]*(1+misCalibStepsEE_[index]);
            //---set et to zero if out of range [e_thr, et_thr+1]
            if(etValues[index]*cosh(eta) < eThresholdsEE_[ring] || etValues[index] > etCutsEE_[ring])
                etValues[index] = 0;
        }
        //---set et to zero if out of range [e_thr, et_thr+1]
        if(energy < eThresholdsEE_[ring] || etValues[0] > etCutsEE_[ring])
            etValues[0] = 0;
        else
            ++totHitsEE;
        //---update the rechHit sumEt
        recHitCollEE_->at(eeHit.denseIndex()).AddHit(etValues,
        					     laser.product()->getLaserCorrection(recHit.id(), evtTimeStamp));
     
        //---fill the plain tree
        if(makeSpectraTreeEE_)
        {
            outEETree_.iring = ring;
            outEETree_.ix = eeHit.ix();
            outEETree_.iy = eeHit.iy();
            outEETree_.et = recHit.energy()/cosh(eta);
            outEETree_.GetTTreePtr()->Fill();
        }
    }

    //---update the beamspot    
    lumiInfo_->back().Update(recoBeamSpotHandle_.product(), totHitsEB, totHitsEE);
}

DEFINE_FWK_MODULE(PhiSymProducer);

#endif

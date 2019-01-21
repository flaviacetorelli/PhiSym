#ifndef _PHISYM_MERGER_
#define _PHISYM_MERGER_

/* ***************************************************************************************
 * kSides=2; 0->EB+, 1->EB-;
 *           0->EE+, 1->EE-;
 * **************************************************************************************/

#include <fstream>

#include "TGraphErrors.h"
#include "TF1.h"

#include "FWCore/Utilities/interface/BranchType.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Common/interface/Provenance.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "Calibration/Tools/interface/EcalRingCalibrationTools.h"

#include "DataFormats/Provenance/interface/ProductIDToBranchID.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "PhiSym/EcalCalibAlgos/interface/utils.h"
#include "PhiSym/EcalCalibDataFormats/interface/PhiSymRecHit.h"
#include "PhiSym/EcalCalibDataFormats/interface/PhiSymInfo.h"
#include "PhiSym/EcalCalibDataFormats/interface/CalibrationFile.h"

using namespace std;

class PhiSymMerger : public edm::EDAnalyzer
{
public:
    explicit PhiSymMerger(const edm::ParameterSet& pSet);
    ~PhiSymMerger() {};

    //---utils
    void         FillOutput();
    void         SearchLumiIOV(PhiSymRunLumi current);
    
    //---methods
    virtual void beginJob() override;
    virtual void beginLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& setup) override;
    virtual void endLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& setup) override;
    virtual void endRun(edm::Run const& run, edm::EventSetup const& setup) override;
    virtual void endJob() override;
    virtual void analyze(edm::Event const&, edm::EventSetup const&) override {};
    
private:
    //---inputs
    edm::Handle<PhiSymInfoCollection> infoHandle_;
    edm::Handle<PhiSymRecHitCollection> recHitEBHandle_;
    edm::Handle<PhiSymRecHitCollection> recHitEEHandle_;
    edm::EDGetTokenT<PhiSymInfoCollection> infoToken_;
    edm::EDGetTokenT<PhiSymRecHitCollection> recHitEBToken_;
    edm::EDGetTokenT<PhiSymRecHitCollection> recHitEEToken_;
    int    blocksToSum_;
    int    nSummedLumis_;
    int    nMisCalib_;
    int    nEvents_;
    int    nBlocks_;
    vector<float> misCalibValuesEB_;
    vector<float> misCalibValuesEE_;

    //---IOV splitting
    int                   IOV_;
    string                iovFile_;
    vector<PhiSymRunLumi> IOVBegins_, IOVEnds_;
    vector<double>        IOVTimes_;
    PhiSymRunLumi         nextIOVMark_;
    
    //---phisym objects
    EcalRingCalibrationTools calibRing_;
    static const short kNRingsEB = EcalRingCalibrationTools::N_RING_BARREL;
    static const short kNRingsEE = EcalRingCalibrationTools::N_RING_ENDCAP;
    PhiSymRecHit ebXstals_[EBDetId::kSizeForDenseIndexing];
    PhiSymRecHit ebXstalsEven_[EBDetId::kSizeForDenseIndexing];
    PhiSymRecHit ebXstalsOdd_[EBDetId::kSizeForDenseIndexing];
    int chStatusEB_[EBDetId::kSizeForDenseIndexing];
    PhiSymRecHit eeXstals_[EEDetId::kSizeForDenseIndexing];
    PhiSymRecHit eeXstalsEven_[EEDetId::kSizeForDenseIndexing];
    PhiSymRecHit eeXstalsOdd_[EEDetId::kSizeForDenseIndexing];
    int chStatusEE_[EEDetId::kSizeForDenseIndexing];

    //---outputs
    auto_ptr<CalibrationFile> outFile_;
    edm::Service<TFileService> fs_;
    bool firstRound_;
};

PhiSymMerger::PhiSymMerger(const edm::ParameterSet& pSet):    
    infoToken_(consumes<PhiSymInfoCollection, edm::BranchType::InLumi>(pSet.getUntrackedParameter<edm::InputTag>("infoTag"))),
    recHitEBToken_(consumes<PhiSymRecHitCollection, edm::BranchType::InLumi>(pSet.getUntrackedParameter<edm::InputTag>("recHitEBTag"))),
    recHitEEToken_(consumes<PhiSymRecHitCollection, edm::BranchType::InLumi>(pSet.getUntrackedParameter<edm::InputTag>("recHitEETag"))),
    blocksToSum_(pSet.getUntrackedParameter<int>("blocksToSum")),
    nSummedLumis_(1),
    nMisCalib_(-1),
    nEvents_(0),
    nBlocks_(0),
    IOV_(0),
    iovFile_(pSet.getUntrackedParameter<string>("IOVfile")),
    nextIOVMark_(PhiSymRunLumi()),
    firstRound_(true)
{}

void PhiSymMerger::beginJob()
{
    //---read list of IOV and add underflow and overflow IOVs if file exist
    //---otherwise make only one IOV (data will still be split in runs).
    TFile* file = TFile::Open(iovFile_.c_str(), "READ");
    if(file->IsOpen())
    {
        TTree* map = (TTree*)file->Get("iov_map");
        int firstRun, lastRun;
        int firstLumi, lastLumi;
        double avg_time;
        map->SetBranchAddress("firstRun", &firstRun);
        map->SetBranchAddress("lastRun", &lastRun);
        map->SetBranchAddress("firstLumi", &firstLumi);
        map->SetBranchAddress("lastLumi", &lastLumi);
        map->SetBranchAddress("unixTimeMean", &avg_time);
        for(int iEntry=0; iEntry<map->GetEntriesFast(); ++iEntry)
        {
            map->GetEntry(iEntry);
            IOVBegins_.push_back(PhiSymRunLumi(firstRun, firstLumi));
            IOVEnds_.push_back(PhiSymRunLumi(lastRun, lastLumi));
            IOVTimes_.push_back(avg_time);
        }
        file->Close();    
        IOVBegins_.insert(IOVBegins_.begin(), PhiSymRunLumi(0, 0));
        IOVEnds_.insert(IOVEnds_.begin(), PhiSymRunLumi(IOVBegins_[1].run, IOVBegins_[1].lumi-1));    
        IOVBegins_.push_back(PhiSymRunLumi(IOVEnds_.back().run, IOVEnds_.back().lumi+1));
        IOVEnds_.push_back(PhiSymRunLumi(10000000, 10000000));
        IOVTimes_.insert(IOVTimes_.begin(), 0);
        IOVTimes_.push_back(-1);
    }
    else
    {
        IOVBegins_.push_back(PhiSymRunLumi(0, 0));
        IOVEnds_.push_back(PhiSymRunLumi(10000000, 10000000));
        IOVTimes_.push_back(-1);
    }
    
    //---output file
    outFile_ = auto_ptr<CalibrationFile>(new CalibrationFile(&fs_->file()));
}

void PhiSymMerger::endJob()
{
    //---collect spare lumis
    if(nBlocks_!=0)
    {
        //---increment block output trees counters
        outFile_->eb_xstals.avg_time = IOVTimes_[IOV_];
        outFile_->eb_xstals.block++;
        outFile_->eb_xstals.n_lumis = nBlocks_*nSummedLumis_;
        outFile_->ee_xstals.block++;
        outFile_->ee_xstals.n_lumis = nBlocks_*nSummedLumis_;

        FillOutput();
    }
    //---finalize outputs
    outFile_->cd();
    outFile_->eb_xstals.Write("eb_xstals");
    outFile_->eb_xstals_even.Write("eb_even");
    outFile_->eb_xstals_odd.Write("eb_odd");
    outFile_->ee_xstals.Write("ee_xstals");
    outFile_->ee_xstals_even.Write("ee_even");
    outFile_->ee_xstals_odd.Write("ee_odd");
}

void PhiSymMerger::endRun(edm::Run const& run, edm::EventSetup const& setup)
{
    //---stop sums at run end
    //---increment block output trees counters
    outFile_->eb_xstals.avg_time = IOVTimes_[IOV_];
    outFile_->eb_xstals.block++;
    outFile_->eb_xstals.n_lumis = nBlocks_*nSummedLumis_;
    outFile_->ee_xstals.block++;
    outFile_->ee_xstals.n_lumis = nBlocks_*nSummedLumis_;

    FillOutput();
}

void PhiSymMerger::beginLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& setup)
{
    if(firstRound_)
    {
        //---get the correct IOV for the first lumi
        PhiSymRunLumi thisRunLumi(lumi.luminosityBlockAuxiliary().run(),
                                  lumi.luminosityBlockAuxiliary().luminosityBlock());
        SearchLumiIOV(thisRunLumi);
    }        
}

void PhiSymMerger::endLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& setup)
{
    //---get the ecal ring geometry
    edm::ESHandle<CaloGeometry> geoHandle;
    setup.get<CaloGeometryRecord>().get(geoHandle);
    calibRing_.setCaloGeometry(&(*geoHandle));
    
    //---get lumi and run number of this lumi
    PhiSymRunLumi thisRunLumi(lumi.luminosityBlockAuxiliary().run(),
                              lumi.luminosityBlockAuxiliary().luminosityBlock());
    
    //---call the output filler if current run do not belong to current IOV
    //---this effectevely sum lumi blocks up to the previous one
    if(thisRunLumi < IOVBegins_[IOV_] || thisRunLumi > IOVEnds_[IOV_])
    {
        //---increment output trees counters and set end run/lumi infos
        outFile_->eb_xstals.avg_time = IOVTimes_[IOV_];
        outFile_->eb_xstals.block++;
        outFile_->eb_xstals.n_lumis = nBlocks_*nSummedLumis_;
        outFile_->ee_xstals.block++;
        outFile_->ee_xstals.n_lumis = nBlocks_*nSummedLumis_;
        FillOutput();

        //---search the correct IOV for the current run
        SearchLumiIOV(thisRunLumi);
    }
    
    //---get PHISYM collections (skip void lumis)
    lumi.getByToken(infoToken_, infoHandle_);
    lumi.getByToken(recHitEBToken_, recHitEBHandle_);
    lumi.getByToken(recHitEEToken_, recHitEEHandle_);
    if(!infoHandle_.isValid())
        return;
    //---if good block count it and record first run
    if(nBlocks_ == 0)
    {
        outFile_->eb_xstals.begin[0] = thisRunLumi.run;
        outFile_->eb_xstals.begin[1] = thisRunLumi.lumi;
        outFile_->ee_xstals.begin[0] = thisRunLumi.run;
        outFile_->ee_xstals.begin[1] = thisRunLumi.lumi;
    }
    nEvents_ += infoHandle_.product()->back().GetNEvents();
    ++nBlocks_;

    //---get mis calib values from metadata
    if(firstRound_)
    {
        nSummedLumis_ = edm::parameterSet(*recHitEBHandle_.provenance()).getParameter<int>("lumisToSum");
        nMisCalib_ = edm::parameterSet(*recHitEBHandle_.provenance()).getParameter<int>("nMisCalib");
        vector<double> misCalibRangeEB =
            edm::parameterSet(*recHitEBHandle_.provenance()).getParameter<vector<double> >("misCalibRangeEB");
        vector<double> misCalibRangeEE =
            edm::parameterSet(*recHitEEHandle_.provenance()).getParameter<vector<double> >("misCalibRangeEE");

        float misCalibStepEB = fabs(misCalibRangeEB[1]-misCalibRangeEB[0])/nMisCalib_;
        float misCalibStepEE = fabs(misCalibRangeEE[1]-misCalibRangeEE[0])/nMisCalib_;
        misCalibValuesEB_.resize(nMisCalib_+1);
        misCalibValuesEE_.resize(nMisCalib_+1);
        for(int iMis=-nMisCalib_/2; iMis<=nMisCalib_/2; ++iMis)
        {
            //--- 0 -> 0; -i -> [1...n/2]; +i -> [n/2+1...n]
            int index = iMis > 0 ? iMis+nMisCalib_/2 : iMis == 0 ? 0 : iMis+nMisCalib_/2+1;
            misCalibValuesEB_[index] = 1+misCalibStepEB*iMis;
            misCalibValuesEE_[index] = 1+misCalibStepEE*iMis;
        }
        outFile_->StoreMisCalibs(misCalibValuesEB_, misCalibValuesEE_);

        //---do this only the first time
        firstRound_ = false;
    }

    const map<uint32_t, short>* badChMap = infoHandle_.product()->back().GetBadChannels();
    const map<uint32_t, short>* goodChMap = infoHandle_.product()->back().GetGoodChannels();

    //---EB---
    //---fill the rings Et sum
    for(auto& recHit : *recHitEBHandle_.product())
    {
        if(badChMap->find(recHit.GetRawId()) != badChMap->end())
            continue;
        EBDetId ebXstal(recHit.GetRawId());
        ebXstals_[ebXstal.denseIndex()] += recHit;
        chStatusEB_[ebXstal.denseIndex()] = (int)goodChMap->find(recHit.GetRawId())->second;
        if(thisRunLumi.lumi % 2 == 0)
            ebXstalsEven_[ebXstal.denseIndex()] += recHit;
        else
            ebXstalsOdd_[ebXstal.denseIndex()] += recHit;
    }
    //---EE---
    //---fill the rings Et sum
    for(auto& recHit : *recHitEEHandle_.product())
    {
        if(badChMap->find(recHit.GetRawId()) != badChMap->end())
            continue;

        EEDetId eeXstal(recHit.GetRawId());
        eeXstals_[eeXstal.denseIndex()] += recHit;
        chStatusEE_[eeXstal.denseIndex()] =  (int)goodChMap->find(recHit.GetRawId())->second;
        if(thisRunLumi.lumi % 2 == 0)
            eeXstalsEven_[eeXstal.denseIndex()] += recHit;
        else
            eeXstalsOdd_[eeXstal.denseIndex()] += recHit;
    }
    //---BeamSpot info---
    //---update current sum, each bs lumi value is weighted by the number of events.
    //---the means will be computed in FillOutput()
    outFile_->eb_xstals.mean_bs_x += infoHandle_.product()->back().GetSum('X');
    outFile_->eb_xstals.mean_bs_sigmax += infoHandle_.product()->back().GetSumSigma('X');
    outFile_->eb_xstals.mean_bs_y += infoHandle_.product()->back().GetSum('Y');
    outFile_->eb_xstals.mean_bs_sigmay += infoHandle_.product()->back().GetSumSigma('Y');
    outFile_->eb_xstals.mean_bs_z += infoHandle_.product()->back().GetSum('Z');
    outFile_->eb_xstals.mean_bs_sigmaz += infoHandle_.product()->back().GetSumSigma('Z');
    
    //---keep track of the current lumi as last lumi of the block
    outFile_->eb_xstals.end[0] = thisRunLumi.run;
    outFile_->eb_xstals.end[1] = thisRunLumi.lumi;
    outFile_->ee_xstals.end[0] = thisRunLumi.run;
    outFile_->ee_xstals.end[1] = thisRunLumi.lumi;
}

//---do not compute k-factors and IC, just add up lumis
void PhiSymMerger::FillOutput()
{
    //---compute beamspot position and spread means
    outFile_->eb_xstals.mean_bs_x /= nEvents_;
    outFile_->eb_xstals.mean_bs_sigmax /= nEvents_;
    outFile_->eb_xstals.mean_bs_y /= nEvents_;
    outFile_->eb_xstals.mean_bs_sigmay /= nEvents_;
    outFile_->eb_xstals.mean_bs_z /= nEvents_;
    outFile_->eb_xstals.mean_bs_sigmaz /= nEvents_;
    outFile_->ee_xstals.mean_bs_x = outFile_->eb_xstals.mean_bs_x;
    outFile_->ee_xstals.mean_bs_sigmax = outFile_->eb_xstals.mean_bs_sigmax;
    outFile_->ee_xstals.mean_bs_y = outFile_->eb_xstals.mean_bs_y;
    outFile_->ee_xstals.mean_bs_sigmay = outFile_->eb_xstals.mean_bs_sigmay;
    outFile_->ee_xstals.mean_bs_z = outFile_->eb_xstals.mean_bs_z;
    outFile_->ee_xstals.mean_bs_sigmaz = outFile_->eb_xstals.mean_bs_sigmaz;

    //---loop over the EB channels and store summed rechits
    for(uint32_t index=0; index<EBDetId::kSizeForDenseIndexing; ++index)
    {
        EBDetId ebXstal = EBDetId::detIdFromDenseIndex(index);
        //---fill the output trees
        //---global
        outFile_->eb_xstals.n_events = nEvents_;
        outFile_->eb_xstals.rec_hit = &ebXstals_[index];
        outFile_->eb_xstals.ieta = ebXstal.ieta();
        outFile_->eb_xstals.iphi = ebXstal.iphi();
        outFile_->eb_xstals.ch_status = chStatusEB_[ebXstal.denseIndex()];
        outFile_->eb_xstals.Fill();
        //---even lumis (or block of lumis)
        outFile_->eb_xstals_even.rec_hit = &ebXstalsEven_[index];
        outFile_->eb_xstals_even.ieta = ebXstal.ieta();
        outFile_->eb_xstals_even.iphi = ebXstal.iphi();
        outFile_->eb_xstals_even.ch_status = chStatusEB_[ebXstal.denseIndex()];
        outFile_->eb_xstals_even.Fill();
        //---odd lumis (or block of lumis)
        outFile_->eb_xstals_odd.rec_hit = &ebXstalsOdd_[index];
        outFile_->eb_xstals_odd.ieta = ebXstal.ieta();
        outFile_->eb_xstals_odd.iphi = ebXstal.iphi();
        outFile_->eb_xstals_odd.ch_status = chStatusEB_[ebXstal.denseIndex()];
        outFile_->eb_xstals_odd.Fill();
        
        //---reset channel status and sum
        ebXstals_[index].Reset();
        ebXstalsEven_[index].Reset();
        ebXstalsOdd_[index].Reset();
    }
    //---loop over the EE channels and store summed rechits
    for(uint32_t index=0; index<EEDetId::kSizeForDenseIndexing; ++index)
    {
        EEDetId eeXstal = EEDetId::detIdFromDenseIndex(index);
        int currentRing=calibRing_.getRingIndex(eeXstal)-kNRingsEB;            
        //---fill the output tree
        outFile_->ee_xstals.n_events = nEvents_;
        outFile_->ee_xstals.rec_hit = &eeXstals_[index];
        outFile_->ee_xstals.iring = currentRing<kNRingsEE/2 ? currentRing-kNRingsEE/2 : currentRing-kNRingsEE/2 + 1;
        outFile_->ee_xstals.ix = eeXstal.ix();
        outFile_->ee_xstals.iy = eeXstal.iy();
        outFile_->ee_xstals.ch_status = chStatusEE_[eeXstal.denseIndex()];
        outFile_->ee_xstals.Fill();            
        //---even lumis (or block of lumis)
        outFile_->ee_xstals_even.rec_hit = &eeXstalsEven_[index];
        outFile_->ee_xstals_even.iring = outFile_->ee_xstals.iring;
        outFile_->ee_xstals_even.ix = eeXstal.ix();
        outFile_->ee_xstals_even.iy = eeXstal.iy();
        outFile_->ee_xstals_even.ch_status = chStatusEE_[eeXstal.denseIndex()];
        outFile_->ee_xstals_even.Fill();
        //---odd lumis (or block of lumis)
        outFile_->ee_xstals_odd.rec_hit = &eeXstalsOdd_[index];
        outFile_->ee_xstals_odd.iring = outFile_->ee_xstals.iring;
        outFile_->ee_xstals_odd.ix = eeXstal.ix();
        outFile_->ee_xstals_odd.iy = eeXstal.iy();
        outFile_->ee_xstals_odd.ch_status = chStatusEE_[eeXstal.denseIndex()];
        outFile_->ee_xstals_odd.Fill();

        //---reset channel status and sum
        eeXstals_[index].Reset();
        eeXstalsEven_[index].Reset();
        eeXstalsOdd_[index].Reset();
    }
    
    //---reset global variables
    nEvents_=0;
    nBlocks_=0;
    outFile_->eb_xstals.mean_bs_x = 0;
    outFile_->eb_xstals.mean_bs_sigmax = 0;
    outFile_->eb_xstals.mean_bs_y = 0;
    outFile_->eb_xstals.mean_bs_sigmay = 0;
    outFile_->eb_xstals.mean_bs_z = 0;
    outFile_->eb_xstals.mean_bs_sigmaz = 0;
}

void PhiSymMerger::SearchLumiIOV(PhiSymRunLumi current)
{
    for(unsigned int iIOV=0; iIOV<IOVEnds_.size(); ++iIOV)
        if(current >= IOVBegins_[iIOV] && current <= IOVEnds_[iIOV])
        {
            IOV_ = iIOV;
            return;
        }

    IOV_ = 0;
    return;
}    

DEFINE_FWK_MODULE(PhiSymMerger);

#endif

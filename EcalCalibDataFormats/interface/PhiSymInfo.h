#ifndef DATAFORMAT_PHISYMINFO_H
#define DATAFORMAT_PHISYMINFO_H

#include <vector>
#include <map>

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Provenance/interface/LuminosityBlockID.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"

//****************************************************************************************
class PhiSymRunLumi
{
public:
    //---ctors---
    PhiSymRunLumi(){run=-1; lumi=-1;};
    PhiSymRunLumi(int r, int l){run=r; lumi=l;};
    
    //---dtor---
    ~PhiSymRunLumi(){};

    //---operators---
    friend bool operator==(const PhiSymRunLumi& lx, const PhiSymRunLumi& rx);
    friend bool operator!=(const PhiSymRunLumi& lx, const PhiSymRunLumi& rx){return !(lx==rx);};
    friend bool operator<(const PhiSymRunLumi& lx, const PhiSymRunLumi& rx);
    friend bool operator>(const PhiSymRunLumi& lx, const PhiSymRunLumi& rx){return rx<lx;};
    friend bool operator>=(const PhiSymRunLumi& lx, const PhiSymRunLumi& rx){return lx>rx || lx==rx;};
    friend bool operator<=(const PhiSymRunLumi& lx, const PhiSymRunLumi& rx){return lx<rx || lx==rx;};
    
    int run;
    int lumi;
};

//****************************************************************************************
class PhiSymInfo
{
public:
    //---ctors---
    PhiSymInfo();
    
    //---dtor---
    ~PhiSymInfo();

    //---setters---
    inline void     SetBadChannel(uint32_t id, short status) {badChannels_[id] = status;};
    inline void     SetGoodChannel(uint32_t id, short status) {goodChannels_[id] = status;};
    void            SetStartLumi(edm::LuminosityBlock const& lumi);
    void            SetEndLumi(edm::LuminosityBlock const& lumi);

    //---getters---
    inline uint64_t GetTotHitsEB()       const {return totHitsEB_;};
    inline uint64_t GetTotHitsEE()       const {return totHitsEE_;};
    inline uint32_t GetNEvents()         const {return nEvents_;};
    float           GetMean(char k)      const;
    float           GetMeanSigma(char k) const;
    float           GetSum(char k)       const;
    float           GetSumSigma(char k)  const;
    inline edm::LuminosityBlockID getStartLumi() const {return startLumi_;};
    inline edm::LuminosityBlockID getEndLumi() const {return endLumi_;};
    inline const std::map<uint32_t, short>* GetBadChannels() const {return &badChannels_;}
    inline const std::map<uint32_t, short>* GetGoodChannels() const {return &goodChannels_;}

    //---utils---
    void            Update(const reco::BeamSpot* bs, uint64_t& nEB, uint64_t& nEE);

    //---operators---
    friend std::ostream& operator<<(std::ostream& out, const PhiSymInfo& obj);

private:
    edm::LuminosityBlockID startLumi_;
    edm::LuminosityBlockID endLumi_;
    std::map<uint32_t, short>  badChannels_;
    std::map<uint32_t, short>  goodChannels_;
    
    uint64_t totHitsEB_;
    uint64_t totHitsEE_;
    uint32_t nEvents_;
    
    float    sumX_;
    float    sumSigmaX_;
    float    sumY_;
    float    sumSigmaY_;
    float    sumZ_;
    float    sumSigmaZ_;
};

typedef std::vector<PhiSymInfo> PhiSymInfoCollection;

#endif

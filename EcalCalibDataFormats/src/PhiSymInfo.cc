#include "PhiSym/EcalCalibDataFormats/interface/PhiSymInfo.h"

//**********PhiSymRunLumi*****************************************************************
bool operator==(const PhiSymRunLumi& lx, const PhiSymRunLumi& rx)
{
    return lx.run == rx.run && lx.lumi == rx.lumi;
}

bool operator<(const PhiSymRunLumi& lx, const PhiSymRunLumi& rx)
{
    return lx.run < rx.run || (lx.run == rx.run && lx.lumi < rx.lumi);
}

//**********PhiSymInfo********************************************************************
//**********constructors******************************************************************
PhiSymInfo::PhiSymInfo():
    totHitsEB_(0), totHitsEE_(0), nEvents_(0),
    sumX_(0), sumSigmaX_(0), sumY_(0), sumSigmaY_(0), sumZ_(0), sumSigmaZ_(0)
{}

//**********destructor********************************************************************
PhiSymInfo::~PhiSymInfo()
{}

//**********getters***********************************************************************
float PhiSymInfo::GetMean(char k) const
{
    if(k == 'X')
        return sumX_/nEvents_;
    if(k == 'Y')
        return sumY_/nEvents_;
    if(k == 'Z')
        return sumZ_/nEvents_;

    return -999;
}

float PhiSymInfo::GetMeanSigma(char k) const
{
    if(k == 'X')
        return sumSigmaX_/nEvents_;
    if(k == 'Y')
        return sumSigmaY_/nEvents_;
    if(k == 'Z')
        return sumSigmaZ_/nEvents_;

    return -999;
}

void PhiSymInfo::SetStartLumi(edm::LuminosityBlock const& lumi)
{
  startLumi_=lumi.luminosityBlockAuxiliary().id();
}

void PhiSymInfo::SetEndLumi(edm::LuminosityBlock const& lumi)
{
  endLumi_=lumi.luminosityBlockAuxiliary().id();
}

//**********utils*************************************************************************
void PhiSymInfo::Update(const reco::BeamSpot* bs, uint64_t& nEB, uint64_t& nEE)
{
    sumX_ += bs->x0();
    sumY_ += bs->y0();
    sumZ_ += bs->z0();

    sumSigmaX_ += bs->BeamWidthX();
    sumSigmaY_ += bs->BeamWidthY();
    sumSigmaZ_ += bs->sigmaZ();
    
    totHitsEB_ += nEB;
    totHitsEE_ += nEE;

    ++nEvents_;
}

//**********operators*********************************************************************

std::ostream& operator<<(std::ostream& out, const PhiSymInfo& obj)
{
    //---dump all the informations
    out << std::endl;
    out << std::setw(20) << "# event: " << obj.GetNEvents() << std::endl;
    out << std::setw(20) << "hits EB: " << obj.GetTotHitsEB() << std::endl;
    out << std::setw(20) << "hits EE: " << obj.GetTotHitsEE() << std::endl;
    out << std::setw(20) << "X: mean | meanSigma: " << obj.GetMean('X') << " | " << obj.GetMeanSigma('X') <<std::endl;
    out << std::setw(20) << "Y: mean | meanSigma: " << obj.GetMean('Y') << " | " << obj.GetMeanSigma('Y') <<std::endl;
    out << std::setw(20) << "Z: mean | meanSigma: " << obj.GetMean('Z') << " | " << obj.GetMeanSigma('Z') <<std::endl;
    out << std::endl;
    
    return out;
}

    



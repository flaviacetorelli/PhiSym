#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TMath.h"
#include "TRandom3.h"

#include <iostream>
#include <algorithm>
#ifdef __ROOTCINT__
#include <omp.h>
#endif

using namespace std;

class AlphaFitter
{
public:
    AlphaFitter() {};
    AlphaFitter(int n_points, const double* eflow_ratio, const double* eflow_ratio_err, const double* lc, float alpha_db)
        {
            n_points_ = n_points;
            eflow_ratio_ = eflow_ratio;
            eflow_ratio_err_ = eflow_ratio_err;
            lc_ = lc;
            alpha_db_ = alpha_db;
            best_fit_alpha_ = -999;
            best_fit_alpha_err_ = -1;
            best_fit_chi2_ = -1;
        };
    AlphaFitter(AlphaFitter& other)
        {
            alpha_db_ = other.alpha_db_;
            n_points_ = other.n_points_;
            eflow_ratio_ = other.eflow_ratio_;
            eflow_ratio_err_ = other.eflow_ratio_err_;
            lc_ = other.lc_;
            best_fit_alpha_ = other.best_fit_alpha_;
            best_fit_alpha_err_ = other.best_fit_alpha_err_;
            best_fit_scale_ = other.best_fit_scale_;
            best_fit_scale_err_ = other.best_fit_scale_err_;
            best_fit_chi2_ = other.best_fit_chi2_;
        }            
    AlphaFitter& operator=(AlphaFitter& other)
        {
            alpha_db_ = other.alpha_db_;
            n_points_ = other.n_points_;
            eflow_ratio_ = other.eflow_ratio_;
            eflow_ratio_err_ = other.eflow_ratio_err_;
            lc_ = other.lc_;
            best_fit_alpha_ = other.best_fit_alpha_;
            best_fit_alpha_err_ = other.best_fit_alpha_err_;
            best_fit_scale_ = other.best_fit_scale_;
            best_fit_scale_err_ = other.best_fit_scale_err_;
            best_fit_chi2_ = other.best_fit_chi2_;

            return *this;
        }            

    ~AlphaFitter() {};

    //---Getters---
    float GetAlpha()      {return best_fit_alpha_;};
    float GetAlphaError() {return best_fit_alpha_err_;};
    float GetScale()      {return best_fit_scale_;};
    float GetScaleError() {return best_fit_scale_err_;};
    float GetChi2()       {return best_fit_chi2_;};
    
    //---Main method
    float Fit()
        {
            ROOT::Math::Functor fChi2(this, &AlphaFitter::chi2, 2);
            ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");            
            min->SetMaxFunctionCalls(1000000);
            min->SetMaxIterations(100000);
            min->SetTolerance(0.0005);
            min->SetPrintLevel(0);

            min->SetFunction(fChi2);
            min->SetLimitedVariable(0, "delta_alpha", 0., 1e-3, -2., 2.);
            min->SetLimitedVariable(1, "scale", 1., 1e-3, 0.9, 1.1);
            
            min->Minimize();

            auto delta_alpha = min->X()[0];
            best_fit_chi2_ = min->MinValue()/(n_points_-2);
            best_fit_alpha_ = best_fit_chi2_>0 ? alpha_db_+delta_alpha : -1;
            best_fit_alpha_err_ = best_fit_chi2_>0 ? min->Errors()[0] : -1;
            best_fit_scale_ = min->X()[1];
            best_fit_scale_err_ = min->Errors()[1];

            delete min;
            return delta_alpha;
        };

protected:
    double chi2(const double* par)
        {
            double chi2=0;
#pragma omp parallel for reduction(+:chi2) 
            for(int point=0; point<n_points_; ++point)
            {
                if(lc_[point] > 0)
                    chi2 += pow((eflow_ratio_[point]*correction(lc_[point]/lc_[0], par[0])-par[1])/(eflow_ratio_err_[point]), 2);
            }
            
            return chi2;
        };

    double correction(double lc, const double& delta_alpha)
        {
            double transparencyLoss = TMath::Power(lc, 1./alpha_db_);
            if (transparencyLoss<0.01 || transparencyLoss>2. )
            {
                //std::cout << "***** RETURNING 1. TO AVOID MORE ISSUES" << std::endl;
                return 1.;
            }
            return TMath::Power(transparencyLoss, delta_alpha);
        };
        
private:
    float alpha_db_;
    int n_points_;
    const double* eflow_ratio_;
    const double* eflow_ratio_err_;
    const double* lc_;
    float best_fit_alpha_;
    float best_fit_alpha_err_;
    float best_fit_scale_;
    float best_fit_scale_err_;
    float best_fit_chi2_;
    
};
    

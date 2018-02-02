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

class AlphaAndPNFitter
{
public:
    AlphaAndPNFitter() {};
    AlphaAndPNFitter(int n_points,  vector<vector<double> >& eflow_ratios, vector<vector<double> >& eflow_ratios_err,
                     vector<vector<double> >& lc, vector<int>& avg_times, vector<float> alphas_db)
        {
            alphas_db_ = alphas_db;
            n_points_ = n_points;
            n_crystals_ = eflow_ratios.size();
            eflow_ratios_ = eflow_ratios;
            eflow_ratios_err_ = eflow_ratios_err;
            lc_ = lc;
            avg_times_ = avg_times;
            best_fit_alphas_ = vector<float>(n_crystals_, -999);
            best_fit_alphas_err_ = vector<float>(n_crystals_, -1);
            best_fit_scales_ = vector<float>(n_crystals_, -999);
            best_fit_scales_err_ = vector<float>(n_crystals_, -1);
            best_fit_pn_corr_ = vector<float>(n_points_, -999);
            best_fit_pn_corr_err_ = vector<float>(n_points_, -1);            
            best_fit_chi2_ = -1;
        };
    AlphaAndPNFitter(AlphaAndPNFitter& other)
        {
            alphas_db_ = other.alphas_db_;
            n_points_ = other.n_points_;
            n_crystals_ = other.n_crystals_;
            eflow_ratios_ = other.eflow_ratios_;
            eflow_ratios_err_ = other.eflow_ratios_err_;
            lc_ = other.lc_;
            avg_times_ = other.avg_times_;
            best_fit_alphas_ = other.best_fit_alphas_;
            best_fit_alphas_err_ = other.best_fit_alphas_err_;
            best_fit_scales_ = other.best_fit_scales_;
            best_fit_scales_err_ = other.best_fit_scales_err_;
            best_fit_pn_corr_ = other.best_fit_pn_corr_;
            best_fit_pn_corr_err_ = other.best_fit_pn_corr_err_;
            best_fit_chi2_ = other.best_fit_chi2_;
        }            
    AlphaAndPNFitter& operator=(AlphaAndPNFitter& other)
        {
            alphas_db_ = other.alphas_db_;
            n_points_ = other.n_points_;
            n_crystals_ = other.n_crystals_;            
            eflow_ratios_ = other.eflow_ratios_;
            eflow_ratios_err_ = other.eflow_ratios_err_;
            lc_ = other.lc_;
            avg_times_ = other.avg_times_;            
            best_fit_alphas_ = other.best_fit_alphas_;
            best_fit_alphas_err_ = other.best_fit_alphas_err_;
            best_fit_scales_ = other.best_fit_scales_;
            best_fit_scales_err_ = other.best_fit_scales_err_;
            best_fit_pn_corr_ = other.best_fit_pn_corr_;
            best_fit_pn_corr_err_ = other.best_fit_pn_corr_err_;
            best_fit_chi2_ = other.best_fit_chi2_;

            return *this;
        }            

    ~AlphaAndPNFitter() {};

    //---Getters---
    float GetAlphaDB(int ixstal)      {return alphas_db_[ixstal];};
    float GetAlpha(int ixstal)        {return best_fit_alphas_[ixstal];};
    float GetAlphaError(int ixstal)   {return best_fit_alphas_err_[ixstal];};
    float GetScale(int ixstal)        {return best_fit_scales_[ixstal];};
    float GetScaleError(int ixstal)   {return best_fit_scales_err_[ixstal];};
    float GetPNCorrection(int ixstal) {return best_fit_pn_corr_[ixstal];};
    float GetPNCorrError(int ixstal)  {return best_fit_pn_corr_err_[ixstal];};
    vector<float>& GetAlphasDB()      {return alphas_db_;};
    vector<float>& GetAlphas()        {return best_fit_alphas_;};
    vector<float>& GetAlphaErrors()   {return best_fit_alphas_err_;};
    vector<float>& GetScales()        {return best_fit_scales_;};
    vector<float>& GetScaleErrors()   {return best_fit_scales_err_;};
    vector<float>& GetPNCorrections() {return best_fit_pn_corr_;};
    vector<float>& GetPNCorrErrors()  {return best_fit_pn_corr_err_;};
    float GetChi2()                   {return best_fit_chi2_;};
    int   GetNCrystals()              {return n_crystals_;};
    
    //---Main method
    float Fit()
        {
            //ROOT::Math::Functor fChi2(this, &AlphaAndPNFitter::chi2, n_crystals_+n_points_);
            ROOT::Math::Functor fChi2(this, &AlphaAndPNFitter::chi2, n_crystals_*2+1);
            ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");            
            min->SetMaxFunctionCalls(1000000);
            min->SetMaxIterations(100000);
            min->SetTolerance(0.0005);
            min->SetPrintLevel(0);

            min->SetFunction(fChi2);
            for(int ixstal=0; ixstal<n_crystals_; ++ixstal)
            {
                string xstal_num = to_string(ixstal);
                min->SetLimitedVariable(ixstal*2, ("delta_alpha_"+xstal_num).c_str(), 0., 1e-3, -1.5, 2.);
                min->SetLimitedVariable(ixstal*2+1, ("scale_"+xstal_num).c_str(), 1., 1e-3, 0.95, 1.05);
            }
            // for(int point=n_crystals_+1; point<n_crystals_+n_points_; ++point)
            //     min->SetLimitedVariable(point, (string("pn_corr")+to_string(point)).c_str(), 0., 1e-3, -0.01, 0.01);
            min->SetLimitedVariable(n_crystals_*2, "pn_corr", 0., 1e-3, -0.1, 0.1);
            
            min->Minimize();

            best_fit_chi2_ = min->MinValue()/(n_points_*(n_crystals_-1)-n_crystals_);
            for(int ixstal=0; ixstal<n_crystals_; ++ixstal)
            {
                best_fit_alphas_[ixstal] = best_fit_chi2_>0 ? alphas_db_[ixstal]+min->X()[ixstal*2] : -1;
                best_fit_alphas_err_[ixstal] = best_fit_chi2_>0 ? min->Errors()[ixstal*2] : -1;
                best_fit_scales_[ixstal] = min->X()[ixstal*2+1];
                best_fit_scales_err_[ixstal] = min->Errors()[ixstal*2+1];
            }
            for(int point=0; point<1/*n_points_*/; ++point)
            {
                int index = n_crystals_*2+point;
                best_fit_pn_corr_[point] = min->X()[index];
                best_fit_pn_corr_err_[point] = min->Errors()[index];
            }

            delete min;
            return best_fit_chi2_;
        };

protected:
    double chi2(const double* par)
        {
            double chi2=0;
            for(int ixstal=0; ixstal<n_crystals_; ++ixstal)
            {
                auto eflow_ratio = eflow_ratios_[ixstal];
                auto eflow_ratio_err = eflow_ratios_err_[ixstal];
                auto lc = lc_[ixstal];
#pragma omp parallel for reduction(+:chi2) 
                for(int point=0; point<n_points_; ++point)
                {
                    if(lc[point] > 0)                        
                        chi2 += pow((eflow_ratio[point]*
                                     correction(ixstal, point, lc[point]/lc[0], par[n_crystals_*2], par[ixstal*2])-par[ixstal*2+1])/
                                    (eflow_ratio_err[point]), 2);
                }
                //---alpha_fit penalty
                chi2 += pow((par[ixstal*2]+alphas_db_[ixstal]-1.52)/0.25, 2);
            }
            //---PN correction penalty
            chi2 += pow(par[n_crystals_*2]/0.005, 2);
            
            return chi2;
        };

    double correction(int ixstal, int iov, double lc, const double pn_corr, const double& delta_alpha)
        {
            double transparencyLoss = TMath::Power(lc, 1./alphas_db_[ixstal]);
            if (transparencyLoss<0.01 || transparencyLoss>2. )
            {
                //std::cout << "***** RETURNING 1. TO AVOID MORE ISSUES" << std::endl;
                return 1.;
            }
            return (pow(transparencyLoss, delta_alpha)*
                    pow(1+pn_corr*double(avg_times_[iov]-avg_times_[0])/15e6, alphas_db_[ixstal]+delta_alpha));
        };
        
private:
    vector<float> alphas_db_;
    int n_points_;
    int n_crystals_;
    vector<vector<double> > eflow_ratios_;
    vector<vector<double> > eflow_ratios_err_;
    vector<vector<double> > lc_;
    vector<int>             avg_times_;
    vector<float> best_fit_alphas_;
    vector<float> best_fit_alphas_err_;
    vector<float> best_fit_scales_;
    vector<float> best_fit_scales_err_;
    vector<float> best_fit_pn_corr_;
    vector<float> best_fit_pn_corr_err_;
    float best_fit_chi2_;
    
};
    

#ifndef _UTILS_H_
#define _UTILS_H_

using namespace std;

namespace PhiSym
{   
//----------compute std::vector averages in a selected range------------------------------
    pair<float, float> VectorMeanRMS(vector<float>& data, int low, int high)
    {
        float sum=0;
        float sum2=0;
        int n=0;
        for(int i=low; i<high; ++i)
        {
            sum += data[i];
            sum2 += data[i]*data[i];
            ++n;
        }

        return make_pair(sum/n, sqrt(sum2/n - sum*sum/(n*n)));
    }

//----------compute std::vector averages in a selected range------------------------------
    pair<float, float> VectorWeightsMeanRMS(vector<pair<float, float> >& data, int low, int high)
    {
        float sum=0;
        float sum2=0;
        float w_sum=0;
        int n=0;
        for(int i=low; i<high; ++i)
        {
            float value=data[i].first*data[i].second;
            sum += value;
            sum2 += value*value;
            w_sum += data[i].second;
            ++n;
        }

        return make_pair(sum/w_sum, sqrt(sum2/w_sum - sum*sum/(w_sum*w_sum)));
    }

    
//----------compute std::vector outliers iterative removal--------------------------------
    pair<float, float> IterativeCut(vector<float>& data, int low, int high, double eps)
    {
        if(low >= high)
            return make_pair(-1, -1);
    
        pair<float, float> current = VectorMeanRMS(data, low, high);
        if(fabs(current.first - VectorMeanRMS(data, low+1, high-1).first) < eps)
            return current;
        else 
            return IterativeCut(data, low+1, high-1, eps);
    }

//----------compute std::vector outliers iterative removal--------------------------------
    pair<float, float> IterativeCutWeights(vector<pair<float, float> >& data, int low, int high, double eps)
    {
        if(low >= high)
            return make_pair(-1, -1);
    
        pair<float, float> current = VectorWeightsMeanRMS(data, low, high);
        if(fabs(current.first - VectorWeightsMeanRMS(data, low+1, high-1).first) < eps)
            return current;
        else 
            return IterativeCutWeights(data, low+1, high-1, eps);
    }
    
//----------compute the TH1 effctive sigma------------------------------------------------
    float EffectiveSigma(TH1* histo)
    {
        int rmsBins[2]={0, 0};       
        int meanBin = histo->FindBin(histo->GetMean());
        int range = min(meanBin-1, histo->GetNbinsX()-meanBin);
        float binWidth = histo->GetBinWidth(1);
        float grossIntegral=0;
        for(int iBin=1; iBin<range; ++iBin)
        {
            if(histo->Integral(meanBin-iBin, meanBin+iBin) > 0.683*histo->Integral())
            {
                grossIntegral = histo->Integral(meanBin-iBin, meanBin+iBin);
                rmsBins[0] = meanBin-iBin;
                rmsBins[1] = meanBin+iBin;
                break;
            }
        }
        for(int i=0; i<100; ++i)
        {
            if(grossIntegral-histo->GetBinContent(rmsBins[0])*i-histo->GetBinContent(rmsBins[1])*i < 0.683*histo->Integral())
                return histo->GetBinCenter(meanBin)-(histo->GetBinCenter(rmsBins[0])+(i-50)*binWidth/100);
        }
        return -1;
    }
}

#endif

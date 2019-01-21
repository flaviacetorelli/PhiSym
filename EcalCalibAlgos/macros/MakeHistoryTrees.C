#include <map>
#include <vector>
#include "TTree.h"

template<class INT, class OUTT>
void MakeHistoryTreesEB(std::vector<INT> eb_iovs, std::vector<INT> eb_even_iovs, std::vector<INT> eb_odd_iovs,
                        OUTT eb_ratio, vector<double>& etSumEB, map<int, vector<double> >& etSumPN)
{
    etSumEB.resize(eb_iovs.size(), 0);
    for(int iV=0; iV<eb_iovs.size(); ++iV)
    { 
        while(eb_iovs[iV]->NextEntry())
        { 
            if(etSumPN.find(WhichPN(eb_iovs[iV]->ieta, eb_iovs[iV]->iphi)) == etSumPN.end())
                etSumPN[WhichPN(eb_iovs[iV]->ieta, eb_iovs[iV]->iphi)] = vector<double>(eb_iovs.size(), 0);
            etSumPN[WhichPN(eb_iovs[iV]->ieta, eb_iovs[iV]->iphi)][iV] += eb_iovs[iV]->rec_hit->GetSumEt();
            etSumEB[iV] += eb_iovs[iV]->rec_hit->GetSumEt();            
        }
    }
    for(int i=0; i<61200; ++i)
    {
        if(i % 10000 == 0)
            cout << "Crystal #" << i << endl;
        for(int iV=0; iV<eb_iovs.size(); ++iV)
        { 
            eb_iovs[iV]->NextEntry(); 
            eb_even_iovs[iV]->NextEntry(); 
            eb_odd_iovs[iV]->NextEntry();
        } 
        int ref=0; 
        for(int iV=ref; iV<eb_iovs.size(); ++iV)
        { 
            if(iV==ref) eb_ratio->status = eb_iovs[iV]->rec_hit->GetSumEt() == 0 ? 0 : 1; 
            eb_ratio->n_events[iV-ref] = eb_iovs[iV]->n_events; 
            eb_ratio->iov[iV-ref] = iV; 
            eb_ratio->avg_time[iV-ref] = eb_iovs[iV]->avg_time; 
            eb_ratio->hashId = i; 
            eb_ratio->ieta = eb_iovs[ref]->ieta; 
            eb_ratio->iphi = eb_iovs[ref]->iphi; 
            eb_ratio->pn = WhichPN(eb_ratio->ieta, eb_ratio->iphi); 
            eb_ratio->sm = eb_ratio->ieta>0 ? (eb_ratio->iphi-1)/20+1 : -(eb_ratio->iphi-1)/20-1; 
            eb_ratio->mean_bs_x[iV-ref] = eb_iovs[iV]->mean_bs_x; 
            eb_ratio->mean_bs_sigmax[iV-ref] = eb_iovs[iV]->mean_bs_sigmax; 
            eb_ratio->mean_bs_y[iV-ref] = eb_iovs[iV]->mean_bs_y; 
            eb_ratio->mean_bs_sigmay[iV-ref] = eb_iovs[iV]->mean_bs_sigmay; 
            eb_ratio->mean_bs_z[iV-ref] = eb_iovs[iV]->mean_bs_z; 
	    //            eb_ratio->mean_bs_sigmaz[iV-ref] = eb_iovs[iV]->mean_bs_sigmaz; 
            eb_ratio->ring_average[iV-ref] = eb_iovs[iV]->ring_average; 
            eb_ratio->ring_average_uncleaned[iV-ref] = eb_iovs[iV]->ring_average_uncleaned; 
            eb_ratio->ic_ratio_abs[iV-ref] = eb_iovs[iV]->ic_ch/eb_iovs[ref]->ic_ch; 
            eb_ratio->ic_ratio_rel[iV-ref] = eb_iovs[iV]->ic_ch/eb_iovs[iV==0?0:iV-1]->ic_ch; 
            eb_ratio->ic_precision[iV-ref] = eb_even_iovs[iV]->ic_ch/eb_odd_iovs[iV]->ic_ch; 
            eb_ratio->k_ratio_abs[iV-ref] = eb_iovs[iV]->k_ch/eb_iovs[ref]->k_ch; 
            eb_ratio->k_ratio_rel[iV-ref] = eb_iovs[iV]->k_ch/eb_iovs[iV==0?0:iV-1]->k_ch; 
            eb_ratio->mean_et_ratio[iV-ref] = (eb_iovs[iV]->rec_hit->GetSumEt()/eb_iovs[iV]->rec_hit->GetNhits())/(eb_iovs[ref]->rec_hit->GetSumEt()/eb_iovs[ref]->rec_hit->GetNhits()); 
            eb_ratio->eflow_pn_abs[iV-ref] = 1+((eb_iovs[iV]->rec_hit->GetSumEt()/etSumPN[eb_ratio->pn][iV])/(eb_iovs[ref]->rec_hit->GetSumEt()/etSumPN[eb_ratio->pn][ref])-1)/eb_iovs[iV]->k_ch; 
            eb_ratio->eflow_abs[iV-ref] = 1+((eb_iovs[iV]->rec_hit->GetSumEt()/etSumEB[iV])/(eb_iovs[ref]->rec_hit->GetSumEt()/etSumEB[ref])-1)/eb_iovs[iV]->k_ch; 
            //eb_ratio->eflow_norm[iV-ref] = 1+((eb_iovs[iV]->rec_hit->GetSumEt()/eb_iovs[iV]->eflow_norm)/(eb_iovs[ref]->rec_hit->GetSumEt()/eb_iovs[ref]->eflow_norm)-1)/eb_iovs[iV]->k_ch; 
            //eb_ratio->lc[iV-ref] = eb_iovs[iV]->rec_hit->GetNhits()>0 ? eb_iovs[iV]->rec_hit->GetLCSum()/eb_iovs[iV]->rec_hit->GetNhits() : 0;
        } 
        eb_ratio->n_iovs = eb_iovs.size()-ref; 
        eb_ratio->GetTTreePtr()->Fill();
    }
    eb_ratio->GetTTreePtr()->BuildIndex("hashId"); 
    eb_ratio->GetTTreePtr()->AddFriend("EB", "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/data/sorted_crystals_info_EB.root"); 
}

template<class INT, class OUTT>
void MakeHistoryTreesEE(std::vector<INT> ee_iovs, std::vector<INT> ee_even_iovs, std::vector<INT> ee_odd_iovs,
                        OUTT ee_ratio, vector<double>& etSumEB, map<int, vector<double> >& etSumPN)
{
    for(int i=0; i<ee_iovs[0]->GetTTreePtr()->GetEntriesFast(); ++i)
    { 
        if(i % 10000 == 0)
            cout << "Crystal #" << i << endl;
        for(int iV=0; iV<ee_iovs.size(); ++iV)
        { 
            ee_iovs[iV]->NextEntry(); 
            ee_even_iovs[iV]->NextEntry(); 
            ee_odd_iovs[iV]->NextEntry();
        } 
        int ref=1; 
        for(int iV=ref; iV<ee_iovs.size(); ++iV)
        { 
            ee_ratio->n_events[iV-ref] = ee_iovs[iV]->n_events; 
            ee_ratio->avg_time[iV-ref] = ee_iovs[iV]->avg_time; 
            ee_ratio->iov[iV-ref] = iV; 
            ee_ratio->hashId = i;        
            ee_ratio->ix = ee_iovs[ref]->ix; 
            ee_ratio->iy = ee_iovs[ref]->iy; 
            ee_ratio->iring = ee_iovs[ref]->iring; 
            ee_ratio->mean_bs_x[iV-ref] = ee_iovs[iV]->mean_bs_x; 
            ee_ratio->mean_bs_sigmax[iV-ref] = ee_iovs[iV]->mean_bs_sigmax; 
            ee_ratio->mean_bs_y[iV-ref] = ee_iovs[iV]->mean_bs_y; 
            ee_ratio->mean_bs_sigmay[iV-ref] = ee_iovs[iV]->mean_bs_sigmay; 
            ee_ratio->mean_bs_z[iV-ref] = ee_iovs[iV]->mean_bs_z; 
            //ee_ratio->mean_bs_sigmaz[iV-ref] = ee_iovs[iV]->mean_bs_sigmaz; 
            ee_ratio->ring_average[iV-ref] = ee_iovs[iV]->ring_average; 
            ee_ratio->ring_average_uncleaned[iV-ref] = ee_iovs[iV]->ring_average_uncleaned; 
            ee_ratio->ic_ratio_abs[iV-ref] = ee_iovs[iV]->ic_ch/ee_iovs[ref]->ic_ch; 
            ee_ratio->ic_ratio_rel[iV-ref] = ee_iovs[iV]->ic_ch/ee_iovs[iV==0?0:iV-1]->ic_ch; 
            ee_ratio->ic_precision[iV-ref] = ee_even_iovs[iV]->ic_ch/ee_odd_iovs[iV]->ic_ch;   
            ee_ratio->k_ratio_abs[iV-ref] = ee_iovs[iV]->k_ch/ee_iovs[ref]->k_ch; 
            ee_ratio->k_ratio_rel[iV-ref] = ee_iovs[iV]->k_ch/ee_iovs[iV==0?0:iV-1]->k_ch; 
            ee_ratio->eflow_abs[iV-ref] = 1+((ee_iovs[iV]->rec_hit->GetSumEt()/etSumEB[iV])/(ee_iovs[ref]->rec_hit->GetSumEt()/etSumEB[ref])-1)/ee_iovs[iV]->k_ch; 
            //ee_ratio->eflow_norm[iV-ref] = 1+((ee_iovs[iV]->rec_hit->GetSumEt()/ee_iovs[iV]->eflow_norm)/(ee_iovs[ref]->rec_hit->GetSumEt()/ee_iovs[ref]->eflow_norm)-1)/ee_iovs[iV]->k_ch; 
            //ee_ratio->lc[iV-ref] = ee_iovs[iV]->rec_hit->GetNhits()>0 ? ee_iovs[iV]->rec_hit->GetLCSum()/ee_iovs[iV]->rec_hit->GetNhits() : 0;
        } 
        ee_ratio->n_iovs = ee_iovs.size()-ref; 
        ee_ratio->GetTTreePtr()->Fill();
    } 
    ee_ratio->GetTTreePtr()->AddFriend("EE", "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/data/sorted_crystals_info_EE.root"); 
}

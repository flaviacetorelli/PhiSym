{
    string file_name = "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_newGT_v5/summed_251562_251562.root";
    
    gSystem->Load("libFWCoreFWLite.so"); 
    AutoLibraryLoader::enable();
    gSystem->Load("libDataFormatsFWLite.so");
    gSystem->Load("libDataFormatsEcalDetId.so");
    gSystem->Load("libPhiSymEcalCalibDataFormats.so");

    float sm_ic_mean[18][2]={0};
    float sm_ic_rms[18][2]={0};    
    float sm_ic_sum[18][2]={0};
    float sm_ic_sum2[18][2]={0};    
    int sm_n_alive[18][2]={0};
    float iphi_ic_sum[361][2]={0};
    float iphi_ic_sum2[361][2]={0};
    float iphi_n_alive[361][2]={0};
    float sm_ic_mean_b65[18][2]={0};
    float sm_ic_rms_b65[18][2]={0};    
    float sm_ic_sum_b65[18][2]={0};
    float sm_ic_sum2_b65[18][2]={0};    
    int sm_n_alive_b65[18][2]={0};
    float iphi_ic_sum_b65[361][2]={0};
    float iphi_ic_sum2_b65[361][2]={0};
    float iphi_n_alive_b65[361][2]={0};
    int n_hits[EBDetId::kSizeForDenseIndexing]={0};
    bool is_good[EBDetId::kSizeForDenseIndexing]={0};
    float ic_uncorr[EBDetId::kSizeForDenseIndexing];
    float corrections[361][2]={0};
    for(int iPhi=1; iPhi<=360; ++iPhi)
    {
        corrections[iPhi][0]=1;
        corrections[iPhi][1]=1;
    }
    pair<int, int> ebMap[EBDetId::kSizeForDenseIndexing];
        
    TGraphErrors* gr_uncorr_EBp = new TGraphErrors();
    TGraphErrors* gr_uncorr_EBm = new TGraphErrors();
    TGraphErrors* gr_sm_sub_EBp = new TGraphErrors();
    TGraphErrors* gr_sm_sub_EBm = new TGraphErrors();
    TGraphErrors* gr_uncorr_b65_EBp = new TGraphErrors();
    TGraphErrors* gr_uncorr_b65_EBm = new TGraphErrors();
    TGraphErrors* gr_sm_sub_b65_EBp = new TGraphErrors();
    TGraphErrors* gr_sm_sub_b65_EBm = new TGraphErrors();
    TH1F* sm_mean_EBm = new TH1F("sb_mean_EBm", "", 360, 1, 360);
    TH1F* sm_mean_EBp = new TH1F("sb_mean_EBp", "", 360, 1, 360);
    TH2F* map_ic_uncorr = new TH2F("map_ic_uncorr", "PhiSym uncorrected ICs map;#it{i#phi};#it{i#eta}", 360, 1, 360, 171, -85, 85);
    TH2F* map_ic_corr = new TH2F("map_ic_corr", "PhiSym corrected ICs map;#it{i#phi};#it{i#eta}", 360, 1, 360, 171, -85, 85);
    TH2F* map_corrections = new TH2F("map_corrections", "Corrections map;#it{i#phi};#it{i#eta}", 360, 1, 360, 171, -85, 85);
    
    TFile* file = TFile::Open(file_name.c_str(), "READ");
    CrystalsEBTree ebTree((TTree*)file->Get("eb_xstals"));

    while(ebTree.NextEntry())
    {
        int index = EBDetId(ebTree.ieta, ebTree.iphi).hashedIndex();
        if(abs(ebTree.ieta) > 65)
            continue;
        is_good[index] = ebTree.rec_hit->GetSumEt()>ebTree.bounds[0] && ebTree.rec_hit->GetSumEt()<ebTree.bounds[1];
        n_hits[index] = ebTree.rec_hit->GetNhits();
        ic_uncorr[index] = ebTree.ic_ch;
        ebMap[index] = make_pair(ebTree.ieta, ebTree.iphi);
        //---skip a couple of bad TT (probably recovered from 2012 --- TO BE CHECKED)
        if(ebTree.ieta > 35 && ebTree.ieta < 41 && ebTree.iphi > 310 && ebTree.iphi < 316)
            n_hits[index] = 0;
        if(ebTree.ieta > -50 && ebTree.ieta < -44 && ebTree.iphi > 10 && ebTree.iphi < 16)
            n_hits[index] = 0;
    }
    
    for(int index=0; index<EBDetId::kSizeForDenseIndexing; ++index)
    {
        if(n_hits[index] == 0 || ic_uncorr[index]<0.7 || ic_uncorr[index]>1.5)
            continue;
        //---inside Tracker Barrel
        if(abs(ebMap[index].first) < 64)
        {
            int sm = (ebMap[index].second-1) / 20;
            int side = ebMap[index].first < 0 ? 0 : 1;
            iphi_ic_sum[ebMap[index].second][side] += ic_uncorr[index];
            iphi_ic_sum2[ebMap[index].second][side] += ic_uncorr[index]*ic_uncorr[index];
            ++iphi_n_alive[ebMap[index].second][side];
            //---outliers rejection 
            if(!is_good[index])
                continue;
            //---SM boundries rejection
            if((ebMap[index].second % 20 == 0 && ebMap[index].first > 0) ||
               (ebMap[index].second % 20 == 1 && ebMap[index].first < 0))
                continue;
            sm_ic_sum[sm][side] += ic_uncorr[index];
            sm_ic_sum2[sm][side] += ic_uncorr[index]*ic_uncorr[index];
            ++sm_n_alive[sm][side];
        }
        //---outside Tracker Barrel
        else
        {            
            int sm = (ebMap[index].second-1) / 20;
            int side = ebMap[index].first < 0 ? 0 : 1;
            iphi_ic_sum_b65[ebMap[index].second][side] += ic_uncorr[index];
            iphi_ic_sum2_b65[ebMap[index].second][side] += ic_uncorr[index]*ic_uncorr[index];
            ++iphi_n_alive_b65[ebMap[index].second][side];
            //---outliers rejection 
            if(!is_good[index])
                continue;
            //---SM boundries rejection
            if((ebMap[index].second % 20 == 0 && ebMap[index].first > 0) ||
               (ebMap[index].second % 20 == 1 && ebMap[index].first < 0))
                continue;
            sm_ic_sum_b65[sm][side] += ic_uncorr[index];
            sm_ic_sum2_b65[sm][side] += ic_uncorr[index]*ic_uncorr[index];
            ++sm_n_alive_b65[sm][side];
        }
    }
    
    for(int iSM=0; iSM<18; ++iSM)
    {
        //---with TB
        sm_ic_mean[iSM][0] = sm_ic_sum[iSM][0]/sm_n_alive[iSM][0];
        sm_ic_mean[iSM][1] = sm_ic_sum[iSM][1]/sm_n_alive[iSM][1];        
        sm_ic_rms[iSM][0] = sqrt(sm_ic_sum2[iSM][0]/sm_n_alive[iSM][0]-pow(sm_ic_sum[iSM][0]/sm_n_alive[iSM][0], 2));
        sm_ic_rms[iSM][1] = sqrt(sm_ic_sum2[iSM][1]/sm_n_alive[iSM][1]-pow(sm_ic_sum[iSM][1]/sm_n_alive[iSM][1], 2));
        //---no TB
        sm_ic_mean_b65[iSM][0] = sm_ic_sum_b65[iSM][0]/sm_n_alive_b65[iSM][0];
        sm_ic_mean_b65[iSM][1] = sm_ic_sum_b65[iSM][1]/sm_n_alive_b65[iSM][1];        
        sm_ic_rms_b65[iSM][0] = sqrt(sm_ic_sum2_b65[iSM][0]/sm_n_alive_b65[iSM][0]-
                                     pow(sm_ic_sum_b65[iSM][0]/sm_n_alive_b65[iSM][0], 2));
        sm_ic_rms_b65[iSM][1] = sqrt(sm_ic_sum2_b65[iSM][1]/sm_n_alive_b65[iSM][1]-
                                     pow(sm_ic_sum_b65[iSM][1]/sm_n_alive_b65[iSM][1], 2));
    }

    for(int iPhi=1; iPhi<=360; ++iPhi)
    {
        float point_unc[2], error_unc[2];
        float point_sub[2], error_sub[2];
        //---with TB
        point_unc[0] = iphi_ic_sum[iPhi][0]/iphi_n_alive[iPhi][0];
        point_unc[1] = iphi_ic_sum[iPhi][1]/iphi_n_alive[iPhi][1];
        error_unc[0] = sqrt(iphi_ic_sum2[iPhi][0]/iphi_n_alive[iPhi][0] -
                        pow(iphi_ic_sum[iPhi][0]/iphi_n_alive[iPhi][0], 2))/sqrt(iphi_n_alive[iPhi][1]);
        error_unc[1] = sqrt(iphi_ic_sum2[iPhi][1]/iphi_n_alive[iPhi][1] -
                        pow(iphi_ic_sum[iPhi][1]/iphi_n_alive[iPhi][1], 2))/sqrt(iphi_n_alive[iPhi][1]);
        point_sub[0] = point_unc[0]/sm_ic_mean[(iPhi-1)/20][0];
        point_sub[1] = point_unc[1]/sm_ic_mean[(iPhi-1)/20][1];
        error_sub[0] = point_sub[0]*sqrt(pow(error_unc[0]/point_unc[0], 2) +
                                         pow(sm_ic_rms[(iPhi-1)/20][0]/sm_n_alive[(iPhi-1)/20][0]/sm_ic_mean[(iPhi-1)/20][0], 2));
        error_sub[1] = point_sub[1]*sqrt(pow(error_unc[1]/point_unc[1], 2) +
                                         pow(sm_ic_rms[(iPhi-1)/20][1]/sm_n_alive[(iPhi-1)/20][1]/sm_ic_mean[(iPhi-1)/20][1], 2));
        //---uncorr
        gr_uncorr_EBm->SetPoint(iPhi-1, iPhi, point_unc[0]);
        gr_uncorr_EBm->SetPointError(iPhi-1, 0, error_unc[0]);
        gr_uncorr_EBp->SetPoint(iPhi-1, iPhi, point_unc[1]);
        gr_uncorr_EBp->SetPointError(iPhi-1, 0, error_unc[1]);
        //---sm averages subtracted
        gr_sm_sub_EBm->SetPoint(iPhi-1, iPhi, point_sub[0]);
        gr_sm_sub_EBm->SetPointError(iPhi-1, 0, error_sub[0]);
        gr_sm_sub_EBp->SetPoint(iPhi-1, iPhi, point_sub[1]);
        gr_sm_sub_EBp->SetPointError(iPhi-1, 0, error_sub[1]);
        //---sm mean and rms
        sm_mean_EBm->SetBinContent(iPhi, 1);
        sm_mean_EBm->SetBinError(iPhi, sm_ic_rms[(iPhi-1)/20][0]);
        sm_mean_EBp->SetBinContent(iPhi, 1);
        sm_mean_EBp->SetBinError(iPhi, sm_ic_rms[(iPhi-1)/20][1]);

        //---no TB
        point_unc[0] = iphi_ic_sum_b65[iPhi][0]/iphi_n_alive_b65[iPhi][0];
        point_unc[1] = iphi_ic_sum_b65[iPhi][1]/iphi_n_alive_b65[iPhi][1];
        error_unc[0] = sqrt(iphi_ic_sum2_b65[iPhi][0]/iphi_n_alive_b65[iPhi][0] -
                        pow(iphi_ic_sum_b65[iPhi][0]/iphi_n_alive_b65[iPhi][0], 2))/sqrt(iphi_n_alive_b65[iPhi][1]);
        error_unc[1] = sqrt(iphi_ic_sum2_b65[iPhi][1]/iphi_n_alive_b65[iPhi][1] -
                        pow(iphi_ic_sum_b65[iPhi][1]/iphi_n_alive_b65[iPhi][1], 2))/sqrt(iphi_n_alive_b65[iPhi][1]);
        point_sub[0] = point_unc[0]/sm_ic_mean_b65[(iPhi-1)/20][0];
        point_sub[1] = point_unc[1]/sm_ic_mean_b65[(iPhi-1)/20][1];
        error_sub[0] = point_sub[0]*sqrt(pow(error_unc[0]/point_unc[0], 2) +
                                         pow(sm_ic_rms_b65[(iPhi-1)/20][0]/
                                             sm_n_alive_b65[(iPhi-1)/20][0]/sm_ic_mean_b65[(iPhi-1)/20][0], 2));
        error_sub[1] = point_sub[1]*sqrt(pow(error_unc[1]/point_unc[1], 2) +
                                         pow(sm_ic_rms_b65[(iPhi-1)/20][1]/
                                             sm_n_alive_b65[(iPhi-1)/20][1]/sm_ic_mean_b65[(iPhi-1)/20][1], 2));
        //---uncorr
        gr_uncorr_b65_EBm->SetPoint(iPhi-1, iPhi, point_unc[0]);
        gr_uncorr_b65_EBm->SetPointError(iPhi-1, 0, error_unc[0]);
        gr_uncorr_b65_EBp->SetPoint(iPhi-1, iPhi, point_unc[1]);
        gr_uncorr_b65_EBp->SetPointError(iPhi-1, 0, error_unc[1]);
        //---sm averages subtracted
        gr_sm_sub_b65_EBm->SetPoint(iPhi-1, iPhi, point_sub[0]);
        gr_sm_sub_b65_EBm->SetPointError(iPhi-1, 0, error_sub[0]);
        gr_sm_sub_b65_EBp->SetPoint(iPhi-1, iPhi, point_sub[1]);
        gr_sm_sub_b65_EBp->SetPointError(iPhi-1, 0, error_sub[1]);
    }

    //---EB+ / EB- comparison
    Double_t* ebm_corr=gr_sm_sub_EBm->GetY();
    Double_t* ebm_corr_err=gr_sm_sub_EBm->GetEY();
    Double_t* ebp_corr=gr_sm_sub_EBp->GetY();
    Double_t* ebp_corr_err=gr_sm_sub_EBp->GetEY();

    TH1F* corr_syst=new TH1F("corr_syst", "Correction syst", 200, -0.1, 0.1);
    TH1F* corr_syst_err=new TH1F("corr_syst_err", "Correction syst error", 200, 0, 0.02);
    TH1F* corr_syst_pull=new TH1F("corr_syst_pull", "Correction pull", 40, -10, 10);
  
    Double_t ebm_corr_diff[360];
    Double_t ebm_corr_diff_err[360];
    
    for(int iPhi=0; iPhi<360; ++iPhi)
    {
        //---compute the bloody correction
        corrections[iPhi+1][0] = sm_ic_mean[iPhi/20][0]/ebm_corr[iPhi];
        corrections[iPhi+1][1] = sm_ic_mean[iPhi/20][1]/ebp_corr[iPhi];
        if((iPhi+1) % 20 == 1)
            corrections[iPhi+1][0] = sm_ic_mean[iPhi/20][0]/ebm_corr[iPhi];
        if((iPhi+1) % 20 == 0)
            corrections[iPhi+1][1] = sm_ic_mean[iPhi/20][1]/ebp_corr[iPhi];
        
        if((iPhi%20==0) || (iPhi%20==1) || (iPhi%20==19))
            continue;
        ebm_corr_diff[iPhi]=(*(ebm_corr+iPhi))-(*(ebp_corr+iPhi));
        ebm_corr_diff_err[iPhi]=TMath::Sqrt((*(ebm_corr_err+iPhi))*(*(ebm_corr_err+iPhi))+
                                            (*(ebp_corr_err+iPhi))*(*(ebp_corr_err+iPhi)));
        corr_syst->Fill(ebm_corr_diff[iPhi]);
        corr_syst_err->Fill(ebm_corr_diff_err[iPhi]);
        corr_syst_pull->Fill(ebm_corr_diff[iPhi]/ebm_corr_diff_err[iPhi]);            
    }    

    for(int index=0; index<EBDetId::kSizeForDenseIndexing; ++index)
    {
        double correction=corrections[ebMap[index].second][ebMap[index].first<0 ? 0 : 1];
        map_ic_uncorr->Fill(ebMap[index].second, ebMap[index].first, ic_uncorr[index]);
        map_ic_corr->Fill(ebMap[index].second, ebMap[index].first, ic_uncorr[index]*correction);
        map_corrections->Fill(ebMap[index].second, ebMap[index].first, correction);
    }

    //---output plots
    TFile* outFile = TFile::Open("geo_and_meterial_corr.root", "RECREATE");
    outFile->cd();
    //---style
    gr_uncorr_EBm->SetTitle("PhiSym IC vs phi - uncorrectred;#it{i#phi};#it{IC}");
    gr_uncorr_EBm->SetMarkerColor(kBlue+1);
    gr_uncorr_EBm->SetLineColor(kBlue+1);
    gr_uncorr_EBp->SetTitle("PhiSym IC vs phi - uncorrectred;#it{i#phi};#it{IC}");
    gr_uncorr_EBp->SetMarkerColor(kRed+1);
    gr_uncorr_EBp->SetLineColor(kRed+1);
    gr_sm_sub_EBm->SetTitle("PhiSym IC vs phi - uncorrectred (SM averages subtracted);#it{i#phi};#it{IC}");
    gr_sm_sub_EBm->SetMarkerColor(kBlue+1);
    gr_sm_sub_EBm->SetLineColor(kBlue+1);
    gr_sm_sub_EBp->SetTitle("PhiSym IC vs phi - uncorrectred (SM averages subtracted);#it{i#phi};#it{IC}");
    gr_sm_sub_EBp->SetMarkerColor(kRed+1);
    gr_sm_sub_EBp->SetLineColor(kRed+1);
    sm_mean_EBm->SetFillStyle(3004);
    sm_mean_EBm->SetFillColor(kBlue+2);
    sm_mean_EBp->SetFillStyle(3005);
    sm_mean_EBp->SetFillColor(kRed+2);
    gr_uncorr_b65_EBm->SetTitle("PhiSym IC vs phi - uncorrectred;#it{i#phi};#it{IC}");
    gr_uncorr_b65_EBm->SetMarkerColor(kBlue+1);
    gr_uncorr_b65_EBm->SetLineColor(kBlue+1);
    gr_uncorr_b65_EBp->SetTitle("PhiSym IC vs phi - uncorrectred;#it{i#phi};#it{IC}");
    gr_uncorr_b65_EBp->SetMarkerColor(kRed+1);
    gr_uncorr_b65_EBp->SetLineColor(kRed+1);
    gr_sm_sub_b65_EBm->SetTitle("PhiSym IC vs phi - uncorrectred (SM averages subtracted);#it{i#phi};#it{IC}");
    gr_sm_sub_b65_EBm->SetMarkerColor(kBlue+1);
    gr_sm_sub_b65_EBm->SetLineColor(kBlue+1);
    gr_sm_sub_b65_EBp->SetTitle("PhiSym IC vs phi - uncorrectred (SM averages subtracted);#it{i#phi};#it{IC}");
    gr_sm_sub_b65_EBp->SetMarkerColor(kRed+1);
    gr_sm_sub_b65_EBp->SetLineColor(kRed+1);
    map_ic_uncorr->SetContour(100);
    map_ic_uncorr->SetAxisRange(0.95, 1.05, "Z");    
    map_ic_corr->SetContour(100);
    map_ic_corr->SetAxisRange(0.95, 1.05, "Z");
    map_corrections->SetContour(100);
    map_corrections->SetAxisRange(0.9, 1.1, "Z");
    //---write histos
    gr_uncorr_EBm->Write("ic_uncorr_EBm");
    gr_uncorr_EBp->Write("ic_uncorr_EBp");
    gr_sm_sub_EBm->Write("ic_sm_sub_EBm");
    gr_sm_sub_EBp->Write("ic_sm_sub_EBp");
    sm_mean_EBm->Write("sm_mean_EBm");
    sm_mean_EBp->Write("sm_mean_EBp");
    corr_syst->Write("corr_syst");
    corr_syst_err->Write("corr_syst_err");
    corr_syst_pull->Write("corr_syst_pull");
    map_ic_uncorr->Write("map_ic_uncorr");
    map_ic_corr->Write("map_ic_corr");
    map_corrections->Write("map_corrections");
}
    

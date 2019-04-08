#ifndef __MATERIAL_CORRECTION__
#define __MATERIAL_CORRECTION__

#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm> 

#include "TSystem.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TH2F.h"
#include "TGraphErrors.h"

#include "FWCore/FWLite/interface/FWLiteEnabler.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

#include "Calibration/Tools/interface/EcalRingCalibrationTools.h"

#include "PhiSym/EcalCalibDataFormats/interface/PhiSymRecHit.h"
#include "PhiSym/EcalCalibDataFormats/interface/CalibrationFile.h"
#include "PhiSym/EcalCalibAlgos/interface/utils.h"

using namespace std;

//**********MAIN**************************************************************************
int main(int argc, char *argv[])
{
    FWLiteEnabler::enable();
    gSystem->Load("libFWCoreFWLite.so"); 
    gSystem->Load("libDataFormatsFWLite.so");
    gSystem->Load("libDataFormatsEcalDetId.so");
    gSystem->Load("libPhiSymEcalCalibDataFormats.so");

    //---Setup---
    if(argc < 2)
    {
        cout << "Usage : " << argv[0] << " [parameters.py]" << endl;
        return 0;
    }
    if(!edm::boost_python::readPSetsFrom(argv[1])->existsAs<edm::ParameterSet>("process"))
    {
        cout << " ERROR: ParametersSet 'process' is missing in your configuration file"
             << endl;
        return 0;
    }

    //---get the python configuration
    const edm::ParameterSet &process = edm::boost_python::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("process");
    const edm::ParameterSet &filesOpt = process.getParameter<edm::ParameterSet>("ioFilesOpt");
    
    bool absoluteICs = process.getParameter<bool>("absoluteICs");
    bool userOutputName = filesOpt.getParameter<bool>("userOutputName");
    string outputFileBase = filesOpt.getParameter<string>("outputFileBase");
    vector<string> outputFiles = filesOpt.getParameter<vector<string> >("outputFiles");
    vector<string> inputFiles = filesOpt.getParameter<vector<string> >("inputFiles");

    //---loop over input files (one input -> one output)
    int iFile=-1;
    for(auto& inputFile: inputFiles)
    {
        ++iFile;
        string runsRange(inputFile.begin()+inputFile.find_first_of("_", inputFile.find_last_of("/"))+1,
                         inputFile.end()-5);

        //---init
        float sm_ic_mean[18][2]={};
        float sm_ic_rms[18][2]={};    
        int sm_n_alive[18][2]={};
        vector<float> sm_ics[18][2];
        vector<float> iphi_ic[361][2];
        float sm_ic_mean_b60[18][2]={};
        float sm_ic_rms_b60[18][2]={};
        int sm_n_alive_b60[18][2]={};
        vector<float> sm_ics_b60[18][2];
        vector<float> iphi_ic_b60[361][2];
        int n_hits[EBDetId::kSizeForDenseIndexing]={0};
        // bool is_good[EBDetId::kSizeForDenseIndexing]={0};
        float ic_uncorr[EBDetId::kSizeForDenseIndexing]={0};
        float corrections[EBDetId::kSizeForDenseIndexing]={1};   
        pair<int, int> ebMap[EBDetId::kSizeForDenseIndexing];

        TGraphErrors* gr_uncorr_EBp = new TGraphErrors();
        TGraphErrors* gr_uncorr_EBm = new TGraphErrors();
        TGraphErrors* gr_sm_sub_EBp = new TGraphErrors();
        TGraphErrors* gr_sm_sub_EBm = new TGraphErrors();
        TGraphErrors* gr_uncorr_b60_EBp = new TGraphErrors();
        TGraphErrors* gr_uncorr_b60_EBm = new TGraphErrors();
        TGraphErrors* gr_sm_sub_b60_EBp = new TGraphErrors();
	TGraphErrors* gr_sm_sub_b60_EBm = new TGraphErrors();

        TH1F* sm_mean_EBm = new TH1F("sb_mean_EBm", "", 360, 0.5, 360.5);
        TH1F* sm_mean_EBp = new TH1F("sb_mean_EBp", "", 360, 0.5, 360.5);
        TH1F* h_ic_uncorr = new TH1F("h_ic_uncorr", "PhiSym uncorrected ICs ;#it{IC_{uncorr}};", 100, 0.95, 1.05);
        TH1F* h_ic_corr = new TH1F("h_ic_corr", "PhiSym corrected ICs ;#it{IC_{corr}};", 100, 0.95, 1.05);
        TH2F* map_ic_uncorr = new TH2F("map_ic_uncorr", "PhiSym uncorrected ICs map;#it{i#phi};#it{i#eta}",
                                       360, 0.5, 360.5, 171, -85.5, 85.5);
        TH2F* map_ic_corr = new TH2F("map_ic_corr", "PhiSym corrected ICs map;#it{i#phi};#it{i#eta}",
                                     360, 0.5, 360.5, 171, -85.5, 85.5);
        TH2F* map_ic_corr_b60 = new TH2F("map_ic_corr_b60", "PhiSym corrected ICs map;#it{i#phi};#it{i#eta}",
                                     360, 0.5, 360.5, 121, -60.5, 60.5);
        TH2F* map_corrections = new TH2F("map_corrections", "Corrections map;#it{i#phi};#it{i#eta}",
                                         360, 0.5, 360.5, 171, -85.5, 85.5);
        TH2F* map_EEp = new TH2F("map_EEp", "IC map EE+;#it{ix};#it{iy}",
                                         101, -0.5, 100.5, 100, 0.5, 100.5);
        TH2F* map_EEm = new TH2F("map_EEm", "IC map EE-;#it{ix};#it{iy}",
                                         101, -0.5, 100.5, 100, 0.5, 100.5);
        TH2F* map_lc_EB = new TH2F("map_lc_EB", "LC map EB;#it{i#phi};#it{i#eta}",
                                         360, 0.5, 360.5, 171, -85.5, 85.5);
        TH2F* map_lc_EEp = new TH2F("map_lc_EEp", "LC map EE+;#it{ix};#it{iy}",
                                         101, -0.5, 100.5, 100, 0.5, 100.5);
        TH2F* map_lc_EEm = new TH2F("map_lc_EEm", "LC map EE-;#it{ix};#it{iy}",
                                         101, -0.5, 100.5, 100, 0.5, 100.5);

       
        TFile* file = TFile::Open(inputFile.c_str(), "READ");
        CrystalsEBTree ebTree((TTree*)file->Get("eb_xstals"));
        CrystalsEETree eeTree((TTree*)file->Get("ee_xstals"));

        for(int iPhi=1; iPhi<=360; ++iPhi)
        {
            iphi_ic[iPhi][0].clear();
            iphi_ic[iPhi][1].clear();
            iphi_ic_b60[iPhi][0].clear();
            iphi_ic_b60[iPhi][1].clear();
        }
        
        while(ebTree.NextEntry())
        {
            int index = EBDetId(ebTree.ieta, ebTree.iphi).hashedIndex();
            // is_good[index] = ebTree.rec_hit->GetSumEt() > ebTree.bounds[0] &&
            //     ebTree.rec_hit->GetSumEt() < ebTree.bounds[1];
            n_hits[index] = ebTree.rec_hit->GetNhits();            
            if(absoluteICs)
                ic_uncorr[index] = ebTree.ic_abs*ebTree.ic_ch;
            else
                ic_uncorr[index] = ebTree.ic_ch;
            ebMap[index] = make_pair(ebTree.ieta, ebTree.iphi);
            //---skip a couple of bad TT (probably recovered from 2012 --- TO BE CHECKED)
            if(ebTree.ieta > 35 && ebTree.ieta < 41 && ebTree.iphi > 310 && ebTree.iphi < 316)
                n_hits[index] = 0;
            if(ebTree.ieta > -50 && ebTree.ieta < -44 && ebTree.iphi > 10 && ebTree.iphi < 16)
                n_hits[index] = 0;

	    map_lc_EB -> Fill(ebTree.iphi, ebTree.ieta, ebTree.rec_hit->GetLCSum()/ebTree.rec_hit->GetNhits());
        }

        while(eeTree.NextEntry())
        {
            float ic_EE = eeTree.ic_ch;
            int ix = eeTree.ix;
	    int iy = eeTree.iy;
	    int iring = eeTree.iring;

	    if(iring>0)
	    {  map_EEp->Fill(ix, iy, ic_EE);
	       map_lc_EEp->Fill(ix, iy,  eeTree.rec_hit->GetLCSum()/eeTree.rec_hit->GetNhits());
            }
            else
	    {  map_EEm->Fill(ix, iy, ic_EE);
	       map_lc_EEm->Fill(ix, iy,  eeTree.rec_hit->GetLCSum()/eeTree.rec_hit->GetNhits());
            }

        }

        //---SM averages
        for(int index=0; index<EBDetId::kSizeForDenseIndexing; ++index)
        {
            if(n_hits[index] == 0)
                continue;
            int sm = (ebMap[index].second-1) / 20;
            int side = ebMap[index].first < 0 ? 0 : 1;
            if(abs(ebMap[index].first) < 60)
                iphi_ic[ebMap[index].second][side].push_back(ic_uncorr[index]);
            else
                iphi_ic_b60[ebMap[index].second][side].push_back(ic_uncorr[index]);
            //---SM boundries rejection
            if(!(ebMap[index].second % 20 == 0 && ebMap[index].first > 0) &&
               !(ebMap[index].second % 20 == 1 && ebMap[index].first < 0))
            {
                //---inside Tracker Barrel
                if(abs(ebMap[index].first) < 60)
                    sm_ics[sm][side].push_back(ic_uncorr[index]);
                //---outside Tracker Barrel
                else
                    sm_ics_b60[sm][side].push_back(ic_uncorr[index]);
            }
        }
    
        for(int iSM=0; iSM<18; ++iSM)
        {
            for(int iSide=0; iSide<2; ++iSide)
            {
                //---with TB
                sort(sm_ics[iSM][iSide].begin(), sm_ics[iSM][iSide].end());
                pair<float, float> itCutResults = PhiSym::IterativeCut(sm_ics[iSM][iSide], 0, sm_ics[iSM][iSide].size(), 0.001);
                for(auto& ic : sm_ics[iSM][iSide])
                    if(ic > itCutResults.first-2*itCutResults.second && ic < itCutResults.first+2*itCutResults.second)
                        ++sm_n_alive[iSM][iSide];

                sm_ic_mean[iSM][iSide] = itCutResults.first;
                sm_ic_rms[iSM][iSide] = itCutResults.second;
                //---no TB
                sort(sm_ics_b60[iSM][iSide].begin(), sm_ics_b60[iSM][iSide].end());
                itCutResults = PhiSym::IterativeCut(sm_ics_b60[iSM][iSide], 0, sm_ics_b60[iSM][iSide].size(), 0.001);
                for(auto& ic : sm_ics_b60[iSM][iSide])
                    if(ic > itCutResults.first-2*itCutResults.second && ic < itCutResults.first+2*itCutResults.second)
                        ++sm_n_alive_b60[iSM][iSide];
                sm_ic_mean_b60[iSM][iSide] = itCutResults.first;
                sm_ic_rms_b60[iSM][iSide] = itCutResults.second;
            }
        }

        for(int iPhi=1; iPhi<=360; ++iPhi)
        {
            float point_unc[2], error_unc[2];
            float point_sub[2], error_sub[2];
            //---with TB
            sort(iphi_ic[iPhi][0].begin(), iphi_ic[iPhi][0].end());
            sort(iphi_ic[iPhi][1].begin(), iphi_ic[iPhi][1].end());
            pair<float, float> tmp_ebm = PhiSym::IterativeCut(iphi_ic[iPhi][0], 0, iphi_ic[iPhi][0].size(), 0.001);
            pair<float, float> tmp_ebp = PhiSym::IterativeCut(iphi_ic[iPhi][1], 0, iphi_ic[iPhi][1].size(), 0.001);
            point_unc[0] = tmp_ebm.first;
            point_unc[1] = tmp_ebp.first;
            error_unc[0] = tmp_ebm.second;
            error_unc[1] = tmp_ebp.second;
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
            sort(iphi_ic_b60[iPhi][0].begin(), iphi_ic_b60[iPhi][0].end());
            sort(iphi_ic_b60[iPhi][1].begin(), iphi_ic_b60[iPhi][1].end());
            tmp_ebm = PhiSym::IterativeCut(iphi_ic_b60[iPhi][0], 0, iphi_ic_b60[iPhi][0].size(), 0.005);
            tmp_ebp = PhiSym::IterativeCut(iphi_ic_b60[iPhi][1], 0, iphi_ic_b60[iPhi][1].size(), 0.005);
            point_unc[0] = tmp_ebm.first;
            point_unc[1] = tmp_ebp.first;
            error_unc[0] = tmp_ebm.second;
            error_unc[1] = tmp_ebp.second;
            point_sub[0] = point_unc[0]/sm_ic_mean_b60[(iPhi-1)/20][0];
            point_sub[1] = point_unc[1]/sm_ic_mean_b60[(iPhi-1)/20][1];
            error_sub[0] = point_sub[0]*sqrt(pow(error_unc[0]/point_unc[0], 2) +
                                             pow(sm_ic_rms_b60[(iPhi-1)/20][0]/
                                                 sm_n_alive_b60[(iPhi-1)/20][0]/sm_ic_mean_b60[(iPhi-1)/20][0], 2));
            error_sub[1] = point_sub[1]*sqrt(pow(error_unc[1]/point_unc[1], 2) +
                                             pow(sm_ic_rms_b60[(iPhi-1)/20][1]/
                                                 sm_n_alive_b60[(iPhi-1)/20][1]/sm_ic_mean_b60[(iPhi-1)/20][1], 2));
            //---uncorr
            gr_uncorr_b60_EBm->SetPoint(iPhi-1, iPhi, point_unc[0]);
            gr_uncorr_b60_EBm->SetPointError(iPhi-1, 0, error_unc[0]);
            gr_uncorr_b60_EBp->SetPoint(iPhi-1, iPhi, point_unc[1]);
            gr_uncorr_b60_EBp->SetPointError(iPhi-1, 0, error_unc[1]);
            //---sm averages subtracted
            gr_sm_sub_b60_EBm->SetPoint(iPhi-1, iPhi, point_sub[0]);
            gr_sm_sub_b60_EBm->SetPointError(iPhi-1, 0, error_sub[0]);
            gr_sm_sub_b60_EBp->SetPoint(iPhi-1, iPhi, point_sub[1]);
            gr_sm_sub_b60_EBp->SetPointError(iPhi-1, 0, error_sub[1]);
        }

        //---EB+ / EB- comparison
        Double_t* ebm_corr=gr_sm_sub_EBm->GetY();
        Double_t* ebm_corr_err=gr_sm_sub_EBm->GetEY();
        Double_t* ebp_corr=gr_sm_sub_EBp->GetY();
        Double_t* ebp_corr_err=gr_sm_sub_EBp->GetEY();
        Double_t* ebm_corr_b60=gr_sm_sub_b60_EBm->GetY();
        //Double_t* ebm_corr_err_b60=gr_sm_sub_b60_EBm->GetEY();
        Double_t* ebp_corr_b60=gr_sm_sub_b60_EBp->GetY();
        //Double_t* ebp_corr_err_b60=gr_sm_sub_b60_EBp->GetEY();

        TH1F* corr_syst=new TH1F("corr_syst", "Correction syst", 200, -0.1, 0.1);
        TH1F* corr_syst_err=new TH1F("corr_syst_err", "Correction syst error", 200, 0, 0.02);
        TH1F* corr_syst_pull=new TH1F("corr_syst_pull", "Correction pull", 40, -10, 10);
  
        Double_t ebm_corr_diff[360];
        Double_t ebm_corr_diff_err[360];
    
        for(int index=0; index<EBDetId::kSizeForDenseIndexing; ++index)
        {
            int ieta = ebMap[index].first;
            int iphi = ebMap[index].second;
            //---compute the bloody correction
            //---with TB
            if(abs(ieta) < 60)
            {
                if(ieta < 0)
                    corrections[index] = 1/ebm_corr[iphi-1];
                if(ieta > 0)
                    corrections[index] = 1/ebp_corr[iphi-1];
            }
            //---no TB
            else
            {
                if(ieta < 0)
                    corrections[index] = 1/ebm_corr_b60[iphi-1];
                if(ieta > 0)
                    corrections[index] = 1/ebp_corr_b60[iphi-1];
            }
        
            if((iphi-1)%20 == 0 || (iphi-1)%20 == 1)
                continue;
            ebm_corr_diff[iphi-1]=(*(ebm_corr+iphi-1))-(*(ebp_corr+iphi-1));
            ebm_corr_diff_err[iphi-1]=TMath::Sqrt((*(ebm_corr_err+iphi-1))*(*(ebm_corr_err+iphi-1))+
                                                  (*(ebp_corr_err+iphi-1))*(*(ebp_corr_err+iphi-1)));
            corr_syst->Fill(ebm_corr_diff[iphi-1]);
            corr_syst_err->Fill(ebm_corr_diff_err[iphi-1]);
            corr_syst_pull->Fill(ebm_corr_diff[iphi-1]/ebm_corr_diff_err[iphi-1]);            
        }    

        //---generate txt corrections file
        ofstream outTxtFile((outputFileBase+runsRange+".txt").c_str(), ios::out);
        for(int index=0; index<EBDetId::kSizeForDenseIndexing; ++index)
        {
            double correction=corrections[index];
            h_ic_uncorr->Fill(ic_uncorr[index]);
            h_ic_corr->Fill(ic_uncorr[index]*correction);
            map_ic_uncorr->Fill(ebMap[index].second, ebMap[index].first, ic_uncorr[index]);
            map_ic_corr->Fill(ebMap[index].second, ebMap[index].first, ic_uncorr[index]*correction);
	    if(abs(ebMap[index].first)<60)
	            map_ic_corr_b60->Fill(ebMap[index].second, ebMap[index].first, ic_uncorr[index]*correction);
            map_corrections->Fill(ebMap[index].second, ebMap[index].first, correction);
            outTxtFile << ebMap[index].first << "   " << ebMap[index].second << "   0   " << correction << endl;

        }

        outTxtFile.close();

	TProfile* test  = map_ic_corr -> ProfileX();
	TProfile* test2 = map_ic_corr -> ProfileX();
	TProfile* test3 = map_ic_corr -> ProfileX();


	for(int i=1; i<test2->GetNbinsX(); i++)
	{
		float ic = test2 -> GetBinContent(i);

		float sm_avg = sm_ic_mean[(i-1)/20][0];
		//cout << i << " " << sm_avg << endl;

		test2 -> SetBinContent(i, ic/sm_avg);
	}

	test3 -> Divide(test2);


	test -> SetMarkerStyle(20);
	test -> SetMarkerColor(kBlue);
	test -> GetYaxis() -> SetRangeUser(0.95, 1.05);

	test2 -> SetMarkerStyle(20);
	test2 -> SetMarkerColor(kBlue);
	test2 -> GetYaxis() -> SetRangeUser(0.95, 1.05);

	test3 -> SetMarkerStyle(20);
	test3 -> SetMarkerColor(kBlue);
	test3 -> GetYaxis() -> SetRangeUser(0.95, 1.05);

        //---output plots
        TFile* outFile;
        if(!userOutputName)
            outFile = TFile::Open((outputFileBase+runsRange+".root").c_str(), "RECREATE");
        else
            outFile = TFile::Open((outputFileBase+outputFiles[iFile]).c_str(), "RECREATE");
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
        gr_uncorr_b60_EBm->SetTitle("PhiSym IC vs phi - uncorrectred;#it{i#phi};#it{IC}");
        gr_uncorr_b60_EBm->SetMarkerColor(kBlue+1);
        gr_uncorr_b60_EBm->SetLineColor(kBlue+1);
        gr_uncorr_b60_EBp->SetTitle("PhiSym IC vs phi - uncorrectred;#it{i#phi};#it{IC}");
        gr_uncorr_b60_EBp->SetMarkerColor(kRed+1);
        gr_uncorr_b60_EBp->SetLineColor(kRed+1);
        gr_sm_sub_b60_EBm->SetTitle("PhiSym IC vs phi - uncorrectred (SM averages subtracted);#it{i#phi};#it{IC}");
        gr_sm_sub_b60_EBm->SetMarkerColor(kBlue+1);
        gr_sm_sub_b60_EBm->SetLineColor(kBlue+1);
        gr_sm_sub_b60_EBp->SetTitle("PhiSym IC vs phi - uncorrectred (SM averages subtracted);#it{i#phi};#it{IC}");
        gr_sm_sub_b60_EBp->SetMarkerColor(kRed+1);
        gr_sm_sub_b60_EBp->SetLineColor(kRed+1);
        map_ic_uncorr->SetContour(100000);
        map_ic_uncorr->SetAxisRange(0.95, 1.05, "Z");    
        map_ic_corr->SetContour(100000);
        map_ic_corr->SetAxisRange(0.95, 1.05, "Z");
        map_ic_corr_b60->SetContour(100000);
        map_ic_corr_b60->SetAxisRange(0.95, 1.05, "Z");
        map_corrections->SetContour(100000);
        map_corrections->SetAxisRange(0.9, 1.1, "Z");
        //---write histos
        gr_uncorr_EBm->Write("ic_uncorr_EBm");
        gr_uncorr_EBp->Write("ic_uncorr_EBp");
        gr_sm_sub_EBm->Write("ic_sm_sub_EBm");
        gr_sm_sub_EBp->Write("ic_sm_sub_EBp");
        gr_uncorr_b60_EBm->Write("ic_uncorr_b60_EBm");
        gr_uncorr_b60_EBp->Write("ic_uncorr_b60_EBp");
        gr_sm_sub_b60_EBm->Write("ic_sm_sub_b60_EBm");
        gr_sm_sub_b60_EBp->Write("ic_sm_sub_b60_EBp");	
	test -> Write("test");
	test2 -> Write("test2");
	test3 -> Write("test3");
        sm_mean_EBm->Write("sm_mean_EBm");
        sm_mean_EBp->Write("sm_mean_EBp");
        corr_syst->Write("corr_syst");
        corr_syst_err->Write("corr_syst_err");
        corr_syst_pull->Write("corr_syst_pull");
        h_ic_uncorr->Write("h_ic_uncorr");
        h_ic_corr->Write("h_ic_corr");
        map_ic_uncorr->Write("map_ic_uncorr");
        map_ic_corr->Write("map_ic_corr");
        map_ic_corr_b60->Write("map_ic_corr_b60");
        map_EEp->Write("map_ic_EEp");
        map_EEm->Write("map_ic_EEm");
        map_corrections->Write("map_corrections");
	map_lc_EB->Write("map_lc_EB");
	map_lc_EEm->Write("map_lc_EEm");
	map_lc_EEp->Write("map_lc_EEp");
        
        file->Close();
    }
}
    
#endif

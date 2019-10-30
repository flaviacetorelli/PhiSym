/*
 g++ -Wall -o plotter_prova `root-config --cflags --glibs` -L $ROOTSYS/lib  plotter_prova.cpp
 */


#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"

#include "TProfile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TChain.h"
#include "TVirtualFitter.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TMath.h"
#include "TCut.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace std;

void MakevsTimePlotsRing(TChain* t,TString var1_var2,  float y_min, float y_max, TString Ytitle, TString title)
{
    TCanvas c;
    gStyle->SetOptStat(0);
     TString cut = "ieta==10";
    t->Draw(var1_var2,cut,"prof");
    TH2F *gr = (TH2F*)gPad->GetPrimitive("htemp");
    gr->SetTitle("; Time(day/month); ");
    gr->GetXaxis()->SetTimeDisplay(1);
    gr->GetXaxis()->SetTimeFormat("%d/%m%F1970-01-01 00:00:00");
    gr->GetXaxis()->SetRangeUser(1.4975e+09,1505012800);
    gr->GetYaxis()->SetTitle(Ytitle);
    gr->GetYaxis()->SetRangeUser(y_min,y_max);
    gPad->Update();
    for ( int i=2; i<9; i++)
    {
    cut=("ieta==" + to_string(i*10)).c_str();
    t->SetMarkerColor(2+i);
    if (i==8) t->SetMarkerColor(11);
    t->Draw(var1_var2,cut,"prof same");
    }
    t->SetMarkerStyle(23);
    for ( int i=1; i<9; i++)
    {
    cut=("ieta==" + to_string(-i*10)).c_str();
    t->SetMarkerColor(1+i);
    if (i==8) t->SetMarkerColor(11);
    t->Draw(var1_var2,cut,"prof same");
    }
    
    c.SaveAs(title + ".png");
    c.SaveAs(title + ".pdf");
    c.SaveAs(title + ".root");



}


void MakevsTimePlots(TChain* t,TString var1_var2, TString cut,  float y_min, float y_max, TString Ytitle, TString title)
{
    gStyle->SetOptStat(0);
    TCanvas c;
    c.SetGrid();
    t->Draw(var1_var2,cut,"prof");
    TH2F *gr = (TH2F*)gPad->GetPrimitive("htemp");
    gr->SetTitle("; Time(day/month); ");
    gr->GetXaxis()->SetTimeDisplay(1);
    gr->GetXaxis()->SetTimeFormat("%d/%m%F1970-01-01 00:00:00");
    gr->GetXaxis()->SetRangeUser(1.49698e+09,1511490000);
    gr->GetYaxis()->SetTitle(Ytitle);
    gr->GetYaxis()->SetRangeUser(y_min,y_max);
    gPad->Update();
    c.SaveAs(title + ".png");
    c.SaveAs(title + ".pdf");



}

void MakevsTimePlots(TChain* t,TString var1_var2,  float y_min, float y_max, TString Ytitle, TString title)
{
    TCanvas c;
    t->Draw(var1_var2,"","prof");
    TH2F *gr = (TH2F*)gPad->GetPrimitive("htemp");
    gr->SetTitle("; Time(day/month); ");
    gr->GetXaxis()->SetTimeDisplay(1);
    gr->GetXaxis()->SetTimeFormat("%d/%m%F1970-01-01 00:00:00");
    gr->GetXaxis()->SetRangeUser(1.4975e+09,1505012800 );
    gr->GetYaxis()->SetTitle(Ytitle);
    gr->GetYaxis()->SetRangeUser(y_min,y_max);
    gPad->Update();
    c.SaveAs(title + ".png");



}



void MakevsTimePlots(TChain* t,TString var1_var2, TString Ytitle, TString title)
{
    TCanvas c;
    t->Draw(var1_var2,"","prof");
    TH2F *gr = (TH2F*)gPad->GetPrimitive("htemp");
    gr->SetTitle("; Time(day/month); ");
    gr->GetXaxis()->SetTimeDisplay(1);
    gr->GetXaxis()->SetTimeFormat("%d/%m%F1970-01-01 00:00:00");
    gr->GetXaxis()->SetRangeUser(1.4975e+09,1505012800);
    gr->GetYaxis()->SetTitle(Ytitle);
    gPad->Update();
    c.SaveAs(title + ".png");



}

void MakeHistoIC(TChain* t,float min, float max, int sup , TString title)
{
    TCanvas c;
    float ic_ratio_abs; 
    float ic_ratio_eflow; 
    int ieta; 
    int sm;
    TH1F* h = new TH1F("","",100, 0.8,1.2); 
    TH1F* h1 = new TH1F("","",100, 0.8,1.2); 
    t->SetBranchAddress("ic_ratio_abs",&ic_ratio_abs);
    t->SetBranchAddress("ic_ratio_eflow",&ic_ratio_eflow);
    t->SetBranchAddress("ieta",&ieta);
    t->SetBranchAddress("sm",&sm);
    for (int i=0; i<t->GetEntries(); i++)
    {
      t->GetEntry(i);

     if (float(ieta) > min && float(ieta) < max  && sm == sup)
     {           std::cout << "ciao4"<< std::endl;
       //std::cout << ic_eflow << std::endl;
       h->Fill(ic_ratio_eflow);
       h1->Fill(ic_ratio_abs);

     }
    }

    h1->SetTitle(";IC ; events ");
    h1->SetLineColor(kBlue);
    h->SetLineColor(kRed);
    h1->Draw("histo"); 
    h->Draw("histo same"); 
        TLegend leg;
    leg.AddEntry(h1, "ic ring", "l");
    leg.AddEntry(h, "ic eflow", "l");
    leg.Draw("same");
    c.SaveAs("IC_SM_1_" + title +".png");
    
    
 
    //c.SaveAs("IC_eflow.pdf");
    //c.SaveAs("IC_eflow.root");

}




int main(int argc, char *argv[])
{

    TChain* ZeroBias;
    TString fold = "/afs/cern.ch/work/f/fcetorel/private/work2/prova_slc6/CMSSW_9_4_0/src/PhiSym/EcalCalibAlgos/";
    //ZeroBias = new TChain("eb_xstals");
    //ZeroBias -> Add(fold + "ntuple_2017B/summed*.root" );
    //ZeroBias -> Add(fold + "ntuple_2017C/summed*.root" );
    ZeroBias = new TChain("eb");
    ZeroBias -> Add(fold + "macros/history_eflow_2017.root" );
    ZeroBias -> SetMarkerStyle(20);
    ZeroBias -> SetMarkerColor(kRed);


   // MakevsTimePlots(ZeroBias,"eflow_wnorm:avg_time","Energy(GeV)" , "ESumW");
   //MakevsTimePlotsRing(ZeroBias, "ring_average/eflow_wnorm:avg_time" , 0.9, 1.1," <E_{ring}>/<E_{EB}>", "ringObarr");

/*
    MakeHistoIC(ZeroBias,0,26,1,  "p_mod1");
    MakeHistoIC(ZeroBias,25,46, 1,  "p_mod2");
    MakeHistoIC(ZeroBias,45, 66, 1, "p_mod3");
    MakeHistoIC(ZeroBias,65, 86, 1, "p_mod4");
    MakeHistoIC(ZeroBias,-26,5, 1, "m_mod1");
    MakeHistoIC(ZeroBias,-46,-25, 1,  "m_mod2");
    MakeHistoIC(ZeroBias,-66, -45, 1, "m_mod3");
    MakeHistoIC(ZeroBias,-86, -65, 1, "m_mod4");*/
     MakevsTimePlots(ZeroBias,"ic_ratio_eflow:avg_time", "n_events>110e6 && ic_ratio_eflow > 0 && ic_ratio_eflow < 100",0.98, 1.01, "IC_{n}/IC_{0}", "IC_all");
   for (int i=1; i<19 ; i++)
    {
    MakevsTimePlots(ZeroBias,"ic_ratio_eflow:avg_time", ("n_events>110e6 && ic_ratio_eflow > 0 && ic_ratio_eflow < 100 && ieta>0 && ieta<26 && sm=="+to_string(i)).c_str(), 0.98, 1.01, "IC_{n}/IC_{0}", ("IC_SM_p"+to_string(i)+"_module1").c_str());
    MakevsTimePlots(ZeroBias,"ic_ratio_eflow:avg_time", ("n_events>110e6 && ic_ratio_eflow > 0 && ic_ratio_eflow < 100 && ieta<0 && ieta>-26 && sm==-"+to_string(i)).c_str(),0.98, 1.01, "IC_{n}/IC_{0}", ("IC_SM_m"+to_string(i)+"_module1").c_str())  ;
    MakevsTimePlots(ZeroBias,"ic_ratio_eflow:avg_time", ("n_events>110e6 && ic_ratio_eflow > 0 && ic_ratio_eflow < 100 && ieta>25 && ieta<46 && sm=="+to_string(i)).c_str(), 0.98, 1.01, "IC_{n}/IC_{0}", ("IC_SM_p"+to_string(i)+"_module2").c_str());
    MakevsTimePlots(ZeroBias,"ic_ratio_eflow:avg_time", ("n_events>110e6 && ic_ratio_eflow > 0 && ic_ratio_eflow < 100 && ieta<-25 && ieta>-46 && sm==-"+to_string(i)).c_str(),0.98, 1.01, "IC_{n}/IC_{0}", ("IC_SM_m"+to_string(i)+"_module2").c_str())  ;
    MakevsTimePlots(ZeroBias,"ic_ratio_eflow:avg_time", ("n_events>110e6 && ic_ratio_eflow > 0 && ic_ratio_eflow < 100 && ieta>45 && ieta<66 && sm=="+to_string(i)).c_str(), 0.98, 1.01, "IC_{n}/IC_{0}", ("IC_SM_p"+to_string(i)+"_module3").c_str());
    MakevsTimePlots(ZeroBias,"ic_ratio_eflow:avg_time", ("n_events>110e6 && ic_ratio_eflow > 0 && ic_ratio_eflow < 100 && ieta<-45 && ieta>-66 && sm==-"+to_string(i)).c_str(), 0.98, 1.01, "IC_{n}/IC_{0}", ("IC_SM_m"+to_string(i)+"_module3").c_str())  ;
    MakevsTimePlots(ZeroBias,"ic_ratio_eflow:avg_time", ("n_events>110e6 && ic_ratio_eflow > 0 && ic_ratio_eflow < 100 && ieta>65 && ieta<86 && sm=="+to_string(i)).c_str(), 0.98, 1.01, "IC_{n}/IC_{0}", ("IC_SM_p"+to_string(i)+"_module4").c_str());
    MakevsTimePlots(ZeroBias,"ic_ratio_eflow:avg_time", ("n_events>110e6 && ic_ratio_eflow > 0 && ic_ratio_eflow < 100 && ieta<-65 && ieta>-86 && sm==-"+to_string(i)).c_str(), 0.98, 1.01, "IC_{n}/IC_{0}", ("IC_SM_m"+to_string(i)+"_module4").c_str())  ;
  
    }
  

    
    
    
    
    
    
    
    
    
    
    
    
    /*
    TCanvas c;
    float ic_eflow; 

    TH1F* h = new TH1F("","",100, 0.8,1.2); 
    ZeroBias->SetBranchAddress("ic_eflow",&ic_eflow);
    for (int i=0; i<ZeroBias->GetEntries(); i++)
    {
     ZeroBias->GetEntry(i);
     //std::cout << ic_eflow << std::endl;
     h->Fill(ic_eflow);
    
    }
    h->SetTitle(";IC eflow ; events ");

    h->Draw("histo"); 
    c.SaveAs("IC_eflow_v3.png");
    c.SaveAs("IC_eflow_v3.pdf");
    c.SaveAs("IC_eflow_v3.root");
     */
 
     system("mv *.png /eos/user/f/fcetorel/www/PhiSym/eflow");
     system("mv *.pdf /eos/user/f/fcetorel/www/PhiSym/eflow");
   //  system("mv *.root /eos/user/f/fcetorel/www/PhiSym/eflow");


}

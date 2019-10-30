/*
 g++ -Wall -o Eop `root-config --cflags --glibs` -L $ROOTSYS/lib  Eop.cpp
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




void MakeCfrPlots( TString fold, TString cut, TString title)
{
    TChain* ZeroBias;
    ZeroBias = new TChain("eb");
    ZeroBias -> Add("/afs/cern.ch/work/f/fcetorel/private/work2/prova_slc6/CMSSW_9_4_0/src/PhiSym/EcalCalibAlgos/macros/history_eflow_2017.root" );
    ZeroBias -> SetMarkerStyle(20);
    ZeroBias -> SetMarkerColor(kBlue+2);
    
    TChain* EoPi;
    EoPi = new TChain("Eop_monitoring");
    EoPi -> Add("/eos/user/f/fmonti/www/EOverpMonitoring/Run2017ULRereco/history/harness/" + fold );
    EoPi -> SetMarkerStyle(20);
    EoPi -> SetMarkerColor(kRed+2);    
    gStyle->SetOptStat(0);
    TCanvas c;
    c.SetGrid();
    TGraph *gr = new TGraph();
    double t; 
    double Eop_tempfit; 
    double ref; 
    double Eop_norm;
    int n; 
    EoPi->SetBranchAddress("t",&t);
    EoPi->SetBranchAddress("Eop_tempfit",&Eop_tempfit);
    gr -> SetMarkerStyle(20);
    gr -> SetMarkerColor(kRed+2);  
   /* for (int i=0; i< EoPi->GetEntries(); i++)
    {
    EoPi->GetEntry(i);
    std::cout << "t  " << t << std::endl;
    if (t > 1000) 
    {ref = Eop_tempfit; 
    n=i;
    std::cout << "n " << n << std::endl;
    break;
    }
    }*/
    EoPi->GetEntry(1);
     
    ref = Eop_tempfit; 
    //std::cout<< "ref  " << ref << std::endl;
    for (int i=1; i< EoPi->GetEntries(); i++)
    {
   // std::cout << "i " << i << std::endl;
    
    EoPi -> GetEntry(i);
    Eop_norm = Eop_tempfit / ref; 
    //std::cout << "Eop tempfit  " << Eop_tempfit << " / " << "ref " << ref << std::endl;
    //std::cout << "Eop norm  " << Eop_norm << std::endl;
    //std::cout << "t  " << t << std::endl;
    Eop_norm = 1 / Eop_norm; 
    //std::cout << "Entries  " << EoPi->GetEntries() << std::endl;
    if (t > 1000 ) 
    {
    gr->SetPoint(i,t,Eop_norm);
    }
    
    }

    gPad->Update();
    TString cut1 = "& n_events > 100e6";
    ZeroBias->Draw("ic_ratio_eflow_hr:avg_time",cut + cut1,"prof");
    TH2F *gr1 = (TH2F*)gPad->GetPrimitive("htemp");
    gr1->SetTitle("; Time(day/month); ");

    gr1->GetXaxis()->SetTimeDisplay(1);
    gr1->GetXaxis()->SetTimeFormat("%d/%m%F1970-01-01 00:00:00");
   // gr->GetXaxis()->SetRangeUser(1.49698e+09,1511490000);
    gr1->GetYaxis()->SetRangeUser(0.985,1.04);
    gPad->Update();
    gr->Draw("P same");
    
    TLegend leg;
    leg.AddEntry(gr, "1 / EoP", "p");
    leg.AddEntry(gr1, "IC eflow", "p");
    leg.Draw("same");
    c.SaveAs(title + "HRnew.png");
    c.SaveAs(title + "HRnew.pdf");



}

//t>1497932000 && t<1498006000 primo bin di EoP




int main(int argc, char *argv[])
{
     MakeCfrPlots( "IEta_6_25_IPhi_221_230__/IEta_6_25_IPhi_221_230___histos.root", " ieta > 6 && ieta < 25 && iphi > 221 && ieta < 230", "Time_Ev_IEta_6_25_IPhi_221_230");
     MakeCfrPlots( "IEta_66_85_IPhi_221_230__/IEta_66_85_IPhi_221_230___histos.root", " ieta > 66 && ieta < 85 && iphi > 221 && iphi < 230", "Time_Ev_IEta_66_85_IPhi_221_230");
     MakeCfrPlots( "IEta_26_45_IPhi_221_230__/IEta_26_45_IPhi_221_230___histos.root", " ieta > 26 && ieta < 45 && iphi > 221 && iphi < 230", "Time_Ev_IEta_26_45_IPhi_221_230");
     MakeCfrPlots( "IEta_46_65_IPhi_221_230__/IEta_46_65_IPhi_221_230___histos.root", " ieta > 46 && ieta < 65 && iphi > 221 && iphi < 230", "Time_Ev_IEta_46_65_IPhi_221_230");
     MakeCfrPlots( "IEta_6_25_IPhi_11_20__/IEta_6_25_IPhi_11_20___histos.root", " ieta > 6 && ieta < 25 && iphi >11 && iphi < 20", "Time_Ev_IEta_6_25_IPhi_11_20");
     MakeCfrPlots( "IEta_66_85_IPhi_11_20__/IEta_66_85_IPhi_11_20___histos.root", " ieta > 66 && ieta < 85 && iphi > 11 && iphi < 20", "Time_Ev_IEta_66_85_IPhi_11_20");
     MakeCfrPlots( "IEta_26_45_IPhi_11_20__/IEta_26_45_IPhi_11_20___histos.root", " ieta > 26 && ieta < 45 && iphi > 11 && iphi < 20", "Time_Ev_IEta_26_45_IPhi_11_20");
     MakeCfrPlots( "IEta_46_65_IPhi_11_20__/IEta_46_65_IPhi_11_20___histos.root", " ieta > 46 && ieta < 65 && iphi > 11 && iphi < 20", "Time_Ev_IEta_46_65_IPhi_11_20");
     MakeCfrPlots( "IEta_6_25_IPhi_101_110__/IEta_6_25_IPhi_101_110___histos.root", " ieta > 6 && ieta < 25 && iphi >101 && iphi < 110", "Time_Ev_IEta_6_25_IPhi_101_110");
     MakeCfrPlots( "IEta_66_85_IPhi_101_110__/IEta_66_85_IPhi_101_110___histos.root", " ieta > 66 && ieta < 85 && iphi > 101 && iphi < 110", "Time_Ev_IEta_66_85_IPhi_101_110");
     MakeCfrPlots( "IEta_26_45_IPhi_101_110__/IEta_26_45_IPhi_101_110___histos.root", " ieta > 26 && ieta < 45 && iphi > 101 && iphi < 110", "Time_Ev_IEta_26_45_IPhi_101_110");
     MakeCfrPlots( "IEta_46_65_IPhi_101_110__/IEta_46_65_IPhi_101_110___histos.root", " ieta > 46 && ieta < 65 && iphi > 101 && iphi < 110", "Time_Ev_IEta_46_65_IPhi_101_110");

       
    
    
 
     system("mv *.png /eos/user/f/fcetorel/www/PhiSym/eflow/cfr_Eop");
     system("mv *.pdf /eos/user/f/fcetorel/www/PhiSym/eflow/cfr_Eop");


}

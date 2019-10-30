/*
 g++ -Wall -o reweight `root-config --cflags --glibs` -L $ROOTSYS/lib -lRooFitCore -lFoam -lMinuit -lMathMore reweight.cpp
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

int main(int argc, char *argv[])
{
	vector<float> w_etaring;
	TChain* DoubleEG;
	TChain* ZeroBias;
	TString fold_EP = "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/ecalelf/ntuples/13TeV/ALCARERECO/PromptReco2017_pedv1_ps_ICv1_laserv3_LC_Alpha4";
	TString fold_PS = "/afs/cern.ch/work/a/abeschi/public/4Flavia/PhiSYmNtu2017";
TString ntuple_EP[] = {"/DoubleEG-Run2017B-ZSkim-Prompt-v1/297046-297723/pedNoise/DoubleEG-Run2017B-ZSkim-Prompt-v1-297046-297723.root","/DoubleEG-Run2017B-ZSkim-Prompt-v2/298678-299329/pedNoise/DoubleEG-Run2017B-ZSkim-Prompt-v2-298678-299329.root"};
		TString ntuple_PS[] ={"/ntuples_2017B_HarnessCorrection/summed*.root"};
	DoubleEG = new TChain("selected");
	for (int i=0; i<2; i++)
	{
		DoubleEG -> Add(fold_EP + ntuple_EP[i]);
		
	}

	ZeroBias = new TChain("eb_xstals");
	for (int i=0; i<1; i++)
	{
		ZeroBias -> Add(fold_PS + ntuple_PS[i]);
	}
	


	

	Short_t chargeEle[3];
	float etaEle[3];
	
	TH1F* etaEle_histo[2];
	TH1F* etaEle_corr_histo = new TH1F("eta_corr", ";#eta; Counts", 171, -85,85);
	TH1F* weight_histo = new TH1F("weight_histo", ";#eta(Ep/Phi)*norm; Counts", 171, -85,85);

	int ieta = 0.;
	float ring_average = 0.;

	for (int i=0; i<2; i++)
	{
	etaEle_histo[i] = new TH1F(("eta"+ std::to_string(i)).c_str(), ";#eta; Counts", 171, -85, 85);
	}	
	

	DoubleEG-> SetBranchAddress("chargeEle", chargeEle);
	DoubleEG-> SetBranchAddress("etaEle", etaEle);
	ZeroBias-> SetBranchAddress("ieta", &ieta);
	ZeroBias-> SetBranchAddress("ring_average", &ring_average);
	

	int nentries= ZeroBias-> GetEntries();
	for (int i=0; i<nentries; i++)
	{

		ZeroBias->GetEntry(i);
		etaEle_histo[1] -> Fill(ieta, ring_average);

		
		
	}
	nentries= DoubleEG-> GetEntries();
	for (int i=0; i<nentries; i++)
	{
		DoubleEG->GetEntry(i);
		
		if (abs(chargeEle[1])==1 && abs(chargeEle[0])==1 && abs(etaEle[0]) < 1.556)
		{
			etaEle_histo[0] -> Fill(int(etaEle[0]/0.0175));
		}

		
		
	}
	 	




	TCanvas* c1 = new TCanvas();
	c1->cd();
	for (int i=0; i<2; i++) // distribuzioni in eta di Ep e PhiSym 

	{
		etaEle_histo[i] -> Draw("histo");
		c1->SaveAs(("EtaEle"+ std::to_string(i) + ".png").c_str());
		c1->Update();
	}

        // calcolo del peso

 
		float norm_Ep = etaEle_histo[0] ->Integral();
		float norm_phiSym = etaEle_histo[1] ->Integral();
		etaEle_histo[0] -> Divide(etaEle_histo[1]);  
	
        //distribuzione dei pesi
		etaEle_histo[0] -> Draw("histo");
		c1->SaveAs("Pesi_ele.png");
		c1->Update();
		

	for (int i=1; i<172; i++)
	{

	float w=(etaEle_histo[0]->GetBinContent(i) )* (norm_phiSym/norm_Ep);
	etaEle_corr_histo -> SetBinContent(i,w*(etaEle_histo[1]->GetBinContent(i)));
	weight_histo -> SetBinContent(i,w);

	}
		etaEle_corr_histo -> Draw("histo");
		c1->SaveAs("Etaele_corr.png");
		c1->Update();
		weight_histo -> Draw("histo");
		c1->SaveAs("weight.png");
		c1->Update();



	TFile *MyFile = new TFile("weight.root","RECREATE");
	weight_histo -> Write("weight_histo");
	MyFile -> Write();
	
	MyFile -> Close();
		

cout << "# eventi prima delle correzione  " << etaEle_histo[1] -> Integral() ; 
cout << "# eventi dopo delle correzione  " << etaEle_corr_histo -> Integral() << endl; 

system("mv *.png /eos/user/f/fcetorel/www/PhiSym");











}

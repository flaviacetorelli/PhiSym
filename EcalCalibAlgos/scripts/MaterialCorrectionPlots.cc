#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm> 
#include <iostream>


#include "TSystem.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TStyle.h"




using namespace std;


void Make1DPlot(TH1F* h, TString name)
{

	h -> SetLineColor(kAzure);
	h -> SetLineWidth(3);
	h -> SetFillStyle(0);

	TCanvas* c = new TCanvas();
	c -> cd();

	h -> Draw("histo");
	
	c -> SaveAs(name + ".pdf");
	c -> SaveAs(name + ".png");

	return;
}

void Make2DPlot(TH2F* h, TString name, float zRangeLow=0.9, float zRangeHigh=1.1)
{
	TCanvas* c = new TCanvas();
	c -> cd();

	h -> GetZaxis() -> SetRangeUser(zRangeLow, zRangeHigh);
	h -> Draw("COLZ");
	
	c -> SaveAs(name + ".pdf");
	c -> SaveAs(name + ".png");

	return;
}



int main(int argc, char *argv[])
{	
	if(argc<2)
	{
		std::cerr << "No input file provided, please specify the input root file" << endl;
		return -1;
	}

	gStyle -> SetOptStat(0);
	gStyle -> SetPalette(kRainBow);

	TString fileName = argv[1];
	TFile* f = new TFile(fileName, "READ");
	
	TH2F* uncorEBMap = (TH2F*)f->Get("map_ic_uncorr");
	TH2F* corEBMap = (TH2F*)f->Get("map_ic_corr");
	TH2F* EEpMap = (TH2F*)f->Get("map_ic_EEp");
	TH2F* EEmMap = (TH2F*)f->Get("map_ic_EEm");
	TH2F* corrections = (TH2F*)f->Get("map_corrections");
	TH2F* map_lc_EB = (TH2F*)f->Get("map_lc_EB");
	TH2F* map_lc_EEm = (TH2F*)f->Get("map_lc_EEm");
	TH2F* map_lc_EEp = (TH2F*)f->Get("map_lc_EEp");
	TH1F* IC = (TH1F*)f->Get("h_ic_corr");

	Make2DPlot(uncorEBMap, "IC_EB_uncorrected");
	Make2DPlot(corEBMap, "IC_EB_corrected");
	Make2DPlot(corEBMap, "IC_EB_correctedUnzoomed", 0.5, 1.5);
	Make2DPlot(EEpMap, "IC_EEpMap");
	Make2DPlot(EEmMap, "IC_EEmMap");
	Make2DPlot(EEpMap, "IC_EEpMapUnzoomed", 0.5, 1.5);
	Make2DPlot(EEmMap, "IC_EEmMapUnzoomed", 0.5, 1.5);
	Make2DPlot(corrections, "map_corrections");
	Make2DPlot(map_lc_EB, "map_lc_EB", 1, 1.5);
	Make2DPlot(map_lc_EEm, "map_lc_EEm", 1, 5);
	Make2DPlot(map_lc_EEp, "map_lc_EEp", 1, 5);
	Make1DPlot(IC, "ic");

	if(argc<3)
		cout << "No output folder specified, plots will be produced in the execution folder" << endl;
	else
	{
		string outFolder = argv[2];
		cout << "Plots moved in folder " << outFolder << endl;
		system(("mv *.pdf " + outFolder).c_str());
		system(("mv *.png " + outFolder).c_str());
	}

	return 0;
}










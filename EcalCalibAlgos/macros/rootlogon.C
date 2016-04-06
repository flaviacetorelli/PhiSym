{
	gSystem->Load("libFWCoreFWLite.so"); 
	AutoLibraryLoader::enable();
	gSystem->Load("libDataFormatsFWLite.so");
	gSystem->Load("libDataFormatsPatCandidates.so");
        gSystem->Load("libPhiSymEcalCalibDataFormats.so");        

        Int_t MyPalette[200];
        Double_t red[9]   = {  0./255.,   5./255.,  15./255.,  35./255., 102./255., 196./255., 208./255., 199./255., 110./255.};
        Double_t green[9] = {  0./255.,  48./255., 124./255., 192./255., 206./255., 226./255.,  97./255.,  16./255.,   0./255.};
        Double_t blue[9]  = { 99./255., 142./255., 198./255., 201./255.,  90./255.,  22./255.,  13./255.,   8./255.,   2./255.};
        Double_t stops[9] = { 0.0000, 0.2500, 0.3750, 0.4375, 0.5000, 0.5635, 0.6250, 0.7500, 1.0000};
        Int_t FI = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
        for(int i=0;i<200;i++)
            MyPalette[i] = FI+i;        
}

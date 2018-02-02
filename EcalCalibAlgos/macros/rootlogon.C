{
    gSystem->Load("DynamicTTreeDict.so");

    static Int_t  colors[2];
    static Bool_t initialized = kFALSE;
    Double_t Red[2]    = { 1.00, 0.00};
    Double_t Green[2]  = { 0.00, 1.00};
    Double_t Blue[2]   = { 0.00, 0.00};
    Double_t Length[2] = { 0.00, 1.00};
    if(!initialized){
        Int_t FI = TColor::CreateGradientColorTable(2,Length,Red,Green,Blue,2);
        for (int i=0; i<2; i++) colors[i] = FI+i;
        initialized = kTRUE;
        return;
    }
    gStyle->SetPalette(2,colors);
}


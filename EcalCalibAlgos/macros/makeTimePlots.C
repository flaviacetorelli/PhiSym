{
    string type_name;
    int type=0; //IC
    // int type=1; // LC
    // int type=2; // SumEt
    // int type=3; //Nhits/Nhits_tot    
    // int type=4; //k-factors
    // int type=5; //corrections
    bool applyCorr=true;
    
    gSystem->Load("libFWCoreFWLite.so"); 
    AutoLibraryLoader::enable();
    gSystem->Load("libDataFormatsFWLite.so");
    gSystem->Load("libDataFormatsEcalDetId.so");
    gSystem->Load("libPhiSymEcalCalibDataFormats.so");

    vector<string> files={        
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012C_v1/summed_200781_200798.root",
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015A_v2/summed_250866_250866.root"

        "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/newThr_2012D/summed_208538_208686.root",
        "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_newGT_v5/summed_251562_251562.root"

        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015A_v2/summed_248003_248003.root",
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015A_v2/summed_248005_248005.root",
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015A_v2/summed_248006_248006.root",
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015A_v2/summed_248009_248009.root",
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015A_v2/summed_248025_248025.root",
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015A_v2/summed_248026_248026.root",
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015A_v2/summed_248028_248028.root",
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015A_v2/summed_248029_248029.root",
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015A_v2/summed_248030_248030.root",
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015A_v2/summed_248031_248031.root",
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015A_v2/summed_248033_248033.root",
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015A_v2/summed_248036_248036.root",
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015A_v2/summed_248038_248038.root",
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015A_v2/summed_250866_250866.root",
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015A_v2/summed_250867_250867.root",
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015A_v2/summed_250869_250869.root",
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015A_v2/summed_250886_250886.root",
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015A_v2/summed_250890_250890.root",
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015A_v2/summed_250893_250893.root",
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015A_v2/summed_250896_250896.root",
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015A_v2/summed_250897_250897.root",
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015A_v2/summed_250930_250930.root",
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015A_v2/summed_250932_250932.root"
        
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015C_Boff_v1/summed_254292_254292.root",
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015C_Boff_v1/summed_254293_254293.root",
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015C_Boff_v1/summed_254294_254294.root",
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015C_Boff_v1/summed_254306_254306.root",
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015C_Boff_v1/summed_254307_254307.root",
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015C_Boff_v1/summed_254319_254319.root",
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015C_Boff_v1/summed_254341_254341.root",
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015C_Boff_v1/summed_254342_254342.root",
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015C_Boff_v1/summed_254349_254349.root"
        
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_newGT_v5/summed_251244_251244.root",
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_newGT_v5/summed_251251_251251.root",
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_newGT_v5/summed_251252_251252.root",
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_newGT_v5/summed_251521_251521.root",
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_newGT_v5/summed_251522_251522.root",
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_newGT_v5/summed_251548_251548.root",
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_newGT_v5/summed_251559_251559.root",
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_newGT_v5/summed_251560_251560.root",
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_newGT_v5/summed_251561_251561.root",
        // "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_newGT_v5/summed_251562_251562.root"
    };

    vector<string> corrections_files={
        //---Boff
        // "/afs/cern.ch/user/s/spigazzi/work/ECAL/CMSSW_7_4_6_patch6/src/PhiSym/EcalCalibAlgos/ntuples/2012C_v1/geo_and_material_corr.txt",
        // "/afs/cern.ch/user/s/spigazzi/work/ECAL/CMSSW_7_4_6_patch6/src/PhiSym/EcalCalibAlgos/ntuples/2015A_v2/geo_and_material_corr.txt"

        //---Bon
        "/afs/cern.ch/user/s/spigazzi/work/ECAL/CMSSW_7_4_6_patch6/src/PhiSym/EcalCalibAlgos/ntuples/newThr_2012D/corrections_208538_208686.txt",
        "/afs/cern.ch/user/s/spigazzi/work/ECAL/CMSSW_7_4_6_patch6/src/PhiSym/EcalCalibAlgos/ntuples/2015B_newGT_v5/corrections_251562_251562.txt"
    };

    float ebVar[70][EBDetId::kSizeForDenseIndexing]={};
    float eeVar[70][EEDetId::kSizeForDenseIndexing]={};
    float ebCorr[70][EBDetId::kSizeForDenseIndexing]={};
    pair<int, int> ebMap[EBDetId::kSizeForDenseIndexing];
    int eeMap[EEDetId::kSizeForDenseIndexing];
    pair<int, int> eeMap2[EEDetId::kSizeForDenseIndexing];

    for(int iFile=0; iFile<files.size(); ++iFile)
    {
        TFile* file = TFile::Open(files[iFile].c_str(), "READ");
        if(!file)
            continue;
        CrystalsEBTree ebTree((TTree*)file->Get("eb_xstals"));
        CrystalsEETree eeTree((TTree*)file->Get("ee_xstals"));
        
        int index=-1;
        long int tot_hits_EB=0;
        long int tot_hits_EE=0;
        while(ebTree.NextEntry())
        {
            if(ebTree.rec_hit->GetNhits() > 0)
                tot_hits_EB += ebTree.rec_hit->GetNhits();
            index=EBDetId(ebTree.ieta, ebTree.iphi).hashedIndex();;
            switch(type)
            {
            case 0:
                ebVar[iFile][index] = ebTree.ic_abs*ebTree.ic_ch;
                type_name = "IC";
                break;
            case 1:                
                ebVar[iFile][index] = ebTree.rec_hit->GetLCSum()/ebTree.rec_hit->GetNhits();
                type_name = "LC";
                break;
            case 2:
                ebVar[iFile][index] = ebTree.rec_hit->GetSumEt(0)/ebTree.rec_hit->GetNhits();
                type_name = "SumEt";
                break;
            case 3:
                ebVar[iFile][index] = ebTree.rec_hit->GetNhits();
                type_name = "Nhits";
                break;
            case 4:
                ebVar[iFile][index] = ebTree.k_ch;
                type_name = "Kfact";
                break;
            case 5:
                ebVar[iFile][index] = 1;
                type_name = "Corr";
                break;
            }
            if(iFile==0)
                ebMap[index] = make_pair(ebTree.ieta, ebTree.iphi);
        }

        while(eeTree.NextEntry())
        {
            if(eeTree.rec_hit->GetNhits() > 0)
                tot_hits_EE += eeTree.rec_hit->GetNhits();
            index=EEDetId(eeTree.ix, eeTree.iy, eeTree.iring>0?1:-1).hashedIndex();
            switch(type)
            {
            case 0:
                eeVar[iFile][index] = eeTree.ic_ch;
                type_name = "IC";
                break;
            case 1:                
                eeVar[iFile][index] = eeTree.rec_hit->GetLCSum()/eeTree.rec_hit->GetNhits();;
                type_name = "LC";
                break;
            case 2:
                eeVar[iFile][index] = eeTree.rec_hit->GetSumEt(0)/eeTree.rec_hit->GetNhits();
                type_name = "SumEt";
                break;                
            case 3:
                eeVar[iFile][index] = eeTree.rec_hit->GetNhits();
                type_name = "Nhits";
                break;
            case 4:
                eeVar[iFile][index] = eeTree.k_ch;
                type_name = "Kfact";
                break;
            case 5:
                eeVar[iFile][index] = 1;
                type_name = "Corr";
                break;
            }
            if(iFile==0)
            {
                eeMap[index] = eeTree.iring<0? eeTree.iring+39 : eeTree.iring+38;
                eeMap2[index] = make_pair(eeTree.iring<0 ? -eeTree.ix : eeTree.ix, eeTree.iy);
            }
        }
        if(type == 3)
        {
            float sum=0;
            for(int index=0; index<EBDetId::kSizeForDenseIndexing; ++index)
                sum += ebVar[iFile][index];
            for(int index=0; index<EBDetId::kSizeForDenseIndexing; ++index)
                ebVar[iFile][index] = ebVar[iFile][index]/sum;
            sum=0;
            for(int index=0; index<EEDetId::kSizeForDenseIndexing; ++index)
                sum += eeVar[iFile][index];
            for(int index=0; index<EEDetId::kSizeForDenseIndexing; ++index)
                eeVar[iFile][index] = eeVar[iFile][index]/sum;
        }
        cout << files[iFile] << " nhits/crystal (EB/EE): " << tot_hits_EB/71200. << "  " << tot_hits_EE/14000. << endl;
        file->Close();
        
        //---get geo&material correction if type == ic
        if((type == 0 && applyCorr) || type==5)
        {
            int ieta, iphi, side;
            float corr;
            
            ifstream corrections(corrections_files[iFile], ios::in);
            while(corrections.good())
            {
                corrections >> ieta >> iphi >> side >> corr;
                ebVar[iFile][EBDetId(ieta, iphi).hashedIndex()] *= corr;
            }
            corrections.close();
        }
    }

    //---compute & draw
    TFile* outFile = new TFile(string(type_name+"_ratios_vs_time.root").c_str(), "RECREATE");
    outFile->cd();
    //---EB
    TH1F* hAbsEB = new TH1F("hAbsEB", "Absolute #sigma_{ratio} spread --- barrel", 250, 0, 0.1);
    TH1F* hRelEB = new TH1F("hRelEB", "Relative #sigma_{ratio} spread --- barrel", 250, 0, 0.1);
    TH1F* prAbsEB[172];
    TH1F* prRelEB[172];
    TGraphErrors* grAbsRingEB[172];
    TGraphErrors* grRelRingEB[172];
    TGraphErrors* grAbsEB = new TGraphErrors();
    TGraphErrors* grRelEB = new TGraphErrors();
    TF1* fitFuncAbsEB = new TF1("fitFuncAbsEB", "gaus", 0.9, 1.1);
    TF1* fitFuncRelEB = new TF1("fitFuncRelEB", "gaus", 0.98, 1.02);
    //---EE
    TH1F* hAbsEE = new TH1F("hAbsEE", "Absolute #sigma_{ratio} spread --- endcaps", 250, 0, 0.5);
    TH1F* hRelEE = new TH1F("hRelEE", "Relative #sigma_{ratio} spread --- endcaps", 250, 0, 0.5);
    TH1F* prAbsEE[80];
    TH1F* prRelEE[80];
    TGraphErrors* grAbsRingEE[80];
    TGraphErrors* grRelRingEE[80];
    TGraphErrors* grAbsEE = new TGraphErrors();
    TGraphErrors* grRelEE = new TGraphErrors();
    TF1* fitFuncAbsEE = new TF1("fitFuncAbsEE", "gaus", 0.8, 1.2);
    TF1* fitFuncRelEE = new TF1("fitFuncRelEE", "gaus", 0.8, 1.2);

    float mapAbsEB_range[2];
    float mapRelEB_range[2];
    float mapAbsEE_range[2];
    float mapRelEE_range[2];
    
    //---init eta histos
    for(int iRing=0; iRing<172; ++iRing)
    {
        prAbsEB[iRing] = new TH1F(to_string(iRing).c_str(), "", 1000, 0.8, 1.2);
        prRelEB[iRing] = new TH1F(to_string(iRing+200).c_str(), "", 1000, 0.9, 1.1);
        grAbsRingEB[iRing] = new TGraphErrors();
        grRelRingEB[iRing] = new TGraphErrors();
    }
    for(int iRing=0; iRing<80; ++iRing)
    {
        prAbsEE[iRing] = new TH1F(to_string(iRing+400).c_str(), "", 1000, 0.7, 1.3);
        prRelEE[iRing] = new TH1F(to_string(iRing+600).c_str(), "", 1000, 0.7, 1.3);
        grAbsRingEE[iRing] = new TGraphErrors();
        grRelRingEE[iRing] = new TGraphErrors();
    }
    //---fill histos
    for(int iFile=1; iFile<files.size(); ++iFile)
    {
        TH1F* tmpAbsEB = new TH1F("tmpAbsEB", "", 2000, 0.5, 1.5);
        TH1F* tmpRelEB = new TH1F("tmpRelEB", "", 2000, 0.5, 1.5);
        TH2F* mapAbsEB = new TH2F("mapAbsEB", "", 360, 0.5, 360.5, 171, -85, 85);        
        TH2F* mapRelEB = new TH2F("mapRelEB", "", 360, 0.5, 360.5, 171, -85, 85);
        TH2F* mapAbsEE = new TH2F("mapAbsEE", "", 201, -100.5, 100.5, 101, 0.5, 100.5);
        TH2F* mapRelEE = new TH2F("mapRelEE", "", 201, -100.5, 100.5, 101, 0.5, 100.5);
        mapAbsEB->SetContour(100);
        mapRelEB->SetContour(100);
        mapAbsEE->SetContour(100);
        mapRelEE->SetContour(100);
        //---EB
        for(int index=0; index<EBDetId::kSizeForDenseIndexing; ++index)
        {
            if(ebVar[0][index] != 0 && ebVar[iFile][index] != 0)
            {
                tmpAbsEB->Fill(ebVar[iFile][index]/ebVar[0][index]);
                mapAbsEB->Fill(ebMap[index].second, ebMap[index].first, ebVar[iFile][index]/ebVar[0][index]);
            }
            if(ebVar[iFile-1][index] != 0 && ebVar[iFile][index] != 0)
            {
                tmpRelEB->Fill(ebVar[iFile][index]/ebVar[iFile-1][index]);
                mapRelEB->Fill(ebMap[index].second, ebMap[index].first, ebVar[iFile][index]/ebVar[iFile-1][index]);
            }
        }
        //---EE
        TH1F* tmpAbsEE = new TH1F("tmpAbsEE", "", 1000, 0.5, 1.5);
        TH1F* tmpRelEE = new TH1F("tmpRelEE", "", 1000, 0.5, 1.5);
        for(int index=0; index<EEDetId::kSizeForDenseIndexing; ++index)
        {
            if(eeVar[0][index] != 0 && eeVar[iFile][index] != 0)
            {
                tmpAbsEE->Fill(eeVar[iFile][index]/eeVar[0][index]);
                mapAbsEE->Fill(eeMap2[index].first, eeMap2[index].second, eeVar[iFile][index]/eeVar[0][index]);
                prAbsEE[eeMap[index]]->Fill(eeVar[iFile][index]/eeVar[0][index]);
            }
            if(eeVar[iFile-1][index] != 0 && eeVar[iFile][index] != 0)
            {
                tmpRelEE->Fill(eeVar[iFile][index]/eeVar[iFile-1][index]);
                mapRelEE->Fill(eeMap2[index].first, eeMap2[index].second, eeVar[iFile][index]/eeVar[iFile-1][index]);
                prRelEE[eeMap[index]]->Fill(eeVar[iFile][index]/eeVar[iFile-1][index]);
            }       
        }
        // if(iFile==1)
        // {
        mapAbsEB_range[0] = tmpAbsEB->GetMean()-2*tmpAbsEB->GetRMS();
        mapAbsEB_range[1] = tmpAbsEB->GetMean()+2*tmpAbsEB->GetRMS();
        mapRelEB_range[0] = tmpRelEB->GetMean()-2*tmpRelEB->GetRMS();
        mapRelEB_range[1] = tmpRelEB->GetMean()+2*tmpRelEB->GetRMS();
        mapAbsEE_range[0] = tmpAbsEE->GetMean()-2*tmpAbsEE->GetRMS();
        mapAbsEE_range[1] = tmpAbsEE->GetMean()+2*tmpAbsEE->GetRMS();
        mapRelEE_range[0] = tmpRelEE->GetMean()-2*tmpRelEE->GetRMS();
        mapRelEE_range[1] = tmpRelEE->GetMean()+2*tmpRelEE->GetRMS();
            //}
        
        //---EB                     
        hAbsEB->Fill(fitFuncAbsEB->GetParameter(2));
        grAbsEB->SetPoint(iFile-1, iFile, tmpAbsEB->GetRMS());
        grAbsEB->SetPointError(iFile-1, 0, tmpAbsEB->GetRMSError());
        hRelEB->Fill(fitFuncRelEB->GetParameter(2));
        grRelEB->SetPoint(iFile-1, iFile, tmpRelEB->GetRMS());
        grRelEB->SetPointError(iFile-1, 0, tmpRelEB->GetRMSError());
        tmpAbsEB->Write(string("AbsEB_"+to_string(iFile)).c_str());
        tmpRelEB->Write(string("RelEB_"+to_string(iFile)).c_str());        
        mapAbsEB->SetAxisRange(mapAbsEB_range[0], mapAbsEB_range[1], "Z");
        mapRelEB->SetAxisRange(mapRelEB_range[0], mapRelEB_range[1], "Z");
        mapAbsEB->Write(string("mapAbsEB_"+to_string(iFile)).c_str());
        mapRelEB->Write(string("mapRelEB_"+to_string(iFile)).c_str());
        tmpAbsEB->Delete();
        tmpRelEB->Delete();
        mapAbsEB->Delete();
        mapRelEB->Delete();
        //---EE
        hAbsEE->Fill(fitFuncAbsEE->GetParameter(2));
        grAbsEE->SetPoint(iFile-1, iFile, tmpAbsEE->GetRMS());
        grAbsEE->SetPointError(iFile-1, 0, tmpAbsEE->GetRMSError());
        hRelEE->Fill(fitFuncRelEE->GetParameter(2));
        grRelEE->SetPoint(iFile-1, iFile, tmpRelEE->GetRMS());
        grRelEE->SetPointError(iFile-1, 0, tmpRelEE->GetRMSError());
        tmpAbsEE->Write(string("AbsEE_"+to_string(iFile)).c_str());
        tmpRelEE->Write(string("RelEE_"+to_string(iFile)).c_str()); 
        mapAbsEE->SetAxisRange(mapAbsEE_range[0], mapAbsEE_range[1], "Z");
        mapRelEE->SetAxisRange(mapRelEE_range[0], mapRelEE_range[1], "Z");
        mapAbsEE->Write(string("mapAbsEE_"+to_string(iFile)).c_str());
        mapRelEE->Write(string("mapRelEE_"+to_string(iFile)).c_str());
        tmpAbsEE->Delete();
        tmpRelEE->Delete();
        mapAbsEE->Delete();
        mapRelEE->Delete();

        //---reset
        for(int iRing=0; iRing<171; ++iRing)
        {
            grAbsRingEB[iRing]->SetPoint(iFile-1, iFile, prAbsEB[iRing]->GetRMS());
            grAbsRingEB[iRing]->SetPointError(iFile-1, 0, prAbsEB[iRing]->GetRMSError());
            grRelRingEB[iRing]->SetPoint(iFile-1, iFile, prRelEB[iRing]->GetRMS());
            grRelRingEB[iRing]->SetPointError(iFile-1, 0, prRelEB[iRing]->GetRMSError());
        }
        for(int iRing=0; iRing<79; ++iRing)
        {
            grAbsRingEE[iRing]->SetPoint(iFile-1, iFile, prAbsEE[iRing]->GetRMS());
            grAbsRingEE[iRing]->SetPointError(iFile-1, 0, prAbsEE[iRing]->GetRMSError());
            grRelRingEE[iRing]->SetPoint(iFile-1, iFile, prRelEE[iRing]->GetRMS());
            grRelRingEE[iRing]->SetPointError(iFile-1, 0, prRelEE[iRing]->GetRMSError());
            prAbsEE[iRing]->Reset();
            prRelEE[iRing]->Reset();
        }
    }
    //---EB
    hAbsEB->SetFillColor(kBlue-4);
    hRelEB->SetFillColor(kRed-4);
    hAbsEB->Write("hAbsEB");
    hRelEB->Write("hRelEB");
    grAbsEB->SetMarkerColor(kBlue);
    grAbsEB->SetMarkerStyle(20);
    grAbsEB->SetMarkerSize(0.7);
    grRelEB->SetMarkerColor(kRed);
    grRelEB->SetMarkerStyle(20);
    grRelEB->SetMarkerSize(0.7);
    grAbsEB->Write("grAbsEB");
    grRelEB->Write("grRelEB");
    // mapAbsEB->Write("mapAbsEB");
    // mapRelEB->Write("mapRelEB");
    //---EE
    hAbsEE->SetFillColor(kBlue-4);
    hRelEE->SetFillColor(kRed-4);
    hAbsEE->Write("hAbsEE");
    hRelEE->Write("hRelEE");
    grAbsEE->SetMarkerColor(kBlue);
    grAbsEE->SetMarkerStyle(20);
    grAbsEE->SetMarkerSize(0.7);
    grRelEE->SetMarkerColor(kRed);
    grRelEE->SetMarkerStyle(20);
    grRelEE->SetMarkerSize(0.7);
    grAbsEE->Write("grAbsEE");
    grRelEE->Write("grRelEE");
    // mapAbsEE->Write("mapAbsEE");
    // mapRelEE->Write("mapRelEE");

    grAbsRingEB[5]->Write("grAbs_ieta_80");
    grAbsRingEB[80]->Write("grAbs_ieta_5");
    grAbsRingEB[5]->Write("grAbs_ieta_165");
    grAbsRingEB[80]->Write("grAbs_ieta_90");        
    grRelRingEB[5]->Write("grRel_ieta_80");
    grRelRingEB[80]->Write("grRel_ieta_5");
    grRelRingEB[5]->Write("grRel_ieta_165");
    grRelRingEB[80]->Write("grRel_ieta_90");

    grAbsRingEE[4]->Write("grAbs_iring_5m");
    grAbsRingEE[14]->Write("grAbs_iring_15m");
    grAbsRingEE[24]->Write("grAbs_iring_25m");
    grAbsRingEE[34]->Write("grAbs_iring_35m");        
    grAbsRingEE[4+39]->Write("grAbs_iring_5p");
    grAbsRingEE[14+39]->Write("grAbs_iring_15p");
    grAbsRingEE[24+38]->Write("grAbs_iring_25p");
    grAbsRingEE[34+39]->Write("grAbs_iring_35p");
    grRelRingEE[4]->Write("grRel_iring_5m");
    grRelRingEE[14]->Write("grRel_iring_15m");
    grRelRingEE[24]->Write("grRel_iring_25m");
    grRelRingEE[34]->Write("grRel_iring_35m");        
    grRelRingEE[4+39]->Write("grRel_iring_5p");
    grRelRingEE[14+39]->Write("grRel_iring_15p");
    grRelRingEE[24+38]->Write("grRel_iring_25p");
    grRelRingEE[34+39]->Write("grRel_iring_35p");
    
    outFile->Close();
}

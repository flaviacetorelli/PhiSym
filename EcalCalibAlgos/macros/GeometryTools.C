vector<int> PNList()
{
    vector<int> pn_regions;
    for(int ism=1; ism<=18; ++ism)
    {
        for(int ipn=0; ipn<9; ++ipn)
        {
            pn_regions.push_back(ism*9+ipn);
            pn_regions.push_back(-(ism*9+ipn));
        }
    }
    return pn_regions;
}

int WhichPN(int ieta, int iphi)
{
    int sm = ieta>0 ? (iphi-1)/20+1 : -(iphi-1)/20-1;
    int local_pn=0;
    if(abs(ieta) > 5 && abs(ieta)<26)
        local_pn = (iphi-1) % 20 >= 10 ? 2 : 1;
    else if(abs(ieta) > 25 && abs(ieta)<46)
        local_pn = (iphi-1) % 20 >= 10 ? 4 : 3;
    else if(abs(ieta) > 45 && abs(ieta)<66)
        local_pn = (iphi-1) % 20 >= 10 ? 6 : 5;    
    else if(abs(ieta) > 65)
        local_pn = (iphi-1) % 20 >= 10 ? 8 : 7;
    
    int pn = sm/(abs(sm))*(abs(sm)*9 + local_pn);

    return pn;
}

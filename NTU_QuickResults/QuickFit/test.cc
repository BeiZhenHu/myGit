

void test(){
    TH1::AddDirectory(kFALSE);

    TH1F *h_Ped_H[2][130];
    TH1F *h_Phy_H[2][130];
    TH1F *h_Ped_L[2][130];
    TH1F *h_Phy_L[2][130];
    TFile fPed( TString::Format("dryrun_42_GU8162.root") );
    
    for (int iChn = 0; iChn<128; iChn++) {
        h_Ped_H[0][iChn] = (TH1F*)fPed.Get(TString::Format("/HighGain/Ping_Chn%1d_HG",iChn));
        h_Ped_L[0][iChn] = (TH1F*)fPed.Get(TString::Format("/LowGain/Ping_Chn%1d_LG",iChn));
        
    }
    
}

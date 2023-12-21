// Fit Histogram v3
//
// by Bei-Zhen Hu 2020.09.14
//
// USAGE: ./main datapath runName
// [datapath]: data path must wiht "/" in the end
//             ex: ~/path/to/data/
// 2021.1.27 for IHEP
// 2022.9.4 for onsite test
// 2022.9.6 onsite quick results
// 2022.9.18 add output files pdf

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <iomanip>
#include "TH1F.h"
#include "TF1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TMath.h"
#include "TLegend.h"

using namespace std;
float unitPico = 1.0*1e-12;
float uniQ = 1.6*1e-19;

struct CalibBrief {
    float hPing_P0, hPing_P0_Err;
    float hPing_P1, hPing_P1_Err;
    float hPong_P0, hPong_P0_Err;
    float hPong_P1, hPong_P1_Err;
    float lPing_P1, lPong_P1;
    
};

struct ChnMap {
    
    int AxonConnChnID; // [0]
    int AxonConnPosID; // [1]
    string HVS_ChanGp; // [2]

    int SpltChnID;   //[3]
    int SamTecPinNo; //[4]
    int GCU_HV_Gp;   //[5]
    int asicID;      //[6]
    int asicChnID;   //[7]
    int ABC_ChnID;   //[8]
    
};

struct PmtTable {
    int SerialNo;
    int PMT_ID;
    int HvLv;
    int CableL;
    
    string Hzc_ConnID;
    string Hzc_ConnID1;
    string Hzc_PlugID;
    
};

struct UwbTable {
    int ConnID;  // connector position on the UWB
    int gcuID;
    int abcID;
    int CableL;
    float gxuHV;
    
    string UWB_PID;  //UWB position
    string Axion_SN; //Zxxx
    
};

struct UwbGcuAbcTable {
    int gcuID;
    int abcID;
    
};



map<int /* ChnN */, vector< CalibBrief > > GainInfo;
map<int /* AsicNo */, vector< ChnMap > > ChannelMap;
map<string /* ConNo */, vector< PmtTable > > PottingMap;
map<int /* GCU No */, vector< UwbTable > > UwbMapping;
map<int /* GCU No */, vector< UwbGcuAbcTable > > GcuAbcMapping;

TH1F *h_Ped_H[2][130];
TH1F *h_Phy_H[2][130];
TH1F *h_Ped_L[2][130];
TH1F *h_Phy_L[2][130];
TH1F *h_PedMean[2];
TH1F *h_PhyMean[2];
TH1F *h_PedSigma[2];
TH1F *h_PhySigma[2];
TH1F *h_QMean[2];
TH1F *h_QSigma[2];
TH1F *h_Gain[2];
TH1F *h_EvtNum[2];
TH1F *h_EvtNumSig[2];
TH1F *h_calbMap_p0[2];
TH1F *h_calbMap_p1[2];
TH1F *h_runT[2];
TH1F *h_ref[2];
TH1F *h_compare[2];
TH1F *h_dkRate;
TH1F *h_dkRef;
TH1F *h_dkCompare;
TLegend* legPing;
TLegend* legPong;
TLegend* legDK;
int BeginJob();
int ResetHist();

int main(int argc,char** argv) {
    //gStyle->SetOptStat(kFALSE);
    gStyle->SetOptFit(kTRUE);
    gStyle->SetPadRightMargin(0.06);
    gStyle->SetPadLeftMargin(0.13);

    if( argc != 3 ) {
        
        std::cerr << "# of arguments error." << std::endl;
        std::cerr << "USAGE: ./main datapath runName" << std::endl;
        return 1;
    }
    
    string rootPath = argv[1];
    string runName   = argv[2];
    int runNo = atoi( runName.substr(7,2).c_str() );
    int gcu_ID = atoi( runName.substr(12,4).c_str() );
   

    cout<<"RunNO: "<<runName<<endl;
    cout<<"GCU: "<<gcu_ID<<endl;
    //cout<<"Physics Run: "<<physRun<<endl;

    ////////// UWB Connector Mapping ///////////////
    ///
    ifstream Data5;
    Data5.open(TString::Format("../Tables/spmtUWBconInfo_0617.txt"));
    string inbuf5;
    getline(Data5, inbuf5, '\n');
    cout<<" test5: "<<inbuf5<<endl;
    string inputTable5[11];

    for (int i=0; i<208; i++) {
       
        for (int iColumn=0; iColumn<11; iColumn++){
            Data5>>inputTable5[iColumn];

        }
        cout<<"55 TEST--- "<<i<<"\t"<<inputTable5[7]<<" "<<stoi(inputTable5[4])<<endl;
        int GcuNo = atoi( inputTable5[3].substr(2,4).c_str() );
        cout<<"GCU: "<<GcuNo<<endl;
        UwbTable uwbMapping;
        uwbMapping.ConnID = stoi(inputTable5[2]);
        uwbMapping.abcID  = stoi(inputTable5[4]);
        uwbMapping.CableL = stoi(inputTable5[6]);
        uwbMapping.gxuHV  = stof(inputTable5[8]);
        
        uwbMapping.UWB_PID  = inputTable5[1];
        uwbMapping.Axion_SN = inputTable5[7];

        UwbMapping[ GcuNo ].push_back( uwbMapping );
        
    }
    Data5.close();
    
    ////////// GCU ABC Mapping ///////////////
    ///
    ifstream Data6;
    Data6.open(TString::Format("../Tables/uwbGcuAbcMap.txt"));
    string inbuf6;
    getline(Data6, inbuf6, '\n');
    cout<<" test6: "<<inbuf6<<endl;
    string inputTable6[4];

    for (int i=0; i<206; i++) {
       
        for (int iColumn=0; iColumn<4; iColumn++){
            Data6>>inputTable6[iColumn];

        }
        cout<<"55 TEST--- "<<i<<"\t"<<inputTable6[3]<<" "<<stoi(inputTable6[2])<<endl;
        int GcuNo = atoi( inputTable6[3].substr(2,4).c_str() );
        cout<<"GCU: "<<GcuNo<<endl;
        UwbGcuAbcTable gcuAbcMapping;
        gcuAbcMapping.abcID  = stoi(inputTable6[2]);
        GcuAbcMapping[ GcuNo ].push_back( gcuAbcMapping );
        
    }
    Data6.close();
    
    int abc_ID = GcuAbcMapping[ gcu_ID ][0].abcID;
    cout<<"ABC: "<<abc_ID<<endl;
    /////////
    
    ifstream Data1;
    Data1.open( TString::Format("../Tables/database_all/ABC_database_%1d/calib_info.txt",abc_ID) );
    cout<<"Check ABC_no: "<<TString::Format("../Tables/database_all/ABC_database_%1d/calib_info.txt",abc_ID)<<endl;
    string inbuf;
    getline(Data1, inbuf, '\n');
    getline(Data1, inbuf, '\n');
    
    double inputTable[15];
    while (!Data1.eof()) {
        for (int iColumn=0; iColumn<15; iColumn++)Data1>>inputTable[iColumn];
        int   ChnNb = inputTable[0];
            
        CalibBrief gRatio;
        gRatio.hPing_P1     = inputTable[1];
        gRatio.hPing_P1_Err = inputTable[2];
        gRatio.hPing_P0     = inputTable[3];
        gRatio.hPing_P0_Err = inputTable[4];
            
        gRatio.hPong_P1     = inputTable[6];
        gRatio.hPong_P1_Err = inputTable[7];
        gRatio.hPong_P0     = inputTable[8];
        gRatio.hPong_P0_Err = inputTable[9];
            
        gRatio.lPing_P1     = inputTable[12];
        gRatio.lPong_P1     = inputTable[14];
            
        GainInfo[ ChnNb ].push_back( gRatio );

    }
    Data1.close();
    
    //////
    ////////// A Connector Mapping ///////////////
    /// Cable Mapping on ABC + splitter and AXON connectors
    
    ifstream Data2;
    Data2.open(TString::Format("../Tables/Onsite_CableMapping.txt"));
    string inbuf2;
    getline(Data2, inbuf2, '\n');
    cout<<" test2: "<<inbuf2<<endl;
    
    //string inputTable2[6];
    string inputTable2[9];
    for (int i=0; i<128; i++) {
        for (int iColumn=0; iColumn<9; iColumn++)Data2>>inputTable2[iColumn];
        cout<<"TEST--- "<<inputTable2[7]<<" "<<inputTable2[1]<<endl;
        //string ConNo =
        int ConnPosition = stoi(inputTable2[1]);
        
        ChnMap AsicConnMapping;

        AsicConnMapping.AxonConnChnID = stoi(inputTable2[0]);
       
        AsicConnMapping.HVS_ChanGp    = inputTable2[2];
        AsicConnMapping.SpltChnID     = stoi(inputTable2[3]);
        AsicConnMapping.SamTecPinNo   = stoi(inputTable2[4]);
        AsicConnMapping.GCU_HV_Gp     = stoi(inputTable2[5]);
        AsicConnMapping.asicID        = stoi(inputTable2[6]);
        AsicConnMapping.asicChnID     = stoi(inputTable2[7]);
        AsicConnMapping.ABC_ChnID     = stoi(inputTable2[8]);
        ChannelMap[ ConnPosition ].push_back( AsicConnMapping ); //ChannelMap[ConnPosition][ConnectorID]
        
    }
    Data2.close();
    //////////
    
    ////////// PMT Connector Mapping ///////////////
    /// mapping for pmt potting
    ifstream Data3;
    //Data3.open(TString::Format("../Tables/Juno_spmtMap14_20220725.txt"));
    Data3.open(TString::Format("../Tables/JunoAndNTU_spmtMap.txt"));
    
    string inbuf3;
    getline(Data3, inbuf3, '\n');
    cout<<" test3: "<<inbuf3<<endl;
    string inputTable3[9];

    //for (int i=0; i<25936; i++) {
    for (int i=0; i<26064; i++) {
        for (int iColumn=0; iColumn<8; iColumn++)Data3>>inputTable3[iColumn];
        cout<<"TEST--- "<<inputTable3[7]<<" "<<inputTable3[1]<<endl;
      
        string ConNo = inputTable3[1];
        PmtTable pmtMapping;
        pmtMapping.PMT_ID = stoi(inputTable3[2]);
        pmtMapping.SerialNo = stoi(inputTable3[3]);
        int Lth = stoi(inputTable3[0].substr(5,3).c_str());
        if (Lth==91) pmtMapping.CableL= 5;
        if (Lth==92) pmtMapping.CableL= 10;
        pmtMapping.HvLv= stoi( inputTable3[5] );

        PottingMap[ ConNo ].push_back( pmtMapping );
        
    }
    Data3.close();
    
    ////////// DarkRoom PMT ///////////////
    int howManyRuns = 1;
    ifstream Data4;
    Data4.open(TString::Format("../Tables/darkRoomPmts.txt"));
    string inbuf4;
    getline(Data4, inbuf4, '\n');
    cout<<" test4: "<<inbuf4<<endl;
    
    
    //////
    int PmtConnPos[howManyRuns][8];
    string PmtConName[howManyRuns][8];
    int shvHV[howManyRuns][8];
    
    for (int k=0; k<howManyRuns; k++) {
        string RunInfo[28];
        for (int iColumn=0; iColumn<24; iColumn++) Data4>>RunInfo[iColumn];
    
        for (int j=0; j<8; j++) {
            PmtConnPos[k][j] = stoi(RunInfo[j]);
            PmtConName[k][j] = RunInfo[8+j];
            shvHV[k][j]      = stoi(RunInfo[16+j]);
            
            cout<<"CHECK: "<<"\t"<<PmtConnPos[k][j]<<"\t"<<PmtConName[k][j]<<"\t"<<shvHV[k][j]<<endl;
        }
        
    }
    Data3.close();
    
    /////////////
    ofstream opt1;
    opt1.open( TString::Format("quickTB_%1d.txt", abc_ID) );
    //opt1<<"runNo\trunType\tMode\tSN\tPlugID\twConnID\tPmtID\tMn_ADC\tMnErr_ADC\tSg_ADC\tSgErr_ADC\tMn_Q\tMnErr_Q\tSg_Q\tSgErr_Q\tappHV"<<endl;
    
    
    ///////
    BeginJob();
    

    TFile *ofile = new TFile(TString::Format("quickLook_%s.root", runName.c_str() ) ,"recreate");
    TFile fPed( TString::Format("%s%s.root",rootPath.c_str(), runName.c_str()) );
    
    
    ofile->cd();
    ofile->mkdir("PedestalRun");

    
    for (int iSign = 0; iSign<=1; iSign++) { // ping:0, pong:1

        float calbP0 = 0;
        float calbP1 = 0;
        float calbP0err = 0;
        float calbP1err = 0;
        
        for (int iChn = 0; iChn<128; iChn++) {
            //if (iChn==60) continue;
            int nEvtH_pi, nEvtH_po, nEvtL_pi, nEvtL_po;
            if (iSign==0) {
                
                calbP0 = GainInfo[iChn][0].hPing_P0;
                calbP0err = GainInfo[iChn][0].hPing_P0_Err;
                
                calbP1 = GainInfo[iChn][0].hPing_P1;
                calbP1err = GainInfo[iChn][0].hPing_P1_Err;
                
                h_Ped_H[iSign][iChn] = (TH1F*)fPed.Get(TString::Format("/HighGain/Ping_Chn%1d_HG",iChn));
                h_Ped_L[iSign][iChn] = (TH1F*)fPed.Get(TString::Format("/LowGain/Ping_Chn%1d_LG",iChn));
                
            }
            if(iSign==1) {
                
                calbP0 = GainInfo[iChn][0].hPong_P0;
                calbP0err = GainInfo[iChn][0].hPong_P0_Err;
                
                calbP1 = GainInfo[iChn][0].hPong_P1;
                calbP1err = GainInfo[iChn][0].hPong_P1_Err;
                
                h_Ped_H[iSign][iChn] = (TH1F*)fPed.Get(TString::Format("/HighGain/Pong_Chn%1d_HG",iChn));
                h_Ped_L[iSign][iChn] = (TH1F*)fPed.Get(TString::Format("/LowGain/Pong_Chn%1d_LG",iChn));
                
            }
            h_calbMap_p0[iSign]->SetBinContent(iChn+1,calbP0);
            h_calbMap_p0[iSign]->SetBinError(iChn+1,calbP0err);
            h_calbMap_p1[iSign]->SetBinContent(iChn+1,calbP1);
            h_calbMap_p1[iSign]->SetBinError(iChn+1,calbP1err);
            
            h_Ped_H[iSign][iChn]->GetXaxis()->SetRange(100,400);
            //h_Ped_H[iSign][iChn]->GetXaxis()->SetRange(50,150);
            
            int nEntries  = h_Ped_H[iSign][iChn]->GetEntries();
            int histMeanH = h_Ped_H[iSign][iChn]->GetMean();
            int histRMSH  = h_Ped_H[iSign][iChn]->GetStdDev();
            int fitL1 = histMeanH-3*histRMSH;
            int fitR1 = histMeanH+3*histRMSH;
            cout<<"====== TEST ======= "<<iChn<<"\t"<<nEntries<<endl;
            if (nEntries<100) continue;
            
            
            TF1 *funcPed;
            funcPed = new TF1("funcPed","gaus",fitL1,fitR1);
            funcPed->SetParameters(nEntries, histMeanH, histRMSH);
            funcPed->SetParLimits(1,histMeanH-2,histMeanH+2);
            h_Ped_H[iSign][iChn]->Fit("funcPed","R","",fitL1,fitR1);
            
            
            float pedMean  = funcPed->GetParameter(1);
            float pedSigma = funcPed->GetParameter(2);
            float pedMeanErr  = funcPed->GetParError(1);
            float pedSigmaErr = funcPed->GetParError(2);
            h_PedMean[iSign] ->SetBinContent(iChn+1,pedMean);
            h_PedSigma[iSign]->SetBinContent(iChn+1,pedSigma);
            //h_PedMean[iSign] ->SetBinError(iChn+1,pedMeanErr);
            //h_PedSigma[iSign]->SetBinError(iChn+1,pedSigmaErr);
            opt1<<runName<<"\t"<<"\t"<<iSign<<"\t"<<iChn<<"\t"<<pedMean<<"\t"<<pedSigma<<endl;
            cout<<runName<<"\t"<<"\t"<<iSign<<"\t"<<iChn<<"\t"<<pedMean<<"\t"<<pedSigma<<endl;
            ofile->cd();
            ofile->cd("PedestalRun");
            h_Ped_H[iSign][iChn]->Write();
            ofile->cd();
            //delete funcPed;
            //delete funcPhy;
        }
       
        h_PedMean[iSign]->Write();
        h_PedSigma[iSign]->Write();

        h_calbMap_p0[iSign]->Write();
        h_calbMap_p1[iSign]->Write();

        //////
        for (int iplug=0; iplug<8; iplug++) {
            TCanvas *cPed = new TCanvas(TString::Format("cPed_%1d",iplug),TString::Format("cPed_%1d",iplug),1500,900);
            cPed->Divide(4,4,0,0);
            TCanvas *cPhy = new TCanvas(TString::Format("cPhy_%1d",iplug),TString::Format("cPhy_%1d",iplug),1500,900);
            cPhy->Divide(4,4,0,0);
            string PmtConnName = PmtConName[0][iplug];
            int plugPosition   = PmtConnPos[0][iplug];
            for (int lp=0;lp<16;lp++){
                
                ///////////
                //string ConnectorName = ChannelMap[ ConnPosition ][lp].ConnName; //AXON connector name: E, F,...
                int pmtConID   = lp+1;
                int asicNo     = ChannelMap[ plugPosition ][lp].asicID;//
                int asicChnNo  = ChannelMap[ plugPosition ][lp].asicChnID;//
                int abcDataChn = ChannelMap[ plugPosition ][lp].ABC_ChnID;
                //if (abcDataChn==60) continue;
                if(PottingMap.count(PmtConnName)!=0){
                    int pID = PottingMap[ PmtConnName ][lp].PMT_ID;
                    int pSN = PottingMap[ PmtConnName ][lp].SerialNo;
                    cout<<"TEST Channel: "<<iplug<<"\t"<<lp<<"\t"<<pmtConID<<"\t"<<PmtConnName<<"\t"<<pID<<"\t"<<abcDataChn<<endl;
                    //if(PmtInfo.count(pSN)==0) {
                    //    cout<<"Chcek SN! "<<PmtConnName<<"\t"<<pmtConID<<"\t"<<pSN<<" Not in the list!"<<endl;
                    //    continue;
                    //}
                   
                    cPed->cd(lp+1);
                    h_Ped_H[iSign][abcDataChn]->Draw();
                    h_Ped_H[iSign][abcDataChn]->Draw();
                    h_Ped_H[iSign][abcDataChn]->GetXaxis()->SetRangeUser(0,200);
                    
                    if(iSign==0) h_Ped_H[iSign][abcDataChn]->SetTitle(TString::Format("Ping: %s_sn%d_Chn%d", PmtConnName.c_str(),pSN,pID));
                    if(iSign==1) h_Ped_H[iSign][abcDataChn]->SetTitle(TString::Format("Pong: %s_sn%d_Chn%d", PmtConnName.c_str(),pSN,pID));
                     
                    TF1 *pedFunc = h_Ped_H[iSign][abcDataChn]->GetFunction("funcPed");
                    pedFunc->SetLineColor(2);
                    
                     
                    /*
                    if(iSign==0) h_Phy_H[iSign][abcDataChn]->SetTitle(TString::Format("Ping: %s_sn%d_Chn%d", PmtConnName.c_str(),pSN,pID));
                    if(iSign==1) h_Phy_H[iSign][abcDataChn]->SetTitle(TString::Format("Pong: %s_sn%d_Chn%d", PmtConnName.c_str(),pSN,pID));
                    
                    TF1 *speFunc = h_Phy_H[iSign][abcDataChn]->GetFunction("funcPhy");
                    speFunc->SetLineColor(2);
                     */
                    
                }
            }
            if (iSign==0 && iplug==0 ) {
                //cPhy->Print(TString::Format("Fit_%s.pdf(",physRun.c_str()),"pdf");
                cPed->Print(TString::Format("Fit_%s.pdf(",runName.c_str()),"pdf");
            }else if (iSign==1 && iplug==7 ) {
                //cPhy->Print(TString::Format("Fit_%s.pdf)",physRun.c_str()),"pdf");
                cPed->Print(TString::Format("Fit_%s.pdf)",runName.c_str()),"pdf");
            }else {
                //cPhy->Print(TString::Format("Fit_%s.pdf",physRun.c_str()),"pdf");
                cPed->Print(TString::Format("Fit_%s.pdf",runName.c_str()),"pdf");
            }
            
            //c1->Print(("./result/"+SubRunName+markedF+".pdf(").c_str(),"pdf");
            
        }
        
    }
    

    ofile->Close();
    
   /////////////
    /*
    TCanvas *c13 = new TCanvas("c13","c13",1400,900);
    c13->Divide(2,2,0.001,0);
    c13->cd(1);
    
    legPing->AddEntry(h_ref[0], TString::Format("Ref_190 values"),"lep");
    //legPing->AddEntry(h_Gain[0], TString::Format("Run: %s", physRun.c_str()),"lep");
    
    h_ref[0]->Draw();
    h_ref[0]->SetLineColor(1);
    h_Gain[0]->Draw("same");
    h_ref[0]->SetStats(0);
    h_Gain[0]->SetStats(0);
    h_ref[0]->SetTitle("Ping");
    h_ref[0]->GetYaxis()->SetLabelFont( 63  );
    h_ref[0]->GetYaxis()->SetLabelSize( 18 );
    h_ref[0]->GetYaxis()->SetRangeUser(1.0*1e+6,6.5*1e+6);
    h_ref[0]->GetYaxis()->SetTitle("Gain");
    h_ref[0]->GetYaxis()->SetTitleOffset(0.8);
    h_ref[0]->GetYaxis()->SetTitleSize(0.06);
    legPing->Draw();
    
    c13->cd(2);
    legPong->AddEntry(h_ref[1], TString::Format("Ref_190 values"),"lep");
    //legPong->AddEntry(h_Gain[1], TString::Format("Run: %s", physRun.c_str()),"lep");
    
    
    h_ref[1]->Draw();
    h_ref[1]->SetLineColor(2);
    h_Gain[1]->Draw("same");
    h_ref[1]->SetStats(0);
    h_Gain[1]->SetStats(0);
    h_ref[1]->SetTitle("Pong");
    h_ref[1]->GetYaxis()->SetLabelFont( 63  );
    h_ref[1]->GetYaxis()->SetLabelSize( 18 );
    h_ref[1]->GetYaxis()->SetRangeUser(1.0*1e+6,6.5*1e+6);
    h_ref[1]->GetYaxis()->SetTitle("Gain");
    h_ref[1]->GetYaxis()->SetTitleOffset(0.8);
    h_ref[1]->GetYaxis()->SetTitleSize(0.06);
    legPong->Draw();
    
    c13->cd(3);
    h_compare[0]->Draw();
    h_compare[0]->SetStats(0);
    h_compare[0]->SetTitle("Comparision");
    h_compare[0]->GetYaxis()->SetTitle("Diff: (Ref-ThisTest)/Ref ");
    h_compare[0]->GetYaxis()->SetTitleOffset(0.8);
    h_compare[0]->GetYaxis()->SetTitleSize(0.06);
    h_compare[0]->GetYaxis()->SetLabelSize( 0.05 );
    h_compare[0]->GetXaxis()->SetTitle("ASIC Channel");
    h_compare[0]->GetXaxis()->SetTitleOffset(0.8);
    h_compare[0]->GetXaxis()->SetTitleSize(0.06);
    h_compare[0]->GetYaxis()->SetRangeUser(-1.0,1.0);
    
    c13->cd(4);
    h_compare[1]->Draw();
    h_compare[1]->SetStats(0);
    h_compare[1]->SetTitle("Comparision");
    h_compare[1]->GetYaxis()->SetTitle("Diff: (Ref-ThisTest)/Ref ");
    h_compare[1]->GetYaxis()->SetTitleOffset(0.8);
    h_compare[1]->GetYaxis()->SetTitleSize(0.06);
    //h_compare[1]->GetYaxis()->SetLabelSize( 18 );
    h_compare[1]->GetXaxis()->SetTitle("ASIC Channel");
    h_compare[1]->GetXaxis()->SetTitleOffset(0.8);
    h_compare[1]->GetXaxis()->SetTitleSize(0.06);
    h_compare[1]->GetYaxis()->SetRangeUser(-1.0,1.0);

    //c13->Print( TString::Format("UwbTestResults_%s.pdf(",physRun.c_str()),"pdf");
    
    TCanvas *c11 = new TCanvas("c11","c11",1000,900);
    c11->Divide(1,2,0.001,0);
    c11->cd(1);
    legDK->AddEntry(h_dkRef, TString::Format("Ref. values"),"lep");
    //legDK->AddEntry(h_dkRate, TString::Format("Run: %s", physRun.c_str()),"lep");
    
    h_dkRate->GetXaxis()->SetTitle("Channel ID");
    h_dkRate->GetYaxis()->SetTitle("DkRate");
    h_dkRate->GetYaxis()->SetLabelFont( 63  );
    h_dkRate->GetYaxis()->SetLabelSize( 18 );
    //h_dkRate->GetYaxis()->SetRangeUser(8*1e+5,10.0*1e+5);
    h_dkRate->Draw();
    h_dkRate->GetYaxis()->SetRangeUser(0,2.0*1e+3);
    h_dkRate->GetYaxis()->SetTitle("DarkRate [totalEntries/RunningTime]");
    h_dkRate->GetYaxis()->SetTitleOffset(0.8);
    h_dkRate->GetYaxis()->SetTitleSize(0.06);
    h_dkRef->Draw("same");
    h_dkRef->SetLineColor(1);
    h_dkRate->SetStats(0);
    h_dkRef->SetStats(0);
    legDK->Draw();
    c11->cd(2);
    
    h_dkCompare->Draw();
    h_dkCompare->SetStats(0);
    h_dkCompare->GetYaxis()->SetTitle("Diff: (Ref-ThisTest)/Ref ");
    h_dkCompare->GetYaxis()->SetTitleOffset(0.8);
    h_dkCompare->GetYaxis()->SetTitleSize(0.06);
    h_dkCompare->GetXaxis()->SetTitle("ASIC Channel");
    h_dkCompare->GetXaxis()->SetTitleOffset(0.8);
    h_dkCompare->GetXaxis()->SetTitleSize(0.06);

    //c11->SaveAs(TString::Format("DK_%s.eps",physRun.c_str()));
    c11->Print( TString::Format("UwbTestResults_%s.pdf)",runName.c_str()),"pdf");
    */
    
    ResetHist();
    
}

int BeginJob()
{

  ///////////////////////////////////
  legPing = new TLegend(0.6,0.7,1,1);
  legPong = new TLegend(0.6,0.7,1,1);
  legDK = new TLegend(0.6,0.7,1,1);
  for (int iSign = 0; iSign<=1; iSign++) {
      h_PedMean[iSign] = new TH1F(TString::Format("P%1d_h_PedMean",iSign),
                                  TString::Format("P%1d_h_PedMean",iSign),
                                  130, 0, 130);
      h_PhyMean[iSign] = new TH1F(TString::Format("P%1d_h_PhyMean",iSign),
                                  TString::Format("P%1d_h_PhyMean",iSign),
                                  130, 0, 130);
      
      h_PedSigma[iSign] = new TH1F(TString::Format("P%1d_h_PedSigma",iSign),
                                   TString::Format("P%1d_h_PedSigma",iSign),
                                   130, 0, 130);
      h_PhySigma[iSign] = new TH1F(TString::Format("P%1d_h_PhySigma",iSign),
                                   TString::Format("P%1d_h_PhySigma",iSign),
                                   130, 0, 130);
      
      h_QMean[iSign] = new TH1F(TString::Format("P%1d_h_QMean",iSign),
                                  TString::Format("P%1d_h_QMean",iSign),
                                  130, 0, 130);
      h_QSigma[iSign] = new TH1F(TString::Format("P%1d_h_QSigma",iSign),
                                   TString::Format("P%1d_h_QSigma",iSign),
                                   130, 0, 130);
      
      h_Gain[iSign] = new TH1F(TString::Format("P%1d_h_Gain",iSign),
                               TString::Format("P%1d_h_Gain",iSign),
                               130, 0, 130);
      h_Gain[iSign]->SetLineColor(iSign+3);
      h_Gain[iSign]->SetMarkerColor(iSign+3);
      
      h_EvtNum[iSign] = new TH1F(TString::Format("P%1d_h_EvtNum",iSign),
                                 TString::Format("P%1d_h_EvtNum",iSign),
                                 130, 0, 130);
      h_EvtNumSig[iSign] = new TH1F(TString::Format("P%1d_h_EvtNumSig",iSign),
                                 TString::Format("P%1d_h_EvtNumSig",iSign),
                                 130, 0, 130);
      
      h_calbMap_p0[iSign] = new TH1F(TString::Format("P%1d_h_calbMap_p0",iSign),
                                 TString::Format("P%1d_h_calbMap_p0",iSign),
                                 130, 0, 130);
      h_calbMap_p1[iSign] = new TH1F(TString::Format("P%1d_h_calbMap_p1",iSign),
                                 TString::Format("P%1d_h_calbMap_p1",iSign),
                                 130, 0, 130);
      
      h_runT[iSign] = new TH1F(TString::Format("P%1d_h_runT",iSign),
                                     TString::Format("P%1d_h_runT",iSign),
                                     130, 0, 130);
      
      h_compare[iSign] = new TH1F(TString::Format("P%1d_Compare",iSign),
                                  TString::Format("P%1d_Compare",iSign),
                                  130, 0, 130);
      h_ref[iSign] = new TH1F(TString::Format("P%1d_Ref",iSign),
                                  TString::Format("P%1d_Ref",iSign),
                                  130, 0, 130);
      
  }
    h_dkRate = new TH1F(TString::Format("dkRate"),
                               TString::Format("dkRate"),
                               130, 0, 130);
    h_dkRate->SetLineColor(4);
    h_dkRef = new TH1F(TString::Format("h_dkRef"),
                               TString::Format("h_dkRef"),
                               130, 0, 130);

  return 1;
}


int ResetHist(){

  ///////////////////////////////////
  
  for (int iSign = 0; iSign<=1; iSign++) {
      delete h_PedMean[iSign];
      delete h_PhyMean[iSign];
      delete h_PedSigma[iSign];
      delete h_PhySigma[iSign];
      delete h_QMean[iSign];
      delete h_QSigma[iSign];
      delete h_Gain[iSign];
      delete h_EvtNum[iSign];
      delete h_EvtNumSig[iSign];
      delete h_calbMap_p0[iSign];
      delete h_calbMap_p1[iSign];
      delete h_runT[iSign];
      delete h_compare[iSign];
      
  }
    delete h_dkRate;
    delete h_dkRef;
  return 1;
  
}



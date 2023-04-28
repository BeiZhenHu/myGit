// RawData to Root file for firmware 2.0
//
// by Bei-Zhen Hu 2020.03.23
//
// USAGE: ./main datapath datalist.txt GCU_no ABC_no
// [datapath]: data path must with "/" in the end
//             ex: ~/path/to/data/
// 2020.11.12 update: cyNB = cycN[ iChnNo ][thisMode]; ==>> cyNB = cycN[ iChnNo ][thisMode]-1;
// 2021.1.26 Duplicated for IHEP DAVIS
// 2022.9.8  for JUNO onsite quick analysis
// 2022.9.17 update output filename format
// 2023.4.26 for firmware 2.0

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include "TH1F.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>
#include <TMath.h>
#include "TStyle.h"

#include "dataDcode_v2.h"
//#include <TApplication.h>
using namespace std;

#define lenOfHderPart1 4
#define lenOfHderPart2 4
#define lenOfEvtChnData 8

unsigned int isHeader = 3138532539;
unsigned int isTrailer= 3141738939;
unsigned int testHeaderID;
unsigned int Header_EvtTimeStmp;  // 5Byte

struct EventBrief {
    int runNB, chnNo;
    int rCT, rFT;
    int qADC, pipoM, lhGain;
    int nbOfRChn;
};

map<int /* RunNo */, vector< EventBrief > > PrevEvtPing;
map<int /* RunNo */, vector< EventBrief > > PrevEvtPong;

unsigned int headerDcode(char *charName, int len){
    testHeaderID = 0;
    int k = len;
    for (int i=0; i<len; i++) {
        k--;
        int shiftBit = 8*k;
        testHeaderID += ( (charName[i] & 0xFFu) << shiftBit );
    }
    return testHeaderID;
}



int main(int argc,char** argv) {
    
    if( argc != 3 ) {

        std::cerr << "# of arguments error." << std::endl;
        std::cerr << "USAGE: ./main datapath datalist.txt" << std::endl;
        return 1;
    }

    string rawDataPath = argv[1];
    string runName[50];
    ifstream rawDataList;
    rawDataList.open(argv[2]);
    
    
    int i=0;
    while (!rawDataList.eof()) {
        rawDataList>>runName[i];
        i++;
    }
    int fileNum = i-1;
    
    cout<<fileNum<<"\t"<<i<<" "<<runName[fileNum-1]<<endl;
    for (int lpFile=0; lpFile<fileNum; lpFile++) {
        int run_nb = atoi(  TString::Format("%s",  runName[lpFile].substr(17,4).c_str())  );
        int gcuNo  = atoi(  TString::Format("%s",  runName[lpFile].substr(0,4).c_str())  );
        cout<<"RunNo: "<<run_nb<<endl;
        cout<<"gcuNo: "<<gcuNo<<endl;
        ///////////////////////////////////////////////
        
        string FileName;
        FileName = TString::Format("%s",runName[lpFile].c_str());
        cout<<FileName<<endl;
        
        ifstream myFile ((rawDataPath+FileName).c_str(), ios::in | ios::binary);
        cout<<"test path: "<<(rawDataPath+FileName).c_str()<<endl;
        ///////////////////////////////////
        
        string NewRunName;
        NewRunName = TString::Format("dryrun_%d_GU%d.root",run_nb, gcuNo);
        cout<<"Test Name formate: "<<NewRunName.c_str()<<endl;
        
        
        TFile *ofile = new TFile(NewRunName.c_str(),"recreate");
        
        TTree *aEvent = new TTree("Event","Event Tree");    /// create tree
        int RunNo = 0;
        int ChnN = 0;
        int ADCU = 0;
        int FT   = 0;
        int CT   = 0;
        int mode = 0;
        int gLH  = 0;
        int cyNB = 0;
        
        aEvent->Branch("RunNo",    &RunNo, "RunNo/I");
        aEvent->Branch("ChnN",     &ChnN,  "ChnN/I");
        aEvent->Branch("ADCU",     &ADCU,  "ADCU/I");
        aEvent->Branch("FT",       &FT  ,  "FT/I");
        aEvent->Branch("CT",       &CT  ,  "CT/I");
        aEvent->Branch("mode",     &mode,  "mode/I");
        aEvent->Branch("gLH",      &gLH,   "gLH/I");
        aEvent->Branch("cyNB",     &cyNB,  "cyNB/I");
        
        ///////////////////////////////////////////////
        
        ofile->cd();
        ofile->mkdir("HighGain");
        ofile->cd("HighGain");
        TH1F *h_Adc_H[2][130];
        for (int iSign = 0; iSign<=1; iSign++) {    // ping:0, pong:1
            for (int iChn=0; iChn<128; iChn++) {
                if (iSign==0) {
                    h_Adc_H[iSign][iChn] = new TH1F(TString::Format("Ping_Chn%1d_HG",iChn),
                                                    TString::Format("Ping_Chn%1d_HG",iChn),
                                                    1024, 0, 1024);
                }
                if (iSign==1) {
                    h_Adc_H[iSign][iChn] = new TH1F(TString::Format("Pong_Chn%1d_HG",iChn),
                                                    TString::Format("Pong_Chn%1d_HG",iChn),
                                                    1024, 0, 1024);
                }
                h_Adc_H[iSign][iChn]->GetXaxis()->SetTitle("ADCu");
            }
        }
        
        ofile->cd();
        ofile->mkdir("LowGain");
        ofile->cd("LowGain");
        TH1F *h_Adc_L[2][130];
        for (int iSign = 0; iSign<=1; iSign++) {
            for (int iChn=0; iChn<128; iChn++) {
                if (iSign==0) {
                    h_Adc_L[iSign][iChn] = new TH1F(TString::Format("Ping_Chn%1d_LG",iChn),
                                                    TString::Format("Ping_Chn%1d_LG",iChn),
                                                    1024, 0, 1024);
                }
                if (iSign==1) {
                    h_Adc_L[iSign][iChn] = new TH1F(TString::Format("Pong_Chn%1d_LG",iChn),
                                                    TString::Format("Pong_Chn%1d_LG",iChn),
                                                    1024, 0, 1024);
                }
                
                h_Adc_L[iSign][iChn]->GetXaxis()->SetTitle("ADCu");
            }
        }
        
        ///////////////////////////////////////////////
        ///////////////////////////////////////////////
        
        DataInfo testInfo;
        //struct stat results;
        char    ABC_header1[8];
        int     loopTotalN = 0;
        int     BlockID = 0;
        int     cycN[500][2] = {0};
        unsigned int HeaderID = 0;
        
        while (!myFile.eof()) {

            myFile.read (ABC_header1, 8);
            HeaderID = headerDcode(ABC_header1, lenOfHderPart1);
           
            if (HeaderID == isHeader) {
                //cout<<"1----------is HeaderID: "<<HeaderID<<endl;
                BlockID = ((ABC_header1[4] & 0xF0u)>>4);
                //cout<<"BlockID: "<< BlockID <<endl;
            }else if(HeaderID == isTrailer){
                //cout<<"2---------- isTrailer: "<<HeaderID<<endl;
                //cout<<"N of channels on block: "<< (ABC_header1[4] & 0xFFu)<<endl;
                
            }else{
                //cout<<"3----------- not Trailer or HeaderID: "<<HeaderID<<endl;
               
                testInfo.inputData(ABC_header1,lenOfEvtChnData);
                    
                int evtType = testInfo.EventType();
                
                if (evtType==0 || evtType==1 ) {
                
                    int iRunNo = run_nb;
                    int iChnNo = testInfo.ChnNum();
                    int abcChn = BlockID*8+iChnNo;
                    int thisFT = testInfo.FineTime();
                    int thisCT = testInfo.CoasTime();
                    int thisADC= testInfo.Charge();
                    int gLHMode= testInfo.GMode();
                    int evtN   = testInfo.EvtCounter();
                    int pipo   = evtN%2;
                    //cout<<"---- test event type: "<<evtType<<endl;
                    //cout<<"-- ChannelID: "<<abcChn<<" CT: "<<thisCT<<" gLH: "<<gLHMode<<"\t ping/pong: "<<pipo<<" Charge: "<<thisADC<<" testQ: "<<testQ<<endl;
                    
                   ////////////////////
                    
                    if (pipo==0){ //PING
                        int nEvt = PrevEvtPing[ abcChn ].size();
                        
                        if (nEvt==1) {
                            int preCT = PrevEvtPing[ abcChn ][0].rCT;
                            
                            cyNB = cycN[ abcChn ][pipo];
                            RunNo= PrevEvtPing[ abcChn ][0].runNB;
                            ChnN = PrevEvtPing[ abcChn ][0].chnNo;
                            gLH  = PrevEvtPong[ abcChn ][0].lhGain;
                            ADCU = PrevEvtPing[ abcChn ][0].qADC;
                            FT   = PrevEvtPing[ abcChn ][0].rFT;
                            CT   = PrevEvtPing[ abcChn ][0].rCT;
                            mode = PrevEvtPing[ abcChn ][0].pipoM;
                           
                            aEvent->Fill();
                            
                            int dT = thisCT-preCT;
                            if (dT<-60000000) cycN[ abcChn ][pipo]++;// cycle counter.
                            vector<EventBrief>::iterator k =  PrevEvtPing[ abcChn ].begin();
                            PrevEvtPing[ abcChn ].erase(k);
                        }
                    }
                    if (pipo==1){ //PONG
                        int nEvt = PrevEvtPong[ abcChn ].size();
                        if (nEvt==1) {
                            int preCT = PrevEvtPong[ abcChn ][0].rCT;
                            
                            cyNB = cycN[ abcChn ][pipo];
                            RunNo= PrevEvtPong[ abcChn ][0].runNB;
                            ChnN = PrevEvtPong[ abcChn ][0].chnNo;
                            gLH  = PrevEvtPong[ abcChn ][0].lhGain;
                            ADCU = PrevEvtPong[ abcChn ][0].qADC;
                            FT   = PrevEvtPong[ abcChn ][0].rFT;
                            CT   = PrevEvtPong[ abcChn ][0].rCT;
                            mode = PrevEvtPong[ abcChn ][0].pipoM;
                            
                            aEvent->Fill();
                            int dT = thisCT-preCT;
                            if (dT<-60000000) cycN[ abcChn ][pipo]++;// cycle counter.
                            vector<EventBrief>::iterator k =  PrevEvtPong[ abcChn ].begin();
                            PrevEvtPong[ abcChn ].erase(k);
                        }
                    }
                   ////////////////////
                    EventBrief preEvt;
                    preEvt.runNB = iRunNo;
                    preEvt.chnNo = abcChn;
                    preEvt.rCT   = thisCT;
                    preEvt.rFT   = thisFT;
                    preEvt.qADC  = thisADC;
                    preEvt.pipoM = pipo;
                    preEvt.lhGain= gLHMode;
                    if (gLHMode==1)   {
                        ofile->cd();
                        ofile->cd("HighGain");
                        h_Adc_H[pipo][abcChn]->Fill(thisADC);
                    }
                    else{
                       
                        ofile->cd();
                        ofile->cd("LowGain");
                        h_Adc_L[pipo][abcChn]->Fill(thisADC);
                    }
                    
                    if (pipo==0) PrevEvtPing[ abcChn ].push_back( preEvt );
                    if (pipo==1) PrevEvtPong[ abcChn ].push_back( preEvt );
                
                }
            }
            
        }//while loop
        /////////
        PrevEvtPing.erase(PrevEvtPing.begin(), PrevEvtPing.end());
        PrevEvtPong.erase(PrevEvtPong.begin(), PrevEvtPong.end());
        myFile.close();
        
        ofile->cd();
        ofile->Write();
        aEvent->Write();
        ofile->Close();
    }
    return 0;
}


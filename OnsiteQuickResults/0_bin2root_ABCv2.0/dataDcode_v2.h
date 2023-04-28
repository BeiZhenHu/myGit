// ABC Firmware 2.1 ref.DocDB-9065

class DataInfo{

public:
    void inputData( char*, int );
    int  EventType();
    int  ChnNum();
    int  CoasTime();
    int  GMode();
    int  EvtCounter();
    int  FineTime();
    int  Charge();
    
    int  DDS_ChnNum();
    int  DDS_EvtCounter();
    int  DDS_TrigWT();
    int  DDS_TrigFET();
    
    


private:
    char inputChar[10];
};

void DataInfo::inputData(char *charName, int len) {
    for (int i=0; i<len; i++)   inputChar[i] = charName[i];
}

int DataInfo::EventType(){
    int shiftBit = 6;
    return  (inputChar[0] & 0xC0u)>>shiftBit;
}

//      0        1        2        3        4        5        6        7
//      |ttChn#<-|---------- Coarse Time--->|G<----EvtN-->|<-Charge->|<---FT--->|
//      |ttnnnn--|--------|--------|--------|G-------|--------|--------|--------|

int DataInfo::ChnNum(){
    int shiftBit = 2;
    return  (inputChar[0] & 0x3Cu)>>shiftBit;
}

int DataInfo::CoasTime(){
    
    unsigned int CoarseT = 0;
    int k = 4;
    for (int i=0; i<4; i++) {
        k--;
        int shiftBit = 8*k;
        
        if (i==0)   CoarseT += ( (inputChar[i] & 0x03u) << shiftBit );
        else        CoarseT += ( (inputChar[i] & 0xFFu) << shiftBit );
    }
    return  CoarseT;
}
int DataInfo::GMode(){
    int shiftBit = 7;
    return  (inputChar[4] & 0x80u)>>shiftBit;
}
int DataInfo::EvtCounter(){
    
    unsigned int EvtN = 0;
    int k = 2;
    for (int i=4; i<6; i++) {
        k--;
        int shiftBit = 8*k;
        
        if (i==4)   EvtN += ( (inputChar[i] & 0x7Fu) << shiftBit ); //01111111
        else        EvtN += ( (inputChar[i] & 0xF0u) << shiftBit ); //11110000
    }
    return  EvtN>>4;
}


int DataInfo::Charge(){
    unsigned int vQ = 0;
    int k = 2;
    for (int i=5; i<7; i++) {
        k--;
        int shiftBit = 8*k;
        
        if (i==5)   vQ += ( (inputChar[i] & 0x0Fu) << shiftBit ); //00001111
        else        vQ += ( (inputChar[i] & 0xFCu) << shiftBit ); //11111100
    }
    
    return  vQ>>2;
}


int DataInfo::FineTime(){
    unsigned int FineT = 0;
    int k = 2;
    for (int i=6; i<8; i++) {
        k--;
        int shiftBit = 8*k;
        
        if (i==6)   FineT += ( (inputChar[i] & 0x03u) << shiftBit );//00000011
        else        FineT += ( (inputChar[i] & 0xFFu) << shiftBit );//11111111
    }
    return  FineT;
}

// Following for firmware V2.1
//      0        1        2        3        4        5        6        7
//      |ttnnnn--|--------|--------|--------|G-------|--------|--------|--------|
//      |tt0nnnn|<------EvtN------->|<-TrgW>|<------------TrigFET-------------->|

int DataInfo::DDS_ChnNum(){
    int shiftBit = 1;
    return  (inputChar[0] & 0x1Eu)>>shiftBit; //00011110
}

int DataInfo::DDS_EvtCounter(){
    unsigned int ddsEvtN = 0;
    int k = 4;
    for (int i=0; i<4; i++) {
        k--;
        int shiftBit = 8*k;
        
        if (i==0)   ddsEvtN += ( (inputChar[i] & 0x01u) << shiftBit ); //00000001
        if (i==3)   ddsEvtN += ( (inputChar[i] & 0x80u) << shiftBit ); //10000000
        else        ddsEvtN += ( (inputChar[i] & 0xFFu) << shiftBit );
    }
    return  ddsEvtN>>7;
}


int DataInfo::DDS_TrigWT(){
    
    unsigned int ddsTrigWT = 0;
    int k = 1;
    for (int i=3; i<4; i++) {
        k--;
        int shiftBit = 8*k;
        
        ddsTrigWT += ( (inputChar[i] & 0x7Fu) << shiftBit ); //01111111
    }
    return  ddsTrigWT;
}


int DataInfo::DDS_TrigFET(){
    unsigned int ddsTrigFET = 0;
    int k = 4;
    for (int i=4; i<8; i++) {
        k--;
        int shiftBit = 8*k;
        ddsTrigFET += ( (inputChar[i] & 0xFFu) << shiftBit ); //11111111
        
    }
    return  ddsTrigFET;
}


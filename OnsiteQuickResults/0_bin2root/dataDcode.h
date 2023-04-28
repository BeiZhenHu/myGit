//      
//      0        1        2        3        4        5        6        7        8
//      |<-ChnN->|<--------- Coarse Time ----------->|<---- Charge --->|<-- Fine Time -->|
//      |--------|xxxxxx--|--------|--------|--------|sxxxxx--|--------|xxxxxx--|--------|
//
//

class DataInfo{

public:
    void inputData( char*, int );
    int  ChnNum();
    int  CoasTime();
    int  FineTime();
    int  Charge();
    int  mode();

private:
    char inputChar[10];
};

void DataInfo::inputData(char *charName, int len) {
    for (int i=0; i<len; i++)   inputChar[i] = charName[i];
}

int DataInfo::ChnNum(){
    return  (inputChar[0] & 0xFFu);
}

//      0        1        2        3        4        5
//      |<-ChnN->|<--------- Coarse Time ----------->|<----
//      |--------|xxxxxx--|--------|--------|--------|sxxxxx--

int DataInfo::CoasTime(){
    
    unsigned int CoarseT = 0;
    int k = 4;
    for (int i=1; i<5; i++) {
        k--;
        int shiftBit = 8*k;
        
        if (i==1)   CoarseT += ( (inputChar[i] & 0x03u) << shiftBit );
        else        CoarseT += ( (inputChar[i] & 0xFFu) << shiftBit );
    }
    return  CoarseT;
}

int DataInfo::FineTime(){
    
    unsigned int FineT = 0;
    int k = 2;
    for (int i=7; i<8; i++) {
        k--;
        int shiftBit = 8*k;
        
        if (i==7)   FineT += ( (inputChar[i] & 0x03u) << shiftBit );
        else        FineT += ( (inputChar[i] & 0xFFu) << shiftBit );
    }
    
    return  FineT;
}

//             5        6        7        8
//      ------>|<---- Charge --->|<-- Fine Time -->|
//      -------|sxxxxx--|--------|xxxxxx--|--------|

int DataInfo::Charge(){
    unsigned int vQ = 0;
    int k = 2;
    for (int i=5; i<7; i++) {
        k--;
        int shiftBit = 8*k;
        
        if (i==5)   vQ += ( (inputChar[i] & 0x03u) << shiftBit );
        else        vQ += ( (inputChar[i] & 0xFFu) << shiftBit );
    }
    return  vQ;
}

int DataInfo::mode(){
    int shiftBit = 7;
    return  (inputChar[5] & 0x80u)>>shiftBit;
}



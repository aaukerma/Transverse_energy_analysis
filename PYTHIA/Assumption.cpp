/******************************************************************************
This program takes specified output runs and makes graphs to show the validity
of the assumptions made. These assumptions include the pion0 assumption, the
kaon0 (long and short) assumptions, the neutron and antineutron assumptions,
as well as the eta and omega assumptions.
******************************************************************************/

#include "Riostream.h"
#include <cstdio>
#include <vector>
#include <string>
#include <fstream>
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1F.h"
#include "TCanvas.h"
#include <iostream>
using namespace std;

struct BIN {
  Double_t EoP;
  Int_t BinCount;
};

struct Data {
  BIN ET;
  BIN PT;
};

int Assumption(){
  Int_t IDNUM=1;
  string a = "=======================================";
  char* filename = Form("hhAnalysisOfRun%i",IDNUM);
  vector<Int_t> SNN {7,11,19,27,39,62,130,200,900,2760,5020,7000,8000,13000,14000};
  vector<Float_t> SNNACTUAL {7.7,11.6,19.6,27,39,62.4,130,200,900,2760,5020,7000,8000,13000,14000};
  TFile* file = TFile::Open(filename, "RECREATE");
  ofstream DUMP;
  DUMP.open("./PythiaData.txt");
  DUMP<<"Pythia Data"<<endl<<a<<endl<<"run number: "<<IDNUM<<endl;
  for (int x;x<14;x++){
    char* temp = Form("hh%iGeVoutfile%i.root",SNN[x],IDNUM);
    DUMP<<endl<<"SNN=="<<SNN[x]<<endl<<a<<endl;
    TFile f(temp);
    if (f.IsZombie()) {
        cout<<"error opening file"<<endl;
        exit(-1);
    }
    TH1F* hNEvents= (TH1F*)f.Get("hNEvents")->Clone();
    TH1F* hETAll= (TH1F*)f.Get("hETAll")->Clone();
    TH1F* hETPiPlus= (TH1F*)f.Get("hETPiPlus")->Clone();
    TH1F* hETPiMinus= (TH1F*)f.Get("hETPiMinus")->Clone();
    TH1F* hETPi0= (TH1F*)f.Get("hETPi0")->Clone();
    TH1F* hETKPlus= (TH1F*)f.Get("hETKPlus")->Clone();
    TH1F* hETKMinus= (TH1F*)f.Get("hETKMinus")->Clone();
    TH1F* hETKL= (TH1F*)f.Get("hETKL")->Clone();
    TH1F* hETKS= (TH1F*)f.Get("hETKS")->Clone();
    TH1F* hETEta= (TH1F*)f.Get("hETEta")->Clone();
    TH1F* hETomega= (TH1F*)f.Get("hETomega")->Clone();
    TH1F* hETLambda0= (TH1F*)f.Get("hETLambda0")->Clone();
    TH1F* hETLambdaBar0= (TH1F*)f.Get("hETLambdaBar0")->Clone();
    TH1F* hETp= (TH1F*)f.Get("hETp")->Clone();
    TH1F* hETn= (TH1F*)f.Get("hETn")->Clone();
    TH1F* hETpBar= (TH1F*)f.Get("hETpBar")->Clone();
    TH1F* hETnBar= (TH1F*)f.Get("hETnBar")->Clone();
    TH1F* hPTAll= (TH1F*)f.Get("hPTAll")->Clone();
    TH1F* hptPiPlus= (TH1F*)f.Get("hptPiPlus")->Clone();
    TH1F* hptPiMinus= (TH1F*)f.Get("hptPiMinus")->Clone();
    TH1F* hptPi0= (TH1F*)f.Get("hptPi0")->Clone();
    TH1F* hptKPlus= (TH1F*)f.Get("hptKPlus")->Clone();
    TH1F* hptKMinus= (TH1F*)f.Get("hptKMinus")->Clone();
    TH1F* hptKL= (TH1F*)f.Get("hptKL")->Clone();
    TH1F* hptKS= (TH1F*)f.Get("hptKS")->Clone();
    TH1F* hptEta= (TH1F*)f.Get("hptEta")->Clone();
    TH1F* hptomega= (TH1F*)f.Get("hptomega")->Clone();
    TH1F* hptLambda0= (TH1F*)f.Get("hptLambda0")->Clone();
    TH1F* hptLambdaBar0= (TH1F*)f.Get("hptLambdaBar0")->Clone();
    TH1F* hptp= (TH1F*)f.Get("hptp")->Clone();
    TH1F* hptn= (TH1F*)f.Get("hptn")->Clone();
    TH1F* hptpBar= (TH1F*)f.Get("hptpBar")->Clone();
    TH1F* hptnBar= (TH1F*)f.Get("hptnBar")->Clone();
    f.Close();
    cout<<"line 89 is confused"<<endl;
    //Int_t t = hNEvents->GetNbinsX();
    //cout<<t<<endl;
    //for (int i; i<t;i++){

      //Double_t TEMP = hNEvents->GetBinContent(hNEvents->FindBin(0));
      cout<<"line 91."<<x<<endl;
      cout<<hNEvents<<endl;
      //DUMP<<"Bin Contents: "<<TEMP<<endl;
      cout<<endl;
    //}
  }




































return 0;
}

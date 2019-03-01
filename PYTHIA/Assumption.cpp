/******************************************************************************
This program takes specified output runs and makes graphs to show the validity
of the assumptions made. These assumptions include the pion0 assumption, the
kaon0 (long and short) assumptions, the neutron and antineutron assumptions,
as well as the eta and omega assumptions.

File outputs: ETBinContent, ETBinLowEdge, PTBinContent, PTBinLowEdge
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
  Double_t BinCount;
};

struct Data {
  BIN ET;
  BIN PT;
};

int Assumption(){
  Int_t IDNUM=1;
  string a = "=======================================";
  string b = "---------------------------------------";
  char* filename = Form("hhAnalysisOfRun%i",IDNUM);
  vector<Int_t> SNN {7,11,19,27,39,62,130,200,900,2760,5020,7000,8000,13000,14000};
  vector<double> SNNACTUAL {7.7,11.6,19.6,27,39,62.4,130,200,900,2760,5020,7000,8000,13000,14000};
  BIN datum;
  Data set;

  TFile* file = TFile::Open(filename, "RECREATE");
  ofstream DUMP;
  DUMP.open("./PythiaData.txt");
  DUMP<<"Pythia Data"<<endl<<"run number: "<<IDNUM<<endl;
  for (int x=0;x<14;x++){
    char* temp = Form("hh%iGeVoutfile%i.root",SNN[x],IDNUM);
    DUMP<<endl<<"SNN=="<<SNN[x]<<endl<<a<<endl;
    TFile f(temp);
    f.Open(temp);
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
    Int_t t = hNEvents->GetNbinsX();
    for (int i=0; i<t;i++){
      Double_t TEMP = hNEvents->GetBinContent(hNEvents->FindBin(0));
      DUMP<<"Number of Events: "<<TEMP<<endl;
    }
    DUMP<<"ET Total || PT Total"<<endl<<b<<endl;
    DUMP<<"# \t\t"<<"ET \t\t"<<"# \t\t"<<"PT"<<endl<<b<<endl;
    Int_t w = hETAll->GetNbinsX();
    for (int i=0; i<w;i++){
      set.ET.EoP = hETAll->GetBinContent(hETAll->FindBin(i));
      set.ET.BinCount = (hETAll->GetBinLowEdge(i));
      set.PT.EoP = hPTAll->GetBinContent(hPTAll->FindBin(i));
      set.PT.BinCount = (hPTAll->GetBinLowEdge(i));
      DUMP<<set.ET.EoP<<"\t\t"<<set.ET.BinCount<<"\t\t"<<set.PT.EoP<<"\t\t"<<set.PT.BinCount<<endl;
    }
    //pi+
    DUMP<<b<<endl<<"pi+"<<endl;
    DUMP<<"# \t\t"<<"ET \t\t"<<"# \t\t"<<"PT"<<endl<<b<<endl;
    for (int i=0; i<w;i++){
      set.ET.EoP = hETPiPlus->GetBinContent(hETPiPlus->FindBin(i));
      set.ET.BinCount = (hETPiPlus->GetBinLowEdge(i));
      set.PT.EoP = hptPiPlus->GetBinContent(hptPiPlus->FindBin(i));
      set.PT.BinCount = (hptPiPlus->GetBinLowEdge(i));
      DUMP<<set.ET.EoP<<"\t\t"<<set.ET.BinCount<<"\t\t"<<set.PT.EoP<<"\t\t"<<set.PT.BinCount<<endl;
    }
    //pi-
    DUMP<<b<<endl<<"pi-"<<endl;
    DUMP<<"# \t\t"<<"ET \t\t"<<"# \t\t"<<"PT"<<endl<<b<<endl;
    for (int i=0; i<w;i++){
      set.ET.EoP = hETPiMinus->GetBinContent(hETPiMinus->FindBin(i));
      set.ET.BinCount = (hETPiMinus->GetBinLowEdge(i));
      set.PT.EoP = hptPiMinus->GetBinContent(hptPiMinus->FindBin(i));
      set.PT.BinCount = (hptPiMinus->GetBinLowEdge(i));
      DUMP<<set.ET.EoP<<"\t\t"<<set.ET.BinCount<<"\t\t"<<set.PT.EoP<<"\t\t"<<set.PT.BinCount<<endl;
    }
    //pi0
    DUMP<<b<<endl<<"pi0"<<endl;
    DUMP<<"# \t\t"<<"ET \t\t"<<"# \t\t"<<"PT"<<endl<<b<<endl;
    for (int i=0; i<w;i++){
      set.ET.EoP = hETPi0->GetBinContent(hETPi0->FindBin(i));
      set.ET.BinCount = (hETPi0->GetBinLowEdge(i));
      set.PT.EoP = hptPi0->GetBinContent(hptPi0->FindBin(i));
      set.PT.BinCount = (hptPi0->GetBinLowEdge(i));
      DUMP<<set.ET.EoP<<"\t\t"<<set.ET.BinCount<<"\t\t"<<set.PT.EoP<<"\t\t"<<set.PT.BinCount<<endl;
    }
    //K+
    DUMP<<b<<endl<<"K+"<<endl;
    DUMP<<"# \t\t"<<"ET \t\t"<<"# \t\t"<<"PT"<<endl<<b<<endl;
    for (int i=0; i<w;i++){
      set.ET.EoP = hETKPlus->GetBinContent(hETKPlus->FindBin(i));
      set.ET.BinCount = (hETKPlus->GetBinLowEdge(i));
      set.PT.EoP = hptKPlus->GetBinContent(hptKPlus->FindBin(i));
      set.PT.BinCount = (hptKPlus->GetBinLowEdge(i));
      DUMP<<set.ET.EoP<<"\t\t"<<set.ET.BinCount<<"\t\t"<<set.PT.EoP<<"\t\t"<<set.PT.BinCount<<endl;
    }
    //K-
    DUMP<<b<<endl<<"K-"<<endl;
    DUMP<<"# \t\t"<<"ET \t\t"<<"# \t\t"<<"PT"<<endl<<b<<endl;
    for (int i=0; i<w;i++){
      set.ET.EoP = hETKMinus->GetBinContent(hETKMinus->FindBin(i));
      set.ET.BinCount = (hETKMinus->GetBinLowEdge(i));
      set.PT.EoP = hptKMinus->GetBinContent(hptKMinus->FindBin(i));
      set.PT.BinCount = (hptKMinus->GetBinLowEdge(i));
      DUMP<<set.ET.EoP<<"\t\t"<<set.ET.BinCount<<"\t\t"<<set.PT.EoP<<"\t\t"<<set.PT.BinCount<<endl;
    }
    //K0L
    DUMP<<b<<endl<<"K0L"<<endl;
    DUMP<<"# \t\t"<<"ET \t\t"<<"# \t\t"<<"PT"<<endl<<b<<endl;
    for (int i=0; i<w;i++){
      set.ET.EoP = hETKL->GetBinContent(hETKL->FindBin(i));
      set.ET.BinCount = (hETKL->GetBinLowEdge(i));
      set.PT.EoP = hptKL->GetBinContent(hptKL->FindBin(i));
      set.PT.BinCount = (hptKL->GetBinLowEdge(i));
      DUMP<<set.ET.EoP<<"\t\t"<<set.ET.BinCount<<"\t\t"<<set.PT.EoP<<"\t\t"<<set.PT.BinCount<<endl;
    }
    //K0S
    DUMP<<b<<endl<<"K0S"<<endl;
    DUMP<<"# \t\t"<<"ET \t\t"<<"# \t\t"<<"PT"<<endl<<b<<endl;
    for (int i=0; i<w;i++){
      set.ET.EoP = hETKS->GetBinContent(hETKS->FindBin(i));
      set.ET.BinCount = (hETKS->GetBinLowEdge(i));
      set.PT.EoP = hptKS->GetBinContent(hptKS->FindBin(i));
      set.PT.BinCount = (hptKS->GetBinLowEdge(i));
      DUMP<<set.ET.EoP<<"\t\t"<<set.ET.BinCount<<"\t\t"<<set.PT.EoP<<"\t\t"<<set.PT.BinCount<<endl;
    }
    //eta
    DUMP<<b<<endl<<"eta"<<endl;
    DUMP<<"# \t\t"<<"ET \t\t"<<"# \t\t"<<"PT"<<endl<<b<<endl;
    for (int i=0; i<w;i++){
      set.ET.EoP = hETEta->GetBinContent(hETEta->FindBin(i));
      set.ET.BinCount = (hETEta->GetBinLowEdge(i));
      set.PT.EoP = hptEta->GetBinContent(hptEta->FindBin(i));
      set.PT.BinCount = (hptEta->GetBinLowEdge(i));
      DUMP<<set.ET.EoP<<"\t\t"<<set.ET.BinCount<<"\t\t"<<set.PT.EoP<<"\t\t"<<set.PT.BinCount<<endl;
    }
    //omega
    DUMP<<b<<endl<<"omega"<<endl;
    DUMP<<"# \t\t"<<"ET \t\t"<<"# \t\t"<<"PT"<<endl<<b<<endl;
    for (int i=0; i<w;i++){
      set.ET.EoP = hETomega->GetBinContent(hETomega->FindBin(i));
      set.ET.BinCount = (hETomega->GetBinLowEdge(i));
      set.PT.EoP = hptomega->GetBinContent(hptomega->FindBin(i));
      set.PT.BinCount = (hptomega->GetBinLowEdge(i));
      DUMP<<set.ET.EoP<<"\t\t"<<set.ET.BinCount<<"\t\t"<<set.PT.EoP<<"\t\t"<<set.PT.BinCount<<endl;
    }
    //Lambda0
    DUMP<<b<<endl<<"Lambda0"<<endl;
    DUMP<<"# \t\t"<<"ET \t\t"<<"# \t\t"<<"PT"<<endl<<b<<endl;
    for (int i=0; i<w;i++){
      set.ET.EoP = hETLambda0->GetBinContent(hETLambda0->FindBin(i));
      set.ET.BinCount = (hETLambda0->GetBinLowEdge(i));
      set.PT.EoP = hptLambda0->GetBinContent(hptLambda0->FindBin(i));
      set.PT.BinCount = (hptLambda0->GetBinLowEdge(i));
      DUMP<<set.ET.EoP<<"\t\t"<<set.ET.BinCount<<"\t\t"<<set.PT.EoP<<"\t\t"<<set.PT.BinCount<<endl;
    }
    //LambdaBar0
    DUMP<<b<<endl<<"LambdaBar0"<<endl;
    DUMP<<"# \t\t"<<"ET \t\t"<<"# \t\t"<<"PT"<<endl<<b<<endl;
    for (int i=0; i<w;i++){
      set.ET.EoP = hETLambdaBar0->GetBinContent(hETLambdaBar0->FindBin(i));
      set.ET.BinCount = (hETLambdaBar0->GetBinLowEdge(i));
      set.PT.EoP = hptLambdaBar0->GetBinContent(hptLambdaBar0->FindBin(i));
      set.PT.BinCount = (hptLambdaBar0->GetBinLowEdge(i));
      DUMP<<set.ET.EoP<<"\t\t"<<set.ET.BinCount<<"\t\t"<<set.PT.EoP<<"\t\t"<<set.PT.BinCount<<endl;
    }
    //p
    DUMP<<b<<endl<<"proton"<<endl;
    DUMP<<"# \t\t"<<"ET \t\t"<<"# \t\t"<<"PT"<<endl<<b<<endl;
    for (int i=0; i<w;i++){
      set.ET.EoP = hETp->GetBinContent(hETp->FindBin(i));
      set.ET.BinCount = (hETp->GetBinLowEdge(i));
      set.PT.EoP = hptp->GetBinContent(hptp->FindBin(i));
      set.PT.BinCount = (hptp->GetBinLowEdge(i));
      DUMP<<set.ET.EoP<<"\t\t"<<set.ET.BinCount<<"\t\t"<<set.PT.EoP<<"\t\t"<<set.PT.BinCount<<endl;
    }
    //n
    DUMP<<b<<endl<<"neutron"<<endl;
    DUMP<<"# \t\t"<<"ET \t\t"<<"# \t\t"<<"PT"<<endl<<b<<endl;
    for (int i=0; i<w;i++){
      set.ET.EoP = hETn->GetBinContent(hETn->FindBin(i));
      set.ET.BinCount = (hETn->GetBinLowEdge(i));
      set.PT.EoP = hptn->GetBinContent(hptn->FindBin(i));
      set.PT.BinCount = (hptn->GetBinLowEdge(i));
      DUMP<<set.ET.EoP<<"\t\t"<<set.ET.BinCount<<"\t\t"<<set.PT.EoP<<"\t\t"<<set.PT.BinCount<<endl;
    }
    //pBar
    DUMP<<b<<endl<<"antiproton"<<endl;
    DUMP<<"# \t\t"<<"ET \t\t"<<"# \t\t"<<"PT"<<endl<<b<<endl;
    for (int i=0; i<w;i++){
      set.ET.EoP = hETpBar->GetBinContent(hETpBar->FindBin(i));
      set.ET.BinCount = (hETpBar->GetBinLowEdge(i));
      set.PT.EoP = hptpBar->GetBinContent(hptpBar->FindBin(i));
      set.PT.BinCount = (hptpBar->GetBinLowEdge(i));
      DUMP<<set.ET.EoP<<"\t\t"<<set.ET.BinCount<<"\t\t"<<set.PT.EoP<<"\t\t"<<set.PT.BinCount<<endl;
    }
    //nBar
    DUMP<<b<<endl<<"antineutron"<<endl;
    DUMP<<"# \t\t"<<"ET \t\t"<<"# \t\t"<<"PT"<<endl<<b<<endl;
    for (int i=0; i<w;i++){
      set.ET.EoP = hETnBar->GetBinContent(hETnBar->FindBin(i));
      set.ET.BinCount = (hETnBar->GetBinLowEdge(i));
      set.PT.EoP = hptnBar->GetBinContent(hptnBar->FindBin(i));
      set.PT.BinCount = (hptnBar->GetBinLowEdge(i));
      DUMP<<set.ET.EoP<<"\t\t"<<set.ET.BinCount<<"\t\t"<<set.PT.EoP<<"\t\t"<<set.PT.BinCount<<endl;
    }
  }
return 0;
}

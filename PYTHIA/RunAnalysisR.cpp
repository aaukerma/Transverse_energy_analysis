/******************************************************************************
This program takes data file "PythiaData.txt". outputs the VALUEs in a
stacked graph

File outputs: ETBinContent, ETBinLowEdge, PTBinContent, PTBinLowEdge
******************************************************************************/


#include "Riostream.h"
#include <cstdio>
#include <vector>
#include <string>
#include <fstream>
#include "TString.h"
#include "TFile.h"
#include "TGraph.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1F.h"
#include "TCanvas.h"
#include <unistd.h>
#include <iostream>
using namespace std;

struct dat {
  Double_t VALUE;
  Double_t stddev;
  Double_t entries;
};
struct d{
    double V;
    double s;
    int n;
    double vsTOT;
};

struct txtdat{
  int NUM;
  double ET;
  double PT;
  double pET;
  double pPT;
};

double GetValues(TH1F* histo){
  double total=0;
  Double_t BinCent=0;
  Double_t BinCont=0;
  double maths=0;
  Int_t binCount=histo->GetNbinsX();
  for (Int_t i=0;i<binCount;i++){
    BinCent=histo->GetBinCenter(i);
    BinCont=histo->GetBinContent(i);
    maths=BinCent*BinCont;
    total+=maths;
  }
  return total;
}

int RunAnalysisR(){
  Int_t IDNUM=1;
  Int_t a=0;//type count 0-1
  Int_t b=0;//part count 0-15
  Int_t c=0;//info count 0-2
  vector<Int_t> SNN {7,11,19,27,39,62,130,200,2760,5020,7000,8000,13000,14000};
  vector<double> SNNACTUAL {7.7,11.5,19.6,27,39,62.4,130,200,2760,5020,7000,8000,13000,14000};
  /*****************************************************************************
  THE FOLLOWING is a reference for counters a, b ,and c
  a{"ET","PT"};
  b{"ALL","PiP","PiM","Pi0","KaP","KaM","K0L","K0S","ETA","OMG","La0","LB0","PRO","NEU","PRB","NEB","OMm","XI0","XIM","SIP","SIM","SI0"};
  c{"AVG","S","N"};
  *****************************************************************************/
  vector<vector<vector<dat>>> Stuff(10, vector<vector<dat>>(2, vector<dat>(16)));
  vector<vector<vector<d>>> StuffTXT(10, vector<vector<d>>(2, vector<d>(24)));
  vector<vector<vector<d>>> TOTs(10, vector<vector<d>>(2, vector<d>(9)));

  //VALUEs for each histo
  Double_t ETALLAVG;
  Double_t ETPiPAVG;
  Double_t ETPiMAVG;
  Double_t ETPi0AVG;
  Double_t ETKaPAVG;
  Double_t ETKaMAVG;
  Double_t ETK0LAVG;
  Double_t ETK0SAVG;
  Double_t ETETAAVG;
  Double_t ETOMGAVG;
  Double_t ETLa0AVG;
  Double_t ETLB0AVG;
  Double_t ETPROAVG;
  Double_t ETNEUAVG;
  Double_t ETPRBAVG;
  Double_t ETNEBAVG;

  Double_t PTALLAVG;
  Double_t PTPiPAVG;
  Double_t PTPiMAVG;
  Double_t PTPi0AVG;
  Double_t PTKaPAVG;
  Double_t PTKaMAVG;
  Double_t PTK0LAVG;
  Double_t PTK0SAVG;
  Double_t PTETAAVG;
  Double_t PTOMGAVG;
  Double_t PTLa0AVG;
  Double_t PTLB0AVG;
  Double_t PTPROAVG;
  Double_t PTNEUAVG;
  Double_t PTPRBAVG;
  Double_t PTNEBAVG;

  //stddev each histo
  Double_t ETALLS;
  Double_t ETPiPS;
  Double_t ETPiMS;
  Double_t ETPi0S;
  Double_t ETKaPS;
  Double_t ETKaMS;
  Double_t ETK0LS;
  Double_t ETK0SS;
  Double_t ETETAS;
  Double_t ETOMGS;
  Double_t ETLa0S;
  Double_t ETLB0S;
  Double_t ETPROS;
  Double_t ETNEUS;
  Double_t ETPRBS;
  Double_t ETNEBS;

  Double_t PTALLS;
  Double_t PTPiPS;
  Double_t PTPiMS;
  Double_t PTPi0S;
  Double_t PTKaPS;
  Double_t PTKaMS;
  Double_t PTK0LS;
  Double_t PTK0SS;
  Double_t PTETAS;
  Double_t PTOMGS;
  Double_t PTLa0S;
  Double_t PTLB0S;
  Double_t PTPROS;
  Double_t PTNEUS;
  Double_t PTPRBS;
  Double_t PTNEBS;

  //number of entries each histo
  Double_t ETALLN;
  Double_t ETPiPN;
  Double_t ETPiMN;
  Double_t ETPi0N;
  Double_t ETKaPN;
  Double_t ETKaMN;
  Double_t ETK0LN;
  Double_t ETK0SN;
  Double_t ETETAN;
  Double_t ETOMGN;
  Double_t ETLa0N;
  Double_t ETLB0N;
  Double_t ETPRON;
  Double_t ETNEUN;
  Double_t ETPRBN;
  Double_t ETNEBN;

  Double_t PTALLN;
  Double_t PTPiPN;
  Double_t PTPiMN;
  Double_t PTPi0N;
  Double_t PTKaPN;
  Double_t PTKaMN;
  Double_t PTK0LN;
  Double_t PTK0SN;
  Double_t PTETAN;
  Double_t PTOMGN;
  Double_t PTLa0N;
  Double_t PTLB0N;
  Double_t PTPRON;
  Double_t PTNEUN;
  Double_t PTPRBN;
  Double_t PTNEBN;

  Double_t BinCent;
  Double_t BinCont;
  Double_t BinErr;
  double e=0;
  double d=0;
  double STDDEVfh=0;

  string TRASH;
  int JOBIDD;
  double NUMBER;
  char ty;

  cout<<"Select Run Number (int): "<<endl;
  cin>>IDNUM;
  cout<<"pp or pA?:"<<endl;
  cin>>ty;
  cout<<"preparing data for run: "<<IDNUM<<endl<<endl;
  for (int x=0;x<10;x++){
    double TET;
    double TPT;
    int num;
    double et;
    double pt;
    double pet;
    double ppt;
    char* temp = Form("p%c%i_%iGeVoutfile.root",ty,IDNUM,SNN[x]);
    char* temptxt = Form("p%c%i_%iGeVoutfile.txt",ty,IDNUM,SNN[x]);
    TFile f(temp);
    f.Open(temp);
    ifstream f2;
    f2.open(temptxt);
    if (f.IsZombie()) {
        cout<<"error opening file"<<endl;
        exit(-1);
    }
    f2>>TRASH>>TRASH>>TRASH>>TRASH>>TRASH>>JOBIDD;
    //cout<<"JOB: "<<JOBIDD<<"\t";
    //cout<<TRASH<<endl;
    f2>>TRASH>>TRASH>>NUMBER;
    cout<<"No.Events: "<<NUMBER<<"\t";
    f2>>TRASH>>NUMBER;
    f2>>TRASH>>TRASH;
    f2>>TRASH>>TRASH>>TRASH;
    f2>>TRASH>>TRASH>>TRASH>>TRASH>>TRASH;
    f2>>TRASH>>NUMBER;
    f2>>TRASH>>TRASH>>TRASH>>TRASH;
    f2>>TRASH>>TRASH>>TRASH;
    f2>>TRASH>>TRASH>>TRASH>>NUMBER;
    TET=NUMBER;
    cout<<"Total ET: "<<TET<<"\t";
    f2>>TRASH>>TRASH>>TRASH>>NUMBER;
    TPT=NUMBER;
    cout<<"Total PT: "<<TPT<<endl;
    f2>>TRASH>>TRASH>>TRASH>>TRASH>>TRASH>>TRASH;
    //start DATA
    //b=pi+,pi-,pi0,K+,K-,K0L,K0S,Lam,Lam_,p,pBAR,n,nBAR,eta,[eta],omega,[omega],Om-,Xi0,Xi-,Sig+,Sig-,Sig0,other
    for (int b=0;b<21;b++){
      f2>>TRASH>>num>>et>>pt>>pet>>ppt;
      //cout<<TRASH<<" "<<num<<" "<<et<<" "<<pt<<" "<<pet<<" "<<ppt<<endl;
      StuffTXT[x][0][b].V=et;
      StuffTXT[x][1][b].V=pt;
      StuffTXT[x][0][b].n=num;
      StuffTXT[x][1][b].n=num;
      StuffTXT[x][0][b].vsTOT=pet;
      StuffTXT[x][1][b].vsTOT=ppt;
      double sd= .1/sqrt(num);
      StuffTXT[x][0][b].s=sd;
      StuffTXT[x][1][b].s=sd;
      //cout<<TRASH<<endl;
    }
    //start TOTALS
    //b=pions,kaons,lams,p/pbar,n/nbar,eta/ome,other
    for (int b=0;b<8;b++){
      //cout<<TRASH<<" "<<num<<" "<<et<<" "<<pt<<" "<<pet<<" "<<ppt<<endl;
      f2>>TRASH>>num>>et>>pt>>pet>>ppt;
      TOTs[x][0][b].V=et;
      TOTs[x][1][b].V=pt;
      TOTs[x][0][b].n=num;
      TOTs[x][1][b].n=num;
      TOTs[x][0][b].vsTOT=pet;
      TOTs[x][1][b].vsTOT=ppt;
      double sd= .1/sqrt(num);
      TOTs[x][0][b].s=sd;
      TOTs[x][1][b].s=sd;
      //cout<<TRASH<<endl;
    }
    cout<<"I made it here Line 243\n";
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
    TH1F* hETXi0= (TH1F*)f.Get("hETXi0")->Clone();
    TH1F* hETXim= (TH1F*)f.Get("hETXim")->Clone();
    TH1F* hETSigmap= (TH1F*)f.Get("hETSigmap")->Clone();
    TH1F* hETSigmam= (TH1F*)f.Get("hETSigmam")->Clone();
    TH1F* hETSigma0= (TH1F*)f.Get("hETSigma0")->Clone();
    TH1F* hETOMEGAm= (TH1F*)f.Get("hETOMEGAm")->Clone();
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
    TH1F* hptXi0= (TH1F*)f.Get("hptXi0")->Clone();
    TH1F* hptXim= (TH1F*)f.Get("hptXim")->Clone();
    TH1F* hptSigmap= (TH1F*)f.Get("hptSigmap")->Clone();
    TH1F* hptSigmam= (TH1F*)f.Get("hptSigmam")->Clone();
    TH1F* hptSigma0= (TH1F*)f.Get("hptSigma0")->Clone();
    TH1F* hptOMEGAm= (TH1F*)f.Get("hptOMEGAm")->Clone();
    TH1F* hptp= (TH1F*)f.Get("hptp")->Clone();
    TH1F* hptn= (TH1F*)f.Get("hptn")->Clone();
    TH1F* hptpBar= (TH1F*)f.Get("hptpBar")->Clone();
    TH1F* hptnBar= (TH1F*)f.Get("hptnBar")->Clone();
    cout<<"Histos Loaded\n";
    f.Close();

    a=0;
    b=0;
    d=0;
    STDDEVfh=0;
    for (Int_t i=0;i<1000;i++){
      BinCent=hETAll->GetBinCenter(i);
      BinCont=hETAll->GetBinContent(i);
      BinErr=sqrt(hETAll->GetBinError(i));
      e=BinCent*BinCont;
      d+=e;
      e=BinErr*BinCont;
      STDDEVfh+=e;
    }
    Stuff[x][a][b].VALUE=d;
    Stuff[x][a][b].stddev=STDDEVfh;
    Stuff[x][a][b].entries=hETAll->GetEntries();

    b++;
    d=0;
    STDDEVfh=0;
    for (Int_t i=0;i<1000;i++){
      BinCent=hETPiPlus->GetBinCenter(i);
      BinCont=hETPiPlus->GetBinContent(i);
      BinErr=sqrt(hETPiPlus->GetBinError(i));
      e=BinCent*BinCont;
      d+=e;
      e=BinErr*BinCont;
      STDDEVfh+=e;
    }
    Stuff[x][a][b].VALUE=d;
    Stuff[x][a][b].stddev=STDDEVfh;
    Stuff[x][a][b].entries=hETPiPlus->GetEntries();

    b++;
    d=0;
    STDDEVfh=0;
    for (Int_t i=0;i<1000;i++){
      BinCent=hETPiMinus->GetBinCenter(i);
      BinCont=hETPiMinus->GetBinContent(i);
      BinErr=sqrt(hETPiMinus->GetBinError(i));
      e=BinCent*BinCont;
      d+=e;
      e=BinErr*BinCont;
      STDDEVfh+=e;
    }
    Stuff[x][a][b].VALUE=d;
    Stuff[x][a][b].stddev=STDDEVfh;
    Stuff[x][a][b].entries=hETPiMinus->GetEntries();

    b++;
    d=0;
    STDDEVfh=0;
    for (Int_t i=0;i<1000;i++){
      BinCent=hETPi0->GetBinCenter(i);
      BinCont=hETPi0->GetBinContent(i);
      BinErr=sqrt(hETPi0->GetBinError(i));
      e=BinCent*BinCont;
      d+=e;
      e=BinErr*BinCont;
      STDDEVfh+=e;
    }
    Stuff[x][a][b].VALUE=d;
    Stuff[x][a][b].stddev=STDDEVfh;
    Stuff[x][a][b].entries=hETPi0->GetEntries();

    b++;
    d=0;
    STDDEVfh=0;
    for (Int_t i=0;i<1000;i++){
      BinCent=hETKPlus->GetBinCenter(i);
      BinCont=hETKPlus->GetBinContent(i);
      BinErr=sqrt(hETKPlus->GetBinError(i));
      e=BinCent*BinCont;
      d+=e;
      e=BinErr*BinCont;
      STDDEVfh+=e;
    }
    Stuff[x][a][b].VALUE=d;
    Stuff[x][a][b].stddev=STDDEVfh;
    Stuff[x][a][b].entries=hETKPlus->GetEntries();

    b++;
    d=0;
    STDDEVfh=0;
    for (Int_t i=0;i<1000;i++){
      BinCent=hETKMinus->GetBinCenter(i);
      BinCont=hETKMinus->GetBinContent(i);
      BinErr=sqrt(hETKMinus->GetBinError(i));
      e=BinCent*BinCont;
      d+=e;
      e=BinErr*BinCont;
      STDDEVfh+=e;
    }
    Stuff[x][a][b].VALUE=d;
    Stuff[x][a][b].stddev=STDDEVfh;
    Stuff[x][a][b].entries=hETKMinus->GetEntries();

    b++;
    d=0;
    STDDEVfh=0;
    for (Int_t i=0;i<1000;i++){
      BinCent=hETKL->GetBinCenter(i);
      BinCont=hETKL->GetBinContent(i);
      BinErr=sqrt(hETKL->GetBinError(i));
      e=BinCent*BinCont;
      d+=e;
      e=BinErr*BinCont;
      STDDEVfh+=e;
    }
    Stuff[x][a][b].VALUE=d;
    Stuff[x][a][b].stddev=STDDEVfh;
    Stuff[x][a][b].entries=hETKL->GetEntries();

    b++;
    d=0;
    STDDEVfh=0;
    for (Int_t i=0;i<1000;i++){
      BinCent=hETKS->GetBinCenter(i);
      BinCont=hETKS->GetBinContent(i);
      BinErr=sqrt(hETKS->GetBinError(i));
      e=BinCent*BinCont;
      d+=e;
      e=BinErr*BinCont;
      STDDEVfh+=e;
    }
    Stuff[x][a][b].VALUE=d;
    Stuff[x][a][b].stddev=STDDEVfh;
    Stuff[x][a][b].entries=hETKS->GetEntries();

    b++;
    d=0;
    STDDEVfh=0;
    for (Int_t i=0;i<1000;i++){
      BinCent=hETEta->GetBinCenter(i);
      BinCont=hETEta->GetBinContent(i);
      BinErr=sqrt(hETEta->GetBinError(i));
      e=BinCent*BinCont;
      d+=e;
      e=BinErr*BinCont;
      STDDEVfh+=e;
    }
    Stuff[x][a][b].VALUE=d;
    Stuff[x][a][b].stddev=STDDEVfh;
    Stuff[x][a][b].entries=hETEta->GetEntries();

    b++;
    d=0;
    STDDEVfh=0;
    for (Int_t i=0;i<1000;i++){
      BinCent=hETomega->GetBinCenter(i);
      BinCont=hETomega->GetBinContent(i);
      BinErr=sqrt(hETomega->GetBinError(i));
      e=BinCent*BinCont;
      d+=e;
      e=BinErr*BinCont;
      STDDEVfh+=e;
    }
    Stuff[x][a][b].VALUE=d;
    Stuff[x][a][b].stddev=STDDEVfh;
    Stuff[x][a][b].entries=hETomega->GetEntries();

    b++;
    d=0;
    STDDEVfh=0;
    for (Int_t i=0;i<1000;i++){
      BinCent=hETLambda0->GetBinCenter(i);
      BinCont=hETLambda0->GetBinContent(i);
      BinErr=sqrt(hETLambda0->GetBinError(i));
      e=BinCent*BinCont;
      d+=e;
      e=BinErr*BinCont;
      STDDEVfh+=e;
    }
    Stuff[x][a][b].VALUE=d;
    Stuff[x][a][b].stddev=STDDEVfh;
    Stuff[x][a][b].entries=hETLambda0->GetEntries();

    b++;
    d=0;
    STDDEVfh=0;
    for (Int_t i=0;i<1000;i++){
      BinCent=hETLambdaBar0->GetBinCenter(i);
      BinCont=hETLambdaBar0->GetBinContent(i);
      BinErr=sqrt(hETLambdaBar0->GetBinError(i));
      e=BinCent*BinCont;
      d+=e;
      e=BinErr*BinCont;
      STDDEVfh+=e;
    }
    Stuff[x][a][b].VALUE=d;
    Stuff[x][a][b].stddev=STDDEVfh;
    Stuff[x][a][b].entries=hETLambdaBar0->GetEntries();

    b++;
    d=0;
    STDDEVfh=0;
    for (Int_t i=0;i<1000;i++){
      BinCent=hETp->GetBinCenter(i);
      BinCont=hETp->GetBinContent(i);
      BinErr=sqrt(hETp->GetBinError(i));
      e=BinCent*BinCont;
      d+=e;
      e=BinErr*BinCont;
      STDDEVfh+=e;
    }
    Stuff[x][a][b].VALUE=d;
    Stuff[x][a][b].stddev=STDDEVfh;
    Stuff[x][a][b].entries=hETp->GetEntries();

    b++;
    d=0;
    STDDEVfh=0;
    for (Int_t i=0;i<1000;i++){
      BinCent=hETn->GetBinCenter(i);
      BinCont=hETn->GetBinContent(i);
      BinErr=sqrt(hETn->GetBinError(i));
      e=BinCent*BinCont;
      d+=e;
      e=BinErr*BinCont;
      STDDEVfh+=e;
    }
    Stuff[x][a][b].VALUE=d;
    Stuff[x][a][b].stddev=STDDEVfh;
    Stuff[x][a][b].entries=hETn->GetEntries();

    b++;
    d=0;
    STDDEVfh=0;
    for (Int_t i=0;i<1000;i++){
      BinCent=hETpBar->GetBinCenter(i);
      BinCont=hETpBar->GetBinContent(i);
      BinErr=sqrt(hETpBar->GetBinError(i));
      e=BinCent*BinCont;
      d+=e;
      e=BinErr*BinCont;
      STDDEVfh+=e;
    }
    Stuff[x][a][b].VALUE=d;
    Stuff[x][a][b].stddev=STDDEVfh;
    Stuff[x][a][b].entries=hETpBar->GetEntries();

    b++;
    d=0;
    STDDEVfh=0;
    for (Int_t i=0;i<1000;i++){
      BinCent=hETnBar->GetBinCenter(i);
      BinCont=hETnBar->GetBinContent(i);
      BinErr=sqrt(hETnBar->GetBinError(i));
      e=BinCent*BinCont;
      d+=e;
      e=BinErr*BinCont;
      STDDEVfh+=e;
    }
    Stuff[x][a][b].VALUE=d;
    Stuff[x][a][b].stddev=STDDEVfh;
    Stuff[x][a][b].entries=hETnBar->GetEntries();

    b=0;


    a++;
    d=0;
    STDDEVfh=0;
    for (Int_t i=0;i<1000;i++){
      BinCent=hPTAll->GetBinCenter(i);
      BinCont=hPTAll->GetBinContent(i);
      BinErr=sqrt(hPTAll->GetBinError(i));
      e=BinCent*BinCont;
      d+=e;
      e=BinErr*BinCont;
      STDDEVfh+=e;
    }
    Stuff[x][a][b].VALUE=d;
    Stuff[x][a][b].stddev=STDDEVfh;
    Stuff[x][a][b].entries=hPTAll->GetEntries();

    b++;
    d=0;
    STDDEVfh=0;
    for (Int_t i=0;i<1000;i++){
      BinCent=hptPiPlus->GetBinCenter(i);
      BinCont=hptPiPlus->GetBinContent(i);
      BinErr=sqrt(hptPiPlus->GetBinError(i));
      e=BinCent*BinCont;
      d+=e;
      e=BinErr*BinCont;
      STDDEVfh+=e;
    }
    Stuff[x][a][b].VALUE=d;
    Stuff[x][a][b].stddev=STDDEVfh;
    Stuff[x][a][b].entries=hptPiPlus->GetEntries();

    b++;
    d=0;
    STDDEVfh=0;
    for (Int_t i=0;i<1000;i++){
      BinCent=hptPiMinus->GetBinCenter(i);
      BinCont=hptPiMinus->GetBinContent(i);
      BinErr=sqrt(hptPiMinus->GetBinError(i));
      e=BinCent*BinCont;
      d+=e;
      e=BinErr*BinCont;
      STDDEVfh+=e;
    }
    Stuff[x][a][b].VALUE=d;
    Stuff[x][a][b].stddev=STDDEVfh;
    Stuff[x][a][b].entries=hptPiMinus->GetEntries();

    b++;
    d=0;
    STDDEVfh=0;
    for (Int_t i=0;i<1000;i++){
      BinCent=hptPi0->GetBinCenter(i);
      BinCont=hptPi0->GetBinContent(i);
      BinErr=sqrt(hptPi0->GetBinError(i));
      e=BinCent*BinCont;
      d+=e;
      e=BinErr*BinCont;
      STDDEVfh+=e;
    }
    Stuff[x][a][b].VALUE=d;
    Stuff[x][a][b].stddev=STDDEVfh;
    Stuff[x][a][b].entries=hptPi0->GetEntries();

    b++;
    d=0;
    STDDEVfh=0;
    for (Int_t i=0;i<1000;i++){
      BinCent=hptKPlus->GetBinCenter(i);
      BinCont=hptKPlus->GetBinContent(i);
      BinErr=sqrt(hptKPlus->GetBinError(i));
      e=BinCent*BinCont;
      d+=e;
      e=BinErr*BinCont;
      STDDEVfh+=e;
    }
    Stuff[x][a][b].VALUE=d;
    Stuff[x][a][b].stddev=STDDEVfh;
    Stuff[x][a][b].entries=hptKPlus->GetEntries();

    b++;
    d=0;
    STDDEVfh=0;
    for (Int_t i=0;i<1000;i++){
      BinCent=hptKMinus->GetBinCenter(i);
      BinCont=hptKMinus->GetBinContent(i);
      BinErr=sqrt(hptKMinus->GetBinError(i));
      e=BinCent*BinCont;
      d+=e;
      e=BinErr*BinCont;
      STDDEVfh+=e;
    }
    Stuff[x][a][b].VALUE=d;
    Stuff[x][a][b].stddev=STDDEVfh;
    Stuff[x][a][b].entries=hptKMinus->GetEntries();

    b++;
    d=0;
    STDDEVfh=0;
    for (Int_t i=0;i<1000;i++){
      BinCent=hptKL->GetBinCenter(i);
      BinCont=hptKL->GetBinContent(i);
      BinErr=sqrt(hptKL->GetBinError(i));
      e=BinCent*BinCont;
      d+=e;
      e=BinErr*BinCont;
      STDDEVfh+=e;
    }
    Stuff[x][a][b].VALUE=d;
    Stuff[x][a][b].stddev=STDDEVfh;
    Stuff[x][a][b].entries=hptKL->GetEntries();

    b++;
    d=0;
    STDDEVfh=0;
    for (Int_t i=0;i<1000;i++){
      BinCent=hptKS->GetBinCenter(i);
      BinCont=hptKS->GetBinContent(i);
      BinErr=sqrt(hptKS->GetBinError(i));
      e=BinCent*BinCont;
      d+=e;
      e=BinErr*BinCont;
      STDDEVfh+=e;
    }
    Stuff[x][a][b].VALUE=d;
    Stuff[x][a][b].stddev=STDDEVfh;
    Stuff[x][a][b].entries=hptKS->GetEntries();

    b++;
    d=0;
    STDDEVfh=0;
    for (Int_t i=0;i<1000;i++){
      BinCent=hptEta->GetBinCenter(i);
      BinCont=hptEta->GetBinContent(i);
      BinErr=sqrt(hptEta->GetBinError(i));
      e=BinCent*BinCont;
      d+=e;
      e=BinErr*BinCont;
      STDDEVfh+=e;
    }
    Stuff[x][a][b].VALUE=d;
    Stuff[x][a][b].stddev=STDDEVfh;
    Stuff[x][a][b].entries=hptEta->GetEntries();

    b++;
    d=0;
    STDDEVfh=0;
    for (Int_t i=0;i<1000;i++){
      BinCent=hptomega->GetBinCenter(i);
      BinCont=hptomega->GetBinContent(i);
      BinErr=sqrt(hptomega->GetBinError(i));
      e=BinCent*BinCont;
      d+=e;
      e=BinErr*BinCont;
      STDDEVfh+=e;
    }
    Stuff[x][a][b].VALUE=d;
    Stuff[x][a][b].stddev=STDDEVfh;
    Stuff[x][a][b].entries=hptomega->GetEntries();

    b++;
    d=0;
    STDDEVfh=0;
    for (Int_t i=0;i<1000;i++){
      BinCent=hptLambda0->GetBinCenter(i);
      BinCont=hptLambda0->GetBinContent(i);
      BinErr=sqrt(hptLambda0->GetBinError(i));
      e=BinCent*BinCont;
      d+=e;
      e=BinErr*BinCont;
      STDDEVfh+=e;
    }
    Stuff[x][a][b].VALUE=d;
    Stuff[x][a][b].stddev=STDDEVfh;
    Stuff[x][a][b].entries=hptLambda0->GetEntries();

    b++;
    d=0;
    STDDEVfh=0;
    for (Int_t i=0;i<1000;i++){
      BinCent=hptLambdaBar0->GetBinCenter(i);
      BinCont=hptLambdaBar0->GetBinContent(i);
      BinErr=sqrt(hptLambdaBar0->GetBinError(i));
      e=BinCent*BinCont;
      d+=e;
      e=BinErr*BinCont;
      STDDEVfh+=e;
    }
    Stuff[x][a][b].VALUE=d;
    Stuff[x][a][b].stddev=STDDEVfh;
    Stuff[x][a][b].entries=hptLambdaBar0->GetEntries();

    b++;
    d=0;
    STDDEVfh=0;
    for (Int_t i=0;i<1000;i++){
      BinCent=hptp->GetBinCenter(i);
      BinCont=hptp->GetBinContent(i);
      BinErr=sqrt(hptp->GetBinError(i));
      e=BinCent*BinCont;
      d+=e;
      e=BinErr*BinCont;
      STDDEVfh+=e;
    }
    Stuff[x][a][b].VALUE=d;
    Stuff[x][a][b].stddev=STDDEVfh;
    Stuff[x][a][b].entries=hptp->GetEntries();

    b++;
    d=0;
    STDDEVfh=0;
    for (Int_t i=0;i<1000;i++){
      BinCent=hptn->GetBinCenter(i);
      BinCont=hptn->GetBinContent(i);
      BinErr=sqrt(hptn->GetBinError(i));
      e=BinCent*BinCont;
      d+=e;
      e=BinErr*BinCont;
      STDDEVfh+=e;
    }
    Stuff[x][a][b].VALUE=d;
    Stuff[x][a][b].stddev=STDDEVfh;
    Stuff[x][a][b].entries=hptn->GetEntries();

    b++;
    d=0;
    STDDEVfh=0;
    for (Int_t i=0;i<1000;i++){
      BinCent=hptpBar->GetBinCenter(i);
      BinCont=hptpBar->GetBinContent(i);
      BinErr=sqrt(hptpBar->GetBinError(i));
      e=BinCent*BinCont;
      d+=e;
      e=BinErr*BinCont;
      STDDEVfh+=e;
    }
    Stuff[x][a][b].VALUE=d;
    Stuff[x][a][b].stddev=STDDEVfh;
    Stuff[x][a][b].entries=hptpBar->GetEntries();

    b++;
    d=0;
    STDDEVfh=0;
    for (Int_t i=0;i<1000;i++){
      BinCent=hptnBar->GetBinCenter(i);
      BinCont=hptnBar->GetBinContent(i);
      BinErr=sqrt(hptnBar->GetBinError(i));
      e=BinCent*BinCont;
      d+=e;
      e=BinErr*BinCont;
      STDDEVfh+=e;
    }
    Stuff[x][a][b].VALUE=d;
    Stuff[x][a][b].stddev=STDDEVfh;
    Stuff[x][a][b].entries=hptnBar->GetEntries();

    b=0;


    cout<<"Energy: "<<SNNACTUAL[x]<<" GeV complete"<<endl;
    //cout<<"ET pi+: "<<ETPiPAVG<<"\t"<<ETPiPS<<"\t"<<ETPiPN<<endl;
    //cout<<"PT pi+: "<<PTPiPAVG<<"\t"<<PTPiPS<<"\t"<<PTPiPN<<endl<<endl;
}
//cout<<"14GeV ET ALL"<<Stuff[14][0][0].VALUE<<"\t"<<Stuff[14][0][0].stddev<<"\t"<<Stuff[14][0][0].entries<<endl;
//cout<<"14GeV PT ALL"<<Stuff[14][1][0].VALUE<<"\t"<<Stuff[14][1][0].stddev<<"\t"<<Stuff[14][1][0].entries<<endl;
cout<<"input complete"<<endl<<endl;


/*******************************************************************************
Start graphing portion
*******************************************************************************/
cout<<"preparing all charts"<<endl<<endl;
//TCanvas *c0 =new TCanvas("c0","Total",200,10,700,500);

Int_t n=8;
Double_t x[n],y[n],ex[n],ey[n];
Double_t THINGY;

TLegend* legend=new TLegend(.127507,.8,.256447,.896842);
    legend->SetTextFont(72);
    legend->SetTextSize(.03);
    legend->SetFillColor(0);

TLegend* legend2=new TLegend(.719092,.115278,.895931,.256944);
    legend2->SetTextFont(72);
    legend2->SetTextSize(.03);
    legend2->SetFillColor(0);

TLegend* legend3=new TLegend(.712833,.709722,.895149,.893056);
    legend3->SetTextFont(72);
    legend3->SetTextSize(.03);
    legend3->SetFillColor(0);

TLegend* legend4=new TLegend(.743349,.781944,.897496,.895833);
    legend4->SetTextFont(72);
    legend4->SetTextSize(.03);
    legend4->SetFillColor(0);

TLegend* legend5=new TLegend(.103286,.784722,.322379,.893056);
    legend5->SetTextFont(72);
    legend5->SetTextSize(.03);
    legend5->SetFillColor(0);
{/*
//Total
for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  y[i]=Stuff[i][0][0].VALUE;
  ex[i]=0;
  THINGY=Stuff[i][0][0].entries;
  ey[i]=(Stuff[i][0][0].stddev)/sqrt(THINGY);
}
TGraphErrors* gALLet=new TGraphErrors(n,x,y,ex,ey);
  gALLet->SetTitle("Total ET avg; SNN(GeV); (GeV)");
  gALLet->SetMarkerStyle(20);
  gALLet->SetMarkerSize(.5);
  gALLet->SetMarkerColor(kRed);
  gALLet->Draw("AP");
c0->SetLogx();
c0->Update();
char* file0= Form("Total_Run%i.png",IDNUM);
c0->SaveAs(file0);
//pi+
TCanvas *c1 =new TCanvas("c1","Pi+",200,10,700,500);
for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  y[i]=Stuff[i][0][1].VALUE;
  ex[i]=0;
  THINGY=Stuff[i][0][1].entries;
  ey[i]=(Stuff[i][0][1].stddev)/sqrt(THINGY);
}
TGraphErrors* gpipet=new TGraphErrors(n,x,y,ex,ey);
  gpipet->SetTitle("Pi+ avg; SNN(GeV); (GeV)");
  gpipet->SetMarkerStyle(20);
  gpipet->SetMarkerSize(.5);
  gpipet->SetMarkerColor(kRed);
  legend->AddEntry(gpipet,"ET","p");
  gpipet->Draw("AP");

for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  y[i]=Stuff[i][1][1].VALUE;
  ex[i]=0;
  THINGY=Stuff[i][1][1].entries;
  ey[i]=(Stuff[i][1][1].stddev)/sqrt(THINGY);
}
TGraphErrors* gpippt=new TGraphErrors(n,x,y,ex,ey);
  gpippt->SetMarkerStyle(20);
  gpippt->SetMarkerSize(.5);
  gpippt->SetMarkerColor(kBlue);
  legend->AddEntry(gpippt,"PT","p");
  gpippt->Draw("P");
  legend->Draw();
c1->SetLogx();
c1->Update();
char* file1= Form("PIP_Run%i.png",IDNUM);
c1->SaveAs(file1);

//pi-
TCanvas *c2 =new TCanvas("c2","Pi-",200,10,700,500);
for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  y[i]=Stuff[i][0][2].VALUE;
  ex[i]=0;
  THINGY=Stuff[i][0][2].entries;
  ey[i]=(Stuff[i][0][2].stddev)/sqrt(THINGY);
}
TGraphErrors* gpimet=new TGraphErrors(n,x,y,ex,ey);
  gpimet->SetTitle("Pi- avg; SNN(GeV); (GeV)");
  gpimet->SetMarkerStyle(20);
  gpimet->SetMarkerSize(.5);
  gpimet->SetMarkerColor(kRed);
  gpimet->Draw("AP");

for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  y[i]=Stuff[i][1][2].VALUE;
  ex[i]=0;
  THINGY=Stuff[i][1][2].entries;
  ey[i]=(Stuff[i][1][2].stddev)/sqrt(THINGY);
}
TGraphErrors* gpimpt=new TGraphErrors(n,x,y,ex,ey);
  gpimpt->SetMarkerStyle(20);
  gpimpt->SetMarkerSize(.5);
  gpimpt->SetMarkerColor(kBlue);
  gpimpt->Draw("P");
  legend->Draw();
c2->SetLogx();
c2->Update();
char* file2= Form("PIM_Run%i.png",IDNUM);
c2->SaveAs(file2);

//pi0
TCanvas *c3 =new TCanvas("c3","Pi0",200,10,700,500);
for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  y[i]=Stuff[i][0][3].VALUE;
  ex[i]=0;
  THINGY=Stuff[i][0][3].entries;
  ey[i]=(Stuff[i][0][3].stddev)/sqrt(THINGY);
}
TGraphErrors* gpi0et=new TGraphErrors(n,x,y,ex,ey);
  gpi0et->SetTitle("Pi0 avg; SNN(GeV); (GeV)");
  gpi0et->SetMarkerStyle(20);
  gpi0et->SetMarkerSize(.5);
  gpi0et->SetMarkerColor(kRed);
  gpi0et->Draw("AP");

for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  y[i]=Stuff[i][1][3].VALUE;
  ex[i]=0;
  THINGY=Stuff[i][1][3].entries;
  ey[i]=(Stuff[i][1][3].stddev)/sqrt(THINGY);
}
TGraphErrors* gpi0pt=new TGraphErrors(n,x,y,ex,ey);
  gpi0pt->SetMarkerStyle(20);
  gpi0pt->SetMarkerSize(.5);
  gpi0pt->SetMarkerColor(kBlue);
  gpi0pt->Draw("P");
  legend->Draw();
c3->SetLogx();
c3->Update();
char* file3= Form("PI0_Run%i.png",IDNUM);
c3->SaveAs(file3);

//k+
TCanvas *c4 =new TCanvas("c4","K+",200,10,700,500);
for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  y[i]=Stuff[i][0][4].VALUE;
  ex[i]=0;
  THINGY=Stuff[i][0][4].entries;
  ey[i]=(Stuff[i][0][4].stddev)/sqrt(THINGY);
}
TGraphErrors* gkpet=new TGraphErrors(n,x,y,ex,ey);
  gkpet->SetTitle("K+ avg; SNN(GeV); (GeV)");
  gkpet->SetMarkerStyle(20);
  gkpet->SetMarkerSize(.5);
  gkpet->SetMarkerColor(kRed);
  gkpet->Draw("AP");

for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  y[i]=Stuff[i][1][4].VALUE;
  ex[i]=0;
  THINGY=Stuff[i][1][4].entries;
  ey[i]=(Stuff[i][1][4].stddev)/sqrt(THINGY);
}
TGraphErrors* gkppt=new TGraphErrors(n,x,y,ex,ey);
  gkppt->SetMarkerStyle(20);
  gkppt->SetMarkerSize(.5);
  gkppt->SetMarkerColor(kBlue);
  gkppt->Draw("P");
  legend->Draw();
c4->SetLogx();
c4->Update();
char* file4= Form("KaP_Run%i.png",IDNUM);
c4->SaveAs(file4);

//k-
TCanvas *c5 =new TCanvas("c5","K-",200,10,700,500);
for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  y[i]=Stuff[i][0][5].VALUE;
  ex[i]=0;
  THINGY=Stuff[i][0][5].entries;
  ey[i]=(Stuff[i][0][5].stddev)/sqrt(THINGY);
}
TGraphErrors* gkmet=new TGraphErrors(n,x,y,ex,ey);
  gkmet->SetTitle("K- avg; SNN(GeV); (GeV)");
  gkmet->SetMarkerStyle(20);
  gkmet->SetMarkerSize(.5);
  gkmet->SetMarkerColor(kRed);
  gkmet->Draw("AP");

for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  y[i]=Stuff[i][1][5].VALUE;
  ex[i]=0;
  THINGY=Stuff[i][1][5].entries;
  ey[i]=(Stuff[i][1][5].stddev)/sqrt(THINGY);
}
TGraphErrors* gkmpt=new TGraphErrors(n,x,y,ex,ey);
  gkmpt->SetMarkerStyle(20);
  gkmpt->SetMarkerSize(.5);
  gkmpt->SetMarkerColor(kBlue);
  gkmpt->Draw("P");
  legend->Draw();
c5->SetLogx();
c5->Update();
char* file5= Form("KaM_Run%i.png",IDNUM);
c5->SaveAs(file5);

//k0l
TCanvas *c6 =new TCanvas("c6","K0L",200,10,700,500);
for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  y[i]=Stuff[i][0][6].VALUE;
  ex[i]=0;
  THINGY=Stuff[i][0][6].entries;
  ey[i]=(Stuff[i][0][6].stddev)/sqrt(THINGY);
}
TGraphErrors* gklet=new TGraphErrors(n,x,y,ex,ey);
  gklet->SetTitle("K0L avg; SNN(GeV); (GeV)");
  gklet->SetMarkerStyle(20);
  gklet->SetMarkerSize(.5);
  gklet->SetMarkerColor(kRed);
  gklet->Draw("AP");

for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  y[i]=Stuff[i][1][6].VALUE;
  ex[i]=0;
  THINGY=Stuff[i][1][6].entries;
  ey[i]=(Stuff[i][1][6].stddev)/sqrt(THINGY);
}
TGraphErrors* gklpt=new TGraphErrors(n,x,y,ex,ey);
  gklpt->SetMarkerStyle(20);
  gklpt->SetMarkerSize(.5);
  gklpt->SetMarkerColor(kBlue);
  gklpt->Draw("P");
  legend->Draw();
c6->SetLogx();
c6->Update();
char* file6= Form("K0L_Run%i.png",IDNUM);
c6->SaveAs(file6);

//k0s
TCanvas *c7 =new TCanvas("c7","K0S",200,10,700,500);
for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  y[i]=Stuff[i][0][7].VALUE;
  ex[i]=0;
  THINGY=Stuff[i][0][7].entries;
  ey[i]=(Stuff[i][0][7].stddev)/sqrt(THINGY);
}
TGraphErrors* gkset=new TGraphErrors(n,x,y,ex,ey);
  gkset->SetTitle("K0S avg; SNN(GeV); (GeV)");
  gkset->SetMarkerStyle(20);
  gkset->SetMarkerSize(.5);
  gkset->SetMarkerColor(kRed);
  gkset->Draw("AP");

for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  y[i]=Stuff[i][1][7].VALUE;
  ex[i]=0;
  THINGY=Stuff[i][1][7].entries;
  ey[i]=(Stuff[i][1][7].stddev)/sqrt(THINGY);
}
TGraphErrors* gkspt=new TGraphErrors(n,x,y,ex,ey);
  gkspt->SetMarkerStyle(20);
  gkspt->SetMarkerSize(.5);
  gkspt->SetMarkerColor(kBlue);
  gkspt->Draw("P");
  legend->Draw();
c7->SetLogx();
c7->Update();
char* file7= Form("K0S_Run%i.png",IDNUM);
c7->SaveAs(file7);

//eta
TCanvas *c8 =new TCanvas("c8","eta",200,10,700,500);
for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  y[i]=Stuff[i][0][8].VALUE;
  ex[i]=0;
  THINGY=Stuff[i][0][8].entries;
  ey[i]=(Stuff[i][0][8].stddev)/sqrt(THINGY);
}
TGraphErrors* geet=new TGraphErrors(n,x,y,ex,ey);
  geet->SetTitle("eta avg; SNN(GeV); (GeV)");
  geet->SetMarkerStyle(20);
  geet->SetMarkerSize(.5);
  geet->SetMarkerColor(kRed);
  geet->Draw("AP");

for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  y[i]=Stuff[i][1][8].VALUE;
  ex[i]=0;
  THINGY=Stuff[i][1][8].entries;
  ey[i]=(Stuff[i][1][8].stddev)/sqrt(THINGY);
}
TGraphErrors* gept=new TGraphErrors(n,x,y,ex,ey);
  gept->SetMarkerStyle(20);
  gept->SetMarkerSize(.5);
  gept->SetMarkerColor(kBlue);
  gept->Draw("P");
  legend->Draw();
c8->SetLogx();
c8->Update();
char* file8= Form("eta_Run%i.png",IDNUM);
c8->SaveAs(file8);

//omega
TCanvas *c9 =new TCanvas("c9","omega",200,10,700,500);
for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  y[i]=Stuff[i][0][9].VALUE;
  ex[i]=0;
  THINGY=Stuff[i][0][9].entries;
  ey[i]=(Stuff[i][0][9].stddev)/sqrt(THINGY);
}
TGraphErrors* goet=new TGraphErrors(n,x,y,ex,ey);
  goet->SetTitle("omega avg; SNN(GeV); (GeV)");
  goet->SetMarkerStyle(20);
  goet->SetMarkerSize(.5);
  goet->SetMarkerColor(kRed);
  goet->Draw("AP");

for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  y[i]=Stuff[i][1][9].VALUE;
  ex[i]=0;
  THINGY=Stuff[i][1][9].entries;
  ey[i]=(Stuff[i][1][9].stddev)/sqrt(THINGY);
}
TGraphErrors* gopt=new TGraphErrors(n,x,y,ex,ey);
  gopt->SetMarkerStyle(20);
  gopt->SetMarkerSize(.5);
  gopt->SetMarkerColor(kBlue);
  gopt->Draw("P");
  legend->Draw();
c9->SetLogx();
c9->Update();
char* file9= Form("omega_Run%i.png",IDNUM);
c9->SaveAs(file9);

//l0
TCanvas *c10 =new TCanvas("c10","Lambda 0",200,10,700,500);
for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  y[i]=Stuff[i][0][10].VALUE;
  ex[i]=0;
  THINGY=Stuff[i][0][10].entries;
  ey[i]=(Stuff[i][0][10].stddev)/sqrt(THINGY);
}
TGraphErrors* glet=new TGraphErrors(n,x,y,ex,ey);
  glet->SetTitle("Lambda0 avg; SNN(GeV); (GeV)");
  glet->SetMarkerStyle(20);
  glet->SetMarkerSize(.5);
  glet->SetMarkerColor(kRed);
  glet->Draw("AP");

for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  y[i]=Stuff[i][1][10].VALUE;
  ex[i]=0;
  THINGY=Stuff[i][1][10].entries;
  ey[i]=(Stuff[i][1][10].stddev)/sqrt(THINGY);
}
TGraphErrors* glpt=new TGraphErrors(n,x,y,ex,ey);
  glpt->SetMarkerStyle(20);
  glpt->SetMarkerSize(.5);
  glpt->SetMarkerColor(kBlue);
  glpt->Draw("P");
  legend->Draw();
c10->SetLogx();
c10->Update();
char* file10= Form("l0_Run%i.png",IDNUM);
c10->SaveAs(file10);

//lb0
TCanvas *c11 =new TCanvas("c11","LambdaBar 0",200,10,700,500);
for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  y[i]=Stuff[i][0][11].VALUE;
  ex[i]=0;
  THINGY=Stuff[i][0][11].entries;
  ey[i]=(Stuff[i][0][11].stddev)/sqrt(THINGY);
}
TGraphErrors* glbet=new TGraphErrors(n,x,y,ex,ey);
  glbet->SetTitle("LambdaBar0 avg; SNN(GeV); (GeV)");
  glbet->SetMarkerStyle(20);
  glbet->SetMarkerSize(.5);
  glbet->SetMarkerColor(kRed);
  glbet->Draw("AP");

for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  y[i]=Stuff[i][1][11].VALUE;
  ex[i]=0;
  THINGY=Stuff[i][1][11].entries;
  ey[i]=(Stuff[i][1][11].stddev)/sqrt(THINGY);
}
TGraphErrors* glbpt=new TGraphErrors(n,x,y,ex,ey);
  glbpt->SetMarkerStyle(20);
  glbpt->SetMarkerSize(.5);
  glbpt->SetMarkerColor(kBlue);
  glbpt->Draw("P");
  legend->Draw();
c11->SetLogx();
c11->Update();
char* file11= Form("lb0_Run%i.png",IDNUM);
c11->SaveAs(file11);

//p
TCanvas *c12 =new TCanvas("c12","proton",200,10,700,500);
for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  y[i]=Stuff[i][0][12].VALUE;
  ex[i]=0;
  THINGY=Stuff[i][0][12].entries;
  ey[i]=(Stuff[i][0][12].stddev)/sqrt(THINGY);
}
TGraphErrors* gpet=new TGraphErrors(n,x,y,ex,ey);
  gpet->SetTitle("proton avg; SNN(GeV); (GeV)");
  gpet->SetMarkerStyle(20);
  gpet->SetMarkerSize(.5);
  gpet->SetMarkerColor(kRed);
  gpet->Draw("AP");

for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  y[i]=Stuff[i][1][12].VALUE;
  ex[i]=0;
  THINGY=Stuff[i][1][12].entries;
  ey[i]=(Stuff[i][1][12].stddev)/sqrt(THINGY);
}
TGraphErrors* gppt=new TGraphErrors(n,x,y,ex,ey);
  gppt->SetMarkerStyle(20);
  gppt->SetMarkerSize(.5);
  gppt->SetMarkerColor(kBlue);
  gppt->Draw("P");
  legend->Draw();
c12->SetLogx();
c12->Update();
char* file12= Form("proton_Run%i.png",IDNUM);
c12->SaveAs(file12);

//n
TCanvas *c13 =new TCanvas("c13","neutron",200,10,700,500);
for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  y[i]=Stuff[i][0][13].VALUE;
  ex[i]=0;
  THINGY=Stuff[i][0][13].entries;
  ey[i]=(Stuff[i][0][13].stddev)/sqrt(THINGY);
}
TGraphErrors* gnet=new TGraphErrors(n,x,y,ex,ey);
  gnet->SetTitle("neutron avg; SNN(GeV); (GeV)");
  gnet->SetMarkerStyle(20);
  gnet->SetMarkerSize(.5);
  gnet->SetMarkerColor(kRed);
  gnet->Draw("AP");

for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  y[i]=Stuff[i][1][13].VALUE;
  ex[i]=0;
  THINGY=Stuff[i][1][13].entries;
  ey[i]=(Stuff[i][1][13].stddev)/sqrt(THINGY);
}
TGraphErrors* gnpt=new TGraphErrors(n,x,y,ex,ey);
  gnpt->SetMarkerStyle(20);
  gnpt->SetMarkerSize(.5);
  gnpt->SetMarkerColor(kBlue);
  gnpt->Draw("P");
  legend->Draw();
c13->SetLogx();
c13->Update();
char* file13= Form("neutron_Run%i.png",IDNUM);
c13->SaveAs(file13);

//pb
TCanvas *c14 =new TCanvas("c14","protonBar",200,10,700,500);
for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  y[i]=Stuff[i][0][14].VALUE;
  ex[i]=0;
  THINGY=Stuff[i][0][14].entries;
  ey[i]=(Stuff[i][0][14].stddev)/sqrt(THINGY);
}
TGraphErrors* gpbet=new TGraphErrors(n,x,y,ex,ey);
  gpbet->SetTitle("protonBar avg; SNN(GeV); (GeV)");
  gpbet->SetMarkerStyle(20);
  gpbet->SetMarkerSize(.5);
  gpbet->SetMarkerColor(kRed);
  gpbet->Draw("AP");

for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  y[i]=Stuff[i][1][14].VALUE;
  ex[i]=0;
  THINGY=Stuff[i][1][14].entries;
  ey[i]=(Stuff[i][1][14].stddev)/sqrt(THINGY);
}
TGraphErrors* gpbpt=new TGraphErrors(n,x,y,ex,ey);
  gpbpt->SetMarkerStyle(20);
  gpbpt->SetMarkerSize(.5);
  gpbpt->SetMarkerColor(kBlue);
  gpbpt->Draw("P");
  legend->Draw();
c14->SetLogx();
c14->Update();
char* file14= Form("protonBar_Run%i.png",IDNUM);
c14->SaveAs(file14);

//nb
TCanvas *c15 =new TCanvas("c15","neutronBar",200,10,700,500);
for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  y[i]=Stuff[i][0][15].VALUE;
  ex[i]=0;
  THINGY=Stuff[i][0][15].entries;
  ey[i]=(Stuff[i][0][15].stddev)/sqrt(THINGY);
}
TGraphErrors* gnbet=new TGraphErrors(n,x,y,ex,ey);
  gnbet->SetTitle("neutronBar avg; SNN(GeV); (GeV)");
  gnbet->SetMarkerStyle(20);
  gnbet->SetMarkerSize(.5);
  gnbet->SetMarkerColor(kRed);
  gnbet->Draw("AP");

for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  y[i]=Stuff[i][1][15].VALUE;
  ex[i]=0;
  THINGY=Stuff[i][1][15].entries;
  ey[i]=(Stuff[i][1][15].stddev)/sqrt(THINGY);
}
TGraphErrors* gnbpt=new TGraphErrors(n,x,y,ex,ey);
  gnbpt->SetMarkerStyle(20);
  gnbpt->SetMarkerSize(.5);
  gnbpt->SetMarkerColor(kBlue);
  gnbpt->Draw("P");
  legend->Draw();
c15->SetLogx();
c15->Update();
char* file15= Form("neutronBar_Run%i.png",IDNUM);
c15->SaveAs(file15);
*/} //extra Graphs
cout<<"plotting comparisions"<<endl;
Double_t blarg;

//pi0 test
/******************************************************************************
The following are the checks required our assumptions.
The following applies to all charts that follow:

colors:
Black = pions
Blue = kaons
Red = p/n

markers:
Circle = neutral over charged

Observe the different versions of markers and their purpose
******************************************************************************/
//pi
TCanvas *c16 =new TCanvas("c16","pi0 check",0,645,1280,745);
for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  blarg=((Stuff[i][0][1].VALUE+Stuff[i][0][2].VALUE)/2)/(Stuff[i][0][3].VALUE);
  y[i]=blarg;
  ex[i]=0;
  blarg=((Stuff[i][0][1].stddev+Stuff[i][0][2].stddev)/2)/(Stuff[i][0][3].stddev);
  THINGY=Stuff[i][0][1].entries;
  ey[i]=blarg/sqrt(THINGY);
}
TGraphErrors* gc1=new TGraphErrors(n,x,y,ex,ey);
  gc1->SetTitle("Pi0 check (ET); SNN(GeV); Ratio");
  gc1->SetMarkerStyle(20);
  gc1->SetMarkerSize(1.5);
  gc1->SetMarkerColor(kBlack);
  gc1->SetLineColor(kBlack);
  legend2->AddEntry(gc1,"Mean(Pi+,Pi-)/Pi0","p");
  gc1->Draw("AP");
  gc1->GetXaxis()->SetLimits(6,300);
  gc1->GetYaxis()->SetRangeUser(.7,1.3);
  gc1->Draw("AP");

for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i]*1.033;
  blarg=(Stuff[i][0][1].VALUE)/(Stuff[i][0][3].VALUE);
  y[i]=blarg;
  ex[i]=0;
  blarg=(Stuff[i][0][1].stddev)/(Stuff[i][0][3].stddev);
  THINGY=Stuff[i][0][1].entries;
  ey[i]=blarg/sqrt(THINGY);
}
TGraphErrors* gc10=new TGraphErrors(n,x,y,ex,ey);
  gc10->SetMarkerStyle(41);
  gc10->SetMarkerSize(1.5);
  gc10->SetMarkerColor(kBlack);
  gc10->SetLineColor(kBlack);
  legend2->AddEntry(gc10,"Pi+/Pi0","p");
  gc10->Draw("P");

for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i]*1.066;
  blarg=(Stuff[i][0][2].VALUE)/(Stuff[i][0][3].VALUE);
  y[i]=blarg;
  ex[i]=0;
  blarg=(Stuff[i][0][2].stddev)/(Stuff[i][0][3].stddev);
  THINGY=Stuff[i][0][2].entries;
  ey[i]=blarg/sqrt(THINGY);
}
TGraphErrors* gc15=new TGraphErrors(n,x,y,ex,ey);
  gc15->SetMarkerStyle(22);
  gc15->SetMarkerSize(1.5);
  gc15->SetMarkerColor(kBlack);
  gc15->SetLineColor(kBlack);
  legend2->AddEntry(gc15,"Pi-/Pi0","p");
  gc15->Draw("P");
/*
for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i]*1.1;
  blarg=(Stuff[i][0][1].VALUE+Stuff[i][0][2].VALUE)/(Stuff[i][0][3].VALUE);
  y[i]=blarg;
  ex[i]=0;
  blarg=(Stuff[i][0][1].stddev+Stuff[i][0][2].stddev)/(Stuff[i][0][3].stddev);
  THINGY=Stuff[i][0][1].entries;
  ey[i]=blarg/sqrt(THINGY);
}
TGraphErrors* gc16=new TGraphErrors(n,x,y,ex,ey);
  gc16->SetMarkerStyle(21);
  gc16->SetMarkerSize(1.5);
  gc16->SetMarkerColor(kBlack);
  gc16->SetLineColor(kBlack);
  legend2->AddEntry(gc16,"Pi+,pi-/Pi0","p");
  gc16->Draw("P");
*/
legend2->Draw();
Int_t n1=3;
Int_t x1[3]={0,100,15000};
Int_t y1[3]={1,1,1};
TGraph* gl=new TGraph(n1,x1,y1);
  gl->SetLineColor(kBlack);
  gl->Draw("L");
c16->SetLogx();
c16->Update();
char* file16= Form("p%c%iPi0Check.png",ty,IDNUM);
c16->SaveAs(file16);

//kaon check
TCanvas *c17 =new TCanvas("c17","kaon checks",0,645,1280,745);
for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  y[i]=(Stuff[i][0][4].VALUE)/(Stuff[i][0][5].VALUE);
  ex[i]=0;
  THINGY=Stuff[i][0][4].entries;
  ey[i]=((Stuff[i][0][4].stddev)/(Stuff[i][0][5].stddev))/sqrt(THINGY);
}
TGraphErrors* gc2=new TGraphErrors(n,x,y,ex,ey);
  gc2->SetTitle("Kaon Checks avg; SNN(GeV); Ratio");
  gc2->SetMarkerStyle(22);
  gc2->SetMarkerSize(1.5);
  gc2->SetMarkerColor(kBlue);
  gc2->SetLineColor(kBlack);
  legend3->AddEntry(gc2,"K+/K-","p");
  gc2->Draw("AP");
  gc2->GetXaxis()->SetLimits(6,300);
  gc2->GetYaxis()->SetRangeUser(.6,1.4);
  gc2->Draw("AP");

for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i]*1.033;
  y[i]=(Stuff[i][0][4].VALUE)/(Stuff[i][0][6].VALUE);
  ex[i]=0;
  THINGY=Stuff[i][0][4].entries;
  ey[i]=((Stuff[i][0][4].stddev)/(Stuff[i][0][6].stddev))/sqrt(THINGY);
}
TGraphErrors* gc3=new TGraphErrors(n,x,y,ex,ey);
  gc3->SetMarkerStyle(21);
  gc3->SetMarkerSize(1.5);
  gc3->SetMarkerColor(kBlue);
  gc3->SetLineColor(kBlack);
  legend3->AddEntry(gc3,"K+/K0L","p");
  gc3->Draw("P");

for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i]*1.066;
  y[i]=(Stuff[i][0][4].VALUE)/(Stuff[i][0][7].VALUE);
  ex[i]=0;
  THINGY=Stuff[i][0][4].entries;
  ey[i]=((Stuff[i][0][4].stddev)/(Stuff[i][0][7].stddev))/sqrt(THINGY);
}
TGraphErrors* gc4=new TGraphErrors(n,x,y,ex,ey);
  gc4->SetMarkerStyle(33);
  gc4->SetMarkerSize(2);
  gc4->SetMarkerColor(kBlue);
  gc4->SetLineColor(kBlack);
  legend3->AddEntry(gc4,"K+/K0S","p");
  gc4->Draw("P");

for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i]*1.1;
  y[i]=(Stuff[i][0][4].VALUE+Stuff[i][0][5].VALUE)/(Stuff[i][0][6].VALUE+Stuff[i][0][7].VALUE);
  ex[i]=0;
  THINGY=Stuff[i][0][5].entries;
  ey[i]=((Stuff[i][0][4].stddev+Stuff[i][0][5].stddev)/(Stuff[i][0][6].stddev+Stuff[i][0][7].stddev))/sqrt(THINGY);
}
TGraphErrors* gc5=new TGraphErrors(n,x,y,ex,ey);
  gc5->SetMarkerStyle(20);
  gc5->SetMarkerSize(1.5);
  gc5->SetMarkerColor(kBlue);
  gc5->SetLineColor(kBlack);
  legend3->AddEntry(gc5,"Kp+Km/K0L+K0S","p");
  gc5->Draw("P");

for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i]*1.133;
  y[i]=(Stuff[i][0][6].VALUE)/(Stuff[i][0][7].VALUE);
  ex[i]=0;
  THINGY=Stuff[i][0][6].entries;
  ey[i]=((Stuff[i][0][6].stddev)/(Stuff[i][0][7].stddev))/sqrt(THINGY);
}
TGraphErrors* gc7=new TGraphErrors(n,x,y,ex,ey);
  gc7->SetMarkerStyle(41);
  gc7->SetMarkerSize(1.5);
  gc7->SetMarkerColor(kBlue);
  gc7->SetLineColor(kBlack);
  legend3->AddEntry(gc7,"K0L/K0S","p");
  gc7->Draw("P");
  legend3->Draw();
gl->Draw("L");
c17->SetLogx();
c17->Update();
char* file17= Form("p%c%iKaonChecks.png",ty,IDNUM);
c17->SaveAs(file17);

//pi0 test
TCanvas *c18 =new TCanvas("c18","p-n checks",0,645,1280,745);
for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  blarg=((Stuff[i][0][12].VALUE)/(Stuff[i][0][13].VALUE));
  y[i]=blarg;
  ex[i]=0;
  blarg=((Stuff[i][0][12].stddev)/(Stuff[i][0][13].stddev));
  THINGY=Stuff[i][0][12].entries;
  ey[i]=blarg/sqrt(THINGY);
}
TGraphErrors* gc8=new TGraphErrors(n,x,y,ex,ey);
  gc8->SetTitle("Proton and Neutron Checks; SNN(GeV); Ratio");
  gc8->SetMarkerStyle(22);
  gc8->SetMarkerSize(1.5);
  gc8->SetMarkerColor(kRed);
  gc8->SetLineColor(kBlack);
  legend4->AddEntry(gc8,"p/n","p");
  gc8->Draw("AP");
  gc8->GetXaxis()->SetLimits(6,300);
  gc8->GetYaxis()->SetRangeUser(.6,1.4);
  gc8->Draw("AP");

for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i]*1.033;
  blarg=((Stuff[i][0][14].VALUE)/(Stuff[i][0][15].VALUE));
  y[i]=blarg;
  ex[i]=0;
  blarg=((Stuff[i][0][14].stddev)/(Stuff[i][0][15].stddev));
  THINGY=Stuff[i][0][14].entries;
  ey[i]=blarg/sqrt(THINGY);
}
TGraphErrors* gc9=new TGraphErrors(n,x,y,ex,ey);
  gc9->SetMarkerStyle(41);
  gc9->SetMarkerSize(1.5);
  gc9->SetMarkerColor(kRed);
  gc9->SetLineColor(kBlack);
  legend4->AddEntry(gc9,"pBar/nBar","p");
  gc9->Draw("P");

for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i]*1.066;
  blarg=((Stuff[i][0][12].VALUE+Stuff[i][0][14].VALUE)/(Stuff[i][0][13].VALUE+Stuff[i][0][15].VALUE));
  y[i]=blarg;
  ex[i]=0;
  blarg=((Stuff[i][0][12].stddev+Stuff[i][0][14].stddev)/(Stuff[i][0][13].stddev+Stuff[i][0][15].stddev));
  THINGY=Stuff[i][0][12].entries;
  ey[i]=blarg/sqrt(THINGY);
}
TGraphErrors* gc11=new TGraphErrors(n,x,y,ex,ey);
  gc11->SetMarkerStyle(20);
  gc11->SetMarkerSize(1.5);
  gc11->SetMarkerColor(kRed);
  gc11->SetLineColor(kBlack);
  legend4->AddEntry(gc11,"p+pbar/n+nbar","p");
  gc11->Draw("P");
  legend4->Draw();
gl->Draw("L");
c18->SetLogx();
c18->Update();
char* file18= Form("p%c%ipncheck.png",ty,IDNUM);
c18->SaveAs(file18);

{/*//vsALL
TCanvas *c19 =new TCanvas("c19","Particle Et vs ETAll",0,645,1280,745);
for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  blarg=1/((Stuff[i][0][1].VALUE+Stuff[i][0][2].VALUE+Stuff[i][0][3].VALUE)/(Stuff[i][0][0].VALUE));
  y[i]=blarg;
  ex[i]=0;
  blarg=1/((Stuff[i][0][1].stddev+Stuff[i][0][2].stddev+Stuff[i][0][3].stddev)/(Stuff[i][0][0].stddev));
  THINGY=Stuff[i][0][1].entries;
  ey[i]=blarg/sqrt(THINGY);
}
TGraphErrors* gc12=new TGraphErrors(n,x,y,ex,ey);
  gc12->SetTitle("Species ET vs All; SNN(GeV); Ratio");
  gc12->SetMarkerStyle(23);
  gc12->SetMarkerSize(1.5);
  gc12->SetMarkerColor(kBlack);
  gc12->SetLineColor(kBlack);
  legend5->AddEntry(gc12,"Pion","p");
  gc12->SetMinimum(0);
  gc12->SetMaximum(1.2);
  gc12->Draw("AP");

for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i]*1.033;
  blarg=1/((Stuff[i][0][4].VALUE+Stuff[i][0][5].VALUE+Stuff[i][0][6].VALUE+Stuff[i][0][7].VALUE)/(Stuff[i][0][0].VALUE));
  y[i]=blarg;
  ex[i]=0;
  blarg=1/((Stuff[i][0][4].stddev+Stuff[i][0][5].stddev+Stuff[i][0][6].stddev+Stuff[i][0][7].stddev)/(Stuff[i][0][0].stddev));
  THINGY=Stuff[i][0][1].entries;
  ey[i]=blarg/sqrt(THINGY);
}
TGraphErrors* gc14=new TGraphErrors(n,x,y,ex,ey);
  gc14->SetMarkerStyle(23);
  gc14->SetMarkerSize(1.5);
  gc14->SetMarkerColor(kBlue);
  gc14->SetLineColor(kBlack);
  legend5->AddEntry(gc14,"Kaon","p");
  gc14->Draw("P");

for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i]*1.066;
  blarg=1/((Stuff[i][0][12].VALUE+Stuff[i][0][13].VALUE+Stuff[i][0][14].VALUE+Stuff[i][0][15].VALUE)/(Stuff[i][0][0].VALUE));
  y[i]=blarg;
  ex[i]=0;
  blarg=1/((Stuff[i][0][12].stddev+Stuff[i][0][13].stddev+Stuff[i][0][14].stddev+Stuff[i][0][15].stddev)/(Stuff[i][0][0].stddev));
  THINGY=Stuff[i][0][1].entries;
  ey[i]=blarg/sqrt(THINGY);
}
TGraphErrors* gc13=new TGraphErrors(n,x,y,ex,ey);
  gc13->SetMarkerStyle(23);
  gc13->SetMarkerSize(1.5);
  gc13->SetMarkerColor(kRed);
  gc13->SetLineColor(kBlack);
  legend5->AddEntry(gc13,"protons and neutrons","p");
  gc13->Draw("P");
  legend5->Draw();
  gl->Draw("L");
c19->SetLogx();
c19->Update();
char* file19= Form("ALLCheck_Run%i.png",IDNUM);
c19->SaveAs(file19);
*/}


TLegend* legend6=new TLegend(.719092,.115278,.895931,.256944);
    legend6->SetTextFont(72);
    legend6->SetTextSize(.03);
    legend6->SetFillColor(0);
TLegend* legend7=new TLegend(.712833,.709722,.895149,.893056);
    legend7->SetTextFont(72);
    legend7->SetTextSize(.03);
    legend7->SetFillColor(0);
TLegend* legend8=new TLegend(.743349,.781944,.897496,.895833);
    legend8->SetTextFont(72);
    legend8->SetTextSize(.03);
    legend8->SetFillColor(0);

TCanvas *c80 =new TCanvas("c80","pi0 check txt",0,645,1280,745);
for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  blarg=((StuffTXT[i][0][0].V+StuffTXT[i][0][1].V)/2)/(StuffTXT[i][0][2].V);
  y[i]=blarg;
  ex[i]=0;
  blarg=((StuffTXT[i][0][0].s+StuffTXT[i][0][1].s)/2)/(StuffTXT[i][0][2].s);
  THINGY=StuffTXT[i][0][0].n;
  ey[i]=blarg/sqrt(THINGY);
}
TGraphErrors* gc70=new TGraphErrors(n,x,y,ex,ey);
  gc70->SetTitle("Pi0 check (from .TXT); SNN(GeV); Ratio");
  gc70->SetMarkerStyle(20);
  gc70->SetMarkerSize(1.5);
  gc70->SetMarkerColor(kBlack);
  gc70->SetLineColor(kBlack);
  legend6->AddEntry(gc70,"Mean(Pi+,Pi-)/Pi0","p");
  gc70->Draw("AP");
  gc70->GetXaxis()->SetLimits(6,300);
  gc70->GetYaxis()->SetRangeUser(.7,1.3);
  gc70->Draw("AP");

for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i]*1.033;
  blarg=(StuffTXT[i][0][0].V)/(StuffTXT[i][0][2].V);
  y[i]=blarg;
  ex[i]=0;
  blarg=(StuffTXT[i][0][0].s)/(StuffTXT[i][0][2].s);
  THINGY=StuffTXT[i][0][0].n;
  ey[i]=blarg/sqrt(THINGY);
}
TGraphErrors* gc71=new TGraphErrors(n,x,y,ex,ey);
  gc71->SetMarkerStyle(41);
  gc71->SetMarkerSize(1.5);
  gc71->SetMarkerColor(kBlack);
  gc71->SetLineColor(kBlack);
  legend6->AddEntry(gc71,"Pi+/Pi0","p");
  gc71->Draw("P");

for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i]*1.066;
  blarg=(StuffTXT[i][0][1].V)/(StuffTXT[i][0][2].V);
  y[i]=blarg;
  ex[i]=0;
  blarg=(StuffTXT[i][0][1].s)/(StuffTXT[i][0][2].s);
  THINGY=StuffTXT[i][0][1].n;
  ey[i]=blarg/sqrt(THINGY);
}
TGraphErrors* gc72=new TGraphErrors(n,x,y,ex,ey);
  gc72->SetMarkerStyle(22);
  gc72->SetMarkerSize(1.5);
  gc72->SetMarkerColor(kBlack);
  gc72->SetLineColor(kBlack);
  legend6->AddEntry(gc72,"Pi-/Pi0","p");
  gc72->Draw("P");
legend6->Draw();
TGraph* gl73=new TGraph(n1,x1,y1);
  gl73->SetLineColor(kBlack);
  gl73->Draw("L");
c80->SetLogx();
c80->Update();
char* file80= Form("p%c%iPi0CheckFROMTXT.png",ty,IDNUM);
c80->SaveAs(file80);




/////KAONS
TCanvas *c81 =new TCanvas("c81","kaon checks",0,645,1280,745);
for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  y[i]=(StuffTXT[i][0][3].V)/(StuffTXT[i][0][4].V);
  ex[i]=0;
  THINGY=StuffTXT[i][0][3].n;
  ey[i]=((StuffTXT[i][0][3].s)/(StuffTXT[i][0][4].s))/sqrt(THINGY);
}
TGraphErrors* gc74=new TGraphErrors(n,x,y,ex,ey);
  gc74->SetTitle("Kaon Checks avg (from .TXT); SNN(GeV); Ratio");
  gc74->SetMarkerStyle(22);
  gc74->SetMarkerSize(1.5);
  gc74->SetMarkerColor(kBlue);
  gc74->SetLineColor(kBlack);
  legend7->AddEntry(gc74,"K+/K-","p");
  gc74->Draw("AP");
  gc74->GetXaxis()->SetLimits(6,300);
  gc74->GetYaxis()->SetRangeUser(.6,1.4);
  gc74->Draw("AP");

for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i]*1.033;
  y[i]=(StuffTXT[i][0][3].V)/(StuffTXT[i][0][5].V);
  ex[i]=0;
  THINGY=StuffTXT[i][0][3].n;
  ey[i]=((StuffTXT[i][0][3].s)/(StuffTXT[i][0][5].s))/sqrt(THINGY);
}
TGraphErrors* gc75=new TGraphErrors(n,x,y,ex,ey);
  gc75->SetMarkerStyle(21);
  gc75->SetMarkerSize(1.5);
  gc75->SetMarkerColor(kBlue);
  gc75->SetLineColor(kBlack);
  legend7->AddEntry(gc75,"K+/K0L","p");
  gc75->Draw("P");

for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i]*1.066;
  y[i]=(StuffTXT[i][0][3].V)/(StuffTXT[i][0][6].V);
  ex[i]=0;
  THINGY=StuffTXT[i][0][3].n;
  ey[i]=((StuffTXT[i][0][3].s)/(StuffTXT[i][0][6].s))/sqrt(THINGY);
}
TGraphErrors* gc76=new TGraphErrors(n,x,y,ex,ey);
  gc76->SetMarkerStyle(33);
  gc76->SetMarkerSize(2);
  gc76->SetMarkerColor(kBlue);
  gc76->SetLineColor(kBlack);
  legend7->AddEntry(gc76,"K+/K0S","p");
  gc76->Draw("P");

for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i]*1.1;
  y[i]=(StuffTXT[i][0][3].V+StuffTXT[i][0][4].V)/(StuffTXT[i][0][5].V+StuffTXT[i][0][6].V);
  ex[i]=0;
  THINGY=StuffTXT[i][0][4].n;
  ey[i]=((StuffTXT[i][0][3].s+StuffTXT[i][0][4].s)/(StuffTXT[i][0][5].s+StuffTXT[i][0][6].s))/sqrt(THINGY);
}
TGraphErrors* gc77=new TGraphErrors(n,x,y,ex,ey);
  gc77->SetMarkerStyle(20);
  gc77->SetMarkerSize(1.5);
  gc77->SetMarkerColor(kBlue);
  gc77->SetLineColor(kBlack);
  legend7->AddEntry(gc77,"Kp+Km/K0L+K0S","p");
  gc77->Draw("P");

for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i]*1.133;
  y[i]=(StuffTXT[i][0][5].V)/(StuffTXT[i][0][6].V);
  ex[i]=0;
  THINGY=StuffTXT[i][0][5].n;
  ey[i]=((StuffTXT[i][0][5].s)/(StuffTXT[i][0][6].s))/sqrt(THINGY);
}
TGraphErrors* gc78=new TGraphErrors(n,x,y,ex,ey);
  gc78->SetMarkerStyle(41);
  gc78->SetMarkerSize(1.5);
  gc78->SetMarkerColor(kBlue);
  gc78->SetLineColor(kBlack);
  legend7->AddEntry(gc78,"K0L/K0S","p");
  gc78->Draw("P");
  legend3->Draw();
gl73->Draw("L");
c81->SetLogx();
c81->Update();
char* file81= Form("p%c%iKaonChecksFROMTXT.png",ty,IDNUM);
c81->SaveAs(file81);

//proton/neutron
TCanvas *c82 =new TCanvas("c82","p-n checks",0,645,1280,745);
for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  blarg=((StuffTXT[i][0][9].V)/(StuffTXT[i][0][11].V));
  y[i]=blarg;
  ex[i]=0;
  blarg=((StuffTXT[i][0][9].s)/(StuffTXT[i][0][11].s));
  THINGY=StuffTXT[i][0][9].n;
  ey[i]=blarg/sqrt(THINGY);
}
TGraphErrors* gc79=new TGraphErrors(n,x,y,ex,ey);
  gc79->SetTitle("Proton and Neutron Checks (from .TXT); SNN(GeV); Ratio");
  gc79->SetMarkerStyle(22);
  gc79->SetMarkerSize(1.5);
  gc79->SetMarkerColor(kRed);
  gc79->SetLineColor(kBlack);
  legend8->AddEntry(gc79,"p/n","p");
  gc79->Draw("AP");
  gc79->GetXaxis()->SetLimits(6,300);
  gc79->GetYaxis()->SetRangeUser(.6,1.4);
  gc79->Draw("AP");

for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i]*1.033;
  blarg=((StuffTXT[i][0][10].V)/(StuffTXT[i][0][12].V));
  y[i]=blarg;
  ex[i]=0;
  blarg=((StuffTXT[i][0][10].s)/(StuffTXT[i][0][12].s));
  THINGY=StuffTXT[i][0][10].n;
  ey[i]=blarg/sqrt(THINGY);
}
TGraphErrors* gc80=new TGraphErrors(n,x,y,ex,ey);
  gc80->SetMarkerStyle(41);
  gc80->SetMarkerSize(1.5);
  gc80->SetMarkerColor(kRed);
  gc80->SetLineColor(kBlack);
  legend8->AddEntry(gc80,"pBar/nBar","p");
  gc80->Draw("P");

for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i]*1.066;
  blarg=((StuffTXT[i][0][9].V+StuffTXT[i][0][10].V)/(StuffTXT[i][0][11].V+StuffTXT[i][0][12].V));
  y[i]=blarg;
  ex[i]=0;
  blarg=((StuffTXT[i][0][9].s+StuffTXT[i][0][10].s)/(StuffTXT[i][0][11].s+StuffTXT[i][0][12].s));
  THINGY=StuffTXT[i][0][9].n;
  ey[i]=blarg/sqrt(THINGY);
}
TGraphErrors* gc81=new TGraphErrors(n,x,y,ex,ey);
  gc81->SetMarkerStyle(20);
  gc81->SetMarkerSize(1.5);
  gc81->SetMarkerColor(kRed);
  gc81->SetLineColor(kBlack);
  legend8->AddEntry(gc81,"p+pbar/n+nbar","p");
  gc81->Draw("P");
  legend8->Draw();
gl73->Draw("L");
c82->SetLogx();
c82->Update();
char* file82= Form("p%c%ipncheckFROMTXT.png",ty,IDNUM);
c82->SaveAs(file82);





cout<<"program complete"<<endl<<endl;
return 0;
}

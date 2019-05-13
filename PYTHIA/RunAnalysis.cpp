/******************************************************************************
This program takes data file "PythiaData.txt". outputs the averages in a
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
  Double_t average;
  Double_t stddev;
  Double_t entries;
};

int RunAnalysis(){
  Int_t IDNUM=1;
  Int_t a=0;//type count 0-1
  Int_t b=0;//part count 0-15
  Int_t c=0;//info count 0-2
  vector<Int_t> SNN {7,11,19,27,39,62,130,200,900,2760,5020,7000,8000,13000,14000};
  vector<double> SNNACTUAL {7.7,11.6,19.6,27,39,62.4,130,200,900,2760,5020,7000,8000,13000,14000};
  /*****************************************************************************
  THE FOLLOWING is a reference for counters a, b ,and c
  a{"ET","PT"};
  b{"ALL","PiP","PiM","Pi0","KaP","KaM","K0L","K0S","ETA","OMG","La0","LB0","PRO","NEU","PRB","NEB"};
  c{"AVG","S","N"};
  *****************************************************************************/
  vector<vector<vector<dat>>> Stuff(15, vector<vector<dat>>(2, vector<dat>(16)));

  //averages for each histo
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

  cout<<"Select Run Number (int): "<<endl;
  cin>>IDNUM;

  cout<<"preparing data for run: "<<IDNUM<<endl<<endl;
  for (int x=0;x<15;x++){
    char* temp = Form("hh%iGeVoutfile%i.root",SNN[x],IDNUM);
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

    a=0;
    b=0;
    ETALLAVG=hETAll->GetMean(1);
    Stuff[x][a][b].average=ETALLAVG; //DEBUG
    b++;
    ETPiPAVG=hETPiPlus->GetMean(1);
    Stuff[x][a][b].average=ETPiPAVG;
    b++;
    ETPiMAVG=hETPiMinus->GetMean(1);
    Stuff[x][a][b].average=ETPiMAVG;
    b++;
    ETPi0AVG=hETPi0->GetMean(1);
    Stuff[x][a][b].average=ETPi0AVG;
    b++;
    ETKaPAVG=hETKPlus->GetMean(1);
    Stuff[x][a][b].average=ETKaPAVG;
    b++;
    ETKaMAVG=hETKMinus->GetMean(1);
    Stuff[x][a][b].average=ETKaMAVG;
    b++;
    ETK0LAVG=hETKL->GetMean(1);
    Stuff[x][a][b].average=ETK0LAVG;
    b++;
    ETK0SAVG=hETKS->GetMean(1);
    Stuff[x][a][b].average=ETK0SAVG;
    b++;
    ETETAAVG=hETEta->GetMean(1);
    Stuff[x][a][b].average=ETETAAVG;
    b++;
    ETOMGAVG=hETomega->GetMean(1);
    Stuff[x][a][b].average=ETOMGAVG;
    b++;
    ETLa0AVG=hETLambda0->GetMean(1);
    Stuff[x][a][b].average=ETLa0AVG;
    b++;
    ETLB0AVG=hETLambdaBar0->GetMean(1);
    Stuff[x][a][b].average=ETLB0AVG;
    b++;
    ETPROAVG=hETp->GetMean(1);
    Stuff[x][a][b].average=ETPROAVG;
    b++;
    ETNEUAVG=hETn->GetMean(1);
    Stuff[x][a][b].average=ETNEUAVG;
    b++;
    ETPRBAVG=hETpBar->GetMean(1);
    Stuff[x][a][b].average=ETPRBAVG;
    b++;
    ETNEBAVG=hETnBar->GetMean(1);
    Stuff[x][a][b].average=ETNEBAVG;
    b=0;


    ETALLS=hETAll->GetStdDev(1);
    Stuff[x][a][b].stddev=ETALLS;
    b++;
    ETPiPS=hETPiPlus->GetStdDev(1);
    Stuff[x][a][b].stddev=ETPiPS;
    b++;
    ETPiMS=hETPiMinus->GetStdDev(1);
    Stuff[x][a][b].stddev=ETPiMS;
    b++;
    ETPi0S=hETPi0->GetStdDev(1);
    Stuff[x][a][b].stddev=ETPi0S;
    b++;
    ETKaPS=hETKPlus->GetStdDev(1);
    Stuff[x][a][b].stddev=ETKaPS;
    b++;
    ETKaMS=hETKMinus->GetStdDev(1);
    Stuff[x][a][b].stddev=ETKaMS;
    b++;
    ETK0LS=hETKL->GetStdDev(1);
    Stuff[x][a][b].stddev=ETK0LS;
    b++;
    ETK0SS=hETKS->GetStdDev(1);
    Stuff[x][a][b].stddev=ETK0SS;
    b++;
    ETETAS=hETEta->GetStdDev(1);
    Stuff[x][a][b].stddev=ETETAS;
    b++;
    ETOMGS=hETomega->GetStdDev(1);
    Stuff[x][a][b].stddev=ETOMGS;
    b++;
    ETLa0S=hETLambda0->GetStdDev(1);
    Stuff[x][a][b].stddev=ETLa0S;
    b++;
    ETLB0S=hETLambdaBar0->GetStdDev(1);
    Stuff[x][a][b].stddev=ETLB0S;
    b++;
    ETPROS=hETp->GetStdDev(1);
    Stuff[x][a][b].stddev=ETPROS;
    b++;
    ETNEUS=hETn->GetStdDev(1);
    Stuff[x][a][b].stddev=ETNEUS;
    b++;
    ETPRBS=hETpBar->GetStdDev(1);
    Stuff[x][a][b].stddev=ETPRBS;
    b++;
    ETNEBS=hETnBar->GetStdDev(1);
    Stuff[x][a][b].stddev=ETNEBS;
    b=0;


    ETALLN=hETAll->GetEntries();
    Stuff[x][a][b].entries=ETALLN;
    b++;
    ETPiPN=hETPiPlus->GetEntries();
    Stuff[x][a][b].entries=ETPiPN;
    b++;
    ETPiMN=hETPiMinus->GetEntries();
    Stuff[x][a][b].entries=ETPiMN;
    b++;
    ETPi0N=hETPi0->GetEntries();
    Stuff[x][a][b].entries=ETPi0N;
    b++;
    ETKaPN=hETKPlus->GetEntries();
    Stuff[x][a][b].entries=ETKaPN;
    b++;
    ETKaMN=hETKMinus->GetEntries();
    Stuff[x][a][b].entries=ETKaMN;
    b++;
    ETK0LN=hETKL->GetEntries();
    Stuff[x][a][b].entries=ETK0LN;
    b++;
    ETK0SN=hETKS->GetEntries();
    Stuff[x][a][b].entries=ETK0SN;
    b++;
    ETETAN=hETEta->GetEntries();
    Stuff[x][a][b].entries=ETETAN;
    b++;
    ETOMGN=hETomega->GetEntries();
    Stuff[x][a][b].entries=ETOMGN;
    b++;
    ETLa0N=hETLambda0->GetEntries();
    Stuff[x][a][b].entries=ETLa0N;
    b++;
    ETLB0N=hETLambdaBar0->GetEntries();
    Stuff[x][a][b].entries=ETLB0N;
    b++;
    ETPRON=hETp->GetEntries();
    Stuff[x][a][b].entries=ETPRON;
    b++;
    ETNEUN=hETn->GetEntries();
    Stuff[x][a][b].entries=ETNEUN;
    b++;
    ETPRBN=hETpBar->GetEntries();
    Stuff[x][a][b].entries=ETPRBN;
    b++;
    ETNEBN=hETnBar->GetEntries();
    Stuff[x][a][b].entries=ETNEBN;
    b=0;

    a++;
    PTALLAVG=hPTAll->GetMean(1);
    Stuff[x][a][b].average=PTALLAVG;
    b++;
    PTPiPAVG=hptPiPlus->GetMean(1);
    Stuff[x][a][b].average=PTPiPAVG;
    b++;
    PTPiMAVG=hptPiMinus->GetMean(1);
    Stuff[x][a][b].average=PTPiMAVG;
    b++;
    PTPi0AVG=hptPi0->GetMean(1);
    Stuff[x][a][b].average=PTPi0AVG;
    b++;
    PTKaPAVG=hptKPlus->GetMean(1);
    Stuff[x][a][b].average=PTKaPAVG;
    b++;
    PTKaMAVG=hptKMinus->GetMean(1);
    Stuff[x][a][b].average=PTKaMAVG;
    b++;
    PTK0LAVG=hptKL->GetMean(1);
    Stuff[x][a][b].average=PTK0LAVG;
    b++;
    PTK0SAVG=hptKS->GetMean(1);
    Stuff[x][a][b].average=PTK0SAVG;
    b++;
    PTETAAVG=hptEta->GetMean(1);
    Stuff[x][a][b].average=PTETAAVG;
    b++;
    PTOMGAVG=hptomega->GetMean(1);
    Stuff[x][a][b].average=PTOMGAVG;
    b++;
    PTLa0AVG=hptLambda0->GetMean(1);
    Stuff[x][a][b].average=PTLa0AVG;
    b++;
    PTLB0AVG=hptLambdaBar0->GetMean(1);
    Stuff[x][a][b].average=PTLB0AVG;
    b++;
    PTPROAVG=hptp->GetMean(1);
    Stuff[x][a][b].average=PTPROAVG;
    b++;
    PTNEUAVG=hptn->GetMean(1);
    Stuff[x][a][b].average=PTNEUAVG;
    b++;
    PTPRBAVG=hptpBar->GetMean(1);
    Stuff[x][a][b].average=PTPRBAVG;
    b++;
    PTNEBAVG=hptnBar->GetMean(1);
    Stuff[x][a][b].average=PTNEBAVG;
    b=0;

    PTALLS=hPTAll->GetStdDev(1);
    Stuff[x][a][b].stddev=PTALLS;
    b++;
    PTPiPS=hptPiPlus->GetStdDev(1);
    Stuff[x][a][b].stddev=PTPiPS;
    b++;
    PTPiMS=hptPiMinus->GetStdDev(1);
    Stuff[x][a][b].stddev=PTPiMS;
    b++;
    PTPi0S=hptPi0->GetStdDev(1);
    Stuff[x][a][b].stddev=PTPi0S;
    b++;
    PTKaPS=hptKPlus->GetStdDev(1);
    Stuff[x][a][b].stddev=PTKaPS;
    b++;
    PTKaMS=hptKMinus->GetStdDev(1);
    Stuff[x][a][b].stddev=PTKaMS;
    b++;
    PTK0LS=hptKL->GetStdDev(1);
    Stuff[x][a][b].stddev=PTK0LS;
    b++;
    PTK0SS=hptKS->GetStdDev(1);
    Stuff[x][a][b].stddev=PTK0SS;
    b++;
    PTETAS=hptEta->GetStdDev(1);
    Stuff[x][a][b].stddev=PTETAS;
    b++;
    PTOMGS=hptomega->GetStdDev(1);
    Stuff[x][a][b].stddev=PTOMGS;
    b++;
    PTLa0S=hptLambda0->GetStdDev(1);
    Stuff[x][a][b].stddev=PTLa0S;
    b++;
    PTLB0S=hptLambdaBar0->GetStdDev(1);
    Stuff[x][a][b].stddev=PTLB0S;
    b++;
    PTPROS=hptp->GetStdDev(1);
    Stuff[x][a][b].stddev=PTPROS;
    b++;
    PTNEUS=hptn->GetStdDev(1);
    Stuff[x][a][b].stddev=PTNEUS;
    b++;
    PTPRBS=hptpBar->GetStdDev(1);
    Stuff[x][a][b].stddev=PTPRBS;
    b++;
    PTNEBS=hptnBar->GetStdDev(1);
    Stuff[x][a][b].stddev=PTNEBS;
    b=0;

    PTALLN=hPTAll->GetEntries();
    Stuff[x][a][b].entries=PTALLN;
    b++;
    PTPiPN=hptPiPlus->GetEntries();
    Stuff[x][a][b].entries=PTPiPN;
    b++;
    PTPiMN=hptPiMinus->GetEntries();
    Stuff[x][a][b].entries=PTPiMN;
    b++;
    PTPi0N=hptPi0->GetEntries();
    Stuff[x][a][b].entries=PTPi0N;
    b++;
    PTKaPN=hptKPlus->GetEntries();
    Stuff[x][a][b].entries=PTKaPN;
    b++;
    PTKaMN=hptKMinus->GetEntries();
    Stuff[x][a][b].entries=PTKaMN;
    b++;
    PTK0LN=hptKL->GetEntries();
    Stuff[x][a][b].entries=PTK0LN;
    b++;
    PTK0SN=hptKS->GetEntries();
    Stuff[x][a][b].entries=PTK0SN;
    b++;
    PTETAN=hptEta->GetEntries();
    Stuff[x][a][b].entries=PTETAN;
    b++;
    PTOMGN=hptomega->GetEntries();
    Stuff[x][a][b].entries=PTOMGN;
    b++;
    PTLa0N=hptLambda0->GetEntries();
    Stuff[x][a][b].entries=PTLa0N;
    b++;
    PTLB0N=hptLambdaBar0->GetEntries();
    Stuff[x][a][b].entries=PTLB0N;
    b++;
    PTPRON=hptp->GetEntries();
    Stuff[x][a][b].entries=PTPRON;
    b++;
    PTNEUN=hptn->GetEntries();
    Stuff[x][a][b].entries=PTNEUN;
    b++;
    PTPRBN=hptpBar->GetEntries();
    Stuff[x][a][b].entries=PTPRBN;
    b++;
    PTNEBN=hptnBar->GetEntries();
    Stuff[x][a][b].entries=PTNEBN;
    b=0;


    cout<<"Energy: "<<SNNACTUAL[x]<<" GeV complete"<<endl;
    //cout<<"ET pi+: "<<ETPiPAVG<<"\t"<<ETPiPS<<"\t"<<ETPiPN<<endl;
    //cout<<"PT pi+: "<<PTPiPAVG<<"\t"<<PTPiPS<<"\t"<<PTPiPN<<endl<<endl;
}
//cout<<"14GeV ET ALL"<<Stuff[14][0][0].average<<"\t"<<Stuff[14][0][0].stddev<<"\t"<<Stuff[14][0][0].entries<<endl;
//cout<<"14GeV PT ALL"<<Stuff[14][1][0].average<<"\t"<<Stuff[14][1][0].stddev<<"\t"<<Stuff[14][1][0].entries<<endl;
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
/*
//Total
for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  y[i]=Stuff[i][0][0].average;
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
  y[i]=Stuff[i][0][1].average;
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
  y[i]=Stuff[i][1][1].average;
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
  y[i]=Stuff[i][0][2].average;
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
  y[i]=Stuff[i][1][2].average;
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
  y[i]=Stuff[i][0][3].average;
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
  y[i]=Stuff[i][1][3].average;
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
  y[i]=Stuff[i][0][4].average;
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
  y[i]=Stuff[i][1][4].average;
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
  y[i]=Stuff[i][0][5].average;
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
  y[i]=Stuff[i][1][5].average;
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
  y[i]=Stuff[i][0][6].average;
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
  y[i]=Stuff[i][1][6].average;
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
  y[i]=Stuff[i][0][7].average;
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
  y[i]=Stuff[i][1][7].average;
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
  y[i]=Stuff[i][0][8].average;
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
  y[i]=Stuff[i][1][8].average;
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
  y[i]=Stuff[i][0][9].average;
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
  y[i]=Stuff[i][1][9].average;
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
  y[i]=Stuff[i][0][10].average;
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
  y[i]=Stuff[i][1][10].average;
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
  y[i]=Stuff[i][0][11].average;
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
  y[i]=Stuff[i][1][11].average;
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
  y[i]=Stuff[i][0][12].average;
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
  y[i]=Stuff[i][1][12].average;
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
  y[i]=Stuff[i][0][13].average;
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
  y[i]=Stuff[i][1][13].average;
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
  y[i]=Stuff[i][0][14].average;
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
  y[i]=Stuff[i][1][14].average;
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
  y[i]=Stuff[i][0][15].average;
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
  y[i]=Stuff[i][1][15].average;
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
*/
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
  blarg=((Stuff[i][0][1].average+Stuff[i][0][2].average)/2)/(Stuff[i][0][3].average);
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
  gc1->GetYaxis()->SetRangeUser(.5,2.5);
  gc1->Draw("AP");

for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i]*1.033;
  blarg=(Stuff[i][0][1].average)/(Stuff[i][0][3].average);
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
  blarg=(Stuff[i][0][2].average)/(Stuff[i][0][3].average);
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

for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i]*1.1;
  blarg=(Stuff[i][0][1].average+Stuff[i][0][2].average)/(Stuff[i][0][3].average);
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
  legend2->Draw();

Int_t n1=3;
Int_t x1[3]={0,100,15000};
Int_t y1[3]={1,1,1};
TGraph* gl=new TGraph(n1,x1,y1);
  gl->SetLineColor(kBlack);
  gl->Draw("L");
c16->SetLogx();
c16->Update();
char* file16= Form("Pi0Check_Run%i.png",IDNUM);
c16->SaveAs(file16);

//kaon check
TCanvas *c17 =new TCanvas("c17","kaon checks",0,645,1280,745);
for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  y[i]=(Stuff[i][0][4].average)/(Stuff[i][0][5].average);
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
  gc2->GetYaxis()->SetRangeUser(.5,1.5);
  gc2->Draw("AP");

for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i]*1.033;
  y[i]=(Stuff[i][0][4].average)/(Stuff[i][0][6].average);
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
  y[i]=(Stuff[i][0][4].average)/(Stuff[i][0][7].average);
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
  y[i]=(Stuff[i][0][4].average+Stuff[i][0][5].average)/(Stuff[i][0][6].average+Stuff[i][0][7].average);
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
  y[i]=(Stuff[i][0][6].average)/(Stuff[i][0][7].average);
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
char* file17= Form("KaonChecks_Run%i.png",IDNUM);
c17->SaveAs(file17);

//pi0 test
TCanvas *c18 =new TCanvas("c18","p-n checks",0,645,1280,745);
for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  blarg=((Stuff[i][0][12].average)/(Stuff[i][0][13].average));
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
  gc8->GetYaxis()->SetRangeUser(.5,1.5);
  gc8->Draw("AP");

for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i]*1.033;
  blarg=((Stuff[i][0][14].average)/(Stuff[i][0][15].average));
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
  blarg=((Stuff[i][0][12].average)/(Stuff[i][0][13].average));
  y[i]=blarg;
  ex[i]=0;
  blarg=((Stuff[i][0][12].stddev)/(Stuff[i][0][13].stddev));
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
char* file18= Form("pncheck_Run%i.png",IDNUM);
c18->SaveAs(file18);

//vsALL
TCanvas *c19 =new TCanvas("c19","Particle Et vs ETAll",0,645,1280,745);
for (Int_t i=0;i<n;i++){
  x[i]=SNNACTUAL[i];
  blarg=1/((Stuff[i][0][1].average+Stuff[i][0][2].average+Stuff[i][0][3].average)/(Stuff[i][0][0].average));
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
  blarg=1/((Stuff[i][0][4].average+Stuff[i][0][5].average+Stuff[i][0][6].average+Stuff[i][0][7].average)/(Stuff[i][0][0].average));
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
  blarg=1/((Stuff[i][0][12].average+Stuff[i][0][13].average+Stuff[i][0][14].average+Stuff[i][0][15].average)/(Stuff[i][0][0].average));
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

cout<<"program complete"<<endl<<endl;
return 0;
}

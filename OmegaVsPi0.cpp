/*******************************************************************************
OmegaVsPi0.cpp
this program takes data from fig4 of 'Common Suppression Pattern of [eta] and
[pi0] Mesons at High Transverse Momentum in Au+Au Collisions at [sqrt(sNN)]
=200Gev' and interpolates a fitting curve for each data set (with error) for the
purpose of aiding the calculation of eta in biswas's program
*******************************************************************************/

#include <cstdio>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "TGraph.h"
#include "TSpline.h"
#include "TString.h"
#include "TMath.h"
#include "TH1D.h"
#include "TKey.h"
#include "Riostream.h"
#include "TString.h"
using namespace std;

struct BinData {
  double pTl;
  double pTh;
  double pTSpec;
  double ErrStat;
  double ErrSys;
  float SNN;
  int cent;
};

struct Data{
  float evp;
  float pt;
  float errR;
  float errPt;
};

struct DataBES{
  float evp;
  float pt;
  float errR;
  float errPt;
  float SNN;
};

double ErrorBuilder(double, double);

void OmegaVsPi0(){
  int p=0; //counter for iteration number
  int k=0; //errorcounter
  int r=0;
  int s=0;
  int rsize;
  float SNN;
  double data;
  double temp1;
  double temp2;
  double temp3;
  string type;
  string garbage;
  Data D;
  vector <Data> epluseminus(0); //Au+Au for central collisions
  vector <Data> pi0pippim(0); //Au+Au for semi-cental Collisions
  vector <Data> pi0gam(0); //Au+Au for peripheral collisions
  vector <Data> pi0pippimPRC75(0); //d+Au minimum bias
  vector <Data> pi0gamPRC(0);  //p+p
  vector <Data> ISR(0);
  vector <Data> E706(0);
  vector <BinData> omega(0);
  vector <BinData> pi0(0);
  vector <DataBES> BESDataInclomega(0);
  vector <Data> BESDataInclomega7(0);
  vector <Data> BESDataInclomega11(0);
  vector <Data> BESDataInclomega19(0);
  vector <Data> BESDataInclomega27(0);
  vector <Data> BESDataInclomega39(0);

  ifstream in;
  ifstream in2;
  in.open("./OmegaVsPi0.dat");
  in2.open("BESData_sorted_PlusOmega.txt");

  while(in>>type){
    cout<<endl<<"==================================================="<<endl;
    cout<<"processing: "<<type<<endl;
    in>>garbage;
    in>>garbage;
    in>>garbage;
    in>>garbage;
    if (p==0){
      for (int j=0;j<9;j++){
        epluseminus.push_back(Data());
        in>>epluseminus[j].evp;
        in>>epluseminus[j].pt;
        in>>epluseminus[j].errR;
        in>>epluseminus[j].errPt;
        cout<<epluseminus[j].evp<<"\t"<<epluseminus[j].pt<<"\t"<<epluseminus[j].errR<<"\t"<<epluseminus[j].errPt<<endl;
        in.clear();
      }
      cout<<"Part1 input complete"<<endl;
    }
    if (p==1){
      for (int j=0;j<21;j++){
        pi0pippim.push_back(Data());
        in>>pi0pippim[j].evp;
        in>>pi0pippim[j].pt;
        in>>pi0pippim[j].errR;
        in>>pi0pippim[j].errPt;
        cout<<pi0pippim[j].evp<<"\t"<<pi0pippim[j].pt<<"\t"<<pi0pippim[j].errR<<"\t"<<pi0pippim[j].errPt<<endl;
      }
      cout<<"Part2 input complete"<<endl;
    }
    if (p==2){
      for (int j=0;j<8;j++){
        pi0gam.push_back(Data());
        in>>pi0gam[j].evp;
        in>>pi0gam[j].pt;
        in>>pi0gam[j].errR;
        in>>pi0gam[j].errPt;
        cout<<pi0gam[j].evp<<"\t"<<pi0gam[j].pt<<"\t"<<pi0gam[j].errR<<"\t"<<pi0gam[j].errPt<<endl;
      }
      cout<<"Part3 input complete"<<endl;
    }
    if (p==3){
      for (int j=0;j<12;j++){
        pi0pippimPRC75.push_back(Data());
        in>>pi0pippimPRC75[j].evp;
        in>>pi0pippimPRC75[j].pt;
        in>>pi0pippimPRC75[j].errR;
        in>>pi0pippimPRC75[j].errPt;
        cout<<pi0pippimPRC75[j].evp<<"\t"<<pi0pippimPRC75[j].pt<<"\t"<<pi0pippimPRC75[j].errR<<"\t"<<pi0pippimPRC75[j].errPt<<endl;
      }
      cout<<"Part4 input complete"<<endl;
    }
    if (p==4){
      for (int j=0;j<5;j++){
        pi0gamPRC.push_back(Data());
        in>>pi0gamPRC[j].evp;
        in>>pi0gamPRC[j].pt;
        in>>pi0gamPRC[j].errR;
        in>>pi0gamPRC[j].errPt;
        cout<<pi0gamPRC[j].evp<<"\t"<<pi0gamPRC[j].pt<<"\t"<<pi0gamPRC[j].errR<<"\t"<<pi0gamPRC[j].errPt<<endl;
      }
      cout<<"Part5 input complete"<<endl<<endl;
    }
    if (p==5){
      for (int j=0;j<3;j++){
        ISR.push_back(Data());
        in>>ISR[j].evp;
        in>>ISR[j].pt;
        in>>ISR[j].errR;
        in>>ISR[j].errPt;
        cout<<ISR[j].evp<<"\t"<<ISR[j].pt<<"\t"<<ISR[j].errR<<"\t"<<ISR[j].errPt<<endl;
      }
      cout<<"Part6 input complete"<<endl<<endl;
    }
    if (p==6){
      for (int j=0;j<9;j++){
        E706.push_back(Data());
        in>>E706[j].evp;
        in>>E706[j].pt;
        in>>E706[j].errR;
        in>>E706[j].errPt;
        cout<<E706[j].evp<<"\t"<<E706[j].pt<<"\t"<<E706[j].errR<<"\t"<<E706[j].errPt<<endl;
      }
      cout<<"Part7 input complete"<<endl<<endl;
    }
    if (p==7){
      cout<<"Fatal Error: Feed in overcount"<<endl;
      return;
    }
    p++;
  }
  cout<<"==================================================="<<endl<<"processing BESData_InclEta.txt"<<endl;

  for(int q=0; q<5;q++){
    in2>>garbage;
    while (garbage!="Au+Au"){
      in2>>garbage;
    }
    in2>>SNN;
    in2>>garbage;
    while (garbage!="pi0"){
      in2>>garbage;
    }
    for(int x=0; x<9;x++){
      in2>>garbage;
      while (in2>>data){
        pi0.push_back(BinData());
        pi0[s].pTl=data;
        in2>>pi0[s].pTh;
        in2>>pi0[s].pTSpec;
        in2>>pi0[s].ErrStat;
        in2>>pi0[s].ErrSys;
        pi0[s].SNN=SNN;
        pi0[s].cent=x;
        //cout<<pi0[s].pTl<<" "<<pi0[s].pTh<<" "<<pi0[s].pTSpec<<" "<<pi0[s].ErrStat<<" "<<pi0[s].ErrSys<<" "<<pi0[s].SNN<<endl;
        s++;
      }
      in2.clear();
      in2>>garbage>>garbage;
    }
    while (garbage!="omega"){
      in2>>garbage;
    }
    for(int x=0; x<9;x++){
      in2>>garbage;
      while (in2>>data){
        omega.push_back(BinData());
        omega[r].pTl=data;
        in2>>omega[r].pTh;
        in2>>omega[r].pTSpec;
        in2>>omega[r].ErrStat;
        in2>>omega[r].ErrSys;
        omega[r].SNN=SNN;
        omega[r].cent=x;
        //cout<<omega[r].pTl<<" "<<omega[r].pTh<<" "<<omega[r].pTSpec<<" "<<omega[r].ErrStat<<" "<<omega[r].ErrSys<<" "<<omega[r].SNN<<endl;
        r++;
      }
      in2.clear();
      in2>>garbage>>garbage;
    }
  }
  in.clear();
  in2.clear();
  in.close();
  in2.close();
  rsize=r;
  for(r=0;r<rsize;r++){
      BESDataInclomega.push_back(DataBES());
      BESDataInclomega[r].evp=omega[r].pTSpec/pi0[r].pTSpec;
      BESDataInclomega[r].pt=.5*(omega[r].pTl+omega[r].pTh);
      BESDataInclomega[r].errR=0;//TODO
      BESDataInclomega[r].errPt=0;//TODO
      BESDataInclomega[r].SNN=omega[r].SNN;

/*
      cout<<BESDataInclomega[r].evp<<" ";
      cout<<BESDataInclomega[r].pt<<" ";
      cout<<BESDataInclomega[r].errR<<" ";
      cout<<BESDataInclomega[r].errPt<<" "<<omega[r].SNN<<endl;
*/

  }
  for (r=0;r<rsize;r++){
    if (BESDataInclomega[r].SNN<=8){
        D.evp=BESDataInclomega[r].evp;
        D.pt=BESDataInclomega[r].pt;
        D.errR=BESDataInclomega[r].errR;
        D.errPt=BESDataInclomega[r].errPt;
        BESDataInclomega7.push_back(D);
    }
    if ((BESDataInclomega[r].SNN>=11)&&(BESDataInclomega[r].SNN<=12)){
        D.evp=BESDataInclomega[r].evp;
        D.pt=BESDataInclomega[r].pt;
        D.errR=BESDataInclomega[r].errR;
        D.errPt=BESDataInclomega[r].errPt;
        BESDataInclomega11.push_back(D);
    }
    if ((BESDataInclomega[r].SNN>=19)&&(BESDataInclomega[r].SNN<=20)){
        D.evp=BESDataInclomega[r].evp;
        D.pt=BESDataInclomega[r].pt;
        D.errR=BESDataInclomega[r].errR;
        D.errPt=BESDataInclomega[r].errPt;
        BESDataInclomega19.push_back(D);
    }
    if ((BESDataInclomega[r].SNN>=26)&&(BESDataInclomega[r].SNN<=28)){
        D.evp=BESDataInclomega[r].evp;
        D.pt=BESDataInclomega[r].pt;
        D.errR=BESDataInclomega[r].errR;
        D.errPt=BESDataInclomega[r].errPt;
        BESDataInclomega27.push_back(D);
    }
    if (BESDataInclomega[r].SNN>=38){
        D.evp=BESDataInclomega[r].evp;
        D.pt=BESDataInclomega[r].pt;
        D.errR=BESDataInclomega[r].errR;
        D.errPt=BESDataInclomega[r].errPt;
        BESDataInclomega39.push_back(D);
    }
  }


  cout<<endl<<"Inputs complete: Creating Graphs"<<endl;

  TCanvas *c1 = new TCanvas("c1","omega/pi0", 1215,730,1213,705);
  Int_t n1=9,n2=21,n3=8,n4=12,n5=5,n6=2,n7=2;
  Int_t n8=BESDataInclomega7.size();
  Int_t n9=BESDataInclomega11.size();
  Int_t n10=BESDataInclomega19.size();
  Int_t n11=BESDataInclomega27.size();
  Int_t n12=BESDataInclomega39.size();
  Int_t n13=3, n14=9;
  Double_t x1[n1], y1[n1], ex1[n1], ey1[n1];
  Double_t x2[n2], y2[n2], ex2[n2], ey2[n2];
  Double_t x3[n3], y3[n3], ex3[n3], ey3[n3];
  Double_t x4[n4], y4[n4], ex4[n4], ey4[n4];
  Double_t x5[n5], y5[n5], ex5[n5], ey5[n5];
  Double_t x6[n6], y6[n2], x7[n2], y7[n2];
  Double_t x8[n8], y8[n8], ex8[n8], ey8[n8];
  Double_t x9[n9], y9[n9], ex9[n9], ey9[n9];
  Double_t x10[n10], y10[n10], ex10[n10], ey10[n10];
  Double_t x11[n11], y11[n11], ex11[n11], ey11[n11];
  Double_t x12[n12], y12[n12], ex12[n12], ey12[n12];
  Double_t x13[n13], y13[n13], ex13[n13], ey13[n13];
  Double_t x14[n14], y14[n14], ex14[n14], ey14[n14];

  TLegend* legend=new TLegend(.582164,.589706,.895954,.894118);
    legend->SetTextFont(72);
    legend->SetTextSize(.03);
    legend->SetFillColor(0);

  x7[0]=0;
  x7[1]=15.5;
  y7[0]=0;
  y7[1]=1;
  TGraph* gr7=new TGraph(n6,x7,y7);
    gr7->SetTitle("Omega/Pi0 Ratio at   snn=200; pT; Omega/Pi0");
    gr7->SetMarkerStyle(27);
    gr7->SetMarkerSize(0);
    gr7->SetMarkerColor(kBlack);
    gr7->GetXaxis()->SetRangeUser(0,14);
    gr7->GetYaxis()->SetRangeUser(0,2.5);
    gr7->Draw("Ap");

  for (Int_t l=0;l<n5;l++){;
    x5[l]=pi0gamPRC[l].pt;
    y5[l]=pi0gamPRC[l].evp;
    ex5[l]=pi0gamPRC[l].errPt;
    ey5[l]=pi0gamPRC[l].errR;
  }
  TGraphErrors* gr5= new TGraphErrors(n5,x5,y5,ex5,ey5);
    gr5->SetMarkerStyle(25);
    gr5->SetMarkerSize(1.5);
    gr5->SetMarkerColor(kBlack);
    //gr5->SetLineColor(kGray);
    gr5->Draw("P");

  for (Int_t l=0;l<n1;l++){;
    x1[l]=epluseminus[l].pt;
    y1[l]=epluseminus[l].evp;
    ex1[l]=epluseminus[l].errPt;
    ey1[l]=epluseminus[l].errR;
  }
  TGraphErrors* gr1= new TGraphErrors(n1,x1,y1,ex1,ey1);
    gr1->SetMarkerStyle(22);
    gr1->SetMarkerSize(1.5);
    gr1->SetMarkerColor(kBlack);
    //gr1->SetLineColor(kBlack);
    legend->AddEntry(gr1,"w->e+e-","p");
    gr1->Draw("P");


  for (Int_t l=0;l<n2;l++){;
    x2[l]=pi0pippim[l].pt;
    y2[l]=pi0pippim[l].evp;
    ex2[l]=pi0pippim[l].errPt;
    ey2[l]=pi0pippim[l].errR;
  }
  TGraphErrors* gr2= new TGraphErrors(n2,x2,y2,ex2,ey2);
    gr2->SetMarkerStyle(20);
    gr2->SetMarkerSize(1.5);
    gr2->SetMarkerColor(kBlue);
    //gr2->SetLineColor(kBlue);
    legend->AddEntry(gr2,"w->pi0pi+pi-","p");
    gr2->Draw("P");

  for (Int_t l=0;l<n3;l++){;
    x3[l]=pi0gam[l].pt;
    y3[l]=pi0gam[l].evp;
    ex3[l]=pi0gam[l].errPt;
    ey3[l]=pi0gam[l].errR;
  }
  TGraphErrors* gr3= new TGraphErrors(n3,x3,y3,ex3,ey3);
    gr3->SetMarkerStyle(21);
    gr3->SetMarkerSize(1.5);
    gr3->SetMarkerColor(kRed);
    //gr3->SetLineColor(kRed);
    legend->AddEntry(gr3,"w->pi0gamma","p");
    gr3->Draw("P");

  for (Int_t l=0;l<n4;l++){;
    x4[l]=pi0pippimPRC75[l].pt;
    y4[l]=pi0pippimPRC75[l].evp;
    ex4[l]=pi0pippimPRC75[l].errPt;
    ey4[l]=pi0pippimPRC75[l].errR;
  }
  TGraphErrors* gr4= new TGraphErrors(n4,x4,y4,ex4,ey4);
    gr4->SetMarkerStyle(4);
    gr4->SetMarkerSize(1.5);
    gr4->SetMarkerColor(kBlack);
    //gr4->SetLineColor(kBlack);
    legend->AddEntry(gr4,"w->pi0pi+pi-PRC75","p");
    legend->AddEntry(gr5,"w->pi0gammaPRC75","p");
    gr4->Draw("P");

  for (Int_t l=0;l<n13;l++){;
    x13[l]=ISR[l].pt;
    y13[l]=ISR[l].evp;
    ex13[l]=ISR[l].errPt;
    ey13[l]=ISR[l].errR;
  }
  TGraphErrors* gr13= new TGraphErrors(n13,x13,y13,ex13,ey13);
    gr13->SetMarkerStyle(31);
    gr13->SetMarkerSize(1.5);
    gr13->SetMarkerColor(kViolet);
    //gr4->SetLineColor(kBlack);
    legend->AddEntry(gr13,"ISR (p+p at SNN=62GeV)","p");
    gr13->Draw("P");

  for (Int_t l=0;l<n14;l++){;
    x14[l]=E706[l].pt;
    y14[l]=E706[l].evp;
    ex14[l]=E706[l].errPt;
    ey14[l]=E706[l].errR;
  }
  TGraphErrors* gr14= new TGraphErrors(n14,x14,y14,ex14,ey14);
    gr14->SetMarkerStyle(23);
    gr14->SetMarkerSize(1.5);
    gr14->SetMarkerColor(kViolet);
    //gr4->SetLineColor(kBlack);
    legend->AddEntry(gr14,"E706 (pi- + Be at SNN=31GeV)","p");
    gr14->Draw("P");

  x6[0]=0;
  x6[1]=18;
  y6[0]=1;
  y6[1]=1;
  TGraph* gr6=new TGraph(n6,x6,y6);
    gr6->SetMarkerStyle(27);
    gr6->SetMarkerSize(1.5);
    gr6->SetMarkerColor(kBlack);
    gr6->Draw("L");
    legend->Draw();

  for (Int_t l=0;l<n8;l++){;
    x8[l]=BESDataInclomega7[l].pt;
    y8[l]=BESDataInclomega7[l].evp;
    temp2=x8[l];
    temp3=y8[l];
    temp1=ErrorBuilder(temp2,temp3);
    BESDataInclomega7[l].errR=temp1;
    ex8[l]=BESDataInclomega7[l].errPt;
    ey8[l]=BESDataInclomega7[l].errR;
  }
  TGraphErrors* gr8= new TGraphErrors(n8,x8,y8,ex8,ey8);
    gr8->SetMarkerStyle(47);
    gr8->SetMarkerSize(1);
    //gr8->SetLineColor(kRed);
    gr8->SetMarkerColor(kGray+3);
    legend->AddEntry(gr8,"Au+Au BES 7.7","p");
    gr8->Draw("P");

  for (Int_t l=0;l<n9;l++){;
    x9[l]=BESDataInclomega11[l].pt;
    y9[l]=BESDataInclomega11[l].evp;
    temp2=x9[l];
    temp3=y9[l];
    temp1=ErrorBuilder(temp2,temp3);
    BESDataInclomega11[l].errR=temp1;
    ex9[l]=BESDataInclomega11[l].errPt;
    ey9[l]=BESDataInclomega11[l].errR;
  }
  TGraphErrors* gr9= new TGraphErrors(n9,x9,y9,ex9,ey9);
    gr9->SetMarkerStyle(33);
    gr9->SetMarkerSize(1);
    //gr8->SetLineColor(kRed);
    gr9->SetMarkerColor(kRed+2);
    legend->AddEntry(gr9,"Au+Au BES 11.5","p");
    gr9->Draw("P");

  for (Int_t l=0;l<n10;l++){;
    x10[l]=BESDataInclomega19[l].pt;
    y10[l]=BESDataInclomega19[l].evp;
    temp2=x10[l];
    temp3=y10[l];
    temp1=ErrorBuilder(temp2,temp3);
    BESDataInclomega19[l].errR=temp1;
    ex10[l]=BESDataInclomega19[l].errPt;
    ey10[l]=BESDataInclomega19[l].errR;
  }
  TGraphErrors* gr10= new TGraphErrors(n10,x10,y10,ex10,ey10);
    gr10->SetMarkerStyle(34);
    gr10->SetMarkerSize(1);
    //gr8->SetLineColor(kRed);
    gr10->SetMarkerColor(kBlue+3);
    legend->AddEntry(gr10,"Au+Au BES 19.6","p");
    gr10->Draw("P");

  for (Int_t l=0;l<n11;l++){;
    x11[l]=BESDataInclomega27[l].pt;
    y11[l]=BESDataInclomega27[l].evp;
    temp2=x11[l];
    temp3=y11[l];
    temp1=ErrorBuilder(temp2,temp3);
    BESDataInclomega27[l].errR=temp1;
    ex11[l]=BESDataInclomega27[l].errPt;
    ey11[l]=BESDataInclomega27[l].errR;
  }
  TGraphErrors* gr11= new TGraphErrors(n11,x11,y11,ex11,ey11);
    gr11->SetMarkerStyle(29);
    gr11->SetMarkerSize(1);
    //gr8->SetLineColor(kRed);
    gr11->SetMarkerColor(kOrange+9);
    legend->AddEntry(gr11,"Au+Au BES 27","p");
    gr11->Draw("P");

  for (Int_t l=0;l<n12;l++){;
    x12[l]=BESDataInclomega39[l].pt;
    y12[l]=BESDataInclomega39[l].evp;
    temp2=x12[l];
    temp3=y12[l];
    temp1=ErrorBuilder(temp2,temp3);
    BESDataInclomega39[l].errR=temp1;
    ex12[l]=BESDataInclomega39[l].errPt;
    ey12[l]=BESDataInclomega39[l].errR;
  }
  TGraphErrors* gr12= new TGraphErrors(n12,x12,y12,ex12,ey12);
    gr12->SetMarkerStyle(41);
    gr12->SetMarkerSize(1);
    //gr12->SetLineWidth(5);
    //gr12->SetLineColor(kGray);
    //gr8->SetLineColor(kRed);
    gr12->SetMarkerColor(kViolet-7);
    legend->AddEntry(gr12,"Au+Au BES 39","p");
    gr12->Draw("P");

  TF1* expo2 = new TF1("expo2","exp([0]-[1]/x)",-100000000000,100000000000);
    expo2->SetParName(0,"p0");
    expo2->SetParName(1,"p1");
    expo2->SetParameters(0.2,0.9);
    expo2->SetParLimits(0,0,100000);
    expo2->SetParLimits(0,0,10000);

  //c1->SetLogy();
  //c1->SetLogx();
  c1->Update();
  gStyle->SetOptFit(0);
  cout<<endl<<"interpolating..."<<endl;


  cout<<"w->e+e- fit is: "<<endl;
  gr1->Fit("expo2");
  gr1->GetFunction("expo2")->SetLineColor(kBlack);
  gr1->GetFunction("expo2")->SetLineWidth(1);
  cout<<"==================================================================="<<endl;
  cout<<"w->pi0pi+pi- fit is: "<<endl;
  gr2->Fit("expo2");
  gr2->GetFunction("expo2")->SetLineColor(kBlue);
  gr2->GetFunction("expo2")->SetLineWidth(1);
  cout<<"==================================================================="<<endl;
  cout<<"w->pi0gamma fit is: "<<endl;
  gr3->Fit("expo2");
  gr3->GetFunction("expo2")->SetLineColor(kRed);
  gr3->GetFunction("expo2")->SetLineWidth(1);
  cout<<"==================================================================="<<endl;
  cout<<"w->pi0pi+pi-PRC75 fit is: "<<endl;
  gr4->Fit("expo2");
  gr4->GetFunction("expo2")->SetLineColor(kBlack);
  gr4->GetFunction("expo2")->SetLineWidth(.25);
  cout<<"==================================================================="<<endl;
  cout<<"w->pi0gammaPRC75 fit is: "<<endl;
  gr5->Fit("expo2");
  gr5->GetFunction("expo2")->SetLineColor(kBlack);
  gr5->GetFunction("expo2")->SetLineWidth(.5);
  cout<<"==================================================================="<<endl;
  cout<<"ISR(pp) fit is: "<<endl;
  gr13->Fit("expo2");
  gr13->GetFunction("expo2")->SetLineColor(kViolet);
  gr13->GetFunction("expo2")->SetLineWidth(1);
  cout<<"==================================================================="<<endl;
  cout<<"E706 fit is: "<<endl;
  gr14->Fit("expo2");
  gr14->GetFunction("expo2")->SetLineColor(kMagenta);
  gr14->GetFunction("expo2")->SetLineWidth(1);
  cout<<"==================================================================="<<endl;


  cout<<"Au+Au BES 7.7 fit is: "<<endl;
  gr8->Fit("expo2");
  gr8->GetFunction("expo2")->SetLineColor(kBlack);
  gr8->GetFunction("expo2")->SetLineWidth(1);
  cout<<"==================================================================="<<endl;
  cout<<"Au+Au BES 11.5 fit is: "<<endl;
  gr9->Fit("expo2");
  gr9->GetFunction("expo2")->SetLineColor(kRed);
  gr9->GetFunction("expo2")->SetLineWidth(1);
  cout<<"==================================================================="<<endl;
  cout<<"Au+Au BES 19.6 fit is: "<<endl;
  gr10->Fit("expo2");
  gr10->GetFunction("expo2")->SetLineColor(kBlue);
  gr10->GetFunction("expo2")->SetLineWidth(1);
  cout<<"==================================================================="<<endl;
  cout<<"Au+Au BES 27 fit is: "<<endl;
  gr11->Fit("expo2");
  gr11->GetFunction("expo2")->SetLineColor(kGreen);
  gr11->GetFunction("expo2")->SetLineWidth(1);
  cout<<"==================================================================="<<endl;
  cout<<"Au+Au BES 39 fit is: "<<endl;
  gr12->Fit("expo2");
  gr12->GetFunction("expo2")->SetLineColor(kViolet);
  gr12->GetFunction("expo2")->SetLineWidth(1);
  cout<<"==================================================================="<<endl;

  c1->SaveAs("OmegaPi0_vs_pt.png");

  TCanvas *c2 = new TCanvas("c2","results",200,10,700,500);

  TLegend* legend2=new TLegend(.13,.6,.35,.9);
    legend2->SetTextFont(72);
    legend2->SetTextSize(.03);
    legend2->SetFillColor(0);

  x7[0]=0;
  x7[1]=2.1;
  y7[0]=0;
  y7[1]=0.59;
  TGraph* gr17=new TGraph(n6,x7,y7);
    gr17->SetTitle("omega/Pi0 Ratio at   snn=200; pT; omega/Pi0");
    gr17->GetXaxis()->SetRangeUser(0,2.2);
    gr17->GetYaxis()->SetRangeUser(0,1);
    gr17->SetMarkerStyle(27);
    gr17->SetMarkerSize(0);
    gr17->SetMarkerColor(kBlack);
    gr17->Draw("Ap");

  gr8->SetMarkerStyle(47);
  gr8->SetMarkerSize(1);
  gr8->SetMarkerColor(kGray+3);
  gr9->SetMarkerStyle(33);
  gr9->SetMarkerSize(1);
  gr9->SetMarkerColor(kRed+2);
  gr10->SetMarkerStyle(34);
  gr10->SetMarkerSize(1);
  gr10->SetMarkerColor(kBlue+3);
  gr11->SetMarkerStyle(29);
  gr11->SetMarkerSize(1);
  gr11->SetMarkerColor(kOrange+9);
  gr12->SetMarkerStyle(41);
  gr12->SetMarkerSize(1);
  gr12->SetMarkerColor(kViolet-7);

  gr8->Draw("P");
  gr9->Draw("P");
  gr10->Draw("P");
  gr11->Draw("P");
  gr12->Draw("P");
  legend2->AddEntry(gr8,"Au+Au BES 7.7","p");
  legend2->AddEntry(gr9,"Au+Au BES 11.5","p");
  legend2->AddEntry(gr10,"Au+Au BES 19.6","p");
  legend2->AddEntry(gr11,"Au+Au BES 27","p");
  legend2->AddEntry(gr12,"Au+Au BES 39","p");
  legend2->Draw();
  c2->Update();
  c2->SaveAs("omegaPi0_vs_pt_r.png");

  cout<<"complete";
  return;
}

double ErrorBuilder(double ptA, double ratio){
  vector<double> p0{0.0000000000205391,0.0000000000111022,0.0000000000355271,0,0.00000000522138,0,0};
  vector<double> p1{1.08474, .88624, .918029, .424807, .813936, .813922, .813922};
  vector<double> e0{2.31657e-01, 9.60156e-03, 2.72294e-02, 6.21560e-02, 1.32981e-01, 1.32966e-01, 1.32966e-01};
  vector<double> e1{1.43492e-01, 9.41344e-02, 1.53131e-01, 2.52255e-01, 3.49523e-01, 3.49520e-01, 3.49520e-01};
  double test=0;
  double x0;
  double x1;
  double comp=0;
  double error=0;
  int index=0;

  for (int i=0;i<5;i++){
    x0=abs(p0[i])+e0[i];
    x1=abs(p1[i])+e1[i];
    test=exp((-x0)-x1/ptA);
    comp=abs(test-ratio);
    if (comp>=error){
      error=comp;
    }
  }

  return error;
}

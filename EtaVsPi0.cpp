/*******************************************************************************
EtaVsPi0.cpp
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

void EtaVsPi0(){
  int p=0; //counter for iteration number
  int k=0; //errorcounter
  int r=0;
  int s=0;
  int rsize;
  float SNN;
  double data;
  string type;
  string garbage;
  Data D;
  vector <Data> auaucent(0); //Au+Au for central collisions
  vector <Data> auauScent(0); //Au+Au for semi-cental Collisions
  vector <Data> auauPeri(0); //Au+Au for peripheral collisions
  vector <Data> dau(0); //d+Au minimum bias
  vector <Data> pp(0);  //p+p
  vector <BinData> Eta(0);
  vector <BinData> pi0(0);
  vector <DataBES> BESDataInclEta(0);
  vector <Data> BESDataInclEta7(0);
  vector <Data> BESDataInclEta11(0);
  vector <Data> BESDataInclEta19(0);
  vector <Data> BESDataInclEta27(0);
  vector <Data> BESDataInclEta39(0);

  ifstream in;
  ifstream in2;
  in.open("./EtaVsPi0.dat");
  in2.open("BESData_sorted_PlusEta.txt");

  while(in>>type){
    cout<<endl<<"==================================================="<<endl;
    cout<<"processing: "<<type<<endl;
    in>>garbage;
    in>>garbage;
    in>>garbage;
    in>>garbage;
    if (p==0){
      for (int j=0;j<10;j++){
        auaucent.push_back(Data());
        in>>auaucent[j].evp;
        in>>auaucent[j].pt;
        in>>auaucent[j].errR;
        in>>auaucent[j].errPt;
        cout<<auaucent[j].evp<<"\t"<<auaucent[j].pt<<"\t"<<auaucent[j].errR<<"\t"<<auaucent[j].errPt<<endl;
        in.clear();
      }
      cout<<"Part1 input complete"<<endl;
    }
    if (p==1){
      for (int j=0;j<10;j++){
        auauScent.push_back(Data());
        in>>auauScent[j].evp;
        in>>auauScent[j].pt;
        in>>auauScent[j].errR;
        in>>auauScent[j].errPt;
        cout<<auauScent[j].evp<<"\t"<<auauScent[j].pt<<"\t"<<auauScent[j].errR<<"\t"<<auauScent[j].errPt<<endl;
      }
      cout<<"Part2 input complete"<<endl;
    }
    if (p==2){
      for (int j=0;j<6;j++){
        auauPeri.push_back(Data());
        in>>auauPeri[j].evp;
        in>>auauPeri[j].pt;
        in>>auauPeri[j].errR;
        in>>auauPeri[j].errPt;
        cout<<auauPeri[j].evp<<"\t"<<auauPeri[j].pt<<"\t"<<auauPeri[j].errR<<"\t"<<auauPeri[j].errPt<<endl;
      }
      cout<<"Part3 input complete"<<endl;
    }
    if (p==3){
      for (int j=0;j<13;j++){
        dau.push_back(Data());
        in>>dau[j].evp;
        in>>dau[j].pt;
        in>>dau[j].errR;
        in>>dau[j].errPt;
        cout<<dau[j].evp<<"\t"<<dau[j].pt<<"\t"<<dau[j].errR<<"\t"<<dau[j].errPt<<endl;
      }
      cout<<"Part4 input complete"<<endl;
    }
    if (p==4){
      for (int j=0;j<12;j++){
        pp.push_back(Data());
        in>>pp[j].evp;
        in>>pp[j].pt;
        in>>pp[j].errR;
        in>>pp[j].errPt;
        cout<<pp[j].evp<<"\t"<<pp[j].pt<<"\t"<<pp[j].errR<<"\t"<<pp[j].errPt<<endl;
      }
      cout<<"Part5 input complete"<<endl<<endl;
    }
    if (p==5){
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
    while (garbage!="eta"){
      in2>>garbage;
    }
    for(int x=0; x<9;x++){
      in2>>garbage;
      while (in2>>data){
        Eta.push_back(BinData());
        Eta[r].pTl=data;
        in2>>Eta[r].pTh;
        in2>>Eta[r].pTSpec;
        in2>>Eta[r].ErrStat;
        in2>>Eta[r].ErrSys;
        Eta[r].SNN=SNN;
        Eta[r].cent=x;
        //cout<<Eta[r].pTl<<" "<<Eta[r].pTh<<" "<<Eta[r].pTSpec<<" "<<Eta[r].ErrStat<<" "<<Eta[r].ErrSys<<" "<<Eta[r].SNN<<endl;
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
      BESDataInclEta.push_back(DataBES());
      BESDataInclEta[r].evp=Eta[r].pTSpec/pi0[r].pTSpec;//BUG not matched to pt bins!!
      BESDataInclEta[r].pt=.5*(Eta[r].pTl+Eta[r].pTh);
      BESDataInclEta[r].errR=0;//TODO
      BESDataInclEta[r].errPt=0;//TODO
      BESDataInclEta[r].SNN=Eta[r].SNN;
/*
      cout<<BESDataInclEta[r].evp<<" ";
      cout<<BESDataInclEta[r].pt<<" ";
      cout<<BESDataInclEta[r].errR<<" ";
      cout<<BESDataInclEta[r].errPt<<" "<<Eta[r].SNN<<endl;
*/
  }
  for (r=0;r<rsize;r++){
    if (BESDataInclEta[r].SNN<=7.7){
        D.evp=BESDataInclEta[r].evp;
        D.pt=BESDataInclEta[r].pt;
        D.errR=BESDataInclEta[r].errR;
        D.errPt=BESDataInclEta[r].errPt;
        BESDataInclEta7.push_back(D);
    }
    if (BESDataInclEta[r].SNN==11.5){
        D.evp=BESDataInclEta[r].evp;
        D.pt=BESDataInclEta[r].pt;
        D.errR=BESDataInclEta[r].errR;
        D.errPt=BESDataInclEta[r].errPt;
        BESDataInclEta11.push_back(D);
    }
    if (BESDataInclEta[r].SNN==19.6){
        D.evp=BESDataInclEta[r].evp;
        D.pt=BESDataInclEta[r].pt;
        D.errR=BESDataInclEta[r].errR;
        D.errPt=BESDataInclEta[r].errPt;
        BESDataInclEta19.push_back(D);
    }
    if (BESDataInclEta[r].SNN==27){
        D.evp=BESDataInclEta[r].evp;
        D.pt=BESDataInclEta[r].pt;
        D.errR=BESDataInclEta[r].errR;
        D.errPt=BESDataInclEta[r].errPt;
        BESDataInclEta27.push_back(D);
    }
    if (BESDataInclEta[r].SNN==39){
        D.evp=BESDataInclEta[r].evp;
        D.pt=BESDataInclEta[r].pt;
        D.errR=BESDataInclEta[r].errR;
        D.errPt=BESDataInclEta[r].errPt;
        BESDataInclEta39.push_back(D);
    }
  }


  cout<<endl<<"Inputs complete: Creating Graphs"<<endl;

  TCanvas *c1 = new TCanvas("c1","Au+Au Central", 200,10,700,500);
  Int_t n1=10,n2=10,n3=6,n4=13,n5=12,n6=2,n7=2;
  Int_t n8=BESDataInclEta7.size();
  Int_t n9=BESDataInclEta11.size();
  Int_t n10=BESDataInclEta19.size();
  Int_t n11=BESDataInclEta27.size();
  Int_t n12=BESDataInclEta39.size();
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

  TLegend* legend=new TLegend(.13,.6,.35,.9);
    legend->SetTextFont(72);
    legend->SetTextSize(.03);
    legend->SetFillColor(0);

  x7[0]=0;
  x7[1]=11.5;
  y7[0]=0;
  y7[1]=1.2;
  TGraph* gr7=new TGraph(n6,x7,y7);
    gr7->SetTitle("Eta/Pi0 Ratio at   snn=200; pT; Eta/Pi0");
    gr7->SetMarkerStyle(27);
    gr7->SetMarkerSize(0);
    gr7->SetMarkerColor(kBlack);
    gr7->Draw("Ap");

  for (Int_t l=0;l<n5;l++){;
    x5[l]=pp[l].pt;
    y5[l]=pp[l].evp;
    ex5[l]=pp[l].errPt;
    ey5[l]=pp[l].errR;
  }
  TGraphErrors* gr5= new TGraphErrors(n5,x5,y5,ex5,ey5);
    gr5->SetMarkerStyle(28);
    gr5->SetMarkerSize(1.5);
    gr5->SetMarkerColor(kBlack);
    //gr5->SetLineColor(kGray);
    gr5->Draw("P");

  for (Int_t l=0;l<n1;l++){;
    x1[l]=auaucent[l].pt;
    y1[l]=auaucent[l].evp;
    ex1[l]=auaucent[l].errPt;
    ey1[l]=auaucent[l].errR;
  }
  TGraphErrors* gr1= new TGraphErrors(n1,x1,y1,ex1,ey1);
    gr1->SetMarkerStyle(20);
    gr1->SetMarkerSize(1.5);
    gr1->SetMarkerColor(kRed);
    //gr1->SetLineColor(kRed);
    legend->AddEntry(gr1,"Au+Au Cent","p");
    gr1->Draw("P");


  for (Int_t l=0;l<n2;l++){;
    x2[l]=auauScent[l].pt;
    y2[l]=auauScent[l].evp;
    ex2[l]=auauScent[l].errPt;
    ey2[l]=auauScent[l].errR;
  }
  TGraphErrors* gr2= new TGraphErrors(n2,x2,y2,ex2,ey2);
    gr2->SetMarkerStyle(21);
    gr2->SetMarkerSize(1.5);
    gr2->SetMarkerColor(kOrange);
    //gr2->SetLineColor(kOrange);
    legend->AddEntry(gr2,"Au+Au Semi-Cent","p");
    gr2->Draw("P");

  for (Int_t l=0;l<n3;l++){;
    x3[l]=auauPeri[l].pt;
    y3[l]=auauPeri[l].evp;
    ex3[l]=auauPeri[l].errPt;
    ey3[l]=auauPeri[l].errR;
  }
  TGraphErrors* gr3= new TGraphErrors(n3,x3,y3,ex3,ey3);
    gr3->SetMarkerStyle(22);
    gr3->SetMarkerSize(1.5);
    gr3->SetMarkerColor(kBlue);
    //gr3->SetLineColor(kBlue);
    legend->AddEntry(gr3,"Au+Au Periph","p");
    gr3->Draw("P");

  for (Int_t l=0;l<n4;l++){;
    x4[l]=dau[l].pt;
    y4[l]=dau[l].evp;
    ex4[l]=dau[l].errPt;
    ey4[l]=dau[l].errR;
  }
  TGraphErrors* gr4= new TGraphErrors(n4,x4,y4,ex4,ey4);
    gr4->SetMarkerStyle(27);
    gr4->SetMarkerSize(1.5);
    gr4->SetMarkerColor(kViolet);
    //gr4->SetLineColor(kBlack);
    legend->AddEntry(gr4,"d+Au MinBias","p");
    legend->AddEntry(gr5,"p+p","p");
    gr4->Draw("P");

  x6[0]=0;
  x6[1]=13;
  y6[0]=1;
  y6[1]=1;
  TGraph* gr6=new TGraph(n6,x6,y6);
    gr6->SetMarkerStyle(27);
    gr6->SetMarkerSize(1.5);
    gr6->SetMarkerColor(kBlack);
    gr6->Draw("L");
    legend->Draw();

  for (Int_t l=0;l<n8;l++){;
    x8[l]=BESDataInclEta7[l].pt;
    y8[l]=BESDataInclEta7[l].evp;
    ex8[l]=BESDataInclEta7[l].errPt;
    ey8[l]=BESDataInclEta7[l].errR;
  }
  TGraphErrors* gr8= new TGraphErrors(n8,x8,y8,ex8,ey8);
    gr8->SetMarkerStyle(20);
    gr8->SetMarkerSize(1);
    //gr8->SetLineColor(kRed);
    gr8->SetMarkerColor(kBlack);
    legend->AddEntry(gr8,"Au+Au BES 7.7","p");
    gr8->Draw("P");

  for (Int_t l=0;l<n9;l++){;
    x9[l]=BESDataInclEta11[l].pt;
    y9[l]=BESDataInclEta11[l].evp;
    ex9[l]=BESDataInclEta11[l].errPt;
    ey9[l]=BESDataInclEta11[l].errR;
  }
  TGraphErrors* gr9= new TGraphErrors(n9,x9,y9,ex9,ey9);
    gr9->SetMarkerStyle(20);
    gr9->SetMarkerSize(1);
    //gr8->SetLineColor(kRed);
    gr9->SetMarkerColor(kRed);
    legend->AddEntry(gr9,"Au+Au BES 11.5","p");
    gr9->Draw("P");

  for (Int_t l=0;l<n10;l++){;
    x10[l]=BESDataInclEta19[l].pt;
    y10[l]=BESDataInclEta19[l].evp;
    ex10[l]=BESDataInclEta19[l].errPt;
    ey10[l]=BESDataInclEta19[l].errR;
  }
  TGraphErrors* gr10= new TGraphErrors(n10,x10,y10,ex10,ey10);
    gr10->SetMarkerStyle(20);
    gr10->SetMarkerSize(1);
    //gr8->SetLineColor(kRed);
    gr10->SetMarkerColor(kBlue);
    legend->AddEntry(gr10,"Au+Au BES 19.6","p");
    gr10->Draw("P");

  for (Int_t l=0;l<n11;l++){;
    x11[l]=BESDataInclEta27[l].pt;
    y11[l]=BESDataInclEta27[l].evp;
    ex11[l]=BESDataInclEta27[l].errPt;
    ey11[l]=BESDataInclEta27[l].errR;
  }
  TGraphErrors* gr11= new TGraphErrors(n11,x11,y11,ex11,ey11);
    gr11->SetMarkerStyle(20);
    gr11->SetMarkerSize(1);
    //gr8->SetLineColor(kRed);
    gr11->SetMarkerColor(kGray);
    legend->AddEntry(gr11,"Au+Au BES 27","p");
    gr11->Draw("P");

  for (Int_t l=0;l<n12;l++){;
    x12[l]=BESDataInclEta39[l].pt;
    y12[l]=BESDataInclEta39[l].evp;
    ex12[l]=BESDataInclEta39[l].errPt;
    ey12[l]=BESDataInclEta39[l].errR;
  }
  TGraphErrors* gr12= new TGraphErrors(n12,x12,y12,ex12,ey12);
    gr12->SetMarkerStyle(20);
    gr12->SetMarkerSize(1);
    //gr8->SetLineColor(kRed);
    gr12->SetMarkerColor(kViolet);
    legend->AddEntry(gr12,"Au+Au BES 39","p");
    gr12->Draw("P");

  TF1* expo2 = new TF1("expo2","exp([0]-[1]/x)",-100000000000,100000000000);
    expo2->SetParName(0,"p0");
    expo2->SetParName(1,"p1");
    expo2->SetParameters(0,12);

  //c1->SetLogy();
  //c1->SetLogx();
  c1->Update();
  gr5->GetXaxis()->SetRangeUser(.001,12);
  gr5->GetYaxis()->SetRangeUser(.001,1.2);
  gStyle->SetOptFit(0);
  cout<<endl<<"interpolating..."<<endl;


  cout<<"Au+Au Central fit is: "<<endl;
  gr1->Fit("expo2");
  gr1->GetFunction("expo2")->SetLineColor(kRed);
  gr1->GetFunction("expo2")->SetLineWidth(1);
  cout<<"==================================================================="<<endl;
  cout<<"Au+Au Semi-Central fit is: "<<endl;
  gr2->Fit("expo2");
  gr2->GetFunction("expo2")->SetLineColor(kOrange);
  gr2->GetFunction("expo2")->SetLineWidth(1);
  cout<<"==================================================================="<<endl;
  cout<<"Au+Au Peripheral fit is: "<<endl;
  gr3->Fit("expo2");
  gr3->GetFunction("expo2")->SetLineColor(kBlue);
  gr3->GetFunction("expo2")->SetLineWidth(1);
  cout<<"==================================================================="<<endl;
  cout<<"d+Au Min Bias fit is: "<<endl;
  gr4->Fit("expo2");
  gr4->GetFunction("expo2")->SetLineColor(kViolet);
  gr4->GetFunction("expo2")->SetLineWidth(1);
  cout<<"==================================================================="<<endl;
  cout<<"p+p fit is: "<<endl;
  gr5->Fit("expo2");
  gr5->GetFunction("expo2")->SetLineColor(kBlack);
  gr5->GetFunction("expo2")->SetLineWidth(1);
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
  //gr10->GetFunction("expo2")->SetLineColor(kBlue);
  //gr10->GetFunction("expo2")->SetLineWidth(1);
  cout<<"==================================================================="<<endl;
  cout<<"Au+Au BES 27 fit is: "<<endl;
  gr11->Fit("expo2");
  gr11->GetFunction("expo2")->SetLineColor(kGray);
  gr11->GetFunction("expo2")->SetLineWidth(1);
  cout<<"==================================================================="<<endl;
  cout<<"Au+Au BES 39 fit is: "<<endl;
  gr12->Fit("expo2");
  gr12->GetFunction("expo2")->SetLineColor(kViolet);
  gr12->GetFunction("expo2")->SetLineWidth(1);
  cout<<"==================================================================="<<endl;
  c1->SaveAs("EtaPi0_vs_pt.png");
  cout<<"complete";
  return;
}

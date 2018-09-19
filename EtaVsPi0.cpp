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

struct Data{
  float evp;
  float pt;
  float errR;
  float errPt;
};

void EtaVsPi0(){
  int p=0; //counter for iteration number
  int k=0; //errorcounter
  string type;
  string garbage;
  vector <Data> auaucent(0); //Au+Au for central collisions
  vector <Data> auauScent(0); //Au+Au for semi-cental Collisions
  vector <Data> auauPeri(0); //Au+Au for peripheral collisions
  vector <Data> dau(0); //d+Au minimum bias
  vector <Data> pp(0);  //p+p

  ifstream in;
  in.open("./EtaVsPi0.dat");

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


  cout<<endl<<"Inputs complete: Creating Graphs"<<endl<<endl<<endl;

  TCanvas *c1 = new TCanvas("c1","Au+Au Central", 200,10,700,500);
  Double_t x1[11], y1[11], ex1[11], ey1[11], iy1[11];
  Double_t xi1[11], yi1[11];
  Double_t x2[11], y2[11], ex2[11], ey2[11], iy2[11];
  Double_t xi2[11], yi2[11];
  Double_t x3[7], y3[7], ex3[7], ey3[7], iy3[7];
  Double_t xi3[7], yi3[7];
  Double_t x4[14], y4[14], ex4[14], ey4[14], iy4[14];
  Double_t xi4[14], yi4[14];
  Double_t x5[13], y5[13], ex5[13], ey5[13], iy5[13];
  Double_t xi5[13], yi5[13];
  Double_t x6[2], y6[2],x7[2],y7[2];
  Int_t n1=10,n2=10,n3=6,n4=13,n5=12,n6=2;

  TLegend* legend=new TLegend(.13,.7,.35,.9);
    legend->SetTextFont(72);
    legend->SetTextSize(.03);
    legend->SetFillColor(0);

  x7[0]=0;
  x7[1]=11.5;
  y7[0]=0;
  y7[1]=1.2;
  TGraph* gr7=new TGraph(n6,x7,y7);
    gr7->SetTitle("Eta/Pi0 Ratio at snn=200; pT; Eta/Pi0");
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
    gr5->SetMarkerSize(2);
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
    gr1->SetMarkerSize(2);
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
    gr2->SetMarkerSize(2);
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
    gr3->SetMarkerSize(2);
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
    gr4->SetMarkerSize(2);
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
    gr6->SetMarkerSize(2);
    gr6->SetMarkerColor(kBlack);
    gr6->Draw("L");
    legend->Draw();

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
  c1->SaveAs("EtaPi0_vs_pt.png");
  cout<<"complete";
  return;
}

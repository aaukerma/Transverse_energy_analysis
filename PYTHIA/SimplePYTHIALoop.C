//____________________________________________________________________
//
// Using Pythia6 with ROOT
// To make an event sample (of size 100) do
//
//    shell> root
//    root [0] .L pythiaExample.C
//    root [1] makeEventSample(1000)
//
// To start the tree view on the generated tree, do
//
//    shell> root
//    root [0] .L pythiaExample.C
//    root [1] showEventSample()
//
//
// The following session:
//    shell> root
//    root [0] .x pythiaExample.C(500)
// will execute makeEventSample(500) and showEventSample()
//
// Alternatively, you can compile this to a program
// and then generate 1000 events with
//
//    ./pythiaExample 1000
//
// To use the program to start the viewer, do
//
//    ./pythiaExample -1
//
// NOTE 1: To run this example, you must have a version of ROOT
// compiled with the Pythia6 version enabled and have Pythia6 installed.
// The statement gSystem->Load("$HOME/pythia6/libPythia6");  (see below)
// assumes that the directory containing the Pythia6 library
// is in the pythia6 subdirectory of your $HOME.  Locations
// that can specify this, are:
//
//  Root.DynamicPath resource in your ROOT configuration file
//    (/etc/root/system.rootrc or ~/.rootrc).
//  Runtime load paths set on the executable (Using GNU ld,
//    specified with flag `-rpath').
//  Dynamic loader search path as specified in the loaders
//    configuration file (On GNU/Linux this file is
//    etc/ld.so.conf).
//  For Un*x: Any directory mentioned in LD_LIBRARY_PATH
//  For Windows: Any directory mentioned in PATH
//
// NOTE 2: The example can also be run with ACLIC:
//  root > gSystem->Load("libEG");
//  root > gSystem->Load("$ROOTSYS/../pythia6/libPythia6"); //change to your setup
//  root > gSystem->Load("libEGPythia6");
//  root > .x pythiaExample.C+
//
//
//____________________________________________________________________
//
// Author: Christian Holm Christensen <cholm@hilux15.nbi.dk>
// Update: 2002-08-16 16:40:27+0200
// Copyright: 2002 (C) Christian Holm Christensen
// Copyright (C) 2006, Rene Brun and Fons Rademakers.
// For the licensing terms see $ROOTSYS/LICENSE.
//
#ifndef __CINT__
#include "TApplication.h"
#include "TPythia6.h"
#include "TFile.h"
#include "TError.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "Riostream.h"
#include <cstdlib>
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TMath.h"
#include "TMCParticle.h"
#include "THnSparse.h"
#include <iostream>
#include <vector>
using namespace std;
#endif

#define FILENAME   "pythia.root"
#define TREENAME   "tree"
#define BRANCHNAME "particles"
#define HISTNAME   "ptSpectra"
#define PDGNUMBER  211
class TH3F;
class TH1F;
class TPythia6;

struct KF_Code {
  Int_t KFc;
  TString name;
};

// This function just load the needed libraries if we're executing from
// an interactive session.
int test(){
return 0;
}

TH3F *CreateHistogram(char *name){
  //assoc pt, trig pt, dphi
  TH3F *histo = new TH3F(name,name,20,1.0,6.0,20,2.0,7.0,144,-TMath::Pi(),TMath::Pi());
  histo->GetXaxis()->SetTitle("p_{T}^{assoc}");
  histo->GetYaxis()->SetTitle("p_{T}^{trig}");
  histo->GetZaxis()->SetTitle("#Delta#phi");
  return histo;
}
TH1F *CreateTriggerHistogram(char *name){
  //assoc pt, trig pt, dphi
  TH1F *histo = new TH1F(name,name,20,2.0,7.0);
  histo->GetYaxis()->SetTitle("N_{trig}");
  histo->GetXaxis()->SetTitle("p_{T}^{trig}");
  return histo;
}
TH1F *CreateSpeciesHistogram(char *name,int bin,int low,int hi){
  TH1F *histo = new TH1F(name,name,bin,low,hi);
  histo->GetYaxis()->SetTitle("N_{trig}");
  histo->GetXaxis()->SetTitle("KFID");
  return histo;
}
TH3F *CreateParentHistogram(char *name){
  TH3F *histo = new TH3F(name,name,246,0,246,246,0,246,400,0,400);
  histo->GetYaxis()->SetTitle("cPKF");
  histo->GetXaxis()->SetTitle("cKF");
  histo->GetZaxis()->SetTitle("index");
  return histo;
}
TH2F *CreateEnergyHistogram(char *name,Int_t nEvents){
  TH2F *histo = new TH2F(name,name,1000,0,2,nEvents,0,nEvents);
  histo->GetYaxis()->SetTitle("number of entries");
  histo->GetXaxis()->SetTitle("Energy_Tpart");
  return histo;
}
TH1F *CreateEnergy1Histogram(char *name,Int_t nEvents){
  TH1F *histo = new TH1F(name,name,10000,0,15);
  histo->GetYaxis()->SetTitle("number of entries");
  histo->GetXaxis()->SetTitle("Energy_Tpart");
  return histo;
}
TH3F *CreateEventHistogram(char *name,Int_t nEvents){
  TH3F *histo = new TH3F(name,name,246,0,246,nEvents,0,nEvents,250,0,250);
  histo->GetYaxis()->SetTitle("event");
  histo->GetXaxis()->SetTitle("cKF");
  histo->GetZaxis()->SetTitle("index");
  return histo;
}
Int_t GetKFConversion(const Int_t kfc, const vector<KF_Code>& partname){
  Int_t i=0;
  Int_t j=0;
  Int_t k=0;
  Int_t n=1;
  Int_t KF;
  Int_t comp;
  Int_t comp1;
  TString name;
  KF=kfc;

  if (KF<0){
    KF=KF*-1;
    n=-1;
  }
  for (Int_t i=0;i<248;i++){
    comp=partname[i].KFc;
    if(comp==KF){
      j=comp;
      k=i;
      break;
    }
  }
  //NOTE: to include if -KF matters
  if (n==-1){
//    name=partname[k].name+" BAR";
    k=k*(-1);
  }
  //else
//    name=partname[k].name;
  return k;
}

Double_t dPhi(Double_t phi1, Double_t phi2) {
  Double_t deltaPhi;
  deltaPhi = phi1 - phi2;
  if (deltaPhi>(1.5*TMath::Pi())) deltaPhi-=2*(TMath::Pi());
  if (deltaPhi<(-0.5*TMath::Pi())) deltaPhi+=2*(TMath::Pi());
  return deltaPhi;
}

Double_t dEta(Double_t eta1, Double_t eta2) {
  Double_t deltaEta;
  deltaEta = eta1 - eta2;
  return deltaEta;
}

Double_t RelativeEPPhi(Double_t jetAng){ // function to calculate angle between jet and EP in the 1st quadrant (0,Pi/2)
  Double_t dphi =  jetAng;


  if( dphi<-1*TMath::Pi() ){
    dphi = dphi + 1*TMath::Pi();
  }

  if( dphi>1*TMath::Pi() ){
    dphi = dphi - 1*TMath::Pi();
  }

  if( (dphi>0) && (dphi<1*TMath::Pi()/2) ){
    // Do nothing! we are in quadrant 1
  }else if( (dphi>1*TMath::Pi()/2) && (dphi<1*TMath::Pi()) ){
    dphi = 1*TMath::Pi() - dphi;
  }else if( (dphi<0) && (dphi>-1*TMath::Pi()/2) ){
    dphi = fabs(dphi);
  }else if( (dphi<-1*TMath::Pi()/2) && (dphi>-1*TMath::Pi()) ){
    dphi = dphi + 1*TMath::Pi();
  }

  // test
  if( dphi < 0 || dphi > TMath::Pi()/2 ) cout<<" dPHI not in range [0, 0.5*Pi]!"<<endl;

  return dphi;   // dphi in [0, Pi/2]
}


void loadLibraries()
{
#ifdef __CINT__
  // Load the Event Generator abstraction library, Pythia 6
  // library, and the Pythia 6 interface library.
  gSystem->Load("libEG");
gSystem->Load("/data/rhip/alice/cnattras/pythia6/libPythia6");

//  gSystem->Load("$ROOTSYS/../pythia6/libPythia6"); //change to your setup
  gSystem->Load("libEGPythia6");
#endif
}

// nEvents is how many events we want.
int makeEventSample(Int_t nEvents, Int_t jobID, Int_t tune, Float_t SNN,Float_t trigEtaMax = 0.5, Float_t assocEtaMax = 0.9)
{
  Int_t l=0;
  vector <KF_Code> partname(0);
  /**************************************************************************
  choose energy
  {11,19,27,39,120,200,300,400,500,600,700,800,900,1000,1500,2000,2760,3000};
  **************************************************************************/
//    float SNN=7.7;
  int modSNN;
  if (SNN==7.7){
    modSNN=7;
  }
  else if (SNN==11.5){
    modSNN=11;
  }
  else if (SNN==19.6){
    modSNN=19;
  }
  else {
    modSNN=SNN;
  }
  cout<<modSNN<<endl;
  cout<<"I made it here, running "<<nEvents<<" events, job id "<<jobID<<", tune "<<tune<<", at "<<SNN<<"GeV"<<endl;
  char *filename = Form("hh%iGeV_outfile%i.root",modSNN,jobID);
  ifstream in;
  in.open("KF_Code.dat");
  while (in.good()){
    partname.push_back(KF_Code());
    in>>partname[l].KFc>>partname[l].name;
    l++;
  }
  //gSystem->Load("$PYTHIA6/libPythia6"); //change to your setup
  //gSystem->Load("libEGPythia6");
  // Load needed libraries
  //loadLibraries();


//-------------------------Creating my THNsparse---------------------------
  //1:	dPhi
  //2:	dEta
  //3:	pT trigger
  //4:  pT assoc
  //5:  dPhi trig

  Int_t nDim=5;
  Int_t *nbins= new Int_t[nDim];
  Double_t *xmin = new Double_t[nDim];
  Double_t *xmax = new Double_t[nDim];
  for (Int_t i=0; i<nDim; i++){
    nbins[i]=96;
    xmin[i]=0;
  }
  nbins[1]=10;//10 eta bins
  xmin[0]=-0.5*TMath::Pi();
  xmax[0]=1.5*TMath::Pi();
  xmin[1]=-(trigEtaMax+assocEtaMax);
  xmax[1]=(trigEtaMax+assocEtaMax);
  xmax[2]=20;
  nbins[2]=40;
  xmax[3]=10;
  nbins[3]=20;
  xmin[4]=0;//-0.5*TMath::Pi();
  xmax[4]=0.5*TMath::Pi();
  nbins[4]=24;
  Double_t *fill = new Double_t[nDim];
//NOTE:this is the THnSparse (DO NOT NEED)
  //THnSparseF *hhCorr = new THnSparseF("hh","hh correlations", nDim, nbins,xmin,xmax);
  //hhCorr->GetAxis(0)->SetTitle("#Delta#phi");
  //hhCorr->GetAxis(1)->SetTitle("#Delta#eta");
  //hhCorr->GetAxis(2)->SetTitle("p_{T}^{trigger}");
  //hhCorr->GetAxis(3)->SetTitle("p_{T}^{assoc}");
  //hhCorr->GetAxis(4)->SetTitle("#phi^{trigger}-#psi");

//Done

  TH1F * numTriggers = new TH1F("trigs", "Number of Triggers", 3, 0, 2);

//----------------Equation for Reaction Plane Dependence---------------
  TF1* trigPhi = new TF1("phiTrig","1+(2*TMath::Power(0.073086,2)*TMath::Cos(2*x))+(2*TMath::Power(0.092,2)*TMath::Cos(3*x))",-.5*TMath::Pi(),1.5*TMath::Pi());

  cout<<"I made it here line 144"<<endl;
  TString piplus = "pi+";
  TString piminus = "pi-";
  TString kplus = "K+";
  TString kminus = "K-";
  TString pplus = "p+";
  TString pminus = "pbar-";
  TString kshort = "K_S0";
  TString lambda = "Lambda0";
  TString antilambda = "Lambdabar0";

  // Create an instance of the Pythia event generator ...
  TPythia6* pythia = new TPythia6();
  cout<<"I made it here line 157"<<endl;
  //turn on decays
  //pythia->SetMSTJ(21,2);
  //pythia->SetMSTJ(22,1);
  //add switches for tunes
  //tune A/home/bs/Transverse_energy_analysis/PYTHIA/SimplePYTHIALoop.C
 /*      pythia->SetPARP(67,4.0);           // Regulates Initial State Radiation (value from best fit to D0 dijet analysis)
       pythia->SetMSTP(81,1);             // Double Gaussian Model
       pythia->SetMSTP(82,4);             // Double Gaussian Model
       pythia->SetPARP(82,2.0);           // [GeV]    PT_min at Ref. energy
       pythia->SetPARP(83,0.5);           // Core radius
       pythia->SetPARP(84,0.4);           // Core radius
       pythia->SetPARP(85,0.90) ;         // Regulates gluon prod. mechanism
       pythia->SetPARP(86,0.95);          // Regulates gluon prod. mechanism
       pythia->SetPARP(89,1800.);         // [GeV]   Ref. energy
       pythia->SetPARP(90,0.25);          // 2*epsilon (exponent in power law)
   */    //pythia->SetPARP(91,2.0);          // 2*epsilon (exponent in power law)
  pythia->SetMSTP(5,tune);//TuneA


//   //adds diffractive events
//   pythia->SetMSEL(0);
//   pythia->SetMSUB(92,1);
//   pythia->SetMSUB(93,1);
//   pythia->SetMSUB(94,1);
//   pythia->SetMSUB(95,1);
//may be important but unlike


//BEN THIS LINE IS IMPORTANT
  pythia->SetMSTJ(22,2);//decays unstable particles

  UInt_t seed = (jobID+1)*17;
  //UInt_t seed = rand()%900000000;
  //cout<<seed<<endl; //gives seed
  if( (seed>=0) && (seed<=900000000) ) {
    pythia->SetMRPY(1, seed);                   // set seed
    pythia->SetMRPY(2, 0);                      // use new seed
    pythia->SetPARP(2,7.0);                     // lowest collision energy
    cout<<"Random Seed : "<<seed<<endl;
  } else {cout << "error: time " << seed << " is not valid" << endl; exit(2);}


  // ... and initialise it to run p+p at sqrt(200) GeV in CMS
  pythia->Initialize("cms", "p", "p", SNN); //NOTE:in GeV
  //pythia->Dump();
  // Open an output file

  cout<<"I made it here line 180"<<endl;
  TFile* file = TFile::Open(filename, "RECREATE");
  if (!file || !file->IsOpen()) {
    Error("makeEventSample", "Couldn;t open file %s", filename);
    return 1;
  }

  TFile *outfile = file;//new TFile(filename,"RECREATE");
    //TH3F *hUnidentifiedCorrelations = CreateHistogram("hUnidentifiedCorrelations");
    //TH3F *hK0Correlations = CreateHistogram("hK0Correlations");
    //TH3F *hLambdaCorrelations = CreateHistogram("hLambdaCorrelations");
    //TH3F *hK0AssocCorrelations = CreateHistogram("hK0AssocCorrelations");
    //TH3F *hLambdaAssocCorrelations = CreateHistogram("hLambdaAssocCorrelations");


    //TH1F *hUnidentifiedTriggers = CreateTriggerHistogram("hUnidentifiedTriggers");
    //TH1F *hK0Triggers = CreateTriggerHistogram("hK0Triggers");
    //TH1F *hLambdaTriggers = CreateTriggerHistogram("hLambdaTriggers");


    //TH3F *hPiCorrelations = CreateHistogram("hPiCorrelations");
    //TH3F *hPiAssocCorrelations = CreateHistogram("hPiAssocCorrelations");
    //TH1F *hPiTriggers = CreateTriggerHistogram("hPiTriggers");
    //TH3F *hProtonCorrelations = CreateHistogram("hProtonCorrelations");
    //TH3F *hProtonAssocCorrelations = CreateHistogram("hProtonAssocCorrelations");
    //TH1F *hProtonTriggers = CreateTriggerHistogram("hProtonTriggers");
    //TH3F *hKCorrelations = CreateHistogram("hKCorrelations");
    //TH3F *hKAssocCorrelations = CreateHistogram("hKAssocCorrelations");
    //TH1F *hKTriggers = CreateTriggerHistogram("hKTriggers");

    //TH1F *hSPALL =CreateSpeciesHistogram("hSPALL",19804440,-9902220,9902220);
/*
    TH3F *hPar = CreateParentHistogram("hPar");
    TH2F *hEnergy = CreateEnergyHistogram("hEnergy",nEvents);
    TH2F *hETAll = CreateEnergyHistogram("hETAll",nEvents);
    TH2F *hETpip = CreateEnergyHistogram("h ET_Pi+",nEvents);
    TH2F *hETpim = CreateEnergyHistogram("h ET_Pi-",nEvents);
    */
    TH2F *hETpi0 = CreateEnergyHistogram("h ET_Pi0",nEvents);
    /*
    TH2F *hETKp = CreateEnergyHistogram("h ET_K+",nEvents);
    TH2F *hETKm = CreateEnergyHistogram("h ET_K-",nEvents);
    TH2F *hETK0 = CreateEnergyHistogram("h ET_K0",nEvents);
    TH2F *hETKL = CreateEnergyHistogram("h ET_KL",nEvents);
    TH2F *hETKS = CreateEnergyHistogram("h ET_KS",nEvents);
    TH2F *hETEta = CreateEnergyHistogram("h ET_Eta",nEvents);
    */
    TH2F *hETOmega = CreateEnergyHistogram("h ET_omega",nEvents);
    /*
    TH2F *hETOmegaP = CreateEnergyHistogram("h ET_Omega+",nEvents);
    TH2F *hETOmegaM = CreateEnergyHistogram("h ET_Omega-",nEvents);
    TH2F *hETLambda0 = CreateEnergyHistogram("h ET_Lambda0",nEvents);
    TH2F *hETLambdaBar0 = CreateEnergyHistogram("h ET_LambdaBar0",nEvents);
    TH2F *hETRhoP = CreateEnergyHistogram("h ET_Rho+",nEvents);
    TH2F *hETRho0 = CreateEnergyHistogram("h ET_Rho0",nEvents);
    TH2F *hETSigmaP = CreateEnergyHistogram("h ET_Sigma+",nEvents);
    TH2F *hETSigmaM = CreateEnergyHistogram("h ET_Simga-",nEvents);
    TH2F *hETSigma0 = CreateEnergyHistogram("h ET_Sigma0",nEvents);
    TH2F *hETDeltaPP = CreateEnergyHistogram("h ET_Delta++",nEvents);
    TH2F *hETDeltaP = CreateEnergyHistogram("h ET_Delta+",nEvents);
    TH2F *hETDeltaM = CreateEnergyHistogram("h ET_Delta-",nEvents);
    TH2F *hETDelta0 = CreateEnergyHistogram("h ET_Delta0",nEvents);
    TH2F *hETXiP = CreateEnergyHistogram("h ET_Xi+",nEvents);
    TH2F *hETXiM = CreateEnergyHistogram("h ET_Xi-",nEvents);
    TH2F *hETXi0 = CreateEnergyHistogram("h ET_Xi0",nEvents);
    TH2F *hETgamma = CreateEnergyHistogram("h ET_gamma",nEvents);
    TH2F *hETp = CreateEnergyHistogram("h ET_p",nEvents);
    TH2F *hETn = CreateEnergyHistogram("h ET_n",nEvents);
    TH2F *hETp_ = CreateEnergyHistogram("h ET_pBar",nEvents);
    TH2F *hETn_ = CreateEnergyHistogram("h ET_nBar",nEvents);
    TH2F *hETmu = CreateEnergyHistogram("h ET_mu",nEvents);
    TH2F *hETmu_ = CreateEnergyHistogram("h ET_muBar",nEvents);
    TH2F *hETe = CreateEnergyHistogram("h ET_e",nEvents);
    TH2F *hETe_ = CreateEnergyHistogram("h ET_eBar",nEvents);
    */
    TH1F *homegaVSpi0 = CreateEnergy1Histogram("h ET_omega/pi0",nEvents);
    //TH3F *hEve = CreateEventHistogram("hEve",nEvents);
    THStack hs("hs","test");

//NOTE: VERY IMPORTANT
    TH1F *hNEvents = new TH1F("hNEvents","Number of events",1,0,1.0);
    hNEvents->GetYaxis()->SetTitle("N_{events}");
    hNEvents->GetXaxis()->SetTitle("no title");


    Float_t lamTrigPt[3][25] = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
    Float_t k0TrigPt[3][25] = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
    Float_t hTrigPt[3][25] = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
    Float_t piTrigPt[3][25] = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
    Float_t kTrigPt[3][25] = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
    Float_t pTrigPt[3][25] = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
    Int_t hPartNumber[25] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    Int_t hPiPartNumber[25] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    Int_t hKPartNumber[25] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    Int_t hPPartNumber[25] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    Int_t nLamTrig, nK0Trig, nHTrig, nPiTrig, nKTrig, nPTrig;
    Int_t maxNtrig = 25;



  // ... and register a the cache of pythia on a branch (It's a
  // TClonesArray of TMCParticle objects. )
  TClonesArray* particles = (TClonesArray*)pythia->GetListOfParticles();
  //tree->Branch(BRANCHNAME, &particles);
  cout<<"I made it here line 235"<<endl;
  // Now we make some events
  //counters
  Int_t cpip=0;
  Int_t cpim=0;
  Int_t cpi0=0;
  Int_t cKp=0;
  Int_t cKm=0;
  Int_t cKL=0;
  Int_t cKS=0;
  Int_t cEta=0;
  Int_t cOmega=0; //THIS IS omega not Omega
  Int_t cOmegaP=0;
  Int_t cOmegaM=0;
  Int_t cLambda0=0;
  Int_t cLambdaBar0=0;
  Int_t cRhoP=0;
  Int_t cRho0=0;
  Int_t cSigmaP=0;
  Int_t cSigmaM=0;
  Int_t cSigma0=0;
  Int_t cDeltaPP=0;
  Int_t cDeltaP=0;
  Int_t cDeltaM=0;
  Int_t cDelta0=0;
  Int_t cXiM=0;
  Int_t cXiP=0;
  Int_t cXi0=0;
  Int_t cgamma=0;
  Int_t cp=0;
  Int_t cn=0;
  Int_t cp_=0;
  Int_t cn_=0;
  Int_t cmu=0;
  Int_t cmu_=0;
  Int_t ce_=0;
  Int_t ce=0;



  for (Int_t event = 0; event < nEvents; event++) {
    Double_t ETAll=0;
    Double_t ETpip=0;
    Double_t ETpim=0;
    Double_t ETpi0=0;
    Double_t ptpi0=0;
    Double_t ETKp=0;
    Double_t ETKm=0;
    Double_t ETK0=0;
    Double_t ETKL=0;
    Double_t ETKS=0;
    Double_t ETEta=0;
    Double_t ETOmega=0; //THIS IS omega not Omega
    Double_t ptome=0;
    Double_t ETOmegaP=0;
    Double_t ETOmegaM=0;
    Double_t ETLambda0=0;
    Double_t ETLambdaBar0=0;
    Double_t ETRhoP=0;
    Double_t ETRho0=0;
    Double_t ETSigmaP=0;
    Double_t ETSigmaM=0;
    Double_t ETSigma0=0;
    Double_t ETDeltaPP=0;
    Double_t ETDeltaP=0;
    Double_t ETDeltaM=0;
    Double_t ETDelta0=0;
    Double_t ETXiM=0;
    Double_t ETXiP=0;
    Double_t ETXi0=0;
    Double_t ETgamma=0;
    Double_t ETp=0;
    Double_t ETn=0;
    Double_t ETp_=0;
    Double_t ETn_=0;
    Double_t ETmu=0;
    Double_t ETmu_=0;
    Double_t ETe_=0;
    Double_t ETe=0;

    //ETpart/ETALL
    Double_t pETpip=0;
    Double_t pETpim=0;
    Double_t pETpi0=0;
    Double_t pETKp=0;
    Double_t pETKm=0;
    Double_t pETKL=0;
    Double_t pETKS=0;
    Double_t pETEta=0;
    Double_t pETOmega=0; //THIS IS omega not Omega
    Double_t pETOmegaP=0;
    Double_t pETOmegaM=0;
    Double_t pETLambda0=0;
    Double_t pETLambdaBar0=0;
    Double_t pETRhoP=0;
    Double_t pETRho0=0;
    Double_t pETSigmaP=0;
    Double_t pETSigmaM=0;
    Double_t pETSigma0=0;
    Double_t pETDeltaPP=0;
    Double_t pETDeltaP=0;
    Double_t pETDeltaM=0;
    Double_t pETDelta0=0;
    Double_t pETXiM=0;
    Double_t pETXiP=0;
    Double_t pETXi0=0;
    Double_t pETgamma=0;
    Double_t pETp=0;
    Double_t pETn=0;
    Double_t pETp_=0;
    Double_t pETn_=0;
    Double_t pETmu=0;
    Double_t pETmu_=0;
    Double_t pETe_=0;
    Double_t pETe=0;

    // Show how far we got every 100'th event.
    if (event % 100 == 0)
      cout <<"Event # " << event <<endl;
    // Make one event.
    pythia->GenerateEvent();

      hNEvents->Fill(0.5);
      nLamTrig = 0;
      nK0Trig = 0;
      nHTrig = 0;
      nPiTrig = 0;
      nKTrig = 0;
      nPTrig = 0;
      for(int i=0;i<25;i++){
	for(int j=0;j<3;j++){
	  lamTrigPt[j][i] = 0;
	  k0TrigPt[j][i] = 0;
	  hTrigPt[j][i] = 0;
	  piTrigPt[j][i] = 0;
	  kTrigPt[j][i] = 0;
	  pTrigPt[j][i] = 0;
	}
      }

      //pythia->Pylist(1); //DEBUG helper
      Int_t npart = particles->GetEntries();
      //printf("Analyse %d Particles\n", npart);
      for (Int_t part=0; part<npart; part++) {
	//TObject *object = particles->At(part);
	//cout<<"I am a "<<object->ClassName()<<endl;
	TMCParticle *MPart = (TMCParticle *) particles->At(part);
	Double_t pt= TMath::Sqrt(MPart->GetPx() * MPart->GetPx() + MPart->GetPy() * MPart->GetPy());
	Double_t phi = TMath::Pi()+TMath::ATan2(-MPart->GetPy(),-MPart->GetPx());
	Double_t p =  TMath::Sqrt(MPart->GetPx() * MPart->GetPx() + MPart->GetPy() * MPart->GetPy() +  MPart->GetPz() * MPart->GetPz());
	Double_t eta = 0.5*TMath::Log((p+MPart->GetPz())/(p-MPart->GetPz()));
	//cout<<"phi "<<phi<<" pt "<<pt<<" p "<<p<<" eta "<<eta<<endl;
	//if(MPart && MPart->GetPDG()) cout<<MPart->GetPDG()->GetName()<<endl;

  Int_t Ckf = MPart->GetKF();//convereted kf code to index
  Ckf=GetKFConversion(Ckf,partname);
  Int_t ind = part+1; //index
  //hEve->Fill(Ckf,event,ind);
  char N ='n';
  Float_t E= MPart ->GetEnergy();
  Float_t partE;
  Float_t px=MPart->GetPx();
  Float_t py=MPart->GetPy();
  Float_t pz=MPart->GetPz();
  if (pz<0){
    pz*=(-1);
  }
  Float_t pT= sqrt(pow(px,2)+pow(py,2));
  Float_t Theta = atan(pT/pz);
  Float_t p_tot=(pT)/sin(Theta);
  Float_t m=MPart->GetMass();
  Float_t Ei_TOT= sqrt(pow(p_tot,2)+pow(m,2));
  Float_t pseudorapidity=-log(tan(Theta/2));
  Float_t rapidity=(0.5)*log((E+pz)/(E-pz));
  //cout<<pseudorapidity<<" "<<rapidity<<endl;
  if (rapidity<0.1){
  //check part type for ET
  if ((Ckf==162)||(Ckf==135)||(Ckf==134)||(Ckf==136)||(Ckf==137)||(Ckf==138)||(Ckf==139)||(Ckf==133)||(Ckf==132)){
    partE = (Ei_TOT-m) * sin(Theta);
  }
  else if ((Ckf==-162)||(Ckf==-135)||(Ckf==-138)){
    partE = (Ei_TOT+m) * sin(Theta);
  }
  else {
    partE = Ei_TOT * sin(Theta);
  }

  Int_t mpartD = MPart->GetFirstChild();
  Int_t mpP = MPart->GetParent();
  Float_t mpL = MPart->GetLifetime();

  if (Ckf==58){ //pi+
    if ((mpartD==0)&&(mpP!=221)&&(mpP!=223)){ //counted if final particle AND not daughter of eta or omega NOTE: not subtracting pions from KS or KL
      ETpip+=partE;
      ETAll+=partE;
      cpip++;
  }}
  if (Ckf==-58){ //pi-
    if ((mpartD==0)&&(mpP!=221)&&(mpP!=223)){//counted if final particle AND not daughter of eta or omega NOTE: not subtracting pions from KS or KL
    ETpim+=partE;
    ETAll+=partE;
    cpim++;
  }}
  if (Ckf==68){ //pi0
    if ((mpP!=221)&&(mpP!=223)){//counted if not daughter of eta or omega NOTE: not subtracting pions from KS or KL (mpP!=73) (mpP!=74)
    ETpi0+=partE;
    ETAll+=partE;
    cpi0++;
    ptpi0+=p_tot;
  }}
  if (Ckf==60){ //K+
    if (mpartD==0){//counted if final particle
    ETKp+=partE;
    ETAll+=partE;
    cKp++;
  }}
  if (Ckf==-60){ //K-
    if (mpartD==0){//counted if final particle
    ETKm+=partE;
    ETAll+=partE;
    cKm++;
  }}
  if (Ckf==73){// KL
    if ((mpL>=100)||(mpartD==0)){//counted if final particle or lifetime > 1cm
    ETKL+=partE;
    ETAll+=partE;
    cKL++;
  }}
  if (Ckf==74){ //KS
    if ((mpL>=100)||(mpartD==0)){//counted if final particle or lifetime > 1cm
    ETKS+=partE;
    ETAll+=partE;
    cKS++;
  }}
  if (Ckf==69){ //eta
    ETEta+=partE;
    ETAll+=partE;
    cEta++;
  }
  if (Ckf==86){ // omega
    ETOmega+=partE;
    ETAll+=partE;
    cOmega++;
    ptome+=p_tot;
  }
  if (Ckf==-162){ // Omega+
    if (((mpL>=100)||(mpartD==0))&&(mpartD!=3122)&&(mpartD!=-3122)){ //included if: 1cm lifetime, final state particle, and does not produce a lambda0 or antilambda0
    ETOmegaP+=partE;
    ETAll+=partE;
    cOmegaP++;
  }}
  if (Ckf==162){ // Omega-
    if (((mpL>=100)||(mpartD==0))&&(mpartD!=3122)&&(mpartD!=-3122)){
    ETOmegaM+=partE;
    ETAll+=partE;
    cOmegaM++;
  }}
  if (Ckf==135){ //Lambda0
    if ((mpL>=100)||(mpartD==0)){
    ETLambda0+=partE;
    ETAll+=partE;
    cLambda0++;
  }}
  if (Ckf==-135){ //LambdaBar0
    if ((mpL>=100)||(mpartD==0)){
    ETLambdaBar0+=partE;
    ETAll+=partE;
    cLambdaBar0++;
  }}
  if (Ckf==137){ //Sigma+
    if (((mpL>=100)||(mpartD==0))&&(mpartD!=3122)&&(mpartD!=-3122)){
    ETSigmaP+=partE;
    ETAll+=partE;
    cSigmaP++;
  }}
  if (Ckf==134){ //Sigma-
    if (((mpL>=100)||(mpartD==0))&&(mpartD!=3122)&&(mpartD!=-3122)){
    ETSigmaM+=partE;
    ETAll+=partE;
    cSigmaM++;
  }}
  if (Ckf==136){ //Sigma0
    if (((mpL>=100)||(mpartD==0))&&(mpartD!=3122)&&(mpartD!=-3122)){
    ETSigma0+=partE;
    ETAll+=partE;
    cSigma0++;
  }}
  if (Ckf==138){ //Xi-
    if (((mpL>=100)||(mpartD==0))&&(mpartD!=3122)&&(mpartD!=-3122)){
    ETXiM+=partE;
    ETAll+=partE;
    cXiM++;
  }}
  if (Ckf==-138){ //Xi+
    if (((mpL>=100)||(mpartD==0))&&(mpartD!=3122)&&(mpartD!=-3122)){
    ETXiP+=partE;
    ETAll+=partE;
    cXiP++;
  }}
  if (Ckf==139){ //Xi0
    if (((mpL>=100)||(mpartD==0))&&(mpartD!=3122)&&(mpartD!=-3122)){
    ETXi0+=partE;
    ETAll+=partE;
    cXi0++;
  }}
  if (Ckf==17){ //gamma
    if ((mpartD==0)&&(mpP!=68)&&(mpP!=86)){
    ETgamma+=partE;
    ETAll+=partE;
    cgamma++;
  }}
  if (Ckf==133){ //proton
    if (mpartD==0){
    if ((part!=0)&&(part!=1)){
      ETp+=partE;
      ETAll+=partE;
      cp++;
    }}
  }
  if (Ckf==132){ //neutron
    if ((mpartD==0)&&(mpP!=2112)){
    ETn+=partE;
    ETAll+=partE;
    cn++;
  }}
  if (Ckf==-133){ //antiproton
    if (mpartD==0){
    if ((part!=0)&&(part!=1)){
      ETp_+=partE;
      ETAll+=partE;
      cp_++;
    }}
  }
  if (Ckf==-132){ //antineutron
    if ((mpartD==0)&&(mpP!=2112)){
    ETn_+=partE;
    ETAll+=partE;
    cn_++;
  }}
  if (Ckf==10){ //mu
    if (mpartD==0){
    ETmu+=partE;
    ETAll+=partE;
    cmu++;
  }}
  if (Ckf==11){ //muBar
    if (mpartD==0){
    ETmu_+=partE;
    ETAll+=partE;
    cmu_++;
  }}
  if (Ckf==9){ //positron
    if (mpartD==0){
    ETe_+=partE;
    ETAll+=partE;
    ce_++;
  }}
  if (Ckf==8){ //electron
    if (mpartD==0){
    ETe+=partE;
    ETAll+=partE;
    ce++;
  }}

} //end pseudorapidity cut
  if(pt>2.0 /*&& TMath::Abs(eta)<trigEtaMax*/){
	  TString mpart = MPart->GetName();
    Int_t mpartPKF=0;
    Int_t mpartDKF=0;
    Int_t mpartD2KF=0;
    Int_t mpartKF = MPart->GetKF();
    Int_t mpartKFCON=GetKFConversion(mpartKF,partname);
    Int_t mpartPKFCON;
    Float_t mpartE = MPart ->GetEnergy();
    Int_t mpartP = MPart ->GetParent();
    //hEnergy->Fill(mpartE);
    if (mpartP!=0){
      TMCParticle* parent = (TMCParticle *) particles->At(mpartP-1);
      mpartPKF = parent ->GetKF();
      if(mpartPKF<=0){
        mpartPKF=mpartPKF*(-1);
      }
      mpartPKFCON=GetKFConversion(mpartPKF,partname);
    //  hPar->Fill(mpartKFCON,mpartPKFCON,mpartP);
    }
    //cout<<mpartKF<<endl;
	  //Int_t mpartPDG  = MPart->GetPdgCode();
	  //cout<<"Part ID "<<mpart;
	  //cout<<" MPart name "<<MPart->GetName();
	  //if(MPart->GetPDG())cout<<" code "<<MPart->GetPDG()->PdgCode();
	  //cout<<endl;
    //cout<<mpartKF<<"\t "<<mpartP<<"\t "<<mpartPKF<<"\t "<<mpartE<<endl;
    //hSPALL->Fill(mpartKF);

	  if(mpart==kshort && nK0Trig<maxNtrig){//K0S
	    k0TrigPt[0][nK0Trig] = pt;
	    k0TrigPt[1][nK0Trig] = phi;
	    k0TrigPt[2][nK0Trig] = eta;
	    //cout<<"Kaon pt "<<pt<<" "<< k0TrigPt[0][nK0Trig] <<endl;
	    nK0Trig++;
	    //hK0Triggers->Fill(pt);
	    //printf("Particle %d\n", mpart);
	  }
	  if((mpart==lambda || mpart==antilambda) && nLamTrig<maxNtrig){//Lambda
	    //printf("Particle %d\n", mpart);
	    lamTrigPt[0][nLamTrig] =pt;
	    lamTrigPt[1][nLamTrig] = phi;
	    lamTrigPt[2][nLamTrig] = eta;
	    nLamTrig++;
	    //hLambdaTriggers->Fill(pt);
	  }
	  if((mpart==piplus || mpart==piminus || mpart==kplus || mpart==kminus || mpart==pplus || mpart==pminus)  && nHTrig<maxNtrig){//pi+- or p/pbar or K+-
	    hTrigPt[0][nHTrig] =pt;
	    hTrigPt[1][nHTrig] = phi;
	    hTrigPt[2][nHTrig] = eta;
	    hPartNumber[nHTrig] = part;
	    nHTrig++;
//	    cout<<"HADRON TRIG"<<endl;
	    //hUnidentifiedTriggers->Fill(pt);
	    if((mpart==piplus || mpart==piminus) && nPiTrig<maxNtrig){//pi
	      piTrigPt[0][nPiTrig] =pt;
	      piTrigPt[1][nPiTrig] = phi;
	      piTrigPt[2][nPiTrig] = eta;
	      hPiPartNumber[nPiTrig] = part;
	      nPiTrig++;
	      //hPiTriggers->Fill(pt);
	    }
	    if((mpart==pplus || mpart==pminus) && nPTrig<maxNtrig){//p/pbar
	      pTrigPt[0][nPTrig] =pt;
	      pTrigPt[1][nPTrig] = phi;
	      pTrigPt[2][nPTrig] = eta;
	      hPPartNumber[nPTrig] = part;
	      nPTrig++;
	      //hProtonTriggers->Fill(pt);
	    }
	    if((mpart==kplus || mpart==kminus) && nKTrig<maxNtrig){//K+-
	      kTrigPt[0][nKTrig] =pt;
	      kTrigPt[1][nKTrig] = phi;
	      kTrigPt[2][nKTrig] = eta;
	      hKPartNumber[nKTrig] = part;
	      nKTrig++;
	      //hKTriggers->Fill(pt);
	    }
	  }
	}
      }//end particle loop

  numTriggers->Fill(1,nPiTrig);
  numTriggers->Fill(1,nPTrig);
  numTriggers->Fill(1,nKTrig);

      if(nLamTrig>0 || nK0Trig>0 || nHTrig>0){//if any of the trigger particles have triggers, then and only then do we associate with particles
	for (Int_t part=0; part<npart; part++) {
	  TMCParticle *MPart = (TMCParticle *) particles->At(part);
	  Double_t pt= TMath::Sqrt(MPart->GetPx() * MPart->GetPx() + MPart->GetPy() * MPart->GetPy());
	  Double_t phi = TMath::Pi()+TMath::ATan2(-MPart->GetPy(),-MPart->GetPx());
	  Double_t p =  TMath::Sqrt(MPart->GetPx() * MPart->GetPx() + MPart->GetPy() * MPart->GetPy() +  MPart->GetPz() * MPart->GetPz());
	  Double_t eta = 0.5*TMath::Log((p+MPart->GetPz())/(p-MPart->GetPz()));
	  if(pt>1.0 && TMath::Abs(eta)<assocEtaMax){
	    //Int_t mpart  = MPart->GetPdgCode();
	    TString mpart = MPart->GetName();
	    if((mpart==piplus || mpart==piminus || mpart==kplus || mpart==kminus || mpart==pplus || mpart==pminus)){//for all charged pi/K/p
/*	      for(int i=0;i<nK0Trig;i++){
		if(k0TrigPt[0][i]>pt){
			//Need these
//	          cout<<"K0 TRIG"<<endl;
		  Double_t *fill = new Double_t[nDim];
                  fill[0]=dPhi(k0TrigPt[1][i],phi);
		  fill[1]=dEta(k0TrigPt[2][i],eta);
		  fill[2]=k0TrigPt[0][i];
		  fill[3]=pt;
		  hhCorr->Fill(fill);
		  hK0Correlations->Fill(pt,k0TrigPt[0][i],dPhi(k0TrigPt[1][i],phi));
		}
	      }
	      for(int i=0;i<nLamTrig;i++){
		if(lamTrigPt[0][i]>pt){
			//Need these
//	          cout<<"LAMBDA TRIG"<<endl;
		  Double_t *fill = new Double_t[nDim];
                  fill[0]=dPhi(lamTrigPt[1][i],phi);
		  fill[1]=dEta(lamTrigPt[2][i],eta);
		  fill[2]=lamTrigPt[0][i];
		  fill[3]=pt;
		  hhCorr->Fill(fill);
		  hLambdaCorrelations->Fill(pt,lamTrigPt[0][i],dPhi(lamTrigPt[1][i],phi));
		}
	      }
	      for(int i=0;i<nHTrig;i++){
	//	if(hTrigPt[0][i]>pt){
		//	cout<<"ENTERED HTRIG LOOP"<<endl;
		  if(part!=hPartNumber[i]){//don't correlate with itself.  This is only necessary for h-h correlations
		    hUnidentifiedCorrelations->Fill(pt,hTrigPt[0][i],dPhi(hTrigPt[1][i],phi));
		//	cout<<"NOT CORRELATING WITH ITSELF"<<endl;
		    if(mpart==piplus || mpart==piminus){//pi assoc
			//Need these
	          //    cout<<"PI TRIG"<<endl;
		      Double_t *fill = new Double_t[nDim];
                      fill[0]=dPhi(piTrigPt[1][i],phi);
      		      fill[1]=dEta(piTrigPt[2][i],eta);
		      fill[2]=piTrigPt[0][i];
		      fill[3]=pt;
		      hhCorr->Fill(fill);
		      hPiAssocCorrelations->Fill(pt,hTrigPt[0][i],dPhi(hTrigPt[1][i],phi));
		    }
		    if(mpart==pplus || mpart==pminus){//p/pbar assoc
			//Need these
	      	    //  cout<<"P TRIG"<<endl;
		      Double_t *fill = new Double_t[nDim];
                      fill[0]=dPhi(pTrigPt[1][i],phi);
      		      fill[1]=dEta(pTrigPt[2][i],eta);
		      fill[2]=pTrigPt[0][i];
		      fill[3]=pt;
		      hhCorr->Fill(fill);
		      hProtonAssocCorrelations->Fill(pt,hTrigPt[0][i],dPhi(hTrigPt[1][i],phi));
		    }
		    if(mpart==kplus || mpart==kminus){//K+- assoc
			//Need these
	     	      //cout<<"K TRIG"<<endl;
		      Double_t *fill = new Double_t[nDim];
                      fill[0]=dPhi(kTrigPt[1][i],phi);
      		      fill[1]=dEta(kTrigPt[2][i],eta);
		      fill[2]=kTrigPt[0][i];
		      fill[3]=pt;
		      hhCorr->Fill(fill);
		      hProtonAssocCorrelations->Fill(pt,hTrigPt[0][i],dPhi(hTrigPt[1][i],phi));
		      hKAssocCorrelations->Fill(pt,hTrigPt[0][i],dPhi(hTrigPt[1][i],phi));
		    }
		  }
		//}
	      }
*/	      for(int i=0;i<nPiTrig;i++){
//TRIGGER PHI CALCULATION HERE
  Double_t PhiT=RelativeEPPhi(trigPhi->GetRandom());
		if(piTrigPt[0][i]>pt){
		  if(part!=hPiPartNumber[i]){//don't correlate with itself.  This is only necessary for h-h correlations
			//Need these
                      fill[0]=dPhi(piTrigPt[1][i],phi);
      		      fill[1]=dEta(piTrigPt[2][i],eta);
		      fill[2]=piTrigPt[0][i];
		      fill[3]=pt;
		      fill[4]=PhiT;
		      //hhCorr->Fill(fill);
		      //hPiCorrelations->Fill(pt,piTrigPt[0][i],dPhi(piTrigPt[1][i],phi));
		  }
		}
	      }
	      for(int i=0;i<nPTrig;i++){
//TRIGGER PHI CALCULATION HERE
		Double_t PhiT=RelativeEPPhi(trigPhi->GetRandom());
		if(pTrigPt[0][i]>pt){
		  if(part!=hPPartNumber[i]){//don't correlate with itself.  This is only necessary for h-h correlations
			//Need these
		      Double_t *fill = new Double_t[nDim];
                      fill[0]=dPhi(pTrigPt[1][i],phi);
      		      fill[1]=dEta(pTrigPt[2][i],eta);
		      fill[2]=pTrigPt[0][i];
		      fill[3]=pt;
		      fill[4]=PhiT;
		      //hhCorr->Fill(fill);
		      //hProtonCorrelations->Fill(pt,pTrigPt[0][i],dPhi(pTrigPt[1][i],phi));
		  }
		}
	      }
	      for(int i=0;i<nKTrig;i++){
//TRIGGER PHI CALCULATION HERE
		Double_t PhiT=RelativeEPPhi(trigPhi->GetRandom());
		if(kTrigPt[0][i]>pt){
		  if(part!=hKPartNumber[i]){//don't correlate with itself.  This is only necessary for h-h correlations
			//Need these
		      Double_t *fill = new Double_t[nDim];
                      fill[0]=dPhi(kTrigPt[1][i],phi);
      		      fill[1]=dEta(kTrigPt[2][i],eta);
		      fill[2]=kTrigPt[0][i];
		      fill[3]=pt;
		      fill[4]=PhiT;
		      //hhCorr->Fill(fill);
		  }
		}
	      }
	    }
	  }

	}//end particle loop
      }
      Double_t omegaET=0;
      Double_t pi0ET=0;
      Double_t pi0pt=0;
      Double_t omept=0;
      Double_t ptomevptpi0=0;
      if ((ETpi0!=0)&&(ETOmega!=0)){
        hETpi0->Fill(ETpi0,cpi0);
        hETOmega->Fill(ETOmega,cOmega);
        pi0pt=ptpi0/cpi0;
        omept=ptome/cOmega;
        ptomevptpi0=(omept)/(pi0pt);
        homegaVSpi0->Fill(ptomevptpi0);
      }
      /*
      hETAll->Fill(ETAll,cpip);
      //cout<<"TOTAL: "<<ETAll<<endl;
      if (ETpip!=0){
      pETpip=ETpip/ETAll;
      hETpip->Fill(pETpip,cpip);

        //cout<<"pip "<<ETpip<<endl;
      }
      if (ETpim!=0){
      pETpim=ETpim/ETAll;
      hETpim->Fill(pETpim,cpim);
        //cout<<"pim "<<ETpim<<endl;
      }
      if (ETpi0!=0){
      pETpi0=ETpi0/ETAll;
      hETpi0->Fill(pETpi0,cpi0);
        //cout<<"pi0 "<<ETpi0<<endl;
      }
      if (ETKp!=0){
      pETKp=ETKp/ETAll;
      hETKp->Fill(pETKp,cKp);
        //cout<<"Kp "<<ETKp<<endl;
      }
      if (ETKm!=0){
      pETKm=ETKm/ETAll;
      hETKm->Fill(pETKm,cKm);
        //cout<<"Kp "<<ETKp<<endl;
      }
      if (ETKL!=0){
      pETKL=ETKL/ETAll;
      hETKL->Fill(pETKL,cKL);
        //cout<<"KL "<<ETKL<<endl;
      }
      if (ETKS!=0){
      pETKS=ETKS/ETAll;
      hETKS->Fill(pETKS,cKS);
        //cout<<"KS "<<ETKS<<endl;
      }
      if (ETEta!=0){
      pETEta=ETEta/ETAll;
      hETEta->Fill(pETEta,cEta);
        //cout<<"Eta "<<ETEta<<endl;
      }
      if (ETOmega!=0){
      pETOmega=ETOmega/ETAll;
      hETOmega->Fill(pETOmega,cOmega);
        //cout<<"Omega "<<ETOmega<<endl;
      }
      if (ETOmegaP!=0){
      pETOmegaP=ETOmegaP/ETAll;
      hETOmegaP->Fill(pETOmegaP,cOmegaP);
        //cout<<"OmegaM "<<ETOmegaM<<endl;
      }
      if (ETOmegaM!=0){
      pETOmegaM=ETOmegaM/ETAll;
      hETOmegaM->Fill(pETOmegaM,cOmegaM);
        //cout<<"OmegaM "<<ETOmegaM<<endl;
      }
      if (ETLambda0!=0){
      pETLambda0=ETLambda0/ETAll;
      hETLambda0->Fill(pETLambda0,cLambda0);
        //cout<<"Lambda0 "<<ETLambda0<<endl;
      }
      if (ETLambdaBar0!=0){
      pETLambdaBar0=ETLambdaBar0/ETAll;
      hETLambdaBar0->Fill(pETLambdaBar0,cLambdaBar0);
        //cout<<"Lambda0 "<<ETLambda0<<endl;
      }
      if (ETSigmaP!=0){
      pETSigmaP=ETSigmaP/ETAll;
      hETSigmaP->Fill(pETSigmaP,cSigmaP);
        //cout<<"SigmaP "<<ETSigmaP<<endl;
      }
      if (ETSigmaM!=0){
      pETSigmaM=ETSigmaM/ETAll;
      hETSigmaM->Fill(pETSigmaM,cSigmaM);
        //cout<<"SigmaM "<<ETSigmaM<<endl;
      }
      if (ETSigma0!=0){
      pETSigma0=ETSigma0/ETAll;
      hETSigma0->Fill(pETSigma0,cSigma0);
        //cout<<"Sigma0 "<<ETSigma0<<endl;
      }
      if (ETXiP!=0){
      pETXiP=ETXiP/ETAll;
      hETXiP->Fill(pETXiP,cXiP);
        //cout<<"Xi- "<<ETXiM<<endl;
      }
      if (ETXiM!=0){
      pETXiM=ETXiM/ETAll;
      hETXiM->Fill(pETXiM,cXiM);
        //cout<<"Xi- "<<ETXiM<<endl;
      }
      if (ETXi0!=0){
      pETXi0=ETXi0/ETAll;
      hETXi0->Fill(pETXi0,cXi0);
        //cout<<"Xi0 "<<ETXiM<<endl;
      }
      if (ETgamma!=0){
      pETgamma=ETgamma/ETAll;
      hETgamma->Fill(pETgamma,cgamma);
        //cout<<"gamma "<<ETgamma<<endl;
      }
      if (ETp!=0){
      pETp=ETp/ETAll;
      hETp->Fill(pETp,cp);
        //cout<<"p "<<ETp<<endl;
      }
      if (ETn!=0){
      pETn=ETn/ETAll;
      hETn->Fill(pETn,cn);
        //cout<<"n "<<ETn<<endl;
      }
      if (ETp_!=0){
      pETp_=ETp_/ETAll;
      hETp_->Fill(pETp_,cp_);
        //cout<<"p "<<ETp<<endl;
      }
      if (ETn_!=0){
      pETn_=ETn_/ETAll;
      hETn_->Fill(pETn_,cn_);
        //cout<<"n "<<ETn<<endl;
      }
      if (ETmu!=0){
      pETmu=ETmu/ETAll;
      hETmu->Fill(pETmu,cmu);
        //cout<<"mu "<<ETmu<<endl;
      }
      if (ETmu_!=0){
      pETmu_=ETmu_/ETAll;
      hETmu_->Fill(pETmu_,cmu_);//set marker color/style
        //cout<<"muBAR "<<ETmu_<<endl;
      }
      if (ETe!=0){
      pETe=ETe/ETAll;
      hETe->Fill(pETe,ce);
        //cout<<"e "<<ETe<<endl;
      }
      if (ETe_!=0){
      pETe_=ETe_/ETAll;
      hETe_->Fill(pETe_,ce_);
        //cout<<"eBAR "<<ETe_<<endl;
      }*/

  }/*
  Double_t ETAll=0;
  Double_t ETpip=0;
  Double_t ETpim=0;*/
  Double_t ETpi0=0;/*
  Double_t ETKp=0;
  Double_t ETKm=0;
  Double_t ETK0=0;
  Double_t ETKL=0;
  Double_t ETKS=0;*/
  Double_t ETEta=0;
  Double_t ETOmega=0;
  Double_t ptmean=0; //THIS IS omega not Omega
  /*Double_t ETOmegaP=0;
  Double_t ETOmegaM=0;
  Double_t ETLambda0=0;
  Double_t ETLambdaBar0=0;
  Double_t ETSigmaP=0;
  Double_t ETSigmaM=0;
  Double_t ETSigma0=0;
  Double_t ETXiM=0;
  Double_t ETXiP=0;
  Double_t ETXi0=0;
  Double_t ETgamma=0;
  Double_t ETp=0;
  Double_t ETn=0;
  Double_t ETmu=0;
  Double_t ETmu_=0;
  Double_t ETe_=0;
  Double_t ETe=0;*/
  //NOTE:write histograms to file
    //numTriggers->Write();
    //hhCorr->Write();
    //hUnidentifiedTriggers->Write();
    //hPiTriggers->Write();
    //hKTriggers->Write();
    //hK0Triggers->Write();
    //hLambdaTriggers->Write();
    //hProtonTriggers->Write();
    hNEvents->Write();
    //hSPALL->Write();
    //hPar->Write();
    /*
    ETAll=hETAll->GetMean();
    ETpip=hETpip->GetMean();
    ETpim=hETpim->GetMean();*/
    ETpi0=hETpi0->GetMean();
    /*ETKp=hETKp->GetMean();
    ETKm=hETKm->GetMean();
    ETK0=hETK0->GetMean();
    ETKL=hETKL->GetMean();
    ETKS=hETKS->GetMean();*/
    //ETEta=hETEta->GetMean();
    ETOmega=hETOmega->GetMean();
    /*ETOmegaM=hETOmegaM->GetMean();
    ETOmegaP=hETOmegaP->GetMean();
    ETLambda0=hETLambda0->GetMean();
    ETLambdaBar0=hETLambdaBar0->GetMean();
    ETSigmaP=hETSigmaP->GetMean();
    ETSigmaM=hETSigmaM->GetMean();
    ETSigma0=hETSigma0->GetMean();
    ETXiP=hETXiP->GetMean();
    ETXiM=hETXiM->GetMean();
    ETXi0=hETXi0->GetMean();
    ETgamma=hETgamma->GetMean();
    ETp=hETp->GetMean();
    ETn=hETn->GetMean();
    ETmu=hETmu->GetMean();
    ETmu_=hETmu_->GetMean();
    ETe=hETe->GetMean();
    ETe_=hETe_->GetMean();
    */
    ptmean=homegaVSpi0->GetMean();
    //Float_t Ratio=(ETOmega)/(ETpi0);
    //Float_t Ratio2=(ETEta*cEta)/(ETpi0*cpi0);
    // found average to be .889669
    //cout<<"ET Eta to Pi0 :"<<Ratio2<<endl;
    cout<<"PT Omega to Pi0 :"<<ptmean<<endl;
    //hETAll->Write();
    //hETpip->Write();
    //hETpim->Write();
    hETpi0->Write();
    //hETKp->Write();
    //hETKm->Write();
    //hETKL->Write();
    //hETKS->Write();
    //hETEta->Write();
    hETOmega->Write();
    homegaVSpi0->Write();
    //hETOmegaM->Write();
    //hETOmegaP->Write();
    //hETLambda0->Write();
    //hETLambdaBar0->Write();
    //hETSigmaP->Write();
    //hETSigmaM->Write();
    //hETSigma0->Write();
    //hETXiP->Write();
    //hETXiM->Write();
    //hETXi0->Write();
    //hETgamma->Write();
    //hETp->Write();
    //hETn->Write();
    //hETp_->Write();
    ///hETn_->Write();
    //hETmu->Write();
    //hETmu_->Write();
    //hETe->Write();
    //hETe_->Write();

    //hs.Write();
    //hEve->Write();
    outfile->Close();

  return 0;
}

// Show the Pt spectra, and start the tree viewer.
int showEventSample()
{
  // Load needed libraries
  loadLibraries();

  // Open the file
  TFile* file = TFile::Open(FILENAME, "READ");
  if (!file || !file->IsOpen()) {
    Error("showEventSample", "Couldn;t open file %s", FILENAME);
    return 1;
  }

  // Get the tree
  TTree* tree = (TTree*)file->Get(TREENAME);
  if (!tree) {
    Error("showEventSample", "couldn't get TTree %s", TREENAME);
    return 2;
  }

  // Start the viewer.
  tree->StartViewer();

  // Get the histogram
  TH1D* hist = (TH1D*)file->Get(HISTNAME);
  if (!hist) {
    Error("showEventSample", "couldn't get TH1D %s", HISTNAME);
    return 4;
  }

  // Draw the histogram in a canvas
  gStyle->SetOptStat(1);
  TCanvas* canvas = new TCanvas("canvas", "canvas");
  canvas->SetLogy();
  hist->Draw("e1");
  TF1* func = hist->GetFunction("expo");

  char expression[64];
  sprintf(expression,"T #approx %5.1f", -1000 / func->GetParameter(1));
  TLatex* latex = new TLatex(1.5, 1e-4, expression);
  latex->SetTextSize(.1);
  latex->SetTextColor(4);
  latex->Draw();

  return 0;
}

void SimplePYTHIALoop(Int_t n=1000, Int_t jobID=0, Int_t tune = 350,Float_t SNN = 2760,Float_t trigEtaMax = 0.5, Float_t assocEtaMax = 0.9) {
  makeEventSample(n,jobID,tune,SNN, trigEtaMax, assocEtaMax);
}

#ifndef __CINT__
int main(int argc, char** argv)
{
  TApplication app("app", &argc, argv);

  Int_t n = 100;
  if (argc > 1)
    n = strtol(argv[1], NULL, 0);

  int retVal = 0;
  if (n > 0)
    retVal = makeEventSample(n,0,0,7.7,0.5,0.9);
  else {
    retVal = showEventSample();
    app.Run();
  }

  return retVal;
}
#endif

//____________________________________________________________________
//
// EOF
//

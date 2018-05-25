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
int makeEventSample(Int_t nEvents, Int_t jobID, Int_t tune, Float_t trigEtaMax = 0.5, Float_t assocEtaMax = 0.9)
{
  cout<<"I made it here, running "<<nEvents<<" events, job id "<<jobID<<", tune "<<tune<<endl;
  char *filename = Form("hh_outfile%i.root",jobID);
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

  THnSparseF *hhCorr = new THnSparseF("hh","hh correlations", nDim, nbins,xmin,xmax);
  hhCorr->GetAxis(0)->SetTitle("#Delta#phi");
  hhCorr->GetAxis(1)->SetTitle("#Delta#eta");
  hhCorr->GetAxis(2)->SetTitle("p_{T}^{trigger}");
  hhCorr->GetAxis(3)->SetTitle("p_{T}^{assoc}");
  hhCorr->GetAxis(4)->SetTitle("#phi^{trigger}-#psi");

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
  //tune A
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

//BEN THIS LINE IS IMPORTANT
  pythia->SetMSTJ(22,2);//decays unstable particles


  UInt_t seed = (jobID+1)*17;
  if( (seed>=0) && (seed<=900000000) ) {
    pythia->SetMRPY(1, seed);                   // set seed
    pythia->SetMRPY(2, 0);                      // use new seed
    cout<<"Random Seed : "<<seed<<endl;
  } else {cout << "error: time " << seed << " is not valid" << endl; exit(2);}

  // ... and initialise it to run p+p at sqrt(200) GeV in CMS
  pythia->Initialize("cms", "p", "p", 2760);
  //pythia->Dump();
  // Open an output file

  cout<<"I made it here line 180"<<endl;
  TFile* file = TFile::Open(filename, "RECREATE");
  if (!file || !file->IsOpen()) {
    Error("makeEventSample", "Couldn;t open file %s", filename);
    return 1;
  }


  TFile *outfile = file;//new TFile(filename,"RECREATE");
    TH3F *hUnidentifiedCorrelations = CreateHistogram("hUnidentifiedCorrelations");
    TH3F *hK0Correlations = CreateHistogram("hK0Correlations");
    TH3F *hLambdaCorrelations = CreateHistogram("hLambdaCorrelations");
    TH3F *hK0AssocCorrelations = CreateHistogram("hK0AssocCorrelations");
    TH3F *hLambdaAssocCorrelations = CreateHistogram("hLambdaAssocCorrelations");


    TH1F *hUnidentifiedTriggers = CreateTriggerHistogram("hUnidentifiedTriggers");
    TH1F *hK0Triggers = CreateTriggerHistogram("hK0Triggers");
    TH1F *hLambdaTriggers = CreateTriggerHistogram("hLambdaTriggers");


    TH3F *hPiCorrelations = CreateHistogram("hPiCorrelations");
    TH3F *hPiAssocCorrelations = CreateHistogram("hPiAssocCorrelations");
    TH1F *hPiTriggers = CreateTriggerHistogram("hPiTriggers");
    TH3F *hProtonCorrelations = CreateHistogram("hProtonCorrelations");
    TH3F *hProtonAssocCorrelations = CreateHistogram("hProtonAssocCorrelations");
    TH1F *hProtonTriggers = CreateTriggerHistogram("hProtonTriggers");
    TH3F *hKCorrelations = CreateHistogram("hKCorrelations");
    TH3F *hKAssocCorrelations = CreateHistogram("hKAssocCorrelations");
    TH1F *hKTriggers = CreateTriggerHistogram("hKTriggers");



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
  for (Int_t event = 0; event < nEvents; event++) {
    // Show how far we got every 100'th event.
    if (event % 100 == 0)
      cout << "Event # " << event << endl;
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
	if(pt>2.0 && TMath::Abs(eta)<trigEtaMax){
	  TString mpart = MPart->GetName();
	  //Int_t mpart  = MPart->GetPdgCode();
	  //cout<<"Part ID "<<mpart;
	  //cout<<" MPart name "<<MPart->GetName();
	  //if(MPart->GetPDG())cout<<" code "<<MPart->GetPDG()->PdgCode();
	  //cout<<endl;
	  if(mpart==kshort && nK0Trig<maxNtrig){//K0S
	    k0TrigPt[0][nK0Trig] = pt;
	    k0TrigPt[1][nK0Trig] = phi;
	    k0TrigPt[2][nK0Trig] = eta;
	    //cout<<"Kaon pt "<<pt<<" "<< k0TrigPt[0][nK0Trig] <<endl;
	    nK0Trig++;
	    hK0Triggers->Fill(pt);
	    //printf("Particle %d\n", mpart);
	  }
	  if((mpart==lambda || mpart==antilambda) && nLamTrig<maxNtrig){//Lambda
	    //printf("Particle %d\n", mpart);
	    lamTrigPt[0][nLamTrig] =pt;
	    lamTrigPt[1][nLamTrig] = phi;
	    lamTrigPt[2][nLamTrig] = eta;
	    nLamTrig++;
	    hLambdaTriggers->Fill(pt);
	  }
	  if((mpart==piplus || mpart==piminus || mpart==kplus || mpart==kminus || mpart==pplus || mpart==pminus)  && nHTrig<maxNtrig){//pi+- or p/pbar or K+-
	    hTrigPt[0][nHTrig] =pt;
	    hTrigPt[1][nHTrig] = phi;
	    hTrigPt[2][nHTrig] = eta;
	    hPartNumber[nHTrig] = part;
	    nHTrig++;
//	    cout<<"HADRON TRIG"<<endl;
	    hUnidentifiedTriggers->Fill(pt);
	    if((mpart==piplus || mpart==piminus) && nPiTrig<maxNtrig){//pi
	      piTrigPt[0][nPiTrig] =pt;
	      piTrigPt[1][nPiTrig] = phi;
	      piTrigPt[2][nPiTrig] = eta;
	      hPiPartNumber[nPiTrig] = part;
	      nPiTrig++;
	      hPiTriggers->Fill(pt);
	    }
	    if((mpart==pplus || mpart==pminus) && nPTrig<maxNtrig){//p/pbar
	      pTrigPt[0][nPTrig] =pt;
	      pTrigPt[1][nPTrig] = phi;
	      pTrigPt[2][nPTrig] = eta;
	      hPPartNumber[nPTrig] = part;
	      nPTrig++;
	      hProtonTriggers->Fill(pt);
	    }
	    if((mpart==kplus || mpart==kminus) && nKTrig<maxNtrig){//K+-
	      kTrigPt[0][nKTrig] =pt;
	      kTrigPt[1][nKTrig] = phi;
	      kTrigPt[2][nKTrig] = eta;
	      hKPartNumber[nKTrig] = part;
	      nKTrig++;
	      hKTriggers->Fill(pt);
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
		      hhCorr->Fill(fill);
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
		      hhCorr->Fill(fill);
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
		      hhCorr->Fill(fill);
		  }
		}
	      }
	    }
	  }
	}//end particle loop
      }

  }

  //write histograms to file
    numTriggers->Write();
    hhCorr->Write();
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

void SimplePYTHIALoop(Int_t n=1000, Int_t jobID=0, Int_t tune = 350) {
  makeEventSample(n,jobID,tune);
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
    retVal = makeEventSample(n,0,0);
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


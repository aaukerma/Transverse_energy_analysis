/***************************************************************************
 *  Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                         *
 *  Author: ALICE OFFLICE.                                                 *
 *  Contributors are mentioned in the code where appropriate.              *
 *                                                                         *
 *  Permission to use, copy, modify and distribute this software and its   *
 *  documentation strictly for non-commercial purposes is hereby granted   *
 *  without fee, provided that the above copyright notice appears in all   *
 *  copies and that both the copyright notice and this permission notice   *
 *  appear in the supporting documentation. The authors make no claims     *
 *  about the suitability of this software for any purpose. It is          *
 *  provided "as is" without express or implied warranty.                  *
 *                                                                         *
 ***************************************************************************
     $Satyajit Jena || alien:sjena Sun Apr 21 14:05:19 CEST 2013$

       Sterring macro for fast production of MC events - Followed from
       the original macro of fastGenAmpt.C

     Implemented Generators: (Version 1.0: Sun Apr 21 14:05:19 CEST 2013)
     ----------------------------------------------------------------
      kPythia6,            kPythia8,               kPythia6D6T,
      kPythiaPerugia0,     kPythia6ATLAS,
      kPythiaJets,
      kPhojet,
      kDPMjet,             kDPMjet_pA,
      kAmptDefault,        kAmptStringMelting,          kAmptStringMeltingNoART,
      kAmptpA,                kAmptFlow,
      kAmptReducedFlow,

     FIXME:
     kPythia6ATLAS_Flat,
     kHijing,             kHijing2000,            kHijing_pA,


 ***************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TStopwatch.h>
#include <TDatime.h>
#include <TRandom.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include "TNamed.h"
#include <TParticle.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <TRandom3.h>
#include <TChain.h>
#include <TVirtualMC.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliHeader.h"
#include "AliStack.h"
#include "AliPDG.h"
#include "AliGenAmpt.h"
#include "TAmpt.h"
#include "AliDecayerPythia.h"
#include "AliConfig.h"
#include "AliGenerator.h"
#include "AliLog.h"
#include "AliGenHIJINGpara.h"
#include "AliGenHijing.h"
#include "AliGenCocktail.h"
#include "AliGenSlowNucleons.h"
#include "AliSlowNucleonModelExp.h"
#include "AliGenParam.h"
#include "AliGenMUONlib.h"
#include "AliGenSTRANGElib.h"
#include "AliGenMUONCocktail.h"
#include "AliGenCocktail.h"
#include "AliGenGeVSim.h"
#include "AliGeVSimParticle.h"
#include "AliGenPythia.h"
#include "AliGenPythiaPlus.h"
#include "AliPythia8.h"
//#include "AliGenDPMjet.h"
#endif

//___________________________________________________________________________
enum PDC06Proc_t {
  kPythia6,
  kPythia8,
  kPythia6D6T,
  kPythiaPerugia0,
  kPythia6ATLAS,
  kPythia6ATLAS_Flat,
  kPythiaJets,
  kPhojet,
  kHijing,
  kHijing2000,
  kHijing_pA,
  kDPMjet,
  kDPMjet_pA,
  kAmptDefault,
  kAmpt,
  kAmptpA,
  kAmptFlowStringMelting,
  kAmptStringMeltingNoART,
  kAmptFlow,
  kAmptReducedFlow,
  kRunMax
};

//___________________________________________________________________________
const char * pprRunName[] = {
  "kPythia6",
  "kPythia8",
  "kPythia6D6T",
  "kPythiaPerugia0",
  "kPythia6ATLAS",
  "kPythia6ATLAS_Flat",
  "kPythiaJets",
  "kPhojet",
  "kHijing",
  "kHijing2000",
  "kHijing_pA",
  "kDPMjet",
  "kDPMjet_pA",
  "kAmpt",
  "kAmptpA",
  "kAmptFlow",
  "kAmptReducedFlow"
};

enum PprTrigConf_t {kDefaultPPTrig, kDefaultPbPbTrig };
const char * pprTrigConfName[] = {"p-p","Pb-Pb"};

//___________________________________________________________________________
void ProcessEnvironmentVars();
class AliGenPythia;
AliGenerator *MbPythia();
AliGenerator *MyPythia8();
AliGenerator *MbPythiaTuneD6T();
AliGenerator *MbPythiaTunePerugia0();
AliGenerator *MbPythiaTuneATLAS();
AliGenerator *MbPythiaTuneATLAS_Flat();
AliGenerator *PythiaJets();
AliGenerator *MbPhojet();
AliGenerator *Hijing();
AliGenerator *Hijing2000();
AliGenerator *Hijing_pA(Bool_t kSlowN);
AliGenerator *DPMjet();
AliGenerator *DPMjet_pA(Bool_t fragments);
AliGenerator *Ampt();
AliGenerator *AmptpA();
AliGenerator* AmptFlow();
AliGenerator *AmptReducedFlow();
AliGenerator* AmptDefault();
AliGenerator* AmptStringMelting();
AliGenerator* AmptStringMeltingNoART();

//_________________________________________________________________________
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// Geterator, field, beam energy

static Double_t    pBeamEnergy = 4000.0;  // Used during pA runs
//static Double_t  energy  = 2.*pBeamEnergy*TMath::Sqrt(82./208.); //energy in CMS

static Double_t      energy    = 39.;
static PDC06Proc_t   proc      = kAmptDefault;
static Float_t       bMin      = 0.;
static Float_t       bMax      = 100.;
static PprTrigConf_t strig     = kDefaultPPTrig; // default pp trigger configuration


static Double_t  JpsiPol      = 0; // Jpsi polarisation
static Bool_t    JpsiHarderPt = kFALSE; // Jpsi harder pt spectrum (8.8 TeV)
static TString comment;
//static PprTrigConf_t strig    = kDefaultPbPbTrig; // default pp trigger configuration
TDatime dt;
static UInt_t seed    = dt.Get();

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//_________________________________________________________________________

TH3F *CreateHistogram(char *name){
  //assoc pt, trig pt, dphi
  TH3F *histo = new TH3F(name,name,20,1.0,6.0,20,2.0,7.0,144,-TMath::Pi(),TMath::Pi());
  histo->GetXaxis()->SetTitle("p_{T}^{assoc}");
  histo->GetYaxis()->SetTitle("p_{T}^{trig}");
  histo->GetZaxis()->SetTitle("#Delta#phi");
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
TH1F *CreateEnergyHistogram(char *name){
  TH1F *histo = new TH1F(name,name,1000,0,40);
  histo->GetYaxis()->SetTitle("Ratio");
  histo->GetXaxis()->SetTitle("Et");
  return histo;
}
TH2F *CreateEnergyHistogram2(char *name){
  TH2F *histo = new TH2F(name,name,1000,0,10,1000,0,1.1);
  histo->GetYaxis()->SetTitle("EtCh/Et");
  histo->GetXaxis()->SetTitle("EtCh");
  return histo;
}
TH3F *CreateEventHistogram(char *name,Int_t nEvents){
  TH3F *histo = new TH3F(name,name,246,0,246,nEvents,0,nEvents,250,0,250);
  histo->GetYaxis()->SetTitle("event");
  histo->GetXaxis()->SetTitle("cKF");
  histo->GetZaxis()->SetTitle("index");
  return histo;
}
TH1F *CreatePTHistogram(char *name){
  TH1F *histo = new TH1F(name,name,100,0,10);
  histo->GetYaxis()->SetTitle("number of entries");
  histo->GetXaxis()->SetTitle("Energy_Tpart/ETAll");
  return histo;
}
TH2F *CreatePTHistogram2(char *name){
  TH2F *histo = new TH2F(name,name,1000,0,100,1000,0,100);
  histo->GetYaxis()->SetTitle("assumption ratio");
  histo->GetXaxis()->SetTitle("pT");
  return histo;
}


void fastMcProduction39(Int_t nev = 300) {
  ProcessEnvironmentVars();

  gRandom->SetSeed(seed);
  cerr<<" Seed for random number generation = "<<seed<<endl;


/******************************************************************************
Require to insert all the histo creation parameters here. include all variables!!
******************************************************************************/

#if defined(__CINT__)
  gSystem->Load("liblhapdf");
  gSystem->Load("libEGPythia6");

  if (proc == kPythia6 || proc == kPhojet || proc == kDPMjet || proc==kDPMjet_pA) {
    gSystem->Load("libpythia6");        // Pythia 6.2
    gSystem->Load("libAliPythia6");     // ALICE specific implementations
  }

  if (proc == kHijing || proc == kHijing2000 || proc == kHijing_pA ) {
    gSystem->Load("libHIJING");
    gSystem->Load("libTHijing");
  }

  else if ( proc == kDPMjet || proc== kDPMjet_pA ) {
    gSystem->Load("libDPMJET");
    gSystem->Load("libTDPMjet");
  }

  else if (proc == kAmptDefault || kAmptFlowStringMelting || proc ==  kAmptStringMeltingNoART || proc == kAmptpA || proc == kAmptReducedFlow) {
    gSystem->Load("libampt");
    gSystem->Load("libTAmpt");
    gSystem->Load("libpythia6");
    gSystem->Load("libAliPythia6");
  }

  if (proc == kPythia8) {
    gSystem->Load("libpythia8");
    gSystem->Load("libAliPythia8");
    gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8145/xmldoc"));
    gSystem->Setenv("LHAPDF",      gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF"));
    gSystem->Setenv("LHAPATH",     gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF/PDFsets"));
  }
#endif


 AliGenerator* gener = 0x0;

 cout<<"Run type set to ------------- "<<pprRunName[proc]<<"   " << proc << "    " << kDPMjet_pA<< endl;

 if (proc == kPythia6) {
   Printf("<<<<<<<<<<<<<<<<<<< Processing Mb. PYTHIA >>>>>>>>>>>>>>>>>>>>");
   gener = MbPythia();
 }

 else if (proc == kPythia8) {
   Printf("<<<<<<<<<<<<<<<<<<< Processing Mb. Pythia 8 >>>>>>>>>>>>>>>>>>>>");
   gener = MyPythia8();
 }

 else if (proc == kPythia6D6T) {
   Printf("<<<<<<<<<<<<<<<<<<< Processing Mb. PYTHIA D6T >>>>>>>>>>>>>>>>>>>>");
   gener = MbPythiaTuneD6T();
 }

 else if (proc == kPythiaPerugia0) {
   Printf("<<<<<<<<<<<<<<<<<<< Processing Mb. PYTHIA Perugia0 >>>>>>>>>>>>>>>>>>>>");
   gener = MbPythiaTunePerugia0();
 }

 else if (proc == kPythia6ATLAS) {
   Printf("<<<<<<<<<<<<<<<<<<< Processing Mb. PYTHIA ATLAS>>>>>>>>>>>>>>>>>>>>");
   gener = MbPythiaTuneATLAS();
 }

 else if (proc == kPythia6ATLAS_Flat) {
   Printf("<<<<<<<<<<<<<<<<<<< Processing Mb. PYTHIA ATLAS_FLAT >>>>>>>>>>>>>>>>>>>>");
   gener = MbPythiaTuneATLAS_Flat();
 }

 else if (proc == kPythiaJets ) {
   Printf("<<<<<<<<<<<<<<<<<<< Processing Pythia Jets >>>>>>>>>>>>>>>>>>>>");
   gener = PythiaJets();
 }

 else if (proc == kPhojet) {
   Printf("<<<<<<<<<<<<<<<<<<< Processing Mb. PHOJET >>>>>>>>>>>>>>>>>>>>");
   gener = MbPhojet();
 }

 else if (proc == kHijing) {
   Printf("<<<<<<<<<<<<<<<<<<< Processing Mb. HIJING >>>>>>>>>>>>>>>>>>>>");
   gener = Hijing();
 }

 else if (proc == kHijing2000) {
   Printf("<<<<<<<<<<<<<<<<<<< Processing Mb. HIJING 2000 >>>>>>>>>>>>>>>>>>>>");
   gener = Hijing2000();
 }

 else if (proc ==kHijing_pA) {
   Printf("<<<<<<<<<<<<<<<<<<< Processing Mb. pA Hijing >>>>>>>>>>>>>>>>>>>>");
   gener = Hijing_pA(kTRUE);
 }

 else if (proc == kDPMjet) {
   Printf("<<<<<<<<<<<<<<<<<<< Processing  DMPJet  >>>>>>>>>>>>>>>>>>>>");
   gener = DPMjet();
 }

 else if (proc == kDPMjet_pA) {
   Printf("<<<<<<<<<<<<<<<<<<< Processing  DMPJet pA >>>>>>>>>>>>>>>>>>>>");
   gener = DPMjet_pA(kFALSE);
 }

 else if (proc == kAmptDefault) {
   Printf("<<<<<<<<<<<<<<<<<<< Processing  AMPT Default >>>>>>>>>>>>>>>>>>>>");
   gener = AmptDefault();
 }

 else if (proc == kAmptFlowStringMelting) {
     Printf("<<<<<<<<<<<<<<<<<<< Processing  AMPT With Flow  >>>>>>>>>>>>>>>>>>>>");
     gener = AmptStringMelting();
 }

 else if (proc == kAmptStringMeltingNoART) {
     Printf("<<<<<<<<<<<<<<<<<<< Processing  AMPT With Flow  >>>>>>>>>>>>>>>>>>>>");
     gener = AmptStringMeltingNoART();
 }

 else if (proc == kAmptpA) {
   Printf("<<<<<<<<<<<<<<<<<<< Processing  AMPT pA  >>>>>>>>>>>>>>>>>>>>");
   gener = AmptpA();
 }


 else if (proc == kAmptReducedFlow) {
   // Specific Fastgen
   Printf("<<<<<<<<<<<<<<<<<<< Processing  AMPT With Reduced Flow >>>>>>>>>>>>>>>>>>>>");
   gener = AmptReducedFlow();
 }

 else {
   cout << "ERROR : Wrong Procss Selcted !!!" << endl;
   return;
 }

//ALL HISTOS SAME AS PYTHIALoop2.C

    TH1F *hEnergy = CreateEnergyHistogram("hEnergy");
    TH1F *hETAll = CreateEnergyHistogram("hETAll");
    TH1F *hETcharge = CreateEnergyHistogram("hETcharge");
    TH2F *hETChOverAll = CreateEnergyHistogram2("hETChOverAll");
    TH1F *hETpip = CreateEnergyHistogram("hETPiPlus");
    TH1F *hETpim = CreateEnergyHistogram("hETPiMinus");
    TH1F *hETpi0 = CreateEnergyHistogram("hETPi0");
    TH1F *hETKp = CreateEnergyHistogram("hETKPlus");
    TH1F *hETKm = CreateEnergyHistogram("hETKMinus");
    TH1F *hETKL = CreateEnergyHistogram("hETKL");
    TH1F *hETKS = CreateEnergyHistogram("hETKS");
    TH1F *hETEta = CreateEnergyHistogram("hETEta");
    TH1F *hETOmega = CreateEnergyHistogram("hETomega");
    TH1F *hETLambda0 = CreateEnergyHistogram("hETLambda0");
    TH1F *hETLambdaBar0 = CreateEnergyHistogram("hETLambdaBar0");
    TH1F *hETp = CreateEnergyHistogram("hETp");
    TH1F *hETn = CreateEnergyHistogram("hETn");
    TH1F *hETp_ = CreateEnergyHistogram("hETpBar");
    TH1F *hETn_ = CreateEnergyHistogram("hETnBar");

    TH1F *homegaVSpi0 = CreateEnergyHistogram("hETomegaOverpi0");
    TH1F *hEtaVSpi0 = CreateEnergyHistogram("hETEtaOverpi0");
    TH1F *hETPi0Check1 = CreateEnergyHistogram("hETPi0Check1");
    TH1F *hETPi0Check2 = CreateEnergyHistogram("hETPi0Check2");
    TH1F *hETK0LCheck = CreateEnergyHistogram("hETK0LCheck");
    TH1F *hETK0SCheck = CreateEnergyHistogram("hETK0SCheck");
    TH1F *hETKMCheck = CreateEnergyHistogram("hETKMCheck");
    TH1F *hETpCheck = CreateEnergyHistogram("hETpCheck");
    TH1F *hETpbarCheck = CreateEnergyHistogram("hETpbarCheck");


    TH1F *hPTAll = CreatePTHistogram("hPTAll");
    TH1F *hPTcharge = CreatePTHistogram("hPTcharge");
    TH1F *hPTpip = CreatePTHistogram("hptPiPlus");
    TH1F *hPTpim = CreatePTHistogram("hptPiMinus");
    TH1F *hPTpi0 = CreatePTHistogram("hptPi0");
    TH1F *hPTKp = CreatePTHistogram("hptKPlus");
    TH1F *hPTKm = CreatePTHistogram("hptKMinus");
    TH1F *hPTK0 = CreatePTHistogram("hptK0");
    TH1F *hPTKL = CreatePTHistogram("hptKL");
    TH1F *hPTKS = CreatePTHistogram("hptKS");
    TH1F *hPTEta = CreatePTHistogram("hptEta");
    TH1F *hPTOmega = CreatePTHistogram("hptomega");
    TH1F *hPTLambda0 = CreatePTHistogram("hptLambda0");
    TH1F *hPTLambdaBar0 = CreatePTHistogram("hptLambdaBar0");
    TH1F *hPTp = CreatePTHistogram("hptp");
    TH1F *hPTn = CreatePTHistogram("hptn");
    TH1F *hPTp_ = CreatePTHistogram("hptpBar");
    TH1F *hPTn_ = CreatePTHistogram("hptnBar");

//Assumption checks
    TH2F *hpi0check= CreatePTHistogram2("hpi0check"); // pi0=1/2(pip+pim)
    TH2F *hpncheck= CreatePTHistogram2("hpncheck"); // p=n
    TH2F *hpbarnbarcheck= CreatePTHistogram2("hpbarnbarcheck"); // pbar=nbar
    TH2F *hkcheck1= CreatePTHistogram2("hkcheck1"); // k+=k0s
    TH2F *hkcheck2= CreatePTHistogram2("hkcheck2"); // k+=k0L
    TH2F *hkcheck3= CreatePTHistogram2("hkcheck3"); // k-=k0s
    TH2F *hkcheck4= CreatePTHistogram2("hkcheck4"); // k-=k0L

//Assumption Feed
    TH2F *hpi0etpt = CreatePTHistogram2("hpi0etpt");
    TH2F *hpipetpt = CreatePTHistogram2("hpipetpt");
    TH2F *hpimetpt = CreatePTHistogram2("hpimetpt");

    TH2F *hpetpt = CreatePTHistogram2("hpetpt");
    TH2F *hnetpt = CreatePTHistogram2("hnetpt");

    TH2F *hpBARetpt = CreatePTHistogram2("hpBARetpt");
    TH2F *hnBARetpt = CreatePTHistogram2("hnBARetpt");

//NOTE: VERY IMPORTANT
    TH1F *hNEvents = new TH1F("hNEvents","Number of events",1,0,1.0);
    hNEvents->GetYaxis()->SetTitle("N_{events}");
    hNEvents->GetXaxis()->SetTitle("no title");


 AliPDG::AddParticlesToPdgDataBase();
 TDatabasePDG::Instance();
 Int_t SNN = gener->GetEnergyCMS();
 const char* filename = Form("AMPT_%iGeV.root",SNN);
 AliRunLoader* rl = AliRunLoader::Open(filename,"FASTRUN","recreate");
 char* filenameA = "GARBAGE.root";
 TFile* file =TFile::Open(filenameA, "RECREATE");
 if (!file || !file->IsOpen()){
   cout<<"Fatal Error:\n Unable to open file";
   return;
 }
 rl->SetCompressionLevel(2);
 rl->SetNumberOfEventsPerFile(nev);
 rl->LoadKinematics("RECREATE");
 rl->MakeTree("E");
 gAlice->SetRunLoader(rl);
 rl->MakeStack();
 AliStack* stack = rl->Stack();

 AliHeader* header = rl->GetHeader();
 /*
   Float_t sigmaz  = 5.4 / TMath::Sqrt(2.); // [cm]
   Float_t betast  = 3.5;                      // beta* [m]
   Float_t eps     = 3.75e-6;                   // emittance [m]
   Float_t gamma   = energy / 2.0 / 0.938272;  // relativistic gamma [1]
   Float_t sigmaxy = TMath::Sqrt(eps * betast / gamma) / TMath::Sqrt(2.) * 100.;  // [cm]
   printf("\n \n Diamond size x-y: %10.3e z: %10.3e\n \n", sigmaxy, sigmaz);
   gener->SetSigma(sigmaxy, sigmaxy, sigmaz);      // Sigma in (X,Y,Z) (cm) on IP position
   gener->SetVertexSmear(kPerEvent);
 */

 Int_t cpip=0;
 Int_t cpim=0;
 Int_t cpi0=0;
 Int_t cKp=0;
 Int_t cKm=0;
 Int_t cKL=0;
 Int_t cKS=0;
 Int_t cEta=0;
 Int_t cOmega=0; //THIS IS omega not Omega
 Int_t cLambda0=0;
 Int_t cLambdaBar0=0;
 Int_t cp=0;
 Int_t cn=0;
 Int_t cp_=0;
 Int_t cn_=0;

 gener->Init();
 gener->SetStack(stack);

 rl->CdGAFile();

 TStopwatch timer;
 timer.Start();
 for (Int_t iev = 0; iev < nev; iev++) {
   cout <<"============================================= Event number "<< iev << endl;
   //  Initialize event

   Double_t ETAll=0;
   Double_t ETcharge=0;
   Double_t ETpip=0;
   Double_t ETpim=0;
   Double_t ETpi0=0;
   Double_t ETKp=0;
   Double_t ETKm=0;
   Double_t ETK0=0;
   Double_t ETKL=0;
   Double_t ETKS=0;
   Double_t ETEta=0;
   Double_t ETOmega=0; //THIS IS omega not Omega
   Double_t ETLambda0=0;
   Double_t ETLambdaBar0=0;
   Double_t ETp=0;
   Double_t ETn=0;
   Double_t ETp_=0;
   Double_t ETn_=0;

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
   Double_t pETLambda0=0;
   Double_t pETLambdaBar0=0;
   Double_t pETp=0;
   Double_t pETn=0;
   Double_t pETp_=0;
   Double_t pETn_=0;

   Double_t PTAll=0;
   Double_t PTcharge=0;
   Double_t PTpip=0;
   Double_t PTpim=0;
   Double_t PTpi0=0;
   Double_t PTKp=0;
   Double_t PTKm=0;
   Double_t PTK0=0;
   Double_t PTKL=0;
   Double_t PTKS=0;
   Double_t PTEta=0;
   Double_t PTOmega=0; //THIS IS omega not Omega
   Double_t PTLambda0=0;
   Double_t PTLambdaBar0=0;
   Double_t PTp=0;
   Double_t PTn=0;
   Double_t PTp_=0;
   Double_t PTn_=0;

   //ETpart/ETALL
   Double_t pPTpip=0;
   Double_t pPTpim=0;
   Double_t pPTpi0=0;
   Double_t pPTKp=0;
   Double_t pPTKm=0;
   Double_t pPTKL=0;
   Double_t pPTKS=0;
   Double_t pPTEta=0;
   Double_t pPTOmega=0; //THIS IS omega not Omega
   Double_t pPTLambda0=0;
   Double_t pPTLambdaBar0=0;
   Double_t pPTp=0;
   Double_t pPTn=0;
   Double_t pPTp_=0;
   Double_t pPTn_=0;



   header->Reset(0,iev);
   rl->SetEventNumber(iev);
   stack->Reset();
   //rl->MakeTree("K");
   Int_t nprim = 0;
   Int_t ntrial = 0;
   //  Int_t ndstar = 0;
   stack->Reset();
   //stack->ConnectTree(rl->TreeK());
   cout<<"start\n";
   gener->Generate();
   const TObjArray* ALLPART = stack->Particles();
   ntrial++;
   nprim = stack->GetNprimary(); //use TParticle here
/*******************************************************************************
THIS IS THE PART THAT ACTUALLY DOES A THING THAT IS USEFUL
again thanks to andy.

Need to find and eliminate the trees! this will save lots of trouble with memory
therefor potential speed increase
*******************************************************************************/
   cout<<"Number of particles: "<<nprim<<endl;
   cout<<"Number of events: "<<nev<<endl;
   for (int partnum=0;partnum<=(nprim-1);partnum++){
     TParticle *part = (TParticle*) ALLPART->At(partnum);
     //if (partnum % 100 == 0) cout<<"Particle "<<partnum<<endl;
     if (!part){
       cout<<"fail\n";
       continue;
     }

     Int_t PDG=part->GetPdgCode();
     const char* name=part->GetName();
     //if (partnum % 100 == 0) {cout<<"particle no. "<<partnum<<" is "<<PDG<<" or "<<name<<endl;}
     Double_t pT=part->Pt();
     Double_t Pz=part->Pz();
     Double_t E=part->Energy();
     Double_t M=part->GetMass();
     Double_t Theta=atan(pT/Pz);
     Double_t Eta=part->Eta();
     Double_t rapidity=part->Y();
     Double_t P=pT/sin(Theta);
     Double_t Ei=sqrt(pow(P,2)+pow(M,2));
     Int_t mpartD=part->GetFirstDaughter();
     Int_t mpartD2=part->GetLastDaughter();
     Int_t mpP = part->GetFirstMother();
     Double_t mpL = 0; //TODO: get Time somehow!!
     Double_t partE;
     Int_t nobby=15;

     if (rapidity<.1){
/******************************************************************************
To account for  pT properly, E_T is dependant on type of particle.
For baryons and antibaryons, note that rest mass is included (see below)
******************************************************************************/
       if ((PDG==3122)||(PDG==2212)||(PDG==2112)){
         partE = (Ei-M) * sin(Theta);
       }
       else if ((PDG==-3122)||(PDG==-2212)||(PDG==-2112)){
         partE = (Ei+M) * sin(Theta);
       }
       else {
         partE = Ei * sin(Theta);
       }
       if (PDG==211){ //pi+ all included
         ETpip+=partE;
         ETAll+=partE;
         ETcharge+=partE;
         PTpip+=pT;
         PTAll+=pT;
         PTcharge=+pT;
         cpip++;
       }
       if (PDG==-211){ //pi- all included
         ETpim+=partE;
         ETAll+=partE;
         ETcharge+=partE;
         PTpim+=pT;
         PTAll+=pT;
         PTcharge=+pT;
         cpim++;
       }
       if (PDG==111){ //pi0 // BUG double counts for each omega counted!!!!!
         ETpi0+=partE;
         ETAll+=partE;
         PTpi0+=pT;
         PTAll+=pT;
         cpi0++;
       }
       if (PDG==321){ //K+ //all included
         ETKp+=partE;
         ETAll+=partE;
         ETcharge+=partE;
         PTKp+=pT;
         PTAll+=pT;
         PTcharge=+pT;
         cKp++;
       }
       if (PDG==-321){ //K- //all included
         ETKm+=partE;
         ETAll+=partE;
         ETcharge+=partE;
         PTKm+=pT;
         PTAll+=pT;
         PTcharge=+pT;
         cKm++;
       }
       if (PDG==130){// KL
         if ((mpL>=100)||(mpartD==0)){//counted if final particle or lifetime > 1cm
           ETKL+=partE;
           ETAll+=partE;
           PTKL+=pT;
           PTAll+=pT;
           cKL++;
       }}
       if (PDG==310){ //KS
         if ((mpL>=100)||(mpartD==0)){//counted if final particle or lifetime > 1cm
           ETKS+=partE;
           ETAll+=partE;
           PTKS+=pT;
           PTAll+=pT;
           cKS++;
       }}
       if (PDG==221){ //eta
         //cout<<ind<<"= particle number (check eta) running func:"<<endl;
         //nobby=GetDaughterCheck(pylis,ind,221,22);
         //cout<<nobby<<endl;
         if (nobby==1){
           ETEta+=partE;
           ETAll+=partE;
           PTEta+=pT;
           PTAll+=pT;
           cEta++;
       }}
       if (PDG==223){ // omega
         //cout<<ind<<"= particle number (check omega) running func:"<<endl;
         //nobby=GetDaughterCheck(pylis,ind,223,22);
         //cout<<nobby<<endl;
         if (nobby==1){
           ETOmega+=partE;
           ETAll+=partE;
           PTOmega+=pT;
           PTAll+=pT;
           cOmega++;
       }}
       if (PDG==3122){ //Lambda0
         if ((mpL>=100)||(mpartD==0)){
           ETLambda0+=partE;
           ETAll+=partE;
           PTLambda0+=pT;
           PTAll+=pT;
           cLambda0++;
       }}
       if (PDG==-3122){ //LambdaBar0
         if ((mpL>=100)||(mpartD==0)){
           ETLambdaBar0+=partE;
           ETAll+=partE;
           PTLambdaBar0+=pT;
           PTAll+=pT;
           cLambdaBar0++;
       }}
       if (PDG==2212){ //proton
         if (mpP!=0){
           ETp+=partE;
           ETAll+=partE;
           ETcharge+=partE;
           PTp+=pT;
           PTAll+=pT;
           PTcharge=+pT;
           cp++;
       }}
       if (PDG==2112){ //neutron
         ETn+=partE;
         ETAll+=partE;
         PTn+=pT;
         PTAll+=pT;
         cn++;
       }
       if (PDG==-2212){ //antiproton
         ETp_+=partE;
         ETAll+=partE;
         ETcharge+=partE;
         PTp_+=pT;
         PTAll+=pT;
         PTcharge=+pT;
         cp_++;
       }
       if (PDG==-2112){ //antineutron
         ETn_+=partE;
         ETAll+=partE;
         PTn_+=pT;
         PTAll+=pT;
         cn_++;
       }

     }
      }// end particle
   Double_t omegaET=0;
   Double_t pi0ET=0;
   Double_t ChAllRat=0;
   hETAll->Fill(ETAll);
   hETcharge->Fill(ETcharge);
   ChAllRat=ETcharge/ETAll;
   hETChOverAll->Fill(ETcharge,ChAllRat);
   //cout<<"TOTAL: "<<ETAll<<endl;
   if (ETpip!=0){
     pETpip=ETpip/ETAll;
     hETpip->Fill(ETpip);
     hPTpip->Fill(PTpip);
     hpipetpt->Fill(PTpip,ETpip); //for test purposes
     //cout<<"pip "<<ETpip<<endl;
   }
   if (ETpim!=0){
     pETpim=ETpim/ETAll;
     hETpim->Fill(ETpim);
     hPTpim->Fill(PTpim);
     hpimetpt->Fill(PTpim,ETpim); //for test purposes
     //cout<<"pim "<<ETpim<<endl;
   }
   if (ETpi0!=0){
     pETpi0=ETpi0/ETAll;
     hETpi0->Fill(ETpi0);
     hPTpi0->Fill(PTpi0);
     hpi0etpt->Fill(PTpi0,ETpi0); //for test purposes
     //cout<<"pi0 "<<ETpi0<<endl;
   }
   if (ETKp!=0){
     pETKp=ETKp/ETAll;
     hETKp->Fill(ETKp);
     hPTKp->Fill(PTKp);
     //cout<<"Kp "<<ETKp<<endl;
   }
   if (ETKm!=0){
     pETKm=ETKm/ETAll;
     hETKm->Fill(ETKm);
     hPTKm->Fill(PTKm);
     //cout<<"Kp "<<ETKp<<endl;
   }
   if (ETKL!=0){
     pETKL=ETKL/ETAll;
     hETKL->Fill(ETKL);
     hPTKL->Fill(PTKL);
     //cout<<"KL "<<ETKL<<endl;
   }
   if (ETKS!=0){
     pETKS=ETKS/ETAll;
     hETKS->Fill(ETKS);
     hPTKS->Fill(PTKS);
     //cout<<"KS "<<ETKS<<endl;
   }
   if (ETEta!=0){
     pETEta=ETEta/ETAll;
     hETEta->Fill(ETEta);
     hPTEta->Fill(PTEta);
     //cout<<"Eta "<<ETEta<<endl;
   }
   if (ETOmega!=0){
     pETOmega=ETOmega/ETAll;
     hETOmega->Fill(ETOmega);
     hPTOmega->Fill(PTOmega);
     //cout<<"Omega "<<ETOmega<<endl;
   }
   if (ETLambda0!=0){
     pETLambda0=ETLambda0/ETAll;
     hETLambda0->Fill(ETLambda0);
     hPTLambda0->Fill(PTLambda0);
     //cout<<"Lambda0 "<<ETLambda0<<endl;
   }
   if (ETLambdaBar0!=0){
     pETLambdaBar0=ETLambdaBar0/ETAll;
     hETLambdaBar0->Fill(ETLambdaBar0);
     hPTLambdaBar0->Fill(PTLambdaBar0);
     //cout<<"Lambda0 "<<ETLambda0<<endl;
   }
   if (ETp!=0){
     pETp=ETp/ETAll;
     hETp->Fill(ETp);
     hPTp->Fill(PTp);
     hpetpt->Fill(PTp,ETp);
     //cout<<"p "<<ETp<<endl;
   }
   if (ETn!=0){
     pETn=ETn/ETAll;
     hETn->Fill(ETn);
     hPTn->Fill(PTn);
     hnetpt->Fill(PTn,ETn);
     //cout<<"n "<<ETn<<endl;
   }
   if (ETp_!=0){
     pETp_=ETp_/ETAll;
     hETp_->Fill(ETp_);
     hPTp_->Fill(PTp_);
     hpBARetpt->Fill(PTp_,ETp_);
     //cout<<"p "<<ETp<<endl;
   }
   if (ETn_!=0){
     pETn_=ETn_/ETAll;
     hETn_->Fill(ETn_);
     hPTn_->Fill(PTn_);
     hnBARetpt->Fill(PTn_,ETn_);
     //cout<<"n "<<ETn<<endl;
   }
   stack->Reset();
   cout << "Number of particles " << nprim << endl;
   cout << "Number of trials " << ntrial << endl;
   //header->SetNprimary(stack->GetNprimary());
   //header->SetNtrack(stack->GetNtrack());
   stack->FinishEvent();
    //header->SetStack(stack);
    //rl->TreeE()->Fill();
    //rl->WriteKinematics("OVERWRITE");
 } // event loop
 timer.Stop();
 timer.Print();
 gener->FinishRun();
 rl->WriteHeader("OVERWRITE");
 gener->Write();
 rl->Write();
 rl->Clear();
 /**
 rl->Close(); //maybe
 rl->Open(); //maybe
 **/

 hETPi0Check2->Add(hETpip,hETpim,.5,.5);
 hETPi0Check1->Divide(hETPi0Check2,hETpi0);
 hETPi0Check1->Write();
 hETK0LCheck->Divide(hETKp,hETKL);
 hETK0LCheck->Write();
 hETK0SCheck->Divide(hETKp,hETKS);
 hETK0SCheck->Write();
 hETKMCheck->Divide(hETKp,hETKm);
 hETKMCheck->Write();
 hETpCheck->Divide(hETp,hETn);
 hETpCheck->Write();
 hETpbarCheck->Divide(hETp_,hETn_);
 hETpbarCheck->Write();
 hETAll->Write();
 hETChOverAll->Write();
 hETpip->Write();
 hETpim->Write();
 hETpi0->Write();
 hETKp->Write();
 hETKm->Write();
 hETKL->Write();
 hETKS->Write();
 hETEta->Write();
 hETOmega->Write();
 hETLambda0->Write();
 hETLambdaBar0->Write();
 hETp->Write();
 hETn->Write();
 hETp_->Write();
 hETn_->Write();
 hPTAll->Write();
 hPTpip->Write();
 hPTpim->Write();
 hPTpi0->Write();
 hPTKp->Write();
 hPTKm->Write();
 hPTKL->Write();
 hPTKS->Write();
 hPTEta->Write();
 hPTOmega->Write();
 hPTLambda0->Write();
 hPTLambdaBar0->Write();
 hPTp->Write();
 hPTn->Write();
 hPTp_->Write();
 hPTn_->Write();
 file->Write();
 file->Close();
}

//___________________________________________________//
void ProcessEnvironmentVars() {
    // Run type
    if (gSystem->Getenv("CONFIG_RUN_TYPE")) {
      for (Int_t iRun = 0; iRun < kRunMax; iRun++) {
        if (strcmp(gSystem->Getenv("CONFIG_RUN_TYPE"), pprRunName[iRun]) == 0) {
          proc = (PDC06Proc_t)iRun;
          cout<<"Run type set to "<<pprRunName[iRun]<<endl;
        }
      }
    }


    // Energy
    if (gSystem->Getenv("CONFIG_ENERGY")) {
      energy = atoi(gSystem->Getenv("CONFIG_ENERGY"));
      cout<<"Energy set to "<<energy<<" GeV"<<endl;
    }

    // Random Number seed
    if (gSystem->Getenv("CONFIG_SEED")) {
      seed = atoi(gSystem->Getenv("CONFIG_SEED"));
    }

    // Impact param
    if (gSystem->Getenv("CONFIG_BMIN")) {
      bMin = atof(gSystem->Getenv("CONFIG_BMIN"));
    }

    if (gSystem->Getenv("CONFIG_BMAX")) {
      bMax = atof(gSystem->Getenv("CONFIG_BMAX"));
    }
    cout<<"Impact parameter in ["<<bMin<<","<<bMax<<"]"<<endl;
}



//______________________________________________________________________
AliGenerator* MbPythia() // Mb Pythia
{
      comment = comment.Append(" pp: Pythia low-pt");
      AliGenPythia* pythia = new AliGenPythia(-1);
      /* pythia->SetMomentumRange(0, 999999.);
      pythia->SetThetaRange(0., 180.);
      pythia->SetYRange(-12.,12.);
      pythia->SetPtRange(0,1000.);*/
      pythia->SetProcess(kPyMb);
      pythia->SetEnergyCMS(energy);

      return pythia;
}


//______________________________________________________________________
AliGenerator* MyPythia8()
{
  AliGenPythiaPlus* gener = new AliGenPythiaPlus(AliPythia8::Instance());
  gener->SetProcess(kPyMbDefault);
  gener->SetEnergyCMS(energy);
  gener->SetEventListRange(-1, 2);
  return gener;
}



//______________________________________________________________________
AliGenerator* MbPythiaTuneD6T()
{
      comment = comment.Append(" pp: Pythia low-pt");
      AliGenPythia* pythia = new AliGenPythia(-1);
      pythia->SetMomentumRange(0, 999999.);
      pythia->SetThetaRange(0., 180.);
      pythia->SetYRange(-12.,12.);
      pythia->SetPtRange(0,1000.);
      pythia->SetProcess(kPyMb);
      pythia->SetEnergyCMS(energy);
//    Tune
//    109     D6T : Rick Field's CDF Tune D6T (NB: needs CTEQ6L pdfs externally)
      pythia->SetTune(109); // F I X
      pythia->SetStrucFunc(kCTEQ6l);
//
      return pythia;
}

//______________________________________________________________________
AliGenerator* MbPythiaTunePerugia0()
{
      comment = comment.Append(" pp: Pythia low-pt (Perugia0)");
//
//    Pythia
      AliGenPythia* pythia = new AliGenPythia(-1);
      /* pythia->SetMomentumRange(0, 999999.);
      pythia->SetThetaRange(0., 180.);
      pythia->SetYRange(-12.,12.);
      pythia->SetPtRange(0,1000.);*/
      pythia->SetProcess(kPyMb);
      pythia->SetEnergyCMS(energy);
//    Tune
//    320     Perugia 0
      pythia->SetTune(320);
      pythia->UseNewMultipleInteractionsScenario();
//
      return pythia;
}

//______________________________________________________________________
AliGenerator* MbPythiaTuneATLAS()
{
      comment = comment.Append(" pp: Pythia low-pt");
//
//    Pythia
      AliGenPythia* pythia = new AliGenPythia(-1);
      /*   pythia->SetMomentumRange(0, 999999.);
      pythia->SetThetaRange(0., 180.);
      pythia->SetYRange(-12.,12.);
      pythia->SetPtRange(0,1000.);*/
      pythia->SetProcess(kPyMb);
      pythia->SetEnergyCMS(energy);
//    Tune
//    C   306 ATLAS-CSC: Arthur Moraes' (new) ATLAS tune (needs CTEQ6L externally)
      pythia->SetTune(306);
      pythia->SetStrucFunc(kCTEQ6l);
//
      return pythia;
}


//______________________________________________________________________
AliGenerator* MbPythiaTuneATLAS_Flat()
{
  AliGenPythia* pythia = (AliGenPythia*)MbPythiaTuneATLAS();

  comment = comment.Append("; flat multiplicity distribution");

  // set high multiplicity trigger
  // this weight achieves a flat multiplicity distribution
  TH1 *weight = new TH1D("weight","weight",201,-0.5,200.5);
  weight->SetBinContent(1,5.49443);
  weight->SetBinContent(2,8.770816);
  weight->SetBinContent(6,0.4568624);
  weight->SetBinContent(7,0.2919915);
  weight->SetBinContent(8,0.6674189);
  weight->SetBinContent(9,0.364737);
  weight->SetBinContent(10,0.8818444);
  weight->SetBinContent(11,0.531885);
  weight->SetBinContent(12,1.035197);
  weight->SetBinContent(13,0.9394057);
  weight->SetBinContent(14,0.9643193);
  weight->SetBinContent(15,0.94543);
  weight->SetBinContent(16,0.9426507);
  weight->SetBinContent(17,0.9423649);
  weight->SetBinContent(18,0.789456);
  weight->SetBinContent(19,1.149026);
  weight->SetBinContent(20,1.100491);
  weight->SetBinContent(21,0.6350525);
  weight->SetBinContent(22,1.351941);
  weight->SetBinContent(23,0.03233504);
  weight->SetBinContent(24,0.9574557);
  weight->SetBinContent(25,0.868133);
  weight->SetBinContent(26,1.030998);
  weight->SetBinContent(27,1.08897);
  weight->SetBinContent(28,1.251382);
  weight->SetBinContent(29,0.1391099);
  weight->SetBinContent(30,1.192876);
  weight->SetBinContent(31,0.448944);
  weight->SetBinContent(32,1);
  weight->SetBinContent(33,1);
  weight->SetBinContent(34,1);
  weight->SetBinContent(35,1);
  weight->SetBinContent(36,0.9999997);
  weight->SetBinContent(37,0.9999997);
  weight->SetBinContent(38,0.9999996);
  weight->SetBinContent(39,0.9999996);
  weight->SetBinContent(40,0.9999995);
  weight->SetBinContent(41,0.9999993);
  weight->SetBinContent(42,1);
  weight->SetBinContent(43,1);
  weight->SetBinContent(44,1);
  weight->SetBinContent(45,1);
  weight->SetBinContent(46,1);
  weight->SetBinContent(47,0.9999999);
  weight->SetBinContent(48,0.9999998);
  weight->SetBinContent(49,0.9999998);
  weight->SetBinContent(50,0.9999999);
  weight->SetBinContent(51,0.9999999);
  weight->SetBinContent(52,0.9999999);
  weight->SetBinContent(53,0.9999999);
  weight->SetBinContent(54,0.9999998);
  weight->SetBinContent(55,0.9999998);
  weight->SetBinContent(56,0.9999998);
  weight->SetBinContent(57,0.9999997);
  weight->SetBinContent(58,0.9999996);
  weight->SetBinContent(59,0.9999995);
  weight->SetBinContent(60,1);
  weight->SetBinContent(61,1);
  weight->SetBinContent(62,1);
  weight->SetBinContent(63,1);
  weight->SetBinContent(64,1);
  weight->SetBinContent(65,0.9999999);
  weight->SetBinContent(66,0.9999998);
  weight->SetBinContent(67,0.9999998);
  weight->SetBinContent(68,0.9999999);
  weight->SetBinContent(69,1);
  weight->SetBinContent(70,1);
  weight->SetBinContent(71,0.9999997);
  weight->SetBinContent(72,0.9999995);
  weight->SetBinContent(73,0.9999994);
  weight->SetBinContent(74,1);
  weight->SetBinContent(75,1);
  weight->SetBinContent(76,1);
  weight->SetBinContent(77,1);
  weight->SetBinContent(78,0.9999999);
  weight->SetBinContent(79,1);
  weight->SetBinContent(80,1);
  weight->SetEntries(526);
  Int_t limit = weight->GetRandom();
  pythia->SetTriggerChargedMultiplicity(limit, 1.4);
  comment = comment.Append(Form("; multiplicity threshold set to %d in |eta| < 1.4", limit));
  return pythia;
}

//______________________________________________________________________
AliGenerator* PythiaJets()
{
      comment = comment.Append(" pp: Pythia low-pt");
//
//    Pythia
      AliGenPythia* pythia = new AliGenPythia(-1);
      /*   pythia->SetMomentumRange(0, 999999.);
      pythia->SetThetaRange(0., 180.);
      pythia->SetYRange(-12.,12.);
      pythia->SetPtRange(0,1000.);*/
      pythia->SetProcess(kPyJets);
      pythia->SetEnergyCMS(energy);
      pythia->SetStrucFunc(kCTEQ6l);
      //  pythia->SetPtHard(50, 1000);
//
      return pythia;
}



//___________________________________________________________________
AliGenerator* MbPhojet()
{
    comment = comment.Append(" pp: Pythia low-pt");
    //
    //    Pythia
    AliGenPythia* dpmjet = new AliGenPythia(-1);
    /*   pythia->SetMomentumRange(0, 999999.);
     pythia->SetThetaRange(0., 180.);
     pythia->SetYRange(-12.,12.);
     pythia->SetPtRange(0,1000.);*/
    dpmjet->SetProcess(kPyJets);
    dpmjet->SetEnergyCMS(energy);
    dpmjet->SetStrucFunc(kCTEQ6l);
    //  pythia->SetPtHard(50, 1000);
    //
    return dpmjet;
}


//__________________________________________________________________
AliGenerator* Hijing()
{
    AliGenHijing *gener = new AliGenHijing(-1);
    // centre of mass energy
    gener->SetEnergyCMS(energy);
    gener->SetImpactParameterRange(bMin, bMax);
    // reference frame
    gener->SetReferenceFrame("CMS");
    // projectile
    gener->SetProjectile("A", 208, 82);
    gener->SetTarget    ("A", 208, 82);
    // tell hijing to keep the full parent child chain
    gener->KeepFullEvent();
    // enable jet quenching
    gener->SetJetQuenching(1);
    // enable shadowing
    gener->SetDecaysOff(0);
    gener->SetShadowing(1);
    // Don't track spectators
    gener->SetSpectators(0);
    // kinematic selection
    gener->SetSelectAll(0);
    return gener;
}


//__________________________________________________________________
AliGenerator* Hijing2000()
{
  AliGenHijing *gener = (AliGenHijing*) Hijing();
  gener->SetJetQuenching(0);
  gener->SetPtHardMin (2.3);
  return gener;
}



//_____________________________________________________________
AliGenerator* Hijing_pA(Bool_t kSlowN) {
  AliGenCocktail *gener = 0x0;
  if (kSlowN) {
    gener  = new AliGenCocktail();
    gener->SetProjectile("A", 208, 82);
    gener->SetTarget    ("P",   1,  1);
    gener->SetEnergyCMS(energy);
  }

  AliGenHijing   *hijing = new AliGenHijing(-1);
  // centre of mass energy
  hijing->SetEnergyCMS(energy);
  // impact parameter range
  hijing->SetImpactParameterRange(0., 100.);
  // reference frame
  hijing->SetReferenceFrame("CMS");
  hijing->SetBoostLHC(1);
  // projectile
  hijing->SetProjectile("A", 208, 82);
  hijing->SetTarget    ("P", 1, 1);
  // tell hijing to keep the full parent child chain
  hijing->KeepFullEvent();
  // enable jet quenching
  hijing->SetJetQuenching(0);
  // enable shadowing
  hijing->SetShadowing(1);
  // kinematic selection
  hijing->SetSelectAll(0);

  if (!kSlowN) {
    // DO track spectators
    hijing->SetSpectators(1);
    return hijing;
  }
  else {
    // Cocktail with slow nucleons generator
    // DO NOT track spectators
    hijing->SetSpectators(0);
    AliGenSlowNucleons* gray   = new AliGenSlowNucleons(1);
    AliCollisionGeometry* coll = hijing->CollisionGeometry();
    AliSlowNucleonModelExp* model = new AliSlowNucleonModelExp();
    //  Not yet in the release...
    //	  model->SetSaturation(kTRUE);
    gray->SetSlowNucleonModel(model);
    gray->SetTarget(208,82);
    gray->SetThetaDist(1);
    gray->SetProtonDirection(1);
    //	  gray->SetDebug(1);
    gray->SetNominalCmsEnergy(2.*pBeamEnergy);
    gray->NeedsCollisionGeometry();
    gray->SetCollisionGeometry(coll);

    gener->AddGenerator(hijing, "Hijing pPb", 1);
    gener->AddGenerator(gray, "Gray Particles", 1);

    return gener;
  }
}



//__________________________________________________________________
AliGenerator* DPMjet()
{
    comment = comment.Append(" pp: Pythia low-pt");
    //
    //    Pythia
    AliGenPythia* dpmjet = new AliGenPythia(-1);
    /*   pythia->SetMomentumRange(0, 999999.);
     pythia->SetThetaRange(0., 180.);
     pythia->SetYRange(-12.,12.);
     pythia->SetPtRange(0,1000.);*/
    dpmjet->SetProcess(kPyJets);
    dpmjet->SetEnergyCMS(energy);
    dpmjet->SetStrucFunc(kCTEQ6l);
    //  pythia->SetPtHard(50, 1000);
    //
    return dpmjet;
}


//____________________________________________
AliGenerator* DPMjet_pA(Bool_t fragments)
{
    comment = comment.Append(" pp: Pythia low-pt");
    //
    //    Pythia
    AliGenPythia* dpmjet = new AliGenPythia(-1);
    /*   pythia->SetMomentumRange(0, 999999.);
     pythia->SetThetaRange(0., 180.);
     pythia->SetYRange(-12.,12.);
     pythia->SetPtRange(0,1000.);*/
    dpmjet->SetProcess(kPyJets);
    dpmjet->SetEnergyCMS(energy);
    dpmjet->SetStrucFunc(kCTEQ6l);
    //  pythia->SetPtHard(50, 1000);
    //
    return dpmjet;

}

//__________________________________________________________________
AliGenerator* AmptDefault()
{
  AliGenAmpt *genHi = new AliGenAmpt(-1);
  //=========================================================================
  // THE DECAYER
  AliDecayer *decayer = new AliDecayerPythia();
  cout << "*****************************************" << endl;
  genHi->SetForceDecay( kHadronicD );
  genHi->SetDecayer( decayer );
  //=========================================================================
  genHi->SetEnergyCMS(energy);
  genHi->SetReferenceFrame("CMS");
  genHi->SetProjectile("A",208,82);
  genHi->SetTarget("A",208,82);

  genHi->SetIsoft(1); //1: defaul - 4: string melting
  genHi->SetStringFrag(0.5,0.9); //Lund string frgmentation parameters
  genHi->SetPtHardMin (3);
  //genHi->SetImpactParameterRange(9.,9.5);
  genHi->SetImpactParameterRange(bMin,bMax);

  // Xmu = 3.2 fm^-1 and as = 0.33 ==> sigma_{partonic} = 1.5mb
  // Ntmax = 150
  // v_2{2} = 0.102105 +/- 0.000266894
  // v_2{4} = 0.0829477 +/- 0.00106158

  genHi->SetNtMax(150);        // NTMAX: number of timesteps (D=150)
  genHi->SetXmu(3.2264);        // parton screening mass in fm^(-1) (D=3.2264d0)

  genHi->SetJetQuenching(0);  // enable jet quenching
  genHi->SetShadowing(1);     // enable shadowing
  genHi->SetDecaysOff(1);     // neutral pion and heavy particle decays switched off
  genHi->SetSpectators(0);    // track spectators
  //Boost into LHC lab frame
  genHi->SetBoostLHC(1);
  //  genHi->Init();
  genHi->SetRandomReactionPlane(kTRUE);

  return genHi;
}

//__________________________________________________________________
AliGenerator* AmptStringMelting()
{
  AliGenAmpt *genHi = new AliGenAmpt(-1);
  //=========================================================================
  // THE DECAYER
  AliDecayer *decayer = new AliDecayerPythia();
  cout << "*****************************************" << endl;
  genHi->SetForceDecay( kHadronicD );
  genHi->SetDecayer( decayer );
  //=========================================================================
  genHi->SetEnergyCMS(energy);
  genHi->SetReferenceFrame("CMS");
  genHi->SetProjectile("A",208,82);
  genHi->SetTarget("A",208,82);

  genHi->SetIsoft(4); //1: defaul - 4: string melting
  genHi->SetStringFrag(0.5,0.9); //Lund string frgmentation parameters
  genHi->SetPtHardMin (3);
  //genHi->SetImpactParameterRange(9.,9.5);
  genHi->SetImpactParameterRange(bMin,bMax);

  // Xmu = 3.2 fm^-1 and as = 0.33 ==> sigma_{partonic} = 1.5mb
  // Ntmax = 150
  // v_2{2} = 0.102105 +/- 0.000266894
  // v_2{4} = 0.0829477 +/- 0.00106158

  genHi->SetNtMax(150);        // NTMAX: number of timesteps (D=150)
  genHi->SetXmu(3.2264);        // parton screening mass in fm^(-1) (D=3.2264d0)

  genHi->SetJetQuenching(0);  // enable jet quenching
  genHi->SetShadowing(1);     // enable shadowing
  genHi->SetDecaysOff(1);     // neutral pion and heavy particle decays switched off
  genHi->SetSpectators(0);    // track spectators
  //Boost into LHC lab frame
  genHi->SetBoostLHC(1);
//  genHi->Init();
  genHi->SetRandomReactionPlane(kTRUE);
  return genHi;

}

//__________________________________________________________________
AliGenerator* AmptStringMeltingNoART()
{
    AliGenAmpt *genHi = new AliGenAmpt(-1);
    //=========================================================================
    // THE DECAYER
    AliDecayer *decayer = new AliDecayerPythia();
    cout << "*****************************************" << endl;
    genHi->SetForceDecay( kHadronicD );
    genHi->SetDecayer( decayer );
    //=========================================================================
    genHi->SetEnergyCMS(energy);
    genHi->SetReferenceFrame("CMS");
    genHi->SetProjectile("A",208,82);
    genHi->SetTarget("A",208,82);

    genHi->SetIsoft(4); //1: defaul - 4: string melting
    genHi->SetStringFrag(0.5,0.9); //Lund string frgmentation parameters
    genHi->SetPtHardMin (3);
    //genHi->SetImpactParameterRange(9.,9.5);
    genHi->SetImpactParameterRange(bMin,bMax);

    // Xmu = 3.2 fm^-1 and as = 0.33 ==> sigma_{partonic} = 1.5mb
    // Ntmax = 150
    // v_2{2} = 0.102105 +/- 0.000266894
    // v_2{4} = 0.0829477 +/- 0.00106158

    genHi->SetNtMax(3);        // NTMAX: number of timesteps (D=150)
    genHi->SetXmu(3.2264);        // parton screening mass in fm^(-1) (D=3.2264d0)

    genHi->SetJetQuenching(0);  // enable jet quenching
    genHi->SetShadowing(1);     // enable shadowing
    genHi->SetDecaysOff(1);     // neutral pion and heavy particle decays switched off
    genHi->SetSpectators(0);    // track spectators
    //Boost into LHC lab frame
    genHi->SetBoostLHC(1);
    //  genHi->Init();
    genHi->SetRandomReactionPlane(kTRUE);
    return genHi;

}


//__________________________________________________________________
AliGenerator* AmptReducedFlow()
{
  AliGenAmpt *genHi = new AliGenAmpt(-1);
  //=========================================================================
  // THE DECAYER
  AliDecayer *decayer = new AliDecayerPythia();
  cout << "*****************************************" << endl;
  genHi->SetForceDecay( kHadronicD );
  genHi->SetDecayer( decayer );
  //=========================================================================
  genHi->SetEnergyCMS(energy);
  genHi->SetReferenceFrame("CMS");
  genHi->SetProjectile("A",208,82);
  genHi->SetTarget("A",208,82);

  genHi->SetIsoft(4); //1: defaul - 4: string melting
  genHi->SetStringFrag(0.5,0.9); //Lund string frgmentation parameters
  genHi->SetPtHardMin (3);
  //genHi->SetImpactParameterRange(9.,9.5);
  genHi->SetImpactParameterRange(bMin,bMax);

  // Xmu = 12.4 fm^-1 and as = 0.33 ==> sigma_{partonic} = 0.1mb
  // Ntmax = 20
  // flow estimates from Q-cumulants
  // (POI, without weights)
  // v_2{2} = 0.0549735 +/- 0.000270249
  // v_2{4} = 0.0421905 +/- 0.00189449

  genHi->SetNtMax(20);        // NTMAX: number of timesteps (D=150)
  genHi->SetXmu(12.4);        // parton screening mass in fm^(-1) (D=3.2264d0)

  genHi->SetJetQuenching(0);  // enable jet quenching
  genHi->SetShadowing(1);     // enable shadowing
  genHi->SetDecaysOff(1);     // neutral pion and heavy particle decays switched off
  genHi->SetSpectators(0);    // track spectators
  //Boost into LHC lab frame
  genHi->SetBoostLHC(1);
 // genHi->Init();
  genHi->SetRandomReactionPlane(kTRUE);
  return genHi;

}

//__________________________________________________________________
AliGenerator* AmptpA()
{
    AliGenAmpt *genHi = new AliGenAmpt(-1);
    //=========================================================================
    // THE DECAYER
    AliDecayer *decayer = new AliDecayerPythia();
    cout << "*****************************************" << endl;
    genHi->SetForceDecay( kHadronicD );
    genHi->SetDecayer( decayer );
    //=========================================================================
    genHi->SetEnergyCMS(energy);
    genHi->SetReferenceFrame("CMS");
    genHi->SetProjectile("A", 208, 82);
    genHi->SetTarget    ("P", 1, 1);
    genHi->SetIsoft(4); //1: defaul - 4: string melting
    genHi->SetStringFrag(0.5,0.9); //Lund string frgmentation parameters
    genHi->SetPtHardMin (3);
    //genHi->SetImpactParameterRange(9.,9.5);
    genHi->SetImpactParameterRange(bMin,bMax);
    genHi->SetNtMax(1500); //NTMAX: number of timesteps (D=150)
    genHi->SetXmu(3.2264); //parton screening mass in fm^(-1) (D=3.2264d0)
    genHi->SetJetQuenching(0); // enable jet quenching
    genHi->SetShadowing(1);    // enable shadowing
    genHi->SetDecaysOff(1);    // neutral pion and heavy particle decays switched off
    genHi->SetSpectators(0);   // track spectators
    //Boost into LHC lab frame
    genHi->SetBoostLHC(1);
    //  genHi->Init();
    genHi->SetRandomReactionPlane(kTRUE);
    return genHi;

}

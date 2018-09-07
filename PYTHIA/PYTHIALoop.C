//Based on SimplePythiaLoop.C

//See PythiaStartUp.C for more info
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
  TH1F *histo = new TH1F(name,name,1000,0,1.1);
  histo->GetYaxis()->SetTitle("number of entries");
  histo->GetXaxis()->SetTitle("Energy_Tpart/ETAll");
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
Int_t GetKFConversion(const Int_t kfc, const vector<KF_Code>& partname){  //This is a convenient conversion from KF code to an index without spaces
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
  if (n==-1){
    k=k*(-1);
  }
  return k;
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
  char *filename = Form("hh%iGeVoutfile%i.root",modSNN,jobID);
  ifstream in;
  in.open("KF_Code.dat");
  while (in.good()){
    partname.push_back(KF_Code());
    in>>partname[l].KFc>>partname[l].name;
    l++;
  }



  cout<<"I made it here line 144"<<endl;
  // Create an instance of the Pythia event generator ...
  TPythia6* pythia = new TPythia6();
  cout<<"I made it here line 157"<<endl;


  //NOTE: check that the tune is correct
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

//NOTE: set seed as jobID
  UInt_t seed = jobID;
  cout<<seed<<endl; //gives seed
  if( (seed>=0) && (seed<=900000000) ) {
    pythia->SetMRPY(1, seed);                   // set seed
    pythia->SetMRPY(2, 0);                      // use new seed
    pythia->SetPARP(2,7.0);                     // lowest collision energy THIS IS A MOD to allow (add more NOTE)
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

//NOTE:get rid of spaces and underscores and special chars
//NOTE: try sed with simga to fix first
    TH1F *hSPALL =CreateSpeciesHistogram("hSPALL",19804440,-9902220,9902220);
    TH3F *hPar = CreateParentHistogram("hPar");
    TH1F *hEnergy = CreateEnergyHistogram("hEnergy");
    TH1F *hETAll = CreateEnergyHistogram("hETAll");
    TH1F *hETpip = CreateEnergyHistogram("hETPiPlus");
    TH1F *hETpim = CreateEnergyHistogram("hETPiMinus");
    TH1F *hETpi0 = CreateEnergyHistogram("hETPi0");
    TH1F *hETKp = CreateEnergyHistogram("hETKPlus");
    TH1F *hETKm = CreateEnergyHistogram("hETKMinus");
    TH1F *hETKL = CreateEnergyHistogram("hETKL");
    TH1F *hETKS = CreateEnergyHistogram("hETKS");
    TH1F *hETEta = CreateEnergyHistogram("hETEta");
    TH1F *hETOmega = CreateEnergyHistogram("hETomega");
    TH1F *hETOmegaM = CreateEnergyHistogram("hETOmegaMinus");
    TH1F *hETLambda0 = CreateEnergyHistogram("hETLambda0");
    TH1F *hETLambdaBar0 = CreateEnergyHistogram("hETLambdaBar0");
    TH1F *hETSigmaP = CreateEnergyHistogram("hETSigmaPlus");
    TH1F *hETSigmaM = CreateEnergyHistogram("hETSimgaMinus");
    TH1F *hETSigma0 = CreateEnergyHistogram("hETSigma0");
    TH1F *hETXiM = CreateEnergyHistogram("hETXiMinus");
    TH1F *hETXi0 = CreateEnergyHistogram("hETXi0");
    TH1F *hETgamma = CreateEnergyHistogram("hETgamma");
    TH1F *hETp = CreateEnergyHistogram("hETp");
    TH1F *hETn = CreateEnergyHistogram("hETn");
    TH1F *hETp_ = CreateEnergyHistogram("hETpBar");
    TH1F *hETn_ = CreateEnergyHistogram("hETnBar");
    TH1F *hETmu = CreateEnergyHistogram("hETmu");
    TH1F *hETmu_ = CreateEnergyHistogram("hETmuBar");
    TH1F *hETe = CreateEnergyHistogram("hETe");
    TH1F *hETe_ = CreateEnergyHistogram("hETeBar");
    TH1F *homegaVSpi0 = CreateEnergyHistogram("hETomegaOverpi0");
    TH1F *hEtaVSpi0 = CreateEnergyHistogram("hETEtaOverpi0");


    TH1F *hPTAll = CreatePTHistogram("hETAll");
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
    TH1F *hPTOmegaM = CreatePTHistogram("hptOmegaMinus");
    TH1F *hPTLambda0 = CreatePTHistogram("hptLambda0");
    TH1F *hPTLambdaBar0 = CreatePTHistogram("hptLambdaBar0");
    TH1F *hPTSigmaP = CreatePTHistogram("hptSigmaPlus");
    TH1F *hPTSigmaM = CreatePTHistogram("hptSimgaMinus");
    TH1F *hPTSigma0 = CreatePTHistogram("hptSigma0");
    TH1F *hPTXiM = CreatePTHistogram("hptXiMinus");
    TH1F *hPTXi0 = CreatePTHistogram("hptXi0");
    TH1F *hPTgamma = CreatePTHistogram("hptgamma");
    TH1F *hPTp = CreatePTHistogram("hptp");
    TH1F *hPTn = CreatePTHistogram("hptn");
    TH1F *hPTp_ = CreatePTHistogram("hptpBar");
    TH1F *hPTn_ = CreatePTHistogram("hptnBar");
    TH1F *hPTmu = CreatePTHistogram("hptmu");
    TH1F *hPTmu_ = CreatePTHistogram("hptmuBar");
    TH1F *hPTe = CreatePTHistogram("hpte");
    TH1F *hPTe_ = CreatePTHistogram("hpteBar");
//NOTE: VERY IMPORTANT
    TH1F *hNEvents = new TH1F("hNEvents","Number of events",1,0,1.0);
    hNEvents->GetYaxis()->SetTitle("N_{events}");
    hNEvents->GetXaxis()->SetTitle("no title");



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
    Double_t ETKp=0;
    Double_t ETKm=0;
    Double_t ETK0=0;
    Double_t ETKL=0;
    Double_t ETKS=0;
    Double_t ETEta=0;
    Double_t ETOmega=0; //THIS IS omega not Omega
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


      Double_t PTAll=0;
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
      Double_t PTOmegaM=0;
      Double_t PTLambda0=0;
      Double_t PTLambdaBar0=0;
      Double_t PTRhoP=0;
      Double_t PTRho0=0;
      Double_t PTSigmaP=0;
      Double_t PTSigmaM=0;
      Double_t PTSigma0=0;
      Double_t PTDeltaPP=0;
      Double_t PTDeltaP=0;
      Double_t PTDeltaM=0;
      Double_t PTDelta0=0;
      Double_t PTXiM=0;
      Double_t PTXi0=0;
      Double_t PTgamma=0;
      Double_t PTp=0;
      Double_t PTn=0;
      Double_t PTp_=0;
      Double_t PTn_=0;
      Double_t PTmu=0;
      Double_t PTmu_=0;
      Double_t PTe_=0;
      Double_t PTe=0;

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
      Double_t pPTOmegaM=0;
      Double_t pPTLambda0=0;
      Double_t pPTLambdaBar0=0;
      Double_t pPTRhoP=0;
      Double_t pPTRho0=0;
      Double_t pPTSigmaP=0;
      Double_t pPTSigmaM=0;
      Double_t pPTSigma0=0;
      Double_t pPTDeltaPP=0;
      Double_t pPTDeltaP=0;
      Double_t pPTDeltaM=0;
      Double_t pPTDelta0=0;
      Double_t pPTXiM=0;
      Double_t pPTXi0=0;
      Double_t pPTgamma=0;
      Double_t pPTp=0;
      Double_t pPTn=0;
      Double_t pPTp_=0;
      Double_t pPTn_=0;
      Double_t pPTmu=0;
      Double_t pPTmu_=0;
      Double_t pPTe_=0;
      Double_t pPTe=0;

    // Show how far we got every 100'th event.
    if (event % 100 == 0)
      cout <<"Event # " << event <<endl;
    // Make one event.
    pythia->GenerateEvent();
      hNEvents->Fill(0.5);
/*****************************************************************************
pythia->Pylist(1); when run with a single event will display the main structure
of data from the pythia event. For help with understanding the structure of this
code, uncomment the below line
*****************************************************************************/
      //pythia->Pylist(1); //DEBUG helper
      Int_t npart = particles->GetEntries();
      //printf("Analyse %d Particles\n", npart);
      for (Int_t part=0; part<npart; part++) {
	//TObject *object = particles->At(part);
	//cout<<"I am a "<<object->ClassName()<<endl;
	TMCParticle *MPart = (TMCParticle *) particles->At(part);
  //NOTE:might consider using these instead of below  [preferred for high energy world]
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
    pz*=(-1); //useful to keep numbers positive while calculating
    N='y'; //not used but can be used to make pz negative again later
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
/******************************************************************************
To account for  pT properly, E_T is dependant on type of particle.
For baryons and antibaryons, note that rest mass is included (see below)
******************************************************************************/
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

/****************************************************
The following section counts particles based on the below conditions to ensure no double counting

Particle:        inclusion requirements:
pion+-           final state AND not daughter of decay from eta, omega, k_0L, k_0S, lambda0,lambdaBar0, Omega-,Sigma_+-0, Xi_0-
pion0            not daughter of decay from eta, omega, k_0L, k_0S, lambda0,lambdaBar0, Omega-,Sigma_+-0, Xi_0-
kaon+            ALL INCLUDED
kaon-            not daughter of Omega-
kaon0_L          lifespan>=1cm
kaon0_S          lifespan>=1cm
eta              ALL INCLUDED
omega            ALL INCLUDED
Omega-           lifespan>=1cm
Lambda0          lifespan>1cm AND not daughter of Sigma+-0, Omega-, or Xi0-
LambdaBar0       lifespan>1cm AND not daughter of Sigma+-0, Omega-, or Xi0-
Sigma+-0         lifespan>=1cm
Xi0-             lifespan>=1cm
gamma            final state AND not daughter of pion0, eta, omega
proton           NOT COLLIDING PARTICLES AND not daughter of Sigma+-0
neutron          not daughter of Sigma+-0
muon/bar         ALL INCLUDED
electron/bar     ALL INCLUDED

For reference, mpP is in KF format, Ckf is converted to index (from KF_Code.dat)
****************************************************/
  if (Ckf==58){ //pi+
    if ((mpartD==0)&&(mpP!=221)&&(mpP!=223)&&(mpP!=130)&&(mpP!=310)&&(mpP!=3334)&&(mpP!=3122)&&(mpP!=3122)&&(mpP!=3222)&&(mpP!=-3222)&&(mpP!=3212)&&(mpP!=3312)&&(mpP!=3322)){
      ETpip+=partE;
      ETAll+=partE;
      PTpip+=pT;
      PTAll+=pT;
      cpip++;
  }}
  if (Ckf==-58){ //pi-
    if ((mpartD==0)&&(mpP!=221)&&(mpP!=223)&&(mpP!=130)&&(mpP!=310)&&(mpP!=3334)&&(mpP!=3122)&&(mpP!=3122)&&(mpP!=3222)&&(mpP!=-3222)&&(mpP!=3212)&&(mpP!=3312)&&(mpP!=3322)){
    ETpim+=partE;
    ETAll+=partE;
    PTpim+=pT;
    PTAll+=pT;
    cpim++;
  }}
  if (Ckf==68){ //pi0
    if ((mpP!=221)&&(mpP!=223)&&(mpP!=130)&&(mpP!=310)&&(mpP!=3334)&&(mpP!=3122)&&(mpP!=3122)&&(mpP!=3222)&&(mpP!=-3222)&&(mpP!=3212)&&(mpP!=3312)&&(mpP!=3322)){
    ETpi0+=partE;
    ETAll+=partE;
    PTpi0+=pT;
    PTAll+=pT;
    cpi0++;
  }}
  if (Ckf==60){ //K+
    if (mpartD==0){//counted if final particle
    ETKp+=partE;
    ETAll+=partE;
    PTKp+=pT;
    PTAll+=pT;
    cKp++;
  }}
  if (Ckf==-60){ //K-
    if ((mpartD==0)&&(mpP!=3334)){//counted if final particle/ not daughter of Omega-
    ETKm+=partE;
    ETAll+=partE;
    PTKm+=pT;
    PTAll+=pT;
    cKm++;
  }}
  if (Ckf==73){// KL
    if ((mpL>=100)||(mpartD==0)){//counted if final particle or lifetime > 1cm
    ETKL+=partE;
    ETAll+=partE;
    PTKL+=pT;
    PTAll+=pT;
    cKL++;
  }}
  if (Ckf==74){ //KS
    if ((mpL>=100)||(mpartD==0)){//counted if final particle or lifetime > 1cm
    ETKS+=partE;
    ETAll+=partE;
    PTKS+=pT;
    PTAll+=pT;
    cKS++;
  }}
  if (Ckf==69){ //eta // include ALLLLLL of them!!!
    ETEta+=partE;
    ETAll+=partE;
    PTEta+=pT;
    PTAll+=pT;
    cEta++;
  }
  if (Ckf==86){ // omega include all of them.
    ETOmega+=partE;
    ETAll+=partE;
    PTOmega+=pT;
    PTAll+=pT;
    cOmega++;
  }
  if (Ckf==162){ // Omega-
    if (((mpL>=100)||(mpartD==0))){
    ETOmegaM+=partE;
    ETAll+=partE;
    PTOmegaM+=pT;
    PTAll+=pT;
    cOmegaM++;
  }}
  if (Ckf==135){ //Lambda0
    if (((mpL>=100)||(mpartD==0))&&(mpP!=3334)&&(mpP!=3122)&&(mpP!=3122)&&(mpP!=3222)&&(mpP!=-3222)&&(mpP!=3212)&&(mpP!=3312)&&(mpP!=3322)){
    ETLambda0+=partE;
    ETAll+=partE;
    PTLambda0+=pT;
    PTAll+=pT;
    cLambda0++;
  }}
  if (Ckf==-135){ //LambdaBar0
    if ((mpL>=100)||(mpartD==0)&&(mpP!=3334)&&(mpP!=3122)&&(mpP!=3122)&&(mpP!=3222)&&(mpP!=-3222)&&(mpP!=3212)&&(mpP!=3312)&&(mpP!=3322)){
    ETLambdaBar0+=partE;
    ETAll+=partE;
    PTLambdaBar0+=pT;
    PTAll+=pT;
    cLambdaBar0++;
  }}
  if (Ckf==137){ //Sigma+
    if (((mpL>=100)||(mpartD==0))){
    ETSigmaP+=partE;
    ETAll+=partE;
    PTSigmaP+=pT;
    PTAll+=pT;
    cSigmaP++;
  }}
  if (Ckf==134){ //Sigma-
    if (((mpL>=100)||(mpartD==0))){
    ETSigmaM+=partE;
    ETAll+=partE;
    PTSigmaM+=pT;
    PTAll+=pT;
    cSigmaM++;
  }}
  if (Ckf==136){ //Sigma0
    if (((mpL>=100)||(mpartD==0))){
    ETSigma0+=partE;
    ETAll+=partE;
    PTSigma0+=pT;
    PTAll+=pT;
    cSigma0++;
  }}
  if (Ckf==138){ //Xi-
    if (((mpL>=100)||(mpartD==0))){
    ETXiM+=partE;
    ETAll+=partE;
    PTXiM+=pT;
    PTAll+=pT;
    cXiM++;
  }}
  if (Ckf==139){ //Xi0
    if (((mpL>=100)||(mpartD==0))){
    ETXi0+=partE;
    ETAll+=partE;
    PTXi0+=pT;
    PTAll+=pT;
    cXi0++;
  }}
  if (Ckf==17){ //gamma
    if ((mpartD==0)&&(mpP!=111)&&(mpP!=221)&&(mpP!=223)){ //record gamma only if not produced by pi0, eta, omega
    ETgamma+=partE;
    ETAll+=partE;
    PTgamma+=pT;
    PTAll+=pT;
    cgamma++;
  }}
  if (Ckf==133){ //proton
    if (mpartD==0){
    if ((part!=0)&&(part!=1)&&(mpP!=3222)&&(mpP!=-3222)&&(mpP!=3212)){
      ETp+=partE;
      ETAll+=partE;
      PTp+=pT;
      PTAll+=pT;
      cp++;
    }}
  }
  if (Ckf==132){ //neutron
    if ((mpartD==0)&&(mpP!=2112)&&(mpP!=3222)&&(mpP!=-3222)&&(mpP!=3212)){
    ETn+=partE;
    ETAll+=partE;
    PTn+=pT;
    PTAll+=pT;
    cn++;
  }}
  if (Ckf==-133){ //antiproton
    if (mpartD==0){
    if ((part!=0)&&(part!=1)&&(mpP!=3222)&&(mpP!=-3222)&&(mpP!=3212)){
      ETp_+=partE;
      ETAll+=partE;
      PTp_+=pT;
      PTAll+=pT;
      cp_++;
    }}
  }
  if (Ckf==-132){ //antineutron
    if ((mpartD==0)&&(mpP!=2112)&&(mpP!=3222)&&(mpP!=-3222)&&(mpP!=3212)){
    ETn_+=partE;
    ETAll+=partE;
    PTn_+=pT;
    PTAll+=pT;
    cn_++;
  }}
  if (Ckf==10){ //mu
    if (mpartD==0){
    ETmu+=partE;
    ETAll+=partE;
    PTmu+=pT;
    PTAll+=pT;
    cmu++;
  }}
  if (Ckf==11){ //muBar
    if (mpartD==0){
    ETmu_+=partE;
    ETAll+=partE;
    PTmu_+=pT;
    PTAll+=pT;
    cmu_++;
  }}
  if (Ckf==9){ //positron
    if (mpartD==0){
    ETe_+=partE;
    ETAll+=partE;
    PTe_+=pT;
    PTAll+=pT;
    ce_++;
  }}
  if (Ckf==8){ //electron
    if (mpartD==0){
    ETe+=partE;
    ETAll+=partE;
    PTe+=pT;
    PTAll+=pT;
    ce++;
  }}

} //end pseudorapidity cut

      }//end particle loop
/******************************************************************************
Histogarms are filled here, the histograms are of the format:
hET<particleName> xaxis is ET particle/ Total ET
                  yaxis is number of events at that fraction
*******************************************************************************/
      Double_t omegaET=0;
      Double_t pi0ET=0;
      hETAll->Fill(ETAll);
      //cout<<"TOTAL: "<<ETAll<<endl;
      if (ETpip!=0){
      pETpip=ETpip/ETAll;
      hETpip->Fill(pETpip);
      hPTpip->Fill(PTpip);

        //cout<<"pip "<<ETpip<<endl;
      }
      if (ETpim!=0){
      pETpim=ETpim/ETAll;
      hETpim->Fill(pETpim);
      hPTpim->Fill(PTpim);
        //cout<<"pim "<<ETpim<<endl;
      }
      if (ETpi0!=0){
      pETpi0=ETpi0/ETAll;
      hETpi0->Fill(pETpi0);
      hPTpi0->Fill(PTpi0);
        //cout<<"pi0 "<<ETpi0<<endl;
      }
      if (ETKp!=0){
      pETKp=ETKp/ETAll;
      hETKp->Fill(pETKp);
      hPTKp->Fill(PTKp);
        //cout<<"Kp "<<ETKp<<endl;
      }
      if (ETKm!=0){
      pETKm=ETKm/ETAll;
      hETKm->Fill(pETKm);
      hPTKm->Fill(PTKm);
        //cout<<"Kp "<<ETKp<<endl;
      }
      if (ETKL!=0){
      pETKL=ETKL/ETAll;
      hETKL->Fill(pETKL);
      hPTKL->Fill(PTKL);
        //cout<<"KL "<<ETKL<<endl;
      }
      if (ETKS!=0){
      pETKS=ETKS/ETAll;
      hETKS->Fill(pETKS);
      hPTKS->Fill(PTKS);
        //cout<<"KS "<<ETKS<<endl;
      }
      if (ETEta!=0){
      pETEta=ETEta/ETAll;
      hETEta->Fill(pETEta);
      hPTEta->Fill(PTEta);
        //cout<<"Eta "<<ETEta<<endl;
      }
      if (ETOmega!=0){
      pETOmega=ETOmega/ETAll;
      hETOmega->Fill(pETOmega);
      hPTOmega->Fill(PTOmega);
        //cout<<"Omega "<<ETOmega<<endl;
      }
      if (ETOmegaM!=0){
      pETOmegaM=ETOmegaM/ETAll;
      hETOmegaM->Fill(pETOmegaM);
      hPTOmegaM->Fill(PTOmegaM);
        //cout<<"OmegaM "<<ETOmegaM<<endl;
      }
      if (ETLambda0!=0){
      pETLambda0=ETLambda0/ETAll;
      hETLambda0->Fill(pETLambda0);
      hPTLambda0->Fill(PTLambda0);
        //cout<<"Lambda0 "<<ETLambda0<<endl;
      }
      if (ETLambdaBar0!=0){
      pETLambdaBar0=ETLambdaBar0/ETAll;
      hETLambdaBar0->Fill(pETLambdaBar0);
      hPTLambdaBar0->Fill(PTLambdaBar0);
        //cout<<"Lambda0 "<<ETLambda0<<endl;
      }
      if (ETSigmaP!=0){
      pETSigmaP=ETSigmaP/ETAll;
      hETSigmaP->Fill(pETSigmaP);
      hPTSigmaP->Fill(PTSigmaP);
        //cout<<"SigmaP "<<ETSigmaP<<endl;
      }
      if (ETSigmaM!=0){
      pETSigmaM=ETSigmaM/ETAll;
      hETSigmaM->Fill(pETSigmaM);
      hPTSigmaM->Fill(PTSigmaM);
        //cout<<"SigmaM "<<ETSigmaM<<endl;
      }
      if (ETSigma0!=0){
      pETSigma0=ETSigma0/ETAll;
      hETSigma0->Fill(pETSigma0);
      hPTSigma0->Fill(PTSigma0);
        //cout<<"Sigma0 "<<ETSigma0<<endl;
      }
      if (ETXiM!=0){
      pETXiM=ETXiM/ETAll;
      hETXiM->Fill(pETXiM);
      hPTXiM->Fill(PTXiM);
        //cout<<"Xi- "<<ETXiM<<endl;
      }
      if (ETXi0!=0){
      pETXi0=ETXi0/ETAll;
      hETXi0->Fill(pETXi0);
      hPTXi0->Fill(PTXi0);
        //cout<<"Xi0 "<<ETXiM<<endl;
      }
      if (ETgamma!=0){
      pETgamma=ETgamma/ETAll;
      hETgamma->Fill(pETgamma);
      hPTgamma->Fill(PTgamma);
        //cout<<"gamma "<<ETgamma<<endl;
      }
      if (ETp!=0){
      pETp=ETp/ETAll;
      hETp->Fill(pETp);
      hPTp->Fill(PTp);
        //cout<<"p "<<ETp<<endl;
      }
      if (ETn!=0){
      pETn=ETn/ETAll;
      hETn->Fill(pETn);
      hPTn->Fill(PTn);
        //cout<<"n "<<ETn<<endl;
      }
      if (ETp_!=0){
      pETp_=ETp_/ETAll;
      hETp_->Fill(pETp_);
      hPTp_->Fill(PTp_);
        //cout<<"p "<<ETp<<endl;
      }
      if (ETn_!=0){
      pETn_=ETn_/ETAll;
      hETn_->Fill(pETn_);
      hPTn_->Fill(PTn_);
        //cout<<"n "<<ETn<<endl;
      }
      if (ETmu!=0){
      pETmu=ETmu/ETAll;
      hETmu->Fill(pETmu);
      hPTmu->Fill(PTmu);
        //cout<<"mu "<<ETmu<<endl;
      }
      if (ETmu_!=0){
      pETmu_=ETmu_/ETAll;
      hETmu_->Fill(pETmu_);
      hPTmu_->Fill(PTmu_);//set marker color/style
        //cout<<"muBAR "<<ETmu_<<endl;
      }
      if (ETe!=0){
      pETe=ETe/ETAll;
      hETe->Fill(pETe);
      hPTe->Fill(PTe);
        //cout<<"e "<<ETe<<endl;
      }
      if (ETe_!=0){
      pETe_=ETe_/ETAll;
      hETe_->Fill(pETe_);
      hPTe_->Fill(PTe_);
        //cout<<"eBAR "<<ETe_<<endl;
      }
      Float_t Rat1=ETOmega/ETpi0;
      Float_t Rat2=ETEta/ETpi0;
      if (Rat1!=0){

      }

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
  Double_t ETOmega=0; //THIS IS omega not Omega
  /*
  Double_t ETOmegaM=0;
  Double_t ETLambda0=0;
  Double_t ETLambdaBar0=0;
  Double_t ETSigmaP=0;
  Double_t ETSigmaM=0;
  Double_t ETSigma0=0;
  Double_t ETXiM=0;
  Double_t ETXi0=0;
  Double_t ETgamma=0;
  Double_t ETp=0;
  Double_t ETn=0;
  Double_t ETmu=0;
  Double_t ETmu_=0;
  Double_t ETe_=0;
  Double_t ETe=0;*/
//NOTE: turn into pt graphs for pi0 and eta (and omega), normalize by multiplying by 1/2pi divide histos
hPTpip=hPTpip*((1.)/(2*TMath::pi())); //n is bin width as histogram
hPTpim=hPTpim*((1.)/(2*TMath::pi()));
hPTpi0=hPTpi0*((1.)/(2*TMath::pi()));
hPTKp=hPTKp*((1.)/(2*TMath::pi()));
hPTKm=hPTKm*((1.)/(2*TMath::pi()));
hPTKL=hPTKL*((1.)/(2*TMath::pi()));
hPTKS=hPTKS*((1.)/(2*TMath::pi()));
hPTEta=hPTEta*((1.)/(2*TMath::pi()));
hPTOmega=hPTOmega*((1.)/(2*TMath::pi()));
hPTOmegaM=hPTOmegaM*((1.)/(2*TMath::pi()));
hPTLambda0=hPTLambda0*((1.)/(2*TMath::pi()));
hPTLambdaBar0=hPTLambdaBar0*((1.)/(2*TMath::pi()));
hPTSigma0=hPTSigma0*((1.)/(2*TMath::pi()));
hPTSigmaP=hPTSigmaP*((1.)/(2*TMath::pi()));
hPTSigmaM=hPTSigmaM*((1.)/(2*TMath::pi()));
hPTXiM=hPTXiM*((1.)/(2*TMath::pi()));
hPTXi0=hPTXi0*((1.)/(2*TMath::pi()));
hPTgamma=hPTgamma*((1.)/(2*TMath::pi()));
hPTp=hPTp*((1.)/(2*TMath::pi()));
hPTn=hPTn*((1.)/(2*TMath::pi()));
hPTp_=hPTp_*((1.)/(2*TMath::pi()));
hPTn_=hPTn_*((1.)/(2*TMath::pi()));
hPTmu=hPTmu*((1.)/(2*TMath::pi()));
hPTmu_=hPTmu_*((1.)/(2*TMath::pi()));
hPTe=hPTe*((1.)/(2*TMath::pi()));
hPTe_=hPTe_*((1.)/(2*TMath::pi()));
//just use raw number over other, dont need pt spectra

  //NOTE:write histograms to file
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
    ETEta=hETEta->GetMean();
    ETOmega=hETOmega->GetMean();
    /*ETOmegaM=hETOmegaM->GetMean();
    ETLambda0=hETLambda0->GetMean();
    ETLambdaBar0=hETLambdaBar0->GetMean();
    ETSigmaP=hETSigmaP->GetMean();
    ETSigmaM=hETSigmaM->GetMean();
    ETSigma0=hETSigma0->GetMean();
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
    Float_t Ratio=(ETOmega*cOmega)/(ETpi0*cpi0);
    Float_t Ratio2=(ETEta*cEta)/(ETpi0*cpi0);
    cout<<"ET omega to pion0 :"<<Ratio<<endl;
    cout<<"ET Eta to pion0 :"<<Ratio2<<endl;
    hETAll->Write();
    hETpip->Write();
    hETpim->Write();
    hETpi0->Write();
    hETKp->Write();
    hETKm->Write();
    hETKL->Write();
    hETKS->Write();
    hETEta->Write();
    hETOmega->Write();
    hETOmegaM->Write();
    hETLambda0->Write();
    hETLambdaBar0->Write();
    hETSigmaP->Write();
    hETSigmaM->Write();
    hETSigma0->Write();
    hETXiM->Write();
    hETXi0->Write();
    hETgamma->Write();
    hETp->Write();
    hETn->Write();
    hETp_->Write();
    hETn_->Write();
    hETmu->Write();
    hETmu_->Write();
    hETe->Write();
    hETe_->Write();


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
    hPTOmegaM->Write();
    hPTLambda0->Write();
    hPTLambdaBar0->Write();
    hPTSigmaP->Write();
    hPTSigmaM->Write();
    hPTSigma0->Write();
    hPTXiM->Write();
    hPTXi0->Write();
    hPTgamma->Write();
    hPTp->Write();
    hPTn->Write();
    hPTp_->Write();
    hPTn_->Write();
    hPTmu->Write();
    hPTmu_->Write();
    hPTe->Write();
    hPTe_->Write();

    outfile->Close();

  return 0;
}

// Show the Pt spectra, and start the tree viewer.
int showEventSample()
{


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

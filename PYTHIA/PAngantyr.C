// main113.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This test program will generate Pb-Pb collisions at
// sqrt(S_NN)=2.76TeV using the Angantyr model for Heavy Ion
// collisions. The analysis will divide the event in centrality
// classes using the same observable as was used for p-Pb in the ATLAS
// analysis in arXiv:1508.00848 [hep-ex] (see main112.cc). The
// centrality classes are same as in the ALICE analysis in
// arXiv:1012.1657 [nucl-ex] although the actual observable used is
// not the same. Histograms of multiplicity distributions are measured
// for each centrality percentile.

// Note that heavy ion collisions are computationally quite CPU
// intensive and generating a single event will take around a second
// on a reasonable desktop. To get reasonable statistics, this program
// will take a couple of hours to run.

/*
In order to run, ypou may need to adjust PythiaStdlib.h ln26

// Stdlib header file for dynamic library loading.
   #ifndef __CINT__
   #define dlsym __
   #include <dlfcn.h>
   #undef dlsym
   #endif
*/


#ifndef __CINT__
#include "TApplication.h"
#include "Pythia8/Pythia.h"
#include "Pythia8/HeavyIons.h"
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
#include <ctime>
using namespace Pythia8;
#endif

#define FILENAME   "pythia.root"
#define TREENAME   "tree"
#define BRANCHNAME "particles"
#define HISTNAME   "ptSpectra"
#define PDGNUMBER  211
class TH3F;
class TH1F;

struct KF_Code {
  Int_t KFc;
  TString name;
};
struct pylista {
  Int_t inde;
  Int_t KFpart;
  Int_t indePar;
};

int test(){
return 0;
}

TH1F *CreateEnergyHistogram(char *name){
  TH1F *histo = new TH1F(name,name,1000,0,10);
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
TH1F *CreatePTHistogram(char *name){
  TH1F *histo = new TH1F(name,name,1000,0,10);
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

Int_t GetDaughterCheck(vector<pylista> pylis,Int_t ind, Int_t INPUT, Int_t CHECK){
  int s=pylis.size();
  int t=0;
  int R=15;
  vector <int> d(3);
  for (int i=0;i<s;i++){
    if (pylis[i].indePar==ind){
      d[t]=pylis[i].KFpart;
      t++;
      //cout<<pylis[i].inde<<"\t"<<pylis[i].KFpart<<"\t"<<pylis[i].indePar<<endl;
    }
  }
  //cout<<d[0]<<" "<<d[1]<<" "<<d[2]<<endl;
  if (INPUT==111){ //eta
    if ((d[0]==CHECK)||(d[1]==CHECK)){
      R=1;
    }
    else
      R=0;
  }
  if (INPUT==221){ //eta
    if ((d[0]==CHECK)&&(d[1]==CHECK)){
      R=1;
    }
    else
      R=0;
  }
  else if (INPUT==223) { //omega
    if ((d[0]==CHECK)||(d[1]==CHECK)){
      R=1;
    }
    else
      R=0;
  }
  else if (INPUT==3122) { //lambda
    if ((d[0]==CHECK)||(d[1]==CHECK)){
      R=1; //there is a proton
    }
    else
      R=0;
  }
  return R;
}



int makeEventSample(Int_t nEvents, Int_t jobID, Int_t tune, Float_t SNN, char CUT, Float_t yncut, char DTYPE, char DMODE) {
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
  else if (SNN==62.4){
    modSNN=62;
  }
  else {
    modSNN=SNN;
  }
  cout<<modSNN<<endl;
  cout<<"I made it here, running "<<nEvents<<" events, job id "<<jobID<<", tune "<<tune<<", at "<<SNN<<"GeV"<<endl;
  char *filename = Form("%i_%iGeVoutfile.root",jobID,modSNN);
  char *filename2 = Form("%i_%iGeVoutfile.txt",jobID,modSNN);
  ifstream in;
  in.open("KF_Code.dat");
  while (in.good()){
    partname.push_back(KF_Code());
    in>>partname[l].KFc>>partname[l].name;
    l++;
  }
  ofstream out;
  out.open(filename2);
  cout<<"I made it here line 200"<<endl;

  //Create new instance of the Pythia event generator
  string pVERSION;
  Pythia pythia;
  pVERSION="Pythia8 Angantyr";

  cout<<"I made it here line 207"<<endl;

  string temporarything="fail";
  out<<"RUN DATA\n******************************************************\n";
  out<<"JobID: "<<jobID<<"\tNo. Events: "<<nEvents<<"\tSNN: "<<SNN<<"\t tune: "<<tune<<endl<<endl;
  out<<"Simulation Parameters:\nVersion: "<<pVERSION<<"\nType of Cut: ";
  if (CUT=='y'){
    temporarything = "rapidity (y)";
  }
  else if (CUT=='n'){
    temporarything = "pseudorapidity (n)";
  }
  out<<temporarything<<"\nCut: "<<yncut<<"\nCalorimeter (r=regular, c=calorimeter): "<<DTYPE<<endl;
  out<<"Counting method: "<<DMODE<<endl<<endl;

  cout<<"I made it here line 222"<<endl;
  // Setup the beams.
  char* SNNthing=Form("Beams:eCM = %f",SNN);
  pythia.readString("Beams:idA = 1000791970"); //sets gold (Au-197)
  pythia.readString("Beams:idB = 1000791970");
  pythia.readString(SNNthing); //SNN
  pythia.settings.addFlag("Main:writeRoot",false);
  pythia.readString("Main:timesAllowErrors = 10000000");
  pythia.readString("Random:setSeed = on");
  // using PDG codes of the format 100ZZZAAAI: 1000020040 = 4He , 1000030060 = 6Li, 1000060120 = 12C, 1000080160 = 16O, 1000290630 = 63Cu, 1000791970 = 197Au, and 1000822080 = 208Pb.
  cout<<"I made it here line 229"<<endl;
  // Initialize the Angantyr model to fit the total and semi-includive
  // cross sections in Pythia within some tolerance.
  pythia.readString("HeavyIon:SigFitErr = "
                    "0.02,0.02,0.1,0.05,0.05,0.0,0.1,0.0");
  // These parameters are typicall suitable for sqrt(S_NN)=5TeV
  pythia.readString("HeavyIon:SigFitDefPar = "
                    "17.24,2.15,0.33,0.0,0.0,0.0,0.0,0.0");
  // A simple genetic algorithm is run for 20 generations to fit the
  // parameters.
  pythia.readString("HeavyIon:SigFitNGen = 20");
  cout<<"I made it here line 240"<<endl;
  // There will be nine centrality bins based on the sum transverse
  // emergy in a rapidity interval between 3.2 and 4.9 obtained from
  // the borders from the generated transverse energy spectrum. The
  // default settings should give approximately the following:
  double genlim[] = {2979.4, 2400.1, 1587.5, 1028.8, 669.9,
                     397.4, 220.3, 116.3, 54.5};
  // If you change any parameters these should also be changed.

  // The upper edge of the correponding percentiles:
  double pclim[] = {0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};

  // Book the pseudorapidity and multiplicity histograms and get
  // counters for number of events and sum of event weights:
  typedef map<double,int,std::greater<double> > MapIdx;
  MapIdx genetaidx;
  vector<Hist*> etadist(9), lmult(9), hmult(9);
  string etaname("EtadistC"), mname("MultC");
  vector<double> gensumw(9, 0.0), gensumn(9, 0.0);
  for ( int i = 0; i < 9; ++i ) {
    genetaidx[genlim[i]] = i;
    etadist[i] = new Hist(etaname + char('0' + i), 54, -2.7, 2.7);
    hmult[i] = new Hist(mname + 'H' + char('0' + i),
                        75, -0.5, 2999.5);
    lmult[i] = new Hist(mname + 'L' + char('0' + i),
                        75, -0.5, 299.5);
  }

  // Book histogram for the centrality measure.
  Hist sumet("SumETfwd", 200, 0.0, 4000.0);

  // Also make a map of all weight to check the generated centrality
  // classes.
  multimap<double,double> gencent;

  // Book a histogram for the distribution of number of wounded
  // nucleons.
  Hist wounded("Nwounded", 209, -0.5, 417.5);

  // Profile for average central multiplicity and number of wounded
  // nucleons as a function of centrality (with errors).
  vector<double> cmult(9, 0.0), cmult2(9, 0.0);
  vector<double> wound(9, 0.0), wound2(9, 0.0);

  // Sum up the weights of all generated events.
  double sumw = 0.0;
  char* IDthing=Form("Random:seed = %i",jobID);
  pythia.readString(IDthing);
  pythia.readString("Beams:frameType = 1");
  char* Eventthing=Form("Next:numberCount = %i",nEvents);
  pythia.readString(Eventthing);
  pythia.readString("PartonLevel:all = 1");
  pythia.readString("ProcessLevel:all = 1");
  pythia.readString("Tune:ee = 1");
  pythia.readString("Tune:pp = 1");
  pythia.readString("Next:numberShowEvent = 1");

  cout<<"I made it here line 286"<<endl;
  // Initialise pythia.
  pythia.init();

  cout<<"I made it here line 290"<<endl;
  TFile* file = TFile::Open(filename, "RECREATE");
  if (!file || !file->IsOpen()) {
    Error("makeEventSample", "Couldn;t open file %s", filename);
    return 1;
  }

  TFile *outfile = file;//new TFile(filename,"RECREATE");
  // Loop over events.
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
  TH1F *hETOMEGAm = CreateEnergyHistogram("hETOMEGAm"); //3334
  TH1F *hETXi0 = CreateEnergyHistogram("hETXi0"); //3322
  TH1F *hETXim = CreateEnergyHistogram("hETXim"); //3312
  TH1F *hETSigmap = CreateEnergyHistogram("hETSigmap"); //3222
  TH1F *hETSigmam = CreateEnergyHistogram("hETSigmam"); //3112
  TH1F *hETSigma0 = CreateEnergyHistogram("hETSigma0"); //3212
  TH1F *hETp = CreateEnergyHistogram("hETp");
  TH1F *hETn = CreateEnergyHistogram("hETn");
  TH1F *hETp_ = CreateEnergyHistogram("hETpBar");
  TH1F *hETn_ = CreateEnergyHistogram("hETnBar");
  TH1F *hETother = CreateEnergyHistogram("hETother");

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
  TH1F *hPTOMEGAm = CreatePTHistogram("hptOMEGAm"); //3334
  TH1F *hPTXi0 = CreatePTHistogram("hptXi0"); //3322
  TH1F *hPTXim = CreatePTHistogram("hptXim"); //3312
  TH1F *hPTSigmap = CreatePTHistogram("hptSigmap"); //3222
  TH1F *hPTSigmam = CreatePTHistogram("hptSigmam"); //3112
  TH1F *hPTSigma0 = CreatePTHistogram("hptSigma0"); //3212
  TH1F *hPTp = CreatePTHistogram("hptp");
  TH1F *hPTn = CreatePTHistogram("hptn");
  TH1F *hPTp_ = CreatePTHistogram("hptpBar");
  TH1F *hPTn_ = CreatePTHistogram("hptnBar");
  TH1F *hPTother = CreatePTHistogram("hptother");

  //NOTE: VERY IMPORTANT
  TH1F *hNEvents = new TH1F("hNEvents","Number of events",1,0,1.0);
  hNEvents->GetYaxis()->SetTitle("N_{events}");
  hNEvents->GetXaxis()->SetTitle("no title");


  //tree->Branch(BRANCHNAME, &particles);
  cout<<"I made it here line 359"<<endl;
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
  Int_t cLambda0=0;
  Int_t cLambdaBar0=0;
  Int_t cOMEGAm=0;
  Int_t cXi0=0;
  Int_t cXim=0;
  Int_t cSigmap=0;
  Int_t cSigmam=0;
  Int_t cSigma0=0;
  Int_t cp=0;
  Int_t cn=0;
  Int_t cp_=0;
  Int_t cn_=0;
  Int_t cother=0;
  Int_t cETAall=0;
  Int_t cOmegaall=0;

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
  Double_t pETOMEGAm=0;
  Double_t pETXi0=0;
  Double_t pETXim=0;
  Double_t pETSigmap=0;
  Double_t pETSigmam=0;
  Double_t pETSigma0=0;
  Double_t pETp=0;
  Double_t pETn=0;
  Double_t pETp_=0;
  Double_t pETn_=0;
  Double_t pETother=0;
  Double_t pETETAall=0;
  Double_t pETOmegaall=0;

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
  Double_t pPTOMEGAm=0;
  Double_t pPTXi0=0;
  Double_t pPTXim=0;
  Double_t pPTSigmap=0;
  Double_t pPTSigmam=0;
  Double_t pPTSigma0=0;
  Double_t pPTp=0;
  Double_t pPTn=0;
  Double_t pPTp_=0;
  Double_t pPTn_=0;
  Double_t pPTother=0;
  Double_t pPTETAall=0;
  Double_t pPTOmegaall=0;


  for ( int event = 0; event < nEvents; ++event ) {
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
    Double_t ETOMEGAm=0;
    Double_t ETXi0=0;
    Double_t ETXim=0;
    Double_t ETSigmap=0;
    Double_t ETSigmam=0;
    Double_t ETSigma0=0;
    Double_t ETp=0;
    Double_t ETn=0;
    Double_t ETp_=0;
    Double_t ETn_=0;
    Double_t ETother=0;


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
      Double_t PTOMEGAm=0;
      Double_t PTXi0=0;
      Double_t PTXim=0;
      Double_t PTSigmap=0;
      Double_t PTSigmam=0;
      Double_t PTSigma0=0;
      Double_t PTp=0;
      Double_t PTn=0;
      Double_t PTp_=0;
      Double_t PTn_=0;//int LambDecayCount=0;
      Double_t PTother=0;

      vector<pylista> pylis(1);
    // Show how far we got every 100'th event.
    if (event % 10 == 0)
      cout <<"Collision Energy: "<< SNN <<" Event # " << event <<endl;
    // Make one event.
    //pythia.GenerateEvent(); //?

    if ( !pythia.next() ) continue;
      hNEvents->Fill(0.5);
    Event TheEvent = pythia.event;
    Int_t npart = TheEvent.size();
    for (Int_t parta=0; parta<npart; parta++) {
      Particle* particle =& TheEvent.at(parta);
      pylis.push_back(pylista());
      pylis[parta].inde=parta;
      pylis[parta].indePar=particle->mother1();
      pylis[parta].KFpart=TheEvent.at(pylis[parta].indePar).id();
    }
    for (Int_t part=0; part<npart; part++) {
      //TObject *object = particles->At(part);
      //cout<<"I am a "<<object->ClassName()<<endl;
      Particle* MPart =& TheEvent.at(part);
      Int_t KFID = MPart->id();
      Int_t Ckf=GetKFConversion(KFID,partname);
      Int_t ind = part+1; //index
      //hEve->Fill(Ckf,event,ind);
      char N ='n';
      Float_t E= MPart ->e();
      Float_t partE;
      Float_t px=MPart->px();
      Float_t py=MPart->py();
      Float_t pz=MPart->pz();
      if (pz<0){
        pz*=(-1); //useful to keep numbers positive while calculating
        N='y'; //not used but can be used to make pz negative again later
      }
      Float_t pT= sqrt(pow(px,2)+pow(py,2));
      Float_t Theta = atan(pT/pz);
      Float_t p_tot=(pT)/sin(Theta);
      Float_t m=MPart->m();
      Float_t Ei_TOT= sqrt(pow(p_tot,2)+pow(m,2));
      Float_t pseudorapidity=-log(tan(Theta/2));
      Float_t rapidity=(0.5)*log((E+pz)/(E-pz));
      //cout<<pseudorapidity<<" "<<rapidity<<endl;
      Int_t mpartD = MPart->daughter1();
      Int_t mpartD2 = MPart->daughter2();
      Int_t mpP = MPart->mother1();
      Float_t mpL = MPart->tau0();
      Int_t nobby=15;
      Int_t XEM=31;
      Float_t cuttype;
      if (CUT=='y'){
        cuttype=rapidity;
      }
      else if (CUT=='n'){
        cuttype=pseudorapidity;
      }
      if (cuttype<yncut){
/******************************************************************************
To account for  pT properly, E_T is dependant on type of particle.
For baryons and antibaryons, note that rest mass is included (see below)
******************************************************************************/
      if (DTYPE=='r'){
        if ((KFID==3122)||(KFID==2212)||(KFID==2112)){
          partE = (Ei_TOT-m) * sin(Theta);
        }
        else if ((KFID==-3122)||(KFID==-2212)||(KFID==-2112)){
          partE = (Ei_TOT+m) * sin(Theta);
        }
        else {
          partE = Ei_TOT * sin(Theta);
        }
      }
      else if (DTYPE=='c')
      partE=Ei_TOT;
      if (DMODE=='a'){ //ALL INCLUDED
/****************************************************
The following section counts particles based on the below conditions to ensure no double counting

Particle:        inclusion requirements:
pion+-           ALL INCLUDED
pion0            ALL INCLUDED
kaon+            ALL INCLUDED
kaon-            ALL INCLUDED
kaon0_L          ALL INCLUDED
kaon0_S          ALL INCLUDED
eta              ALL INCLUDED
omega            ALL INCLUDED
Lambda0          ALL INCLUDED
LambdaBar0       ALL INCLUDED
Omega-           ALL INCLUDED
Xi0              ALL INCLUDED
Xi-              ALL INCLUDED
Sigma+           ALL INCLUDED
Sigma0           ALL INCLUDED
Sigma-           ALL INCLUDED
proton           ALL INCLUDED
neutron          ALL INCLUDED

****************************************************/

  if (KFID==211){ //pi+ all included
    ETpip+=partE; //repeated for event
    pETpip+=partE; //repeated over total run
    ETAll+=partE;
    ETcharge+=partE;
    PTpip+=pT;
    pPTpip+=pT;
    PTAll+=pT;
    PTcharge=+pT;
    cpip++;  //repeated over total run
  }
  if (KFID==-211){ //pi- all included
    ETpim+=partE;
    pETpim+=partE;
    ETAll+=partE;
    ETcharge+=partE;
    PTpim+=pT;
    pPTpim+=pT;
    PTAll+=pT;
    PTcharge=+pT;
    cpim++;
  }
  if (KFID==111){ //pi0 // BUG double counts for each omega counted!!!!!
    ETpi0+=partE;
    pETpi0+=partE;
    ETAll+=partE;
    PTpi0+=pT;
    pPTpi0+=pT;
    PTAll+=pT;
    cpi0++;
  }
  if (KFID==321){ //K+ //all included
    ETKp+=partE;
    pETKp+=partE;
    ETAll+=partE;
    ETcharge+=partE;
    PTKp+=pT;
    pPTKp+=pT;
    PTAll+=pT;
    PTcharge=+pT;
    cKp++;
  }
  if (KFID==-321){ //K- //all included
    ETKm+=partE;
    pETKm+=partE;
    ETAll+=partE;
    ETcharge+=partE;
    PTKm+=pT;
    pPTKm+=pT;
    PTAll+=pT;
    PTcharge=+pT;
    cKm++;
  }
  if (KFID==130){// KL
    //if ((mpL>=100)||(mpartD==0)){//counted if final particle or lifetime > 1cm
    ETKL+=partE;
    pETKL+=partE;
    ETAll+=partE;
    PTKL+=pT;
    pPTKL+=pT;
    PTAll+=pT;
    cKL++;
  }
  if (KFID==310){ //KS
    //if ((mpL>=100)||(mpartD==0)){//counted if final particle or lifetime > 1cm
    ETKS+=partE;
    pETKS+=partE;
    ETAll+=partE;
    PTKS+=pT;
    pPTKS+=pT;
    PTAll+=pT;
    cKS++;
  }
  if (KFID==221){ //eta
    //cout<<ind<<"= particle number (check eta) running func:"<<endl;
    //nobby=GetDaughterCheck(pylis,ind,221,22);
    //cout<<nobby<<endl;
    //  if (nobby==1){
    ETEta+=partE;
    pETEta+=partE;
    ETAll+=partE;
    PTEta+=pT;
    pPTEta+=pT;
    PTAll+=pT;
    cEta++;
  }
  if (KFID==223){ // omega
    //cout<<ind<<"= particle number (check omega) running func:"<<endl;
    //  nobby=GetDaughterCheck(pylis,ind,223,22);
    //cout<<nobby<<endl;
    //if (nobby==1){
    ETOmega+=partE;
    pETOmega+=partE;
    ETAll+=partE;
    PTOmega+=pT;
    pPTOmega+=pT;
    PTAll+=pT;
    cOmega++;
  }
  if (KFID==3122){ //Lambda0
    //if ((mpL>=100)||(mpartD==0)){
    ETLambda0+=partE;
    pETLambda0+=partE;
    ETAll+=partE;
    PTLambda0+=pT;
    pPTLambda0+=pT;
    PTAll+=pT;
    cLambda0++;
  }
  if (KFID==-3122){ //LambdaBar0
    //if ((mpL>=100)||(mpartD==0)){
    ETLambdaBar0+=partE;
    pETLambdaBar0+=partE;
    ETAll+=partE;
    PTLambdaBar0+=pT;
    pPTLambdaBar0+=pT;
    PTAll+=pT;
    cLambdaBar0++;
  }
  if (KFID==3334){ //Omega-
    //if ((mpL>=100)||(mpartD==0)){
    ETOMEGAm+=partE;
    pETOMEGAm+=partE;
    ETAll+=partE;
    PTOMEGAm+=pT;
    pPTOMEGAm+=pT;
    PTAll+=pT;
    PTcharge=+pT;
    cOMEGAm++;
  }
  if (KFID==3322){ //Xi0
    //if ((mpL>=100)||(mpartD==0)){
    ETXi0+=partE;
    pETXi0+=partE;
    ETAll+=partE;
    PTXi0+=pT;
    pPTXi0+=pT;
    PTAll+=pT;
    cXi0++;
  }
  if (KFID==3312){ //Xi-
    //if ((mpL>=100)||(mpartD==0)){
    ETXim+=partE;
    pETXim+=partE;
    ETAll+=partE;
    PTXim+=pT;
    pPTXim+=pT;
    PTAll+=pT;
    PTcharge=+pT;
    cXim++;
  }
  if (KFID==3222){ //Sigma+
    //if ((mpL>=100)||(mpartD==0)){
    ETSigmap+=partE;
    pETSigmap+=partE;
    ETAll+=partE;
    PTSigmap+=pT;
    pPTSigmap+=pT;
    PTAll+=pT;
    PTcharge=+pT;
    cSigmap++;
  }
  if (KFID==3112){ //sigma-
    //if ((mpL>=100)||(mpartD==0)){
    ETSigmam+=partE;
    pETSigmam+=partE;
    ETAll+=partE;
    PTSigmam+=pT;
    pPTSigmam+=pT;
    PTAll+=pT;
    PTcharge=+pT;
    cSigmam++;
  }
  if (KFID==3212){ //sigma0
    //if ((mpL>=100)||(mpartD==0)){
    ETSigma0+=partE;
    pETSigma0+=partE;
    ETAll+=partE;
    PTSigma0+=pT;
    pPTSigma0+=pT;
    PTAll+=pT;
    PTcharge=+pT;
    cSigma0++;
  }
  if (KFID==2212){ //proton
    //if (mpP!=0){ //THIS IS AT MID RAPIDITY so this is not really required
    ETp+=partE;
    pETp+=partE;
    ETAll+=partE;
    ETcharge+=partE;
    PTp+=pT;
    pPTp+=pT;
    PTAll+=pT;
    PTcharge=+pT;
    cp++;
  }
  if (KFID==2112){ //neutron
    ETn+=partE;
    pETn+=partE;
    ETAll+=partE;
    PTn+=pT;
    pPTn+=pT;
    PTAll+=pT;
    cn++;
  }
  if (KFID==-2212){ //antiproton
    ETp_+=partE;
    pETp_+=partE;
    ETAll+=partE;
    ETcharge+=partE;
    PTp_+=pT;
    pPTp_+=pT;
    PTAll+=pT;
    PTcharge=+pT;
    cp_++;
  }
  if (KFID==-2112){ //antineutron
    ETn_+=partE;
    pETn_+=partE;
    ETAll+=partE;
    PTn_+=pT;
    pPTn_+=pT;
    PTAll+=pT;
    cn_++;
  }
  if ((KFID==411)||(KFID==421)||(KFID==431)||(KFID==511)||(KFID==521)||(KFID==531)||(KFID==541)||(KFID==331)||(KFID==441)||(KFID==551)||(KFID==333)||(KFID==443)||(KFID==553)||(KFID==110)||(KFID==4112)||(KFID==4122)||(KFID==4212)||(KFID==4222)||(KFID==4132)||(KFID==4312)||(KFID==4232)||(KFID==4322)||(KFID==4332)||(KFID==5112)||(KFID==5122)||(KFID==5212)||(KFID==5222)||(KFID==3114)||(KFID==3214)||(KFID==3224)||(KFID==3314)||(KFID==3324)||(KFID==4114)||(KFID==4214)||(KFID==4224)||(KFID==4314)||(KFID==4324)||(KFID==4334)||(KFID==5114)||(KFID==5214)||(KFID==5224)) {
    ETother+=partE;
    pETother+=partE;
    ETAll+=partE;
    PTother+=partE;
    pPTother+=partE;
    PTAll+=partE;
    cother++;
    //out<<event<<"\t"<<KFID<<";\t";
  }
  }//FULL INCLUSION
  if (DMODE=='b'){ //ONLY IF FINAL
/****************************************************
The following section counts particles based on the below conditions to ensure no double counting

Particle:        inclusion requirements:
pion+-           ALL INCLUDED
pion0            ALL INCLUDED
kaon+            Only if final
kaon-            Only if final
kaon0_L          Only if final
kaon0_S          Only if final
eta              Only if final
omega            Only if final
Lambda0          Only if final
LambdaBar0       Only if final
Omega-           Only if final
Xi0              Only if final
Xi-              Only if final
Sigma+           Only if final
Sigma0           Only if final
Sigma-           Only if final
proton           ALL INCLUDED
neutron          ALL INCLUDED
other            Only if final

****************************************************/

  if (KFID==211){ //pi+ all included
    ETpip+=partE; //repeated for event
    pETpip+=partE; //repeated over total run
    ETAll+=partE;
    ETcharge+=partE;
    PTpip+=pT;
    pPTpip+=pT;
    PTAll+=pT;
    PTcharge=+pT;
    cpip++;  //repeated over total run
  }
  if (KFID==-211){ //pi- all included
    ETpim+=partE;
    pETpim+=partE;
    ETAll+=partE;
    ETcharge+=partE;
    PTpim+=pT;
    pPTpim+=pT;
    PTAll+=pT;
    PTcharge=+pT;
    cpim++;
  }
  if (KFID==111){ //pi0
    ETpi0+=partE;
    pETpi0+=partE;
    ETAll+=partE;
    PTpi0+=pT;
    pPTpi0+=pT;
    PTAll+=pT;
    cpi0++;
  }
  if (KFID==321){ //K+ //all included
    if (mpartD==0){
      ETKp+=partE;
      pETKp+=partE;
      ETAll+=partE;
      ETcharge+=partE;
      PTKp+=pT;
      pPTKp+=pT;
      PTAll+=pT;
      PTcharge=+pT;
      cKp++;
    }
  }
  if (KFID==-321){ //K- //all included
    if (mpartD==0){
      ETKm+=partE;
      pETKm+=partE;
      ETAll+=partE;
      ETcharge+=partE;
      PTKm+=pT;
      pPTKm+=pT;
      PTAll+=pT;
      PTcharge=+pT;
      cKm++;
  }}
  if (KFID==130){// KL
    //if ((mpL>=100)||(mpartD==0)){//counted if final particle or lifetime > 1cm
    if (mpartD==0){
      ETKL+=partE;
      pETKL+=partE;
      ETAll+=partE;
      PTKL+=pT;
      pPTKL+=pT;
      PTAll+=pT;
      cKL++;
    }
  }
  if (KFID==310){ //KS
    //if ((mpL>=100)||(mpartD==0)){//counted if final particle or lifetime > 1cm
    if (mpartD==0){
      ETKS+=partE;
      pETKS+=partE;
      ETAll+=partE;
      PTKS+=pT;
      pPTKS+=pT;
      PTAll+=pT;
      cKS++;
    }
  }
  if (KFID==221){ //eta
    //cout<<ind<<"= particle number (check eta) running func:"<<endl;
    //nobby=GetDaughterCheck(pylis,ind,221,22);
    //cout<<nobby<<endl;
    //  if (nobby==1){
    if (mpartD==0){
      ETEta+=partE;
      pETEta+=partE;
      ETAll+=partE;
      PTEta+=pT;
      pPTEta+=pT;
      PTAll+=pT;
      cEta++;
  }}
  if (KFID==223){ // omega
    //cout<<ind<<"= particle number (check omega) running func:"<<endl;
    //  nobby=GetDaughterCheck(pylis,ind,223,22);
    //cout<<nobby<<endl;
    //if (nobby==1){
    if (mpartD==0){
      ETOmega+=partE;
      pETOmega+=partE;
      ETAll+=partE;
      PTOmega+=pT;
      pPTOmega+=pT;
      PTAll+=pT;
      cOmega++;
  }}
  if (KFID==3122){ //Lambda0
    //if ((mpL>=100)||(mpartD==0)){
    if (mpartD==0){
      ETLambda0+=partE;
      pETLambda0+=partE;
      ETAll+=partE;
      PTLambda0+=pT;
      pPTLambda0+=pT;
      PTAll+=pT;
      cLambda0++;
  }}
  if (KFID==-3122){ //LambdaBar0
    //if ((mpL>=100)||(mpartD==0)){
    if (mpartD==0){
      ETLambdaBar0+=partE;
      pETLambdaBar0+=partE;
      ETAll+=partE;
      PTLambdaBar0+=pT;
      pPTLambdaBar0+=pT;
      PTAll+=pT;
      cLambdaBar0++;
  }}
  if (KFID==3334){ //Omega-
    //if ((mpL>=100)||(mpartD==0)){
    if (mpartD==0){
      ETOMEGAm+=partE;
      pETOMEGAm+=partE;
      ETAll+=partE;
      PTOMEGAm+=pT;
      pPTOMEGAm+=pT;
      PTAll+=pT;
      PTcharge=+pT;
      cOMEGAm++;
  }}
  if (KFID==3322){ //Xi0
    //if ((mpL>=100)||(mpartD==0)){
    if (mpartD==0){
      ETXi0+=partE;
      pETXi0+=partE;
      ETAll+=partE;
      PTXi0+=pT;
      pPTXi0+=pT;
      PTAll+=pT;
      cXi0++;
  }}
  if (KFID==3312){ //Xi-
    //if ((mpL>=100)||(mpartD==0)){
    if (mpartD==0){
      ETXim+=partE;
      pETXim+=partE;
      ETAll+=partE;
      PTXim+=pT;
      pPTXim+=pT;
      PTAll+=pT;
      PTcharge=+pT;
      cXim++;
  }}
  if (KFID==3222){ //Sigma+
    //if ((mpL>=100)||(mpartD==0)){
    if (mpartD==0){
      ETSigmap+=partE;
      pETSigmap+=partE;
      ETAll+=partE;
      PTSigmap+=pT;
      pPTSigmap+=pT;
      PTAll+=pT;
      PTcharge=+pT;
      cSigmap++;
  }}
  if (KFID==3112){ //sigma-
    //if ((mpL>=100)||(mpartD==0)){
    if (mpartD==0){
      ETSigmam+=partE;
      pETSigmam+=partE;
      ETAll+=partE;
      PTSigmam+=pT;
      pPTSigmam+=pT;
      PTAll+=pT;
      PTcharge=+pT;
      cSigmam++;
  }}
  if (KFID==3212){ //sigma0
    //if ((mpL>=100)||(mpartD==0)){
    if (mpartD==0){
      ETSigma0+=partE;
      pETSigma0+=partE;
      ETAll+=partE;
      PTSigma0+=pT;
      pPTSigma0+=pT;
      PTAll+=pT;
      PTcharge=+pT;
      cSigma0++;
  }}
  if (KFID==2212){ //proton
    //if (mpP!=0){ //THIS IS AT MID RAPIDITY so this is not really required
    if (mpartD==0){
      ETp+=partE;
      pETp+=partE;
      ETAll+=partE;
      ETcharge+=partE;
      PTp+=pT;
      pPTp+=pT;
      PTAll+=pT;
      PTcharge=+pT;
      cp++;
  }}
  if (KFID==2112){ //neutron
    if (mpartD==0){
      ETn+=partE;
      pETn+=partE;
      ETAll+=partE;
      PTn+=pT;
      pPTn+=pT;
      PTAll+=pT;
      cn++;
  }}
  if (KFID==-2212){ //antiproton
    if (mpartD==0){
      ETp_+=partE;
      pETp_+=partE;
      ETAll+=partE;
      ETcharge+=partE;
      PTp_+=pT;
      pPTp_+=pT;
      PTAll+=pT;
      PTcharge=+pT;
      cp_++;
  }}
  if (KFID==-2112){ //antineutron
    if (mpartD==0){
      ETn_+=partE;
      pETn_+=partE;
      ETAll+=partE;
      PTn_+=pT;
      pPTn_+=pT;
      PTAll+=pT;
      cn_++;
  }}
  if ((KFID==411)||(KFID==421)||(KFID==431)||(KFID==511)||(KFID==521)||(KFID==531)||(KFID==541)||(KFID==331)||(KFID==441)||(KFID==551)||(KFID==333)||(KFID==443)||(KFID==553)||(KFID==110)||(KFID==4112)||(KFID==4122)||(KFID==4212)||(KFID==4222)||(KFID==4132)||(KFID==4312)||(KFID==4232)||(KFID==4322)||(KFID==4332)||(KFID==5112)||(KFID==5122)||(KFID==5212)||(KFID==5222)||(KFID==3114)||(KFID==3214)||(KFID==3224)||(KFID==3314)||(KFID==3324)||(KFID==4114)||(KFID==4214)||(KFID==4224)||(KFID==4314)||(KFID==4324)||(KFID==4334)||(KFID==5114)||(KFID==5214)||(KFID==5224)) {
    if (mpartD==0){
      ETother+=partE;
      pETother+=partE;
      ETAll+=partE;
      PTother+=partE;
      pPTother+=partE;
      PTAll+=partE;
      cother++;
      out<<event<<"\t"<<KFID<<endl;
  }}
  }//Least Inlcusion
  if (DMODE=='c'){ //Decay version 1
/****************************************************
The following section counts particles based on the below conditions to ensure no double counting
FOR ET TOTAL:

Particle:        inclusion requirements:
pion+-           All but KOS, [anti]Lambda daughters
pion0            Most Included (if omega parent is 2nd decay mode or KOS/lambda/antilam daughter, pi0 is not included)
kaon+            ALL INCLUDED
kaon-            ALL INCLUDED
kaon0_L          ALL INCLUDED ct=15.34m
kaon0_S          ALL INCLUDED ct=2.6844cm
eta              Included only if it decays to gamma t=5.02E-19s
omega            Included only if it decays to gamma t=7.75E-23s
Lambda0          ALL INCLUDED ct=7.89cm
LambdaBar0       ALL INCLUDED
Omega-           ALL INCLUDED ct=2.461cm
Xi0              ALL INCLUDED ct=8.71cm
Xi-              ALL INCLUDED ct=4.91cm
Sigma+           ALL INCLUDED ct=2.404cm
Sigma0           ALL INCLUDED ct=2.22x10^(-11)m
Sigma-           ALL INCLUDED ct=4.424cm
proton           All except from [anti]lambda decays
antiproton       All except from [anti]lambda decays
neutron          All except from [anti]lambda decays
Antineutron      All except from [anti]lambda decays
exotics          Final State Particles Only

****************************************************/
  if (KFID==211){ //pi+
    XEM=pylis[mpP].KFpart;
    if ((XEM!=310)&&(XEM!=3122)&&(XEM!=-3122)){ //do not count if daughter of KOS
      ETpip+=partE; //repeated for event
      pETpip+=partE; //repeated over total run
      ETAll+=partE;
      ETcharge+=partE;
      PTpip+=pT;
      pPTpip+=pT;
      PTAll+=pT;
      PTcharge=+pT;
      cpip++;  //repeated over total run
  }}
  if (KFID==-211){ //pi-
    XEM=pylis[mpP].KFpart;
    if ((XEM!=310)&&(XEM!=3122)&&(XEM!=-3122)){ //do not count if daughter of KOS
      ETpim+=partE;
      pETpim+=partE;
      ETAll+=partE;
      ETcharge+=partE;
      PTpim+=pT;
      pPTpim+=pT;
      PTAll+=pT;
      PTcharge=+pT;
      cpim++;
  }}
  if (KFID==111){ //pi0
    XEM=pylis[mpP].KFpart; // is one of parents daughters a gamma? if so not included (this is to exclude certain omega and eta modes of decay)
    nobby=0;
    if (XEM==223){
      nobby=GetDaughterCheck(pylis,mpP,111,22);
    }
    if (nobby==0){
      if ((XEM!=310)&&(XEM!=3122)&&(XEM!=-3122)){ //do not count if daughter of KOS
        ETpi0+=partE;
        pETpi0+=partE;
        ETAll+=partE;
        PTpi0+=pT;
        pPTpi0+=pT;
        PTAll+=pT;
        cpi0++;
  }}}
  if (KFID==321){ //K+ //all included
    ETKp+=partE;
    pETKp+=partE;
    ETAll+=partE;
    ETcharge+=partE;
    PTKp+=pT;
    pPTKp+=pT;
    PTAll+=pT;
    PTcharge=+pT;
    cKp++;
  }
  if (KFID==-321){ //K- //all included
    ETKm+=partE;
    pETKm+=partE;
    ETAll+=partE;
    ETcharge+=partE;
    PTKm+=pT;
    pPTKm+=pT;
    PTAll+=pT;
    PTcharge=+pT;
    cKm++;
  }
  if (KFID==130){// KL
    ETKL+=partE;
    pETKL+=partE;
    ETAll+=partE;
    PTKL+=pT;
    pPTKL+=pT;
    PTAll+=pT;
    cKL++;
  }
  if (KFID==310){ //KS
    ETKS+=partE;
    pETKS+=partE;
    ETAll+=partE;
    PTKS+=pT;
    pPTKS+=pT;
    PTAll+=pT;
    cKS++;
  }
  if (KFID==221){ //eta
    nobby=GetDaughterCheck(pylis,ind,221,22);
    cETAall++;
    pETETAall+=partE;
    pPTETAall+=pT;
    if (nobby==1){
      ETEta+=partE;
      pETEta+=partE;
      ETAll+=partE;
      PTEta+=pT;
      pPTEta+=pT;
      PTAll+=pT;
      cEta++;
  }}
  if (KFID==223){ // omega
    nobby=GetDaughterCheck(pylis,ind,223,22);
    cOmegaall++;
    pETOmegaall+=partE;
    pPTOmegaall+=pT;
    if (nobby==1){
      ETOmega+=partE;
      pETOmega+=partE;
      ETAll+=partE;
      PTOmega+=pT;
      pPTOmega+=pT;
      PTAll+=pT;
      cOmega++;
  }}
  if (KFID==3122){ //Lambda0
    //if ((mpL>=100)||(mpartD==0)){
    ETLambda0+=partE;
    pETLambda0+=partE;
    ETAll+=partE;
    PTLambda0+=pT;
    pPTLambda0+=pT;
    PTAll+=pT;
    cLambda0++;
  }
  if (KFID==-3122){ //LambdaBar0
    //if ((mpL>=100)||(mpartD==0)){
    ETLambdaBar0+=partE;
    pETLambdaBar0+=partE;
    ETAll+=partE;
    PTLambdaBar0+=pT;
    pPTLambdaBar0+=pT;
    PTAll+=pT;
    cLambdaBar0++;
  }
  if (KFID==3334){ //Omega-
    //if ((mpL>=100)||(mpartD==0)){
    ETOMEGAm+=partE;
    pETOMEGAm+=partE;
    ETAll+=partE;
    PTOMEGAm+=pT;
    pPTOMEGAm+=pT;
    PTAll+=pT;
    PTcharge=+pT;
    cOMEGAm++;
  }
  if (KFID==3322){ //Xi0
    //if ((mpL>=100)||(mpartD==0)){
    ETXi0+=partE;
    pETXi0+=partE;
    ETAll+=partE;
    PTXi0+=pT;
    pPTXi0+=pT;
    PTAll+=pT;
    cXi0++;
  }
  if (KFID==3312){ //Xi-
    //if ((mpL>=100)||(mpartD==0)){
    ETXim+=partE;
    pETXim+=partE;
    ETAll+=partE;
    PTXim+=pT;
    pPTXim+=pT;
    PTAll+=pT;
    PTcharge=+pT;
    cXim++;
  }
  if (KFID==3222){ //Sigma+
    if ((mpL>=100)||(mpartD==0)){
      ETSigmap+=partE;
      pETSigmap+=partE;
      ETAll+=partE;
      PTSigmap+=pT;
      pPTSigmap+=pT;
      PTAll+=pT;
      PTcharge=+pT;
      cSigmap++;
  }}
  if (KFID==3112){ //sigma-
    if ((mpL>=100)||(mpartD==0)){
      ETSigmam+=partE;
      pETSigmam+=partE;
      ETAll+=partE;
      PTSigmam+=pT;
      pPTSigmam+=pT;
      PTAll+=pT;
      PTcharge=+pT;
      cSigmam++;
  }}
  if (KFID==3212){ //sigma0
    if ((mpL>=100)||(mpartD==0)){
      ETSigma0+=partE;
      pETSigma0+=partE;
      ETAll+=partE;
      PTSigma0+=pT;
      pPTSigma0+=pT;
      PTAll+=pT;
      PTcharge=+pT;
      cSigma0++;
  }}
  if (KFID==2212){ //proton
    XEM=pylis[mpP].KFpart;
    if ((XEM!=3122)&&(XEM!=-3122)){
      ETp+=partE;
      pETp+=partE;
      ETAll+=partE;
      ETcharge+=partE;
      PTp+=pT;
      pPTp+=pT;
      PTAll+=pT;
      PTcharge=+pT;
      cp++;
  }}
  if (KFID==2112){ //neutron
    XEM=pylis[mpP].KFpart;
    if ((XEM!=3122)&&(XEM!=-3122)){
      ETn+=partE;
      pETn+=partE;
      ETAll+=partE;
      PTn+=pT;
      pPTn+=pT;
      PTAll+=pT;
      cn++;
  }}
  if (KFID==-2212){ //antiproton
    XEM=pylis[mpP].KFpart;
    if ((XEM!=3122)&&(XEM!=-3122)){
      ETp_+=partE;
      pETp_+=partE;
      ETAll+=partE;
      ETcharge+=partE;
      PTp_+=pT;
      pPTp_+=pT;
      PTAll+=pT;
      PTcharge=+pT;
      cp_++;
  }}
  if (KFID==-2112){ //antineutron
    XEM=pylis[mpP].KFpart;
    if ((XEM!=3122)&&(XEM!=-3122)){
      ETn_+=partE;
      pETn_+=partE;
      ETAll+=partE;
      PTn_+=pT;
      pPTn_+=pT;
      PTAll+=pT;
      cn_++;
  }}
  if ((KFID==411)||(KFID==421)||(KFID==431)||(KFID==511)||(KFID==521)||(KFID==531)||(KFID==541)||(KFID==331)||(KFID==441)||(KFID==551)||(KFID==333)||(KFID==443)||(KFID==553)||(KFID==110)||(KFID==4112)||(KFID==4122)||(KFID==4212)||(KFID==4222)||(KFID==4132)||(KFID==4312)||(KFID==4232)||(KFID==4322)||(KFID==4332)||(KFID==5112)||(KFID==5122)||(KFID==5212)||(KFID==5222)||(KFID==3114)||(KFID==3214)||(KFID==3224)||(KFID==3314)||(KFID==3324)||(KFID==4114)||(KFID==4214)||(KFID==4224)||(KFID==4314)||(KFID==4324)||(KFID==4334)||(KFID==5114)||(KFID==5214)||(KFID==5224)) {
    if ((mpartD==22)||mpartD2==22){
      ETother+=partE;
      pETother+=partE;
      ETAll+=partE;
      PTother+=partE;
      pPTother+=partE;
      PTAll+=partE;
      cother++;
      //out<<event<<"\t"<<KFID<<endl;
  }}
  }//FULL INCLUSION
  } //end pseudorapidity cut
  }//end particle loop
  /******************************************************************************
  Histogarms are filled here, the histograms are of the format:
  hET<particleName> xaxis is ET particle/ Total ET
                    yaxis is number of events at that fraction
  *******************************************************************************/
        Double_t omegaET=0;
        Double_t pi0ET=0;
        Double_t ChAllRat=0;
        hETAll->Fill(ETAll);
        hETcharge->Fill(ETcharge);
        ChAllRat=ETcharge/ETAll;
        hETChOverAll->Fill(ETcharge,ChAllRat);
        //cout<<"TOTAL: "<<ETAll<<endl;
        if (ETpip!=0){
        //pETpip=ETpip/ETAll;
        hETpip->Fill(ETpip);
        hPTpip->Fill(PTpip);

          //cout<<"pip "<<ETpip<<endl;
        }
        if (ETpim!=0){
        //pETpim=ETpim/ETAll;
        hETpim->Fill(ETpim);
        hPTpim->Fill(PTpim);
          //cout<<"pim "<<ETpim<<endl;
        }
        if (ETpi0!=0){
        //pETpi0=ETpi0/ETAll;
        hETpi0->Fill(ETpi0);
        hPTpi0->Fill(PTpi0);
          //cout<<"pi0 "<<ETpi0<<endl;
        }
        if (ETKp!=0){
        //pETKp=ETKp/ETAll;
        hETKp->Fill(ETKp);
        hPTKp->Fill(PTKp);
          //cout<<"Kp "<<ETKp<<endl;
        }
        if (ETKm!=0){
        //pETKm=ETKm/ETAll;
        hETKm->Fill(ETKm);
        hPTKm->Fill(PTKm);
          //cout<<"Kp "<<ETKp<<endl;
        }
        if (ETKL!=0){
        //pETKL=ETKL/ETAll;
        hETKL->Fill(ETKL);
        hPTKL->Fill(PTKL);
          //cout<<"KL "<<ETKL<<endl;
        }
        if (ETKS!=0){
        //pETKS=ETKS/ETAll;
        hETKS->Fill(ETKS);
        hPTKS->Fill(PTKS);
          //cout<<"KS "<<ETKS<<endl;
        }
        if (ETEta!=0){
        //pETEta=ETEta/ETAll;
        hETEta->Fill(ETEta);
        hPTEta->Fill(PTEta);
          //cout<<"Eta "<<ETEta<<endl;
        }
        if (ETOmega!=0){
        //pETOmega=ETOmega/ETAll;
        hETOmega->Fill(ETOmega);
        hPTOmega->Fill(PTOmega);
          //cout<<"omega "<<ETOmega<<endl;
        }
        if (ETLambda0!=0){
        //pETLambda0=ETLambda0/ETAll;
        hETLambda0->Fill(ETLambda0);
        hPTLambda0->Fill(PTLambda0);
          //cout<<"Lambda0 "<<ETLambda0<<endl;
        }
        if (ETLambdaBar0!=0){
        //pETLambdaBar0=ETLambdaBar0/ETAll;
        hETLambdaBar0->Fill(ETLambdaBar0);
        hPTLambdaBar0->Fill(PTLambdaBar0);
          //cout<<"Lambda0 "<<ETLambda0<<endl;
        }
        if (ETOMEGAm!=0){
        //pETOMEGAm=ETOMEGAm/ETAll;
        hETOMEGAm->Fill(ETOMEGAm);
        hPTOMEGAm->Fill(PTOMEGAm);
          //cout<<"Omega "<<ETOMEGAm<<endl;
        }
        if (ETXi0!=0){
        //pETXi0=ETXi0/ETAll;
        hETXi0->Fill(ETXi0);
        hPTXi0->Fill(PTXi0);
          //cout<<"Xi0 "<<ETXi0<<endl;
        }
        if (ETXim!=0){
        //pETXim=ETXim/ETAll;
        hETXim->Fill(ETXim);
        hPTXim->Fill(PTXim);
          //cout<<"Xim "<<ETXim<<endl;
        }
        if (ETSigmap!=0){
        //pETSigmap=ETSigmap/ETAll;
        hETSigmap->Fill(ETSigmap);
        hPTSigmap->Fill(PTSigmap);
          //cout<<"Sigma+ "<<ETSigmap<<endl;
        }
        if (ETSigmam!=0){
        //pETSigmam=ETSigmam/ETAll;
        hETSigmam->Fill(ETSigmam);
        hPTSigmam->Fill(PTSigmam);
          //cout<<"Sigma- "<<ETSigmam<<endl;
        }
        if (ETSigma0!=0){
        //pETSigma0=ETSigma0/ETAll;
        hETSigma0->Fill(ETSigma0);
        hPTSigma0->Fill(PTSigma0);
          //cout<<"Sigma0 "<<ETSigma0<<endl;
        }
        if (ETp!=0){
        //pETp=ETp/ETAll;
        hETp->Fill(ETp);
        hPTp->Fill(PTp);
          //cout<<"p "<<ETp<<endl;
        }
        if (ETn!=0){
        //pETn=ETn/ETAll;
        hETn->Fill(ETn);
        hPTn->Fill(PTn);
          //cout<<"n "<<ETn<<endl;
        }
        if (ETp_!=0){
        //pETp_=ETp_/ETAll;
        hETp_->Fill(ETp_);
        hPTp_->Fill(PTp_);
          //cout<<"p "<<ETp<<endl;
        }
        if (ETn_!=0){
        //pETn_=ETn_/ETAll;
        hETn_->Fill(ETn_);
        hPTn_->Fill(PTn_);
          //cout<<"n "<<ETn<<endl;
        }
        if (ETother!=0){
        //pETother=ETother/ETAll;
        hETother->Fill(ETother);
        hPTother->Fill(PTother);
          //cout<<"other "<<ETother<<endl;
        }
    // First sum up transverse energy for centrality measure and also
    // check that the trigger requiring ar least one charged particle
    // forward and backward.
    /*
    double etfwd = 0.0;
    bool trigfwd = false;
    bool trigbwd = false;
    int nc = 0;
    for (int i = 0; i < pythia.event.size(); ++i) {
      Particle & p = pythia.event[i];
      if ( p.isFinal() ) {
        double eta = p.eta();
        if ( p.isCharged() && p.pT() > 0.1 && eta < -2.09 && eta > -3.84 )
          trigfwd = true;
        if ( p.isCharged() && p.pT() > 0.1 && eta > 2.09 && eta < 3.84 )
          trigbwd = true;
        if ( p.pT() > 0.1 && abs(eta) > 3.2 && abs(eta) < 4.9 )
          etfwd += p.eT();
        if ( p.isCharged() && p.pT() > 0.1 && abs(eta) < 0.5 ) ++nc;
      }
    }
    // Skip if not triggered
    if ( !(trigfwd && trigbwd) ) continue;

    // Keep track of the sum of waights
    double weight = pythia.info.weight();
    sumw += weight;

    // Histogram and save the summed Et.
    sumet.fill(etfwd, weight);
    gencent.insert(make_pair(etfwd, weight));

    // Also fill the number of (absorptively and diffractively)
    // wounded nucleaons.
    int nw = pythia.info.hiinfo->nAbsTarg() +
             pythia.info.hiinfo->nDiffTarg() +
             pythia.info.hiinfo->nAbsProj() +
             pythia.info.hiinfo->nDiffProj();
    wounded.fill(nw, weight);

    // Find the correct centrality histograms.
    MapIdx::iterator genit = genetaidx.upper_bound(etfwd);
    int genidx = genit== genetaidx.end()? -1: genit->second;

    // Sum the weights in the centrality classes, skip if not in a class.
    if ( genidx < 0 ) continue;
    gensumw[genidx] += weight;
    hmult[genidx]->fill(nc, weight);
    lmult[genidx]->fill(nc, weight);
    gensumn[genidx] += 1.0;
    cmult[genidx] += nc*weight;
    cmult2[genidx] += nc*nc*weight;
    wound[genidx] += nw*weight;
    wound2[genidx] += nw*nw*weight;

    // Go through the event again and fill the eta distributions.
    for (int i = 0; i < pythia.event.size(); ++i) {
      Particle & p = pythia.event[i];
      if ( p.isFinal() && p.isCharged() &&
           abs(p.eta()) < 2.7 && p.pT() > 0.1 ) {
         etadist[genidx]->fill(p.eta(), weight);
      }
    }
    */
  }//end event loop

    hNEvents->Write();

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
    hETOMEGAm->Write();
    hETXi0->Write();
    hETXim->Write();
    hETSigmap->Write();
    hETSigmam->Write();
    hETSigma0->Write();
    hETp->Write();
    hETn->Write();
    hETp_->Write();
    hETn_->Write();
    hETother->Write();

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
    hPTOMEGAm->Write();
    hPTXi0->Write();
    hPTXim->Write();
    hPTSigmap->Write();
    hPTSigmam->Write();
    hPTSigma0->Write();
    hPTp->Write();
    hPTn->Write();
    hPTp_->Write();
    hPTn_->Write();
    hPTother->Write();

    double TET=pETpip+pETpim+pETpi0+pETKp+pETKm+pETKL+pETKS+pETLambda0+pETLambdaBar0+pETp+pETn+pETp_+pETn_+pETEta+pETOmega+pETOMEGAm+pETXi0+pETXim+pETSigmap+pETSigmam+pETSigma0;
    double TPT=pPTpip+pPTpim+pPTpi0+pPTKp+pPTKm+pPTKL+pPTKS+pPTLambda0+pPTLambdaBar0+pPTp+pPTn+pPTp_+pPTn_+pPTEta+pPTOmega+pPTOMEGAm+pPTXi0+pPTXim+pPTSigmap+pPTSigmam+pPTSigma0;
    out<<"Total Transverse Energy: "<<TET<<"\tTotal Transverse Momentum: "<<TPT<<"\n\n";

    out<<"Particle\tNumber\tEt\tPt\t%EtVsTotal\t%PtVsTotal\n";
    out<<"Pi+:\t\t"<<cpip<<"\t"<<pETpip<<"\t"<<pPTpip<<"\t"<<(pETpip/TET)*100<<"\t"<<(pPTpip/TPT)*100<<endl;
    out<<"Pi-:\t\t"<<cpim<<"\t"<<pETpim<<"\t"<<pPTpim<<"\t"<<(pETpim/TET)*100<<"\t"<<(pPTpim/TPT)*100<<endl;
    out<<"Pi0:\t\t"<<cpi0<<"\t"<<pETpi0<<"\t"<<pPTpi0<<"\t"<<(pETpi0/TET)*100<<"\t"<<(pPTpi0/TPT)*100<<endl;
    out<<"K+:\t\t"<<cKp<<"\t"<<pETKp<<"\t"<<pPTKp<<"\t"<<(pETKp/TET)*100<<"\t"<<(pPTKp/TPT)*100<<endl;
    out<<"K-:\t\t"<<cKm<<"\t"<<pETKm<<"\t"<<pPTKm<<"\t"<<(pETKm/TET)*100<<"\t"<<(pPTKm/TPT)*100<<endl;
    out<<"K0L:\t\t"<<cKL<<"\t"<<pETKL<<"\t"<<pPTKL<<"\t"<<(pETKL/TET)*100<<"\t"<<(pPTKL/TPT)*100<<endl;
    out<<"K0S:\t\t"<<cKS<<"\t"<<pETKS<<"\t"<<pPTKS<<"\t"<<(pETKS/TET)*100<<"\t"<<(pPTKS/TPT)*100<<endl;
    out<<"Lambda0:\t"<<cLambda0<<"\t"<<pETLambda0<<"\t"<<pPTLambda0<<"\t"<<(pETLambda0/TET)*100<<"\t"<<(pPTLambda0/TPT)*100<<endl;
    out<<"LambdaBar0:\t"<<cLambdaBar0<<"\t"<<pETLambdaBar0<<"\t"<<pPTLambdaBar0<<"\t"<<(pETLambdaBar0/TET)*100<<"\t"<<(pPTLambdaBar0/TPT)*100<<endl;
    out<<"Proton:\t"<<cp<<"\t"<<pETp<<"\t"<<pPTp<<"\t"<<(pETp/TET)*100<<"\t"<<(pPTp/TPT)*100<<endl;
    out<<"Antiproton:\t"<<cp_<<"\t"<<pETp_<<"\t"<<pPTp_<<"\t"<<(pETp_/TET)*100<<"\t"<<(pPTp_/TPT)*100<<endl;
    out<<"Neutron:\t"<<cn<<"\t"<<pETn<<"\t"<<pPTn<<"\t"<<(pETn/TET)*100<<"\t"<<(pPTn/TPT)*100<<endl;
    out<<"Antineutron:\t"<<cn_<<"\t"<<pETn_<<"\t"<<pPTn_<<"\t"<<(pETn_/TET)*100<<"\t"<<(pPTn_/TPT)*100<<endl;
    out<<"eta:\t"<<cEta<<"\t"<<pETEta<<"\t"<<pPTEta<<"\t"<<(pETEta/TET)*100<<"\t"<<(pPTEta/TPT)*100<<endl;
    out<<"[eta]:\t"<<cETAall<<"\t"<<pETETAall<<"\t"<<pPTETAall<<"\t"<<(pETETAall/TET)*100<<"\t"<<(pPTETAall/TPT)*100<<endl;
    out<<"omega:\t"<<cOmega<<"\t"<<pETOmega<<"\t"<<pPTOmega<<"\t"<<(pETOmega/TET)*100<<"\t"<<(pPTOmega/TPT)*100<<endl;
    out<<"[omega]:\t"<<cOmegaall<<"\t"<<pETOmegaall<<"\t"<<pPTOmegaall<<"\t"<<(pETOmegaall/TET)*100<<"\t"<<(pPTOmegaall/TPT)*100<<endl;
    out<<"Omega-:\t"<<cOMEGAm<<"\t"<<pETOMEGAm<<"\t"<<pPTOMEGAm<<"\t"<<(pETOMEGAm/TET)*100<<"\t"<<(pPTOMEGAm/TPT)*100<<endl;
    out<<"Xi0:\t"<<cXi0<<"\t"<<pETXi0<<"\t"<<pPTXi0<<"\t"<<(pETXi0/TET)*100<<"\t"<<(pPTXi0/TPT)*100<<endl;
    out<<"Xi-:\t"<<cXim<<"\t"<<pETXim<<"\t"<<pPTXim<<"\t"<<(pETXim/TET)*100<<"\t"<<(pPTXim/TPT)*100<<endl;
    out<<"Sigma+:\t"<<cSigmap<<"\t"<<pETSigmap<<"\t"<<pPTSigmap<<"\t"<<(pETSigmap/TET)*100<<"\t"<<(pPTSigmap/TPT)*100<<endl;
    out<<"Sigma-:\t"<<cSigmam<<"\t"<<pETSigmam<<"\t"<<pPTSigmam<<"\t"<<(pETSigmam/TET)*100<<"\t"<<(pPTSigmam/TPT)*100<<endl;
    out<<"Sigma0:\t"<<cSigma0<<"\t"<<pETSigma0<<"\t"<<pPTSigma0<<"\t"<<(pETSigma0/TET)*100<<"\t"<<(pPTSigma0/TPT)*100<<endl;
    out<<endl;
    out<<"[Pions]:\t"<<(cpip+cpim+cpi0)<<"\t"<<(pETpip+pETpim+pETpi0)<<"\t"<<(pPTpip+pPTpim+pPTpi0)<<"\t"<<((pETpip+pETpim+pETpi0)/TET)*100<<"\t"<<((pPTpip+pPTpim+pPTpi0)/TPT)*100<<endl;
    out<<"[Kaons]:\t"<<(cKp+cKm+cKL+cKS)<<"\t"<<(pETKp+pETKm+pETKL+pETKS)<<"\t"<<(pPTKp+pPTKm+pPTKL+pPTKS)<<"\t"<<((pETKp+pETKm+pETKL+pETKS)/TET)*100<<"\t"<<((pPTKp+pPTKm+pPTKL+pPTKS)/TPT)*100<<endl;
    out<<"[Lambdas]:\t"<<(cLambda0+cLambdaBar0)<<"\t"<<(pETLambda0+pETLambdaBar0)<<"\t"<<(pPTLambda0+pPTLambdaBar0)<<"\t"<<((pETLambda0+pETLambdaBar0)/TET)*100<<"\t"<<((pPTLambda0+pPTLambdaBar0)/TPT)*100<<endl;
    out<<"[p,pbar]:\t"<<(cp+cp_)<<"\t"<<(pETp+pETp_)<<"\t"<<(pPTp+pPTp_)<<"\t"<<((pETp+pETp_)/TET)*100<<"\t"<<((pPTp+pPTp_)/TPT)*100<<endl;
    out<<"[n,nbar]:\t"<<(cn+cn_)<<"\t"<<(pETn+pETn_)<<"\t"<<(pPTn+pPTn_)<<"\t"<<((pETn+pETn_)/TET)*100<<"\t"<<((pPTn+pPTn_)/TPT)*100<<endl;
    out<<"[eta,ome]:\t"<<(cEta+cOmega)<<"\t"<<(pETEta+pETOmega)<<"\t"<<(pPTEta+pPTOmega)<<"\t"<<((pETEta+pETOmega)/TET)*100<<"\t"<<((pPTEta+pPTOmega)/TPT)*100<<endl;
    out<<"[exotics]:\t"<<(cOMEGAm+cXi0+cXim+cSigmap+cSigmam+cSigma0)<<"\t"<<(pETOMEGAm+pETXi0+pETXim+pETSigmap+pETSigmam+pETSigma0)<<"\t"<<(pPTOMEGAm+pPTXi0+pPTXim+pPTSigmap+pPTSigmam+pPTSigma0)<<"\t"<<((pETOMEGAm+pETXi0+pETXim+pETSigmap+pETSigmam+pETSigma0)/TET)*100<<"\t"<<((pPTOMEGAm+pPTXi0+pPTXim+pPTSigmap+pPTSigmam+pPTSigma0)/TPT)*100<<endl;
    out<<"[other]:\t"<<cother<<"\t"<<pETother<<"\t"<<pPTother<<"\t"<<(pETother/TET)*100<<"\t"<<(pPTother/TPT)*100<<endl;

    out.close();
    outfile->Close();
    cout<<"I made it here line 1716"<<endl;
  // The run is over, so we write out some statistics.


  // Now, we just have to normalize and prtint out the histograms. We
  // choose to print the histograms to a file that can be read by
  // eg. gnuplot.
  char* nameOfile= Form("%i%fPAngatyr.dat",jobID,SNN);
  ofstream ofs(nameOfile);

  sumet /= sumw*2.0;
  ofs << "# " << sumet.getTitle() << endl;
  sumet.table(ofs);

  wounded /= sumw*2.0;
  ofs << "\n# " << wounded.getTitle() << endl;
  wounded.table(ofs);

  // Print out the centrality binned eta distributions and delete the
  // heap-allocate histograms.
  for ( int idx = 0; idx < 9; ++idx ) {
    *hmult[idx] /= gensumw[idx]*40.0;
    ofs << "\n# " << hmult[idx]->getTitle() << endl;
    hmult[idx]->table(ofs);
    delete hmult[idx];
    *lmult[idx] /= gensumw[idx]*4.0;
    ofs << "\n# " << lmult[idx]->getTitle() << endl;
    lmult[idx]->table(ofs);
    delete lmult[idx];
    *etadist[idx] /= gensumw[idx]*0.1;
    ofs << "\n# " << etadist[idx]->getTitle() << endl;
    etadist[idx]->table(ofs);
    delete etadist[idx];
  }

  // Print out average central charged multiplicity as a function of
  // centrality.
  ofs << "\n# Nch0\n";
  for ( int idx = 0; idx < 9; ++idx ) {
    double Nch = cmult[idx]/gensumw[idx];
    cmult2[idx] = (cmult2[idx]/gensumw[idx] - pow2(Nch))/gensumn[idx];
    ofs << setprecision(2) << setw(4) << int(pclim[idx]*100.0 + 0.5)
        << setw(10) << Nch << setw(10) << sqrt(cmult2[idx]) <<endl;
  }
  ofs << "\n# Nwc\n";
  for ( int idx = 0; idx < 9; ++idx ) {
    double Nw = wound[idx]/gensumw[idx];
    wound2[idx] = (wound2[idx]/gensumw[idx] - pow2(Nw))/gensumn[idx];
    ofs << setprecision(2) << setw(4) << int(pclim[idx]*100.0 + 0.5)
        << setw(10) << Nw << setw(10) << sqrt(wound2[idx]) <<endl;
  }

  // Befor we end we print out some statistics. Also, we want to check
  // that our generated centrality classes were the same as we
  // guessed.
  pythia.stat();
  double curr = 0.0;
  double prev = 0.0;
  double acc = 0.0;
  int idxa = 8;
  double lim = sumw*(1.0 - pclim[idxa]);
  vector<double> newlim(9);
  for ( multimap<double, double>::iterator it = gencent.begin();
        it != gencent.end(); ++it ) {
    prev = curr;
    curr = it->first;
    double w = it->second;
    if ( acc < lim && acc + w >= lim ) {
      newlim[idxa--] = prev + (curr - prev)*(lim - acc)/w;
      if ( idxa < 0 ) break;
      lim = sumw*(1.0 - pclim[idxa]);
    }
    acc += w;
  }

  cout << "The generated limits between centrality classes in this run:\n"
       << "   %   assumed    actual\n";
  for ( int idx = 0; idx < 9; ++idx )
    cout << setw(4) << int(pclim[idx]*100.0 + 0.5)
         << setw(10) << fixed << setprecision(1) << genlim[idx]
         << setw(10) << fixed << setprecision(1) << newlim[idx] << endl;

  // And we're done!
  return 0;
}






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

void SimplePYTHIALoop(Int_t n=1000, Int_t jobID=0, Int_t tune = 350,Float_t SNN = 2760,char CUT='y',Float_t yncut=0.1,char DTYPE='r',char DMODE='a') {
  makeEventSample(n,jobID,tune,SNN,CUT,yncut,DTYPE,DMODE);
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
    retVal = makeEventSample(n,0,0,7.7,'y',0.1,'r','a');
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

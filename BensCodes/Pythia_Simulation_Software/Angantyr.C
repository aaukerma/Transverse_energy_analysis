/*
TODO: do it the old way for comparison: pion contribution added from omegas and etas
*/
#ifndef __CINT__
#include "Pythia8/Pythia.h"
#include "Pythia8/HeavyIons.h"
#include "TFile.h"
#include "TError.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TMath.h"
#include <iostream>
#include <vector>
#include <ctime>
#include <assert.h>
using namespace Pythia8;
//using namespace std;
#endif

#define FILENAME   "pythia.root"
#define TREENAME   "tree"
#define BRANCHNAME "particles"
#define HISTNAME   "ptSpectra"
#define PDGNUMBER  211
class TH3F;
class TH1F;
class TPythia8;

struct KF_Code {
  Int_t KFc;
  TString name;
};
struct pylista {
  Int_t inde;
  Int_t KFpart;
  Int_t indePar;
};
struct DEC {
  Int_t g;
  Double_t ET_a;
  Double_t ET_b;
  Double_t ET_c;
  Double_t ET_d;
  Double_t ET_e;
  Int_t a;
  Int_t b;
  Int_t c;
  Int_t d;
  Int_t e;
};

int test(){
return 0;
}

TH1F *CreateEnergyHistogram(char *name){
  TH1F *histo = new TH1F(name,name,1000,-0.9,99.9);
  histo->GetYaxis()->SetTitle("number of entries");
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
  TH1F *histo = new TH1F(name,name,1000,0,100);
  histo->GetYaxis()->SetTitle("number of entries");
  histo->GetXaxis()->SetTitle("Pt");
  return histo;
}
TH1F *CreateMultHistogram(char *name){
  TH1F *histo = new TH1F(name,name,1000,0,100);
  histo->GetYaxis()->SetTitle("number of entries");
  histo->GetXaxis()->SetTitle("Number of Particles");
  return histo;
}
TH1F *CreateJHistogram(char *name){
  TH1F *histo = new TH1F(name,name,1000,.999,1.5);
  histo->GetYaxis()->SetTitle("number of entries");
  histo->GetXaxis()->SetTitle("J");
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

DEC GetEtaDecayMode(vector<Particle>& DAU2){
  int x = 5;// DAU2.size();
  DEC ret;
  //following inputs daughter IDS
  int a=abs(DAU2[0].id());
  int b=abs(DAU2[1].id());
  int c=abs(DAU2[2].id());
  int d=abs(DAU2[3].id());
  int e=abs(DAU2[4].id());
  int f=0;
  int g=-1;
  //following modifies electron/positron code for unique combinations
  if(a==11) a=93;
  if(b==11) b=93;
  if(c==11) c=93;
  if(d==11) d=93;
  if(e==11) e=93;
  //add all particles IDs for unique decay mode code
  f=a+b+c+d+e;
  switch (f){
    case 44:{   //2gamma
      g=0;
      break;
    }
    case 333:{    //3pi0
      g=1;
      break;
    }
    case 533:{    //pi+,pi-,pi0
      g=2;
      break;
    }
    case 444:{    //pi+,pi-,gamma
      g=3;
      break;
    }
    case 208:{    //e+,e-,gamma
      g=4;
      break;
    }
    case 48:{    //mu+,mu-,gamma
      g=5;
      break;
    }
    case 155:{    //pi0,2gamma
      g=6;
      break;
    }
    case 266:{    //2pi0,2gamma
      g=7;
      break;
    }
    case 88:{    //4gamma
      g=8;
      break;
    }
    case 186:{    //e+,e-
      g=9;
      break;
    }
    case 26:{    //mu+,mu-
      g=10;
      break;
    }
    case 372:{    //2e+,2e-
      g=11;
      break;
    }
    case 630:{    //pi+,pi-,e+,e-,(gamma)
      g=12;
      break;
    }
    case 212:{    //e+,e-,mu+,mu-
      g=13;
      break;
    }
    case 52:{    //2mu+,2mu-
      g=14;
      break;
    }
    default:{    //other
      g=15;
      break;
    }
  }
  ret.g = g;
  ret.a = a;
  ret.b = b;
  ret.c = c;
  ret.d = d;
  ret.e = e;
  //cout<<" eta decay! mode: "<<g<<" "<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<e<<endl;
  for (int i=0;i<5;i++){
    Particle p=DAU2[i];
    Float_t partE;
    Float_t Theta = p.theta();
    Float_t pT= sqrt(pow(p.px(),2)+pow(p.py(),2));
    Float_t p_tot=(pT)/sin(Theta);
    Float_t Ei_TOT= sqrt(pow(p_tot,2)+pow(p.m(),2));
    if ((p.id()==3122)||(p.id()==2212)||(p.id()==2112)){
      partE = (Ei_TOT-p.m()) * sin(Theta);
    }
    else if ((p.id()==-3122)||(p.id()==-2212)||(p.id()==-2112)){
      partE = (Ei_TOT+p.m()) * sin(Theta);
    }
    else {
      partE = Ei_TOT * sin(Theta);
    }
    switch (i){
      case 0:{
        ret.ET_a = partE;
        break;
      }
      case 1:{
        ret.ET_b = partE;
        break;
      }
      case 2:{
        ret.ET_c = partE;
        break;
      }
      case 3:{
        ret.ET_d = partE;
        break;
      }
      case 4:{
        ret.ET_e = partE;
        break;
      }
    }
  }
  return ret;
}

DEC GetOmegaDecayMode(vector<Particle>& DAU2){
  int x = 5;// DAU2.size();
  DEC ret;
  //following inputs daughter IDS
  int a=abs(DAU2[0].id());
  int b=abs(DAU2[1].id());
  int c=abs(DAU2[2].id());
  int d=abs(DAU2[3].id());
  int e=abs(DAU2[4].id());
  int f=0;
  int g=-1;
  //following modifies electron/positron code for unique combinations
  if(a==11) a=93;
  if(b==11) b=93;
  if(c==11) c=93;
  if(d==11) d=93;
  if(e==11) e=93;
  //add all particles IDs for unique decay mode code
  f=a+b+c+d+e;
  switch (f){
    case 533:{
      g=0;
      break;
    }
    case 133:{
      g=1;
      break;
    }
    case 422:{
      g=2;
      break;
    }
    case 243:{
      g=3;
      break;
    }
    case 297:{
      g=4;
      break;
    }
    case 137:{
      g=5;
      break;
    }
    case 186:{
      g=6;
      break;
    }
    case 644:{
      g=7;
      break;
    }
    case 444:{
      g=8;
      break;
    }
    case 844:{
      g=9;
      break;
    }
    case 244:{
      g=10;
      break;
    }
    case 354:{
      g=11;
      break;
    }
    case 26:{
      g=12;
      break;
    }
    case 66:{
      g=13;
      break;
    }
    default:{
      g=14;
      break;
    }
  }
  ret.g = g;
  ret.a = a;
  ret.b = b;
  ret.c = c;
  ret.d = d;
  ret.e = e;
  //cout<<" omega decay! mode: "<<g<<" "<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<e<<endl;
  for (int i=0;i<5;i++){
    Particle p=DAU2[i];
    Float_t partE;
    Float_t Theta = p.theta();
    Float_t pT= sqrt(pow(p.px(),2)+pow(p.py(),2));
    Float_t p_tot=(pT)/sin(Theta);
    Float_t Ei_TOT= sqrt(pow(p_tot,2)+pow(p.m(),2));
    if ((p.id()==3122)||(p.id()==2212)||(p.id()==2112)){
      partE = (Ei_TOT-p.m()) * sin(Theta);
    }
    else if ((p.id()==-3122)||(p.id()==-2212)||(p.id()==-2112)){
      partE = (Ei_TOT+p.m()) * sin(Theta);
    }
    else {
      partE = Ei_TOT * sin(Theta);
    }
    switch (i){
      case 0:{
        ret.ET_a = partE;
        break;
      }
      case 1:{
        ret.ET_b = partE;
        break;
      }
      case 2:{
        ret.ET_c = partE;
        break;
      }
      case 3:{
        ret.ET_d = partE;
        break;
      }
      case 4:{
        ret.ET_e = partE;
        break;
      }
    }
  }
  return ret;
}

int makeEventSample(Int_t nEvents, Int_t jobID, Float_t SNN, char CUT, Float_t yncut, char DTYPE, char DMODE) {
/********************
Variable Declarations
*********************/
  //modifiable
  int print = 10;                 //Event Number printout frequency

  //DO NOT MODIFY
  int debug =0;                   //Runs debug statements to find bugs
  Int_t l=0;                      //counter for KF_CODE array builder
  int iAbort = 0;                 //counter of pythia.next() failures
  vector <KF_Code> partname(0);   //vector containing KF codes
  int modSNN;                     //integer version of SNN for file naming
  string pVERSION;                //version of pythia
  string temporarything="fail";   //temporary string
  int t=0;
  Float_t BeamHalfSNN;

/************************
Pythia startup and tuning
************************/
  cout<<"Running "<<nEvents<<" events, job id "<<jobID<<", at "<<SNN<<"GeV"<<endl;
  Pythia pythia;
  pythia.checkVersion();
  pVERSION="Pythia8 Angantyr";
  pythia.initPtrs();
  pythia.readFile("Angantyr.cmnd");
  char* SNNthing=Form("Beams:eCM = %f",SNN);
  if (SNN<6){
    pythia.readString("Beams:frameType = 2");
    BeamHalfSNN=SNN/2;
    char* SNNa=Form("Beams:eA = %f",BeamHalfSNN);
    char* SNNb=Form("Beams:eB = %f",BeamHalfSNN);
    pythia.readString(SNNa);
    pythia.readString(SNNb);
  }
  char* IDthing=Form("Random:seed = %i",jobID);
  char* Eventthing=Form("Main:numberOfEvents = %i",nEvents);
  pythia.readString(SNNthing);
  pythia.readString(IDthing);
  pythia.readString(Eventthing);
  int colparA = pythia.mode("Beams:idA");
  int colparB = pythia.mode("Beams:idB");
  pythia.settings.listChanged();
  if (!pythia.init()){
    cout << "something went wrong\n";
  }
  pythia.settings.listChanged();                      //list of changed tune switches only
  //pythia.settings.listAll();                        //Full list of tuning switches
  int nAbort = pythia.mode("Main:timesAllowErrors");  //number of allowed .next() errors
  cout<<"Pythia startup and tuning complete\n";
/**********
Input Files
************/
  ifstream in;
  in.open("KF_Code.dat");
  while (in.good()){                //Loads PDG/KF code indexing vector
    partname.push_back(KF_Code());
    in>>partname[l].KFc>>partname[l].name;
    l++;
  }//END PDG/KF code index loop
  cout<<"PDG and KF code database loaded\n";

/***********
Output Files
************/
  if (SNN==7.7){ //folloing if statements allow for file naming without decimals
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
  char ty;
  if (colparA == 2212 && colparB == 2212) ty='p';
  else if (colparA == 1000791970 && colparB == 1000791970) ty = 'A';
  else if (colparA == 1000020040 && colparB == 1000020040) ty = 'a';
  char *filename = Form("Data/root/p%c%i_%iGeVoutfile.root",ty,jobID,modSNN);
  char *filename2 = Form("Data/outfile/p%c%i_%iGeVoutfile.txt",ty,jobID,modSNN);
  char *filename3 = Form("Data/exotics/p%c%i_%iGeVexotics.txt",ty,jobID,modSNN);
  char *filename4 = Form("Data/decay/p%c%i_%iGeVDecayStats.txt",ty,jobID,modSNN);
  char *filename5 = Form("Data/list/p%c%i_%iGeVParticleList.txt",ty,jobID,modSNN);
  TFile* file = TFile::Open(filename, "RECREATE");
  if (!file || !file->IsOpen()) {
    Error("TEST", "Couldn't open file %s", filename);
    return 1;
  }
  ofstream out;
  ofstream out2;
  ofstream out3;
  ofstream out4;
  out.open(filename2);
  out2.open(filename3);
  out3.open(filename4);
  out4.open(filename5);

  //following loads header to output .txt files
  out<<"Event ET Analysis\n******************************************************\n";
  out<<"JobID: "<<jobID<<"\tNo. Events: "<<nEvents<<"\tSNN: "<<SNN<<endl<<endl;
  out<<"Simulation Parameters:\nVersion: "<<pVERSION<<"\nType of Cut: ";
  out2<<"Exotic Particle Analysis\n******************************************************\n";
  out2<<"JobID: "<<jobID<<"\tNo. Events: "<<nEvents<<"\tSNN: "<<SNN<<endl<<endl;
  out2<<"Simulation Parameters:\nVersion: "<<pVERSION<<"\nType of Cut: ";
  out3<<"Decay (eta and omega) Analysis\n******************************************************\n";
  out3<<"JobID: "<<jobID<<"\tNo. Events: "<<nEvents<<"\tSNN: "<<SNN<<endl<<endl;
  out3<<"Simulation Parameters:\nVersion: "<<pVERSION<<"\nType of Cut: ";
  out4<<"Particle List (1st event only)\n******************************************************\n";
  out4<<"JobID: "<<jobID<<"\tNo. Events: "<<nEvents<<"\tSNN: "<<SNN<<endl<<endl;
  out4<<"Simulation Parameters:\nVersion: "<<pVERSION<<"\nType of Cut: ";
  if (CUT=='y'){
    temporarything = "rapidity (y)";
  }
  else if (CUT=='n'){
    temporarything = "pseudorapidity (n)";
  }
  out<<temporarything<<"\nCut: "<<yncut<<"\nCalorimeter (r=regular, c=calorimeter): "<<DTYPE<<endl;
  out<<"Counting method: "<<DMODE<<endl<<endl;
  out2<<temporarything<<"\nCut: "<<yncut<<"\nCalorimeter (r=regular, c=calorimeter): "<<DTYPE<<endl;
  out2<<"Counting method: "<<DMODE<<endl<<endl;
  out3<<temporarything<<"\nCut: "<<yncut<<"\nCalorimeter (r=regular, c=calorimeter): "<<DTYPE<<endl;
  out3<<"Counting method: "<<DMODE<<endl<<endl;
  out2<<"Event\tParticle\tET\tpT\n";
  out4<<temporarything<<"\nCut: "<<yncut<<"\nCalorimeter (r=regular, c=calorimeter): "<<DTYPE<<endl;
  out4<<"Counting method: "<<DMODE<<endl<<endl;
  cout<<"Output files created\n";

/***********************
Histogram Intialization
************************/

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
  TH1F *hETOmegasub = CreateEnergyHistogram("hETomegasub");
  TH1F *hETLambda0 = CreateEnergyHistogram("hETLambda0");
  TH1F *hETLambdaBar0 = CreateEnergyHistogram("hETLambdaBar0");
  TH1F *hETOMEGAm = CreateEnergyHistogram("hETOMEGAm"); //3334
  TH1F *hETXi0 = CreateEnergyHistogram("hETXi0"); //3322
  TH1F *hETXim = CreateEnergyHistogram("hETXim"); //3312
  TH1F *hETSigmap = CreateEnergyHistogram("hETSigmap"); //3222
  TH1F *hETSigmam = CreateEnergyHistogram("hETSigmam"); //3112
  TH1F *hETSigma0 = CreateEnergyHistogram("hETSigma0"); //3212
  TH1F *hETSigmapBAR = CreateEnergyHistogram("hETSigmapBAR"); //3222
  TH1F *hETSigmamBAR = CreateEnergyHistogram("hETSigmamBAR"); //3112
  TH1F *hETSigma0BAR = CreateEnergyHistogram("hETSigma0BAR"); //3212
  TH1F *hETp = CreateEnergyHistogram("hETp");
  TH1F *hETn = CreateEnergyHistogram("hETn");
  TH1F *hETp_ = CreateEnergyHistogram("hETpBar");
  TH1F *hETn_ = CreateEnergyHistogram("hETnBar");
  TH1F *hETgamma = CreateEnergyHistogram("hETgamma");
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
  TH1F *hPTSigmapBAR = CreateEnergyHistogram("hPTSigmapBAR"); //3222
  TH1F *hPTSigmamBAR = CreateEnergyHistogram("hPTSigmamBAR"); //3112
  TH1F *hPTSigma0BAR = CreateEnergyHistogram("hPTSigma0BAR"); //3212
  TH1F *hPTp = CreatePTHistogram("hptp");
  TH1F *hPTn = CreatePTHistogram("hptn");
  TH1F *hPTp_ = CreatePTHistogram("hptpBar");
  TH1F *hPTn_ = CreatePTHistogram("hptnBar");
  TH1F *hPTgamma = CreatePTHistogram("hptgamma");
  TH1F *hPTother = CreatePTHistogram("hptother");
  TH1F *hNEvents = new TH1F("hNEvents","Number of events",1,0,1.0);
    hNEvents->GetYaxis()->SetTitle("N_{events}");
    hNEvents->GetXaxis()->SetTitle("no title");

  TH1F *hccpip = CreateMultHistogram("hccPiPlus");
  TH1F *hccpim = CreateMultHistogram("hccPiMinus");
  TH1F *hccpi0 = CreateMultHistogram("hccPi0");
  TH1F *hccKp = CreateMultHistogram("hccKPlus");
  TH1F *hccKm = CreateMultHistogram("hccKMinus");
  TH1F *hccKL = CreateMultHistogram("hccKL");
  TH1F *hccKS = CreateMultHistogram("hccKS");
  TH1F *hccEta = CreateMultHistogram("hccEta");
  TH1F *hccOmega = CreateMultHistogram("hccomega");
  TH1F *hccLambda0 = CreateMultHistogram("hccLambda0");
  TH1F *hccLambdaBar0 = CreateMultHistogram("hccLambdaBar0");
  TH1F *hccOMEGAm = CreateMultHistogram("hccOMEGAm"); //3334
  TH1F *hccXi0 = CreateMultHistogram("hccXi0"); //3322
  TH1F *hccXim = CreateMultHistogram("hccXim"); //3312
  TH1F *hccSigmap = CreateMultHistogram("hccSigmap"); //3222
  TH1F *hccSigmam = CreateMultHistogram("hccSigmam"); //3112
  TH1F *hccSigma0 = CreateMultHistogram("hccSigma0"); //3212
  TH1F *hccSigmapBAR = CreateMultHistogram("hccSigmapBAR"); //3222
  TH1F *hccSigmamBAR = CreateMultHistogram("hccSigmamBAR"); //3112
  TH1F *hccSigma0BAR = CreateMultHistogram("hccSigma0BAR"); //3212
  TH1F *hccp = CreateMultHistogram("hccp");
  TH1F *hccn = CreateMultHistogram("hccn");
  TH1F *hccp_ = CreateMultHistogram("hccpBar");
  TH1F *hccn_ = CreateMultHistogram("hccnBar");
  TH1F *hccgamma = CreateMultHistogram("hccgamma");
  TH1F *hccother = CreateMultHistogram("hccother");

  TH1F *hJHisto = CreateJHistogram("hJHisto");

  cout<<"Histogram intialization complete\n";

/***************************
Total Run Counters and Data
***************************/

  // individual particle counters

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
  Int_t cSigmapBAR=0;
  Int_t cSigmamBAR=0;
  Int_t cSigma0BAR=0;
  Int_t cp=0;
  Int_t cn=0;
  Int_t cp_=0;
  Int_t cn_=0;
  Int_t cother=0;
  Int_t cgamma=0;
  Int_t cETAall=0;
  Int_t cOmegaall=0;

  Int_t cetapis=0;
  Int_t comegapis=0;

  //Transverse energy
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
  Double_t pETSigmapBAR=0;
  Double_t pETSigmamBAR=0;
  Double_t pETSigma0BAR=0;
  Double_t pETp=0;
  Double_t pETn=0;
  Double_t pETp_=0;
  Double_t pETn_=0;
  Double_t pETgamma=0;
  Double_t pETother=0;
  Double_t pETETAall=0;
  Double_t pETOmegaall=0;
  Double_t ETCHECK=0;

  Double_t pETALL=0;

  Double_t ETetapis=0;
  Double_t ETomegapis=0;

  //Transverse Momenta
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
  Double_t pPTSigmapBAR=0;
  Double_t pPTSigmamBAR=0;
  Double_t pPTSigma0BAR=0;
  Double_t pPTp=0;
  Double_t pPTn=0;
  Double_t pPTp_=0;
  Double_t pPTn_=0;
  Double_t pPTgamma=0;
  Double_t pPTother=0;
  Double_t pPTETAall=0;
  Double_t pPTOmegaall=0;

  Double_t PTetapis=0;
  Double_t PTomegapis=0;

  //DECAY COUNTERS
  vector <double> ETADECAYS ={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  vector <double> OMEGADECAYS ={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  vector <double> ETADECAYSET ={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  vector <double> OMEGADECAYSET ={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  vector <double> ETADECAYSETa ={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  vector <double> OMEGADECAYSETa ={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  vector <double> ETADECAYSETb ={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  vector <double> OMEGADECAYSETb ={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  vector <double> ETADECAYSETc ={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  vector <double> OMEGADECAYSETc ={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  vector <double> ETADECAYSETd ={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  vector <double> OMEGADECAYSETd ={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  vector <double> ETADECAYSETe ={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  vector <double> OMEGADECAYSETe ={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  cout<<"Beginning Events\n";

/***********************
========================
      EVENT LOOP
========================
***********************/
  for ( int event = 0; event < nEvents; ++event ) {
    if (event % print == 0){
      cout <<"Collision Energy: "<< SNN <<" Event # " << event <<endl;
    }
    /**************************
    In-Event Counters and Data
    **************************/
    int counterA=0;
    int counterB=0;

    int npart = 0;
    //particle counters
    Int_t ccpip=0;
    Int_t ccpim=0;
    Int_t ccpi0=0;
    Int_t ccKp=0;
    Int_t ccKm=0;
    Int_t ccKL=0;
    Int_t ccKS=0;
    Int_t ccEta=0;
    Int_t ccOmega=0; //THIS IS omega not Omega
    Int_t ccLambda0=0;
    Int_t ccLambdaBar0=0;
    Int_t ccOMEGAm=0;
    Int_t ccXi0=0;
    Int_t ccXim=0;
    Int_t ccSigmap=0;
    Int_t ccSigmam=0;
    Int_t ccSigma0=0;
    Int_t ccSigmapBAR=0;
    Int_t ccSigmamBAR=0;
    Int_t ccSigma0BAR=0;
    Int_t ccp=0;
    Int_t ccn=0;
    Int_t ccp_=0;
    Int_t ccn_=0;
    Int_t ccgamma=0;
    Int_t ccother=0;
    Int_t ccETAall=0;
    Int_t ccOmegaall=0;

    //Transverse Energy
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
    Double_t ETSigmapBAR=0;
    Double_t ETSigmamBAR=0;
    Double_t ETSigma0BAR=0;
    Double_t ETp=0;
    Double_t ETn=0;
    Double_t ETp_=0;
    Double_t ETn_=0;
    Double_t ETgamma=0;
    Double_t ETother=0;
    Double_t etetap=0;
    Double_t etomep=0;

    //Transverse Momenta
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
    Double_t PTSigmapBAR=0;
    Double_t PTSigmamBAR=0;
    Double_t PTSigma0BAR=0;
    Double_t PTp=0;
    Double_t PTn=0;
    Double_t PTp_=0;
    Double_t PTn_=0;//int LambDecayCount=0;
    Double_t PTgamma=0;
    Double_t PTother=0;
    Double_t ptetap=0;
    Double_t ptomep=0;

    vector<pylista> pylis(0);            //particle list (akin to pylist())

    /**************************
    Event Generation Sequence
    **************************/
    if(debug)cout << "1 ";             //DEBUG STATEMENT
    Int_t NEXT =0;
    NEXT=pythia.next();
    if(debug)cout << "2 ";//DEBUG STATEMENT
    if (!NEXT) {
      iAbort++;
      cout<<" Critical Error: Failure "<<iAbort<<" of pythia.next()\n";
      if (iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }
    if(debug)cout << "3 ";             //DEBUG STATEMENT
    //cout<<"A ";
    Event* TheEvent =& pythia.event;
    npart = TheEvent->size();
    //cout<<"B\n";
    if (debug)cout<< "45";
    //cout<<npart<<endl;

    /****************************************
    Particle Parent and Decay Chain Structure
    *****************************************/
    for (int parta=0; parta<npart; parta++) {
      //cout<<"-";
      Particle particle = TheEvent->at(parta);
      pylis.push_back(pylista());
      pylis[parta].inde=parta;
      pylis[parta].indePar=particle.mother1();
      pylis[parta].KFpart=TheEvent->at(pylis[parta].indePar).id();
      //cout<<"=="<<pylis.size()<<"==";
      //cout<<pylis[parta].inde<<"\t"<<pylis[parta].indePar<<"\t"<<pylis[parta].KFpart<<endl;
    }
    if (event==0){
      out4<<"I\tparticle/jet\tstatus\tKF\torig\tp_x\tp_y\tp_z\tE\tm\n";
      out4<<"=========================================================\n";
    }
    /************************
    =========================
         PARTICLE LOOP
    =========================
    ************************/
    for (int part=0; part<npart; part++) {
      //cout<<".";
      Particle MPart = TheEvent->at(part);
      Int_t KFID = MPart.id();
      Int_t Ckf=GetKFConversion(KFID,partname);
      Int_t ind = part; //index   +1
      char N ='n';
      Float_t E= MPart.e();
      string NameOfParticle = MPart.name();
      Float_t partE;
      Float_t px=MPart.px();
      Float_t py=MPart.py();
      Float_t pz=MPart.pz();
      if (pz<0){
        pz*=(-1); //useful to keep numbers positive while calculating
        N='y';
      }
      Float_t p=sqrt(px*px+py*py+pz*pz);
      Float_t pT= sqrt(pow(px,2)+pow(py,2));
      Float_t Theta1 = atan(pT/pz);
      Float_t Theta = MPart.theta();
      Float_t p_tot=(pT)/sin(Theta);
      Float_t m=MPart.m();
      Float_t Ei_TOT= sqrt(pow(p_tot,2)+pow(m,2));
      //Float_t pseudorapidity=-log(tan(Theta/2));
//      cout<<-log(tan(Theta/2))<<" vs "<<0.5*log((p+pz)/(p-pz))<<"      and      "<<Theta<<" vs "<<Theta1<<endl;
      Float_t pseudorapidity = 0.5*log((p+pz)/(p-pz));
      Float_t rapidity=(0.5)*log((E+pz)/(E-pz));
      //cout<<pseudorapidity<<" "<<rapidity<<endl;
      Int_t mpartD = MPart.daughter1();
      Int_t mpartD2 = MPart.daughter2();
      vector<int> DAU = MPart.daughterList();
      int ds = DAU.size();
      vector<Particle> DAU2={0,0,0,0,0};
      Int_t mpP = MPart.mother1();
      Float_t mpL = MPart.tau0();
      bool finalp = MPart.isFinalPartonLevel();
      int STATUS = MPart.status();
      Int_t nobby=15;
      Int_t XEM=31;
      Int_t XEM2=31;
      DEC decayed;
      Float_t cuttype;
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
      else if (DTYPE=='c') partE=Ei_TOT;
      if (event==0){
        if (N =='y') pz*=(-1);
        out4<<ind<<"\t"<<NameOfParticle<<"\t"<<STATUS<<"\t"<<KFID<<"\t"<<mpP<<"\t"<<px<<"\t"<<py<<"\t"<<pz<<"\t"<<E<<"\t"<<m<<"\t"<<partE<<endl;
        if (N =='y') pz*=(-1);
      }
      if (CUT=='y'){
        cuttype=rapidity;
      }
      else if (CUT=='n'){
        cuttype=pseudorapidity;
      }
      if (cuttype<yncut){
        Double_t J = (pT*TMath::CosH(pseudorapidity))/sqrt(m*m+pT*pT*TMath::CosH(pseudorapidity));
        Double_t Ja = (pT*TMath::CosH(0))/sqrt(m*m+pT*pT*TMath::CosH(0));
        J=J/Ja;
        hJHisto->Fill(J);
        //if ((STATUS>=81)&&(STATUS<=86)){pETALL+=partE;} //all particles that are primary
        if (DMODE=='d'){ //COUTNING METHOD: Christine's preferred method
          if (KFID==211){ //pi+
            XEM=pylis[ind].KFpart;
            if ((XEM!=310)&&(XEM!=3122)&&(XEM!=-3122)&&(XEM!=3222)&&(XEM!=3112)&&(XEM!=3312)&&(XEM!=3322)&&(XEM!=3334)){ //not daughter of k0s,lam,lambar,sig+,sig-,xi0,xi-
              ETpip+=partE; //repeated for event
              pETpip+=partE; //repeated over total run
              ETAll+=partE;
              ETcharge+=partE;
              PTpip+=pT;
              pPTpip+=pT;
              PTAll+=pT;
              PTcharge=+pT;
              cpip++;  //repeated over total run
              ccpip++;
            }
            if (XEM==223){//if this pion is a daughter of an omega decay it will be subtracted from total omega ET/pt etc
              ETomegapis+=partE;
              etomep+=partE;
              PTomegapis+=pT;
              ptomep+=pT;
              comegapis++;
            }
            if (XEM==221){//if this pion is a daughter of an eta decay it will be subtracted from total omega ET/pt etc
              ETetapis+=partE;
              etetap+=partE;
              PTetapis+=pT;
              ptetap+=pT;
              cetapis++;
          }}
          if (KFID==-211){ //pi-
            XEM=pylis[ind].KFpart;
            if ((XEM!=310)&&(XEM!=3122)&&(XEM!=-3122)&&(XEM!=3222)&&(XEM!=3112)&&(XEM!=3312)&&(XEM!=3322)&&(XEM!=3334)){ //do not count if daughter of KOS
              ETpim+=partE;
              pETpim+=partE;
              ETAll+=partE;
              ETcharge+=partE;
              PTpim+=pT;
              pPTpim+=pT;
              PTAll+=pT;
              PTcharge=+pT;
              cpim++;
              ccpim++;
            }
            if (XEM==223){//if this pion is a daughter of an omega decay it will be subtracted from total omega ET/pt etc
              ETomegapis+=partE;
              etomep+=partE;
              PTomegapis+=pT;
              ptomep+=pT;
              comegapis++;
            }
            if (XEM==221){//if this pion is a daughter of an eta decay it will be subtracted from total omega ET/pt etc
              ETetapis+=partE;
              etetap+=partE;
              PTetapis+=pT;
              ptetap+=pT;
              cetapis++;
          }}
          if (KFID==111){ //pi0
            XEM=pylis[ind].KFpart; // is one of parents daughters a gamma? if so not included (this is to exclude certain omega and eta modes of decay)
            if ((XEM!=310)&&(XEM!=3122)&&(XEM!=-3122)&&(XEM!=3222)&&(XEM!=3112)&&(XEM!=3312)&&(XEM!=3322)&&(XEM!=3334)){ //do not count if daughter of KOS
              ETpi0+=partE;
              pETpi0+=partE;
              ETAll+=partE;
              PTpi0+=pT;
              pPTpi0+=pT;
              PTAll+=pT;
              cpi0++;
              ccpi0++;
            }
            if (XEM==223){//if this pion is a daughter of an omega decay it will be subtracted from total omega ET/pt etc
              ETomegapis+=partE;
              etomep+=partE;
              PTomegapis+=pT;
              ptomep+=pT;
              comegapis++;
            }
            if (XEM==221){//if this pion is a daughter of an eta decay it will be subtracted from total omega ET/pt etc
              ETetapis+=partE;
              etetap+=partE;
              PTetapis+=pT;
              ptetap+=pT;
              cetapis++;
            }
          }
          if (KFID==321){ //K+ //all included
            XEM=pylis[ind].KFpart;
            if ((XEM!=310)&&(XEM!=3122)&&(XEM!=-3122)&&(XEM!=3222)&&(XEM!=3112)&&(XEM!=3312)&&(XEM!=3322)&&(XEM!=3334)){ //not daughter of k0s,lam,lambar,sig+,sig-,xi0,xi-
              ETKp+=partE;
              pETKp+=partE;
              ETAll+=partE;
              ETcharge+=partE;
              PTKp+=pT;
              pPTKp+=pT;
              PTAll+=pT;
              PTcharge=+pT;
              cKp++;
              ccKp++;
          }}
          if (KFID==-321){ //K- //all included
            XEM=pylis[ind].KFpart;
            if ((XEM!=310)&&(XEM!=3122)&&(XEM!=-3122)&&(XEM!=3222)&&(XEM!=3112)&&(XEM!=3312)&&(XEM!=3322)&&(XEM!=3334)){ //not daughter of k0s,lam,lambar,sig+,sig-,xi0,xi-
              ETKm+=partE;
              pETKm+=partE;
              ETAll+=partE;
              ETcharge+=partE;
              PTKm+=pT;
              pPTKm+=pT;
              PTAll+=pT;
              PTcharge=+pT;
              cKm++;
              ccKm++;
          }}
          if (KFID==130){// KL
            ETKL+=partE;
            pETKL+=partE;
            ETAll+=partE;
            PTKL+=pT;
            pPTKL+=pT;
            PTAll+=pT;
            cKL++;
            ccKL++;
          }
          if (KFID==310){ //KS
            //THIS IS MEASURING ALL OF THEM!!!!!
            ETKS+=partE;
            pETKS+=partE;
            ETAll+=partE;
            PTKS+=pT;
            pPTKS+=pT;
            PTAll+=pT;
            cKS++;
            ccKS++;
          }
          if (KFID==221){ //eta
            for (int hxq=0;hxq<ds;hxq++){
              Particle XE = TheEvent->at(DAU[hxq]);
              DAU2[hxq]=XE;
            }
            decayed = GetEtaDecayMode(DAU2);
            ETADECAYS[decayed.g]+=1;
            ETADECAYSET[decayed.g]+=partE;
            ETADECAYSETa[decayed.g]+=decayed.ET_a;
            ETADECAYSETb[decayed.g]+=decayed.ET_b;
            ETADECAYSETc[decayed.g]+=decayed.ET_c;
            ETADECAYSETd[decayed.g]+=decayed.ET_d;
            ETADECAYSETe[decayed.g]+=decayed.ET_e;
            cETAall++;
            ETEta+=partE;
            pETEta+=partE;
            ETAll+=partE;
            PTEta+=pT;
            pPTEta+=pT;
            PTAll+=pT;
            cEta++;
            ccEta++;
          }
          if (KFID==223){ // omega
            nobby=0;
            for (int hxd=0;hxd<ds;hxd++){
              Particle XE = TheEvent->at(DAU[hxd]);
              DAU2[hxd]=XE;
            }
            decayed = GetOmegaDecayMode(DAU2);
            OMEGADECAYS[decayed.g]+=1;
            OMEGADECAYSET[decayed.g]+=partE;
            OMEGADECAYSETa[decayed.g]+=decayed.ET_a;
            OMEGADECAYSETb[decayed.g]+=decayed.ET_b;
            OMEGADECAYSETc[decayed.g]+=decayed.ET_c;
            OMEGADECAYSETd[decayed.g]+=decayed.ET_d;
            OMEGADECAYSETe[decayed.g]+=decayed.ET_e;
            cOmegaall++;
            ETOmega+=partE;
            pETOmega+=partE;
            ETAll+=partE;
            PTOmega+=pT;
            pPTOmega+=pT;
            PTAll+=pT;
            cOmega++;
            ccOmega++;
          }
          if (KFID==3122){ //Lambda0
            XEM=pylis[ind].KFpart;
            if ((XEM!=310)&&(XEM!=3222)&&(XEM!=3112)&&(XEM!=3312)&&(XEM!=3322)&&(XEM!=3334)){ //not daughter of k0s,lam,lambar,sig+,sig-,xi0,xi-
              ETLambda0+=partE;
              pETLambda0+=partE;
              ETAll+=partE;
              PTLambda0+=pT;
              pPTLambda0+=pT;
              PTAll+=pT;
              cLambda0++;
              ccLambda0++;
          }}
          if (KFID==-3122){ //LambdaBar0
            XEM=pylis[ind].KFpart;
            if ((XEM!=310)&&(XEM!=3222)&&(XEM!=3112)&&(XEM!=3312)&&(XEM!=3322)&&(XEM!=3334)){ //not daughter of k0s,lam,lambar,sig+,sig-,xi0,xi-
              ETLambdaBar0+=partE;
              pETLambdaBar0+=partE;
              ETAll+=partE;
              PTLambdaBar0+=pT;
              pPTLambdaBar0+=pT;
              PTAll+=pT;
              cLambdaBar0++;
              ccLambdaBar0++;
          }}
          if (KFID==3334){ //Omega-
            ETOMEGAm+=partE;
            pETOMEGAm+=partE;
            ETAll+=partE;
            PTOMEGAm+=pT;
            pPTOMEGAm+=pT;
            PTAll+=pT;
            PTcharge=+pT;
            cOMEGAm++;
            ccOMEGAm++;
          }
          if (KFID==3322){ //Xi0
            ETXi0+=partE;
            pETXi0+=partE;
            ETAll+=partE;
            PTXi0+=pT;
            pPTXi0+=pT;
            PTAll+=pT;
            cXi0++;
            ccXi0++;
          }
          if (KFID==3312){ //Xi-
            ETXim+=partE;
            pETXim+=partE;
            ETAll+=partE;
            PTXim+=pT;
            pPTXim+=pT;
            PTAll+=pT;
            PTcharge=+pT;
            cXim++;
            ccXim++;
          }
          if (KFID==3222){ //Sigma+
              ETSigmap+=partE;
              pETSigmap+=partE;
              ETAll+=partE;
              PTSigmap+=pT;
              pPTSigmap+=pT;
              PTAll+=pT;
              PTcharge=+pT;
              cSigmap++;
              ccSigmap++;
          }
          if (KFID==3112){ //sigma-
              ETSigmam+=partE;
              pETSigmam+=partE;
              ETAll+=partE;
              PTSigmam+=pT;
              pPTSigmam+=pT;
              PTAll+=pT;
              PTcharge=+pT;
              cSigmam++;
              ccSigmam++;
          }
          if (KFID==3212){ //sigma0
              ETSigma0+=partE;
              pETSigma0+=partE;
              ETAll+=partE;
              PTSigma0+=pT;
              pPTSigma0+=pT;
              PTAll+=pT;
              PTcharge=+pT;
              cSigma0++;
              ccSigma0++;
          }
          if (KFID==-3222){ //Sigma+BAR
              ETSigmapBAR+=partE;
              pETSigmapBAR+=partE;
              ETAll+=partE;
              PTSigmapBAR+=pT;
              pPTSigmapBAR+=pT;
              PTAll+=pT;
              PTcharge=+pT;
              cSigmapBAR++;
              ccSigmapBAR++;
          }
          if (KFID==-3112){ //sigma-BAR
              ETSigmamBAR+=partE;
              pETSigmamBAR+=partE;
              ETAll+=partE;
              PTSigmamBAR+=pT;
              pPTSigmamBAR+=pT;
              PTAll+=pT;
              PTcharge=+pT;
              cSigmamBAR++;
              ccSigmamBAR++;
          }
          if (KFID==-3212){ //sigma0BAR
              ETSigma0BAR+=partE;
              pETSigma0BAR+=partE;
              ETAll+=partE;
              PTSigma0BAR+=pT;
              pPTSigma0BAR+=pT;
              PTAll+=pT;
              PTcharge=+pT;
              cSigma0BAR++;
              ccSigma0BAR++;
          }
          if (KFID==2212){ //proton
            XEM=pylis[ind].KFpart;
            if ((XEM!=310)&&(XEM!=3122)&&(XEM!=-3122)&&(XEM!=3222)&&(XEM!=3112)&&(XEM!=3312)&&(XEM!=3322)&&(XEM!=3334)){ //not daughter of k0s,lam,lambar,sig+,sig-,xi0,xi-
              //are we measuing these right!!
              ETp+=partE;
              pETp+=partE;
              ETAll+=partE;
              ETcharge+=partE;
              PTp+=pT;
              pPTp+=pT;
              PTAll+=pT;
              PTcharge=+pT;
              cp++;
              ccp++;
          }}
          if (KFID==2112){ //neutron
            XEM=pylis[ind].KFpart;
            if ((XEM!=310)&&(XEM!=3122)&&(XEM!=-3122)&&(XEM!=3222)&&(XEM!=3112)&&(XEM!=3312)&&(XEM!=3322)&&(XEM!=3334)){ //not daughter of k0s,lam,lambar,sig+,sig-,xi0,xi-
              ETn+=partE;
              pETn+=partE;
              ETAll+=partE;
              PTn+=pT;
              pPTn+=pT;
              PTAll+=pT;
              cn++;
              ccn++;
          }}
          if (KFID==-2212){ //antiproton
            XEM=pylis[ind].KFpart;
            counterA++;
            if ((XEM!=310)&&(XEM!=3122)&&(XEM!=-3122)&&(XEM!=3222)&&(XEM!=3112)&&(XEM!=3312)&&(XEM!=3322)&&(XEM!=3334)){ //not daughter of k0s,lam,lambar,sig+,sig-,xi0,xi-
              ETp_+=partE;
                counterB++;
              pETp_+=partE;
              ETAll+=partE;
              ETcharge+=partE;
              PTp_+=pT;
              pPTp_+=pT;
              PTAll+=pT;
              PTcharge=+pT;
              cp_++;
              ccp_++;
          }}
          if (KFID==-2112){ //antineutron
            XEM=pylis[ind].KFpart;
            if ((XEM!=310)&&(XEM!=3122)&&(XEM!=-3122)&&(XEM!=3222)&&(XEM!=3112)&&(XEM!=3312)&&(XEM!=3322)&&(XEM!=3334)){ //not daughter of k0s,lam,lambar,sig+,sig-,xi0,xi-
              ETn_+=partE;
              pETn_+=partE;
              ETAll+=partE;
              PTn_+=pT;
              pPTn_+=pT;
              PTAll+=pT;
              cn_++;
              ccn_++;
          }}
          if (KFID==22){ //gamma
            XEM=pylis[ind].KFpart;
            if ((XEM!=111)&&(XEM!=221)&&(XEM!=223)&&(XEM!=3212)){
              ETgamma+=partE;
              pETgamma+=partE;
              ETAll+=partE;
              PTgamma+=pT;
              pPTgamma+=pT;
              PTAll+=pT;
              cgamma++;
              ccgamma++;
              out2<<event<<"\t"<<KFID<<"\t"<<partE<<"\t"<<pT<<"\t"<<STATUS<<endl;
            }
          }
          if ((KFID!=211)&&(KFID!=-211)&&(KFID!=111)&&(KFID!=321)&&(KFID!=-321)&&(KFID!=130)&&(KFID!=310)&&(KFID!=221)&&(KFID!=223)&&(KFID!=3122)&&(KFID!=-3122)&&(KFID!=3334)&&(KFID!=3322)&&(KFID!=3312)&&(KFID!=3222)&&(KFID!=3112)&&(KFID!=3212)&&(KFID!=2212)&&(KFID!=2112)&&(KFID!=-2212)&&(KFID!=-2112)&&(KFID!=22)&&((KFID>=100)||(KFID==11)||(KFID==12)||(KFID==13)||(KFID==14)||(KFID==15)||(KFID==16)||(KFID==17)||(KFID==18)||(KFID==-11)||(KFID==-12)||(KFID==-13)||(KFID==-14)||(KFID==-15)||(KFID==-16)||(KFID==-17)||(KFID==-18))) {
            if (mpartD==0){
              XEM=pylis[ind].KFpart;
              if ((XEM!=221)&&(XEM!=223)){
              ETother+=partE;
              pETother+=partE;
              ETAll+=partE;
              PTother+=pT;
              pPTother+=pT;
              PTAll+=pT;
              cother++;
              ccother++;
              out2<<event<<"\t"<<KFID<<"\t"<<partE<<"\t"<<pT<<"\t"<<STATUS<<endl;
          }}}
        }//END DMODE 'd'
      }//end rap cut
    }//end part loop
    /*****************************************
    Histogram Filling Sequence
    *****************************************/
    if(debug)cout << "4\n";             //DEBUG STATEMENT

    Double_t omegaET=0;
    Double_t pi0ET=0;
    Double_t ChAllRat=0;
    hETAll->Fill(ETAll);
    hETcharge->Fill(ETcharge);
    ChAllRat=ETcharge/ETAll;
    hETChOverAll->Fill(ETcharge,ChAllRat);
    pETALL+=ETAll;
    //cout<<"TOTAL: "<<ETAll<<endl;
    { //ALL HISTOGRAM FILL
    if (ETpip!=0){ //pi+ hist
      ETCHECK+=ETpip;
      hETpip->Fill(ETpip);
      hPTpip->Fill(PTpip);
      hccpip->Fill(ccpip);
    }
    if (ETpim!=0){ //pi- histo
      hETpim->Fill(ETpim);
      hPTpim->Fill(PTpim);
      hccpim->Fill(ccpim);
    }
    if (ETpi0!=0){ //pi0 histo
      hETpi0->Fill(ETpi0);
      hPTpi0->Fill(PTpi0);
      hccpi0->Fill(ccpi0);
    }
    if (ETKp!=0){ //K+ histo
      hETKp->Fill(ETKp);
      hPTKp->Fill(PTKp);
      hccKp->Fill(ccKp);
    }
    if (ETKm!=0){ //K- histo
      hETKm->Fill(ETKm);
      hPTKm->Fill(PTKm);
      hccKm->Fill(ccKm);
    }
    if (ETKL!=0){ //=K0L histo
      hETKL->Fill(ETKL);
      hPTKL->Fill(PTKL);
      hccKL->Fill(ccKL);
    }
    if (ETKS!=0){ //K0S histo
      hETKS->Fill(ETKS);
      hPTKS->Fill(PTKS);
      hccKS->Fill(ccKS);
    }

    //NOTE: The eta ET, PT, and Total count is subtracted by any ET, PT, count made
    //from daughter pions! same for omega

    if (ETEta!=0){ //eta histo
      hETEta->Fill(ETEta-etetap);
      hPTEta->Fill(PTEta-ptetap);
      hccEta->Fill(ccEta-cetapis);
    }
    double kick = ETOmega-etomep;
    if (kick<0) kick=0;
    if (ETOmega!=0){ //omega histo
      hETOmega->Fill(ETOmega);
      hETOmegasub->Fill(kick);
      hPTOmega->Fill(PTOmega-ptomep);
      hccOmega->Fill(ccOmega-comegapis);
    }

    if (ETLambda0!=0){ //Lambda0 histo
      hETLambda0->Fill(ETLambda0);
      hPTLambda0->Fill(PTLambda0);
      hccLambda0->Fill(ccLambda0);
    }
    if (ETLambdaBar0!=0){ //LambdaBar histo
      hETLambdaBar0->Fill(ETLambdaBar0);
      hPTLambdaBar0->Fill(PTLambdaBar0);
      hccLambdaBar0->Fill(ccLambdaBar0);
    }
    if (ETOMEGAm!=0){ //Omega histo
      hETOMEGAm->Fill(ETOMEGAm);
      hPTOMEGAm->Fill(PTOMEGAm);
      hccOMEGAm->Fill(ccOMEGAm);
    }
    if (ETXi0!=0){ //Xi0 histo
      hETXi0->Fill(ETXi0);
      hPTXi0->Fill(PTXi0);
      hccXi0->Fill(ccXi0);
    }
    if (ETXim!=0){ //Xi- histo
      hETXim->Fill(ETXim);
      hPTXim->Fill(PTXim);
      hccXim->Fill(ccXim);
    }
    if (ETSigmap!=0){ //Sigma+ histo
      hETSigmap->Fill(ETSigmap);
      hPTSigmap->Fill(PTSigmap);
      hccSigmap->Fill(ccSigmap);
    }
    if (ETSigmam!=0){ //Sigma- histo
      hETSigmam->Fill(ETSigmam);
      hPTSigmam->Fill(PTSigmam);
      hccSigmam->Fill(ccSigmam);
    }
    if (ETSigma0!=0){ //Sigma0 histo
      hETSigma0->Fill(ETSigma0);
      hPTSigma0->Fill(PTSigma0);
      hccSigma0->Fill(ccSigma0);
    }
    if (ETSigmapBAR!=0){ //SigmaBar+ histo
      hETSigmapBAR->Fill(ETSigmapBAR);
      hPTSigmapBAR->Fill(PTSigmapBAR);
      hccSigmapBAR->Fill(ccSigmapBAR);
    }
    if (ETSigmamBAR!=0){ //SigmaBar- histo
      hETSigmamBAR->Fill(ETSigmamBAR);
      hPTSigmamBAR->Fill(PTSigmamBAR);
      hccSigmamBAR->Fill(ccSigmamBAR);
    }
    if (ETSigma0BAR!=0){ //SigmaBar0 histo
      hETSigma0BAR->Fill(ETSigma0BAR);
      hPTSigma0BAR->Fill(PTSigma0BAR);
      hccSigma0BAR->Fill(ccSigma0BAR);
    }
    if (ETp!=0){ //proton histo
      hETp->Fill(ETp);
      hPTp->Fill(PTp);
      hccp->Fill(ccp);
    }
    if (ETn!=0){ //neutron histo
      hETn->Fill(ETn);
      hPTn->Fill(PTn);
      hccn->Fill(ccn);
    }
    if (ETp_!=0){ //antiproton histo
      hETp_->Fill(ETp_);
      hPTp_->Fill(PTp_);
      hccp_->Fill(ccp_);
    }
    if (ETn_!=0){ //antineutron histo
      hETn_->Fill(ETn_);
      hPTn_->Fill(PTn_);
      hccn_->Fill(ccn_);
    }
    if (ETgamma!=0){ //gamma histo
      hETgamma->Fill(ETgamma);
      hPTgamma->Fill(PTgamma);
      hccgamma->Fill(ccgamma);
    }
    if (ETother!=0){ //other/exotics histo
      hETother->Fill(ETother);
      hPTother->Fill(PTother);
      hccother->Fill(ccother);
    }
    }


    //cout<<counterA<<" "<<counterB<<endl;
    //TheEvent->free();
    TheEvent->reset();
  }//end event loop
  cout<<"Run complete, Writing to Disk\n";
/**************
Writing to File
**************/
  //to .root
  Double_t BinCent,BinCont,d=0,e,BinErr,STDDEV=0;
  /*
  for (Int_t i=0;i<1000;i++){
    BinCent=hETOmega->GetBinCenter(i);
    BinCont=hETOmega->GetBinContent(i);
    BinErr=sqrt(hETOmega->GetBinError(i));
    e=BinCent*BinCont;
    d+=e;
    e=BinErr*BinCont;
    STDDEV+=e;
  }*/
  Double_t omegaHIGH = hETOmega->Integral()*hETOmega->GetMean();
  Double_t omegaHIGHerr = STDDEV;
  Double_t omega = hETOmegasub->Integral()*hETOmegasub->GetMean();
  Double_t omegaerr = STDDEV;


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
  hETOmegasub->Write();
  hETOmega->Write();
  hETLambda0->Write();
  hETLambdaBar0->Write();
  hETOMEGAm->Write();
  hETXi0->Write();
  hETXim->Write();
  hETSigmap->Write();
  hETSigmam->Write();
  hETSigma0->Write();
  hETSigmapBAR->Write();
  hETSigmamBAR->Write();
  hETSigma0BAR->Write();
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
  hPTSigmapBAR->Write();
  hPTSigmamBAR->Write();
  hPTSigma0BAR->Write();
  hPTp->Write();
  hPTn->Write();
  hPTp_->Write();
  hPTn_->Write();
  hPTother->Write();
  hccpip->Write();
  hccpim->Write();
  hccpi0->Write();
  hccKp->Write();
  hccKm->Write();
  hccKL->Write();
  hccKS->Write();
  hccEta->Write();
  hccOmega->Write();
  hccLambda0->Write();
  hccLambdaBar0->Write();
  hccOMEGAm->Write();
  hccXi0->Write();
  hccXim->Write();
  hccSigmap->Write();
  hccSigmam->Write();
  hccSigma0->Write();
  hccp->Write();
  hccn->Write();
  hccp_->Write();
  hccn_->Write();
  hccother->Write();
  hJHisto->Write();


  //to .txt
  pETEta-=ETetapis;
  pETOmega-=ETomegapis;
  pPTEta-=PTetapis;
  pPTOmega-=PTomegapis;
  double TET=pETpip+pETpim+pETpi0+pETKp+pETKm+pETKL+pETKS+pETLambda0+pETLambdaBar0+pETp+pETn+pETp_+pETn_+pETEta+pETOmega+pETOMEGAm+pETXi0+pETXim+pETSigmap+pETSigmam+pETSigma0+pETgamma+pETother;
  double TPT=pPTpip+pPTpim+pPTpi0+pPTKp+pPTKm+pPTKL+pPTKS+pPTLambda0+pPTLambdaBar0+pPTp+pPTn+pPTp_+pPTn_+pPTEta+pPTOmega+pPTOMEGAm+pPTXi0+pPTXim+pPTSigmap+pPTSigmam+pPTSigma0+pPTgamma+pPTother;
  out<<"Total Transverse Energy: "<<TET<<"\tTotal Transverse Momentum: "<<TPT<<"\n\n";
  pETEta+=ETetapis;
  pETOmega+=ETomegapis;
  pPTEta+=PTetapis;
  pPTOmega+=PTomegapis;
  //out<<ETCHECK<<endl;
  out<<"Particle     Number\tEt\tPt\t%EtVsTotal\t%PtVsTotal\n";
  out<<"Pi+:         "<<cpip<<"\t"<<pETpip<<"\t"<<pPTpip<<"\t"<<(pETpip/TET)*100<<"\t"<<(pPTpip/TPT)*100<<endl;
  out<<"Pi-:         "<<cpim<<"\t"<<pETpim<<"\t"<<pPTpim<<"\t"<<(pETpim/TET)*100<<"\t"<<(pPTpim/TPT)*100<<endl;
  out<<"Pi0:         "<<cpi0<<"\t"<<pETpi0<<"\t"<<pPTpi0<<"\t"<<(pETpi0/TET)*100<<"\t"<<(pPTpi0/TPT)*100<<endl;
  out<<"K+:          "<<cKp<<"\t"<<pETKp<<"\t"<<pPTKp<<"\t"<<(pETKp/TET)*100<<"\t"<<(pPTKp/TPT)*100<<endl;
  out<<"K-:          "<<cKm<<"\t"<<pETKm<<"\t"<<pPTKm<<"\t"<<(pETKm/TET)*100<<"\t"<<(pPTKm/TPT)*100<<endl;
  out<<"K0L:         "<<cKL<<"\t"<<pETKL<<"\t"<<pPTKL<<"\t"<<(pETKL/TET)*100<<"\t"<<(pPTKL/TPT)*100<<endl;
  out<<"K0S:         "<<cKS<<"\t"<<pETKS<<"\t"<<pPTKS<<"\t"<<(pETKS/TET)*100<<"\t"<<(pPTKS/TPT)*100<<endl;
  out<<"Lambda0:     "<<cLambda0<<"\t"<<pETLambda0<<"\t"<<pPTLambda0<<"\t"<<(pETLambda0/TET)*100<<"\t"<<(pPTLambda0/TPT)*100<<endl;
  out<<"LambdaBar0:  "<<cLambdaBar0<<"\t"<<pETLambdaBar0<<"\t"<<pPTLambdaBar0<<"\t"<<(pETLambdaBar0/TET)*100<<"\t"<<(pPTLambdaBar0/TPT)*100<<endl;
  out<<"Proton:      "<<cp<<"\t"<<pETp<<"\t"<<pPTp<<"\t"<<(pETp/TET)*100<<"\t"<<(pPTp/TPT)*100<<endl;
  out<<"Antiproton:  "<<cp_<<"\t"<<pETp_<<"\t"<<pPTp_<<"\t"<<(pETp_/TET)*100<<"\t"<<(pPTp_/TPT)*100<<endl;
  out<<"Neutron:     "<<cn<<"\t"<<pETn<<"\t"<<pPTn<<"\t"<<(pETn/TET)*100<<"\t"<<(pPTn/TPT)*100<<endl;
  out<<"Antineutron: "<<cn_<<"\t"<<pETn_<<"\t"<<pPTn_<<"\t"<<(pETn_/TET)*100<<"\t"<<(pPTn_/TPT)*100<<endl;
  out<<"etaTOT:      "<<cEta<<"\t"<<pETEta<<"\t"<<pPTEta<<"\t"<<(pETEta/TET)*100<<"\t"<<(pPTEta/TPT)*100<<endl;
  out<<"omegaTOT:    "<<cOmega<<"\t"<<pETOmega<<"\t"<<omegaHIGH<<"\t"<<(pETOmega/TET)*100<<"\t"<<(pPTOmega/TPT)*100<<endl;
  pETEta-=ETetapis;
  pETOmega-=ETomegapis;
  pPTEta-=PTetapis;
  pPTOmega-=PTomegapis;
  out<<"eta:         "<<cEta<<"\t"<<pETEta<<"\t"<<pPTEta<<"\t"<<(pETEta/TET)*100<<"\t"<<(pPTEta/TPT)*100<<endl;
  out<<"omega:       "<<cOmega<<"\t"<<pETOmega<<"\t"<<omega<<"\t"<<(pETOmega/TET)*100<<"\t"<<(pPTOmega/TPT)*100<<endl;
  out<<"Omega-:      "<<cOMEGAm<<"\t"<<pETOMEGAm<<"\t"<<pPTOMEGAm<<"\t"<<(pETOMEGAm/TET)*100<<"\t"<<(pPTOMEGAm/TPT)*100<<endl;
  out<<"Xi0:         "<<cXi0<<"\t"<<pETXi0<<"\t"<<pPTXi0<<"\t"<<(pETXi0/TET)*100<<"\t"<<(pPTXi0/TPT)*100<<endl;
  out<<"Xi-:         "<<cXim<<"\t"<<pETXim<<"\t"<<pPTXim<<"\t"<<(pETXim/TET)*100<<"\t"<<(pPTXim/TPT)*100<<endl;
  out<<"Sigma+:      "<<cSigmap<<"\t"<<pETSigmap<<"\t"<<pPTSigmap<<"\t"<<(pETSigmap/TET)*100<<"\t"<<(pPTSigmap/TPT)*100<<endl;
  out<<"Sigma-:      "<<cSigmam<<"\t"<<pETSigmam<<"\t"<<pPTSigmam<<"\t"<<(pETSigmam/TET)*100<<"\t"<<(pPTSigmam/TPT)*100<<endl;
  out<<"Sigma0:      "<<cSigma0<<"\t"<<pETSigma0<<"\t"<<pPTSigma0<<"\t"<<(pETSigma0/TET)*100<<"\t"<<(pPTSigma0/TPT)*100<<endl;
  out<<"Sigma+BAR:   "<<cSigmapBAR<<"\t"<<pETSigmapBAR<<"\t"<<pPTSigmapBAR<<"\t"<<(pETSigmapBAR/TET)*100<<"\t"<<(pPTSigmapBAR/TPT)*100<<endl;
  out<<"Sigma-BAR:   "<<cSigmamBAR<<"\t"<<pETSigmamBAR<<"\t"<<pPTSigmamBAR<<"\t"<<(pETSigmamBAR/TET)*100<<"\t"<<(pPTSigmamBAR/TPT)*100<<endl;
  out<<"Sigma0BAR:   "<<cSigma0BAR<<"\t"<<pETSigma0BAR<<"\t"<<pPTSigma0BAR<<"\t"<<(pETSigma0BAR/TET)*100<<"\t"<<(pPTSigma0BAR/TPT)*100<<endl;
  out<<"gamma:       "<<cgamma<<"\t"<<pETgamma<<"\t"<<pPTgamma<<"\t"<<(pETgamma/TET)*100<<"\t"<<(pPTgamma/TPT)*100<<endl;
  out<<endl;
  out<<"[Pions]:     "<<(cpip+cpim+cpi0)<<"\t"<<(pETpip+pETpim+pETpi0)<<"\t"<<(pPTpip+pPTpim+pPTpi0)<<"\t"<<((pETpip+pETpim+pETpi0)/TET)*100<<"\t"<<((pPTpip+pPTpim+pPTpi0)/TPT)*100<<endl;
  out<<"[Kaons]:     "<<(cKp+cKm+cKL+cKS)<<"\t"<<(pETKp+pETKm+pETKL+pETKS)<<"\t"<<(pPTKp+pPTKm+pPTKL+pPTKS)<<"\t"<<((pETKp+pETKm+pETKL+pETKS)/TET)*100<<"\t"<<((pPTKp+pPTKm+pPTKL+pPTKS)/TPT)*100<<endl;
  out<<"[Lambdas]:   "<<(cLambda0+cLambdaBar0)<<"\t"<<(pETLambda0+pETLambdaBar0)<<"\t"<<(pPTLambda0+pPTLambdaBar0)<<"\t"<<((pETLambda0+pETLambdaBar0)/TET)*100<<"\t"<<((pPTLambda0+pPTLambdaBar0)/TPT)*100<<endl;
  out<<"[p,pbar]:    "<<(cp+cp_)<<"\t"<<(pETp+pETp_)<<"\t"<<(pPTp+pPTp_)<<"\t"<<((pETp+pETp_)/TET)*100<<"\t"<<((pPTp+pPTp_)/TPT)*100<<endl;
  out<<"[n,nbar]:    "<<(cn+cn_)<<"\t"<<(pETn+pETn_)<<"\t"<<(pPTn+pPTn_)<<"\t"<<((pETn+pETn_)/TET)*100<<"\t"<<((pPTn+pPTn_)/TPT)*100<<endl;
  out<<"[eta,ome]:   "<<(cEta+cOmega)<<"\t"<<(pETEta+pETOmega)<<"\t"<<(pPTEta+pPTOmega)<<"\t"<<((pETEta+pETOmega)/TET)*100<<"\t"<<((pPTEta+pPTOmega)/TPT)*100<<endl;
  out<<"[Xi,Om,Sig]: "<<(cOMEGAm+cXi0+cXim+cSigmap+cSigmam+cSigma0)<<"\t"<<(pETOMEGAm+pETXi0+pETXim+pETSigmap+pETSigmam+pETSigma0)<<"\t"<<(pPTOMEGAm+pPTXi0+pPTXim+pPTSigmap+pPTSigmam+pPTSigma0)<<"\t"<<((pETOMEGAm+pETXi0+pETXim+pETSigmap+pETSigmam+pETSigma0)/TET)*100<<"\t"<<((pPTOMEGAm+pPTXi0+pPTXim+pPTSigmap+pPTSigmam+pPTSigma0)/TPT)*100<<endl;
  out<<"[gamma]:     "<<cgamma<<"\t"<<pETgamma<<"\t"<<pPTgamma<<"\t"<<(pETgamma/TET)*100<<"\t"<<(pPTgamma/TPT)*100<<endl;
  out<<"[exotic]:    "<<cother<<"\t"<<pETother<<"\t"<<pPTother<<"\t"<<(pETother/TET)*100<<"\t"<<(pPTother/TPT)*100<<endl;
  out<<endl;
  out<<"counted ET:  "<<TET<<"\tTotal ET:\t"<<pETALL<<"\t(counted/total)%:\t"<<(TET/pETALL)*100<<endl;
  out<<endl;
  out<<"Decay Products vs 'Primary'\nGroup\tET\tPT\tORIGINAL: ETo\tPTo\tET%\tPT%\n";
  out<<"[Pions (eta)]   "<<ETetapis<<"\t"<<PTetapis<<"\t"<<(pETpip+pETpim+pETpi0)<<"\t"<<(pPTpip+pPTpim+pPTpi0)<<"\t"<<(ETetapis/(pETpip+pETpim+pETpi0))*100<<"\t"<<(PTetapis/(pPTpip+pPTpim+pPTpi0))*100<<endl;
  out<<"[Pions (omega)] "<<ETomegapis<<"\t"<<PTomegapis<<"\t"<<(pETpip+pETpim+pETpi0)<<"\t"<<(pPTpip+pPTpim+pPTpi0)<<"\t"<<(ETomegapis/(pETpip+pETpim+pETpi0))*100<<"\t"<<(PTomegapis/(pPTpip+pPTpim+pPTpi0))*100<<endl;

  pETEta+=ETetapis;
  pETOmega+=ETomegapis;
  pPTEta+=PTetapis;
  pPTOmega+=PTomegapis;
  for (int i=0;i<16;i++){
    double t = ETADECAYS[i];
    ETADECAYS[i]=(t/cETAall)*100.;
  }
  cout<<"\n";
  for (int j=0;j<15;j++){
    double t = OMEGADECAYS[j];
    OMEGADECAYS[j]=(t/cOmegaall)*100.;
  }

  out3<<"  Decay Mode       ET \t %of_particle_total \t %ofTotal_ET \t ET_a \t ET_b \t ET_c \t ET_d \t ET_e \t fraction \t expected [Particle Physics Booklet 2018]\n";
  out3<<"--------------eta Decay Branching Ratios-----------------\n";
  out3<<"  2gamma          "<<ETADECAYSET[0]<<"\t"<<(ETADECAYSET[0]/pETEta)*100<<"%\t"<<(ETADECAYSET[0]/TET)*100<<"%\t"<<ETADECAYSETa[0]<<"\t"<<ETADECAYSETb[0]<<"\t"<<ETADECAYSETc[0]<<"\t"<<ETADECAYSETd[0]<<"\t"<<ETADECAYSETe[0]<<"\t"<<ETADECAYS[0]<<"%\t 39.41%\n";
  out3<<"  3pi0            "<<ETADECAYSET[1]<<"\t"<<(ETADECAYSET[1]/pETEta)*100<<"%\t"<<(ETADECAYSET[1]/TET)*100<<"%\t"<<ETADECAYSETa[1]<<"\t"<<ETADECAYSETb[1]<<"\t"<<ETADECAYSETc[1]<<"\t"<<ETADECAYSETd[1]<<"\t"<<ETADECAYSETe[1]<<"\t"<<ETADECAYS[1]<<"%\t 32.68%\n";
  out3<<"  pi+,pi-,pi0     "<<ETADECAYSET[2]<<"\t"<<(ETADECAYSET[2]/pETEta)*100<<"%\t"<<(ETADECAYSET[2]/TET)*100<<"%\t"<<ETADECAYSETa[2]<<"\t"<<ETADECAYSETb[2]<<"\t"<<ETADECAYSETc[2]<<"\t"<<ETADECAYSETd[2]<<"\t"<<ETADECAYSETe[2]<<"\t"<<ETADECAYS[2]<<"%\t 22.92%\n";
  out3<<"  pi+,pi-,gamma   "<<ETADECAYSET[3]<<"\t"<<(ETADECAYSET[3]/pETEta)*100<<"%\t"<<(ETADECAYSET[3]/TET)*100<<"%\t"<<ETADECAYSETa[3]<<"\t"<<ETADECAYSETb[3]<<"\t"<<ETADECAYSETc[3]<<"\t"<<ETADECAYSETd[3]<<"\t"<<ETADECAYSETe[3]<<"\t"<<ETADECAYS[3]<<"%\t 4.22%\n";
  out3<<"  e+,e-,gamma     "<<ETADECAYSET[4]<<"\t"<<(ETADECAYSET[4]/pETEta)*100<<"%\t"<<(ETADECAYSET[4]/TET)*100<<"%\t"<<ETADECAYSETa[4]<<"\t"<<ETADECAYSETb[4]<<"\t"<<ETADECAYSETc[4]<<"\t"<<ETADECAYSETd[4]<<"\t"<<ETADECAYSETe[4]<<"\t"<<ETADECAYS[4]<<"%\t 6.9E-3%\n";
  out3<<"  mu+,mu-,gamma   "<<ETADECAYSET[5]<<"\t"<<(ETADECAYSET[5]/pETEta)*100<<"%\t"<<(ETADECAYSET[5]/TET)*100<<"%\t"<<ETADECAYSETa[5]<<"\t"<<ETADECAYSETb[5]<<"\t"<<ETADECAYSETc[5]<<"\t"<<ETADECAYSETd[5]<<"\t"<<ETADECAYSETe[5]<<"\t"<<ETADECAYS[5]<<"%\t 3.1E-4%\n";
  out3<<"  pi0,2gamma      "<<ETADECAYSET[6]<<"\t"<<(ETADECAYSET[6]/pETEta)*100<<"%\t"<<(ETADECAYSET[6]/TET)*100<<"%\t"<<ETADECAYSETa[6]<<"\t"<<ETADECAYSETb[6]<<"\t"<<ETADECAYSETc[6]<<"\t"<<ETADECAYSETd[6]<<"\t"<<ETADECAYSETe[6]<<"\t"<<ETADECAYS[6]<<"%\t 2.56E-4%\n";
  out3<<"  2pi0,2gamma     "<<ETADECAYSET[7]<<"\t"<<(ETADECAYSET[7]/pETEta)*100<<"%\t"<<(ETADECAYSET[7]/TET)*100<<"%\t"<<ETADECAYSETa[7]<<"\t"<<ETADECAYSETb[7]<<"\t"<<ETADECAYSETc[7]<<"\t"<<ETADECAYSETd[7]<<"\t"<<ETADECAYSETe[7]<<"\t"<<ETADECAYS[7]<<"%\t <1.2E-3%\n";
  out3<<"  4gamma          "<<ETADECAYSET[8]<<"\t"<<(ETADECAYSET[8]/pETEta)*100<<"%\t"<<(ETADECAYSET[8]/TET)*100<<"%\t"<<ETADECAYSETa[8]<<"\t"<<ETADECAYSETb[8]<<"\t"<<ETADECAYSETc[8]<<"\t"<<ETADECAYSETd[8]<<"\t"<<ETADECAYSETe[8]<<"\t"<<ETADECAYS[8]<<"%\t <2.8E-4%\n";
  out3<<"  e+,e-           "<<ETADECAYSET[9]<<"\t"<<(ETADECAYSET[9]/pETEta)*100<<"%\t"<<(ETADECAYSET[9]/TET)*100<<"%\t"<<ETADECAYSETa[9]<<"\t"<<ETADECAYSETb[9]<<"\t"<<ETADECAYSETc[9]<<"\t"<<ETADECAYSETd[9]<<"\t"<<ETADECAYSETe[9]<<"\t"<<ETADECAYS[9]<<"%\t <2.3E-6%\n";
  out3<<"  mu+,mu-         "<<ETADECAYSET[10]<<"\t"<<(ETADECAYSET[10]/pETEta)*100<<"%\t"<<(ETADECAYSET[10]/TET)*100<<"%\t"<<ETADECAYSETa[10]<<"\t"<<ETADECAYSETb[10]<<"\t"<<ETADECAYSETc[10]<<"\t"<<ETADECAYSETd[10]<<"\t"<<ETADECAYSETe[10]<<"\t"<<ETADECAYS[10]<<"%\t 5.8E-6%\n";
  out3<<"  2e+,2e-         "<<ETADECAYSET[11]<<"\t"<<(ETADECAYSET[11]/pETEta)*100<<"%\t"<<(ETADECAYSET[11]/TET)*100<<"%\t"<<ETADECAYSETa[11]<<"\t"<<ETADECAYSETb[11]<<"\t"<<ETADECAYSETc[11]<<"\t"<<ETADECAYSETd[11]<<"\t"<<ETADECAYSETe[11]<<"\t"<<ETADECAYS[11]<<"%\t 2.4E-5%\n";
  out3<<"  pi+,pi-,e+,e-,(gamma) "<<ETADECAYSET[12]<<"\t"<<(ETADECAYSET[12]/pETEta)*100<<"%\t"<<(ETADECAYSET[12]/TET)*100<<"%\t"<<ETADECAYSETa[12]<<"\t"<<ETADECAYSETb[12]<<"\t"<<ETADECAYSETc[12]<<"\t"<<ETADECAYSETd[12]<<"\t"<<ETADECAYSETe[12]<<"\t"<<ETADECAYS[12]<<"%\t 2.68E-4%\n";
  out3<<"  e+,e-,mu+,mu-   "<<ETADECAYSET[13]<<"\t"<<(ETADECAYSET[13]/pETEta)*100<<"%\t"<<(ETADECAYSET[13]/TET)*100<<"%\t"<<ETADECAYSETa[13]<<"\t"<<ETADECAYSETb[13]<<"\t"<<ETADECAYSETc[13]<<"\t"<<ETADECAYSETd[13]<<"\t"<<ETADECAYSETe[13]<<"\t"<<ETADECAYS[13]<<"%\t <1.6E-4%\n";
  out3<<"  2mu+,2mu-       "<<ETADECAYSET[14]<<"\t"<<(ETADECAYSET[14]/pETEta)*100<<"%\t"<<(ETADECAYSET[14]/TET)*100<<"%\t"<<ETADECAYSETa[14]<<"\t"<<ETADECAYSETb[14]<<"\t"<<ETADECAYSETc[14]<<"\t"<<ETADECAYSETd[14]<<"\t"<<ETADECAYSETe[14]<<"\t"<<ETADECAYS[14]<<"%\t <3.6E-4%\n";
  out3<<"  [OTHER]         "<<ETADECAYSET[15]<<"\t"<<(ETADECAYSET[15]/pETEta)*100<<"%\t"<<(ETADECAYSET[15]/TET)*100<<"%\t"<<ETADECAYSETa[15]<<"\t"<<ETADECAYSETb[15]<<"\t"<<ETADECAYSETc[15]<<"\t"<<ETADECAYSETd[15]<<"\t"<<ETADECAYSETe[15]<<"\t"<<ETADECAYS[15]<<"%\t -------%\n";
  out3<<"--------------omega Decay Branching Ratios-----------------\n";
  out3<<"  pi+,pi-,pi0     "<<OMEGADECAYSET[0]<<"\t"<<(OMEGADECAYSET[0]/pETOmega)*100<<"%\t"<<(OMEGADECAYSET[0]/TET)*100<<"%\t"<<OMEGADECAYSETa[0]<<"\t"<<OMEGADECAYSETb[0]<<"\t"<<OMEGADECAYSETc[0]<<"\t"<<OMEGADECAYSETd[0]<<"\t"<<OMEGADECAYSETe[0]<<"\t"<<OMEGADECAYS[0]<<"%\t 89.2%\n";
  out3<<"  pi0,gamma       "<<OMEGADECAYSET[1]<<"\t"<<(OMEGADECAYSET[1]/pETOmega)*100<<"%\t"<<(OMEGADECAYSET[1]/TET)*100<<"%\t"<<OMEGADECAYSETa[1]<<"\t"<<OMEGADECAYSETb[1]<<"\t"<<OMEGADECAYSETc[1]<<"\t"<<OMEGADECAYSETd[1]<<"\t"<<OMEGADECAYSETe[1]<<"\t"<<OMEGADECAYS[1]<<"%\t 8.40%\n";
  out3<<"  pi+,pi-         "<<OMEGADECAYSET[2]<<"\t"<<(OMEGADECAYSET[2]/pETOmega)*100<<"%\t"<<(OMEGADECAYSET[2]/TET)*100<<"%\t"<<OMEGADECAYSETa[2]<<"\t"<<OMEGADECAYSETb[2]<<"\t"<<OMEGADECAYSETc[2]<<"\t"<<OMEGADECAYSETd[2]<<"\t"<<OMEGADECAYSETe[2]<<"\t"<<OMEGADECAYS[2]<<"%\t 1.53%\n";
  out3<<"  eta,gamma       "<<OMEGADECAYSET[3]<<"\t"<<(OMEGADECAYSET[3]/pETOmega)*100<<"%\t"<<(OMEGADECAYSET[3]/TET)*100<<"%\t"<<OMEGADECAYSETa[3]<<"\t"<<OMEGADECAYSETb[3]<<"\t"<<OMEGADECAYSETc[3]<<"\t"<<OMEGADECAYSETd[3]<<"\t"<<OMEGADECAYSETe[3]<<"\t"<<OMEGADECAYS[3]<<"%\t 4.5E-4%\n";
  out3<<"  pi0,e+,e-       "<<OMEGADECAYSET[4]<<"\t"<<(OMEGADECAYSET[4]/pETOmega)*100<<"%\t"<<(OMEGADECAYSET[4]/TET)*100<<"%\t"<<OMEGADECAYSETa[4]<<"\t"<<OMEGADECAYSETb[4]<<"\t"<<OMEGADECAYSETc[4]<<"\t"<<OMEGADECAYSETd[4]<<"\t"<<OMEGADECAYSETe[4]<<"\t"<<OMEGADECAYS[4]<<"%\t 7.7E-4%\n";
  out3<<"  pi0,mu+,mu-     "<<OMEGADECAYSET[5]<<"\t"<<(OMEGADECAYSET[5]/pETOmega)*100<<"%\t"<<(OMEGADECAYSET[5]/TET)*100<<"%\t"<<OMEGADECAYSETa[5]<<"\t"<<OMEGADECAYSETb[5]<<"\t"<<OMEGADECAYSETc[5]<<"\t"<<OMEGADECAYSETd[5]<<"\t"<<OMEGADECAYSETe[5]<<"\t"<<OMEGADECAYS[5]<<"%\t 1.34E-4%\n";
  out3<<"  e+,e-           "<<OMEGADECAYSET[6]<<"\t"<<(OMEGADECAYSET[6]/pETOmega)*100<<"%\t"<<(OMEGADECAYSET[6]/TET)*100<<"%\t"<<OMEGADECAYSETa[6]<<"\t"<<OMEGADECAYSETb[6]<<"\t"<<OMEGADECAYSETc[6]<<"\t"<<OMEGADECAYSETd[6]<<"\t"<<OMEGADECAYSETe[6]<<"\t"<<OMEGADECAYS[6]<<"%\t 7.36E-5%\n";
  out3<<"  pi+,pi-,2pi0    "<<OMEGADECAYSET[7]<<"\t"<<(OMEGADECAYSET[7]/pETOmega)*100<<"%\t"<<(OMEGADECAYSET[7]/TET)*100<<"%\t"<<OMEGADECAYSETa[7]<<"\t"<<OMEGADECAYSETb[7]<<"\t"<<OMEGADECAYSETc[7]<<"\t"<<OMEGADECAYSETd[7]<<"\t"<<OMEGADECAYSETe[7]<<"\t"<<OMEGADECAYS[7]<<"%\t <2E-4%\n";
  out3<<"  pi+,pi-,gamma   "<<OMEGADECAYSET[8]<<"\t"<<(OMEGADECAYSET[8]/pETOmega)*100<<"%\t"<<(OMEGADECAYSET[8]/TET)*100<<"%\t"<<OMEGADECAYSETa[8]<<"\t"<<OMEGADECAYSETb[8]<<"\t"<<OMEGADECAYSETc[8]<<"\t"<<OMEGADECAYSETd[8]<<"\t"<<OMEGADECAYSETe[8]<<"\t"<<OMEGADECAYS[8]<<"%\t <3.6E-3%\n";
  out3<<"  2pi+,2pi-       "<<OMEGADECAYSET[9]<<"\t"<<(OMEGADECAYSET[9]/pETOmega)*100<<"%\t"<<(OMEGADECAYSET[9]/TET)*100<<"%\t"<<OMEGADECAYSETa[9]<<"\t"<<OMEGADECAYSETb[9]<<"\t"<<OMEGADECAYSETc[9]<<"\t"<<OMEGADECAYSETd[9]<<"\t"<<OMEGADECAYSETe[9]<<"\t"<<OMEGADECAYS[9]<<"%\t <1E-3%\n";
  out3<<"  2pi0,gamma      "<<OMEGADECAYSET[10]<<"\t"<<(OMEGADECAYSET[10]/pETOmega)*100<<"%\t"<<(OMEGADECAYSET[10]/TET)*100<<"%\t"<<OMEGADECAYSETa[10]<<"\t"<<OMEGADECAYSETb[10]<<"\t"<<OMEGADECAYSETc[10]<<"\t"<<OMEGADECAYSETd[10]<<"\t"<<OMEGADECAYSETe[10]<<"\t"<<OMEGADECAYS[10]<<"%\t 6.7E-5%\n";
  out3<<"  eta,pi0,gamma   "<<OMEGADECAYSET[11]<<"\t"<<(OMEGADECAYSET[11]/pETOmega)*100<<"%\t"<<(OMEGADECAYSET[11]/TET)*100<<"%\t"<<OMEGADECAYSETa[11]<<"\t"<<OMEGADECAYSETb[11]<<"\t"<<OMEGADECAYSETc[11]<<"\t"<<OMEGADECAYSETd[11]<<"\t"<<OMEGADECAYSETe[11]<<"\t"<<OMEGADECAYS[11]<<"%\t <3.3E-5%\n";
  out3<<"  mu+,mu-         "<<OMEGADECAYSET[12]<<"\t"<<(OMEGADECAYSET[12]/pETOmega)*100<<"%\t"<<(OMEGADECAYSET[12]/TET)*100<<"%\t"<<OMEGADECAYSETa[12]<<"\t"<<OMEGADECAYSETb[12]<<"\t"<<OMEGADECAYSETc[12]<<"\t"<<OMEGADECAYSETd[12]<<"\t"<<OMEGADECAYSETe[12]<<"\t"<<OMEGADECAYS[12]<<"%\t 7.4E-5%\n";
  out3<<"  3gamma          "<<OMEGADECAYSET[13]<<"\t"<<(OMEGADECAYSET[13]/pETOmega)*100<<"%\t"<<(OMEGADECAYSET[13]/TET)*100<<"%\t"<<OMEGADECAYSETa[13]<<"\t"<<OMEGADECAYSETb[13]<<"\t"<<OMEGADECAYSETc[13]<<"\t"<<OMEGADECAYSETd[13]<<"\t"<<OMEGADECAYSETe[13]<<"\t"<<OMEGADECAYS[13]<<"%\t <1.9E-4%\n";
  out3<<"  [OTHER]         "<<OMEGADECAYSET[14]<<"\t"<<(OMEGADECAYSET[14]/pETOmega)*100<<"%\t"<<(OMEGADECAYSET[14]/TET)*100<<"%\t"<<OMEGADECAYSETa[14]<<"\t"<<OMEGADECAYSETb[14]<<"\t"<<OMEGADECAYSETc[14]<<"\t"<<OMEGADECAYSETd[14]<<"\t"<<OMEGADECAYSETe[14]<<"\t"<<OMEGADECAYS[14]<<"%\t -------%\n";

  out.close();
  out2.close();
  out3.close();
  out4.close();
  file->Close();
  cout<<"Program Complete\n";
  return 0;
}//end program

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

void SimplePYTHIALoop(Int_t n=1000, Int_t jobID=0,Float_t SNN = 2760,char CUT='y',Float_t yncut=0.1,char DTYPE='r',char DMODE='d') {
  makeEventSample(n,jobID,SNN,CUT,yncut,DTYPE,DMODE);
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
    retVal = makeEventSample(n,0,7.7,'y',0.1,'r','d');
  else {
    retVal = showEventSample();
    app.Run();
  }

  return retVal;
}
#endif

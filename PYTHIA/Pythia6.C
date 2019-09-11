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
struct pylista {
  Int_t inde;
  Int_t KFpart;
  Int_t indePar;
};

// This function just load the needed libraries if we're executing from
// an interactive session.
int test(){
return 0;
}

TH1F *CreateEnergyHistogram(char *name){
  TH1F *histo = new TH1F(name,name,1024,0,32);
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
  TH1F *histo = new TH1F(name,name,1024,0,32);
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


// nEvents is how many events we want.
int makeEventSample(Int_t nEvents, Int_t jobID, Int_t tune, Float_t SNN, char CUT, Float_t yncut, char DTYPE, char DMODE, Float_t trigEtaMax = 0.5, Float_t assocEtaMax = 0.9){
/********************
Variable Declarations
*********************/
  //modifiable
  int print = 10;                 //Event Number printout frequency

  //DO NOT MODIFY
  Int_t l=0;                      //counter for KF_CODE array builder
  vector <KF_Code> partname(0);   //vector containing KF codes
  int modSNN;                     //integer version of SNN for file naming
  string temporarything="fail";   //temporary string
  string pVERSION;                //version of pythia

/************************
Pythia startup and tuning
************************/

// Create an instance of the Pythia event generator ...
TPythia6* pythia = new TPythia6();
pVERSION="Pythia6";//NOTE: check that the tune is correct

//turn on decays
//pythia->SetMSTJ(21,2);

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
pythia->SetMSTJ(22,1);//decays unstable particles

//NOTE: set seed as jobID
UInt_t seed = jobID;
cout<<seed<<endl; //gives seed
if( (seed>=0) && (seed<=900000000) ) {
  pythia->SetMRPY(1, seed);                   // set seed
  pythia->SetMRPY(2, 0);                      // use new seed
  pythia->SetPARP(2,7.0);                     // lowest collision energy THIS IS A MOD to allow (add more NOTE)
  cout<<"Random Seed : "<<seed<<endl;
} else {cout << "error: time " << seed << " is not valid" << endl; exit(2);}
cout<<"Running "<<nEvents<<" events, job id "<<jobID<<", tune "<<tune<<", at "<<SNN<<"GeV"<<endl;
// ... and initialise it to run p+p at sqrt(200) GeV in CMS
pythia->Initialize("cms", "p", "p", SNN); //NOTE:in GeV
//pythia->Dump();

cout<<"Pythia startup and tuning complete\n";



/**********
Input Files
************/
  ifstream in;
  in.open("KF_Code.dat");
  while (in.good()){
    partname.push_back(KF_Code());
    in>>partname[l].KFc>>partname[l].name;
    l++;
  }//End PDG/KF Code
  cout<<"PDG and KF code database loaded\n";

/***********
Output Files
************/
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
  char *filename = Form("pp%i_%iGeVoutfile.root",jobID,modSNN);
  char *filename2 = Form("pp%i_%iGeVoutfile.txt",jobID,modSNN);
  char *filename3 = Form("pp%i_%iGeVexotics.txt",jobID,modSNN);
  TFile* file = TFile::Open(filename, "RECREATE");
  TFile *outfile = file;//new TFile(filename,"RECREATE");
  if (!file || !file->IsOpen()) {
    Error("makeEventSample", "Couldn;t open file %s", filename);
    return 1;
  }
  ofstream out;
  ofstream out2;
  out.open(filename2);
  out2.open(filename3);

  out<<"RUN DATA\n******************************************************\n";
  out<<"JobID: "<<jobID<<"\tNo. Events: "<<nEvents<<"\tSNN: "<<SNN<<"\t tune: "<<tune<<endl<<endl;
  out<<"Simulation Parameters:\nVersion: "<<pVERSION<<"\nType of Cut: ";
  out2<<"Exotic Particle Analysis\n******************************************************\n";
  out2<<"JobID: "<<jobID<<"\tNo. Events: "<<nEvents<<"\tSNN: "<<SNN<<endl<<endl;
  out2<<"Simulation Parameters:\nVersion: "<<pVERSION<<"\nType of Cut: ";
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
  out2<<"Event\tParticle\tET\tpT\n";
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
    TH1F *hccp = CreateMultHistogram("hccp");
    TH1F *hccn = CreateMultHistogram("hccn");
    TH1F *hccp_ = CreateMultHistogram("hccpBar");
    TH1F *hccn_ = CreateMultHistogram("hccnBar");
    TH1F *hccgamma = CreateMultHistogram("hccgamma");
    TH1F *hccother = CreateMultHistogram("hccother");
    cout<<"Histogram intialization complete\n";

  /***************************
  Total Run Counters and Data
  ***************************/

  TClonesArray* particles = (TClonesArray*)pythia->GetListOfParticles();
  //tree->Branch(BRANCHNAME, &particles);

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
  Int_t cp=0;
  Int_t cn=0;
  Int_t cp_=0;
  Int_t cn_=0;
  Int_t cother=0;
  Int_t cETAall=0;
  Int_t cOmegaall=0;

  Int_t cetapis=0;
  Int_t comegapis=0;

  //Transverse Energy
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
  Double_t pPTp=0;
  Double_t pPTn=0;
  Double_t pPTp_=0;
  Double_t pPTn_=0;
  Double_t pPTother=0;
  Double_t pPTETAall=0;
  Double_t pPTOmegaall=0;

  Double_t PTetapis=0;
  Double_t PTomegapis=0;

/***********************
========================
      EVENT LOOP
========================
***********************/
  for (Int_t event = 0; event < nEvents; event++) {
    if (event % print == 0)
      cout <<"Collision Energy: "<< SNN <<" Event # " << event <<endl;
    /**************************
    In-Event Counters and Data
    **************************/
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
    Int_t ccp=0;
    Int_t ccn=0;
    Int_t ccp_=0;
    Int_t ccn_=0;
    Int_t ccgamma=0;
    Int_t ccother=0;
    Int_t ccETAall=0;
    Int_t ccOmegaall=0;

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
    Double_t ETgamma=0;
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
    Double_t PTgamma=0;
    Double_t PTother=0;

    vector<pylista> pylis(1);            //particle list (akin to pylist())

    /**************************
    Event Generation Sequence
    **************************/

    pythia->GenerateEvent();
      hNEvents->Fill(0.5);
    //pythia->Pylist(1);                  //DEBUG helper (prints full particle list)
    Int_t npart = particles->GetEntries();

    /****************************************
    Particle Parent and Decay Chain Structure
    *****************************************/
    for (Int_t parta=0; parta<npart; parta++) {
      pylis.push_back(pylista());
      pylis[parta].inde=parta;
      pylis[parta].KFpart=pythia->GetK(parta,2);
      pylis[parta].indePar=pythia->GetK(parta,3);
    }

    /************************
    =========================
         PARTICLE LOOP
    =========================
    ************************/

    for (Int_t part=0; part<npart; part++) {
	    TMCParticle *MPart = (TMCParticle *) particles->At(part);
      Int_t KFID = MPart->GetKF();
      Int_t Ckf=GetKFConversion(KFID,partname);
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
      Int_t mpartD = MPart->GetFirstChild();
      Int_t mpartD2 = MPart->GetLastChild();
      Int_t mpP = MPart->GetParent();
      Float_t mpL = MPart->GetLifetime();
      Int_t nobby=15;
      Int_t XEM=pylis[mpP].KFpart;
      Float_t cuttype;
      if (CUT=='y'){
        cuttype=rapidity;
      }
      else if (CUT=='n'){
        cuttype=pseudorapidity;
      }
      if (cuttype<yncut){
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
        if ((XEM==92)||(XEM==91)){pETALL+=partE;}
        if (DMODE=='a'){ //ALL INCLUDED
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
        }//END DMODE 'a'
        if (DMODE=='b'){ //ONLY IF FINAL
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
            if ((mpartD==0)){
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
            if ((mpartD==0)){
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
          }
          if (KFID==130){// KL
            //if ((mpL>=100)||(mpartD==0)){//counted if final particle or lifetime > 1cm
            if ((mpartD==0)){
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
            if ((mpartD==0)){
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
            if ((mpartD==0)){
              ETEta+=partE;
              pETEta+=partE;
              ETAll+=partE;
              PTEta+=pT;
              pPTEta+=pT;
              PTAll+=pT;
              cEta++;
            }
          }
          if (KFID==223){ // omega
            //cout<<ind<<"= particle number (check omega) running func:"<<endl;
            //  nobby=GetDaughterCheck(pylis,ind,223,22);
            //cout<<nobby<<endl;
            //if (nobby==1){
            if ((mpartD==0)){
              ETOmega+=partE;
              pETOmega+=partE;
              ETAll+=partE;
              PTOmega+=pT;
              pPTOmega+=pT;
              PTAll+=pT;
              cOmega++;
            }
          }
          if (KFID==3122){ //Lambda0
            //if ((mpL>=100)||(mpartD==0)){
            if ((mpartD==0)){
              ETLambda0+=partE;
              pETLambda0+=partE;
              ETAll+=partE;
              PTLambda0+=pT;
              pPTLambda0+=pT;
              PTAll+=pT;
              cLambda0++;
            }
          }
          if (KFID==-3122){ //LambdaBar0
            //if ((mpL>=100)||(mpartD==0)){
            if ((mpartD==0)){
              ETLambdaBar0+=partE;
              pETLambdaBar0+=partE;
              ETAll+=partE;
              PTLambdaBar0+=pT;
              pPTLambdaBar0+=pT;
              PTAll+=pT;
              cLambdaBar0++;
            }
          }
          if (KFID==3334){ //Omega-
            //if ((mpL>=100)||(mpartD==0)){
            if ((mpartD==0)){
              ETOMEGAm+=partE;
              pETOMEGAm+=partE;
              ETAll+=partE;
              PTOMEGAm+=pT;
              pPTOMEGAm+=pT;
              PTAll+=pT;
              PTcharge=+pT;
              cOMEGAm++;
            }
          }
          if (KFID==3322){ //Xi0
            //if ((mpL>=100)||(mpartD==0)){
            if ((mpartD==0)){
              ETXi0+=partE;
              pETXi0+=partE;
              ETAll+=partE;
              PTXi0+=pT;
              pPTXi0+=pT;
              PTAll+=pT;
              cXi0++;
            }
          }
          if (KFID==3312){ //Xi-
            //if ((mpL>=100)||(mpartD==0)){
            if ((mpartD==0)){
              ETXim+=partE;
              pETXim+=partE;
              ETAll+=partE;
              PTXim+=pT;
              pPTXim+=pT;
              PTAll+=pT;
              PTcharge=+pT;
              cXim++;
            }
          }
          if (KFID==3222){ //Sigma+
            //if ((mpL>=100)||(mpartD==0)){
            if ((mpartD==0)){
              ETSigmap+=partE;
              pETSigmap+=partE;
              ETAll+=partE;
              PTSigmap+=pT;
              pPTSigmap+=pT;
              PTAll+=pT;
              PTcharge=+pT;
              cSigmap++;
            }
          }
          if (KFID==3112){ //sigma-
            //if ((mpL>=100)||(mpartD==0)){
            if ((mpartD==0)){
              ETSigmam+=partE;
              pETSigmam+=partE;
              ETAll+=partE;
              PTSigmam+=pT;
              pPTSigmam+=pT;
              PTAll+=pT;
              PTcharge=+pT;
              cSigmam++;
            }
          }
          if (KFID==3212){ //sigma0
            //if ((mpL>=100)||(mpartD==0)){
            if ((mpartD==0)){
              ETSigma0+=partE;
              pETSigma0+=partE;
              ETAll+=partE;
              PTSigma0+=pT;
              pPTSigma0+=pT;
              PTAll+=pT;
              PTcharge=+pT;
              cSigma0++;
            }
          }
          if (KFID==2212){ //proton
            //if (mpP!=0){ //THIS IS AT MID RAPIDITY so this is not really required
            if ((mpartD==0)){
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
          }
          if (KFID==2112){ //neutron
            if ((mpartD==0)){
              ETn+=partE;
              pETn+=partE;
              ETAll+=partE;
              PTn+=pT;
              pPTn+=pT;
              PTAll+=pT;
              cn++;
            }
          }
          if (KFID==-2212){ //antiproton
            if ((mpartD==0)){
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
          }
          if (KFID==-2112){ //antineutron
            if ((mpartD==0)){
              ETn_+=partE;
              pETn_+=partE;
              ETAll+=partE;
              PTn_+=pT;
              pPTn_+=pT;
              PTAll+=pT;
              cn_++;
            }
          }
          if ((KFID==411)||(KFID==421)||(KFID==431)||(KFID==511)||(KFID==521)||(KFID==531)||(KFID==541)||(KFID==331)||(KFID==441)||(KFID==551)||(KFID==333)||(KFID==443)||(KFID==553)||(KFID==110)||(KFID==4112)||(KFID==4122)||(KFID==4212)||(KFID==4222)||(KFID==4132)||(KFID==4312)||(KFID==4232)||(KFID==4322)||(KFID==4332)||(KFID==5112)||(KFID==5122)||(KFID==5212)||(KFID==5222)||(KFID==3114)||(KFID==3214)||(KFID==3224)||(KFID==3314)||(KFID==3324)||(KFID==4114)||(KFID==4214)||(KFID==4224)||(KFID==4314)||(KFID==4324)||(KFID==4334)||(KFID==5114)||(KFID==5214)||(KFID==5224)) {
            if ((mpartD==0)){
              ETother+=partE;
              pETother+=partE;
              ETAll+=partE;
              PTother+=partE;
              pPTother+=partE;
              PTAll+=partE;
              cother++;
              out<<event<<"\t"<<KFID<<endl;
            }
          }
        }//END DMODE 'b'
        if (DMODE=='c'){ //Decay version 1
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
            }
          }
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
            }
          }
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
              }
            }
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
            }
          }
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
            }
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
            if ((mpL>=100)||(mpartD==0)){
              ETSigmap+=partE;
              pETSigmap+=partE;
              ETAll+=partE;
              PTSigmap+=pT;
              pPTSigmap+=pT;
              PTAll+=pT;
              PTcharge=+pT;
              cSigmap++;
            }
          }
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
            }
          }
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
            }
          }
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
            }
          }
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
            }
          }
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
            }
          }
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
            }
          }
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
            }
          }

        }//END DMODE 'c'
        if (DMODE=='d'){ //Christines preferred method
          if (KFID==211){ //pi+
            XEM=pylis[ind].KFpart;
            if ((XEM!=310)&&(XEM!=223)&&(XEM!=221)&&(XEM!=3122)&&(XEM!=-3122)&&(XEM!=3222)&&(XEM!=3112)&&(XEM!=3312)&&(XEM!=3322)&&(XEM!=3334)){ //not daughter of k0s,lam,lambar,sig+,sig-,xi0,xi-
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
            if (XEM==223){
              ETomegapis+=partE;
              PTomegapis+=pT;
              comegapis++;
            }
            if (XEM==221){
              ETetapis+=partE;
              PTetapis+=pT;
              cetapis++;
          }}
          if (KFID==-211){ //pi-
            XEM=pylis[ind].KFpart;
            if ((XEM!=310)&&(XEM!=223)&&(XEM!=221)&&(XEM!=3122)&&(XEM!=-3122)&&(XEM!=3222)&&(XEM!=3112)&&(XEM!=3312)&&(XEM!=3322)&&(XEM!=3334)){ //do not count if daughter of KOS
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
            if (XEM==223){
              ETomegapis+=partE;
              PTomegapis+=pT;
              comegapis++;
            }
            if (XEM==221){
              ETetapis+=partE;
              PTetapis+=pT;
              cetapis++;
          }}
          if (KFID==111){ //pi0
            XEM=pylis[ind].KFpart; // is one of parents daughters a gamma? if so not included (this is to exclude certain omega and eta modes of decay)
            if ((XEM!=310)&&(XEM!=223)&&(XEM!=221)&&(XEM!=3122)&&(XEM!=-3122)&&(XEM!=3222)&&(XEM!=3112)&&(XEM!=3312)&&(XEM!=3322)&&(XEM!=3334)){ //do not count if daughter of KOS
              ETpi0+=partE;
              pETpi0+=partE;
              ETAll+=partE;
              PTpi0+=pT;
              pPTpi0+=pT;
              PTAll+=pT;
              cpi0++;
              ccpi0++;
            }
            if (XEM==223){
              ETomegapis+=partE;
              PTomegapis+=pT;
              comegapis++;
            }
            if (XEM==221){
              ETetapis+=partE;
              PTetapis+=pT;
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
          if (KFID==2212){ //proton
            XEM=pylis[ind].KFpart;
            if ((XEM!=310)&&(XEM!=3122)&&(XEM!=-3122)&&(XEM!=3222)&&(XEM!=3112)&&(XEM!=3312)&&(XEM!=3322)&&(XEM!=3334)){ //not daughter of k0s,lam,lambar,sig+,sig-,xi0,xi-
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
            if ((XEM!=310)&&(XEM!=3122)&&(XEM!=-3122)&&(XEM!=3222)&&(XEM!=3112)&&(XEM!=3312)&&(XEM!=3322)&&(XEM!=3334)){ //not daughter of k0s,lam,lambar,sig+,sig-,xi0,xi-
              ETp_+=partE;
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
      } //end pseudorapidity cut
    }//end particle loop
    /*****************************************
    Histogram Filling Sequence
    *****************************************/
    Double_t omegaET=0;
    Double_t pi0ET=0;
    Double_t ChAllRat=0;
    hETAll->Fill(ETAll);
    hETcharge->Fill(ETcharge);
    ChAllRat=ETcharge/ETAll;
    hETChOverAll->Fill(ETcharge,ChAllRat);
    //cout<<"TOTAL: "<<ETAll<<endl;
    if (ETpip!=0){
      hETpip->Fill(ETpip);
      hPTpip->Fill(PTpip);
      hccpip->Fill(ccpip);
    }
    if (ETpim!=0){
      hETpim->Fill(ETpim);
      hPTpim->Fill(PTpim);
      hccpim->Fill(ccpim);
    }
    if (ETpi0!=0){
      hETpi0->Fill(ETpi0);
      hPTpi0->Fill(PTpi0);
      hccpi0->Fill(ccpi0);
    }
    if (ETKp!=0){
      hETKp->Fill(ETKp);
      hPTKp->Fill(PTKp);
      hccKp->Fill(ccKp);
    }
    if (ETKm!=0){
      hETKm->Fill(ETKm);
      hPTKm->Fill(PTKm);
      hccKm->Fill(ccKm);
    }
    if (ETKL!=0){
      hETKL->Fill(ETKL);
      hPTKL->Fill(PTKL);
      hccKL->Fill(ccKL);
    }
    if (ETKS!=0){
      hETKS->Fill(ETKS);
      hPTKS->Fill(PTKS);
      hccKS->Fill(ccKS);
    }
    if (ETEta!=0){
      hETEta->Fill(ETEta);
      hPTEta->Fill(PTEta);
      hccEta->Fill(ccEta);
    }
    if (ETOmega!=0){
      hETOmega->Fill(ETOmega);
      hPTOmega->Fill(PTOmega);
      hccOmega->Fill(ccOmega);
    }
    if (ETLambda0!=0){
      hETLambda0->Fill(ETLambda0);
      hPTLambda0->Fill(PTLambda0);
      hccLambda0->Fill(ccLambda0);
    }
    if (ETLambdaBar0!=0){
      hETLambdaBar0->Fill(ETLambdaBar0);
      hPTLambdaBar0->Fill(PTLambdaBar0);
      hccLambdaBar0->Fill(ccLambdaBar0);
    }
    if (ETOMEGAm!=0){
      hETOMEGAm->Fill(ETOMEGAm);
      hPTOMEGAm->Fill(PTOMEGAm);
      hccOMEGAm->Fill(ccOMEGAm);
    }
    if (ETXi0!=0){
      hETXi0->Fill(ETXi0);
      hPTXi0->Fill(PTXi0);
      hccXi0->Fill(ccXi0);
    }
    if (ETXim!=0){
      hETXim->Fill(ETXim);
      hPTXim->Fill(PTXim);
      hccXim->Fill(ccXim);
    }
    if (ETSigmap!=0){
      hETSigmap->Fill(ETSigmap);
      hPTSigmap->Fill(PTSigmap);
      hccSigmap->Fill(ccSigmap);
    }
    if (ETSigmam!=0){
      hETSigmam->Fill(ETSigmam);
      hPTSigmam->Fill(PTSigmam);
      hccSigmam->Fill(ccSigmam);
    }
    if (ETSigma0!=0){
      hETSigma0->Fill(ETSigma0);
      hPTSigma0->Fill(PTSigma0);
      hccSigma0->Fill(ccSigma0);
    }
    if (ETp!=0){
      hETp->Fill(ETp);
      hPTp->Fill(PTp);
      hccp->Fill(ccp);
    }
    if (ETn!=0){
      hETn->Fill(ETn);
      hPTn->Fill(PTn);
      hccn->Fill(ccn);
    }
    if (ETp_!=0){
      hETp_->Fill(ETp_);
      hPTp_->Fill(PTp_);
      hccp_->Fill(ccp_);
    }
    if (ETn_!=0){
      hETn_->Fill(ETn_);
      hPTn_->Fill(PTn_);
      hccn_->Fill(ccn_);
    }
    if (ETgamma!=0){
      hETgamma->Fill(ETgamma);
      hPTgamma->Fill(PTgamma);
      hccgamma->Fill(ccgamma);
    }
    if (ETother!=0){
      hETother->Fill(ETother);
      hPTother->Fill(PTother);
      hccother->Fill(ccother);
    }
  }// end Event loop
  cout<<"Run complete, Writing to Disk\n";
/**************
Writing to File
**************/
  //to .root
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

  //to .txt
  double TET=pETpip+pETpim+pETpi0+pETKp+pETKm+pETKL+pETKS+pETLambda0+pETLambdaBar0+pETp+pETn+pETp_+pETn_+pETEta+pETOmega+pETOMEGAm+pETXi0+pETXim+pETSigmap+pETSigmam+pETSigma0+pETgamma+pETother;
  double TPT=pPTpip+pPTpim+pPTpi0+pPTKp+pPTKm+pPTKL+pPTKS+pPTLambda0+pPTLambdaBar0+pPTp+pPTn+pPTp_+pPTn_+pPTEta+pPTOmega+pPTOMEGAm+pPTXi0+pPTXim+pPTSigmap+pPTSigmam+pPTSigma0+pPTgamma+pPTother;
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
  out<<"omega:\t"<<cOmega<<"\t"<<pETOmega<<"\t"<<pPTOmega<<"\t"<<(pETOmega/TET)*100<<"\t"<<(pPTOmega/TPT)*100<<endl;
  out<<"Omega-:\t"<<cOMEGAm<<"\t"<<pETOMEGAm<<"\t"<<pPTOMEGAm<<"\t"<<(pETOMEGAm/TET)*100<<"\t"<<(pPTOMEGAm/TPT)*100<<endl;
  out<<"Xi0:\t"<<cXi0<<"\t"<<pETXi0<<"\t"<<pPTXi0<<"\t"<<(pETXi0/TET)*100<<"\t"<<(pPTXi0/TPT)*100<<endl;
  out<<"Xi-:\t"<<cXim<<"\t"<<pETXim<<"\t"<<pPTXim<<"\t"<<(pETXim/TET)*100<<"\t"<<(pPTXim/TPT)*100<<endl;
  out<<"Sigma+:\t"<<cSigmap<<"\t"<<pETSigmap<<"\t"<<pPTSigmap<<"\t"<<(pETSigmap/TET)*100<<"\t"<<(pPTSigmap/TPT)*100<<endl;
  out<<"Sigma-:\t"<<cSigmam<<"\t"<<pETSigmam<<"\t"<<pPTSigmam<<"\t"<<(pETSigmam/TET)*100<<"\t"<<(pPTSigmam/TPT)*100<<endl;
  out<<"Sigma0:\t"<<cSigma0<<"\t"<<pETSigma0<<"\t"<<pPTSigma0<<"\t"<<(pETSigma0/TET)*100<<"\t"<<(pPTSigma0/TPT)*100<<endl;
  out<<"gamma:\t"<<cgamma<<"\t"<<pETgamma<<"\t"<<pPTgamma<<"\t"<<(pETgamma/TET)*100<<"\t"<<(pPTgamma/TPT)*100<<endl;
  out<<endl;
  out<<"[Pions]:\t"<<(cpip+cpim+cpi0)<<"\t"<<(pETpip+pETpim+pETpi0)<<"\t"<<(pPTpip+pPTpim+pPTpi0)<<"\t"<<((pETpip+pETpim+pETpi0)/TET)*100<<"\t"<<((pPTpip+pPTpim+pPTpi0)/TPT)*100<<endl;
  out<<"[Kaons]:\t"<<(cKp+cKm+cKL+cKS)<<"\t"<<(pETKp+pETKm+pETKL+pETKS)<<"\t"<<(pPTKp+pPTKm+pPTKL+pPTKS)<<"\t"<<((pETKp+pETKm+pETKL+pETKS)/TET)*100<<"\t"<<((pPTKp+pPTKm+pPTKL+pPTKS)/TPT)*100<<endl;
  out<<"[Lambdas]:\t"<<(cLambda0+cLambdaBar0)<<"\t"<<(pETLambda0+pETLambdaBar0)<<"\t"<<(pPTLambda0+pPTLambdaBar0)<<"\t"<<((pETLambda0+pETLambdaBar0)/TET)*100<<"\t"<<((pPTLambda0+pPTLambdaBar0)/TPT)*100<<endl;
  out<<"[p,pbar]:\t"<<(cp+cp_)<<"\t"<<(pETp+pETp_)<<"\t"<<(pPTp+pPTp_)<<"\t"<<((pETp+pETp_)/TET)*100<<"\t"<<((pPTp+pPTp_)/TPT)*100<<endl;
  out<<"[n,nbar]:\t"<<(cn+cn_)<<"\t"<<(pETn+pETn_)<<"\t"<<(pPTn+pPTn_)<<"\t"<<((pETn+pETn_)/TET)*100<<"\t"<<((pPTn+pPTn_)/TPT)*100<<endl;
  out<<"[eta,ome]:\t"<<(cEta+cOmega)<<"\t"<<(pETEta+pETOmega)<<"\t"<<(pPTEta+pPTOmega)<<"\t"<<((pETEta+pETOmega)/TET)*100<<"\t"<<((pPTEta+pPTOmega)/TPT)*100<<endl;
  out<<"[Xi,Om,Sig]:\t"<<(cOMEGAm+cXi0+cXim+cSigmap+cSigmam+cSigma0)<<"\t"<<(pETOMEGAm+pETXi0+pETXim+pETSigmap+pETSigmam+pETSigma0)<<"\t"<<(pPTOMEGAm+pPTXi0+pPTXim+pPTSigmap+pPTSigmam+pPTSigma0)<<"\t"<<((pETOMEGAm+pETXi0+pETXim+pETSigmap+pETSigmam+pETSigma0)/TET)*100<<"\t"<<((pPTOMEGAm+pPTXi0+pPTXim+pPTSigmap+pPTSigmam+pPTSigma0)/TPT)*100<<endl;
  out<<"[gamma]:\t"<<cgamma<<"\t"<<pETgamma<<"\t"<<pPTgamma<<"\t"<<(pETgamma/TET)*100<<"\t"<<(pPTgamma/TPT)*100<<endl;
  out<<"[exotic]:\t"<<cother<<"\t"<<pETother<<"\t"<<pPTother<<"\t"<<(pETother/TET)*100<<"\t"<<(pPTother/TPT)*100<<endl;
  out<<endl;
  out<<"counted ET:\t"<<TET<<"\tTotal ET:\t"<<pETALL<<"\t(counted/total)%:\t"<<(TET/pETALL)*100<<endl;
  out<<endl;
  out<<"Decay Products vs 'Primary'\nGroup\tET\tPT\tORIGINAL: ETo\tPTo\tET%\tPT%\n";
  out<<"[Pions (eta)]\t"<<ETetapis<<"\t"<<PTetapis<<"\t"<<(pETpip+pETpim+pETpi0)<<"\t"<<(pPTpip+pPTpim+pPTpi0)<<"\t"<<(ETetapis/(pETpip+pETpim+pETpi0))*100<<"\t"<<(PTetapis/(pPTpip+pPTpim+pPTpi0))*100<<endl;
  out<<"[Pions (omega)]\t"<<ETomegapis<<"\t"<<PTomegapis<<"\t"<<(pETpip+pETpim+pETpi0)<<"\t"<<(pPTpip+pPTpim+pPTpi0)<<"\t"<<(ETomegapis/(pETpip+pETpim+pETpi0))*100<<"\t"<<(PTomegapis/(pPTpip+pPTpim+pPTpi0))*100<<endl;

  out.close();
  out2.close();
  outfile->Close();
  return 0;
}//end makeEventSample

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

void SimplePYTHIALoop(Int_t n=1000, Int_t jobID=0, Int_t tune = 350,Float_t SNN = 2760,char CUT='y',Float_t yncut=0.1,char DTYPE='r',char DMODE='a',Float_t trigEtaMax = 0.5, Float_t assocEtaMax = 0.9) {
  makeEventSample(n,jobID,tune,SNN,CUT,yncut,DTYPE,DMODE, trigEtaMax, assocEtaMax);
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
    retVal = makeEventSample(n,0,0,7.7,'y',0.1,'r','a',0.5,0.9);
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

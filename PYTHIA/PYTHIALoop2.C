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
//#include <HeavyIons.h>
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
  TH1F *histo = new TH1F(name,name,1000,0,20);
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
int makeEventSample(Int_t nEvents, Int_t jobID, Int_t tune, Float_t SNN, char CUT, Float_t yncut, char DTYPE, Float_t trigEtaMax = 0.5, Float_t assocEtaMax = 0.9)
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
  else if (SNN==62.4){
    modSNN=62;
  }
  else {
    modSNN=SNN;
  }
  cout<<modSNN<<endl;
  cout<<"I made it here, running "<<nEvents<<" events, job id "<<jobID<<", tune "<<tune<<", at "<<SNN<<"GeV"<<endl;
  char *filename = Form("hh%iGeVoutfile%i.root",modSNN,jobID);
  char *filename2 = Form("hh%iGeVoutfile%i.txt",modSNN,jobID);
  ifstream in;
  in.open("KF_Code.dat");
  while (in.good()){
    partname.push_back(KF_Code());
    in>>partname[l].KFc>>partname[l].name;
    l++;
  }
  ofstream out;
  out.open(filename2);

  cout<<"I made it here line 144"<<endl;
  // Create an instance of the Pythia event generator ...
  string pVERSION;
  TPythia6* pythia = new TPythia6();
  pVERSION="Pythia6";
  //TPythia8* pythia = new TPythia8();
  //pVERSION="Pythia8";
  //Angantyr* pythia = new Angantyr();
  //pVERSION="Angantyr";

  cout<<"I made it here line 157"<<endl;

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
  out<<temporarything<<"\nCut: "<<yncut<<"\nCalorimeter (r=regular, c=calorimeter): "<<DTYPE<<endl<<endl;

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
  Int_t cLambda0=0;
  Int_t cLambdaBar0=0;
  Int_t cp=0;
  Int_t cn=0;
  Int_t cp_=0;
  Int_t cn_=0;

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



  for (Int_t event = 0; event < nEvents; event++) {
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

      vector<pylista> pylis(1);
    // Show how far we got every 100'th event.
    if (event % 10 == 0)
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
      for (Int_t parta=0; parta<npart; parta++) {
        pylis.push_back(pylista());
        pylis[parta].inde=parta;
        pylis[parta].KFpart=pythia->GetK(parta,2);
        pylis[parta].indePar=pythia->GetK(parta,3);
      }
      //printf("Analyse %d Particles\n", npart);
      for (Int_t part=0; part<npart; part++) {
	//TObject *object = particles->At(part);
	//cout<<"I am a "<<object->ClassName()<<endl;
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
proton           ALL INCLUDED
neutron          ALL INCLUDED

For reference, mpP is in KF format, Ckf is converted to index (from KF_Code.dat)
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
        //cout<<"Omega "<<ETOmega<<endl;
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
  }
  Double_t ETAll=0;
  Double_t ETpip=0;
  Double_t ETpim=0;
  Double_t ETpi0=0;
  Double_t ETKp=0;
  Double_t ETKm=0;
  Double_t ETKL=0;
  Double_t ETKS=0;
  Double_t ETEta=0;
  Double_t ETOmega=0; //THIS IS omega not Omega
  Double_t ETOmegaM=0;
  Double_t ETLambda0=0;
  Double_t ETLambdaBar0=0;
  Double_t ETp=0;
  Double_t ETn=0;
//NOTE: turn into pt graphs for pi0 and eta (and omega), normalize by multiplying by 1/2pi divide histos

/*hPTpip->Scale((1.)/(2*TMath::Pi())); //n is bin width as histogram
hPTpim->Scale((1.)/(2*TMath::Pi()));
hPTpi0->Scale((1.)/(2*TMath::Pi()));
hPTKp->Scale((1.)/(2*TMath::Pi()));
hPTKm->Scale((1.)/(2*TMath::Pi()));
hPTKL->Scale((1.)/(2*TMath::Pi()));
hPTKS->Scale((1.)/(2*TMath::Pi()));
hPTEta->Scale((1.)/(2*TMath::Pi()));
hPTOmega->Scale((1.)/(2*TMath::Pi()));
hPTLambda0->Scale((1.)/(2*TMath::Pi()));
hPTLambdaBar0->Scale((1.)/(2*TMath::Pi()));
hPTp->Scale((1.)/(2*TMath::Pi()));
hPTn->Scale((1.)/(2*TMath::Pi()));
hPTp_->Scale((1.)/(2*TMath::Pi()));
hPTn_->Scale((1.)/(2*TMath::Pi()));
*/
//just use raw number over other, dont need pt spectra

  //NOTE:write histograms to file
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

    out<<"Pion Species Data:\n**********************\n";
    out<<"Pi+:\n";
    out<<"Number: "<<cpip<<"\tTotal Et: "<<pETpip<<"\tTotal Pt: "<<pPTpip<<endl;
    out<<"Pi-:\n";
    out<<"Number: "<<cpim<<"\tTotal Et: "<<pETpim<<"\tTotal Pt: "<<pPTpim<<endl;
    out<<"Pi0:\n";
    out<<"Number: "<<cpi0<<"\tTotal Et: "<<pETpi0<<"\tTotal Pt: "<<pPTpi0<<endl;
    out<<"\nKaon Species Data:\n**********************\n";
    out<<"K+:\n";
    out<<"Number: "<<cKp<<"\tTotal Et: "<<pETKp<<"\tTotal Pt: "<<pPTKp<<endl;
    out<<"K-:\n";
    out<<"Number: "<<cKm<<"\tTotal Et: "<<pETKm<<"\tTotal Pt: "<<pPTKm<<endl;
    out<<"K0L:\n";
    out<<"Number: "<<cKL<<"\tTotal Et: "<<pETKL<<"\tTotal Pt: "<<pPTKL<<endl;
    out<<"K0S:\n";
    out<<"Number: "<<cKS<<"\tTotal Et: "<<pETKS<<"\tTotal Pt: "<<pPTKS<<endl;
    out<<"\nLambda Species Data:\n**********************\n";
    out<<"Lambda0:\n";
    out<<"Number: "<<cLambda0<<"\tTotal Et: "<<pETLambda0<<"\tTotal Pt: "<<pPTLambda0<<endl;
    out<<"LambdaBar0:\n";
    out<<"Number: "<<cLambdaBar0<<"\tTotal Et: "<<pETLambdaBar0<<"\tTotal Pt: "<<pPTLambdaBar0<<endl;
    out<<"\nProton/Neutron Species Data:\n**********************\n";
    out<<"Proton:\n";
    out<<"Number: "<<cp<<"\tTotal Et: "<<pETp<<"\tTotal Pt: "<<pPTp<<endl;
    out<<"Neutron:\n";
    out<<"Number: "<<cn<<"\tTotal Et: "<<pETn<<"\tTotal Pt: "<<pPTn<<endl;
    out<<"Antiproton:\n";
    out<<"Number: "<<cp_<<"\tTotal Et: "<<pETp_<<"\tTotal Pt: "<<pPTp_<<endl;
    out<<"Antineutron:\n";
    out<<"Number: "<<cn_<<"\tTotal Et: "<<pETn_<<"\tTotal Pt: "<<pPTn_<<endl;
    out<<"\nExotic Species Data:\n**********************\n";
    out<<"eta:\n";
    out<<"Number: "<<cEta<<"\tTotal Et: "<<pETEta<<"\tTotal Pt: "<<pPTEta<<endl;
    out<<"omega:\n";
    out<<"Number: "<<cOmega<<"\tTotal Et: "<<pETOmega<<"\tTotal Pt: "<<pPTOmega<<endl;
    out.close();
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

void SimplePYTHIALoop(Int_t n=1000, Int_t jobID=0, Int_t tune = 350,Float_t SNN = 2760,char CUT='y',Float_t yncut=0.1,char DTYPE='r',Float_t trigEtaMax = 0.5, Float_t assocEtaMax = 0.9) {
  makeEventSample(n,jobID,tune,SNN,CUT,yncut,DTYPE, trigEtaMax, assocEtaMax);
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
    retVal = makeEventSample(n,0,0,7.7,'y',0.1,'r',0.5,0.9);
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

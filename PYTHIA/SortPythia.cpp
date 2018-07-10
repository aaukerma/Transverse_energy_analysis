/********************************************
SortPythia.cpp
designed to syphon data from pythia output files
Start:20180626.1135
Completed:--
**********************************************/

#include "Riostream.h"
#include <cstdio>
#include <vector>
#include <string>
#include <fstream>
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <iostream>
using namespace std;

struct Species {
  Int_t KF;
  Int_t num;
  double frac;
  string PN;
};

struct KF_Code {
  int KFc;
  string name;
};

string GetPName(const int kfc, const vector<KF_Code>& partname);

void SortPythia(){
  int k=0;
  int l=0;
  double partsum=0;
  double KFCode;
  double temp;
  string name;
  Int_t temp1;
  vector <Species> partnum(0);
  vector <KF_Code> partname(0);

  cout<<"start"<<endl;
  ifstream in;
  in.open("KF_Code.dat");
  ofstream file;
  file.open("./PythiaData3.txt");
  TFile f("hh_outfile3.root");
  if (f.IsZombie()) {
    cout<<"error opening file"<<endl;
    exit(-1);
  }
  file<<"KF\t"<<"PName\t"<<"quantity\t"<<"percent"<<endl;
  TH1F *trigs= (TH1F*)f.Get("trigs");
  TH1F *hSPALL= (TH1F*)((TH1F*)f.FindObjectAny("hSPALL"))->Clone();

  while (in.good()){
    partname.push_back(KF_Code());
    in>>partname[l].KFc>>partname[l].name;
    l++;
  }
  for(int i=-9910553;i<9910553;i++){
    KFCode=hSPALL->GetBinContent(hSPALL->FindBin(i));
    if (KFCode!=0){
      partnum.push_back(Species());
      partnum[k].KF=i;
      partnum[k].num=KFCode;
      name=GetPName(i,partname);
      partnum[k].PN=name;
      if ((i>100)||(i<-100)){ //NOTE:may include bad data!! TODO:exclude all bad numbers
      partsum=partsum+partnum[k].num;
      }
      k++;
    }
  }
  for (int i=0;i<k;i++){
    temp1=partnum[i].KF;
    if ((temp1>100)||(temp1<-100)){
      temp=partnum[i].num;
      partnum[i].frac = (temp/partsum)*100;
      file<<partnum[i].KF<<" \t"<<partnum[i].PN<<"\t"<<partnum[i].num<<" \t"<<partnum[i].frac<<endl;
    }
  }
  file<<endl<<endl<<"Total particles counted: "<<partsum<<endl<<endl;

  //from christine:[
  /*
  for all events{
    ETpiplus=0;
    ETall=0;
    for all particles{
      if(piplus) etpiplus+=Etparticle;
      ETall += Etparticle;
    }
    histopiplus->Fill(etpiplus);
    histoall->Fill(ETall);
  }
  */

  //]
  file.close();
  cout<<"done"<<endl;
  return;
}

string GetPName(const int kfc, const vector<KF_Code>& partname){
  int i=0;
  int j=0;
  int k=0;
  int n=1;
  int KF;
  int comp;
  int comp1;
  string name;
  KF=kfc;

  if (KF<0){
    KF=KF*-1;
    n=-1;
  }
  for (int i=0;i<248;i++){
    comp=partname[i].KFc;
    if(comp==KF){
      j=comp;
      k=i;
      break;
    }
  }
  //NOTE: to include if -KF matters
//  if (n==-1){
//    name=partname[k].name+" BAR";
//  }
//  else
    name=partname[k].name;
  return name;
}

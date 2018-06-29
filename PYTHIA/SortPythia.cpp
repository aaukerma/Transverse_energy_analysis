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
#include <iostream>
using namespace std;

//TFile::TFile(const char* hh_outfile2.root,Option_t* option="",const char* "hh_outfile2", int_t compress=0)
int SortPythia(){
  cout<<"start"<<endl;
  ofstream file;
  file.open("./PythiaData.txt");
  TFile f("hh_outfile2.root");
  if (f.IsZombie()) {
    cout<<"error opening file"<<endl;
    exit(-1);
  }
  TH1F *trig= (TH1F*)f.Get("trig");
  //hh->THnSparse::GetNbins();  works in root
  //hh->THnSparse::GetBinContent(0, 0)
  //see THnSparse on interwebs
  TH1F *axis0= (TH1F*)f.Get("axis0");
  TH1F *axis1= (TH1F*)f.Get("axis1");
  TH1F *axis2= (TH1F*)f.Get("axis2");
  TH1F *axis3= (TH1F*)f.Get("axis3");
  TH1F *axis4= (TH1F*)f.Get("axis4");

file.close();
  cout<<"done"<<endl;



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
  return 0;
}

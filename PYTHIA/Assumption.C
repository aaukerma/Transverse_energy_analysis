/******************************************************************************
This program takes specified output runs and makes graphs to show the validity
of the assumptions made. These assumptions include the pion0 assumption, the
kaon0 (long and short) assumptions, the neutron and antineutron assumptions,
as well as the eta and omega assumptions.
******************************************************************************/

#include "Riostream.h"
#include <cstdio>
#include <vector>
#include <string>
#include <fstream>
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1F.h"
#include "TCanvas.h"
#include <iostream>
using namespace std;

int Assumption(){
  Int_t IDNUM=1;
  char* filename = Form("hhAnalysisOfRun%i",IDNUM);
  vector<Int_t> SNN {7,11,19,27,39,62,130,200,900,2760,5020,7000,8000,13000,14000};
  vector<float> SNNACTUAL {7.7,11.6,19.6,27,39,62.4,130,200,900,2760,5020,7000,8000,13000,14000};
  TFile* file = TFile::Open(filename, "RECREATE");
  for (int x;x<14;x++){
    char* temp = Form("hh%iGeVoutfile%i.root",SNN[x],IDNUM);
    TFile f(temp);
    if (f.IsZombie()) {
        cout<<"error opening file"<<endl;
        exit(-1);
    }
    temp=Form("hETPi0Check%i",SNN[x]);
    TH1F* (temp)= (TH1F*)((TH1F*)f.FindObjectAny("hETPi0Check1"))->Clone();

    f.Close();
  }





































}

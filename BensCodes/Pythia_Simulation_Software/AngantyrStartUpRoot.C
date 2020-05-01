void AngantyrStartUpRoot(){
  gSystem->Load("libEG");
  gSystem->Load("libEGPythia8");
  gSystem->Load("libpythia8");
  cout<<"\n\n\nSTART PROGRAM\n\n\n";

  int EVN=100000;                //number of events to run
  int ID=20200430;              //jobID (recommend date in format YYYYMMDD)
  char CUT='n';                 //Type of cut (rapidity 'y' or pseudorapidity 'n')
  float YNCUT=0.1;              //size of above cut
  char DTYPE='r';               //regular or calrimeter settings (see README)
  char DMODE='d';               //Counting method employed, (see README) DO NOT CHANGE!!!!!!

  /*************************
      Collision Energy
  SNN values for RHIP Group:
            7.7*
            11.5
            19.6
            27
            39
            62.4
            130
            200
            2760
            5020
  (program limit: 10GeV and up)
   *7.7GeV returns 10GeV BUG
  ***************************/
  float SNN = 200;

  gROOT->ProcessLine(".L  Angantyr.C");
  char* event1= Form("makeEventSample(%i,%i,%f,'%c',%f,'%c','%c')",EVN,ID,SNN,CUT,YNCUT,DTYPE,DMODE);
  gROOT->ProcessLine(event1);
}

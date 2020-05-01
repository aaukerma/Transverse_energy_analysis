void PythiaStartUpRoot(){
  gSystem->Load("libEG");
  gSystem->Load("libEGPythia6");
  gSystem->Load("libPythia6");
  int EVN=1;
  int ID=20200117;
  int TUNE=350;
  char CUT='y'; //
  float YNCUT=0.1;
  char DTYPE='r';
  char DMODE='d';

  /*************************
      Collision Energy
  SNN values for RHIP Group:
            7.7
            11.5
            19.6
            27
            39
            62.4
            130
            200
            2760
            5020
  ***************************/
  float SNN = 200;


  gROOT->ProcessLine(".L Pythia6.C");
  char* event1= Form("makeEventSample(%i,%i,%i,%f,'%c',%f,'%c','%c')",EVN,ID,TUNE,SNN,CUT,YNCUT,DTYPE,DMODE);
  gROOT->ProcessLine(event1);

}

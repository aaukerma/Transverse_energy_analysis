void PAngantyrStartUpRoot(){
  gSystem->Load("libEG");
  gSystem->Load("libEGPythia8");
  gSystem->Load("libpythia8");
  int EVN=500;
  int ID=20190730;
  char CUT='y'; //
  float YNCUT=0.1;
  char DTYPE='r';
  char DMODE='d';
  gROOT->ProcessLine(".L PAngantyr.C");
  char* event1= Form("makeEventSample(%i,%i,200,'%c',%f,'%c','%c')",EVN,ID,CUT,YNCUT,DTYPE,DMODE);
  gROOT->ProcessLine(event1);
}

void PAngantyrStartUpRoot(){
  gSystem->Load("libEG");
  gSystem->Load("libEGPythia8");
  gSystem->Load("libpythia8");
  int EVN=5;
  int ID=20190723;
  char CUT='y'; //
  float YNCUT=0.1;
  char DTYPE='r';
  char DMODE='c';
  gROOT->ProcessLine(".L PAngantyr.C");
  char* event1= Form("makeEventSample(%i,%i,200,'%c',%f,'%c','%c')",EVN,ID,CUT,YNCUT,DTYPE,DMODE);
  gROOT->ProcessLine(event1);
}

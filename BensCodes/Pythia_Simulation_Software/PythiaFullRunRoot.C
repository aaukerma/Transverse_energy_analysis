void PythiaFullRunRoot(){
  gSystem->Load("libEG");
  gSystem->Load("libEGPythia6");
  gSystem->Load("libPythia6");
  int EVN=100000;
  int HIGHEVN=10000;
  int ID=20190905;
  char CUT='y'; //
  float YNCUT=0.1;
  char DTYPE='r';
  char DMODE='d';
  char* event1= Form("makeEventSample(%i,%i,350,7.7,'%c',%f,'%c','%c')",EVN,ID,CUT,YNCUT,DTYPE,DMODE);
  char* event2= Form("makeEventSample(%i,%i,350,11.6,'%c',%f,'%c','%c')",EVN,ID,CUT,YNCUT,DTYPE,DMODE);
  char* event3= Form("makeEventSample(%i,%i,350,19.5,'%c',%f,'%c','%c')",EVN,ID,CUT,YNCUT,DTYPE,DMODE);
  char* event4= Form("makeEventSample(%i,%i,350,27,'%c',%f,'%c','%c')",EVN,ID,CUT,YNCUT,DTYPE,DMODE);
  char* event5= Form("makeEventSample(%i,%i,350,39,'%c',%f,'%c','%c')",EVN,ID,CUT,YNCUT,DTYPE,DMODE);
  char* event6= Form("makeEventSample(%i,%i,350,62.4,'%c',%f,'%c','%c')",EVN,ID,CUT,YNCUT,DTYPE,DMODE);
  char* event7= Form("makeEventSample(%i,%i,350,130,'%c',%f,'%c','%c')",EVN,ID,CUT,YNCUT,DTYPE,DMODE);
  char* event8= Form("makeEventSample(%i,%i,350,200,'%c',%f,'%c','%c')",EVN,ID,CUT,YNCUT,DTYPE,DMODE);
  char* event10= Form("makeEventSample(%i,%i,350,2760,'%c',%f,'%c','%c')",HIGHEVN,ID,CUT,YNCUT,DTYPE,DMODE);
  char* event11= Form("makeEventSample(%i,%i,350,5020,'%c',%f,'%c','%c')",HIGHEVN,ID,CUT,YNCUT,DTYPE,DMODE);
  gROOT->ProcessLine(".L Pythia6.C")
  gROOT->ProcessLine(event1);
  gROOT->ProcessLine(event2);
  gROOT->ProcessLine(event3);
  gROOT->ProcessLine(event4);
  gROOT->ProcessLine(event5);
  gROOT->ProcessLine(event6);
  gROOT->ProcessLine(event7);
  gROOT->ProcessLine(event8);
  gROOT->ProcessLine(event10);
  gROOT->ProcessLine(event11);
}

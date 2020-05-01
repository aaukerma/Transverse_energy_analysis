void AngantyrFullRunRoot(){
  gSystem->Load("libEG");
  gSystem->Load("libEGPythia8");//gSystem->Load("$HOME/alice/sw/ubuntu1804_x86-64/AliRoot/v5-09-47c-1/lib/libpythia6_4_28");
  gSystem->Load("libpythia8");
  int EVN=100000;
  int HIGHEVN=10000;
  int ID=20190918;
  char CUT='y'; //
  float YNCUT=0.1;
  char DTYPE='r';
  char DMODE='d';
  gROOT->ProcessLine(".L Angantyr_4.C");
  char* event1= Form("makeEventSample(%i,%i,7.7,'%c',%f,'%c','%c')",EVN,ID,CUT,YNCUT,DTYPE,DMODE);
  char* event2= Form("makeEventSample(%i,%i,11.6,'%c',%f,'%c','%c')",EVN,ID,CUT,YNCUT,DTYPE,DMODE);
  char* event3= Form("makeEventSample(%i,%i,19.5,'%c',%f,'%c','%c')",EVN,ID,CUT,YNCUT,DTYPE,DMODE);
  char* event4= Form("makeEventSample(%i,%i,27,'%c',%f,'%c','%c')",EVN,ID,CUT,YNCUT,DTYPE,DMODE);
  char* event5= Form("makeEventSample(%i,%i,39,'%c',%f,'%c','%c')",EVN,ID,CUT,YNCUT,DTYPE,DMODE);
  char* event6= Form("makeEventSample(%i,%i,62.4,'%c',%f,'%c','%c')",EVN,ID,CUT,YNCUT,DTYPE,DMODE);
  char* event7= Form("makeEventSample(%i,%i,130,'%c',%f,'%c','%c')",EVN,ID,CUT,YNCUT,DTYPE,DMODE);
  char* event8= Form("makeEventSample(%i,%i,200,'%c',%f,'%c','%c')",EVN,ID,CUT,YNCUT,DTYPE,DMODE);
  char* event10= Form("makeEventSample(%i,%i,2760,'%c',%f,'%c','%c')",HIGHEVN,ID,CUT,YNCUT,DTYPE,DMODE);
  char* event11= Form("makeEventSample(%i,%i,5020,'%c',%f,'%c','%c')",HIGHEVN,ID,CUT,YNCUT,DTYPE,DMODE);
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

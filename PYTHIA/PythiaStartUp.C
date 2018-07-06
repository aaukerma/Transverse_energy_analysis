void PythiaStartUp(){
  gSystem->Load("libEG");
  gSystem->Load("libEGPythia6");
  gSystem->Load("$HOME/alice/sw/ubuntu1604_x86-64/AliRoot/sw-1/lib/libpythia6");
  //gSystem->Load("libpythia6");
  gSystem->Load("libAliPythia6");
  gROOT->ProcessLine(".L SimplePYTHIALoop.C");
}

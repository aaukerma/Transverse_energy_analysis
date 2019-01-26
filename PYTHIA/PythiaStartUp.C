void PythiaStartUp(){
  gSystem->Load("libEG");
  gSystem->Load("libEGPythia6");
  gSystem->Load("$HOME/alice/sw/ubuntu1804_x86-64/AliRoot/master-1/lib/libpythia6_4_28");
  //gSystem->Load("libpythia6");
  gSystem->Load("libAliPythia6");
  gROOT->ProcessLine(".L PYTHIALoop2.C");
}

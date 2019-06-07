void PythiaStartUpRoot(){
  gSystem->Load("libEG");
  gSystem->Load("libEGPythia6");
  //gSystem->Load("$HOME/alice/sw/ubuntu1804_x86-64/AliRoot/v5-09-47c-1/lib/libpythia6_4_28");
  gSystem->Load("libPythia6");
  //gSystem->Load("libAliPythia6");
  gROOT->ProcessLine(".L PYTHIALoop2.C");
}

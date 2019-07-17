void PAngantyrStartUpRoot(){
  gSystem->Load("libEG");
  gSystem->Load("libEGPythia8");
  gSystem->Load("libpythia8");
  gROOT->ProcessLine(".L PAngantyr.C");
}

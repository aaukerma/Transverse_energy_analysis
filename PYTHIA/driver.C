void driver(Int_t nEvents = 100, Int_t jobID = 0, Int_t tune = 100){
  gSystem->Load("$PYTHIA6/libPythia6"); //change to your setup
  gSystem->Load("libEGPythia6");
  gROOT->ProcessLine(".L SimplePYTHIALoop.C++");
  makeEventSample(nEvents,jobID,tune);

}

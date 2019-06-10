Int_t n = 3; // This makes n "visible" to the compiled macro
void rungen200_list(Int_t nev=3){
  // Simulation and reconstruction
  TStopwatch timer;
  timer.Start();
  gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT -I$ALICE_ROOT/EVGEN");
  gSystem->Load("liblhapdf");      // Parton density functions
  gSystem->Load("libEGPythia6");   // TGenerator interface
  gSystem->Load("libpythia6");     // Pythia
  gSystem->Load("libAliPythia6");  // ALICE specific implementations
  n = nev; // Use the requested number of events
  gROOT->Macro("fastMcProduction200_list.C++(n)");
  timer.Stop();
  timer.Print();
}

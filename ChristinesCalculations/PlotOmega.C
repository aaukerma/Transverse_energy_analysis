int nenergies = 8;
double energies[] = {7.7,11.5,19.6,27,39, 62.4,130,200};
double energiesMaster[] = {7.7,11.5,19.6,27,39, 62.4,130,200};
double energiesUnc[10] = {0,0,0,0,0,0,0,0,0,0};
Int_t markers[] = {20,21,22,23,29,33,34,20,20,20};
// Int_t colors[] = {kViolet-6,kBlue+2,kBlue,kCyan+1,
// 		  kTeal-6,kGreen+2, kOrange+1,kRed,kPink,1,1};
Int_t colors[] = {kRed,kOrange+1,kGreen+2,kTeal-6,kCyan+1,kBlue,kBlue+2,kViolet-6,kViolet+3,kBlue+2,kBlue,kCyan+1,
		  kTeal-6,kGreen+2, kOrange+1,kPink,1,1};

double eta[10][10]={{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0}};
double etaerr[10][10]={{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0}};
double omega[10][10]={{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0}};
double omegaerr[10][10]={{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0}};
int ncb = 9;
TGraphErrors *GetsNNTGraph(Double_t *y,Double_t *ey,Int_t color, Int_t marker,Int_t fillstyle = 0){
  TGraphErrors *gr = new TGraphErrors(8,(Double_t *)energies,(Double_t *)y,(Double_t *)energiesUnc,(Double_t *)ey);
  gr->SetMarkerColor(color);
  gr->SetLineColor(color);
  gr->SetFillColor(color);
  gr->SetMarkerStyle(marker);
  gr->SetMarkerSize(3.0);
  gr->SetFillStyle(fillstyle);
  gr->SetLineWidth(2);
  gr->GetHistogram()->GetXaxis()->SetTitle("#sqrt{s_{NN}}");
  gr->GetHistogram()->GetYaxis()->SetTitle("E_{T}^{#eta,#omega}/E_{T}^{#pi^0}");
  gr->GetHistogram()->GetXaxis()->SetTitleSize(0.07);
  gr->GetHistogram()->GetXaxis()->SetLabelSize(0.07);
  gr->GetHistogram()->GetYaxis()->SetTitleSize(0.07);
  gr->GetHistogram()->GetYaxis()->SetTitleOffset(0.8);
  gr->GetHistogram()->GetYaxis()->SetLabelSize(0.07);
  return gr;
}

void PlotOmega(){
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  char *MYFILENAME = "EtRatios.txt";
  ifstream file(MYFILENAME);
  double energy = 7.7;
  int energyN = 0;
  double myenergy, cent, etafrac, etafracerr, omegafrac, omegafracerr;
  float lowOmega = 1;
  float highOmega = -1;
  float lowEta = 1;
  float highEta = -1;
  while(!file.eof()){
    file>>myenergy>>cent>>etafrac>>etafracerr>>omegafrac>>omegafracerr;
    if(myenergy>energy){
      energyN++;
      energy = myenergy+0.01;
    }
    int cb = (int) cent;
    cout<<"e bin "<<energyN<<" cent "<<cent<<endl;
    eta[cb][energyN] = etafrac;
    etaerr[cb][energyN] = etafracerr;
    omega[cb][energyN] = 1.0*omegafrac;
    omegaerr[cb][energyN] = 1.0*omegafracerr;
    if(etafrac+etafracerr>highEta) highEta = etafrac+etafracerr;
    if(etafrac-etafracerr<lowEta) lowEta = etafrac-etafracerr;
    if(omegafrac+omegafracerr>highOmega) highOmega = omegafrac+omegafracerr;
    if(omegafrac-omegafracerr<lowOmega) lowOmega = omegafrac-omegafracerr;
    if(myenergy>100 && myenergy<150){//130 GeV only has 7 bins so we're just going to make this one off the plot
      eta[cb][energyN] = -1;
      omega[cb][energyN] = -1;
    }
    //for(int j=0;j<ncb;j++){
    //for(int i=0;i<nenergies;i++){
	//filenpart>>npartPHENIX[i][j]>>nparterrPHENIX[i][j];
    //}
      //}
  }
  cout<<"eta "<<lowEta<<" - "<<highEta<<endl;
  cout<<"omega "<<lowOmega<<" - "<<highOmega<<endl;
  float xlow = 5;
  float xhigh = 250;
  TBox *boxEta = new TBox(xlow,lowEta,xhigh,highEta);
  boxEta->SetFillStyle(1001);
  boxEta->SetFillColor(kGray);
  boxEta->SetLineWidth(1);
  boxEta->SetLineStyle(1);
  TLine *lineEta = new TLine(xlow,(lowEta+highEta)/2.0,xhigh,(lowEta+highEta)/2.0);
  lineEta->SetLineWidth(2);
   TCanvas *c1 = new TCanvas("c1", "eta/pi0",1204,928/2.0);
   c1->Range(-1.768938,-0.3125,15.92044,2.8125);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);
   c1->SetLeftMargin(0);//0.168885);//0.1425);
   c1->SetRightMargin(0);//0.015);
   c1->SetTopMargin(0);//0.0321508);
   c1->SetBottomMargin(0);//0.160754);
   c1->SetLogx();
   c1->Divide(2);
   TPad *pad2 =(TPad*) c1->cd(2);
   pad2->SetLogx();
   pad2->SetTopMargin(0.0245826);
   pad2->SetBottomMargin(0.168367);
   pad2->SetLeftMargin(0.005165);
   pad2->SetRightMargin(0.195819);//0.166355);
   TPad *pad1 = (TPad*)c1->cd(1);
   pad1->SetLogx();
   pad1->SetTopMargin(0.0245826);
   pad1->SetBottomMargin(0.168367);
   pad1->SetRightMargin(0.005165);
   pad1->SetLeftMargin(0.166355);
   TH1F *frame = new TH1F("frame","frame",200,xlow,xhigh);
   frame->SetMinimum(0.0);
   frame->SetMaximum(0.25);
  frame->GetXaxis()->SetTitle("#sqrt{s_{NN}}");
  frame->GetYaxis()->SetTitle("E_{T}^{#eta}/E_{T}^{#pi^{0}}");
  frame->GetXaxis()->SetTitleSize(0.07);
  frame->GetXaxis()->SetLabelSize(0.07);
  frame->GetYaxis()->SetTitleSize(0.07);
  frame->GetYaxis()->SetTitleOffset(1.2);
  frame->GetYaxis()->SetLabelSize(0.07);
   TH1F *frame2 = new TH1F("frame2","frame2",200,xlow,xhigh);
   frame2->SetMinimum(0.0);
   frame2->SetMaximum(0.055);
  frame2->GetXaxis()->SetTitle("#sqrt{s_{NN}}");
  frame2->GetYaxis()->SetTitle("E_{T}^{#omega}/E_{T}^{#pi^{0}}");
  frame2->GetXaxis()->SetTitleSize(0.07);
  frame2->GetXaxis()->SetLabelSize(0.07);
  frame2->GetYaxis()->SetTitleSize(0.07);
  frame2->GetYaxis()->SetTitleOffset(1.4);
  frame2->GetYaxis()->SetLabelSize(0.07);
   TGraphErrors *graphEta[9] = {NULL,NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL};
   TGraphErrors *graphOmega[9] = {NULL,NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL};
  for(int cb = 0; cb<9;cb++){
    for(int e=0;e<nenergies;e++) energies[e] = energiesMaster[e]+cb*TMath::Log(0.1*energiesMaster[e]);
    graphEta[cb] = GetsNNTGraph(eta[cb],etaerr[cb],colors[cb],markers[cb]);
    graphOmega[cb] = GetsNNTGraph(omega[cb],omegaerr[cb],colors[cb],markers[cb]);
  }
  graphEta[0]->SetHistogram(frame);
  graphOmega[0]->SetHistogram(frame2);
  for(int cb = 0; cb<9;cb++){
    if(cb==0){
      graphEta[cb]->Draw("AP");
      //boxEta->Draw("lf");
      //lineEta->Draw();
    }
    graphEta[cb]->Draw("P same");
  }
  pad2->cd();
  for(int cb = 0; cb<9;cb++){
    if(cb==0)    graphOmega[cb]->Draw("AP Y+");
    else{     graphOmega[cb]->Draw("P same");}
  }
  TLegend *leg = new TLegend(0.00343178,0.643785,0.272081,0.980056);//0.695285,0.03239940.6272879,0.5570321,0.8943428,0.8925803,NULL,"brNDC");
   leg->SetBorderSize(1);
   leg->SetTextFont(72);
   leg->SetTextSize(0.0443459);
   leg->SetLineColor(0);
   leg->SetLineStyle(0);
   leg->SetLineWidth(0);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   leg->AddEntry(graphOmega[0],"0-5% (0-12%)","P");
   leg->AddEntry(graphOmega[1],"5-10%","P");
   leg->AddEntry(graphOmega[2],"10-20%","P");
   leg->AddEntry(graphOmega[3],"20-30%","P");
   leg->AddEntry(graphOmega[4],"30-40%","P");
   //leg->AddEntry(graphOmega[5],"40-50%","P");
   //leg->AddEntry(graphOmega[6],"50-60%","P");
   //leg->AddEntry(graphOmega[7],"60-70%","P");
   //leg->AddEntry(graphOmega[8],"70-80%","P");
   leg->Draw();
   TLegend *leg2 = new TLegend(0.263415,0.727273,0.532065,0.980056);//0.6272879,0.5570321,0.8943428,0.8925803,NULL,"brNDC");
   leg2->SetBorderSize(1);
   leg2->SetTextFont(72);
   leg2->SetTextSize(0.0443459);
   leg2->SetLineColor(0);
   leg2->SetLineStyle(0);
   leg2->SetLineWidth(0);
   leg2->SetFillColor(0);
   leg2->SetFillStyle(0);
   //leg2->AddEntry(graphOmega[0],"0-5% (0-12%)","P");
   //leg2->AddEntry(graphOmega[1],"5-10%","P");
   //leg2->AddEntry(graphOmega[2],"10-20%","P");
   //leg2->AddEntry(graphOmega[3],"20-30%","P");
   //leg2->AddEntry(graphOmega[4],"30-40%","P");
   leg2->AddEntry(graphOmega[5],"40-50%","P");
   leg2->AddEntry(graphOmega[6],"50-60%","P");
   leg2->AddEntry(graphOmega[7],"60-70%","P");
   leg2->AddEntry(graphOmega[8],"70-80%","P");
   leg2->Draw();
   c1->SaveAs("Ratios.pdf");
}

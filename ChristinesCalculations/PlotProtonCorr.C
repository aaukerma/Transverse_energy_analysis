Int_t nenergies = 8;
double yieldratio[8][9];
double yieldratioerr[8][9];
double etratio[8][9];
double etratioerr[8][9];
//double cb[8][9];
//double cberr[8][9];
double cbbinedges[] = {0,5,10,20,30,40,50,60,70,80};
double energy[8] = {7.7,11.5,19.6,27,39,62.4,130,200};
int color[9] = {kViolet-5,kBlue,kTeal-5,kGreen+2,kOrange-3,kOrange+2,kRed-2,kRed,kPink+10};
int marker[9] = {20,21,22,23,29,33,34,47,43};
int marker2[9] = {24,25,26,32,30,27,28,46,42};
int ncent = 9;
void SetStyles(TGraph *g,int color,int marker){
  g->SetMarkerStyle(marker);
  g->SetMarkerColor(color);
  g->SetLineColor(color);
  g->SetMarkerSize(2.0);

}
void PlotProtonCorr(){
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  char *MYFILENAME = "protoncorrection.txt";
  ifstream file(MYFILENAME);
  cout<<"ncent "<<ncent<<" nenergies "<<nenergies<<endl;
  TGraphErrors *graph[8];
  TGraphErrors *graphyield[8];
  TGraphErrors *graphsNN[9];
  TGraphErrors *graphyieldsNN[9];
  for(int i=0;i<nenergies;i++){
    graph[i] = new TGraphErrors();
    SetStyles(graph[i],color[i],marker[i]);
    graphyield[i] = new TGraphErrors();
    SetStyles(graphyield[i],color[i],marker2[i]);
  }
  for(int i=0;i<ncent;i++){
    graphsNN[i] = new TGraphErrors();
    SetStyles(graphsNN[i],color[i],marker[i]);
    graphyieldsNN[i] = new TGraphErrors();
    SetStyles(graphyieldsNN[i],color[i],marker2[i]);
  }
  double dat0,dat1, dat2, dat3, dat4, dat5;
  float shift = 0.3;
  for(int i=0;i<nenergies;i++){
    for(int j=0;j<ncent;j++){
      //if(!file.eof()){
      file>>dat0>>dat1>>dat2>>dat3>>dat4>>dat5;
      //cout<<"i "<<i<<" j "<<j;
      //cout<<" "<<dat0<<" "<<dat1<<" "<<dat2<<" "<<dat3<<" "<<dat4<<" "<<dat5;
      //cout<<endl;
      double cb = (cbbinedges[j+1]+cbbinedges[j])/2.0;
      double cberr = (cbbinedges[j+1]-cbbinedges[j])/2.0;
      //energy[j] = dat1;
      //cb[i][j] = dat2;
      //cberr[i][j] = dat3;
      etratio[i][j] = dat2;
      etratioerr[i][j] = dat3;
      yieldratio[i][j] = dat4;
      yieldratioerr[i][j] = dat5;
      graph[i]->SetPoint(j,cb+0.3*i,dat2);
      graph[i]->SetPointError(j,cberr,dat3);
      graphyield[i]->SetPoint(j,cb+0.3*i,dat4);
      graphyield[i]->SetPointError(j,cberr,dat5);
      float myshift = exp(j*0.03);
      cout<<"shift "<<myshift<<endl;
      graphsNN[j]->SetPoint(i,energy[i]*myshift,dat2);
      graphsNN[j]->SetPointError(i,0,dat3);
      cout<<"energy "<<energy[i]<<" + "<<myshift<<endl;
      graphyieldsNN[j]->SetPoint(i,energy[i]*myshift,dat4);
      graphyieldsNN[j]->SetPointError(i,0,dat5);
      //if(dat4+dat5>high) high = dat4+dat5;
      //file>>energy[j]>>cb[i][j]>>cberr[i][j]>>kaonratio[i][j]>>kaonratioerr[i][j];
      //}
    }
  }
  TH1F *frame = new TH1F("frame","frame",1,0,80);
   frame->GetXaxis()->SetTitle("centrality");
   frame->GetXaxis()->SetLabelFont(42);
   frame->GetXaxis()->SetLabelSize(0.07);
   frame->GetXaxis()->SetTitleSize(0.07);
   frame->GetXaxis()->SetTitleOffset(1.0);
   frame->GetYaxis()->SetTitleOffset(1.0);
   frame->GetXaxis()->SetTitleFont(42);
   frame->GetYaxis()->SetTitle("f_{p},f_{p,yield}");
   frame->GetYaxis()->SetLabelFont(42);
   frame->GetYaxis()->SetLabelSize(0.07);
   frame->GetYaxis()->SetTitleSize(0.07);
   frame->GetYaxis()->SetTitleFont(42);
   frame->GetZaxis()->SetLabelFont(42);
   frame->GetZaxis()->SetLabelSize(0.07);
   frame->GetZaxis()->SetTitleSize(0.07);
   frame->GetZaxis()->SetTitleOffset(1);
   frame->GetZaxis()->SetTitleFont(42);
   frame->SetMaximum(3.5);
   frame->SetMinimum(1.5);
   //frame->Draw();

   TCanvas *c1 = new TCanvas("c1", "omega/pi0",600,400);
   c1->Range(-1.768938,-0.3125,15.92044,2.8125);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);
   c1->SetLeftMargin(0.1425);
   c1->SetRightMargin(0.0234114);
   c1->SetTopMargin(0.0321508);
   c1->SetBottomMargin(0.160754);
   graph[0]->SetHistogram(frame);
   graph[0]->Draw("AP");
   for(int i=0;i<nenergies;i++){graphyield[i]->Draw("P same");}
   for(int i=0;i<nenergies;i++){graph[i]->Draw("P same");}
   TLegend *leg2 = new TLegend(0.172819,0.789894,0.375839,0.960106);//0.753356,0.168449,0.956376,0.491979);//0.17,0.461197,0.3725,0.868071);//0.461918,1.116760.6125,0.657428,0.815833,0.956763);//0.611481,0.537099,0.896839,0.89701,NULL,"brNDC");
leg2->SetBorderSize(0);
leg2->SetTextFont(72);
leg2->SetTextSize(0.0481283);
leg2->SetLineColor(0);
leg2->SetLineStyle(0);
leg2->SetLineWidth(0);
leg2->SetFillColor(0);
leg2->SetFillStyle(0);
 leg2->AddEntry(graphsNN[0],"0-5%","P");
 leg2->AddEntry(graphsNN[1],"5-10%","P");
 leg2->AddEntry(graphsNN[2],"10-20%","P");
 TLegend *leg3 = new TLegend(0.32107,0.789894,0.525084,0.960106);//0.172819,0.176471,0.375839,0.395722);//0.753356,0.168449,0.956376,0.491979);//0.17,0.461197,0.3725,0.868071);//0.461918,1.116760.6125,0.657428,0.815833,0.956763);//0.611481,0.537099,0.896839,0.89701,NULL,"brNDC");
leg3->SetBorderSize(0);
leg3->SetTextFont(72);
leg3->SetTextSize(0.0481283);
leg3->SetLineColor(0);
leg3->SetLineStyle(0);
leg3->SetLineWidth(0);
leg3->SetFillColor(0);
leg3->SetFillStyle(0);
 leg3->AddEntry(graphsNN[3],"20-30%","P");
 leg3->AddEntry(graphsNN[4],"30-40%","P");
 leg3->AddEntry(graphsNN[5],"40-50%","P");
 TLegend *leg5 = new TLegend(0.481605,0.789894,0.685619,0.960106);//0.172819,0.176471,0.375839,0.395722);//0.753356,0.168449,0.956376,0.491979);//0.17,0.461197,0.3725,0.868071);//0.461918,1.116760.6125,0.657428,0.815833,0.956763);//0.611481,0.537099,0.896839,0.89701,NULL,"brNDC");
leg5->SetBorderSize(0);
leg5->SetTextFont(72);
leg5->SetTextSize(0.0481283);
leg5->SetLineColor(0);
leg5->SetLineStyle(0);
leg5->SetLineWidth(0);
leg5->SetFillColor(0);
leg5->SetFillStyle(0);
 leg5->AddEntry(graphsNN[6],"50-60%","P");
 leg5->AddEntry(graphsNN[7],"60-70%","P");
 leg5->AddEntry(graphsNN[8],"70-80%","P");

 float xlow = 5;
 float xhigh = 270;
 TLine *lineLow = new TLine(xlow,2.0,xhigh,2.0);
 lineLow->SetLineStyle(2);
 lineLow->SetLineWidth(3);
 lineLow->SetLineColor(1);
 TLine *lineHigh = new TLine(xlow,1.0+121.0/79.0,xhigh,1.0+121.0/79.0);
 lineHigh->SetLineStyle(2);
 lineHigh->SetLineWidth(3);
 lineHigh->SetLineColor(1);
  TH1F *frame2 = new TH1F("frame2","frame2",200,xlow,xhigh);
   frame2->GetXaxis()->SetTitle("#sqrt{s_{NN}}");
   frame2->GetXaxis()->SetLabelFont(42);
   frame2->GetXaxis()->SetLabelSize(0.07);
   frame2->GetXaxis()->SetTitleSize(0.07);
   frame2->GetXaxis()->SetTitleOffset(1.0);
   frame2->GetYaxis()->SetTitleOffset(1.0);
   frame2->GetXaxis()->SetTitleFont(42);
   frame2->GetYaxis()->SetTitle("f_{p},f_{p,yield}");
   frame2->GetYaxis()->SetLabelFont(42);
   frame2->GetYaxis()->SetLabelSize(0.07);
   frame2->GetYaxis()->SetTitleSize(0.07);
   frame2->GetYaxis()->SetTitleFont(42);
   frame2->GetZaxis()->SetLabelFont(42);
   frame2->GetZaxis()->SetLabelSize(0.07);
   frame2->GetZaxis()->SetTitleSize(0.07);
   frame2->GetZaxis()->SetTitleOffset(1);
   frame2->GetZaxis()->SetTitleFont(42);
   frame2->SetMaximum(3.5);
   frame2->SetMinimum(1.6);
   //frame2->Draw();

   TCanvas *c2 = new TCanvas("c2", "omega/pi0",600,400);
   c2->Range(-1.768938,-0.3125,15.92044,2.8125);
   c2->SetFillColor(0);
   c2->SetBorderMode(0);
   c2->SetBorderSize(2);
   c2->SetFrameBorderMode(0);
   c2->SetFrameBorderMode(0);
   c2->SetLeftMargin(0.1425);
   c2->SetRightMargin(0.0234114);
   c2->SetTopMargin(0.0321508);
   c2->SetBottomMargin(0.160754);
   c2->SetLogx();
   graphsNN[0]->SetHistogram(frame2);
   graphsNN[0]->Draw("AP");
   lineLow->Draw();
   lineHigh->Draw();
   for(int i=0;i<nenergies;i++){graphyieldsNN[i]->Draw("P same");}
   for(int i=0;i<nenergies;i++){graphsNN[i]->Draw("P same");}
   leg2->Draw();
   leg3->Draw();
   leg5->Draw();

   TLegend *leg4 = new TLegend(0.652174,0.805851,0.856187,0.949468);//0.692308,0.789894,0.896321,0.968085);//0.43311,0.710106,0.635452,0.960106);//0.172819,0.176471,0.375839,0.395722);//0.753356,0.168449,0.956376,0.491979);//0.17,0.461197,0.3725,0.868071);//0.461918,1.116760.6125,0.657428,0.815833,0.956763);//0.611481,0.537099,0.896839,0.89701,NULL,"brNDC");
leg4->SetBorderSize(0);
leg4->SetTextFont(72);
leg4->SetTextSize(0.0664894);
leg4->SetLineColor(0);
leg4->SetLineStyle(0);
leg4->SetLineWidth(0);
leg4->SetFillColor(0);
leg4->SetFillStyle(0);
 leg4->AddEntry(graph[0],"f_{p}","P");
 leg4->AddEntry(graphyield[0],"f_{p,yield}","P");
 leg4->Draw();
 c2->SaveAs("ProtonCorr.pdf");
}

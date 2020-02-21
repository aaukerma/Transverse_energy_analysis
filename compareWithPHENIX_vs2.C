// macro to overlay two graphs on the same canvas for comparison

{
TCanvas *c1 = new TCanvas();
//c1 -> DrawFrame(0.000000000000001, 0.000000000000000001, 2800., 12.);
TFile* file1 = TFile::Open("./crossCheckGraphs_TB.root");
TFile* file2 = TFile::Open("./plt4aAdare.root");
TH1D* g1 = (TH1D*)file1->Get("dETdEtaOverNpartBy2SumCent0");
TH1D* g2 = (TH1D*)file2->Get("plt4aAdare_PHENIX");
///////// power law fit
TF1 *func = new TF1("func","[0]*pow(x,[1])",0,2800);
func->SetParameter(1,0.1);
func->SetParLimits(1,0.1,1.);
func->SetLineStyle(9);
func->SetLineColor(1);
g1->Fit("func");
	  		cout<<"Chi^2 "<<func->GetChisquare()
	  			<<" NDF "<<func->GetNDF()<<" Chi^2/NDF "
	  			<<func->GetChisquare()/func->GetNDF()<<endl;


//g1->GetXaxis()->SetTicks("+");
c1 -> DrawFrame(5., 0.5, 50., 5.);
g1->SetMarkerColor(1); g1->SetMarkerSize(1);
g2 -> SetMarkerStyle(kFullCircle); g2->SetMarkerColor(4); g2->SetMarkerSize(1);
TString ystring="E_{T}/(N_{part}/2)";
TString xstring="#sqrt{s_{NN}}";
g1->Draw("same PE");
g1 -> SetMarkerStyle(kFullStar); /// needs to be called after Draw("A*")
g1->Draw("same PE");
g2->Draw("same PE");
func->GetYaxis()->SetTitle(ystring);
func->GetXaxis()->SetTitle(xstring);
g1 -> SetName("g1");
g2 -> SetName("g2");
func -> Draw("same");

TLegend* leg = new TLegend(0.5,0.15,0.7,0.3);

leg -> Draw();
leg -> AddEntry("g1", "This Analysis 0-5% Au+Au (Spectral)", "p");
leg -> AddEntry("g2", "PHENIX 0-5% Au+Au (Calorimetric)", "p");
leg -> SetFillStyle(0);
leg -> SetFillColor(0);
leg -> SetBorderSize(0);
leg -> SetTextSize(0.0334759);
c1->SetLogx(); c1->SetLogy();

c1 -> Update();
c1->SaveAs("./STAR_PHENIX_vs2.png");// STAR, PHENIX, CMS, ALICE
return c1;
}

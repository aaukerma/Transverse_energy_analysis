/// evolved from fitBESData_4_1.cpp
// added centrality column in output file : DONE
	// corresponding change in code to retreive centrality information from root file: DONE
// changed all functions (esp the fitting function funcBGBW) to have 6 parameters
	// otherwise the cov. matrix in later functions, with 5 parameters, is not consistent
// added // h->SetMaximum(5*(h->GetMaximum()));

#include <iostream>
#include <string>
#include "TKey.h"
#include <sstream>
#include <fstream>
#include "fitBESData10.h"
using namespace std;
// forward declarations for methods in fitBESData.h:


Double_t getdNdptOverptIntegrand(Double_t* rad, Double_t* par);// not used
Double_t getdNdpt(Double_t* pT, Double_t* params);
string concatenateHistoname(string,string,string,string);
Double_t* getIntegralsAndErrorsFromData(TH1D*, Double_t, Double_t);
/////Double_t* getIntegralsAndErrorsFromFit(Double_t* myPt, Double_t* par);
Double_t getdNdpt(Double_t* pT, Double_t* params);
Double_t getdETdEtaIntegrand(Double_t* myPt, Double_t* par);
Double_t getdETdyIntegrand(Double_t* myPt, Double_t* par);
Double_t getdNdEtaIntegrand(Double_t* myPt, Double_t* par);
Double_t getdNdyIntegrand(Double_t* myPt, Double_t* par);
Int_t* getNpartAndErr(Double_t collisionEnergy, string centrality);

//-------------------------HAGE
Double_t getdNdptHAGE(Double_t* pT, Double_t* params);
Double_t getdETdEtaIntegrandHAGE(Double_t* myPt, Double_t* par);
Double_t getdETdyIntegrandHAGE(Double_t* myPt, Double_t* par);
Double_t getdNdEtaIntegrandHAGE(Double_t* myPt, Double_t* par);
Double_t getdNdyIntegrandHAGE(Double_t* myPt, Double_t* par);

//--------------------expo
Double_t getdNdptEXPO(Double_t* pT, Double_t* params);
Double_t getdETdEtaIntegrandEXPO(Double_t* myPt, Double_t* par);
Double_t getdETdyIntegrandEXPO(Double_t* myPt, Double_t* par);
Double_t getdNdEtaIntegrandEXPO(Double_t* myPt, Double_t* par);
Double_t getdNdyIntegrandEXPO(Double_t* myPt, Double_t* par);


//---------------------expo2
Double_t getdNdptEXPO2(Double_t* pT, Double_t* params);
Double_t getdETdEtaIntegrandEXPO2(Double_t* myPt, Double_t* par);
Double_t getdETdyIntegrandEXPO2(Double_t* myPt, Double_t* par);
Double_t getdNdEtaIntegrandEXPO2(Double_t* myPt, Double_t* par);
Double_t getdNdyIntegrandEXPO2(Double_t* myPt, Double_t* par);

// main function:
int fit7pp(){
std::ofstream avg ("avgfitResults3.dat", std::ofstream::app);
TFile* myFile = TFile::Open("HEPData-ins1357424-v1-Table_1.root");

TDirectory* dir = myFile->GetDirectory("Table 1");
TH1D* h=(TH1D* )myFile->Get("Table 1/Hist1D_y3");
Double_t mass = 0.93827; 
Double_t type = 0.0;
Double_t collEn= 7000.0;
Double_t* integralDataPtr;
		// TODO : need to fix what function this should be:
		integralDataPtr = getIntegralsAndErrorsFromData(h,type,mass);
	
TF1* EXPO;
TF1* EXPO2;
TCanvas* c1;
TGraphAsymmErrors *g = (TGraphAsymmErrors *)dir->Get("Graph1D_y3");
g->Draw();
Double_t x,y,xerrlow,xerrhigh,yerrlow,yerrhigh;

for(int i=0;i<g->GetN();i++){
	g->GetPoint(i,x,y);
	xerrlow = g->	GetErrorXlow(i);
	xerrhigh = g->	GetErrorXhigh(i);
	yerrlow = g->	GetErrorYlow(i);
	yerrhigh = g->	GetErrorYhigh(i);
	cout<<"i "<<i<<" x "<<x<<" +/- "<<xerrlow<<" y "<<y<< " +/- "<<yerrlow<<endl;
	//cout<< endl<< " 1/N " << (x*integralDataPtr[7])*(1.0/510.0)<< endl;
	cout<<endl<<"status of histo bin:"<<endl;
	
	cout<<" y "<<h->GetBinContent(i+1)<<endl;
	cout<<" x "<<h->GetXaxis()->GetBinLowEdge(i+1)<<" - "<<h->GetXaxis()->GetBinCenter(i+1)<<" - "<<h->GetXaxis()->GetBinLowEdge(i+2)<<endl;
	h->SetBinError(i+1,yerrlow);
}
h->SetLineColor(kRed);
h->Draw("same");


	
	TF1* dETdEtaIntegrandFuncEXPO;
	TF1* dETdyIntegrandFuncEXPO;
	TF1* dNdEtaIntegrandFuncEXPO;
	TF1* dNdyIntegrandFuncEXPO;

	TF1* dETdEtaIntegrandFuncEXPO2;
	TF1* dETdyIntegrandFuncEXPO2;
	TF1* dNdEtaIntegrandFuncEXPO2;
	TF1* dNdyIntegrandFuncEXPO2;
	string histoName = h->GetName();
	c1 = new TCanvas(); 



	///////
	
	EXPO = new TF1("getdNdptEXPO", getdNdptEXPO,0.000000000001,10.,4);

		EXPO2 = new TF1("getdNdptEXPO2", getdNdptEXPO2,0.000000000001,10.,5);
		
	
		///////
		dETdEtaIntegrandFuncEXPO = new TF1("dETdEtaIntegrandEXPO",
									getdETdEtaIntegrandEXPO,
									0, 10, 4 );// function goes from 0 to 10
										// and has 6 parameters"
										// mass, pt,
		dETdyIntegrandFuncEXPO = new TF1("dETdyIntegrandEXPO",
								  getdETdyIntegrandEXPO,
								  0,10,4);
		dNdEtaIntegrandFuncEXPO = new TF1("dETdyIntegrandEXPO",
								  getdNdEtaIntegrandEXPO,
								  0,10,4); // 5 parameters:m, 6th is type*0
		dNdyIntegrandFuncEXPO = new TF1("dETdyIntegrandEXPO",
								  getdNdEtaIntegrandEXPO,
								  0,10,4);// 5 parameters:m, 6th is type*0

	
		
		//Double_t chi2BGBWE = EXPO->GetChisquare();
		//Double_t nDFBGBW2E = EXPO->GetNDF();
		//Double_t p2E = EXPO->GetParameter(2);
		//Double_t e2E = EXPO->GetParError(2);
		///////
dETdEtaIntegrandFuncEXPO2 = new TF1("dETdEtaIntegrandEXPO2",
									getdETdEtaIntegrandEXPO2,
									0, 10, 5 );// function goes from 0 to 10
										// and has 6 parameters"
										// mass, pt,
		dETdyIntegrandFuncEXPO2 = new TF1("dETdyIntegrandEXPO2",
								  getdETdyIntegrandEXPO2,
								  0,10,5);
		dNdEtaIntegrandFuncEXPO2 = new TF1("dETdyIntegrandEXPO2",
								  getdNdEtaIntegrandEXPO2,
								  0,10,5); // 5 parameters:m, 6th is type*0
		dNdyIntegrandFuncEXPO2 = new TF1("dETdyIntegrandEXPO2",
								  getdNdEtaIntegrandEXPO2,  0,10,5);
								  
				/////////////////
				
//////////////////////
				//gPad->SetLogy();
		gStyle->SetOptFit(1111);// display fit parameters; customizable
		gStyle->SetOptDate();
		gROOT-> SetBatch(kTRUE);// save canvases without displaying them
		c1->Update();				// 5 parameters:m, 6th is type*0
		//Double_t chi2BGBWE2 = EXPO2->GetChisquare();
		//Double_t nDFBGBW2E2 = EXPO2->GetNDF();
		//Double_t p2E2 = EXPO2->GetParameter(2);
		//Double_t e2E2 = EXPO2->GetParError(2);
		///////////////////
		

//string centrality;
//string particleID;
// int proBinNumber_eta_inel0 = h->GetNbinsX();
 /// cout<<"Get the number of bins: ";

  //cout<<proBinNumber_eta_inel0;
Double_t avgET=0.0;
		Double_t avgET_err=0.0;
		Double_t avgN=0.0;
		Double_t avgN_err=0.0;
		Double_t avgNpart=0.0;
		Double_t avgNpart_err=0.0;
		
		
		
		EXPO->SetParameters(10000.,0.01,mass,type);
		EXPO->SetParNames("Ae","Be","mass","type");
		//EXPO->SetParLimits(0,1,2);
		//EXPO->SetParLimits(1,-0.1,0.01);
		EXPO->SetLineColor(kBlack);
		EXPO->FixParameter(3,type);
		EXPO->FixParameter(2,mass);
		
		EXPO2->SetParameters(1.,1.,1.,mass,type);
		EXPO2->SetParNames("A2","B2","C2","mass","type");
			EXPO2->SetParLimits(0,4.0,5.515);//A
		EXPO2->SetParLimits(1,-0.9,-0.8757);//B
		EXPO2->SetParLimits(2,0.15,0.2633);//C
		EXPO2->SetLineColor(kViolet);
		EXPO2->FixParameter(4,type);
		EXPO2->FixParameter(3,mass);
		ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(20000);
		
		
		TFitResultPtr v = h->Fit("getdNdptEXPO","S+","M",0.000000000001,10.);
		Double_t meanpt4= EXPO->GetHistogram()->GetMean();
		
		TFitResultPtr w = h->Fit("getdNdptEXPO2","S+","M",0.000000000001,10.);
		Double_t meanpt5= EXPO2->GetHistogram()->GetMean();
		
		Double_t chi4prob4 = v->Prob();
		Double_t chi5prob5 = w->Prob();
		h->SetMaximum(5*(h->GetMaximum()));
		//h-> GetYaxis()->SetRangeUser(0.,maxY);
		//TMatrixDSym cov = r->GetCovarianceMatrix();
		h-> GetXaxis()->SetRangeUser(0.,10.);
		TString xlabel = "p_{T}";
		TString ylabel = "#frac{d^{2}N}{dydp_{T}}";
		h-> SetXTitle(xlabel);
		h-> SetYTitle(ylabel);
		
		
		
		Double_t Ae = EXPO->GetParameter(0);
		Double_t Be = EXPO->GetParameter(1);
		Double_t AeErr = EXPO->GetParError(0);
		Double_t BeErr = EXPO->GetParError(1);
		
		Double_t A2 = EXPO2->GetParameter(0);
		Double_t B2 = EXPO2->GetParameter(1);
		Double_t C2 = EXPO2->GetParameter(2);
		Double_t A2Err = EXPO2->GetParError(0);
		Double_t B2Err = EXPO2->GetParError(1);
		Double_t C2Err = EXPO2->GetParError(2);
			
		/////////////
		EXPO			 	-> SetParameters(Ae,Be,mass,type);
		dETdEtaIntegrandFuncEXPO 	-> SetParameters(Ae,Be,mass,type);
		dETdEtaIntegrandFuncEXPO	-> FixParameter(3,type);
		dETdEtaIntegrandFuncEXPO	-> FixParameter(2,mass);
		dETdyIntegrandFuncEXPO 		-> SetParameters(Ae,Be,mass,type);
		dETdyIntegrandFuncEXPO		-> FixParameter(3,type);
		dETdyIntegrandFuncEXPO		-> FixParameter(2,mass);
		dNdEtaIntegrandFuncEXPO 	-> SetParameters(Ae,Be,mass,type);
		dNdEtaIntegrandFuncEXPO		-> FixParameter(2,mass);
		dNdEtaIntegrandFuncEXPO		-> FixParameter(3,type);
		dNdyIntegrandFuncEXPO		-> SetParameters(Ae,Be,mass,type);
		dNdyIntegrandFuncEXPO		-> FixParameter(2,mass);
		dNdyIntegrandFuncEXPO		-> FixParameter(3,type);

		Int_t totBinsE 	= h->GetNbinsX();
		Int_t binx1E 	= 0;
		Int_t binx2E 	= totBinsE+1;

		Double_t leftCutE 	= h->GetXaxis()->GetBinLowEdge(binx1E+2);
		Double_t rightCutE 	= h->GetXaxis()->GetBinUpEdge(binx2E-1);

		Double_t dETdEtaLeftE 	= dETdEtaIntegrandFuncEXPO -> Integral(0.,leftCutE);
		Double_t dETdEtaRightE 	= dETdEtaIntegrandFuncEXPO -> Integral(rightCutE,10.);
		Double_t dETdyLeftE 		= dETdyIntegrandFuncEXPO -> Integral(0.,leftCutE);
		Double_t dETdyRightE 	= dETdyIntegrandFuncEXPO -> Integral(rightCutE,10.);
		Double_t dNdEtaLeftE 	= dNdEtaIntegrandFuncEXPO -> Integral(0.,leftCutE);
		Double_t dNdEtaRightE 	= dNdEtaIntegrandFuncEXPO -> Integral(rightCutE,10.);
		Double_t dNdyLeftE 		= dNdyIntegrandFuncEXPO -> Integral(0.,leftCutE);
		Double_t dNdyRightE 		= dNdyIntegrandFuncEXPO -> Integral(rightCutE,10.);
		// Errors:
		Double_t dETdEtaLErrE	=
		dETdEtaIntegrandFuncEXPO->IntegralError(0.,leftCutE,
								v->GetParams(),
								v->GetCovarianceMatrix().GetMatrixArray());
		Double_t dETdEtaRErrE	=
		dETdEtaIntegrandFuncEXPO->IntegralError(rightCutE,10.,
								v->GetParams(),
								v->GetCovarianceMatrix().GetMatrixArray());
		Double_t dETdyLErrE	=
		dETdyIntegrandFuncEXPO->IntegralError(0.,leftCutE,
								v->GetParams(),
								v->GetCovarianceMatrix().GetMatrixArray());
		Double_t dETdyRErrE	=
		dETdyIntegrandFuncEXPO->IntegralError(rightCutE,10.,
								v->GetParams(),
								v->GetCovarianceMatrix().GetMatrixArray());
		Double_t dNdEtaLErrE	=
		dNdEtaIntegrandFuncEXPO->IntegralError(0.,leftCutE,
								v->GetParams(),
								v->GetCovarianceMatrix().GetMatrixArray());
		Double_t dNdEtaRErrE	=
		dETdEtaIntegrandFuncEXPO->IntegralError(rightCutE,10.,
								v->GetParams(),
								v->GetCovarianceMatrix().GetMatrixArray());
		Double_t dNdyLErrE	=
		dNdyIntegrandFuncEXPO->IntegralError(0.,leftCutE,
								v->GetParams(),
								v->GetCovarianceMatrix().GetMatrixArray());
		Double_t dNdyRErrE	=
		dNdyIntegrandFuncEXPO->IntegralError(rightCutE,10.,
								v->GetParams(),
								v->GetCovarianceMatrix().GetMatrixArray());

		Double_t dETdEta_dE= *(integralDataPtr+0);
		Double_t dETdEta_d_errE = *(integralDataPtr+1);
		Double_t dETdEtaTotalE = dETdEtaLeftE+dETdEta_dE+dETdEtaRightE;
		Double_t dETdEtaTErrE = dETdEtaLErrE+dETdEta_d_errE+dETdEtaRErrE;

		Double_t dETdy_dE = *(integralDataPtr+2);
		Double_t dETdy_d_errE = *(integralDataPtr+3);
		Double_t dETdyTotalE = dETdyLeftE+dETdy_dE+dETdyRightE;
		Double_t dETdyTErrE = dETdyLErrE+dETdy_d_errE+dETdyRErrE;

		Double_t dNdEta_dE = *(integralDataPtr+4);
		Double_t dNdEta_d_errE = *(integralDataPtr+5);
		Double_t dNdEtaTotalE = dNdEtaLeftE+dNdEta_dE+dNdEtaRightE;
		Double_t dNdEtaTErrE = dNdEtaLErrE+dNdEta_d_errE+dNdEtaRErrE;

		Double_t dNdy_dE = *(integralDataPtr+6);
		Double_t dNdy_d_errE = *(integralDataPtr+7);
		Double_t dNdyTotalE = dNdyLeftE+dNdy_dE+dNdyRightE;
		Double_t dNdyTErrE = dNdyLErrE+dNdy_d_errE+dNdyRErrE;
		c1->Update();
	
		cout<< " dNdYE " << dNdyTotalE <<endl;
		cout<<" dNdEtaTotalE " << dNdEtaTotalE<<endl;
		Double_t chi2BGBWE = EXPO->GetChisquare();
		Double_t nDFBGBW2E = EXPO->GetNDF();
		Double_t p2E = EXPO->GetParameter(2);
		Double_t e2E = EXPO->GetParError(2);
////////////
	EXPO2			 	-> SetParameters(A2,B2,C2,mass,type);
		dETdEtaIntegrandFuncEXPO2 	-> SetParameters(A2,B2,C2,mass,type);
		dETdEtaIntegrandFuncEXPO2	-> FixParameter(4,type);
		dETdEtaIntegrandFuncEXPO2	-> FixParameter(3,mass);
		dETdyIntegrandFuncEXPO2 		-> SetParameters(A2,B2,C2,mass,type);
		dETdyIntegrandFuncEXPO2		-> FixParameter(4,type);
		dETdyIntegrandFuncEXPO2		-> FixParameter(3,mass);
		dNdEtaIntegrandFuncEXPO2 	-> SetParameters(A2,B2,C2,mass,type);
		dNdEtaIntegrandFuncEXPO2	-> FixParameter(3,mass);
		dNdEtaIntegrandFuncEXPO2		-> FixParameter(4,type);
		dNdyIntegrandFuncEXPO2		-> SetParameters(A2,B2,C2,mass,type);
		dNdyIntegrandFuncEXPO2		-> FixParameter(3,mass);
		dNdyIntegrandFuncEXPO2		-> FixParameter(4,type);

		Int_t totBinsE2 	= h->GetNbinsX();
		Int_t binx1E2 	= 0;
		Int_t binx2E2 	= totBinsE2+1;

		Double_t leftCutE2 	= h->GetXaxis()->GetBinLowEdge(binx1E2+2);
		Double_t rightCutE2 	= h->GetXaxis()->GetBinUpEdge(binx2E2-1);

		Double_t dETdEtaLeftE2 	= dETdEtaIntegrandFuncEXPO2 -> Integral(0.,leftCutE2);
		Double_t dETdEtaRightE2 	= dETdEtaIntegrandFuncEXPO2 -> Integral(rightCutE2,10.);
		Double_t dETdyLeftE2 		= dETdyIntegrandFuncEXPO2 -> Integral(0.,leftCutE2);
		Double_t dETdyRightE2 	= dETdyIntegrandFuncEXPO2 -> Integral(rightCutE2,10.);
		Double_t dNdEtaLeftE2	= dNdEtaIntegrandFuncEXPO2 -> Integral(0.,leftCutE2);
		Double_t dNdEtaRightE2 	= dNdEtaIntegrandFuncEXPO2 -> Integral(rightCutE2,10.);
		Double_t dNdyLeftE2 		= dNdyIntegrandFuncEXPO2 -> Integral(0.,leftCutE2);
		Double_t dNdyRightE2 		= dNdyIntegrandFuncEXPO2 -> Integral(rightCutE2,10.);
		// Errors:
		Double_t dETdEtaLErrE2	=
		dETdEtaIntegrandFuncEXPO2->IntegralError(0.,leftCutE2,
								w->GetParams(),
								w->GetCovarianceMatrix().GetMatrixArray());
		Double_t dETdEtaRErrE2	=
		dETdEtaIntegrandFuncEXPO2->IntegralError(rightCutE2,10.,
								w->GetParams(),
								w->GetCovarianceMatrix().GetMatrixArray());
		Double_t dETdyLErrE2	=
		dETdyIntegrandFuncEXPO2->IntegralError(0.,leftCutE2,
								w->GetParams(),
								w->GetCovarianceMatrix().GetMatrixArray());
		Double_t dETdyRErrE2	=
		dETdyIntegrandFuncEXPO2->IntegralError(rightCutE2,10.,
								w->GetParams(),
								w->GetCovarianceMatrix().GetMatrixArray());
		Double_t dNdEtaLErrE2	=
		dNdEtaIntegrandFuncEXPO2->IntegralError(0.,leftCutE2,
								w->GetParams(),
								w->GetCovarianceMatrix().GetMatrixArray());
		Double_t dNdEtaRErrE2	=
		dETdEtaIntegrandFuncEXPO2->IntegralError(rightCutE2,10.,
								w->GetParams(),
								w->GetCovarianceMatrix().GetMatrixArray());
		Double_t dNdyLErrE2	=
		dNdyIntegrandFuncEXPO2->IntegralError(0.,leftCutE2,w->GetParams(),w->GetCovarianceMatrix().GetMatrixArray());
		Double_t dNdyRErrE2	=
		dNdyIntegrandFuncEXPO2->IntegralError(rightCutE2,10.,w->GetParams(),w->GetCovarianceMatrix().GetMatrixArray());

		Double_t dETdEta_dE2= *(integralDataPtr+0);
		Double_t dETdEta_d_errE2 = *(integralDataPtr+1);
		
	Double_t dETdEtaTotalE2 = dETdEtaLeftE2+dETdEta_dE2+dETdEtaRightE2;
		
Double_t dETdEtaTErrE2 = dETdEtaLErrE2+dETdEta_d_errE2+dETdEtaRErrE2;

		Double_t dETdy_dE2 = *(integralDataPtr+2);
		Double_t dETdy_d_errE2 = *(integralDataPtr+3);
		Double_t dETdyTotalE2 = dETdyLeftE2+dETdy_dE2+dETdyRightE2;
		
	Double_t dETdyTErrE2 = dETdyLErrE2+dETdy_d_errE2+dETdyRErrE2;

		Double_t dNdEta_dE2 = *(integralDataPtr+4);
		Double_t dNdEta_d_errE2 = *(integralDataPtr+5);
		
	Double_t dNdEtaTotalE2 = dNdEtaLeftE2+dNdEta_dE2+dNdEtaRightE2;
		
	Double_t dNdEtaTErrE2 = dNdEtaLErrE2+dNdEta_d_errE2+dNdEtaRErrE2;

		Double_t dNdy_dE2 = *(integralDataPtr+6);
		Double_t dNdy_d_errE2 = *(integralDataPtr+7);
		Double_t dNdyTotalE2 = dNdyLeftE2+dNdy_dE2+dNdyRightE2;
		Double_t dNdyTErrE2 = dNdyLErrE2+dNdy_d_errE2+dNdyRErrE2;
		//cout<< " dNdYE2 " << dNdyTotalE2 <<" ____ " << dNdYTErrE2<< endl;
	//cout<< " dNdYE " << dNdyTotalE << " ____ "<< dNdyTErrE<<endl;
		/*Int_t* NpartAndArrPtrE2;
		Int_t NpartE2;
		Int_t NpartErrE2;
		NpartAndArrPtrE2 = getNpartAndErr(collEn,centrality);
		NpartE2 = *(NpartAndArrPtrE2+0);
		NpartErrE2 = *(NpartAndArrPtrE2+1);*/
		c1->Update();
		Double_t chi2BGBWE2 = EXPO2->GetChisquare();
		Double_t nDFBGBW2E2 = EXPO2->GetNDF();
		Double_t p2E2 = EXPO2->GetParameter(2);
		Double_t e2E2 = EXPO2->GetParError(2);
		
		
		TString hisst = h->GetName();
		c1->Update();
	c1->cd();
	Double_t yloc = EXPO->GetHistogram()->GetMaximum();
	yloc = yloc/10.0;
	
	TF1* hj = new TF1();
	hj =dETdyIntegrandFuncEXPO ;
	hj->DrawCopy();
	hj->Draw();
	hj->GetHistogram()->GetXaxis()->SetRangeUser(0.,3.);
	hj->GetHistogram()->GetYaxis()->SetRangeUser(0.,(11.0*yloc));
	hj->GetHistogram()->SetXTitle(xlabel);
	hj->GetHistogram()->SetYTitle(ylabel);
	hj->SetTitle(hisst);
	hj->SetLineColor(kWhite);
	h->Draw("same");
	h->GetXaxis()->SetRangeUser(0.,10.);
	c1->Update();
	avgET=(dETdEtaTotalE+dETdEtaTotalE2)/2.0;
		avgET_err=(dETdEtaTErrE+dETdEtaTErrE2)/2.0;//averages the error in the points/fit.  Probably OK, especially since it looks like the fits actually pretty much agree and have small uncertainties
		Double_t avgETfitErr=TMath::Abs(dETdEtaTotalE-dETdEtaTotalE2)/2.0;
		Double_t avgETy=(dETdyTotalE+dETdyTotalE2)/2.0;
		Double_t avgETy_err=(dETdyTErrE+dETdyTErrE2)/2.0;//averages the error in the points/fit.  Probably OK, especially since it looks like the fits actually pretty much agree and have small uncertainties
		Double_t avgETyfitErr=TMath::Abs(dETdyTotalE-dETdyTotalE2)/2.0;
		avg <<collEn << "\t"
		    << mass << "\t"
		    << avgET <<  "\t"
		    << avgET_err <<  "\t"
		    << avgETfitErr <<  "\t"
		    << avgETy <<  "\t"
		    << avgETy_err <<  "\t"
		    << avgETyfitErr <<  "\n";
	string imgPathAndName = "./Table1/7Tevpp.png";
				//c1 -> SaveAs("./fittedPlots/trial1.png");
		TImage *png = TImage::Create();// FIXME try to use canvas method instead of png object
		png->FromPad(c1);
		const char* imgPathAndNameConstCharPtr = imgPathAndName.c_str();
		png->WriteImage(imgPathAndNameConstCharPtr);
	gSystem->ProcessEvents();
	
	
delete h;
delete g;

delete EXPO;
delete dETdEtaIntegrandFuncEXPO;
delete dETdyIntegrandFuncEXPO;
delete EXPO2;
delete dETdEtaIntegrandFuncEXPO2;
delete dETdyIntegrandFuncEXPO2;
delete c1;
gObjectTable->Print();
delete myFile;
avg.close();
 
return 0;
}

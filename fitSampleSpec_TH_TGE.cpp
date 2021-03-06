/*
Evolved from fitBESData5.cpp

Used to analyze any TH1D object within a TFile object by customizing:
1. the name of the input TFile object myFile
2. name of the output file (stream) datFile
3. method to read TH1D object from Tfile depending on how it was created:
	3.1. h = (TH1D*)mikey->ReadObj();
			(only) in this case, string histoName = h->GetName(); makes sense
	3.2. h = (TH1D*)myFile->Get(Form("cent%i_proton_plus",0));
			this is the case applicable to SPECTRA_COMB_20120709.root
4. (depending on the histoName) particleID and its dependent variables
5. fit parameters: funcBGBW->SetParameters(mass,0.99,0.30,0.1,1000.,type);
6. ylabel: ylabel = "#frac{d^{2}N}{dydp_{T}}";
7. depending on the ylabel, possible transformation to be applied
	under the associated header file (eg. fitALICE2013Data.h) within functions:
	7.1. getIntegralsAndErrorsFromData(TH1D* hist, Double_t type, Double_t mass)
	7.2. getdNdpt(Double_t* pT, Double_t* params)
*/

#include <iostream>
#include <string>
#include "TKey.h"
#include <fstream>
#include "fitBESData5.h"
using namespace std;

// forward declarations for methods in fitALICE2013Data.h:
Double_t getdNdptOverptIntegrand(Double_t* rad, Double_t* par);// not used
Double_t getdNdpt(Double_t* pT, Double_t* params);
//DONE # change name to concatenatePlotName to include TGraphErrors along with TH1:
string concatenatePlotname(string,string,string,string);
//# create overloaded function with argument TGraphErrors* instead of TH1D*:
Double_t* getIntegralsAndErrorsFromData(TH1D*, Double_t, Double_t);
//# check what exactly the first argument passes, i.e., bin edges or centers:
Double_t getdNdpt(Double_t* pT, Double_t* params);
Double_t getdETdEtaIntegrand(Double_t* myPt, Double_t* par);
Double_t getdETdyIntegrand(Double_t* myPt, Double_t* par);
Double_t getdNdEtaIntegrand(Double_t* myPt, Double_t* par);
Double_t getdNdyIntegrand(Double_t* myPt, Double_t* par);
Int_t* getNpartAndErr(Double_t collisionEnergy, string centrality);

// main function:
int fitSampleSpec_TH_TGE(){
	std::ofstream datFile ("fitResults_la.dat", std::ofstream::out);
	datFile << "CollEn"<< "\t"	
			<< "particle" << "\t"
			<< "centrality" << "\t"
			<< "mass" << "\t"
			<< "beta" <<"\t"
			<< "betaErr" <<"\t"
			<< "temp" <<"\t"
			<< "tempErr" <<"\t"
			<< "n-v_prof" <<"\t"
			<< "nErr" <<"\t"
			<< "norm" <<"\t"
			<< "normErr" <<"\t"
			<< "dETdEta_data" << "\t" //dETdEta_d
			<< "dETdEta_d_err"<< "\t" //dETdEta_d_err
			<< "dETdEtaLeft" << "\t"
			<< "dETdEtaLErr" << "\t"
			<< "dETdEtaRight" << "\t"
			<< "dETdEtaRErr" << "\t"
			<< "dETdEtaTotal" << "\t" // dETdEtaTotal
			<< "dETdEtaTErr" << "\t"
			<< "dETdy_d" << "\t" //dETdy_d
			<< "dETdy_d_err"<< "\t" //dETdy_d_err
			<< "dETdyLeft" << "\t"
			<< "dETdyLErr" << "\t"
			<< "dETdyRight" << "\t"
			<< "dETdyRErr" << "\t"
			<< "dETdyTotal" << "\t" // dETdyTotal
			<< "dETdyTErr" << "\t"
			<< "dNdEta_d" << "\t" //dNdEta_d
			<< "dNdEta_d_err"<< "\t" //dNdEta_d_err
			<< "dNdEtaLeft" << "\t"
			<< "dNdEtaLErr" << "\t"
			<< "dNdEtaRight" << "\t"
			<< "dNdEtaRErr" << "\t"
			<< "dNdEtaTotal" << "\t" // dNdEtaTotal
			<< "dNdEtaTErr" << "\t"
			<< "dNdy_d" << "\t" //dNdy_d
			<< "dNdy_d_err"<< "\t" //dNdy_d_err
			<< "dNdyLeft" << "\t"
			<< "dNdyLErr" << "\t"
			<< "dNdyRight" << "\t"
			<< "dNdyRErr" << "\t"
			<< "dNdyTotal" << "\t" // dNdyTotal
			<< "dNdyTErr" << "\t"
			<< "Npart" << "\t"
			<< "NpartErr" << "\n";
	TFile* myFile = new TFile("BESLambdas.root");
	TIter next(myFile->GetListOfKeys());
	TKey* mikey;
	TH1D* h;
	TGraphErrors* tg;
	TCanvas* c1;
	TClass* class1;
	TF1* funcBGBW;
	TF1* dETdEtaIntegrandFunc;
	TF1* dETdyIntegrandFunc;
	TF1* dNdEtaIntegrandFunc;
	TF1* dNdyIntegrandFunc;
	int breakOutForTesting =0;
	int stop =270; // breakOut after this many iterations (if achieved); default: 140
	///while(1){
	while((mikey=(TKey*)next())){
		class1 = gROOT->GetClass(mikey->GetClassName());
		/////.....
		/*
		if(!class1->InheritsFrom("TH1") || !class1->InheritsFrom("TGraphErrors")){
			cout << "Object not TH1 or TGraphErrors";
			delete class1;
			mikey->DeleteBuffer();
			continue;
		}
		*/	
		c1 = new TCanvas(); // a la Rademakers
		funcBGBW = new TF1("getdNdpt",getdNdpt,0.00000000000001,30.,6); // actually has 5 parameters
												// 6th parameter, type, multiplied by 0 and added
												// for consistency of cov. matrix needed later		
		dETdEtaIntegrandFunc = new TF1("dETdEtaIntegrand", 
									getdETdEtaIntegrand, 
									0, 30, 6 );// function goes from 0 to 10, 30 etc.
										// and has 6 parameters
										// mass, beta, temp, n, norm, type
		dETdyIntegrandFunc = new TF1("dETdyIntegrand",
								  getdETdyIntegrand,
								  0,30,6);
		dNdEtaIntegrandFunc = new TF1("dETdyIntegrand",
								  getdNdEtaIntegrand,
								  0,30,6); // 5 parameters:m,b,t,n,norm, 6th is type*0
		dNdyIntegrandFunc = new TF1("dETdyIntegrand",
								  getdNdEtaIntegrand,
								  0,30,6);// 5 parameters:m,b,t,n,norm, 6th is type*0
		gPad->SetLogy();
		gStyle->SetOptFit(1111);// display fit parameters; customizable
		gStyle->SetOptDate();// display date (at bottom left)
		gROOT-> SetBatch(kTRUE);// save canvases without displaying them
		c1->Update();

		// read histogram object for current iteration of key:
		// was trying to make the method more generic, but due to time constraints,
		// just have the functionality for TGraphErrors objects for the time being
		/*
		if(class1->InheritsFrom("TH1")) h = (TH1D*)mikey->ReadObj();
		else tg = (TGraphErrors*)mikey->ReadObj();
		*/
		tg = (TGraphErrors*)mikey->ReadObj();
		//////////h = (TH1D*)myFile->Get(Form("cent%i_proton_plus",0));
		////////....string histoName = h->GetName();
		string graphName = tg->GetName();
		////////........
		/*
		if(histoName != "cent8_pi+_Au+Au_39"){ // pi-/+, ka-/+, pbar, proton
			breakOutForTesting++;
			///cout << breakOutForTesting << endl;
			if(breakOutForTesting>=stop) break;
			else continue;
		}
		*/
		if(graphName != "cent0_la_Au+Au_7.7"){
			breakOutForTesting++;
			///cout << breakOutForTesting << endl;
			if(breakOutForTesting>=stop) break;
			else continue;
		}
		/////....cout << "test bin content: " << h->GetBinContent(1) << endl;
		Double_t collEn = 0.;// initialize
		//cent8_ka+_Au+Au_7.7 // sample histo name
		if(graphName.substr( graphName.length() - 4 ) == "_7.7") collEn = 7.7;
		else if(graphName.substr( graphName.length() - 4 ) == "11.5") collEn = 11.5;
		else if(graphName.substr( graphName.length() - 4 ) == "19.6") collEn = 19.6;
		else if(graphName.substr( graphName.length() - 4 ) == "u_27") collEn = 27;
		else if(graphName.substr( graphName.length() - 4 ) == "u_39") collEn = 39;
		else if(graphName.substr( graphName.length() - 4 ) == "2.76") collEn = 2760;
		else collEn = 77777.;
		//get first three characters of particle name from histoName:
		string particleID = graphName.substr(6,3);// starting position in array:6, 3 chars total
		string centrality = graphName.substr(4,1);// starting position in array:4, 1 char total
		
		//------------ Assign mass & type to particle -----------------//
		Double_t mass = 0.93827; // in GeV
		
		// type Double_t instead of Int_t 
		 //to use as argument in TF1 method SetParameters()
		Double_t type = -1.;// 0 for mesons, -1 for baryons, 1 for antibaryons
		if		(particleID=="pi-"||particleID=="pi+"||particleID == "pio")
				{mass = 0.13957; type = 0.;}
		else if	(particleID=="ka-"||particleID=="ka+")
				{mass = 0.49368; type = 0.;}
		else if	(particleID=="pro")
				{mass = 0.93827; type = -1.;}
		else if	(particleID=="pba")
				{mass = 0.93827; type = 1.;}
		else if (particleID=="la_")
				{mass = 1.11568; type = -1.;}
		else if (particleID=="ala")
				{mass = 1.11568; type = 1.;}
		else if		(particleID=="pis")
				{mass = 0.13957; type = 0.;}
		else {cout << "Check particle: "
				<< particleID<<endl;/*return 1;*/}
		
		/////.....Double_t* integralDataPtr;

		/////////...........integralDataPtr = getIntegralsAndErrorsFromData(h,type,mass);
							// ^ method verified!!!
		
		
		//------------- Begin BGBW fit --------------------------//

		//////...........
		
		if (	graphName == "cent7_ka-_Au+Au_7.7"
			|| 	graphName == "cent7_ka-_Au+Au_11.5"
			||	graphName == "cent7_pi+_Au+Au_7.7"
			||	graphName == "cent4_pi-_Au+Au_19.6"
			||	graphName == "cent5_ka+_Au+Au_27"
			||	graphName == "cent5_ka-_Au+Au_7.7"
			||	graphName == "cent6_pi+_Au+Au_11.5"
			
			||	graphName == "cent3_pi-_Au+Au_7.7"
			||	graphName == "cent4_pi-_Au+Au_7.7"
			||	graphName == "cent5_pi-_Au+Au_7.7"
			||	graphName == "cent7_pi-_Au+Au_7.7"
			||	graphName == "cent1_pbar_Au+Au_7.7"	// NOT POS-DEF		
			||	graphName == "cent3_pbar_Au+Au_7.7"
			||	graphName == "cent4_pbar_Au+Au_7.7"
			||	graphName == "cent6_pbar_Au+Au_7.7"
			||	graphName == "cent7_pbar_Au+Au_7.7"
			||	graphName == "cent8_ka+_Au+Au_7.7"
			||	graphName == "cent8_ka-_Au+Au_7.7"
			||	graphName == "cent6_pi-_Au+Au_11.5"
			||	graphName == "cent0_pi+_Au+Au_11.5"
			||	graphName == "cent7_pi+_Au+Au_11.5"
			||	graphName == "cent8_proton_Au+Au_19.6"
			||	graphName == "cent5_proton_Au+Au_7.7"
			||	graphName == "cent6_pi+_Au+Au_27"
			||	graphName == "cent2_ka-_Au+Au_27"
			||	graphName == "cent3_ka+_Au+Au_27"
			||	graphName == "cent7_pbar_Au+Au_27"
			||	graphName == "cent8_pi+_Au+Au_39"
			
			||	graphName == "cent0_la_Au+Au_7.7"
						
			// still not working:
			// ||	graphName == "cent8_pbar_Au+Au_7.7" /// ------ FAILED
			)
			{
			cout << "graphname is: " << graphName << endl;
			funcBGBW->SetParameters(mass,0.9,0.03,0.01,10000.,type);
			cout << "alternate init pars: 0.9,0.03,0.01,10000." << endl;
			}
			
			
			
			
		else if(graphName == "cent0_pis_Pb+Pb_2.76"
			||	graphName == "cent8_pbar_Au+Au_7.7"
			)
			{
			cout << "check2" << endl;
			funcBGBW->SetParameters(mass,0.99,0.30,0.1,1000.,type);
			}
			
			
			
		
		/////...... may have to uncomment above code and edit if below...
		// ... code doesn't work		
		else{
			cout << "graphname is: " << graphName << endl;
			funcBGBW->SetParameters(mass,0.95,0.05,0.1,1000000.,type);
		}
		funcBGBW->SetParNames("mass","beta (c)","temp","n","norm","type");
		///.... changed the lower limit of beta from 0.5 to 0.2
		funcBGBW->SetParLimits(1,0.2,0.999999999999999999999);//param 1

		funcBGBW->FixParameter(0,mass);// mass in GeV
		funcBGBW->FixParameter(5,type);
		ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(20000);
		TFitResultPtr r = tg->Fit("getdNdpt","S","",0.00000000000001,5.);
		Double_t chi2Prob = r->Prob();
		cout << "chi-sq prob: " << chi2Prob << endl;
		////////........h->SetMaximum(5*(h->GetMaximum()));
		//h-> GetYaxis()->SetRangeUser(0.,maxY);
		///TMatrixDSym cov = r->GetCovarianceMatrix();
		//////.......h-> GetXaxis()->SetRangeUser(0.,30.);
		TString xlabel = "p_{T}";
		TString ylabel = "#frac{d^{2}N}{dydp_{T}}";
		tg-> GetHistogram()->GetXaxis()->SetTitle(xlabel);/////.....SetXTitle(xlabel);
		tg-> GetHistogram()->GetYaxis()->SetTitle(ylabel);
		Double_t beta 			= funcBGBW->GetParameter(1);
		Double_t temp 			= funcBGBW->GetParameter(2);
		Double_t n	  			= funcBGBW->GetParameter(3);
		Double_t norm 			= funcBGBW->GetParameter(4);
		Double_t betaErr 		= funcBGBW->GetParError(1);
		Double_t tempErr 		= funcBGBW->GetParError(2);
		Double_t nErr 			= funcBGBW->GetParError(3);
		Double_t normErr 		= funcBGBW->GetParError(4);
		tg-> Draw("AP");
		//------------- end BGBW fit ----------------------------
		
		
		//-------- Find integrals left and right of data points -------//
		funcBGBW			 	-> SetParameters(mass,beta,temp,n,norm,type);
		dETdEtaIntegrandFunc 	-> SetParameters(mass,beta,temp,n,norm,type);
		dETdEtaIntegrandFunc	-> FixParameter(5,type);
		dETdEtaIntegrandFunc	-> FixParameter(0,mass);
		dETdyIntegrandFunc 		-> SetParameters(mass,beta,temp,n,norm,type);
		dETdyIntegrandFunc		-> FixParameter(5,type);
		dETdyIntegrandFunc		-> FixParameter(0,mass);
		dNdEtaIntegrandFunc 	-> SetParameters(mass,beta,temp,n,norm,type);
		dNdEtaIntegrandFunc		-> FixParameter(0,mass);
		dNdEtaIntegrandFunc		-> FixParameter(5,type);
		dNdyIntegrandFunc 		-> SetParameters(mass,beta,temp,n,norm,type);
		dNdyIntegrandFunc		-> FixParameter(0,mass);
		dNdyIntegrandFunc		-> FixParameter(5,type);

		///////....Int_t totBins 	= h->GetNbinsX();
		Int_t totPoints = tg->GetN();		
		Int_t binx1 	= 0;
		Int_t binx2 	= totPoints+1;
		
		/////.....Double_t leftCut 	= h->GetXaxis()->GetBinLowEdge(binx1+2); // TODO figure out why +2 with cout
		//Double_t leftCut	= (tg->GetPoint(1, *(tg->GetX()), *(tg->GetY())));//->GetX();
		Double_t leftCut;
		Double_t leftCutY;
		tg->GetPoint(0, leftCut, leftCutY);
		cout << "totPoints: " << totPoints << endl;
		cout << "leftcut: " << leftCut << endl;
		/////.....Double_t rightCut 	= h->GetXaxis()->GetBinUpEdge(binx2-1); // TODO figure out why -1 with cout
		Double_t rightCut;
		Double_t rightCutY;
		tg->GetPoint(totPoints-1, rightCut, rightCutY);
		cout << "rightcut: " << rightCut << endl; 		
		Double_t dETdEtaLeft 	= dETdEtaIntegrandFunc -> Integral(0.,leftCut);
		Double_t dETdEtaRight 	= dETdEtaIntegrandFunc -> Integral(rightCut,30.);
		Double_t dETdyLeft 		= dETdyIntegrandFunc -> Integral(0.,leftCut);
		Double_t dETdyRight 	= dETdyIntegrandFunc -> Integral(rightCut,30.);
		Double_t dNdEtaLeft 	= dNdEtaIntegrandFunc -> Integral(0.,leftCut);
		Double_t dNdEtaRight 	= dNdEtaIntegrandFunc -> Integral(rightCut,30.);
		Double_t dNdyLeft 		= dNdyIntegrandFunc -> Integral(0.,leftCut);
		Double_t dNdyRight 		= dNdyIntegrandFunc -> Integral(rightCut,30.);
		// Errors:
		Double_t dETdEtaLErr	= 
		dETdEtaIntegrandFunc->IntegralError(0.,leftCut,
								r->GetParams(),
								r->GetCovarianceMatrix().GetMatrixArray());
		Double_t dETdEtaRErr	= 
		dETdEtaIntegrandFunc->IntegralError(rightCut,30.,
								r->GetParams(),
								r->GetCovarianceMatrix().GetMatrixArray());
		Double_t dETdyLErr	= 
		dETdyIntegrandFunc->IntegralError(0.,leftCut,
								r->GetParams(),
								r->GetCovarianceMatrix().GetMatrixArray());
		Double_t dETdyRErr	= 
		dETdyIntegrandFunc->IntegralError(rightCut,30.,
								r->GetParams(),
								r->GetCovarianceMatrix().GetMatrixArray());
		Double_t dNdEtaLErr	= 
		dNdEtaIntegrandFunc->IntegralError(0.,leftCut,
								r->GetParams(),
								r->GetCovarianceMatrix().GetMatrixArray());
		Double_t dNdEtaRErr	= 
		dETdEtaIntegrandFunc->IntegralError(rightCut,30.,
								r->GetParams(),
								r->GetCovarianceMatrix().GetMatrixArray());
		Double_t dNdyLErr	= 
		dNdyIntegrandFunc->IntegralError(0.,leftCut,
								r->GetParams(),
								r->GetCovarianceMatrix().GetMatrixArray());
		Double_t dNdyRErr	= 
		dNdyIntegrandFunc->IntegralError(rightCut,30.,
								r->GetParams(),
								r->GetCovarianceMatrix().GetMatrixArray());
								
		//////......Double_t dETdEta_d = *(integralDataPtr+0);
		/////.....
		/*
	Double_t dETdEta_d = *(integralData);
		Double_t dETdEta_d_err = *(integralDataPtr+1);
		Double_t dETdEtaTotal = dETdEtaLeft+dETdEta_d+dETdEtaRight;
		Double_t dETdEtaTErr = dETdEtaLErr+dETdEta_d_err+dETdEtaRErr;
		
		Double_t dETdy_d = *(integralDataPtr+2);
		Double_t dETdy_d_err = *(integralDataPtr+3);
		Double_t dETdyTotal = dETdyLeft+dETdy_d+dETdyRight;
		Double_t dETdyTErr = dETdyLErr+dETdy_d_err+dETdyRErr;
		
		Double_t dNdEta_d = *(integralDataPtr+4);
		Double_t dNdEta_d_err = *(integralDataPtr+5);
		Double_t dNdEtaTotal = dNdEtaLeft+dNdEta_d+dNdEtaRight;
		Double_t dNdEtaTErr = dNdEtaLErr+dNdEta_d_err+dNdEtaRErr;
		
		Double_t dNdy_d = *(integralDataPtr+6);
		Double_t dNdy_d_err = *(integralDataPtr+7);
		Double_t dNdyTotal = dNdyLeft+dNdy_d+dNdyRight;
		Double_t dNdyTErr = dNdyLErr+dNdy_d_err+dNdyRErr;
		
		cout <<"Integral from data for "<<graphName<<": "<<*(integralDataPtr+0)<<endl;// should be 363.7 for pi minus cent 0
		*/		
		cout << "dETdEtaLeft: " << dETdEtaLeft << endl;
		cout << "dETdEtaRight: " << dETdEtaRight << endl;
		cout << "dETdEtaData: " << dETdEtaIntegrandFunc -> Integral(leftCut,rightCut) << endl;
		cout << "dETdEtaTotal: " << dETdEtaLeft+dETdEtaRight+ 
			dETdEtaIntegrandFunc -> Integral(leftCut,rightCut)<< endl;
		cout<<"-----------------------------------"<<endl;				
		//------ end Find integrals left and right of data points ----//
		//------ begin - assign Npart and errors from BES paper -----//
		Int_t* NpartAndArrPtr;
		Int_t Npart;
		Int_t NpartErr;
		NpartAndArrPtr = getNpartAndErr(collEn,centrality);
		Npart = *(NpartAndArrPtr+0);
		NpartErr = *(NpartAndArrPtr+1);
		//------ end - assign Npart and errors from BES paper -------// 
		//-- Output results to file-----------------------------
		datFile << collEn << "\t"	
				<< particleID << "\t"
				<< centrality << "\t"
				<< mass << "\t"
				<< beta <<"\t"
				<< betaErr <<"\t"
				<< temp <<"\t"
				<< tempErr <<"\t"
				<< n <<"\t"
				<< nErr <<"\t"
				<< norm <<"\t"
				<< normErr <<"\t"
				<< 7777 << "\t"//////......<< dETdEta_d << "\t" //dETdEta_d
				<< 7777 << "\t"//////......<< dETdEta_d_err<< "\t" //dETdEta_d_err
				<< dETdEtaLeft << "\t"
				<< dETdEtaLErr << "\t"
				<< dETdEtaRight << "\t"
				<< dETdEtaRErr << "\t"
				<< 7777 << "\t"//////......<< dETdEtaTotal<< "\t" // dETdEtaTotal
				<< 7777 << "\t"//////......<< dETdEtaTErr << "\t"
				<< 7777 << "\t"//////......<< dETdy_d << "\t" //dETdy_d
				<< 7777 << "\t"//////......<< dETdy_d_err<< "\t" //dETdy_d_err
				<< dETdyLeft << "\t"
				<< dETdyLErr << "\t"
				<< dETdyRight << "\t"
				<< dETdyRErr << "\t"
				<< 7777 << "\t"//////......<< dETdyTotal<< "\t" // dETdyTotal
				<< 7777 << "\t"//////......<< dETdyTErr << "\t"
				<< 7777 << "\t"//////......<< dNdEta_d << "\t" //dNdEta_d
				<< 7777 << "\t"//////......<< dNdEta_d_err<< "\t" //dNdEta_d_err
				<< dNdEtaLeft << "\t"
				<< dNdEtaLErr << "\t"
				<< dNdEtaRight << "\t"
				<< dNdEtaRErr << "\t"
				<< 7777 << "\t"//////......<< dNdEtaTotal << "\t" // dNdEtaTotal
				<< 7777 << "\t"//////......<< dNdEtaTErr << "\t"
				<< 7777 << "\t"//////......<< dNdy_d << "\t" //dNdy_d
				<< 7777 << "\t"//////......<< dNdy_d_err<< "\t" //dNdy_d_err
				<< dNdyLeft << "\t"
				<< dNdyLErr << "\t"
				<< dNdyRight << "\t"
				<< dNdyRErr << "\t"
				<< 7777 << "\t"//////......<< dNdyTotal << "\t" // dNdyTotal
				<< 7777 << "\t"//////......<< dNdyTErr << "\t"
				<< Npart << "\t"
				<< NpartErr << "\n";
		
		//-- end- output results to file------------------------
		c1->Update();
		Double_t chi2BGBW = funcBGBW->GetChisquare();
		Double_t nDFBGBW = funcBGBW->GetNDF();
		Double_t p2 = funcBGBW->GetParameter(2);
		Double_t e2 = funcBGBW->GetParError(2);


	
		//cout << "chi2: " << chi2BGBW << "\nndf: " 
			//<< nDFBGBW<< "\nchi2/ndf: " << chi2BGBW/nDFBGBW <<endl;
	
		/* FIXME */
		time_t secFromEpoch;
		string imgPathAndName = "./debugPlots/"+graphName+to_string(time(&secFromEpoch))+".png";
				//c1 -> SaveAs("./fittedPlots/trial1.png");
		TImage *png = TImage::Create();// FIXME try to use canvas method instead of png object
		png->FromPad(c1);
		const char* imgPathAndNameConstCharPtr = imgPathAndName.c_str();
		png->WriteImage(imgPathAndNameConstCharPtr);
		/**/
		//cout << "Draw class here: \n";
		//h-> DrawClass();
		///////h->Delete();// works
		///////////FIXME c1->Clear();// 
		/// sometimes when you delete objects, they stay in the program stack
		//FIXME delete png;
		mikey->DeleteBuffer();// works!
		
		breakOutForTesting++;
		if(breakOutForTesting>=stop) break;
		
		gSystem->ProcessEvents();
		/////.....delete h;
		delete funcBGBW;
		delete dETdEtaIntegrandFunc;
		delete dETdyIntegrandFunc;
		delete tg;
		delete c1;	// Rademakers
		//delete mikey; // FIXME 9 segmentation violation
		//delete class1; // segmentation violation
	}// end of while loop to iterate through every key
	///} // end of while(1)
	/////////////delete c1;
	//delete mikey;
	//delete h;
	//delete class1;
	//delete funcBGBW;
	//delete dETdEtaIntegrandFunc;
	//delete dETdyIntegrandFunc;
	//gObjectTable->Print();
	delete myFile;
	datFile.close();
return 0;
}






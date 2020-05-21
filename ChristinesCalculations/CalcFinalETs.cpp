#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;



struct particle {
  double dETdEta;
  double dETdEta_err;
  double dETdEta_err_extrap;
  double dNdEta;
  double dNdEta_err;
  double dNdEta_err_extrap;
  double npart;
  double npart_err;
  //double dETdEta_exterr;
  //double dNdEta_exterr;
};


int main(){
  //value
  double ETpiplus =1;
  //systematic uncertainty from experimental errors
  double ETpiplusSysExp = 0.1;
  //systematic uncertainty from experimental errors
  double ETpiplusSysExtrap = 0.00;
  //value
  double ETpiminus =1;
  //systematic uncertainty from experimental errors
  double ETpiminusSysExp = 0.1;
  //systematic uncertainty from experimental errors
  double ETpiminusSysExtrap = 0.00;
  //value
  double ETKplus =1;
  //systematic uncertainty from experimental errors
  double ETKplusSysExp = 0.1;
  //systematic uncertainty from experimental errors
  double ETKplusSysExtrap = 0.00;
  //value
  double ETKminus =1;
  //systematic uncertainty from experimental errors
  double ETKminusSysExp = 0.1;
  //systematic uncertainty from experimental errors
  double ETKminusSysExtrap = 0.00;
  //value
  double ETproton =1;
  //systematic uncertainty from experimental errors
  double ETprotonSysExp = 0.1;
  //systematic uncertainty from experimental errors
  double ETprotonSysExtrap = 0.00;
  //value
  double ETantiproton =1;
  //systematic uncertainty from experimental errors
  double ETantiprotonSysExp = 0.1;
  //systematic uncertainty from experimental errors
  double ETantiprotonSysExtrap = 0.00;
  //value
  double ETlambda =1;
  //systematic uncertainty from experimental errors
  double ETlambdaSysExp = 0.1;
  //systematic uncertainty from experimental errors
  double ETlambdaSysExtrap = 0.00;
  //value
  double ETantilambda =1;
  //systematic uncertainty from experimental errors
  double ETantilambdaSysExp = 0.1;
  //systematic uncertainty from experimental errors
  double ETantilambdaSysExtrap = 0;
  //value
  double ETetaomega =.00000001;
  //systematic uncertainty - we're just rolling this into one because the scale factor uncertainty is dominant
  double ETetaomegaSys = 0.1;

  double fpi=1.638;//1.588;
  double fpiErr=0.055;//0.029;
  double fK=1.85; //fk- = 1.5+-0.1        fk+ = 0.5 +- 0.05
  double fKErr=0.15;
  //double fp=1.26; //add in for energy dependance
  //double fpErr=0.26;
  double fLam=1.07;//.4373948
  double fLamErr=0.50;//.12
  double feta=0;//.409; //.111
  double fetaErr=0;//.023237;
  double fome=0;//.0127;//.147
  double fomeErr=0;//.0445;
  double fexotic = 0;//0.0137;
  double fSigma = 1.57; //(ET_sigplus + ET_sigminus + ETsig0)/(ET_piplus + ET_piminus) (difference between sigma species is sigplus~=0.007909 sigminus~=0.008591 sig0~=0.009361)
  double fSigmaErr = 0.07;
  double fp = 2.2658229;
  double fpErr = 0.26582277;
  //vector <double> fp = {2.26,2.25,2.23,2.21,2.18,2.14,2.10,2.05,2.03};
  //vector <double> fpErr = {0.26,0.25,0.23,0.21,0.18,0.14,0.10,0.05,0.03};
  //The last three rows are placeholders because I just copied yields out of the BES spectra paper
  double antiprotonyields[8][9] = {  {0.39,0.32,0.26,0.19,0.14,0.09,0.06,0.033,0.018},
				     {1.5,1.2,0.9,0.7,0.5,0.33,0.21,0.13,0.07},
				     {4.2,3.4,2.7,1.9,1.4,0.95,0.6,0.35,0.18},
				     {6,5.1,4,2.9,2,1.4,0.8,0.49,0.23},
				     {8.5,7.4,5.4,3.9,2.8,1.8,1.2,0.64,0.33},
				     {13.6,11.4,8.77,6.39,4.27,2.77,1.68,0.96,0.464},
				     {20, 15.7, 12.8, 10.5, 8.02, 5.51, 3.33, 1.31, 0},
				     {26.7,21.4,15.7,11.2,7.46,4.93,3.16,1.84,0.915}};
  double antiprotonyieldserr[8][9] = {  {0.05,0.04,0.03,0.02,0.02,0.01,0.01,0.004,0.002},
					{0.2,0.2,0.1,0.1,0.1,0.04,0.03,0.02,0.01},
					{0.5,0.4,0.4,0.3,0.2,0.1,0.1,0.05,0.02},
					{0.7,0.6,0.5,0.3,0.2,0.2,0.1,0.05,0.03},
					{1,0.9,0.7,0.5,0.3,0.2,0.1,0.08,0.04},
					{1.7,1.1,0.78,0.55,0.35,0.19,0.12,0.059,0.047},
					{2.2,1.6,1.6,1,0.81,0.45,0.3,0.09,0.09},
					{3.4,2.5,1.7,1.1,0.72,0.46,0.29,0.16,0.081}};
  double protonyields[8][9] = {{54.9,45.4,33.4,23.2,15.8,9.3,5.4,2.8,1.4},
			       {44,35.2,26.1,17.8,11.8,7.3,4.2,2.1,1},
			       {34.2,29.3,21.9,14.6,9.2,5.8,3.3,1.8,0.8},
			       {31.7,26.5,19.4,12.9,8.9,5.6,3.2,1.7,0.8},
			       {26.5,22.7,17.3,11.9,7.9,4.9,2.9,1.5,0.7},
			       {29,23.8,17.8,12.2,8.08,5.07,2.98,1.6,0.745},
			       {28.2,21.9,17.9,14.4,10.9,7.35,4.38,1.65,1.65},
			       {34.7,28.2,20.1,14.4,9.3,6.17,3.88,2.2,1.09}};
  double protonyieldserr[8][9] = {  {6.1,5,3.7,2.6,1.7,1,0.6,0.3,0.2},
					{5.3,4.2,3.1,2.1,1.4,0.9,0.5,0.3,0.1},
					{4.5,3.8,2.9,1.9,1.2,0.8,0.4,0.2,0.1},
					{3.8,3.2,2.3,1.5,1.1,0.7,0.4,0.2,0.1},
					{2.9,2.5,1.9,1.3,0.9,0.5,0.3,0.2,0.1},
					{3.8,2.4,1.6,1.1,0.67,0.36,0.22,0.12,0.086},
					{3.1,2.3,2.2,1.3,1.1,0.6,0.39,0.11,},
					{4.4,3.3,2.2,1.4,0.89,0.57,0.35,0.2,0.1}};
  double pionminusyield[8][9] = {{100,81.9,62.9,43.3,29.1,18.8,11.8,6.6,3.4},
				 {129.8,102.3,77,52,35.7,22.5,13.6,7.9,4.2},
				 {165.8,133.7,102.1,68.8,46,28.9,17.6,9.7,5.2},
				 {177.1,147.5,111.6,75.9,49.9,31.5,18.9,10.6,5.3},
				 {185.8,155,118.4,80.7,52.9,33.7,20.6,11.3,6.1},
				 {237,192,146,101,67.4,43.7,26.8,14.7,7.43},
				 {280,228,187,140,104,70.9,42.4,16,16},
				 {327.00,261,196,136,89.6,58.9,36.3,21.1,10.9}};
  double pionminusyielderr[8][9] = {{9,7.4,5.7,3.9,2.6,1.7,1.1,0.6,0.3},
				    {13,10.3,7.7,5.2,3.6,2.3,1.4,0.8,0.4},
				    {18.3,14.7,11.3,7.6,5.1,3.2,1.9,1.1,0.6},
				    {19.5,16.3,12.3,8.4,5.5,3.5,2.1,1.2,0.6},
				    {20.5,17.1,13.1,8.9,5.8,3.7,2.2,1.2,0.7},
				    {17,13,11,7,5.2,3.5,2.4,1.3,0.62},
				    {20,16,16,11,8,4.9,3.5,2.1,2.1},
				    {25,20,15,10,6.8,4.5,2.8,1.6,0.8}};
  double pionplusyield[8][9] = {{93.4,76.8,58.7,40.5,26.9,17.6,10.9,6.1,3.1},
				{123.9,97.1,73.4,49.5,33.9,21.3,12.9,7.6,3.9},
				{161.4,130.3,99.3,67.1,44.8,28.1,17.1,9.5,5},
				{172.9,144.3,109.4,74.3,48.8,30.7,18.6,10.4,5.1},
				{182.3,151.4,115.9,78.9,51.8,32.9,20.1,11,5.9},
				{233,191,144,98.9,66.5,43.2,26.5,14.8,7.34},
				{278,228,186,140,103,71.8,42.2,16,16},
				{322,257,194,135,89.2,58.7,36.2,21.1,10.8}};
double pionplusyielderr[8][9] = {{8.4,6.9,5.3,3.7,2.4,1.6,0.9,0.6,0.3},
				 {12.4,9.7,7.4,4.9,3.4,2.1,1.3,0.8,0.4},
				 {17.8,14.4,10.9,7.4,4.9,3.1,1.9,1,0.6},
				 {19.1,15.9,12.1,8.2,5.4,3.4,2,1.1,0.6},
				 {20.1,16.7,12.8,8.7,5.7,3.6,2.2,1.2,0.7},
				 {17,13,11,6.9,5.1,3.5,2.3,1.3,0.62},
				 {20,16,16,11,8,5,3.5,1.9,1.9},
				 {25,20,15,10,6.8,4.5,2.7,1.6,0.8}};
double kaonplusyield[8][9] = {{20.8,17.3,12.4,8.6,5.3,3.2,1.8,0.82,0.33},
			      {25,20.6,14.8,9.6,6.1,3.7,1.9,0.98,0.46},
			      {29.6,24.3,18,12.3,7.8,4.7,2.7,1.3,0.65},
			      {31.1,25.8,19.4,12.9,8.3,5.2,2.9,1.5,0.68},
			      {32,27,20.3,13.6,8.8,5.4,3.2,1.6,0.8},
			      {37.6,31.2,23,15.9,10.4,6.62,3.64,1.95,0.868},
			      {46.3,35.6,29,22.3,16.4,11.2,6.83,2.31,2.31},
			      {51.3,40.8,30,20.5,13.6,8.69,5.4,2.98,1.41}};
double kaonplusyielderr[8][9] = {{1.7,1.4,1,0.7,0.4,0.3,0.1,0.07,0.03},
				 {2.5,2.1,1.5,1,0.6,0.4,0.2,0.09,0.05},
				 {2.9,2.4,1.8,1.2,0.8,0.5,0.3,0.1,0.06},
				 {2.8,2.3,1.8,1.2,0.8,0.5,0.3,0.1,0.06},
				 {2.9,2.4,1.8,1.2,0.8,0.5,0.3,0.1,0.07},
				 {2.7,2.2,1.6,1.1,0.7,0.46,0.25,0.13,0.058},
				 {3,2.6,2.1,1.9,1.4,1,0.48,0.15,0.15},
				 {6.5,4.7,3.2,2,1.3,0.81,0.49,0.27,0.13}};
double kaonpminusyield[8][9] = {{7.7,6.4,4.7,3.2,2.1,1.3,0.71,0.32,0.13},
				{12.3,10.2,7.5,4.9,3.2,1.9,1,0.53,0.25},
				{18.8,15.5,11.6,7.9,5.2,3.2,1.8,0.9,0.45},
				{22.6,18.7,14.5,9.8,6.2,3.9,2.2,1.1,0.51},
				{25,21,15.9,10.7,7,4.4,2.6,1.3,0.7},
				{32.4,27.2,19.8,14,8.89,5.68,3.31,1.74,0.813},
				{42.7,33.1,26.6,20.5,15,10.1,5.81,2.23,2.23},
				{49.5,39.8,28.7,19.7,13.2,8.37,5.19,2.89,1.38}};
double kaonminusyielderr[8][9] = {{0.6,0.5,0.4,0.3,0.2,0.1,0.06,0.03,0.01},
				  {1.2,1,0.7,0.5,0.3,0.2,0.1,0.05,0.03},
				  {1.9,1.6,1.2,0.8,0.5,0.3,0.2,0.1,0.04},
				  {2,1.7,1.3,0.9,0.6,0.3,0.2,0.1,0.05},
				  {2.3,1.9,1.4,1,0.6,0.4,0.2,0.1,0.1},
				  {2.3,1.9,1.4,1,0.62,0.39,0.23,0.12,0.055},
				  {2.8,2.4,1.9,1.8,1.3,0.9,0.41,0.14,0.14},
				  {6.2,4.6,3.1,2,1.3,0.78,0.47,0.26,0.13}};

  double cbedges[] = {0,5,10,20,30,40,50,60,70,80};
  double jacobian = 1.0;
  
  double jacobiancorr = 1.023;//Only applies to 62-200 GeV
  double jacobiancorrerr = 0.001;

  double protonmass = 0.938;

  string GRBG;
  double grbgd;
  ifstream in;
  particle q;
  int numcolumns = 9;
  int a = 9;//number of centralities
  int E  =8;//number of energies
  cout<<"I'm Running!\n";
  vector < vector< particle > > ETpip (E, vector<particle>(a));
  vector < vector< particle > > ETpim (E,vector<particle>(a));
  vector < vector< particle > > ETKp (E,vector<particle>(a));
  vector < vector< particle > > ETKm (E,vector<particle>(a));
  vector < vector< particle > > ETp (E,vector<particle>(a));
  vector < vector< particle > > ETap (E,vector<particle>(a));
  vector < vector< particle > > ETeta (E,vector<particle>(a));
  vector < vector< particle > > ETome (E,vector<particle>(a));
  vector < vector< particle > > ETLa (E,vector<particle>(a));
  vector < vector< particle > > ETLab (E,vector<particle>(a));
  vector <float> SNN = {7.7,11.5,19.6,27,39,62.4,130,200};
  int cent;
  double dat1,dat2,dat3,dat4,dat5,dat6,dat7;
  in.open("avgfitResults.dat");
  for (int i = 0; i<numcolumns;i++){//read in header, which we don't use
    in>>GRBG;
    //cout<<" i "<<i<<" "<<GRBG;
  }
  cout<<endl;
  for (int iIter=0;iIter<E;iIter++){//loop over all energies
    //The 200 and 130 GeV are switched
    int i=iIter;
    if(iIter==E-1){ 
      //cout<<"This is really 130 GeV"<<endl;
      i=E-2;//This is the 130 GeV
    }
    if(iIter==E-2){
      //cout<<"This is really 200 GeV"<<endl;
      i=E-1;//This is the 200 GeV
    }
    for(int j=0;j<a;j++){//loop over all centralities, read in pi-
      if(iIter==E-1 && j==a-1){//this fills the last bin with a negative ET so that it doesn't appear on the plot
	//cout<<"iIter "<<iIter<<" i "<<i<<" j "<<j<<" a-1 "<<a-1<<endl;
	ETpim[i][j].dETdEta=-dat1;
	//cout<<"Only entering here!"<<endl;
      }
      else{
	in>>grbgd>>GRBG>>cent>>dat1>>dat2>>dat3>>dat4>>dat5>>dat6>>dat7;//>>dat8;
	ETpim[i][j].dETdEta=dat1;
	ETpim[i][j].dETdEta_err=dat2;
	ETpim[i][j].dETdEta_err_extrap=dat3;
	ETpim[i][j].dNdEta=dat4;
	ETpim[i][j].dNdEta_err=dat5;
	ETpim[i][j].npart=dat6;
	ETpim[i][j].npart_err=dat7;
      //ETpim[i][j].dETdEta_exterr=dat7;
      //ETpim[i][j].dNdEta_exterr=dat8;
	//if(j==0){
	if(cent!=j) cout<<"Warning! Failure here"<<endl;
	//cout<<"Energy "<<i<<" "<<grbgd<<" "<<GRBG<<" cent "<<cent<<" pion ET "<<dat1<<" "<<ETpim[i][j].dETdEta<<" a "<<a<<" j "<<j<<endl;
	//cout<<" "<<grbgd<<" "<<GRBG<<" "<<cent<<" "<<dat1<<" "<<dat2<<" "<<dat3<<" "<<dat4<<" "<<dat5<<" "<<dat6<<" "<<dat7<<endl;
	  //}
      }
    }
    for(int j=0;j<a;j++){//loop over all centralities, read in pi+
      if(iIter==E-1 && j==a-1){//this fills the last bin with a negative ET so that it doesn't appear on the plot
	//cout<<"Only entering here!"<<endl;
	ETpip[i][j].dETdEta=-dat1;
      }
      else{
	in>>grbgd>>GRBG>>cent>>dat1>>dat2>>dat3>>dat4>>dat5>>dat6>>dat7;//>>dat8;
	ETpip[i][j].dETdEta=dat1;
	ETpip[i][j].dETdEta_err=dat2;
	ETpip[i][j].dETdEta_err_extrap=dat3;
	ETpip[i][j].dNdEta=dat4;
	ETpip[i][j].dNdEta_err=dat5;
	ETpip[i][j].npart=dat6;
	ETpip[i][j].npart_err=dat7;
	//ETpip[i][j].dETdEta_exterr=dat7;
	//ETpip[i][j].dNdEta_exterr=dat8;
	if(cent!=j) cout<<"Warning! Failure here"<<endl;
	//cout<<"Energy "<<i<<" "<<grbgd<<" "<<GRBG<<" cent "<<cent<<" pion+ ET "<<dat1<<" "<<ETpim[i][j].dETdEta<<" a "<<a<<" j "<<j<<endl;
      }
    }
    for(int j=0;j<a;j++){//pi0, not used
      if(iIter==E-1 && j==a-1){//this fills the last bin with a negative ET so that it doesn't appear on the plot
      }
      else{
	in>>grbgd>>GRBG>>grbgd>>grbgd>>grbgd>>grbgd>>grbgd>>grbgd>>grbgd>>grbgd;
      }
      //cout<<GRBG<<endl;
    }
    for(int j=0;j<a;j++){//loop over all centralities, read in K-
      if(iIter==E-1 && j==a-1){//this fills the last bin with a negative ET so that it doesn't appear on the plot
	ETKm[i][j].dETdEta=-dat1;
      }
      else{
	in>>grbgd>>GRBG>>cent>>dat1>>dat2>>dat3>>dat4>>dat5>>dat6>>dat7;//>>dat8;
	ETKm[i][j].dETdEta=dat1;
	ETKm[i][j].dETdEta_err=dat2;
	ETKm[i][j].dETdEta_err_extrap=dat3;
	ETKm[i][j].dNdEta=dat4;
	ETKm[i][j].dNdEta_err=dat5;
	ETKm[i][j].npart=dat6;
	ETKm[i][j].npart_err=dat7;
	if(cent!=j) cout<<"Warning! Failure here"<<endl;
	//cout<<"Energy "<<i<<" "<<grbgd<<" "<<GRBG<<" cent "<<cent<<" Kaon- ET "<<dat1<<" "<<ETpim[i][j].dETdEta<<" a "<<a<<" j "<<j<<endl;
	//ETKm[i][j].dETdEta_exterr=dat7;
	//ETKm[i][j].dNdEta_exterr=dat8;
      }
    }
    for(int j=0;j<a;j++){//loop over all centralities, read in K+
      if(iIter==E-1 && j==a-1){//this fills the last bin with a negative ET so that it doesn't appear on the plot
	ETKp[i][j].dETdEta=-dat1;
      }
      else{
	in>>grbgd>>GRBG>>cent>>dat1>>dat2>>dat3>>dat4>>dat5>>dat6>>dat7;//>>dat8;
	ETKp[i][j].dETdEta=dat1;
	ETKp[i][j].dETdEta_err=dat2;
	ETKp[i][j].dETdEta_err_extrap=dat3;
	ETKp[i][j].dNdEta=dat4;
	ETKp[i][j].dNdEta_err=dat5;
	ETKp[i][j].npart=dat6;
	ETKp[i][j].npart_err=dat7;
	if(cent!=j) cout<<"Warning! Failure here"<<endl;
	//cout<<"Energy "<<i<<" "<<grbgd<<" "<<GRBG<<" cent "<<cent<<" kaon+ ET "<<dat1<<" "<<ETpim[i][j].dETdEta<<" a "<<a<<" j "<<j<<endl;
	//ETKp[i][j].dETdEta_exterr=dat7;
      //ETKp[i][j].dNdEta_exterr=dat8;
      }
    }
    for(int j=0;j<a;j++){//loop over all centralities, read in antiprotons
      if(iIter==E-1 && j==a-1){//this fills the last bin with a negative ET so that it doesn't appear on the plot
	ETap[i][j].dETdEta=-dat1;
      }
      else{
	in>>grbgd>>GRBG>>cent>>dat1>>dat2>>dat3>>dat4>>dat5>>dat6>>dat7;//>>dat8;
	ETap[i][j].dETdEta=dat1;
	ETap[i][j].dETdEta_err=dat2;
	ETap[i][j].dETdEta_err_extrap=dat3;
	ETap[i][j].dNdEta=dat4;
	ETap[i][j].dNdEta_err=dat5;
	ETap[i][j].npart=dat6;
	ETap[i][j].npart_err=dat7;
      //ETap[i][j].dETdEta_exterr=dat7;
      //ETap[i][j].dNdEta_exterr=dat8;
	if(cent!=j) cout<<"Warning! Failure here"<<endl;
	//cout<<"Energy "<<i<<" "<<grbgd<<" "<<GRBG<<" cent "<<cent<<" antiproton ET "<<dat1<<" "<<ETpim[i][j].dETdEta<<" a "<<a<<" j "<<j<<endl;
      }
    }
    for(int j=0;j<a;j++){//loop over all centralities, read in protons
      if(iIter==E-1 && j==a-1){//this fills the last bin with a negative ET so that it doesn't appear on the plot
	ETp[i][j].dETdEta=-dat1;
	//cout<<"Is this getting flagged here? proton"<<endl;
      }
      else{
	in>>grbgd>>GRBG>>cent>>dat1>>dat2>>dat3>>dat4>>dat5>>dat6>>dat7;//>>dat8;
	ETp[i][j].dETdEta=dat1;
	ETp[i][j].dETdEta_err=dat2;
	ETp[i][j].dETdEta_err_extrap=dat3;
	ETp[i][j].dNdEta=dat4;
	ETp[i][j].dNdEta_err=dat5;
	ETp[i][j].npart=dat6;
	ETp[i][j].npart_err=dat7;
	if(cent!=j) cout<<"Warning! Failure here"<<endl;
	//cout<<"Energy "<<i<<" "<<grbgd<<" "<<GRBG<<" cent "<<cent<<" proton ET "<<dat1<<" "<<ETpim[i][j].dETdEta<<" a "<<a<<" j "<<j<<endl;
	//ETp[i][j].dETdEta_exterr=dat7;
	//ETp[i][j].dNdEta_exterr=dat8;
      }
    }
    for(int j=0;j<a;j++){//loop over all centralities, read in eta
      if(iIter==E-1 && j==a-1){//this fills the last bin with a negative ET so that it doesn't appear on the plot
	ETeta[i][j].dETdEta=-dat1;
	//cout<<"Is this getting flagged here? eta"<<endl;
      }
      else{
	in>>grbgd>>GRBG>>cent>>dat1>>dat2>>dat3>>dat4>>dat5>>dat6>>dat7;//>>dat8;
	ETeta[i][j].dETdEta=dat1;
	ETeta[i][j].dETdEta_err=dat2;
	ETeta[i][j].dETdEta_err_extrap=dat3;
	ETeta[i][j].dNdEta=dat4;
	ETeta[i][j].dNdEta_err=dat5;
	ETeta[i][j].npart=dat6;
	ETeta[i][j].npart_err=dat7;
	if(cent!=j) cout<<"Warning! Failure here"<<endl;
	//cout<<"Energy "<<i<<" "<<grbgd<<" "<<GRBG<<" cent "<<cent<<" eta ET "<<dat1<<" "<<ETpim[i][j].dETdEta<<" a "<<a<<" j "<<j<<endl;
	//ETeta[i][j].dETdEta_exterr=dat7;
	//ETeta[i][j].dNdEta_exterr=dat8;
      }
    }
    for(int j=0;j<a;j++){//loop over all centralities, read in omega
      if(iIter==E-1 && j==a-1){//this fills the last bin with a negative ET so that it doesn't appear on the plot
	ETome[i][j].dETdEta=-dat1;
      }
      else{
	in>>grbgd>>GRBG>>cent>>dat1>>dat2>>dat3>>dat4>>dat5>>dat6>>dat7;//>>dat8;
	ETome[i][j].dETdEta=dat1;
	ETome[i][j].dETdEta_err=dat2;
	ETome[i][j].dETdEta_err_extrap=dat3;
	ETome[i][j].dNdEta=dat4;
	ETome[i][j].dNdEta_err=dat5;
	ETome[i][j].npart=dat6;
	ETome[i][j].npart_err=dat7;
	if(cent!=j) cout<<"Warning! Failure here"<<endl;
	//cout<<"OMEGALOOP Energy "<<i<<" "<<grbgd<<" "<<GRBG<<" cent "<<cent<<" omega ET "<<dat1<<" "<<ETpim[i][j].dETdEta<<" a "<<a<<" j "<<j<<endl;
	//ETome[i][j].dETdEta_exterr=dat7;
	//ETome[i][j].dNdEta_exterr=dat8;
      }
    }
  }
  in.close();
  in.open("avgLAMResults.dat");
  for (int i = 0; i<a;i++){
    in>>GRBG;
  }
  a=7;
  for (int i=0;i<E;i++){
    if(i>5) a=5;
    for(int j=0;j<a;j++){
      in>>grbgd>>GRBG>>cent>>dat1>>dat2>>dat3>>dat4>>dat5>>dat6; //>>dat7>>dat8;
      ETLa[i][j].dETdEta=dat1;
      ETLa[i][j].dETdEta_err=dat2;
      ETLa[i][j].dNdEta=dat3;
      ETLa[i][j].dNdEta_err=dat4;
      ETLa[i][j].npart=dat5;
      ETLa[i][j].npart_err=dat6;
      //ETLa[i][j].dETdEta_exterr=dat7;
      //ETLa[i][j].dNdEta_exterr=dat8;
    }
    for(int j=0;j<a;j++){
      in>>grbgd>>GRBG>>cent>>dat1>>dat2>>dat3>>dat4>>dat5>>dat6; //>>dat7>>dat8;
      ETLab[i][j].dETdEta=dat1;
      ETLab[i][j].dETdEta_err=dat2;
      ETLab[i][j].dNdEta=dat3;
      ETLab[i][j].dNdEta_err=dat4;
      ETLab[i][j].npart=dat5;
      ETLab[i][j].npart_err=dat6;
      //ETLab[i][j].dETdEta_exterr=dat7;
      //ETLab[i][j].dNdEta_exterr=dat8;
    }
  }
  in.close();
  ofstream out;
  out.open("EtResults.txt");
  ofstream out2;
  out2.open("EtResultsN.txt");
  ofstream out3;
  out3.open("EtRatios.txt");

  double nominalETdETdEta;
  double nominalETdNdEta;

  double factorVariancedETdEta;
  double factorVariancedNdEta;

  double factorSysUncertaintydETdEta;
  double factorSysUncertaintydNdEta;

  double expVariancedETdEta;
  double expVariancedNdEta;

  double expSysUncertaintydETdEta;
  double expSysUncertaintydNdEta;

  double expSysUncertaintyCorrelateddETdEta;
  double expSysUncertaintyCorrelateddNdEta;

  double extrapVariancedETdEta;
  double extrapVariancedNdEta;

  double extrapSysUncertaintydETdEta;
  double extrapSysUncertaintydNdEta;

  double extrapSysUncertaintyCorrelateddETdEta;
  double extrapSysUncertaintyCorrelateddNdEta;

  double extrapVariance =fpi*(ETpiplusSysExtrap*ETpiplusSysExtrap+ETpiminusSysExtrap*ETpiminusSysExtrap)+ fp*(ETprotonSysExtrap*ETprotonSysExtrap+ETantiprotonSysExtrap*ETantiprotonSysExtrap)+fK*(ETKminusSysExtrap*ETKminusSysExtrap+ETKplusSysExtrap*ETKplusSysExtrap)+fLam*(ETantilambdaSysExtrap*ETantilambdaSysExtrap+ETlambdaSysExtrap*ETlambdaSysExtrap) ;
  double extrapSysUncertainty = sqrt(extrapVariance);
  double extrapSysUncertaintyCorrelated = fpi*(ETpiplusSysExtrap+ETpiminusSysExtrap)+ fp*(ETprotonSysExtrap+ETantiprotonSysExtrap)+fK*(ETKminusSysExtrap+ETKplusSysExtrap)+fLam*(ETantilambdaSysExtrap+ETlambdaSysExtrap);

  //out<<"dETdEta\n";
  //out<<" SNN, bin, corr, nominalET   factorSysUncertainty   expSysUncertainty   extrapSysUncertainty\n";
  //out2<<"dNdEta\n";
  //out2<<" SNN, bin, corr, nominalET   factorSysUncertainty   expSysUncertainty   extrapSysUncertainty\n";
  double beta = 121.0/79.0;
  double jacobianscale = 1.0;
  double jacobianscaleerror = 0.;
  for(int i=0;i<8;i++){
    if(i>=5){
      jacobianscale = jacobiancorr;
      jacobianscaleerror = jacobiancorrerr;
    }
    for(int j=0;j<9;j++){
      double fpALT =( (1.0+beta) +  (3.0-beta)*ETap[i][j].dETdEta/ETp[i][j].dETdEta-(1.0-beta)*2.0*antiprotonyields[i][j]*protonmass/ETp[i][j].dETdEta)/(1.0+ETap[i][j].dETdEta/ETp[i][j].dETdEta);
      double etratioerror = ETap[i][j].dETdEta/ETp[i][j].dETdEta*sqrt((ETp[i][j].dETdEta_err*ETp[i][j].dETdEta_err+ETp[i][j].dETdEta_err_extrap*ETp[i][j].dETdEta_err_extrap)/ETp[i][j].dETdEta/ETp[i][j].dETdEta+(ETap[i][j].dETdEta_err*ETap[i][j].dETdEta_err+ETap[i][j].dETdEta_err_extrap*ETap[i][j].dETdEta_err_extrap)/ETap[i][j].dETdEta/ETap[i][j].dETdEta);
      //I've neglected the correlated error in the denominator but I think that's pretty good.
      double fpalTerr = sqrt((3.0-beta)*(3.0-beta)*etratioerror*etratioerror + (1.0-beta)*2.0*antiprotonyieldserr[i][j]*protonmass/ETp[i][j].dETdEta*(1.0-beta)*2.0*antiprotonyieldserr[i][j]*protonmass/ETp[i][j].dETdEta)/(1.0+ETap[i][j].dETdEta/ETp[i][j].dETdEta);
      //cout<<" & "<<cbedges[j];
      double yieldratio = (1.0+beta + (3.0-beta)*antiprotonyields[i][j]/protonyields[i][j])/(1.0+antiprotonyields[i][j]/protonyields[i][j]);
      double yieldratioerr = yieldratio*sqrt(antiprotonyieldserr[i][j]*antiprotonyieldserr[i][j]/antiprotonyields[i][j]/antiprotonyields[i][j]+protonyieldserr[i][j]*protonyieldserr[i][j]/protonyields[i][j]/protonyields[i][j]);

      //cout<<""<<i<<" "<<j;
      //cout<<" pion: "<<ETpip[i][j].dETdEta<<" "<<ETpim[i][j].dETdEta<<" fpi "<<fpi;
      //cout<<" proton: "<<ETp[i][j].dETdEta<<" "<<ETap[i][j].dETdEta<<" fp "<<fp[i];
      //cout<<" "<<fpALT<<" "<<fpalTerr;
      //cout<<" "<<yieldratio<<" "<<yieldratioerr;
      //cout<<" 1st correction term "<<(3.0-beta)*ETap[i][j].dETdEta/ETp[i][j].dETdEta;
      //cout<<" et ratio "<<ETap[i][j].dETdEta/ETp[i][j].dETdEta<<" +/- "<<etratioerror;
      //cout<<" Npbar "<<antiprotonyields[i][j]<<" second term "<<-(1.0-beta)*2.0*antiprotonyields[i][j]*protonmass/ETp[i][j].dETdEta;
      //cout<<" Np "<<protonyields[i][j]<<" numerator of second term / Np "<<-(1.0-beta)*2.0*antiprotonyields[i][j]*protonmass/protonyields[i][j]<<" Eproton/Nproton "<<ETp[i][j].dETdEta/protonyields[i][j];
      //cout<<" kaon: "<<ETKp[i][j].dETdEta<<" "<<ETKm[i][j].dETdEta<<" fk "<<fK;
      //cout<<" lambda: "<<ETLab[i][j].dETdEta<<" "<<ETLa[i][j].dETdEta<<" fLam "<<fLam;
      nominalETdETdEta=jacobianscale*(fpi*(ETpip[i][j].dETdEta+ETpim[i][j].dETdEta)+ fp*(ETp[i][j].dETdEta+ETap[i][j].dETdEta)+fK*(ETKm[i][j].dETdEta+ETKp[i][j].dETdEta)+fLam*(ETLab[i][j].dETdEta+ETLa[i][j].dETdEta)+feta*ETeta[i][j].dETdEta+fome*ETome[i][j].dETdEta);
      nominalETdNdEta=jacobianscale*(fpi*(ETpip[i][j].dNdEta+ETpim[i][j].dNdEta)+ fp*(ETp[i][j].dNdEta+ETap[i][j].dNdEta)+fK*(ETKm[i][j].dNdEta+ETKp[i][j].dNdEta)+fLam*(ETLab[i][j].dNdEta+ETLa[i][j].dNdEta)+feta*ETeta[i][j].dNdEta+fome*ETome[i][j].dNdEta);
      //cout<<" ET "<<nominalETdETdEta;
      //cout<<endl;

      //1.  Factor uncertainties:
        //Here we are going to treat everything as constant except the factors and we're going to add the uncertainties from the factors as if they are uncorrelated with each other
        //We are going to add in the eta uncertainty because this is largely dominated by the scaling uncertainties
      factorVariancedETdEta = fpiErr*fpiErr*(ETpip[i][j].dETdEta+ETpim[i][j].dETdEta)*(ETpip[i][j].dETdEta+ETpim[i][j].dETdEta)+ fpErr* fpErr*(ETp[i][j].dETdEta+ETap[i][j].dETdEta)*(ETp[i][j].dETdEta+ETap[i][j].dETdEta)+fKErr*fKErr*(ETKm[i][j].dETdEta+ETKp[i][j].dETdEta)*(ETKm[i][j].dETdEta+ETKp[i][j].dETdEta)+fLamErr*fLamErr*(ETLab[i][j].dETdEta+ETLa[i][j].dETdEta)*(ETLab[i][j].dETdEta+ETLa[i][j].dETdEta)+fetaErr*fetaErr*ETeta[i][j].dETdEta+fomeErr*fomeErr*ETome[i][j].dETdEta;
      factorVariancedNdEta = fpiErr*fpiErr*(ETpip[i][j].dNdEta+ETpim[i][j].dNdEta)*(ETpip[i][j].dNdEta+ETpim[i][j].dNdEta)+ fpErr* fpErr*(ETp[i][j].dNdEta+ETap[i][j].dNdEta)*(ETp[i][j].dNdEta+ETap[i][j].dNdEta)+fKErr*fKErr*(ETKm[i][j].dNdEta+ETKp[i][j].dNdEta)*(ETKm[i][j].dNdEta+ETKp[i][j].dNdEta)+fLamErr*fLamErr*(ETLab[i][j].dNdEta+ETLa[i][j].dNdEta)*(ETLab[i][j].dNdEta+ETLa[i][j].dNdEta)+ETetaomegaSys*ETetaomegaSys+fetaErr*fetaErr*ETeta[i][j].dNdEta+fomeErr*fomeErr*ETome[i][j].dNdEta;

      factorSysUncertaintydETdEta = sqrt(jacobianscale*factorVariancedETdEta+jacobianscaleerror*jacobianscaleerror*nominalETdETdEta*nominalETdETdEta);
      factorSysUncertaintydNdEta = sqrt(jacobianscale*factorVariancedNdEta+jacobianscaleerror*jacobianscaleerror*nominalETdNdEta*nominalETdNdEta);

      //2.  Experimental uncertainties
      expVariancedETdEta = fpi*(ETpip[i][j].dETdEta_err*ETpip[i][j].dETdEta_err+ETpim[i][j].dETdEta_err*ETpim[i][j].dETdEta_err)+ fp*(ETp[i][j].dETdEta_err*ETp[i][j].dETdEta_err+ETap[i][j].dETdEta_err*ETap[i][j].dETdEta_err)+fK*(ETKp[i][j].dETdEta_err*ETKp[i][j].dETdEta_err+ETKm[i][j].dETdEta_err*ETKm[i][j].dETdEta_err)+fLam*(ETLa[i][j].dETdEta_err*ETLa[i][j].dETdEta_err+ETLab[i][j].dETdEta_err*ETLab[i][j].dETdEta_err)+feta*ETeta[i][j].dETdEta_err*ETeta[i][j].dETdEta_err+fome*ETome[i][j].dETdEta_err*ETome[i][j].dETdEta_err;
      expVariancedNdEta = fpi*(ETpip[i][j].dNdEta_err*ETpip[i][j].dNdEta_err+ETpim[i][j].dNdEta_err*ETpim[i][j].dNdEta_err)+ fp*(ETp[i][j].dNdEta_err*ETp[i][j].dNdEta_err+ETap[i][j].dNdEta_err*ETap[i][j].dNdEta_err)+fK*(ETKp[i][j].dNdEta_err*ETKp[i][j].dNdEta_err+ETKm[i][j].dNdEta_err*ETKm[i][j].dNdEta_err)+fLam*(ETLa[i][j].dNdEta_err*ETLa[i][j].dNdEta_err+ETLab[i][j].dNdEta_err*ETLab[i][j].dNdEta_err)+feta*ETeta[i][j].dNdEta_err*ETeta[i][j].dNdEta_err+fome*ETome[i][j].dNdEta_err*ETome[i][j].dNdEta_err;

      expSysUncertaintydETdEta = sqrt(expVariancedETdEta);
      expSysUncertaintydNdEta = sqrt(expVariancedNdEta);

      expSysUncertaintyCorrelateddETdEta = fpi*(ETpip[i][j].dETdEta_err+ETpim[i][j].dETdEta_err)+ fp*(ETp[i][j].dETdEta_err+ETap[i][j].dETdEta_err)+fK*(ETKm[i][j].dETdEta_err+ETKp[i][j].dETdEta_err)+fLam*(ETLab[i][j].dETdEta_err+ETLa[i][j].dETdEta_err)+feta*ETeta[i][j].dETdEta_err+fome*ETome[i][j].dETdEta_err;
      expSysUncertaintyCorrelateddNdEta = fpi*(ETpip[i][j].dNdEta_err+ETpim[i][j].dNdEta_err)+ fp*(ETp[i][j].dNdEta_err+ETap[i][j].dNdEta_err)+fK*(ETKm[i][j].dNdEta_err+ETKp[i][j].dNdEta_err)+fLam*(ETLab[i][j].dNdEta_err+ETLa[i][j].dNdEta_err)+feta*ETeta[i][j].dNdEta_err+fome*ETome[i][j].dNdEta_err;
      /*
      //3.  Extrapolation uncertainties
      extrapVariancedETdEta = fpi*(ETpip[i][j].dETdEta_exterr*ETpip[i][j].dETdEta_exterr+ETpim[i][j].dETdEta_exterr*ETpim[i][j].dETdEta_exterr)+ fp*(ETp[i][j].dETdEta_exterr*ETp[i][j].dETdEta_exterr+ETap[i][j].dETdEta_exterr*ETap[i][j].dETdEta_exterr)+fK*(ETKp[i][j].dETdEta_exterr*ETKp[i][j].dETdEta_exterr+ETKm[i][j].dETdEta_exterr*ETKm[i][j].dETdEta_exterr)+fLam*(ETLa[i][j].dETdEta_exterr*ETLa[i][j].dETdEta_exterr+ETLab[i][j].dETdEta_exterr*ETLab[i][j].dETdEta_exterr)+feta*(ETpip[i][j].dETdEta_exterr*ETpip[i][j].dETdEta_exterr+ETpim[i][j].dETdEta_exterr*ETpim[i][j].dETdEta_exterr)+fome*(ETpip[i][j].dETdEta_exterr*ETpip[i][j].dETdEta_exterr+ETpim[i][j].dETdEta_exterr*ETpim[i][j].dETdEta_exterr);
      extrapVariancedNdEta = fpi*(ETpip[i][j].dNdEta_exterr*ETpip[i][j].dNdEta_exterr+ETpim[i][j].dNdEta_exterr*ETpim[i][j].dNdEta_exterr)+ fp*(ETp[i][j].dNdEta_exterr*ETp[i][j].dNdEta_exterr+ETap[i][j].dNdEta_exterr*ETap[i][j].dNdEta_exterr)+fK*(ETKp[i][j].dNdEta_exterr*ETKp[i][j].dNdEta_exterr+ETKm[i][j].dNdEta_exterr*ETKm[i][j].dNdEta_exterr)+fLam*(ETLa[i][j].dNdEta_exterr*ETLa[i][j].dNdEta_exterr+ETLab[i][j].dNdEta_exterr*ETLab[i][j].dNdEta_exterr)+feta*(ETpip[i][j].dNdEta_exterr*ETpip[i][j].dNdEta_exterr+ETpim[i][j].dNdEta_exterr*ETpim[i][j].dNdEta_exterr)+fome*(ETpip[i][j].dNdEta_exterr*ETpip[i][j].dNdEta_exterr+ETpim[i][j].dNdEta_exterr*ETpim[i][j].dNdEta_exterr);

      extrapSysUncertaintydETdEta = sqrt(extrapVariancedETdEta);
      extrapSysUncertaintydNdEta = sqrt(extrapVariancedNdEta);

      extrapSysUncertaintyCorrelateddETdEta = fpi*(ETpip[i][j].dETdEta_exterr+ETpim[i][j].dETdEta_exterr)+ fp*(ETp[i][j].dETdEta_exterr+ETap[i][j].dETdEta_exterr)+fK*(ETKm[i][j].dETdEta_exterr+ETKp[i][j].dETdEta_exterr)+fLam*(ETLab[i][j].dETdEta_exterr+ETLa[i][j].dETdEta_exterr)+feta*(ETpip[i][j].dETdEta_exterr+ETpim[i][j].dETdEta_exterr)+fome*(ETpip[i][j].dETdEta_exterr+ETpim[i][j].dETdEta_exterr);
      extrapSysUncertaintyCorrelateddNdEta = fpi*(ETpip[i][j].dNdEta_exterr+ETpim[i][j].dNdEta_exterr)+ fp*(ETp[i][j].dNdEta_exterr+ETap[i][j].dNdEta_exterr)+fK*(ETKm[i][j].dNdEta_exterr+ETKp[i][j].dNdEta_exterr)+fLam*(ETLab[i][j].dNdEta_exterr+ETLa[i][j].dNdEta_exterr)+feta*(ETpip[i][j].dNdEta_exterr+ETpim[i][j].dNdEta_exterr)+fome*(ETpip[i][j].dNdEta_exterr+ETpim[i][j].dNdEta_exterr);
      */

      //out<<SNN[i]<<"\t"<<j<<" u "<<nominalETdETdEta<<"   "<<factorSysUncertaintydETdEta<<"   "<<expSysUncertaintydETdEta<<"   "<<extrapSysUncertainty<<endl;
      out<<SNN[i]<<"\t"<<j<<" "<<nominalETdETdEta<<" "<<factorSysUncertaintydETdEta<<"   "<<expSysUncertaintyCorrelateddETdEta<<"   "<<extrapSysUncertaintyCorrelated<<endl;
      //out2<<SNN[i]<<"\t"<<j<<" u "<<nominalETdNdEta<<"   "<<factorSysUncertaintydNdEta<<"   "<<expSysUncertaintydNdEta<<"   "<<extrapSysUncertainty<<endl;
      out2<<SNN[i]<<"\t"<<j<<" "<<nominalETdNdEta<<" "<<factorSysUncertaintydNdEta<<"   "<<expSysUncertaintyCorrelateddNdEta<<"   "<<extrapSysUncertaintyCorrelated<<endl;


      //eta and omega fractional errors
      float etafrac = 2.0*ETeta[i][j].dETdEta/(ETpip[i][j].dETdEta+ETpim[i][j].dETdEta);
      float omegafrac = 2.0*ETome[i][j].dETdEta/(ETpip[i][j].dETdEta+ETpim[i][j].dETdEta);
      //the uncertainties faily in the current implementation because they are treated as point to point uncorrelated, which means that the uncertainty is massively overestimated.  70% uncertainty on eta.  35% on omega.
      float etafracerr = etafrac*0.7;//sqrt( (ETpip[i][j].dETdEta_err*ETpip[i][j].dETdEta_err+ETpim[i][j].dETdEta_err*ETpim[i][j].dETdEta_err)/((ETpip[i][j].dETdEta+ETpim[i][j].dETdEta)*(ETpip[i][j].dETdEta+ETpim[i][j].dETdEta)) +ETeta[i][j].dETdEta_err*ETeta[i][j].dETdEta_err/ETeta[i][j].dETdEta/ETeta[i][j].dETdEta);
      float omegafracerr = omegafrac*0.35;//sqrt( (ETpip[i][j].dETdEta_err*ETpip[i][j].dETdEta_err+ETpim[i][j].dETdEta_err*ETpim[i][j].dETdEta_err)/((ETpip[i][j].dETdEta+ETpim[i][j].dETdEta)*(ETpip[i][j].dETdEta+ETpim[i][j].dETdEta)) +ETome[i][j].dETdEta_err*ETome[i][j].dETdEta_err/ETome[i][j].dETdEta/ETome[i][j].dETdEta);
//       if(j==0){
// 	cout<<SNN[i]<<"\t"<<j<<" "<<etafrac<<" "<<etafracerr<<" "<<omegafrac<<" "<<omegafracerr<<" ";
// 	cout<<" eta "<<ETeta[i][j].dETdEta<<" omega "<<ETome[i][j].dETdEta<<" pi+ "<<ETpim[i][j].dETdEta<<" pi- "<<ETpip[i][j].dETdEta<<" avg "<<(ETpip[i][j].dETdEta+ETpim[i][j].dETdEta)/2.0<<endl;
//       }
      //quick fix because eta and omega ET fractions cannot be negative
//       float low = etafrac-etafracerr;
//       float high = etafrac+etafracerr;
//       if(low<0){
// 	low=0;
//       }
//       etafrac = (high+low)/2.0;
//       etafracerr = (high-low)/2.0;
//       low = omegafrac-omegafracerr;
//       high = omegafrac+omegafracerr;
//       if(low<0){
// 	//cout<<"This one fails "<<SNN[i]<<"\t"<<j<<endl;
// 	low=0;
//       }
//       omegafrac = (high+low)/2.0;
//       omegafracerr = (high-low)/2.0;

      //if(j==0){
	//cout<<SNN[i]<<"\t"<<j<<" "<<etafrac<<" "<<etafracerr<<" "<<omegafrac<<" "<<omegafracerr<<" ";
	//cout<<" eta "<<ETeta[i][j].dETdEta<<" omega "<<ETome[i][j].dETdEta<<" pi+ "<<ETpim[i][j].dETdEta<<" pi- "<<ETpip[i][j].dETdEta<<" avg "<<(ETpip[i][j].dETdEta+ETpim[i][j].dETdEta)/2.0<<endl;
	//}
      //if(j==0){
	out3<<SNN[i]<<" "<<j<<" "<<etafrac<<" "<<etafracerr<<" "<<omegafrac<<" "<<omegafracerr<<endl;
	//}
	//else{
	//if(SNN[i]>130 && j==1){
	//out3<<SNN[i]<<"\t"<<j<<" "<<etafrac<<" "<<etafracerr<<" "<<omegafrac<<" "<<omegafracerr<<endl;
	//}
	//}
    }
    //out<<endl;
    //out2<<endl;
  }
  out3.close();
  out.close();
  return 0;
}

//double nominalET = fpi*(ETpiplus+ETpiminus)+ fp*(ETproton+ETantiproton)+fK*(ETKminus+ETKplus)+fLam*(ETantilambda+ETlambda)+ETetaomega;
//We are going to separate the uncertainties into
//1.  Factor uncertainties
//2.  Experimental uncertainties
//3.  Extrapolation uncertainties

//1.  Factor uncertainties:
//Here we are going to treat everything as constant except the factors and we're going to add the uncertainties from the factors as if they are uncorrelated with each other
//We are going to add in the eta uncertainty because this is largely dominated by the scaling uncertainties
//double factorVariance = fpiErr*fpiErr*(ETpiplus+ETpiminus)*(ETpiplus+ETpiminus)+ fpErr* fpErr*(ETproton+ETantiproton)*(ETproton+ETantiproton)+fKErr*fKErr*(ETKminus+ETKplus)*(ETKminus+ETKplus)+fLamErr*fLamErr*(ETantilambda+ETlambda)*(ETantilambda+ETlambda)+ETetaomegaSys*ETetaomegaSys;
//double factorSysUncertainty = sqrt(factorVariance);

//2.  Experimental uncertainties
//Assuming totally uncorrelated
//double expVariance =fpi*(ETpiplusSysExp*ETpiplusSysExp+ETpiminusSysExp*ETpiminusSysExp)+ fp*(ETprotonSysExp*ETprotonSysExp+ETantiprotonSysExp*ETantiprotonSysExp)+fK*(ETKminusSysExp*ETKminusSysExp+ETKplusSysExp*ETKplusSysExp)+fLam*(ETantilambdaSysExp*ETantilambdaSysExp+ETlambdaSysExp*ETlambdaSysExp) ;
//double expSysUncertainty = sqrt(expVariance);
//Assuming totally correlated
//double expSysUncertaintyCorrelated = fpi*(ETpiplusSysExp+ETpiminusSysExp)+ fp*(ETprotonSysExp+ETantiprotonSysExp)+fK*(ETKminusSysExp+ETKplusSysExp)+fLam*(ETantilambdaSysExp+ETlambdaSysExp);

//3.  Extrapolation uncertainties
//double extrapVariance =fpi*(ETpiplusSysExtrap*ETpiplusSysExtrap+ETpiminusSysExtrap*ETpiminusSysExtrap)+ fp*(ETprotonSysExtrap*ETprotonSysExtrap+ETantiprotonSysExtrap*ETantiprotonSysExtrap)+fK*(ETKminusSysExtrap*ETKminusSysExtrap+ETKplusSysExtrap*ETKplusSysExtrap)+fLam*(ETantilambdaSysExtrap*ETantilambdaSysExtrap+ETlambdaSysExtrap*ETlambdaSysExtrap) ;
//double extrapSysUncertainty = sqrt(extrapVariance);
//Assuming totally correlated
//double extrapSysUncertaintyCorrelated = fpi*(ETpiplusSysExtrap+ETpiminusSysExtrap)+ fp*(ETprotonSysExtrap+ETantiprotonSysExtrap)+fK*(ETKminusSysExtrap+ETKplusSysExtrap)+fLam*(ETantilambdaSysExtrap+ETlambdaSysExtrap);

//cout<<"Assuming uncorrelated"<<endl;
//cout<<"ET "<<nominalET<<"   "<<factorSysUncertainty<<"   "<<expSysUncertainty<<"   "<<extrapSysUncertainty<<endl;
//cout<<"Assuming Correlated"<<endl;
//cout<<"ET "<<nominalET<<"   "<<factorSysUncertainty<<"   "<<expSysUncertaintyCorrelated<<"   "<<extrapSysUncertaintyCorrelated<<endl;

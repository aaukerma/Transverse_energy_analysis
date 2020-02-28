#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

//value
Double_t ETpiplus =1;
//systematic uncertainty from experimental errors
Double_t ETpiplusSysExp = 0.1;
//systematic uncertainty from experimental errors
Double_t ETpiplusSysExtrap = 0.05;
//value
Double_t ETpiminus =1;
//systematic uncertainty from experimental errors
Double_t ETpiminusSysExp = 0.1;
//systematic uncertainty from experimental errors
Double_t ETpiminusSysExtrap = 0.05;
//value
Double_t ETKplus =1;
//systematic uncertainty from experimental errors
Double_t ETKplusSysExp = 0.1;
//systematic uncertainty from experimental errors
Double_t ETKplusSysExtrap = 0.05;
//value
Double_t ETKminus =1;
//systematic uncertainty from experimental errors
Double_t ETKminusSysExp = 0.1;
//systematic uncertainty from experimental errors
Double_t ETKminusSysExtrap = 0.05;
//value
Double_t ETproton =1;
//systematic uncertainty from experimental errors
Double_t ETprotonSysExp = 0.1;
//systematic uncertainty from experimental errors
Double_t ETprotonSysExtrap = 0.05;
//value
Double_t ETantiproton =1;
//systematic uncertainty from experimental errors
Double_t ETantiprotonSysExp = 0.1;
//systematic uncertainty from experimental errors
Double_t ETantiprotonSysExtrap = 0.05;
//value
Double_t ETlambda =1;
//systematic uncertainty from experimental errors
Double_t ETlambdaSysExp = 0.1;
//systematic uncertainty from experimental errors
Double_t ETlambdaSysExtrap = 0.05;
//value
Double_t ETantilambda =1;
//systematic uncertainty from experimental errors
Double_t ETantilambdaSysExp = 0.1;
//systematic uncertainty from experimental errors
Double_t ETantilambdaSysExtrap = 0.05;
//value
Double_t ETetaomega =.00000001;
//systematic uncertainty - we're just rolling this into one because the scale factor uncertainty is dominant
Double_t ETetaomegaSys = 0.1;
Double_t fpi=1.588;
Double_t fpiErr=0.029;
Double_t fK=1.8; //fk- = 1.5+-0.1        fk+ = 0.5 +- 0.05
Double_t fKErr=0.1;
Double_t fp=1.26; //add in for energy dependance
Double_t fpErr=0.26;
Double_t fLam=1.08;//.4373948
Double_t fLamErr=0.51;//.12

Double_t EtaMass=1;
Double_t OmeMass=1;

struct particle {
  double dETdEta;
  double dETdEta_err;
  double dNdEta;
  double dNdEta_err;
  double npart;
  double npart_err;
};
void EtaMath(Double_t pt, Double_t bin){
  Double_t et = TMath::Sqrt(pt*pt+EtaMass*EtaMass);
  Double_t J = pt/et;
  Double_t tr = 2. *TMath::Pi() *pt;
  //Double_t dETdEta = bincont*tr*J*et*dx;
  return;
}
void OmegaMath(Double_t pt, Double_t bin){

  return;
}

void CalcFinalETs(){
  string GRBG;
  Double_t grbgd;
  ifstream in;
  particle q;
  int a = 9;
  int E  =8;
  vector < vector< particle > > ETpip (E, vector<particle>(a));
  vector < vector< particle > > ETpim (E,vector<particle>(a));
  vector < vector< particle > > ETKp (E,vector<particle>(a));
  vector < vector< particle > > ETKm (E,vector<particle>(a));
  vector < vector< particle > > ETp (E,vector<particle>(a));
  vector < vector< particle > > ETap (E,vector<particle>(a));
  vector < vector< particle > > ETLa (E,vector<particle>(a));
  vector < vector< particle > > ETLab (E,vector<particle>(a));
  vector <float> SNN = {7.7,11.5,19.6,27,39,62.4,130,200};
  int cent;
  double dat1,dat2,dat3,dat4,dat5,dat6;
  in.open("Averaged_PiKp_ET_results.txt");
  for (int i = 0; i<a;i++){
    in>>GRBG;
  }
  for (int i=0;i<E;i++){
    for(int j=0;j<a;j++){
      in>>grbgd>>GRBG>>cent>>dat1>>dat2>>dat3>>dat4>>dat5>>dat6;
      ETpip[i][j].dETdEta=dat1;
      ETpip[i][j].dETdEta_err=dat2;
      ETpip[i][j].dNdEta=dat3;
      ETpip[i][j].dNdEta_err=dat4;
      ETpip[i][j].npart=dat5;
      ETpip[i][j].npart_err=dat6;
    }
    for(int j=0;j<a;j++){
      in>>grbgd>>GRBG>>cent>>dat1>>dat2>>dat3>>dat4>>dat5>>dat6;
      ETpim[i][j].dETdEta=dat1;
      ETpim[i][j].dETdEta_err=dat2;
      ETpim[i][j].dNdEta=dat3;
      ETpim[i][j].dNdEta_err=dat4;
      ETpim[i][j].npart=dat5;
      ETpim[i][j].npart_err=dat6;
    }
    for(int j=0;j<a;j++){
      in>>grbgd>>GRBG>>cent>>dat1>>dat2>>dat3>>dat4>>dat5>>dat6;
      ETKp[i][j].dETdEta=dat1;
      ETKp[i][j].dETdEta_err=dat2;
      ETKp[i][j].dNdEta=dat3;
      ETKp[i][j].dNdEta_err=dat4;
      ETKp[i][j].npart=dat5;
      ETKp[i][j].npart_err=dat6;
    }
    for(int j=0;j<a;j++){
      in>>grbgd>>GRBG>>cent>>dat1>>dat2>>dat3>>dat4>>dat5>>dat6;
      ETKm[i][j].dETdEta=dat1;
      ETKm[i][j].dETdEta_err=dat2;
      ETKm[i][j].dNdEta=dat3;
      ETKm[i][j].dNdEta_err=dat4;
      ETKm[i][j].npart=dat5;
      ETKm[i][j].npart_err=dat6;
    }
    for(int j=0;j<a;j++){
      in>>grbgd>>GRBG>>cent>>dat1>>dat2>>dat3>>dat4>>dat5>>dat6;
      ETp[i][j].dETdEta=dat1;
      ETp[i][j].dETdEta_err=dat2;
      ETp[i][j].dNdEta=dat3;
      ETp[i][j].dNdEta_err=dat4;
      ETp[i][j].npart=dat5;
      ETp[i][j].npart_err=dat6;
    }
    for(int j=0;j<a;j++){
      in>>grbgd>>GRBG>>cent>>dat1>>dat2>>dat3>>dat4>>dat5>>dat6;
      ETap[i][j].dETdEta=dat1;
      ETap[i][j].dETdEta_err=dat2;
      ETap[i][j].dNdEta=dat3;
      ETap[i][j].dNdEta_err=dat4;
      ETap[i][j].npart=dat5;
      ETap[i][j].npart_err=dat6;
    }
    for(int j=0;j<a;j++){
      in>>grbgd>>GRBG>>grbgd>>grbgd>>grbgd>>grbgd>>grbgd>>grbgd>>grbgd;
    }
    for(int j=0;j<a;j++){
      in>>grbgd>>GRBG>>grbgd>>grbgd>>grbgd>>grbgd>>grbgd>>grbgd>>grbgd;
    }
  }
  in.close();
  in.open("Averaged_Lam_ET_resutls.txt");
  for (int i = 0; i<a;i++){
    in>>GRBG;
  }
  a=7;
  for (int i=0;i<E;i++){
    if(i>5) a=5;
    for(int j=0;j<a;j++){
      in>>grbgd>>GRBG>>cent>>dat1>>dat2>>dat3>>dat4>>dat5>>dat6;
      ETLa[i][j].dETdEta=dat1;
      ETLa[i][j].dETdEta_err=dat2;
      ETLa[i][j].dNdEta=dat3;
      ETLa[i][j].dNdEta_err=dat4;
      ETLa[i][j].npart=dat5;
      ETLa[i][j].npart_err=dat6;
    }
    for(int j=0;j<a;j++){
      in>>grbgd>>GRBG>>cent>>dat1>>dat2>>dat3>>dat4>>dat5>>dat6;
      ETLab[i][j].dETdEta=dat1;
      ETLab[i][j].dETdEta_err=dat2;
      ETLab[i][j].dNdEta=dat3;
      ETLab[i][j].dNdEta_err=dat4;
      ETLab[i][j].npart=dat5;
      ETLab[i][j].npart_err=dat6;
    }
  }
  in.close();
  ofstream out;
  out.open("EtResults.txt");
  ofstream out2;
  out2.open("EtResultsN.txt");
  ofstream out3;
  out3.open("EtResultsnpart.txt");

  Double_t nominalETdETdEta;
  Double_t nominalETdNdEta;
  Double_t nominalETnpart;

  Double_t factorVariancedETdEta;
  Double_t factorVariancedNdEta;
  Double_t factorVariancenpart;

  Double_t factorSysUncertaintydETdEta;
  Double_t factorSysUncertaintydNdEta;
  Double_t factorSysUncertaintynpart;

  Double_t expVariancedETdEta;
  Double_t expVariancedNdEta;
  Double_t expVariancenpart;

  Double_t expSysUncertaintydETdEta;
  Double_t expSysUncertaintydNdEta;
  Double_t expSysUncertaintynpart;

  Double_t expSysUncertaintyCorrelateddETdEta;
  Double_t expSysUncertaintyCorrelateddNdEta;
  Double_t expSysUncertaintyCorrelatednpart;

  Double_t extrapVariance =fpi*(ETpiplusSysExtrap*ETpiplusSysExtrap+ETpiminusSysExtrap*ETpiminusSysExtrap)+fp*(ETprotonSysExtrap*ETprotonSysExtrap+ETantiprotonSysExtrap*ETantiprotonSysExtrap)+fK*(ETKminusSysExtrap*ETKminusSysExtrap+ETKplusSysExtrap*ETKplusSysExtrap)+fLam*(ETantilambdaSysExtrap*ETantilambdaSysExtrap+ETlambdaSysExtrap*ETlambdaSysExtrap) ;
  Double_t extrapSysUncertainty = TMath::Sqrt(extrapVariance);
  Double_t extrapSysUncertaintyCorrelated = fpi*(ETpiplusSysExtrap+ETpiminusSysExtrap)+fp*(ETprotonSysExtrap+ETantiprotonSysExtrap)+fK*(ETKminusSysExtrap+ETKplusSysExtrap)+fLam*(ETantilambdaSysExtrap+ETlambdaSysExtrap);

  out<<"dETdEta\n";
  out<<" SNN, bin, corr, nominalET +/- factorSysUncertainty +/- expSysUncertainty +/- extrapSysUncertainty\n";
  out2<<"dNdEta\n";
  out2<<" SNN, bin, corr, nominalET +/- factorSysUncertainty +/- expSysUncertainty +/- extrapSysUncertainty\n";
  out3<<"npart\n";
  out3<<" SNN, bin, corr, nominalET +/- factorSysUncertainty +/- expSysUncertainty +/- extrapSysUncertainty\n";
  for(int i=0;i<8;i++){
    for(int j=0;j<9;j++){
      nominalETdETdEta=fpi*(ETpip[i][j].dETdEta+ETpim[i][j].dETdEta)+fp*(ETp[i][j].dETdEta+ETap[i][j].dETdEta)+fK*(ETKm[i][j].dETdEta+ETKp[i][j].dETdEta)+fLam*(ETLab[i][j].dETdEta+ETLa[i][j].dETdEta)+ETetaomega;
      nominalETdNdEta=fpi*(ETpip[i][j].dNdEta+ETpim[i][j].dNdEta)+fp*(ETp[i][j].dNdEta+ETap[i][j].dNdEta)+fK*(ETKm[i][j].dNdEta+ETKp[i][j].dNdEta)+fLam*(ETLab[i][j].dNdEta+ETLa[i][j].dNdEta)+ETetaomega;
      nominalETnpart=fpi*(ETpip[i][j].npart+ETpim[i][j].npart)+fp*(ETp[i][j].npart+ETap[i][j].npart)+fK*(ETKm[i][j].npart+ETKp[i][j].npart)+fLam*(ETLab[i][j].npart+ETLa[i][j].npart)+ETetaomega;

      //1.  Factor uncertainties:
        //Here we are going to treat everything as constant except the factors and we're going to add the uncertainties from the factors as if they are uncorrelated with each other
        //We are going to add in the eta uncertainty because this is largely dominated by the scaling uncertainties
      factorVariancedETdEta = fpiErr*fpiErr*(ETpip[i][j].dETdEta+ETpim[i][j].dETdEta)*(ETpip[i][j].dETdEta+ETpim[i][j].dETdEta)+fpErr*fpErr*(ETp[i][j].dETdEta+ETap[i][j].dETdEta)*(ETp[i][j].dETdEta+ETap[i][j].dETdEta)+fKErr*fKErr*(ETKm[i][j].dETdEta+ETKp[i][j].dETdEta)*(ETKm[i][j].dETdEta+ETKp[i][j].dETdEta)+fLamErr*fLamErr*(ETLab[i][j].dETdEta+ETLa[i][j].dETdEta)*(ETLab[i][j].dETdEta+ETLa[i][j].dETdEta)+ETetaomegaSys*ETetaomegaSys;
      factorVariancedNdEta = fpiErr*fpiErr*(ETpip[i][j].dNdEta+ETpim[i][j].dNdEta)*(ETpip[i][j].dNdEta+ETpim[i][j].dNdEta)+fpErr*fpErr*(ETp[i][j].dNdEta+ETap[i][j].dNdEta)*(ETp[i][j].dNdEta+ETap[i][j].dNdEta)+fKErr*fKErr*(ETKm[i][j].dNdEta+ETKp[i][j].dNdEta)*(ETKm[i][j].dNdEta+ETKp[i][j].dNdEta)+fLamErr*fLamErr*(ETLab[i][j].dNdEta+ETLa[i][j].dNdEta)*(ETLab[i][j].dNdEta+ETLa[i][j].dNdEta)+ETetaomegaSys*ETetaomegaSys;
      factorVariancenpart = fpiErr*fpiErr*(ETpip[i][j].npart+ETpim[i][j].npart)*(ETpip[i][j].npart+ETpim[i][j].npart)+fpErr*fpErr*(ETp[i][j].npart+ETap[i][j].npart)*(ETp[i][j].npart+ETap[i][j].npart)+fKErr*fKErr*(ETKm[i][j].npart+ETKp[i][j].npart)*(ETKm[i][j].npart+ETKp[i][j].npart)+fLamErr*fLamErr*(ETLab[i][j].npart+ETLa[i][j].npart)*(ETLab[i][j].npart+ETLa[i][j].npart)+ETetaomegaSys*ETetaomegaSys;

      factorSysUncertaintydETdEta = TMath::Sqrt(factorVariancedETdEta);
      factorSysUncertaintydNdEta = TMath::Sqrt(factorVariancedNdEta);
      factorSysUncertaintynpart = TMath::Sqrt(factorVariancenpart);

      //2.  Experimental uncertainties
      expVariancedETdEta = fpi*(ETpip[i][j].dETdEta_err*ETpip[i][j].dETdEta_err+ETpim[i][j].dETdEta_err*ETpim[i][j].dETdEta_err)+fp*(ETp[i][j].dETdEta_err*ETp[i][j].dETdEta_err+ETap[i][j].dETdEta_err*ETap[i][j].dETdEta_err)+fK*(ETKp[i][j].dETdEta_err*ETKp[i][j].dETdEta_err+ETKm[i][j].dETdEta_err*ETKm[i][j].dETdEta_err)+fLam*(ETLa[i][j].dETdEta_err*ETLa[i][j].dETdEta_err+ETLab[i][j].dETdEta_err*ETLab[i][j].dETdEta_err);
      expVariancedNdEta = fpi*(ETpip[i][j].dNdEta_err*ETpip[i][j].dNdEta_err+ETpim[i][j].dNdEta_err*ETpim[i][j].dNdEta_err)+fp*(ETp[i][j].dNdEta_err*ETp[i][j].dNdEta_err+ETap[i][j].dNdEta_err*ETap[i][j].dNdEta_err)+fK*(ETKp[i][j].dNdEta_err*ETKp[i][j].dNdEta_err+ETKm[i][j].dNdEta_err*ETKm[i][j].dNdEta_err)+fLam*(ETLa[i][j].dNdEta_err*ETLa[i][j].dNdEta_err+ETLab[i][j].dNdEta_err*ETLab[i][j].dNdEta_err);
      expVariancenpart = fpi*(ETpip[i][j].npart_err*ETpip[i][j].npart_err+ETpim[i][j].npart_err*ETpim[i][j].npart_err)+fp*(ETp[i][j].npart_err*ETp[i][j].npart_err+ETap[i][j].npart_err*ETap[i][j].npart_err)+fK*(ETKp[i][j].npart_err*ETKp[i][j].npart_err+ETKm[i][j].npart_err*ETKm[i][j].npart_err)+fLam*(ETLa[i][j].npart_err*ETLa[i][j].npart_err+ETLab[i][j].npart_err*ETLab[i][j].npart_err);

      expSysUncertaintydETdEta = TMath::Sqrt(expVariancedETdEta);
      expSysUncertaintydNdEta = TMath::Sqrt(expVariancedNdEta);
      expSysUncertaintynpart = TMath::Sqrt(expVariancenpart);

      expSysUncertaintyCorrelateddETdEta = fpi*(ETpip[i][j].dETdEta_err+ETpim[i][j].dETdEta_err)+fp*(ETp[i][j].dETdEta_err+ETap[i][j].dETdEta_err)+fK*(ETKm[i][j].dETdEta_err+ETKp[i][j].dETdEta_err)+fLam*(ETLab[i][j].dETdEta_err+ETLa[i][j].dETdEta_err);
      expSysUncertaintyCorrelateddNdEta = fpi*(ETpip[i][j].dNdEta_err+ETpim[i][j].dNdEta_err)+fp*(ETp[i][j].dNdEta_err+ETap[i][j].dNdEta_err)+fK*(ETKm[i][j].dNdEta_err+ETKp[i][j].dNdEta_err)+fLam*(ETLab[i][j].dNdEta_err+ETLa[i][j].dNdEta_err);
      expSysUncertaintyCorrelatednpart = fpi*(ETpip[i][j].npart_err+ETpim[i][j].npart_err)+fp*(ETp[i][j].npart_err+ETap[i][j].npart_err)+fK*(ETKm[i][j].npart_err+ETKp[i][j].npart_err)+fLam*(ETLab[i][j].npart_err+ETLa[i][j].npart_err);

      out<<SNN[i]<<"\t"<<j<<" u "<<nominalETdETdEta<<" +/- "<<factorSysUncertaintydETdEta<<" +/- "<<expSysUncertaintydETdEta<<" +/- "<<extrapSysUncertainty<<endl;
      out<<SNN[i]<<"\t"<<j<<" c "<<nominalETdETdEta<<" +/- "<<factorSysUncertaintydETdEta<<" +/- "<<expSysUncertaintyCorrelateddETdEta<<" +/- "<<extrapSysUncertaintyCorrelated<<endl;
      out2<<SNN[i]<<"\t"<<j<<" u "<<nominalETdNdEta<<" +/- "<<factorSysUncertaintydNdEta<<" +/- "<<expSysUncertaintydNdEta<<" +/- "<<extrapSysUncertainty<<endl;
      out2<<SNN[i]<<"\t"<<j<<" c "<<nominalETdNdEta<<" +/- "<<factorSysUncertaintydNdEta<<" +/- "<<expSysUncertaintyCorrelateddNdEta<<" +/- "<<extrapSysUncertaintyCorrelated<<endl;
      out3<<SNN[i]<<"\t"<<j<<" u "<<nominalETnpart<<" +/- "<<factorSysUncertaintynpart<<" +/- "<<expSysUncertaintynpart<<" +/- "<<extrapSysUncertainty<<endl;
      out3<<SNN[i]<<"\t"<<j<<" c "<<nominalETnpart<<" +/- "<<factorSysUncertaintynpart<<" +/- "<<expSysUncertaintyCorrelatednpart<<" +/- "<<extrapSysUncertaintyCorrelated<<endl;
    }
    out<<endl;
    out2<<endl;
    out3<<endl;
  }
  //Double_t nominalET = fpi*(ETpiplus+ETpiminus)+fp*(ETproton+ETantiproton)+fK*(ETKminus+ETKplus)+fLam*(ETantilambda+ETlambda)+ETetaomega;
  //We are going to separate the uncertainties into
  //1.  Factor uncertainties
  //2.  Experimental uncertainties
  //3.  Extrapolation uncertainties

  //1.  Factor uncertainties:
  //Here we are going to treat everything as constant except the factors and we're going to add the uncertainties from the factors as if they are uncorrelated with each other
  //We are going to add in the eta uncertainty because this is largely dominated by the scaling uncertainties
  //Double_t factorVariance = fpiErr*fpiErr*(ETpiplus+ETpiminus)*(ETpiplus+ETpiminus)+fpErr*fpErr*(ETproton+ETantiproton)*(ETproton+ETantiproton)+fKErr*fKErr*(ETKminus+ETKplus)*(ETKminus+ETKplus)+fLamErr*fLamErr*(ETantilambda+ETlambda)*(ETantilambda+ETlambda)+ETetaomegaSys*ETetaomegaSys;
  //Double_t factorSysUncertainty = TMath::Sqrt(factorVariance);

  //2.  Experimental uncertainties
  //Assuming totally uncorrelated
  //Double_t expVariance =fpi*(ETpiplusSysExp*ETpiplusSysExp+ETpiminusSysExp*ETpiminusSysExp)+fp*(ETprotonSysExp*ETprotonSysExp+ETantiprotonSysExp*ETantiprotonSysExp)+fK*(ETKminusSysExp*ETKminusSysExp+ETKplusSysExp*ETKplusSysExp)+fLam*(ETantilambdaSysExp*ETantilambdaSysExp+ETlambdaSysExp*ETlambdaSysExp) ;
  //Double_t expSysUncertainty = TMath::Sqrt(expVariance);
  //Assuming totally correlated
  //Double_t expSysUncertaintyCorrelated = fpi*(ETpiplusSysExp+ETpiminusSysExp)+fp*(ETprotonSysExp+ETantiprotonSysExp)+fK*(ETKminusSysExp+ETKplusSysExp)+fLam*(ETantilambdaSysExp+ETlambdaSysExp);

  //3.  Extrapolation uncertainties
  //Double_t extrapVariance =fpi*(ETpiplusSysExtrap*ETpiplusSysExtrap+ETpiminusSysExtrap*ETpiminusSysExtrap)+fp*(ETprotonSysExtrap*ETprotonSysExtrap+ETantiprotonSysExtrap*ETantiprotonSysExtrap)+fK*(ETKminusSysExtrap*ETKminusSysExtrap+ETKplusSysExtrap*ETKplusSysExtrap)+fLam*(ETantilambdaSysExtrap*ETantilambdaSysExtrap+ETlambdaSysExtrap*ETlambdaSysExtrap) ;
  //Double_t extrapSysUncertainty = TMath::Sqrt(extrapVariance);
  //Assuming totally correlated
  //Double_t extrapSysUncertaintyCorrelated = fpi*(ETpiplusSysExtrap+ETpiminusSysExtrap)+fp*(ETprotonSysExtrap+ETantiprotonSysExtrap)+fK*(ETKminusSysExtrap+ETKplusSysExtrap)+fLam*(ETantilambdaSysExtrap+ETlambdaSysExtrap);

  //cout<<"Assuming uncorrelated"<<endl;
  //cout<<"ET "<<nominalET<<" +/- "<<factorSysUncertainty<<" +/- "<<expSysUncertainty<<" +/- "<<extrapSysUncertainty<<endl;
  //cout<<"Assuming Correlated"<<endl;
  //cout<<"ET "<<nominalET<<" +/- "<<factorSysUncertainty<<" +/- "<<expSysUncertaintyCorrelated<<" +/- "<<extrapSysUncertaintyCorrelated<<endl;
  out.close();

}

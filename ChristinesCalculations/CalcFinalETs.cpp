#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;



struct particle {
  double dETdEta;
  double dETdEta_err;
  double dNdEta;
  double dNdEta_err;
  double npart;
  double npart_err;
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
  double fpi=1.588;
  double fpiErr=0.029;
  double fK=1.8; //fk- = 1.5+-0.1        fk+ = 0.5 +- 0.05
  double fKErr=0.1;
  double fp=1.26; //add in for energy dependance
  double fpErr=0.26;
  double fLam=1.08;//.4373948
  double fLamErr=0.51;//.12

  double EtaMass=1;
  double OmeMass=1;

  string GRBG;
  double grbgd;
  ifstream in;
  particle q;
  int a = 9;
  int E  =8;
  cout<<"I'm Running!\n";
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

  double nominalETdETdEta;
  double nominalETdNdEta;
  double nominalETnpart;

  double factorVariancedETdEta;
  double factorVariancedNdEta;
  double factorVariancenpart;

  double factorSysUncertaintydETdEta;
  double factorSysUncertaintydNdEta;
  double factorSysUncertaintynpart;

  double expVariancedETdEta;
  double expVariancedNdEta;
  double expVariancenpart;

  double expSysUncertaintydETdEta;
  double expSysUncertaintydNdEta;
  double expSysUncertaintynpart;

  double expSysUncertaintyCorrelateddETdEta;
  double expSysUncertaintyCorrelateddNdEta;
  double expSysUncertaintyCorrelatednpart;

  double extrapVariance =fpi*(ETpiplusSysExtrap*ETpiplusSysExtrap+ETpiminusSysExtrap*ETpiminusSysExtrap)+fp*(ETprotonSysExtrap*ETprotonSysExtrap+ETantiprotonSysExtrap*ETantiprotonSysExtrap)+fK*(ETKminusSysExtrap*ETKminusSysExtrap+ETKplusSysExtrap*ETKplusSysExtrap)+fLam*(ETantilambdaSysExtrap*ETantilambdaSysExtrap+ETlambdaSysExtrap*ETlambdaSysExtrap) ;
  double extrapSysUncertainty = sqrt(extrapVariance);
  double extrapSysUncertaintyCorrelated = fpi*(ETpiplusSysExtrap+ETpiminusSysExtrap)+fp*(ETprotonSysExtrap+ETantiprotonSysExtrap)+fK*(ETKminusSysExtrap+ETKplusSysExtrap)+fLam*(ETantilambdaSysExtrap+ETlambdaSysExtrap);

  out<<"dETdEta\n";
  out<<" SNN, bin, corr, nominalET   factorSysUncertainty   expSysUncertainty   extrapSysUncertainty\n";
  out2<<"dNdEta\n";
  out2<<" SNN, bin, corr, nominalET   factorSysUncertainty   expSysUncertainty   extrapSysUncertainty\n";
  out3<<"npart\n";
  out3<<" SNN, bin, corr, nominalET   factorSysUncertainty   expSysUncertainty   extrapSysUncertainty\n";
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

      factorSysUncertaintydETdEta = sqrt(factorVariancedETdEta);
      factorSysUncertaintydNdEta = sqrt(factorVariancedNdEta);
      factorSysUncertaintynpart = sqrt(factorVariancenpart);

      //2.  Experimental uncertainties
      expVariancedETdEta = fpi*(ETpip[i][j].dETdEta_err*ETpip[i][j].dETdEta_err+ETpim[i][j].dETdEta_err*ETpim[i][j].dETdEta_err)+fp*(ETp[i][j].dETdEta_err*ETp[i][j].dETdEta_err+ETap[i][j].dETdEta_err*ETap[i][j].dETdEta_err)+fK*(ETKp[i][j].dETdEta_err*ETKp[i][j].dETdEta_err+ETKm[i][j].dETdEta_err*ETKm[i][j].dETdEta_err)+fLam*(ETLa[i][j].dETdEta_err*ETLa[i][j].dETdEta_err+ETLab[i][j].dETdEta_err*ETLab[i][j].dETdEta_err);
      expVariancedNdEta = fpi*(ETpip[i][j].dNdEta_err*ETpip[i][j].dNdEta_err+ETpim[i][j].dNdEta_err*ETpim[i][j].dNdEta_err)+fp*(ETp[i][j].dNdEta_err*ETp[i][j].dNdEta_err+ETap[i][j].dNdEta_err*ETap[i][j].dNdEta_err)+fK*(ETKp[i][j].dNdEta_err*ETKp[i][j].dNdEta_err+ETKm[i][j].dNdEta_err*ETKm[i][j].dNdEta_err)+fLam*(ETLa[i][j].dNdEta_err*ETLa[i][j].dNdEta_err+ETLab[i][j].dNdEta_err*ETLab[i][j].dNdEta_err);
      expVariancenpart = fpi*(ETpip[i][j].npart_err*ETpip[i][j].npart_err+ETpim[i][j].npart_err*ETpim[i][j].npart_err)+fp*(ETp[i][j].npart_err*ETp[i][j].npart_err+ETap[i][j].npart_err*ETap[i][j].npart_err)+fK*(ETKp[i][j].npart_err*ETKp[i][j].npart_err+ETKm[i][j].npart_err*ETKm[i][j].npart_err)+fLam*(ETLa[i][j].npart_err*ETLa[i][j].npart_err+ETLab[i][j].npart_err*ETLab[i][j].npart_err);

      expSysUncertaintydETdEta = sqrt(expVariancedETdEta);
      expSysUncertaintydNdEta = sqrt(expVariancedNdEta);
      expSysUncertaintynpart = sqrt(expVariancenpart);

      expSysUncertaintyCorrelateddETdEta = fpi*(ETpip[i][j].dETdEta_err+ETpim[i][j].dETdEta_err)+fp*(ETp[i][j].dETdEta_err+ETap[i][j].dETdEta_err)+fK*(ETKm[i][j].dETdEta_err+ETKp[i][j].dETdEta_err)+fLam*(ETLab[i][j].dETdEta_err+ETLa[i][j].dETdEta_err);
      expSysUncertaintyCorrelateddNdEta = fpi*(ETpip[i][j].dNdEta_err+ETpim[i][j].dNdEta_err)+fp*(ETp[i][j].dNdEta_err+ETap[i][j].dNdEta_err)+fK*(ETKm[i][j].dNdEta_err+ETKp[i][j].dNdEta_err)+fLam*(ETLab[i][j].dNdEta_err+ETLa[i][j].dNdEta_err);
      expSysUncertaintyCorrelatednpart = fpi*(ETpip[i][j].npart_err+ETpim[i][j].npart_err)+fp*(ETp[i][j].npart_err+ETap[i][j].npart_err)+fK*(ETKm[i][j].npart_err+ETKp[i][j].npart_err)+fLam*(ETLab[i][j].npart_err+ETLa[i][j].npart_err);

      //out<<SNN[i]<<"\t"<<j<<" u "<<nominalETdETdEta<<"   "<<factorSysUncertaintydETdEta<<"   "<<expSysUncertaintydETdEta<<"   "<<extrapSysUncertainty<<endl;
      out<<SNN[i]<<"\t"<<j<<" "<<nominalETdETdEta<<" "<<factorSysUncertaintydETdEta<<"   "<<expSysUncertaintyCorrelateddETdEta<<"   "<<extrapSysUncertaintyCorrelated<<endl;
      //out2<<SNN[i]<<"\t"<<j<<" u "<<nominalETdNdEta<<"   "<<factorSysUncertaintydNdEta<<"   "<<expSysUncertaintydNdEta<<"   "<<extrapSysUncertainty<<endl;
      out2<<SNN[i]<<"\t"<<j<<" "<<nominalETdNdEta<<" "<<factorSysUncertaintydNdEta<<"   "<<expSysUncertaintyCorrelateddNdEta<<"   "<<extrapSysUncertaintyCorrelated<<endl;
      //out3<<SNN[i]<<"\t"<<j<<" u "<<nominalETnpart<<"   "<<factorSysUncertaintynpart<<"   "<<expSysUncertaintynpart<<"   "<<extrapSysUncertainty<<endl;
      out3<<SNN[i]<<"\t"<<j<<" "<<nominalETnpart<<"   "<<factorSysUncertaintynpart<<"   "<<expSysUncertaintyCorrelatednpart<<"   "<<extrapSysUncertaintyCorrelated<<endl;
    }
    out<<endl;
    out2<<endl;
    out3<<endl;
  }
  //double nominalET = fpi*(ETpiplus+ETpiminus)+fp*(ETproton+ETantiproton)+fK*(ETKminus+ETKplus)+fLam*(ETantilambda+ETlambda)+ETetaomega;
  //We are going to separate the uncertainties into
  //1.  Factor uncertainties
  //2.  Experimental uncertainties
  //3.  Extrapolation uncertainties

  //1.  Factor uncertainties:
  //Here we are going to treat everything as constant except the factors and we're going to add the uncertainties from the factors as if they are uncorrelated with each other
  //We are going to add in the eta uncertainty because this is largely dominated by the scaling uncertainties
  //double factorVariance = fpiErr*fpiErr*(ETpiplus+ETpiminus)*(ETpiplus+ETpiminus)+fpErr*fpErr*(ETproton+ETantiproton)*(ETproton+ETantiproton)+fKErr*fKErr*(ETKminus+ETKplus)*(ETKminus+ETKplus)+fLamErr*fLamErr*(ETantilambda+ETlambda)*(ETantilambda+ETlambda)+ETetaomegaSys*ETetaomegaSys;
  //double factorSysUncertainty = sqrt(factorVariance);

  //2.  Experimental uncertainties
  //Assuming totally uncorrelated
  //double expVariance =fpi*(ETpiplusSysExp*ETpiplusSysExp+ETpiminusSysExp*ETpiminusSysExp)+fp*(ETprotonSysExp*ETprotonSysExp+ETantiprotonSysExp*ETantiprotonSysExp)+fK*(ETKminusSysExp*ETKminusSysExp+ETKplusSysExp*ETKplusSysExp)+fLam*(ETantilambdaSysExp*ETantilambdaSysExp+ETlambdaSysExp*ETlambdaSysExp) ;
  //double expSysUncertainty = sqrt(expVariance);
  //Assuming totally correlated
  //double expSysUncertaintyCorrelated = fpi*(ETpiplusSysExp+ETpiminusSysExp)+fp*(ETprotonSysExp+ETantiprotonSysExp)+fK*(ETKminusSysExp+ETKplusSysExp)+fLam*(ETantilambdaSysExp+ETlambdaSysExp);

  //3.  Extrapolation uncertainties
  //double extrapVariance =fpi*(ETpiplusSysExtrap*ETpiplusSysExtrap+ETpiminusSysExtrap*ETpiminusSysExtrap)+fp*(ETprotonSysExtrap*ETprotonSysExtrap+ETantiprotonSysExtrap*ETantiprotonSysExtrap)+fK*(ETKminusSysExtrap*ETKminusSysExtrap+ETKplusSysExtrap*ETKplusSysExtrap)+fLam*(ETantilambdaSysExtrap*ETantilambdaSysExtrap+ETlambdaSysExtrap*ETlambdaSysExtrap) ;
  //double extrapSysUncertainty = sqrt(extrapVariance);
  //Assuming totally correlated
  //double extrapSysUncertaintyCorrelated = fpi*(ETpiplusSysExtrap+ETpiminusSysExtrap)+fp*(ETprotonSysExtrap+ETantiprotonSysExtrap)+fK*(ETKminusSysExtrap+ETKplusSysExtrap)+fLam*(ETantilambdaSysExtrap+ETlambdaSysExtrap);

  //cout<<"Assuming uncorrelated"<<endl;
  //cout<<"ET "<<nominalET<<"   "<<factorSysUncertainty<<"   "<<expSysUncertainty<<"   "<<extrapSysUncertainty<<endl;
  //cout<<"Assuming Correlated"<<endl;
  //cout<<"ET "<<nominalET<<"   "<<factorSysUncertainty<<"   "<<expSysUncertaintyCorrelated<<"   "<<extrapSysUncertaintyCorrelated<<endl;
  out.close();
  return 0;
}

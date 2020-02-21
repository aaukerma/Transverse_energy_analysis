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
Double_t ETetaomega =1;
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

void CalcFinalETs(){

  Double_t nominalET = fpi*(ETpiplus+ETpiminus)+fp*(ETproton+ETantiproton)+fK*(ETKminus+ETKplus)+fLam*(ETantilambda+ETlambda)+ETetaomega;
  //We are going to separate the uncertainties into
  //1.  Factor uncertainties
  //2.  Experimental uncertainties
  //3.  Extrapolation uncertainties

  //1.  Factor uncertainties:
  //Here we are going to treat everything as constant except the factors and we're going to add the uncertainties from the factors as if they are uncorrelated with each other
  //We are going to add in the eta uncertainty because this is largely dominated by the scaling uncertainties
  Double_t factorVariance = fpiErr*fpiErr*(ETpiplus+ETpiminus)*(ETpiplus+ETpiminus)+fpErr*fpErr*(ETproton+ETantiproton)*(ETproton+ETantiproton)+fKErr*fKErr*(ETKminus+ETKplus)*(ETKminus+ETKplus)+fLamErr*fLamErr*(ETantilambda+ETlambda)*(ETantilambda+ETlambda)+ETetaomegaSys*ETetaomegaSys;
  Double_t factorSysUncertainty = TMath::Sqrt(factorVariance);

  //2.  Experimental uncertainties
  //Assuming totally uncorrelated
  Double_t expVariance =fpi*(ETpiplusSysExp*ETpiplusSysExp+ETpiminusSysExp*ETpiminusSysExp)+fp*(ETprotonSysExp*ETprotonSysExp+ETantiprotonSysExp*ETantiprotonSysExp)+fK*(ETKminusSysExp*ETKminusSysExp+ETKplusSysExp*ETKplusSysExp)+fLam*(ETantilambdaSysExp*ETantilambdaSysExp+ETlambdaSysExp*ETlambdaSysExp) ;
  Double_t expSysUncertainty = TMath::Sqrt(expVariance);
  //Assuming totally correlated
  Double_t expSysUncertaintyCorrelated = fpi*(ETpiplusSysExp+ETpiminusSysExp)+fp*(ETprotonSysExp+ETantiprotonSysExp)+fK*(ETKminusSysExp+ETKplusSysExp)+fLam*(ETantilambdaSysExp+ETlambdaSysExp);

  //3.  Extrapolation uncertainties
  Double_t extrapVariance =fpi*(ETpiplusSysExtrap*ETpiplusSysExtrap+ETpiminusSysExtrap*ETpiminusSysExtrap)+fp*(ETprotonSysExtrap*ETprotonSysExtrap+ETantiprotonSysExtrap*ETantiprotonSysExtrap)+fK*(ETKminusSysExtrap*ETKminusSysExtrap+ETKplusSysExtrap*ETKplusSysExtrap)+fLam*(ETantilambdaSysExtrap*ETantilambdaSysExtrap+ETlambdaSysExtrap*ETlambdaSysExtrap) ;
  Double_t extrapSysUncertainty = TMath::Sqrt(extrapVariance);
  //Assuming totally correlated
  Double_t extrapSysUncertaintyCorrelated = fpi*(ETpiplusSysExtrap+ETpiminusSysExtrap)+fp*(ETprotonSysExtrap+ETantiprotonSysExtrap)+fK*(ETKminusSysExtrap+ETKplusSysExtrap)+fLam*(ETantilambdaSysExtrap+ETlambdaSysExtrap);

  cout<<"Assuming uncorrelated"<<endl;
  cout<<"ET "<<nominalET<<" +/- "<<factorSysUncertainty<<" +/- "<<expSysUncertainty<<" +/- "<<extrapSysUncertainty<<endl;
  cout<<"Assuming Correlated"<<endl;
  cout<<"ET "<<nominalET<<" +/- "<<factorSysUncertainty<<" +/- "<<expSysUncertaintyCorrelated<<" +/- "<<extrapSysUncertaintyCorrelated<<endl;


}

/****************************
jacobianTest.cpp
Purpose of this program is to test assumtions made in fitBESData5_X.cpp and fitBESData.h
assumption is that eta=0

This code will be inserted into fitBESData5.h for testing in real time
Start:20180629.1330
Complete:20180629.1402

****************************/

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
using namespace std;

int main(){
  double Alert=50; //Threshold of alert
  double eta = 1;
  double pt=0; //average of bins
  double mass = 0.93827; //pi+-
  double J;
  double Ja;
  double comp;

  for (int i=0;i<22;i++){
    pt+=.25;
    J = pt/(sqrt(pt*pt+mass*mass)); //Biswas method
    Ja = (pt*cosh(eta))/(sqrt(pt*pt*cosh(eta)+mass*mass)); //with eta=1
    comp = (Ja-J)/Ja;
    if (comp < 0){
      comp=comp*(-1.0);
    }
    comp=comp*100;
    if (comp>Alert){
      cout<<"JTEST FAIL"<<endl;
      cout<<"Biswas code="<<J<<"|| With eta="<<Ja<<endl;
    }
    cout<<"%DIFF="<<comp<<endl;
  }
  cout<<"done"<<endl;
  return 0;
}

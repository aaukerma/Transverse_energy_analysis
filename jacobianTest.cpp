/****************************
jacobianTest.cpp
Purpose of this program is to test assumtions made in fitBESData5_X.cpp and fitBESData.h
assumption is that eta=0

****************************/

#include <iostream>
#include <string>
#include "TKey.h"
#include <fstream>
#include "fitBESData5.h"
using namespace std;

int jacobianTest(){
  double Alert=10; //Threshold of alert
  double J = pt/(TMath::Sqrt(pt*pt+mass*mass)); //Biswas method
  double Ja = pt/(TMath::Sqrt(pt*pt+mass*mass)); //with eta=1
  double comp = (Ja-J)/Ja;
  if (comp < 0){
    comp=comp*(-1.0);
  }
  comp=comp*100;
  if (comp>Alert){
    cout<<"JTEST FAIL"<<endl;
  }

  return 0;
}

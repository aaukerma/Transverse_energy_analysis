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
#include "TKey.h"
#include <fstream>
#include "fitBESData5.h"
using namespace std;

int jacobianTest(){
  double Alert=10; //Threshold of alert
  double eta = 1;
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

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <vector>
#include "TMath.h"
#include "TH1F.h"
#include "TCanvas.h"
using namespace std;

class particle{
public:
  int num;
  double ET;
  double PT;
  double perET;
  double perPT;
};

int calculateRatio0(){
  int NUM = 8; //CHANGE TO 10 IF USING 2760 and 5020
  vector<int> SNN {7,11,19,27,39,62,130,200,2760,5020};
  char type;
  int date,t;
  double fpi0t=0,fgamt=0,fetat=0,fomet=0,fsigpt=0,fsigmt=0,fsig0t=0,fsigvLamt=0;
  double temp;
  vector <double> fpi0(NUM),fgam(NUM),feta(NUM),fome(NUM),fsigp(NUM),fsigm(NUM),fsig0(NUM),fsigvLam(NUM);
  vector <particle> pip(NUM),pim(NUM),pi0(NUM),eta(NUM),ome(NUM),etaTOT(NUM),omeTOT(NUM),gam(NUM),lam(NUM),lamB(NUM),sigp(NUM),sigm(NUM),sig0(NUM),sigpB(NUM),sigmB(NUM),sig0B(NUM),exo(NUM);
  double NTRASH;
  stringstream filename;
  string fn;
  string TRASH;
  ifstream in;
  printf("Input files to run\nFile Form: char date\n");
  cin>>type>>date;
  for (int i=0;i<NUM;i++){
    filename.flush();
    filename.clear();
    filename.str("");
    filename<<"Data/outfile/p"<<type<<date<<"_"<<SNN[i]<<"GeVoutfile.txt";
    //sprintf(filename,"p%c%d_%dGeVoutfile.txt",type,date,SNN[i]);
    fn=filename.str();
    in.open(fn);
    t=0;
    while(t<5){
      in>>TRASH;
      t++;
    }
    in>>NTRASH>>TRASH>>TRASH>>NTRASH>>TRASH>>NTRASH;
    t=0;
    while(t<11){
      in>>TRASH;
      t++;
    }
    in>>NTRASH;
    t=0;
    while(t<10){
      in>>TRASH;
      t++;
    }
    in>>NTRASH>>TRASH>>TRASH>>TRASH>>NTRASH;
    t=0;
    while(t<6){
      in>>TRASH;
      t++;
    }
    in>>TRASH>>pip[i].num>>pip[i].ET>>pip[i].PT>>pip[i].perET>>pip[i].perPT;
    in>>TRASH>>pim[i].num>>pim[i].ET>>pim[i].PT>>pim[i].perET>>pim[i].perPT;
    in>>TRASH>>pi0[i].num>>pi0[i].ET>>pi0[i].PT>>pi0[i].perET>>pi0[i].perPT;
    t=0;
    while(t<4){
      in>>TRASH>>NTRASH>>NTRASH>>NTRASH>>NTRASH>>NTRASH;
      t++;
    }
    in>>TRASH>>lam[i].num>>lam[i].ET>>lam[i].PT>>lam[i].perET>>lam[i].perPT;
    in>>TRASH>>lamB[i].num>>lamB[i].ET>>lamB[i].PT>>lamB[i].perET>>lamB[i].perPT;
    t=0;
    while(t<4){
      in>>TRASH>>NTRASH>>NTRASH>>NTRASH>>NTRASH>>NTRASH;
      t++;
    }
    in>>TRASH>>etaTOT[i].num>>etaTOT[i].ET>>etaTOT[i].PT>>etaTOT[i].perET>>etaTOT[i].perPT;
    in>>TRASH>>omeTOT[i].num>>omeTOT[i].ET>>omeTOT[i].PT>>omeTOT[i].perET>>omeTOT[i].perPT;
    in>>TRASH>>eta[i].num>>eta[i].ET>>eta[i].PT>>eta[i].perET>>eta[i].perPT;
    in>>TRASH>>ome[i].num>>ome[i].ET>>ome[i].PT>>ome[i].perET>>ome[i].perPT;
    t=0;
    while(t<3){
      in>>TRASH>>NTRASH>>NTRASH>>NTRASH>>NTRASH>>NTRASH;
      t++;
    }
    in>>TRASH>>sigp[i].num>>sigp[i].ET>>sigp[i].PT>>sigp[i].perET>>sigp[i].perPT;
    in>>TRASH>>sigm[i].num>>sigm[i].ET>>sigm[i].PT>>sigm[i].perET>>sigm[i].perPT;
    in>>TRASH>>sig0[i].num>>sig0[i].ET>>sig0[i].PT>>sig0[i].perET>>sig0[i].perPT;
    in>>TRASH>>sigpB[i].num>>sigpB[i].ET>>sigpB[i].PT>>sigpB[i].perET>>sigpB[i].perPT;
    in>>TRASH>>sigmB[i].num>>sigmB[i].ET>>sigmB[i].PT>>sigmB[i].perET>>sigmB[i].perPT;
    in>>TRASH>>sig0B[i].num>>sig0B[i].ET>>sig0B[i].PT>>sig0B[i].perET>>sig0B[i].perPT;
    in>>TRASH>>gam[i].num>>gam[i].ET>>gam[i].PT>>gam[i].perET>>gam[i].perPT;
    t=0;
    while(t<8){
      in>>TRASH>>NTRASH>>NTRASH>>NTRASH>>NTRASH>>NTRASH;
      t++;
    }
    in>>TRASH>>exo[i].num>>exo[i].ET>>exo[i].PT>>exo[i].perET>>exo[i].perPT; //need to changes this to the exotics.
    in.close();
  }
  for (int i=0;i<NUM;i++){
    fpi0[i]=pi0[i].ET/((pip[i].ET+pim[i].ET)/2);
    fgam[i]=(gam[i].ET+exo[i].ET)/(pip[i].ET+pim[i].ET);
    feta[i]=eta[i].ET/etaTOT[i].ET;
    fome[i]=ome[i].ET/omeTOT[i].ET;
    fsigp[i]=sigp[i].ET/(pip[i].ET+pim[i].ET);
    fsigm[i]=sigm[i].ET/(pip[i].ET+pim[i].ET);
    fsig0[i]=sig0[i].ET/(pip[i].ET+pim[i].ET);
    temp=(lam[i].ET+lamB[i].ET+sigp[i].ET+sigm[i].ET+sigpB[i].ET+sigmB[i].ET+sig0[i].ET+sig0B[i].ET)/(lam[i].ET+lamB[i].ET+sig0[i].ET+sig0B[i].ET);
    fsigvLam[i]=temp;
    printf("SNN: %d\tfpi0: %f\tfeta: %f\tfomega: %f\tfexotics+gamma: %f\tfsig: %f\n",SNN[i],fpi0[i],feta[i],fome[i],fgam[i],temp);
    //printf("fpi0= %f\tfgam= %f\tfeta= %f\tfome= %f\n",fpi0,fgam,feta,fome);
  }
  int nobby=0;
  if (type=='A') nobby=1;
  for (int i=nobby;i<NUM;i++){
    fpi0t+=fpi0[i];
    fgamt+=fgam[i];
    fetat+=feta[i];
    fomet+=fome[i];
    fsigpt+=fsigp[i];
    fsigmt+=fsigm[i];
    fsig0t+=fsig0[i];
    fsigvLamt+=fsigvLam[i];
  }
  fpi0t=fpi0t/(NUM-nobby); //
  fgamt=fgamt/(NUM-nobby);
  fetat=fetat/(NUM-nobby);
  fomet=fomet/(NUM-nobby);
  fsigpt=fsigpt/(NUM-nobby);
  fsigmt=fsigmt/(NUM-nobby);
  fsig0t=fsig0t/(NUM-nobby);
  fsigvLamt=fsigvLamt/(NUM-nobby);
  double fsigs = fsigpt+fsigmt+fsig0t;
  printf("fpi0= %f\nf_gamma+exotics= %f\nfeta= %f\nfome= %f\nfsigs= %f\n",fpi0t,fgamt,fetat,fomet,fsigvLamt);

  return 0;
}

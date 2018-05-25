//BESData_InclPi0.cpp
/*This program will take data from BESData_sorted.txt and include
approximated data for Pi0 particle. built from BESDataToRootFile.cpp
start:20180522.11:00
compl:20180525.11:14
*/

#include <cstdio>
#include <vector>
#include <string>
#include <fstream>
#include <iostream> // to use cout for debugging
#include <iomanip>
using namespace std;

struct BinData { //
  double pTl;
  double pTh;
  double pTSpec;
  double ErrStat;
  double ErrSys;
};

struct Bin {
  int BinIndex;
  BinData Dat;
};

void pi0Builder(const vector<vector<Bin>>& vect1, const vector<vector<Bin>>& vect2, vector<vector<Bin>>& Pi0);

int main(){
  string CollType;
  string CollEner;
  string title0 = "------------------------"; //This is manually entered (it could potentially be read in directly)[applies to all "titleX"]
  string title1 = "pTlow [GeV/c]";
  string title2 = "pThigh [GeV/c]";
  string title3 = "d^2N/(2pi*pt*dpt*dy)[(GeV/c)^-2]";
  string title4 = "error(stat.)";
  string title5 = "error(sys.)";
  string title6 = "========================";
  string PartType;
  string CentBin;
  string GARBAGE;
  int c1;
  int c2;
  int j;
  int k;
  int l=0;
  int PiSize;
  double DataPoint;
  vector <Bin> PiM(0);
  vector <Bin> PiP(0);
  vector <vector <Bin>> PiMinus(9);
  vector <vector <Bin>> PiPlus(9);
  vector <vector <Bin>> Pi0(9);

  ifstream in;
	in.open("./BESData_sorted.txt"); //Formatting must match previous versions!!!!
  ofstream out;
  out.open("./BESData_sorted_PlusPi0.txt"); //created file for inclusion of pi0, formatted!

  while(in>>CollType){
    while(CollType!="Au+Au") {in>>CollType;}
    in>>CollEner;
    out<<CollType<<" ";
    out<<CollEner<<" GeV"<<endl;
    for(int i = 0; i<8; i++) {in >> GARBAGE;} //skips up to species
    out <<title0<<endl<<title1<<" "<<title2<<" "<<title3<<" "<<title4<<" "<<title5<<endl;

    for(j = 0; j<9; j++) {
      in>>PartType;
      if (j>=1){ //this is a quick fix, there is a 'title0' that appears and has to be skipped. doesnt work at end of loop for some reason
        in>>PartType;
      }
      else{}
      out<<PartType<<"  ";
      in>>CentBin;
      out<<CentBin<<endl;
      if (PartType == "pi-") {//write in data for vector PiMinus and writes to new file
        k=0;
        while (in>>DataPoint){ //This counts as a (in>>) therefore first bit not in>>
          PiM.push_back(Bin());
          PiM[k].BinIndex = k;
          PiM[k].Dat.pTl = DataPoint;
          out<<fixed<<setprecision(2)<<PiM[k].Dat.pTl<<" \t ";
          in>>PiM[k].Dat.pTh;
          out<<fixed<<setprecision(2)<<PiM[k].Dat.pTh<<" \t ";
          in>>PiM[k].Dat.pTSpec;
          out<<fixed<<setprecision(4)<<PiM[k].Dat.pTSpec<<" \t ";
          in>>PiM[k].Dat.ErrStat;
          out<<fixed<<setprecision(5)<<PiM[k].Dat.ErrStat<<" \t ";
          in>>PiM[k].Dat.ErrSys;
          out<<fixed<<setprecision(5)<<PiM[k].Dat.ErrSys<<endl;
          PiMinus[j].push_back(Bin());
          PiMinus[j][k]=PiM[k];
          k++;
        }
      }
      else if (PartType =="pi+") { //write in data for vector PiPlus and writes to new file
        k=0;
        while (in>>DataPoint){
          PiP.push_back(Bin());
          PiP[k].BinIndex = k;
          PiP[k].Dat.pTl = DataPoint;
          out<<fixed<<setprecision(2)<<PiP[k].Dat.pTl<<" \t ";
          in>>PiP[k].Dat.pTh;
          out<<fixed<<setprecision(2)<<PiP[k].Dat.pTh<<" \t ";
          in>>PiP[k].Dat.pTSpec;
          out<<fixed<<setprecision(4)<<PiP[k].Dat.pTSpec<<" \t ";
          in>>PiP[k].Dat.ErrStat;
          out<<fixed<<setprecision(5)<<PiP[k].Dat.ErrStat<<" \t ";
          in>>PiP[k].Dat.ErrSys;
          out<<fixed<<setprecision(5)<<PiP[k].Dat.ErrSys<<endl;
          PiPlus[j].push_back(Bin());
          PiPlus[j][k]=PiP[k];
          k++;
        }
      }
      else {//writes to new file
        while (in>>DataPoint){
          //in>>DataPoint; this is already done in boolean
          out<<fixed<<setprecision(2)<<DataPoint<<" \t ";
          in>>DataPoint;
          out<<fixed<<setprecision(2)<<DataPoint<<" \t ";
          in>>DataPoint;
          out<<fixed<<setprecision(4)<<DataPoint<<" \t ";
          in>>DataPoint;
          out<<fixed<<setprecision(5)<<DataPoint<<" \t ";
          in>>DataPoint;
          out<<fixed<<setprecision(5)<<DataPoint<<endl;
        }
      }
      out<<title0<<endl;
      in.clear();
    } //end forloop ()
    in>>GARBAGE;
    out<<title6<<endl;
    if (l==1){
      pi0Builder(PiMinus, PiPlus, Pi0); //runs once for each energy group after pi+ is fed
      out<<CollType<<" ";
      out<<CollEner<<" GeV"<<endl;
      out <<title0<<endl<<title1<<" "<<title2<<" "<<title3<<" "<<title4<<" "<<title5<<endl;
      for (int n=0; n<9; n++) { //writes in Pi0
        c1=n*10;
        c2=(n+1)*10;
        out<<"pi0  "<<c1<<"-"<<c2<<'%'<<endl;
        PiSize=Pi0[n].size();
        for (int m=0; m<PiSize;m++) {
          if (Pi0[n][m].Dat.pTl != 0) {
            out<<fixed<<setprecision(2)<<Pi0[n][m].Dat.pTl<<" \t ";
            out<<fixed<<setprecision(2)<<Pi0[n][m].Dat.pTh<<" \t ";
            out<<fixed<<setprecision(4)<<Pi0[n][m].Dat.pTSpec<<" \t ";
            out<<fixed<<setprecision(5)<<Pi0[n][m].Dat.ErrStat<<" \t ";
            out<<fixed<<setprecision(5)<<Pi0[n][m].Dat.ErrSys<<endl;
          }
        }
        out<<title0<<endl;
      }
      out<<title6<<endl;
    }
    l++;
    if (l==6) {l=0;}
    in.clear();
  } //end while (in>>CollType)
cout<<"done";
in.close();
out.close();
return 0;
} //end int 'main'


//pi0Builder passes in vectors consisting of all pi- and pi+ data points.
//pi0Builder takes average of all points and creates new vector with
//all data for pi0 then writes it to file.
void pi0Builder(const vector<vector<Bin>>& vect1, const vector<vector<Bin>>& vect2, vector<vector<Bin>>& Pi0){
  //int j=0; //size counter for vect1
  //int k=0; //size coutner for vect2 (may be different than vect1)
  int size;
  Bin pim;
  Bin pip;
  Bin pizero;
  for (int i=0; i<9; i++){
    /*if (vect1[i][0].Dat.pTl > vect2[i][0].Dat.pTl) {
      while (vect1[i][0].Dat.pTl != vect2[i][k].Dat.pTl) {
        k++;
      }
    }
    else if (vect1[i][0].Dat.pTl < vect2[i][0].Dat.pTl) {
      while (vect1[i][j].Dat.pTl != vect2[i][0].Dat.pTl) {
        j++;
      }
    }
    else {
    }*/  //may need if .Dat.pTl do not match for first value for different data set
    if (vect1[i].size() > vect2[i].size()) {
      size=vect2[i].size();
    }
    else //(vect1[i].size() <= vect2[i].size())
      size=vect1[i].size();
    for (int j=0; j<size;j++){
      pim.Dat.pTSpec=vect1[i][j].Dat.pTSpec;
      pim.Dat.ErrStat=vect1[i][j].Dat.ErrStat;
      pim.Dat.ErrSys=vect1[i][j].Dat.ErrSys;
      pip.Dat.pTSpec=vect2[i][j].Dat.pTSpec;
      pip.Dat.ErrStat=vect2[i][j].Dat.ErrStat;
      pip.Dat.ErrSys=vect2[i][j].Dat.ErrSys;
      pizero.BinIndex=j;
      pizero.Dat.pTl=vect1[i][j].Dat.pTl;
      pizero.Dat.pTh=vect1[i][j].Dat.pTh;
      pizero.Dat.pTSpec= (pim.Dat.pTSpec+pip.Dat.pTSpec)/2;
      pizero.Dat.ErrStat= (pim.Dat.ErrStat+pip.Dat.ErrStat)/2;
      pizero.Dat.ErrSys= (pim.Dat.ErrSys+pip.Dat.ErrSys)/2;
      if (pizero.Dat.pTSpec !=0){
        Pi0[i].push_back(Bin());
        Pi0[i][j]=pizero;
      }
    }
  }
}

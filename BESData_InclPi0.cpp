//BESData_InclPi0.cpp
/*This program will take data from BESData_sorted.txt and include
approximated data for Pi0 particle. built from BESDataToRootFile.cpp
start:20180522.11:00
compl:--
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

void pi0Builder(const vector<vector<Bin>>& vect1, const vector<vector<Bin>>& vect2);

int main(){
  string CollType;
  string CollEner;
  string title0 = "------------------------";
  string title1 = "pTlow [GeV/c]";
  string title2 = "pThigh [GeV/c]";
  string title3 = "d^2N/(2pi*pt*dpt*dy)[(GeV/c)^-2]";
  string title4 = "error(stat.)";
  string title5 = "error(sys.)";
  string title6 = "========================";
  string PartType;
  string CentBin;
  string GARBAGE;
  int BinNum = 0;
  int j;
  int k;
  int l=0;
  double DataPoint;
  vector <Bin> PiM(0);
  vector <Bin> PiP(0);
  vector <vector <Bin>> PiMinus(9);
  vector <vector <Bin>> PiPlus(9);

  ifstream in;
	in.open("./BESData_sorted.txt");
  ofstream out;
  out.open("./BESData_sorted_PlusPi0.txt");

  while(in>>CollType){
    while(CollType!="Au+Au") {in>>CollType;}
    in>>CollEner;
    out<<CollType<<" ";
    out<<CollEner<<" GeV"<<endl;
    for(int i = 0; i<8; i++) {in >> GARBAGE;} //skips up to species
    out <<title0<<endl<<title1<<" "<<title2<<" "<<title3<<" "<<title4<<" "<<title5<<endl;

    for(j = 0; j<9; j++) {
      in>>PartType;
      if (j>=1){
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
        BinNum = 0;
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
        BinNum = 0;
        k=0;
        while (in>>DataPoint){
          BinNum++;
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
    pi0Builder(PiMinus, PiPlus);
    out<<title6<<endl;
    in.clear();
    l++;
  } //end while (in>>CollType)
cout<<"done";
return 0;
} //end int 'main'


//pi0Builder passes in vectors consisting of all pi- and pi+ data points.
//pi0Builder takes average of all points and creates new vector with
//all data for pi0 then writes it to file.
void pi0Builder(const vector<vector<Bin>>& vect1, const vector<vector<Bin>>& vect2){

}

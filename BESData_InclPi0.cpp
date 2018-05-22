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
using namespace std;

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
  struct BinData { //
    double pTl;
    double pTh;
    double pTSpec;
    double ErrStat;
    double ErrSys;
  };
  double DataPoint;
  vector <vector <BinData>> PiMinus(1);
  vector <vector <BinData>> PiPlus(1);
  vector <vector <BinData>> Pi0(1);

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
      out<<PartType<<" ";
      in>>CentBin;
      out<<CentBin<<endl;
      if (PartType == "pi-") {//write in data for vector PiMinus and writes to new file
        BinNum = 0;
        k=0;
        while (in>>DataPoint){ //this will have to have to hold data for either each bin OR have to make j number of vectors
          BinNum++;
          PiMinus[j].push_back(BinData());
            in>>PiMinus[j][k].pTl;
            out<<PiMinus[j][k].pTl<<"\t";
            in>>PiMinus[j][k].pTh;
            out<<PiMinus[j][k].pTh<<"\t";
            in>>PiMinus[j][k].pTSpec;
            out<<PiMinus[j][k].pTSpec<<"\t";
            in>>PiMinus[j][k].ErrStat;
            out<<PiMinus[j][k].ErrStat<<"\t";
            in>>PiMinus[j][k].ErrSys;
            out<<PiMinus[j][k].ErrSys<<endl;
          k++;
        }
      }
      if (PartType =="pi+") { //write in data for vector PiPlus and writes to new file
        BinNum = 0;
        while (in>>DataPoint){
          BinNum++;
          PiPlus[j].push_back(BinData());
            in>>PiPlus[j][k].pTl;
            out<<PiPlus[j][k].pTl;
            in>>PiPlus[j][k].pTh;
            out<<PiPlus[j][k].pTh;
            in>>PiPlus[j][k].pTSpec;
            out<<PiPlus[j][k].pTSpec;
            in>>PiPlus[j][k].ErrStat;
            out<<PiPlus[j][k].ErrStat;
            in>>PiPlus[j][k].ErrSys;
            out<<PiPlus[j][k].ErrSys;
          k++;
        }
      }
      else {//writes to new file
        BinNum = 0;
        while (in>>DataPoint){
          BinNum++;
          out<<DataPoint;
          out<<DataPoint;
          out<<DataPoint;
          out<<DataPoint;
          out<<DataPoint;
          out<<endl;
        }
      }

    } //end forloop ()
    in>>GARBAGE;
    out<<title0<<endl;
    in.clear();
  } //end while (in>>CollType)
cout<<"done";
return 0;
} //end int 'main'

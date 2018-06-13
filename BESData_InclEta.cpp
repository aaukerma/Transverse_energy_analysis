//BESData_InclEta.cpp
/*This program will take data from BESData_sorted_PlusPi0.txt
and estimate the pTSpectra and error for Eta particles. Error will be
assumed high due to the minimal contribution of Eta.
built from BESData_InclPi0.cpp
start:20180529.9:30
compl:
*/

#include <cstdio>
#include <vector>
#include <cmath>
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

void EtaBuilder(const vector<vector<Bin>>& vect1, vector<vector<Bin>>& Eta);

int BESData_InclEta(){
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
  int z=0;
  int EtaSize;
  double DataPoint;
  vector <Bin> Pi0(0);
  vector <vector <Bin>> PiZero(9);
  vector <vector <Bin>> Eta(9);

  ifstream in;
	in.open("./BESData_sorted_PlusPi0.txt"); //Formatting must match previous versions!!!!
  ofstream out;
  out.open("./BESData_sorted_PlusEta.txt"); //created file for inclusion of pi0, formatted!

  while(in>>CollType){
    while(CollType!="Au+Au") {in>>CollType;}
    if (z!=0){
      out<<title0<<endl;
      out<<title6<<endl;
    }
    in>>CollEner;
    out<<CollType<<" ";
    out<<CollEner<<" GeV"<<endl;
    for(int i = 0; i<8; i++) {in >> GARBAGE;} //skips up to species
    out <<title0<<endl<<title1<<" "<<title2<<" "<<title3<<" "<<title4<<" "<<title5<<endl;

    for(j = 0; j<9; j++) {
      if (j!=0) {
        out<<title0<<endl;
      }
      in>>GARBAGE;
      in>>PartType;
      out<<PartType<<"  ";
      in>>CentBin;
      out<<CentBin<<endl;
      if (PartType == "pi0") {//write in data for vector PiMinus and writes to new file
        k=0;
        while (in>>DataPoint){ //This counts as a (in>>) therefore first bit not in>>
          Pi0.push_back(Bin());
          Pi0[k].BinIndex = k;
          Pi0[k].Dat.pTl = DataPoint;
          out<<fixed<<setprecision(2)<<Pi0[k].Dat.pTl<<" \t ";
          in>>Pi0[k].Dat.pTh;
          out<<fixed<<setprecision(2)<<Pi0[k].Dat.pTh<<" \t ";
          in>>Pi0[k].Dat.pTSpec;
          out<<fixed<<setprecision(4)<<Pi0[k].Dat.pTSpec<<" \t ";
          in>>Pi0[k].Dat.ErrStat;
          out<<fixed<<setprecision(5)<<Pi0[k].Dat.ErrStat<<" \t ";
          in>>Pi0[k].Dat.ErrSys;
          out<<fixed<<setprecision(5)<<Pi0[k].Dat.ErrSys<<endl;
          PiZero[j].push_back(Bin());
          PiZero[j][k]=Pi0[k];
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
      in.clear();
    } //end forloop ()
    in>>GARBAGE;
    z++;
    if (l==6){ //determines order in other particles
      EtaBuilder(PiZero, Eta); //runs once for each energy group after pi+ is fed
      out<<title0<<endl<<title6<<endl;
      out<<CollType<<" ";
      out<<CollEner<<" GeV"<<endl;
      out <<title0<<endl<<title1<<" "<<title2<<" "<<title3<<" "<<title4<<" "<<title5<<endl;
      for (int n=0; n<9; n++) { //writes in Pi0
        c1=n*10; //creates bins
        c2=(n+1)*10;
        out<<"eta  "<<c1<<"-"<<c2<<'%'<<endl;
        EtaSize=Eta[n].size();
        for (int m=0; m<EtaSize;m++) {
          if (Eta[n][m].Dat.pTl != 0) {
            out<<fixed<<setprecision(2)<<Eta[n][m].Dat.pTl<<" \t ";
            out<<fixed<<setprecision(2)<<Eta[n][m].Dat.pTh<<" \t ";
            out<<fixed<<setprecision(4)<<Eta[n][m].Dat.pTSpec<<" \t ";
            out<<fixed<<setprecision(5)<<Eta[n][m].Dat.ErrStat<<" \t ";
            out<<fixed<<setprecision(5)<<Eta[n][m].Dat.ErrSys<<endl;
          }
        }
        if (n!=8)
        out<<title0<<endl;
      }
    }
    l++;
    if (l==7) {l=0;}
    in.clear();
  } //end while (in>>CollType)
cout<<"done";
in.close();
out.close();
return 0;
} //end int 'main'


/*EtaBuilder
EtaBuilder takes data of pi0, assumes a scalar relation, and returns data for Eta.
The Scalar value is TODO and relates the pTSpec of pi0 to pTSpec of Eta.
The error is also created here and is TODO
*/
void EtaBuilder(const vector<vector<Bin>>& vect1, vector<vector<Bin>>& Eta){
  int size;
  double pT;
  double c;
  double S=0.482; //scale Factor TODO: find scale factor
  double Se=0.030;
    /*
    S=0.482 based on "Applicability of transverse mass scaling in harmonic collisions at LHC" (page6)
      this has a +- of 0.030 and is based on not pi0 but the average of pi+ and pi- (as was used as estimated pi0 for this program)
      NOTE: this paper was used for the following reasons
                  1: the momentums calculated for were nearest those used in the data sets include (between .025 and 2 GeV/c)
                  2: the data was for higher momentums (where lower momentum data is preferred) in the paper "common suppression pattern of eta and pionzero mesons at high transverse momentum in au+au collisions at snn=200GeV"
                  3: there was no clear correlation presented in "measurement of neutral mesons in p+p collisions as snn=200GeV and scaling properties of hadron production",
                    also I feel it is preferred to use data based in au+au collision data as that is what we are using
    TODO: use PYTHIA v6 to find eta/pi0 vs pT ratio curve (based on "Common Suppression Patterns of [eta] and [pion0] Mesons at High Transverse Momentum in Au+Au collisions at sqrt(s_NN)=200GeV")
      this ratio can be used to calculate S per centrality bin
    come up with solution (defend) [keep it simple]
    another possible assume its same as kaon (incr uncertainty)

    omega pythia based (see alice (bracket to pi0 or eta)) based on ET %
    */
  Bin pi0;
  Bin n;
  for (int i=0; i<9; i++){
    size=vect1[i].size();
    for (int j=0; j<size;j++){
      pi0.Dat.pTSpec=vect1[i][j].Dat.pTSpec;
      pi0.Dat.ErrStat=vect1[i][j].Dat.ErrStat;
      pi0.Dat.ErrSys=vect1[i][j].Dat.ErrSys;
      n.BinIndex=j;
      n.Dat.pTl=vect1[i][j].Dat.pTl;
      n.Dat.pTh=vect1[i][j].Dat.pTh;
      n.Dat.pTSpec= S*(pi0.Dat.pTSpec);
      pT=pow((pi0.Dat.ErrStat*S),2);
      c=pow(Se,2);
      //pTe=pow(pi0.Dat.ErrStat,2);
      //ce=pow(Se,2);
      n.Dat.ErrStat=sqrt(pT+c); //TODO: FIND ERROR STAT rel to pi0 //straight from pion0
      n.Dat.ErrSys=pi0.Dat.ErrSys*S; //TODO: Find ERROR STAT rel to pi0 //straight from pion0
      if (n.Dat.pTSpec !=0){ //solved problem with extra zero entries
        Eta[i].push_back(Bin());
        Eta[i][j]=n; //installs value of eta to vector
      }
    }
  }
}

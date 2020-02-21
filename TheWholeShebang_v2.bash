root -q BESDataToRootFile.cpp
source analyzeAllHistos.bash
root -q BESLambdasToRootFile_TH1.cpp
root -q fitAllHistosInTFile.cpp
root -q interpLaCentTGE.cpp
root -q interpLaCentTGE_y.cpp
root -q finalPlots_TGE_TB.cpp
root -q finalPlots_TGE_y.cpp
root -q stackFinalPlots.C
root -q stackFinalPlots_y.C
root compareWithPHENIX_vs2.C

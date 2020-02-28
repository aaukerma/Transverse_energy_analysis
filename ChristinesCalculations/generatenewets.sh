#!/bin/bash
#Tanner's output is written in Averaged_PiKp_ET_results.txt & Averaged_Lam_ET_resutls.txt
#It creates
#-rw-rw-r-- 1 captainq captainq   7628 Feb 28 14:26 EtResults.txt
#-rw-rw-r-- 1 captainq captainq   7564 Feb 28 14:26 EtResultsN.txt
#-rw-rw-r-- 1 captainq captainq   7526 Feb 28 14:26 EtResultsnpart.txt
g++ -o calc CalcFinalETs.cpp
./calc

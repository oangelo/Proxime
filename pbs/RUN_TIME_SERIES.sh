#!/bin/sh
#PBS -N Chaos
#PBS -M angelo@if.uff.br
#PBS -m bae
#PBS -l nice=16,walltime=1068:00:00 
#PBS -t 0-9
 
cd $PBS_O_WORKDIR
case "$PBS_ARRAYID" in
0)
./chaos time_series Double_Pendulum ab _ns1_ 0.1 0.0 0.0001 50000	
;;
1)
./chaos time_series Double_Pendulum ab _ns2_ 0.2 0.0 0.0001 50000	
;;
2)
./chaos time_series Double_Pendulum ab _ns3_ 0.3 0.0 0.0001 50000	
;;
3)
./chaos time_series Double_Pendulum ab _ns4_ 0.4 0.0 0.0001 50000	
;;
4)
./chaos time_series Double_Pendulum ab _ns5_ 0.5 0.0 0.0001 50000	
;;
5)
./chaos time_series Double_Pendulum ab _s1_ 0.1 0.1 0.0001 50000
;;
6)
./chaos time_series Double_Pendulum ab _s2_ 0.2 0.2 0.0001 50000	
;;
7)
./chaos time_series Double_Pendulum ab _s3_ 0.3 0.3 0.0001 50000	
;;
8)
./chaos time_series Double_Pendulum ab _s4_ 0.4 0.4 0.0001 50000
;;
9)
./chaos time_series Double_Pendulum ab _s5_ 0.5 0.5 0.0001 50000	
;;
esac

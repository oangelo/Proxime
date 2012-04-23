#!/bin/sh
#PBS -N Chaos
#PBS -M angelo@if.uff.br
#PBS -m bae
#PBS -l nice=16,walltime=1068:00:00 
#PBS -t 0-9
 
cd $PBS_O_WORKDIR
case "$PBS_ARRAYID" in
0)
./chaos fft Double_Pendulum rk _s_1_ 0.1 0.1 0.0001 10000
;;
1)
./chaos fft Double_Pendulum rk _s_2_ 0.2 0.2 0.0001 10000	
;;
2)
./chaos fft Double_Pendulum rk _s_3_ 0.3 0.3 0.0001 10000
;;
3)
./chaos fft Double_Pendulum rk _s_4_ 0.4 0.4 0.0001 10000	
;;
4)
./chaos fft Double_Pendulum rk _s_5_ 0.5 0.5 0.0001 10000	
;;
5)
./chaos fft Double_Pendulum rk _ns_1_ 0.1 0.0 0.0001 10000	
;;
6)
./chaos fft Double_Pendulum rk _ns_2_ 0.2 0.0 0.0001 10000	
;;
7)
./chaos fft Double_Pendulum rk _ns_3_ 0.3 0.0 0.0001 10000
;;
8)
./chaos fft Double_Pendulum rk _ns_4_ 0.4 0.0 0.0001 10000
;; 
9)
./chaos fft Double_Pendulum rk _ns_5_ 0.5 0.0 0.0001 10000
;;
esac

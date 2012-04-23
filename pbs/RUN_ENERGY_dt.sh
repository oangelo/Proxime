#!/bin/sh
#PBS -N Chaos
#PBS -M angelo@if.uff.br
#PBS -m bae
#PBS -l nice=16,walltime=1068:00:00 
#PBS -t 0-11
 
cd $PBS_O_WORKDIR
case "$PBS_ARRAYID" in
0)
./chaos energy Double_Pendulum rk _e-5dt_l_ 0.2 0.2 0.00001 10000	
;;
1)
./chaos energy Double_Pendulum rk _e-4dt_l_ 0.2 0.2 0.0001 10000	
;;
2)
./chaos energy Double_Pendulum rk _e-3dt_l_ 0.2 0.2 0.001 10000	
;;
3)
./chaos energy Double_Pendulum ab _e-5dt_l_ 0.2 0.2 0.00001 10000	
;;
4)
./chaos energy Double_Pendulum ab _e-4dt_l_ 0.2 0.2 0.0001 10000	
;;
5)
./chaos energy Double_Pendulum ab _e-3dt_l_ 0.2 0.2 0.001 10000	
;;
6)
./chaos energy Double_Pendulum am _e-5dt_l_ 0.2 0.2 0.00001 10000	
;;
7)
./chaos energy Double_Pendulum am _e-4dt_l_ 0.2 0.2 0.0001 10000	
;;
8)
./chaos energy Double_Pendulum am _e-3dt_l_ 0.2 0.2 0.001 10000	
;;
9)
./chaos energy Double_Pendulum rk _e-6dt_l_ 0.2 0.2 0.000001 10000	
;;
10)
./chaos energy Double_Pendulum ab _e-6dt_l_ 0.2 0.2 0.000001 10000	
;;
11)
./chaos energy Double_Pendulum am _e-6dt_l_ 0.2 0.2 0.000001 10000	
;;

esac

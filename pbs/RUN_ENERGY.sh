#!/bin/sh
#PBS -N Chaos
#PBS -M angelo@if.uff.br
#PBS -m bae
#PBS -l nice=16,walltime=1068:00:00 
#PBS -t 9-26
 
cd $PBS_O_WORKDIR
case "$PBS_ARRAYID" in
0)
./chaos energy Double_Pendulum rk _s1_ 0.1 0.1 0.00001 10000	
;;
1)
./chaos energy Double_Pendulum rk _s2_ 0.2 0.2 0.00001 100000	
;;
2)
./chaos energy Double_Pendulum rk _s3_ 0.3 0.3 0.00001 100000	
;;
3)
./chaos energy Double_Pendulum rk _s4_ 0.4 0.4 0.00001 100000	
;;
4)
./chaos energy Double_Pendulum rk _s5_ 0.5 0.5 0.00001 100000	
;;
5)
./chaos energy Double_Pendulum rk _s6_ 0.6 0.6 0.00001 100000	
;;
6)
./chaos energy Double_Pendulum rk _s7_ 0.7 0.7 0.00001 100000	
;;
7)
./chaos energy Double_Pendulum rk _s8_ 0.8 0.8 0.00001 100000	
;;
8)
./chaos energy Double_Pendulum rk _s9_ 0.9 0.9 0.00001 100000	
;; 

9)
./chaos energy Double_Pendulum ab _s1_ 0.1 0.1 0.00001 10000	
;;
10)
./chaos energy Double_Pendulum ab _s2_ 0.2 0.2 0.00001 100000	
;;
11)
./chaos energy Double_Pendulum ab _s3_ 0.3 0.3 0.00001 100000	
;;
12)
./chaos energy Double_Pendulum ab _s4_ 0.4 0.4 0.00001 100000	
;;
13)
./chaos energy Double_Pendulum ab _s5_ 0.5 0.5 0.00001 100000	
;;
14)
./chaos energy Double_Pendulum ab _s6_ 0.6 0.6 0.00001 100000	
;;
15)
./chaos energy Double_Pendulum ab _s7_ 0.7 0.7 0.00001 100000	
;;
16)
./chaos energy Double_Pendulum ab _s8_ 0.8 0.8 0.00001 100000	
;;
17)
./chaos energy Double_Pendulum ab _s9_ 0.9 0.9 0.00001 100000	
;;

18)
./chaos energy Double_Pendulum am _s1_ 0.1 0.1 0.00001 10000	
;;
19)
./chaos energy Double_Pendulum am _s2_ 0.2 0.2 0.00001 100000	
;;
20)
./chaos energy Double_Pendulum am _s3_ 0.3 0.3 0.00001 100000	
;;
21)
./chaos energy Double_Pendulum am _s4_ 0.4 0.4 0.00001 100000	
;;
22)
./chaos energy Double_Pendulum am _s5_ 0.5 0.5 0.00001 100000	
;;
23)
./chaos energy Double_Pendulum am _s6_ 0.6 0.6 0.00001 100000	
;;
24)
./chaos energy Double_Pendulum am _s7_ 0.7 0.7 0.00001 100000	
;;
25)
./chaos energy Double_Pendulum am _s8_ 0.8 0.8 0.00001 100000	
;;
26)
./chaos energy Double_Pendulum am _s9_ 0.9 0.9 0.00001 100000	
;;
esac

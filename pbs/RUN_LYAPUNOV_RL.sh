#!/bin/sh
#PBS -N Chaos
#PBS -M angelo@if.uff.br
#PBS -m bae
#PBS -l nice=16,walltime=1068:00:00 
#PBS -t 0-5
 
cd $PBS_O_WORKDIR
case "$PBS_ARRAYID" in 
0)
./chaos lyapunov rossler  rk _tab1_ 0.2 0.2 5.7 0.001 9 6	
;;
1)
./chaos lyapunov rossler  rk _tab2_ 0.15 0.2 10 0.001 9 6
;;
2)
./chaos lyapunov rossler  ab _tab1_ 0.2 0.2 5.7 0.001 9 6
;;
3)
./chaos lyapunov rossler  ab _tab2_ 0.15 0.2 10 0.001 9 6
;;
4)
./chaos lyapunov rossler  am _tab1_ 0.2 0.2 5.7 0.001 9 6
;;
5)
./chaos lyapunov rossler  am _tab2_ 0.15 0.2 10 0.001 9 6
;;
6)
./chaos lyapunov lorenz  rk _tab1_  16 40 4 0.001 7 6
;;
7)
./chaos lyapunov lorenz  rk _tab2_  10 45.92 2.666666667 0.001 7 6
;;
8)
./chaos lyapunov lorenz  rk _tab3_  10 28 2.666666667 0.001 7 6
;;
9)
./chaos lyapunov lorenz  rk _tab4_  16 45.92 4 0.001 7 6
;;
10)
./chaos lyapunov lorenz  ab _tab1_  16 40 4 0.001 7 6
;;
11)
./chaos lyapunov lorenz  ab _tab2_  10 45.92 2.666666667 0.001 7 6
;;
12)
./chaos lyapunov lorenz  ab _tab3_  10 28 2.666666667 0.001 7 6
;;
13)
./chaos lyapunov lorenz  ab _tab4_  16 45.92 4 0.001 7 6
;;
14)
./chaos lyapunov lorenz  am _tab1_  16 40 4 0.001 7 6
;;
15)
./chaos lyapunov lorenz  am _tab2_  10 45.92 2.666666667 0.001 7 6
;;
16)
./chaos lyapunov lorenz  am _tab3_  10 28 2.666666667 0.001 7 6
;;
17)
./chaos lyapunov lorenz  am _tab4_  16 45.92 4 0.001 7 6
;;
esac

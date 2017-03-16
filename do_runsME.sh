#!/bin/bash

mpirun -n 1 ./aDGVM 1 1 29 30 -25 5 &
mpirun -n 1 ./aDGVM 2 1 29 30 -25 5 &
mpirun -n 1 ./aDGVM 3 1 29 30 -25 5 &
mpirun -n 1 ./aDGVM 4 1 29 30 -25 5 &
mpirun -n 1 ./aDGVM 5 1 29 30 -25 5 &
mpirun -n 1 ./aDGVM 6 1 29 30 -25 5 &
mpirun -n 1 ./aDGVM 7 1 29 30 -25 5 &
mpirun -n 1 ./aDGVM 8 1 29 30 -25 5 &
mpirun -n 1 ./aDGVM 9 1 29 30 -25 5 &
mpirun -n 1 ./aDGVM 10 1 29 30 -25 5 &


#echo "FIRST 4 RUNS OF SERIES DONE"

#mpirun -n 4 ./aDGVM 145 1 25 27 -23 -21 &
#mpirun -n 4 ./aDGVM 146 1 25 27 -23 -21 &
#mpirun -n 4 ./aDGVM 147 1 25 27 -23 -21 &
#mpirun -n 4 ./aDGVM 148 1 25 27 -23 -21 

#echo "SECOND 4 RUNS OF SERIES DONE, storing files away..."

#cd demand_182.5
#mkdir fire_grazing_freq0.4
#cd ..
#mv pop_14*.nc demand_182.5/fire_grazing_freq0.4
#mv trait_14*.nc demand_182.5/fire_grazing_freq0.4

#echo "DONE!"

#exit

#mv runconfig.cfg runconfig_fq0.4.cfg
#mv runconfig_fq0.2 runconfig.cfg





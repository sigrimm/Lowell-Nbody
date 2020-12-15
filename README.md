# Lowell-Nbody

A N-body integrator for small bodies of the Solar System
Author: Simon Grimm, University of Bern, Switzerland

# Instructions #

## Step 1, get initial conditions ##

- cd TestSet
- run:

python3 jpl_vectors.py <Asteroid_Name>

with the desired asteroid name, e.g. python3 jpl_vectors.py Icarus

This will downlaod the barycentric coordinates for the JPL Horizons database,
and store it in a file name <Aseroid_Name_d>.dat 

-inside jpl_vectors.py, the flag `useHelioCentric`can be used to swith from barycentric to heliocentric coordinates.

## Step 1B, get perturbers file ##
This is only necessary, if the perturbers file needs to be updated.

- use `python3 jpl_vectorsAll.py` to download all perturbers.
- run `python3 merge.py > All_b.dat` to merge all perturbers onto one file
- copy All_b.dat to the directory 'combined'


## Step 2, copy initial conditions ##

-Chose the starting time of the integration, this can not be earlier than 2433286.5
-Copy the corresponding line from <Aseroid_Name_d>.dat to a file named initial.dat
- e.g. 2433286.5 1.046239393196292 -1.497371446392998 -0.4605849986210472 0.002079418803054296 0.005857305828677572 -0.0008158440936451222

-Copy the file initial.dat to the directory 'combined'

## Step 3, run the integration ##

-cd ../combined

- compile with g++ -o RKF45 RKF45.cpp

- ./RKF45 (options)
- options are:
- `-Nsteps` <integer> , default = 4000000, number of time steps.
- `-outInterval` <integer> , default = 1000, interval out output files.
- `-dt` <float> , default = 0.001, time step in days.

This will run the integration and produce an output (barycentric) (every day for default options). The output file names include the number of the time step. 

The output files contain 
- time in days
- index
- mass in Solar masses
- x position in AU
- y position in AU
- z position in AU
- vx velocity in AU/day
- vy velocity in AU/day
- vz velocity in AU/day

## Step 4 compare the results ##
- combine the output files with

cat Out* > out<name>.dat

e.g. cat Out* > outIcarus.dat

-run python3 compare.py > diffIcarus.dat
(The skiprow arguments needs to be adapted to the initial start time)

This produces a file with the difference between the real positions and the integration.
The file contains two columns:
-- time in day
-- difference in meters


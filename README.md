# Lowell-Nbody

A N-body integrator for small bodies of the Solar System
Author: Simon Grimm, University of Bern, Switzerland

# Instructions #

## Step 1A, get initial conditions ##

- cd TestSet1850
- run:

python3 jpl_vectors.py <Asteroid_Name>

with the desired asteroid name, e.g. python3 jpl_vectors.py Icarus

This will download the heliocentric coordinates for the JPL Horizons database,
and store it in a file name <Asteroid_Name_d>_h.dat 

-inside jpl_vectors.py, the flag `useHelioCentric`can be used to swith from heliocentric to barycentric coordinates.

## Step 1B, get perturbers file ##
This is only necessary, if the perturbers file needs to be updated.

- use `python3 jpl_vectorsAll.py` to download all perturbers.
  This will download all perturbers coordinates and store them in files <Perturber_Name>_h.dat
- run `python3 merge.py > All_h.dat` to merge all perturbers onto one file
- copy All_h.dat to the directory 'combined'


## Step 2, copy initial conditions ##
In the TestSet1850 directory

-Chose the starting time of the integration, this can not be earlier than 2396770.5
-Copy the corresponding line from <Asteroid_Name_h>.dat in the TestSet1850 directory to a file named initial.dat
- e.g. for asteroid 105140:

2396770.5 0.752060381543369849133284787968 -0.947384215699796805587595827092 0.300046271431012923081027565786 0.0107856006406296302951863808062 -0.00245802678537547720019618147091 0.00537247903456764940716139378196 0 0 0 0.11126204259999999957 4.6142000000000003013 2.1499999999999999112 5.0929999999999999716 2.8079999999999998295


-Copy the file initial.dat to the directory 'combined'

## Step 3, run the integration ##

-cd ../combined

- compile with g++ -o RKF45 RKF45.cpp

-run: ./RKF45 [options]
- options are:
- `-Nsteps` < integer >, default = 400000, number of time steps.
- `-outInterval` < integer >, default = 1000, interval out output files.
- `-dt` < float >, default = 0.001, time step in days.

This will run the integration and produce an output (heliocentric) (every day for default options). The output file names include the number of the time step. 

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

## Step 4 compare the results with JPL ##
- combine the output files with

cat Out* > out< name >_h.dat

e.g. cat Out* > outIcarus_h.dat

-run python3 compare.py < name > > diff< name >.dat

e.g. python3 compare.py Icarus > diffIcarus.dat

This produces a file with the difference between the real positions from JPL and the integration.
The file contains two columns:
-- time in day
-- difference in meters


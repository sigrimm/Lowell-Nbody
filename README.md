# Lowell-Nbody

A N-body integrator for small bodies of the Solar System
Author: Simon Grimm, ETH Zürich, University of Zürich, Switzerland


# Instructions #

## Step 1, get perturbers file ##

- cd readChebyshev
- Download linux_p1550p2650.440 
- Download header.440
- Download sb441-n16.bsp
- make
- set start time and end time of desired perturbers data in paramChebyshev.dat
- ./readChebyshev


This will create a binary file called 'PerturbersChebyshev.bin'
The structure of the file is:

- The header contains all necessary constants:

  - time0,   fp64, Start time of the entire file, in Julian Days, default = 2450800.5. 
  - time1,   fp64, End time of the entire file, in Julian Days, default = 2459800.5
  - AUtokm,  fp64, Conversion factor AU to km, taken from DE440 file.
  - EM,      fp64, Earth / Moon mass ratio, taken from DE440 file.
  - CLIGHT,  fp64, Speed of light in km/s, taken from DE440 file.
  - RE,      fp64, Earth radius in km, taken from DE440 file.
  - J2E,     fp64, J2 term of Earth (dimensionless), taken from DE440 file.

- After the header follows a list of all perturbers. Containing for each perturber:

  - id,         int32, perturber id
  - nChebyshev, int32, number of Chebyshev coefficients
  - offset0 ,   int32, offset of begiining of data from perturber
  - offset1,    int32, offset of ending of data from perturber
  - GM,         fp64,  G * mass 

- Now follows the data itself. 

  - perturber 1, record 1
  - perturber 1, record 2
  - ...
  - perturber 1, record n
  - perturber 2, record 1
  - perturber 2, record 2
  - ...
  - perturber 2, record n
  - ...

Every record contains

  - start time of record
  - end time of record
  - Chebyshev coefficients, (3 * nChebyshev)


## Step 2, set parameters ##

- cd ../integrator
- make
- modify param.dat


allowed integrators are:

- LF, leapfrog 
- RK4, Runge Kutta 4, fixed time steps
- RKF45, Runge Kutta Fehlberg 45, adaptive time steps
- DP54, Runge Kutta Fehlberg 54, adaptive time steps
- RKF78, Runge Kutta Fehlberg 78, adaptive time steps

## Step 3, set initial conditions ##

- file name specified in param.dat

- x y z vx vy vz A1 A2 A3
- in Barycentric coordinates, AU and AU/day

## Step 4, run the integration ##

- ./integrator

- outputs are in Out.dat

The output files contain 
- time in days
- x position in AU
- y position in AU
- z position in AU
- vx velocity in AU/day
- vy velocity in AU/day
- vz velocity in AU/day
- minimal time step of interval
- in Barycentric coordinates

<!---
## Step 5 compare the results with JPL ##
- combine the output files with

cat Out* > out< name >_h.dat

e.g. cat Out* > outIcarus_h.dat

-run python3 compare.py < name > > diff< name >.dat

e.g. python3 compare.py Icarus > diffIcarus.dat

This produces a file with the difference between the real positions from JPL and the integration.
The file contains two columns:
-- time in day
-- difference in meters
-->


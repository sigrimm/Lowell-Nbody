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
- ./readChebyshev

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


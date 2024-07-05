#include <math.h>
#include "asteroid.h"
#include "force.h"
#include "integrator.h"
#include "Chebyshev.h"

int main(){


	asteroid A;

	//A.infile = fopen("PerturbersChebyshev.dat", "r");
	A.infile = fopen("PerturbersChebyshev.bin", "rb");


	A.Nperturbers = 27;
	A.time_reference = 2451545.0;
	A.timeStart = 8255.5;		//start time of integration
	A.timeEnd = -11744.5;		//end time of integration
	A.dt = -0.01;


	A.allocate();


	// *********************************
	//Set initial conditions

	//Barycentric coordinates of 15TC25   
	A.x[A.Nperturbers + 0] = -0.621994840672447701912517459277;   // x in AU
	A.y[A.Nperturbers + 0] = -0.828032662228602056586623803014;  // y
	A.z[A.Nperturbers + 0] = -0.406943193813317449780697643291;  // z
	A.vx[A.Nperturbers + 0] = 0.0122692503194383149833779356186;   // vx in AU/day
	A.vy[A.Nperturbers + 0] = -0.00859370367531481910150503722434;   // vy
	A.vz[A.Nperturbers + 0] = -0.0046654983615674223973446288482;  // vz

	//Non-grav terms
	A.A1[A.Nperturbers + 0] = 1.599229127169e-10;
	A.A2[A.Nperturbers + 0] = -5.273644346744e-12;
	A.A3[A.Nperturbers + 0] = 0.0;


	// *********************************

	int er = A.loop();	

	fclose(A.infile);
	
}


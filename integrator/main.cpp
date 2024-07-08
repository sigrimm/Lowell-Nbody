#include <math.h>
#include "asteroid.h"
#include "force.h"
#include "RKF.h"
#include "integrator.h"
#include "Chebyshev.h"

int main(){

	asteroid A;

	//A.perturbersFile = fopen("../readChebyshev/PerturbersChebyshev.dat", "r");
	A.perturbersFile = fopen("../readChebyshev/PerturbersChebyshev.bin", "rb");
	if(A.perturbersFile == NULL){
		printf("Error, perturbers file %s not found\n", "../readChebyshev/PerturbersChebyshev.bin");
		return 0;
	}


	A.Nperturbers = 27;

	//set default parameters
	//----------------------------------------------------------
	A.time_reference = 2451545.0;
	A.timeStart = 8255.5;		//start time of integration
	A.timeEnd = -11744.5;		//end time of integration
	A.dt = -0.01;
	A.outputInterval = 10.0;
	A.RKF_atol = 1.0e-16;
	A.RKF_rtol = 1.0e-16;
	A.RKF_fac = 0.84;
	A.RKF_facmin = 0.8;
	A.RKF_facmax = 1.5;
	A.RKFn = 6;
	//----------------------------------------------------------

	int er = 0;

	er = A.readParam();
	if(er <= 0){
		return 0;
	}
	printf("Read param.dat file OK\n");
	printf("Initial condition file = %s\n", A.inputFilename);

	A.timeStart -= A.time_reference;
	A.timeEnd -= A.time_reference;

	//set integrator properties

	A.allocate();
	printf("Allocate memory OK\n");
	
	if(A.RKFn == 4){
		A.setRK4();
	}
	if(A.RKFn == 6){
		A.setRKF45();
	}
	if(A.RKFn == 7){
		A.setDP54();
	}
	if(A.RKFn == 13){
		A.setRKF78();
	}


	// *********************************
	//Set initial conditions
	er = A.readIC();

	if(er <= 0){
		return 0;
	}
	printf("Read initial conditions file OK\n");

	/*
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
	*/
	// *********************************

	er = A.loop();	

	fclose(A.perturbersFile);
	
}


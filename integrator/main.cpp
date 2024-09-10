#include "asteroid.h"

int main(){


#if USEGPU == 1
	printf("Use GPU %d\n", USEGPU);
#endif

	asteroid A;



	A.Nperturbers = 27;


	if(A.Nperturbers > def_NP){
		printf("Error, Number of perturbers is larger than def_NP, increase def_NP. %d %d\n", A.Nperturbers, def_NP);
		return 0;
	}


	//----------------------------------------------------------
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
	sprintf(A.perturbersFilePath, "%s", "../readChebyshev/");
	//----------------------------------------------------------

	int er = 0;

	//----------------------------------------------------------
	//Read param.dat file
	//----------------------------------------------------------
	printf("Read param.dat file\n");
	er = A.readParam();
	if(er <= 0){
		return 0;
	}
	printf("Read param.dat file OK\n");
	//----------------------------------------------------------


	//sprintf(A.perturbersFileName, "%s/PerturbersChebyshev.dat", A.perturbersFilePath);
	//A.perturbersFile = fopen(A.perturbersFileName, "r");
	sprintf(A.perturbersFileName, "%s/PerturbersChebyshev.bin", A.perturbersFilePath);

	printf("Open Perturbers file = %s\n", A.perturbersFileName);

	A.perturbersFile = fopen(A.perturbersFileName, "rb");
	if(A.perturbersFile == NULL){
		printf("Error, perturbers file not found: %s\n",  A.perturbersFileName);
		return 0;
	}
	printf("Open Perturbers file OK\n");


	A.timeStart -= A.time_reference;
	A.timeEnd -= A.time_reference;

	printf("Allocate memory\n");
	A.allocate();
#if USEGPU == 1
	er = A.allocateGPU();
	if(er <= 0){
		return 0;
	}
#endif
	printf("Allocate memory OK\n");

#if USEGPU == 1
	printf("Read data table\n");
	er = A.readData();

	if(er <= 0){
		return 0;
	}

	printf("Read data table OK\n");

#endif
	
	//set integrator properties
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


	//----------------------------------------------------------
	//Read initial conditions
	//----------------------------------------------------------
	printf("Read initial conditions file = %s\n", A.inputFilename);
	er = A.readIC();

	if(er <= 0){
		return 0;
	}

#if USEGPU == 1
	er = A.copyIC();

	if(er <= 0){
		return 0;
	}
#endif
	printf("Read initial conditions file OK\n");
	//----------------------------------------------------------

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


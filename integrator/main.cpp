#include "asteroid.h"
#include <chrono>

int main(int argc, char*argv[]){

	std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

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
	sprintf(A.perturbersFilePath, "%s", "../readChebyshev/");
	sprintf(A.name, "%s", "test");
	//----------------------------------------------------------

	int er = 0;

	//----------------------------------------------------------
	//Read param.dat file
	//----------------------------------------------------------
	printf("Read param.dat file\n");
	er = A.readParam(argc, argv);
	if(er <= 0){
		return 0;
	}
	printf("Read param.dat file OK\n");
	//----------------------------------------------------------

	//----------------------------------------------------------
	//Read Size of initial conditions file
	//----------------------------------------------------------
	if(A.ICformat == 0){
		A.inputFile = fopen(A.inputFilename, "r");
		if(A.inputFile == NULL){
			printf("Error, could not open initial condition file |%s|\n", A.inputFilename);
			return 0;
		}

		printf("Read ICsize\n");
		er = A.readICSize();
		if(er <= 0){
			return 0;
		}
		fclose(A.inputFile);
		printf("Read ICSize OK with %d bodies\n", A.N);
	}
	else{

		A.inputFile = fopen(A.inputFilename, "rb");
		if(A.inputFile == NULL){
			printf("Error, could not open initial condition file |%s|\n", A.inputFilename);
			return 0;
		}
		printf("Read IC file header\n");
		er = A.readHeader();
		if(er <= 0){
			return 0;
		}
		printf("Read IC file header OK with %d bodies\n", A.N);
	}
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

	A.infoFile = fopen(A.infoFilename, "w");
	A.printInfo();


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
	if(A.ICformat == 0){
		A.inputFile = fopen(A.inputFilename, "r");
		er = A.readIC();
	}
	else{
		er = A.readFile();
	}

	if(er <= 0){
		return 0;
	}
	fclose(A.inputFile);

	A.timeStart -= A.time_reference;
	A.timeEnd -= A.time_reference;
	A.time = A.timeStart;

#if USEGPU == 1
	er = A.copyIC();
	if(er <= 0){
		return 0;
	}
	er = A.copyConst();
	if(er <= 0){
		return 0;
	}
#endif
	printf("Read initial conditions file OK\n");
	//----------------------------------------------------------

	er = A.loop();	

	fclose(A.perturbersFile);
	fclose(A.infoFile);

	std::chrono::steady_clock::time_point time_end = std::chrono::steady_clock::now();

	int ms = 	std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_begin).count();
	//int ns = 	std::chrono::duration_cast<std::chrono::nanoseconds> (time_end - time_begin).count();
	printf("Run time in seconds:  %g\n", ms / 1000.0);
	
}


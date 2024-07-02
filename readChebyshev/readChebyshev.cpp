#include "define.h"
#include "planets.h"
#include "perturbers.h"


int main(){


	double time0 = 2450800.5;
	double time1 = 2459800.5;


	planets pl;
	perturbers pert;


	pl.Nplanets = 11;
	pl.Npert = 16;
	pert.Npert = 16;



	pl.alloc();
	pert.alloc();

	pl.id[0] = 0;
	pl.id[1] = 1;
	pl.id[2] = 2;
	pl.id[3] = 3;
	pl.id[4] = 4;
	pl.id[5] = 5;
	pl.id[6] = 6;
	pl.id[7] = 7;
	pl.id[8] = 8;
	pl.id[9] = 9;
	pl.id[10] = 10;


	pert.id[0] = 2000107;	//Camilla
	pert.id[1] = 2000001;	//Ceres
	pert.id[2] = 2000065;	//Cybele
	pert.id[3] = 2000511;	//Davida
	pert.id[4] = 2000015;	//Eunomia
	pert.id[5] = 2000031;	//Euphosyne
	pert.id[6] = 2000052;	//Europa
	pert.id[7] = 2000010;	//Hygiea
	pert.id[8] = 2000704;	//Interamnia
	pert.id[9] = 2000007;	//Iris
	pert.id[10] = 2000003;	//Juno
	pert.id[11] = 2000002;	//Pallas
	pert.id[12] = 2000016;	//Psyche
	pert.id[13] = 2000087;	//Sylvia
	pert.id[14] = 2000088;	//Thisbe
	pert.id[15] = 2000004;	//Vesta


	for(int i = 0; i < pl.Npert; ++i){
		pl.id[pl.Nplanets + i] = pert.id[i] - 2000000;
	}


	char headerFileName[256];
	char planetsFileName[256];
	char perturbersFileName[256];
	char outFileName[256];

	sprintf(headerFileName, "header.440");
	sprintf(planetsFileName, "linux_p1550p2650.440");
	sprintf(perturbersFileName, "sb441-n16.bsp");
	//sprintf(outFileName, "PerturbersChebyshev.dat");
	sprintf(outFileName, "PerturbersChebyshev.bin");

	


	FILE *headerFile;
	FILE *planetsFile;
	FILE *perturbersFile;


	FILE *outFile;


	//outFile = fopen(outFileName, "w");
	outFile = fopen(outFileName, "wb");

	int er;

	// **********************************************
	// Planets
	// **********************************************

	headerFile = fopen(headerFileName, "r");
	planetsFile = fopen(planetsFileName, "rb");

	er = pl.readHeader(headerFile);
	er = pl.readPlanets(planetsFile, outFile, time0, time1);

	fclose(headerFile);
	fclose(planetsFile);

	for(int i = 0; i < pl.Npert; ++i){
		pert.GM[i] = pl.GM[pl.Nplanets + i];
	}

	// **********************************************
	// Perturbers
	// **********************************************


	perturbersFile = fopen(perturbersFileName, "rb");

	er = pert.readPerturbers1(perturbersFile);
	er = pert.readPerturbers2(perturbersFile, outFile, time0, time1, pl.dataSize);

	fclose(perturbersFile);


	er = pl.printPlanets(outFile);
	er = pert.printPerturbers(outFile);

	fclose(outFile);


	return 0;

}


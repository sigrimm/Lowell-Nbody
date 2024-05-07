#include "define.h"
#include "planets.h"
#include "perturbers.h"


int main(){


	double time0 = 2450800.5;
	double time1 = 2459800.5;


	planets pl;
	perturbers pert;


	pl.Nplanets = 11;
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
	
	pert.id[0] = 2000107;
	pert.id[1] = 2000001;
	pert.id[2] = 2000065;
	pert.id[3] = 2000511;
	pert.id[4] = 2000015;
	pert.id[5] = 2000031;
	pert.id[6] = 2000052;
	pert.id[7] = 2000010;
	pert.id[8] = 2000704;
	pert.id[9] = 2000007;
	pert.id[10] = 2000003;
	pert.id[11] = 2000002;
	pert.id[12] = 2000016;
	pert.id[13] = 2000087;
	pert.id[14] = 2000088;
	pert.id[15] = 2000004;



	char headerFileName[256];
	char planetsFileName[256];
	char perturbersFileName[256];
	char outFileName[256];

	sprintf(headerFileName, "header.440");
	sprintf(planetsFileName, "linux_p1550p2650.440");
	sprintf(perturbersFileName, "sb441-n16.bsp");
	sprintf(outFileName, "PerturbersChebyshev.dat");

	


	FILE *headerFile;
	FILE *planetsFile;
	FILE *perturbersFile;


	FILE *outFile;


	headerFile = fopen(headerFileName, "r");
	planetsFile = fopen(planetsFileName, "rb");
	perturbersFile = fopen(perturbersFileName, "rb");
	outFile = fopen(outFileName, "w");



	int er;

	er = pl.readHeader(headerFile);
	er = pl.readPlanets(planetsFile, time0, time1);
	er = pl.printPlanets(outFile);



	er = pert.readPerturbers1(perturbersFile);
	er = pert.readPerturbers2(perturbersFile, time0, time1);
	er = pert.printPerturbers(outFile);


	fclose(headerFile);
	fclose(planetsFile);
	fclose(perturbersFile);
	fclose(outFile);



	return 0;

}


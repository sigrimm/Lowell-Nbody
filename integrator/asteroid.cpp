#include "asteroid.h"


void asteroid::allocate(){
	N = Nperturbers + 1;

	infile = fopen("PerturbersChebyshev.bin", "rb");
	nCm = 0;
	dts = (dt > 0.0) ? 1.0 : -1.0;      //sign of time step
	dt1 = dt;
	stop = 0;

	startTime = (double*)malloc(Nperturbers * sizeof(double));
	endTime = (double*)malloc(Nperturbers * sizeof(double));
	id = (int*)malloc(Nperturbers * sizeof(int));
	nChebyshev = (int*)malloc(Nperturbers * sizeof(int));
	offset0 = (int*)malloc(Nperturbers * sizeof(int));
	offset1 = (int*)malloc(Nperturbers * sizeof(int));
	GM = (double*)malloc(Nperturbers * sizeof(double));

	//read header
	fread(&time0, sizeof(double), 1, infile);
	fread(&time1, sizeof(double), 1, infile);
	fread(&AUtokm, sizeof(double), 1, infile);
	fread(&EM, sizeof(double), 1, infile);
	fread(&CLIGHT, sizeof(double), 1, infile); 
	fread(&RE, sizeof(double), 1, infile);
	fread(&J2E, sizeof(double), 1, infile);

	for(int i = 0; i < Nperturbers; ++i){
		fread(&id[i], sizeof(int), 1, infile);
		fread(&nChebyshev[i], sizeof(int), 1, infile);
		fread(&offset0[i], sizeof(int), 1, infile);
		fread(&offset1[i], sizeof(int), 1, infile);
		fread(&GM[i], sizeof(double), 1, infile);

		nCm = (nCm > nChebyshev[i]) ? nCm : nChebyshev[i];

		offset0[i] += (3 * Nperturbers + 7);    //add size of header 7*double + Nperturbers * (4 int + double)
		offset1[i] += (3 * Nperturbers + 7);    //add size of header

		startTime[i] = 100000000.0;     //large number
		endTime[i] = 0; 
	//printf("%d %d %d %d %.20g\n", id[i], nChebyshev[i], offset0[i], offset1[i], GM[i]);
	}  

	cdata = (double*)malloc(Nperturbers * nCm * 3 * sizeof(double));

	time = timeStart;
	double c = (CLIGHT / AUtokm) * 86400.0;
	c2 = c * c;
	REAU = RE / AUtokm;   //Earth radius in AU


	x = (double*)malloc(N * sizeof(double));
	y = (double*)malloc(N * sizeof(double));
	z = (double*)malloc(N * sizeof(double));

	vx = (double*)malloc(N * sizeof(double));
	vy = (double*)malloc(N * sizeof(double));
	vz = (double*)malloc(N * sizeof(double));

	ax = (double*)malloc(N * sizeof(double));
	ay = (double*)malloc(N * sizeof(double));
	az = (double*)malloc(N * sizeof(double));

	A1 = (double*)malloc(N * sizeof(double));
	A2 = (double*)malloc(N * sizeof(double));
	A3 = (double*)malloc(N * sizeof(double));
}



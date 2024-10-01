#include "asteroid.h"
#include <stdio.h>


int asteroid::readParam(){

	FILE *paramfile;
	paramfile = fopen("param.dat", "r");

	if(paramfile == NULL){
		printf("Error, param.dat file does not exist\n");
		return 0;
	}


	char sp[160];
	int er;
	char *str;	//Needed for return value of fgest, otherwise a compiler warning is generated

	for(int j = 0; j < 1000; ++j){ //loop around all lines in the param.dat file
		int c;
		for(int i = 0; i < 50; ++i){
			c = fgetc(paramfile);
			if(c == EOF){
				break;
			}
			sp[i] = char(c);
			if(c == '=' || c == ':'){
				sp[i + 1] = '\0';
				break;
			}
		}
		if(c == EOF) break;
		if(strcmp(sp, "Initial condition file =") == 0){
			er = fscanf (paramfile, "%s", inputFilename);
			if(er <= 0){
				printf("Error: Initial condition file is not valid!\n");
				return 0;
			}
		str = fgets(sp, 3, paramfile);
		}
		else if(strcmp(sp, "Reference Time =") == 0){
			er = fscanf (paramfile, "%lf", &time_reference);
			if(er <= 0){
				printf("Error: Reference Time value is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
		}
		else if(strcmp(sp, "Start Time =") == 0){
			er = fscanf (paramfile, "%lf", &timeStart);
			if(er <= 0){
				printf("Error: Start Time value is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
		}
		else if(strcmp(sp, "End Time =") == 0){
			er = fscanf (paramfile, "%lf", &timeEnd);
			if(er <= 0){
				printf("Error: End Time value is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
		}
		else if(strcmp(sp, "Output Interval =") == 0){
			er = fscanf (paramfile, "%lf", &outputInterval);
			if(er <= 0){
				printf("Error: Output Interval value is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
		}
		else if(strcmp(sp, "Integrator =") == 0){
			char integratorName[160];
			er = fscanf (paramfile, "%s", integratorName);
			if(er <= 0){
				printf("Error: Integrator Name is not valid!\n");
				return 0;
			}
			if(strcmp(integratorName, "LF") == 0){
				RKFn = 1;
			}
			else if(strcmp(integratorName, "RK4") == 0){
				RKFn = 4;
			}
			else if(strcmp(integratorName, "RKF45") == 0){
				RKFn = 6;
			}
			else if(strcmp(integratorName, "DP54") == 0){
				RKFn = 7;
			}
			else if(strcmp(integratorName, "RKF78") == 0){
				RKFn = 13;
			}
			else{
				printf("Errof, Integrator not valid\n");
				return 0;

			}
			str = fgets(sp, 3, paramfile);
		}
		else if(strcmp(sp, "RKF absolute tolerance =") == 0){
			er = fscanf (paramfile, "%lf", &RKF_atol);
			if(er <= 0){
				printf("Error: RKF absolute tolerance value is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
		}
		else if(strcmp(sp, "RKF relative tolerance =") == 0){
			er = fscanf (paramfile, "%lf", &RKF_rtol);
			if(er <= 0){
				printf("Error: RKF relative tolerance value is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
		}
		else if(strcmp(sp, "Time Step =") == 0){
			er = fscanf (paramfile, "%lf", &dt);
			if(er <= 0){
				printf("Error: dt value is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
		}
		else if(strcmp(sp, "Path to perturbers file =") == 0){
			er = fscanf (paramfile, "%s", perturbersFilePath);
			if(er <= 0){
				printf("Error: Path to perturbers file is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
		}

		else{
			printf("Error: param.dat file is not valid! %s\n", sp);
			return 0;
		}
	}
	fclose(paramfile);

	return 1;

}

//Read the size of the initial conditions file
int asteroid::readICSize(){
	FILE *infile;

	//infile = fopen(inputFilename, "r");
	infile = fopen("initial.dat", "r");
	if(infile == NULL){
		printf("Error, could not open initial condition file |%s|\n", inputFilename);
		return 0;
	}

	N = 0;

	double x, y, z;
	double vx, vy, vz;
	double A1, A2, A3;

	for(int i = 0; i < 1024 * 1024; ++i){
		int er = 0;
		er = fscanf(infile, "%lf", &x);
		er = fscanf(infile, "%lf", &y);
		er = fscanf(infile, "%lf", &z);
		er = fscanf(infile, "%lf", &vx);
		er = fscanf(infile, "%lf", &vy);
		er = fscanf(infile, "%lf", &vz);
		er = fscanf(infile, "%lf", &A1);
		er = fscanf(infile, "%lf", &A2);
		er = fscanf(infile, "%lf", &A3);
		if(er < 0){
			break;
		}
		if(i >= 1024 * 1024 - 1){
			printf("Error, N is too large for scan kernels\n");
			return 0;
		}
		++N;
	//printf("xyz %.40g %.40g %.40g %.40g %.40g %.40g %.40g %.40g %.40g\n", x_h[i], y_h[i], z_h[i], vx_h[i], vy_h[i], vz_h[i], A1_h[i], A2_h[i], A3_h[i]);
	}
	fclose(infile);

	return 1;

}

//Read the initial conditions file
int asteroid::readIC(){
	FILE *infile;

	//infile = fopen(inputFilename, "r");
	infile = fopen("initial.dat", "r");
	if(infile == NULL){
		printf("Error, could not open initial condition file |%s|\n", inputFilename);
		return 0;
	}

	for(int i = 0; i < N; ++i){
		int er = 0;
		er = fscanf(infile, "%lf", &x_h[i]);
		er = fscanf(infile, "%lf", &y_h[i]);
		er = fscanf(infile, "%lf", &z_h[i]);
		er = fscanf(infile, "%lf", &vx_h[i]);
		er = fscanf(infile, "%lf", &vy_h[i]);
		er = fscanf(infile, "%lf", &vz_h[i]);
		er = fscanf(infile, "%lf", &A1_h[i]);
		er = fscanf(infile, "%lf", &A2_h[i]);
		er = fscanf(infile, "%lf", &A3_h[i]);
		if(er < 0){
			printf("Error, reading initial conditions file failed.\n");
			return 0;
		}
	//printf("xyz %.40g %.40g %.40g %.40g %.40g %.40g %.40g %.40g %.40g\n", x_h[i], y_h[i], z_h[i], vx_h[i], vy_h[i], vz_h[i], A1_h[i], A2_h[i], A3_h[i]);
	}
	fclose(infile);

	return 1;

}


void asteroid::allocate(){
	nCm = 0;
	dts = (dt > 0.0) ? 1.0 : -1.0;      //sign of time step
	dt1 = dt;
	stop = 0;

	id_h = (int*)malloc(Nperturbers * sizeof(int));
	nChebyshev_h = (int*)malloc(Nperturbers * sizeof(int));
#if USEGPU == 0
	startTime_h = (double*)malloc(Nperturbers * sizeof(double));
	endTime_h = (double*)malloc(Nperturbers * sizeof(double));
	offset0_h = (int*)malloc(Nperturbers * sizeof(int));
	offset1_h = (int*)malloc(Nperturbers * sizeof(int));
#else
	//When using the GPU, then calcualting the perturbers position can be
	//parallelized also along the number of stages in the Runge-Kutta integrator
	//Therefore we need for every perturber and every stage the offset,
	//Start time and End time of the current record
	startTime_h = (double*)malloc(Nperturbers * RKFn * sizeof(double));
	endTime_h = (double*)malloc(Nperturbers * RKFn * sizeof(double));
	offset0_h = (int*)malloc(Nperturbers * RKFn * sizeof(int));
	offset1_h = (int*)malloc(Nperturbers * RKFn * sizeof(int));
#endif
	GM_h = (double*)malloc(Nperturbers * sizeof(double));

	//read header
	int er;
	er = fread(&time0, sizeof(double), 1, perturbersFile);
	er = fread(&time1, sizeof(double), 1, perturbersFile);
	er = fread(&AUtokm, sizeof(double), 1, perturbersFile);
	er = fread(&EM, sizeof(double), 1, perturbersFile);
	er = fread(&CLIGHT, sizeof(double), 1, perturbersFile); 
	er = fread(&RE, sizeof(double), 1, perturbersFile);
	er = fread(&J2E, sizeof(double), 1, perturbersFile);

#if USEGPU == 0
	for(int i = 0; i < Nperturbers; ++i){
		er = fread(&id_h[i], sizeof(int), 1, perturbersFile);
		er = fread(&nChebyshev_h[i], sizeof(int), 1, perturbersFile);
		er = fread(&offset0_h[i], sizeof(int), 1, perturbersFile);
		er = fread(&offset1_h[i], sizeof(int), 1, perturbersFile);
		er = fread(&GM_h[i], sizeof(double), 1, perturbersFile);

		nCm = (nCm > nChebyshev_h[i]) ? nCm : nChebyshev_h[i];

		offset0_h[i] += (3 * Nperturbers + 7);    //add size of header 7*double + Nperturbers * (4 int + double)
		offset1_h[i] += (3 * Nperturbers + 7);    //add size of header

		startTime_h[i] = 100000000.0;     //large number
		endTime_h[i] = 0.0; 
//printf("%d %d %d %d %.20g\n", id_h[i], nChebyshev_h[i], offset0_h[i], offset1_h[i], GM_h[i]);
	}
#else
	for(int i = 0; i < Nperturbers; ++i){
		er = fread(&id_h[i], sizeof(int), 1, perturbersFile);
		er = fread(&nChebyshev_h[i], sizeof(int), 1, perturbersFile);
		er = fread(&offset0_h[i * RKFn], sizeof(int), 1, perturbersFile);
		er = fread(&offset1_h[i * RKFn], sizeof(int), 1, perturbersFile);
		er = fread(&GM_h[i], sizeof(double), 1, perturbersFile);

		nCm = (nCm > nChebyshev_h[i]) ? nCm : nChebyshev_h[i];

		startTime_h[i * RKFn] = 100000000.0;     //large number
		endTime_h[i * RKFn] = 0.0; 

		for(int j = 0; j < RKFn; ++j){
			//copy values from stage 0 to all other stages
			offset0_h[i * RKFn + j] = offset0_h[i * RKFn];
			offset1_h[i * RKFn + j] = offset1_h[i * RKFn];
			startTime_h[i * RKFn + j] = startTime_h[i * RKFn];
			endTime_h[i * RKFn + j] = endTime_h[i * RKFn];
		}

//printf("%d %d %d %d %.20g\n", id_h[i], nChebyshev_h[i], offset0_h[i], offset1_h[i], GM_h[i]);
	}
#endif
	//printf("nCm %d\n", nCm);

	//Find size of entire data file  

	cdata_h = (double*)malloc(Nperturbers * nCm * 3 * sizeof(double));
	snew_h = (double*)malloc(sizeof(double));
#if USEGPU == 1
	datasize = offset1_h[(Nperturbers - 1) * RKFn] - offset0_h[0 * RKFn];
	printf("size of perturbers data table %d\n", datasize);
	data_h = (double*)malloc(datasize * sizeof(double));
#endif
	time = timeStart;
	double c = (CLIGHT / AUtokm) * 86400.0;
	c2 = c * c;
	REAU = RE / AUtokm;   //Earth radius in AU

	xTable_h = (double*)malloc(Nperturbers * sizeof(double));
	yTable_h = (double*)malloc(Nperturbers * sizeof(double));
	zTable_h = (double*)malloc(Nperturbers * sizeof(double));

	vxTable_h = (double*)malloc(Nperturbers * RKFn * sizeof(double));
	vyTable_h = (double*)malloc(Nperturbers * RKFn * sizeof(double));
	vzTable_h = (double*)malloc(Nperturbers * RKFn * sizeof(double));

	x_h = (double*)malloc(N * sizeof(double));
	y_h = (double*)malloc(N * sizeof(double));
	z_h = (double*)malloc(N * sizeof(double));

	vx_h = (double*)malloc(N * sizeof(double));
	vy_h = (double*)malloc(N * sizeof(double));
	vz_h = (double*)malloc(N * sizeof(double));

	xt_h = (double*)malloc(N * sizeof(double));
	yt_h = (double*)malloc(N * sizeof(double));
	zt_h = (double*)malloc(N * sizeof(double));

	vxt_h = (double*)malloc(N * sizeof(double));
	vyt_h = (double*)malloc(N * sizeof(double));
	vzt_h = (double*)malloc(N * sizeof(double));

	dx_h = (double*)malloc(N * sizeof(double));
	dy_h = (double*)malloc(N * sizeof(double));
	dz_h = (double*)malloc(N * sizeof(double));

	dvx_h = (double*)malloc(N * sizeof(double));
	dvy_h = (double*)malloc(N * sizeof(double));
	dvz_h = (double*)malloc(N * sizeof(double));

	kx_h = (double*)malloc(N * RKFn * sizeof(double));
	ky_h = (double*)malloc(N * RKFn * sizeof(double));
	kz_h = (double*)malloc(N * RKFn * sizeof(double));

	kvx_h = (double*)malloc(N * RKFn * sizeof(double));
	kvy_h = (double*)malloc(N * RKFn * sizeof(double));
	kvz_h = (double*)malloc(N * RKFn * sizeof(double));

	ax_h = (double*)malloc(N * sizeof(double));
	ay_h = (double*)malloc(N * sizeof(double));
	az_h = (double*)malloc(N * sizeof(double));

	A1_h = (double*)malloc(N * sizeof(double));
	A2_h = (double*)malloc(N * sizeof(double));
	A3_h = (double*)malloc(N * sizeof(double));


	RKFa_h = (double*)malloc(RKFn * RKFn * sizeof(double));
	RKFb_h = (double*)malloc(RKFn * sizeof(double));
	RKFbb_h = (double*)malloc(RKFn * sizeof(double));
	RKFc_h = (double*)malloc(RKFn * sizeof(double));

	for(int i = 0; i < RKFn; ++i){
		for(int j = 0; j < RKFn; ++j){
			RKFa_h[i * RKFn + j] = 0.0;
		}
		RKFb_h[i] = 0.0;
		RKFbb_h[i] = 0.0;
		RKFc_h[i] = 0.0;
	}


}



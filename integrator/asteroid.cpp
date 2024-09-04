#include "asteroid.h"


int asteroid::readParam(){

	FILE *paramfile;
	paramfile = fopen("param.dat", "r");

	if(paramfile == NULL){
		printf("Error, param.dat file does not exist\n");
		return 0;
	}


	char sp[160];
	int er;

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
		fgets(sp, 3, paramfile);
		}
		else if(strcmp(sp, "Reference Time =") == 0){
			er = fscanf (paramfile, "%lf", &time_reference);
			if(er <= 0){
				printf("Error: Reference Time value is not valid!\n");
				return 0;
			}
			fgets(sp, 3, paramfile);
		}
		else if(strcmp(sp, "Start Time =") == 0){
			er = fscanf (paramfile, "%lf", &timeStart);
			if(er <= 0){
				printf("Error: Start Time value is not valid!\n");
				return 0;
			}
			fgets(sp, 3, paramfile);
		}
		else if(strcmp(sp, "End Time =") == 0){
			er = fscanf (paramfile, "%lf", &timeEnd);
			if(er <= 0){
				printf("Error: End Time value is not valid!\n");
				return 0;
			}
			fgets(sp, 3, paramfile);
		}
		else if(strcmp(sp, "Output Interval =") == 0){
			er = fscanf (paramfile, "%lf", &outputInterval);
			if(er <= 0){
				printf("Error: Output Interval value is not valid!\n");
				return 0;
			}
			fgets(sp, 3, paramfile);
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
			fgets(sp, 3, paramfile);
		}
		else if(strcmp(sp, "RKF absolute tolerance =") == 0){
			er = fscanf (paramfile, "%lf", &RKF_atol);
			if(er <= 0){
				printf("Error: RKF absolute tolerance value is not valid!\n");
				return 0;
			}
			fgets(sp, 3, paramfile);
		}
		else if(strcmp(sp, "RKF relative tolerance =") == 0){
			er = fscanf (paramfile, "%lf", &RKF_rtol);
			if(er <= 0){
				printf("Error: RKF relative tolerance value is not valid!\n");
				return 0;
			}
			fgets(sp, 3, paramfile);
		}
		else if(strcmp(sp, "Time Step =") == 0){
			er = fscanf (paramfile, "%lf", &dt);
			if(er <= 0){
				printf("Error: dt value is not valid!\n");
				return 0;
			}
			fgets(sp, 3, paramfile);
		}

		else{
			printf("Error: param.dat file is not valid! %s\n", sp);
			return 0;
		}
	}
	fclose(paramfile);

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
	for(int i = Nperturbers; i < N; ++i){
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
	N = Nperturbers + 1;

	perturbersFile = fopen("../readChebyshev/PerturbersChebyshev.bin", "rb");
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
	int er;
	er = fread(&time0, sizeof(double), 1, perturbersFile);
	er = fread(&time1, sizeof(double), 1, perturbersFile);
	er = fread(&AUtokm, sizeof(double), 1, perturbersFile);
	er = fread(&EM, sizeof(double), 1, perturbersFile);
	er = fread(&CLIGHT, sizeof(double), 1, perturbersFile); 
	er = fread(&RE, sizeof(double), 1, perturbersFile);
	er = fread(&J2E, sizeof(double), 1, perturbersFile);

	for(int i = 0; i < Nperturbers; ++i){
		er = fread(&id[i], sizeof(int), 1, perturbersFile);
		er = fread(&nChebyshev[i], sizeof(int), 1, perturbersFile);
		er = fread(&offset0[i], sizeof(int), 1, perturbersFile);
		er = fread(&offset1[i], sizeof(int), 1, perturbersFile);
		er = fread(&GM[i], sizeof(double), 1, perturbersFile);

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


	a_h = (double*)malloc(RKFn * RKFn * sizeof(double));
	b_h = (double*)malloc(RKFn * sizeof(double));
	bb_h = (double*)malloc(RKFn * sizeof(double));
	c_h = (double*)malloc(RKFn * sizeof(double));

	for(int i = 0; i < RKFn; ++i){
		for(int j = 0; j < RKFn; ++j){
			a_h[i * RKFn + j] = 0.0;
		}
		b_h[i] = 0.0;
		bb_h[i] = 0.0;
		c_h[i] = 0.0;
	}


}



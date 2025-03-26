#include "asteroid.h"
#include <stdio.h>


int asteroid::readParam(int argc, char*argv[]){

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
		if(strcmp(sp, "Initial condition file format =") == 0){
			char format[160];
			er = fscanf (paramfile, "%s", format);
			if(er <= 0){
				printf("Error: Initial condition file format is not valid!\n");
				return 0;
			}
			if(strcmp(format, "text") == 0){
				ICformat = 0;
			}
			else if(strcmp(format, "binary") == 0){
				ICformat = 1;
			}
			else{
				printf("Error: Initial condition file format is not valid\n");
				return 0;

			}
			str = fgets(sp, 3, paramfile);
		}
		else if(strcmp(sp, "Initial condition file coordinates =") == 0){
			char format[160];
			er = fscanf (paramfile, " %[^\n]s", format);
			if(er <= 0){
				printf("Error: Initial condition file coordinates is not valid!\n");
				return 0;
			}
			if(strcmp(format, "cartesian heliocentric ecliptic") == 0){
				ICorbital = 0;
				ICecliptic = 1;
				ICheliocentric = 1;
			}
			else if(strcmp(format, "cartesian heliocentric equatorial") == 0){
				ICorbital = 0;
				ICecliptic = 0;
				ICheliocentric = 1;
			}
			else if(strcmp(format, "cartesian barycentric ecliptic") == 0){
				ICorbital = 0;
				ICecliptic = 1;
				ICheliocentric = 0;
			}
			else if(strcmp(format, "cartesian barycentric equatorial") == 0){
				ICorbital = 0;
				ICecliptic = 0;
				ICheliocentric = 0;
			}
			else if(strcmp(format, "orbital heliocentric ecliptic") == 0){
				ICorbital = 1;
				ICecliptic = 1;
				ICheliocentric = 1;
			}
			else{
				printf("Error: Initial condition file coordinates is not valid!\n");
				return 0;

			}
			str = fgets(sp, 3, paramfile);
		}
		else if(strcmp(sp, "Initial condition file =") == 0){
			er = fscanf (paramfile, "%s", inputFilename);
			if(er <= 0){
				printf("Error: Initial condition file is not valid!\n");
				return 0;
			}
		str = fgets(sp, 3, paramfile);
		}
		else if(strcmp(sp, "Output file name =") == 0){
			er = fscanf (paramfile, "%s", name);
			if(er <= 0){
				printf("Error: Output file is not valid!\n");
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
		else if(strcmp(sp, "Use GR correction =") == 0){
			er = fscanf (paramfile, "%d", &useGR);
			if(er <= 0){
				printf("Error: Use GR correction is not valid!\n");
				return 0;
		}
			str = fgets(sp, 3, paramfile);
		}
		else if(strcmp(sp, "Use J2 force =") == 0){
			er = fscanf (paramfile, "%d", &useJ2);
			if(er <= 0){
				printf("Error: Use J2 force is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
		}
		else if(strcmp(sp, "Use non-gravitational force =") == 0){
			er = fscanf (paramfile, "%d", &useNonGrav);
			if(er <= 0){
				printf("Error: Use non-gravitational force is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
		}
		else if(strcmp(sp, "Use comets =") == 0){
			er = fscanf (paramfile, "%d", &cometFlag);
			if(er <= 0){
				printf("Error: Use comets is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
		}
		else if(strcmp(sp, "comet alpha =") == 0){
			er = fscanf (paramfile, "%lf", &nonGrav_alpha);
			if(er <= 0){
				printf("Error: comet alpha is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
		}
		else if(strcmp(sp, "comet nk =") == 0){
			er = fscanf (paramfile, "%lf", &nonGrav_nk);
			if(er <= 0){
				printf("Error: comet nk is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
		}
		else if(strcmp(sp, "comet nm =") == 0){
			er = fscanf (paramfile, "%lf", &nonGrav_nm);
			if(er <= 0){
				printf("Error: comet nm is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
		}
		else if(strcmp(sp, "comet nn =") == 0){
			er = fscanf (paramfile, "%lf", &nonGrav_nn);
			if(er <= 0){
				printf("Error: comet nn is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
		}
		else if(strcmp(sp, "comet r0 =") == 0){
			er = fscanf (paramfile, "%lf", &nonGrav_r0);
			if(er <= 0){
				printf("Error: comet r0 is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
		}
		else if(strcmp(sp, "comet tau =") == 0){
			er = fscanf (paramfile, "%lf", &nonGrav_tau);
			if(er <= 0){
				printf("Error: comet tau is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
		}
		else if(strcmp(sp, "Use binary output format =") == 0){
			er = fscanf (paramfile, "%d", &outBinary);
			if(er <= 0){
				printf("Error: outBinary is not valid!\n");
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

	//read console arguments
	for(int i = 1; i < argc; i += 2){
		if(strcmp(argv[i], "-outInterval") == 0){
			outputInterval = atof(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-dt") == 0){
			dt = atof(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-in") == 0){
			sprintf(inputFilename, "%s", argv[i + 1]);
		}
		else if(strcmp(argv[i], "-out") == 0){
			sprintf(name, "%s", argv[i + 1]);
		}
		else if(strcmp(argv[i], "-mode") == 0){
			GPUMode = atoi(argv[i + 1]);
		}
		else{
			printf("Error, console argument is not valid.\n");
		return 0;
		}
	}





	if(outBinary == 0){
		sprintf(outputFilename, "Out_%s.dat", name);
	}
	else{
		sprintf(outputFilename, "Out_%s.bin", name);
	}

	sprintf(infoFilename, "Info_%s.dat", name);



	if(ICformat > 0 && ICformat != 0){
		printf("Error, combination of initial condition format and initial condition coordinate system is not allowed\n");
		return 0;
	}


	return 1;

}

void asteroid::printInfo(){
	fprintf(infoFile, "Initial condition file format = %d\n", ICformat);
	fprintf(infoFile, "Initial condition file = %s\n", inputFilename);
	fprintf(infoFile, "Output file name = %s\n", outputFilename);
	fprintf(infoFile, "Path to perturbers file = %s\n", perturbersFilePath);
	fprintf(infoFile, "Start Time = %.20g\n", timeStart);
	fprintf(infoFile, "End Time = %.20g\n", timeEnd);
	fprintf(infoFile, "Time Step = %.20g\n", dt);
	fprintf(infoFile, "Output Interval = %g\n", outputInterval);
	fprintf(infoFile, "Integrator RKFn = %d\n", RKFn);
	fprintf(infoFile, "RKF absolute tolerance = %.20g\n", RKF_atol);
	fprintf(infoFile, "RKF relative tolerance = %.20g\n", RKF_rtol);

	fprintf(infoFile, "useGR correction = %d\n", useGR);
	fprintf(infoFile, "useJ2 force = %d\n", useJ2);
	fprintf(infoFile, "use non-gravitational force = %d\n", useNonGrav);
	fprintf(infoFile, "Use GPU = %d\n", USEGPU);
	fprintf(infoFile, "Use binary output format = %d\n", outBinary);
	//fprintf(infoFile, "useFIFO = %d\n", useFIFO);
	//fprintf(infoFile, "InVersion = %g\n", InVersion);

	fprintf(infoFile, "Nperturbers = %d\n", Nperturbers);

	//fprintf(infoFile, "outStart: %.20g, time0: %.20g, time1: %.20g, outInterval: %lld\n", outStart, time0, time1, outInterval);
} 




void asteroid::allocate(){
	nCm = 0;
	dts = (dt > 0.0) ? 1.0 : -1.0;      //sign of time step
	dt1 = dt;
	stop = 0;

	idp_h = (int*)malloc(Nperturbers * sizeof(int));
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
		er = fread(&idp_h[i], sizeof(int), 1, perturbersFile);
		er = fread(&nChebyshev_h[i], sizeof(int), 1, perturbersFile);
		er = fread(&offset0_h[i], sizeof(int), 1, perturbersFile);
		er = fread(&offset1_h[i], sizeof(int), 1, perturbersFile);
		er = fread(&GM_h[i], sizeof(double), 1, perturbersFile);

		nCm = (nCm > nChebyshev_h[i]) ? nCm : nChebyshev_h[i];

		offset0_h[i] += (3 * Nperturbers + 7);    //add size of header 7*double + Nperturbers * (4 int + double)
		offset1_h[i] += (3 * Nperturbers + 7);    //add size of header

		startTime_h[i] = 100000000.0;     //large number
		endTime_h[i] = 0.0; 
//printf("%d %d %d %d %.20g\n", idp_h[i], nChebyshev_h[i], offset0_h[i], offset1_h[i], GM_h[i]);
	}
#else
	for(int i = 0; i < Nperturbers; ++i){
		er = fread(&idp_h[i], sizeof(int), 1, perturbersFile);
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

//printf("%d %d %d %d %.20g\n", idp_h[i], nChebyshev_h[i], offset0_h[i], offset1_h[i], GM_h[i]);
	}
#endif
	//printf("nCm %d\n", nCm);

	//Find size of entire data file  

	cdata_h = (double*)malloc(Nperturbers * nCm * 3 * sizeof(double));
#if USEGPU == 1
	datasize = offset1_h[(Nperturbers - 1) * RKFn] - offset0_h[0 * RKFn];
	printf("size of perturbers data table %d\n", datasize);
	data_h = (double*)malloc(datasize * sizeof(double));
#endif
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

	snew_h = (double*)malloc(N * sizeof(double));
	jd_init_h = (double*)malloc(N * sizeof(double));

	id_h = (long long int*)malloc(N * sizeof(long long int));

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

	Rsave_h = (double*)malloc(Rbuffersize * N * sizeof(double));
	Tsave_h = (double*)malloc(Rbuffersize * sizeof(double));

}

void asteroid::output(double dtmin){
	printf("Reached time %.20g dtmin %.8g\n", time_reference + time, dtmin);

	if(time_reference + time >= outStart){
		for(int i = 0; i < N; ++i){
			if(outBinary == 0){
				fprintf(outputFile, "%.20g %lld %.20g %.20g %.20g %.20g %.20g %.20g %.20g\n", time_reference + time, id_h[i], x_h[i], y_h[i], z_h[i], vx_h[i], vy_h[i], vz_h[i], dtmin);
			}
			else{
				long long int id = id_h[i];
				//unsigned long long int id = __builtin_bswap64 (id_h[i]);
				double tt = time_reference + time;
				double xx = x_h[i];
				double yy = y_h[i];
				double zz = z_h[i];
				double vxx = vx_h[i];
				double vyy = vy_h[i];
				double vzz = vz_h[i];

				fwrite(&tt, sizeof(double), 1, outputFile);
				fwrite(&id, sizeof(long long int), 1, outputFile);
				fwrite(&xx, sizeof(double), 1, outputFile);
				fwrite(&yy, sizeof(double), 1, outputFile);
				fwrite(&zz, sizeof(double), 1, outputFile);
				fwrite(&vxx, sizeof(double), 1, outputFile);
				fwrite(&vyy, sizeof(double), 1, outputFile);
				fwrite(&vzz, sizeof(double), 1, outputFile);
				fwrite(&dtmin, sizeof(double), 1, outputFile);


			}
		}
	}
}


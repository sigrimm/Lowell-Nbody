#include "asteroid.h"


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

	int count[1000];	//counter to check if the param.dat file contains dublicated entries
	for(int i = 0; i < 1000; ++i){
		count[i] = 0;
	}

	for(int i = 0; i < nL; ++i){
		dtlimit[i] = 0.0;
	}

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
			if(c == '\n'){
				//blank line
				i = -1;
				continue;
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
			++count[0];
		}
		else if(strcmp(sp, "Initial condition file coordinates =") == 0){
			char format[160];
			er = fscanf (paramfile, " %159[^\n\r]", format);

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
			++count[1];
		}
		else if(strcmp(sp, "Output coordinates =") == 0){
			char format[160];
			er = fscanf (paramfile, " %159[^\n\r]", format);
			if(er <= 0){
				printf("Error: Output coordinates is not valid!\n");
				return 0;
			}
			if(strcmp(format, "cartesian heliocentric ecliptic") == 0){
				Outecliptic = 1;
				Outheliocentric = 1;
			}
			else if(strcmp(format, "cartesian heliocentric equatorial") == 0){
				Outecliptic = 0;
				Outheliocentric = 1;
			}
			else if(strcmp(format, "cartesian barycentric ecliptic") == 0){
				Outecliptic = 1;
			}
			else if(strcmp(format, "cartesian barycentric equatorial") == 0){
				Outecliptic = 0;
			}
			else if(strcmp(format, "cartesian geocentric ecliptic") == 0){
				Outecliptic = 1;
				Outgeocentric = 1;
			}
			else if(strcmp(format, "cartesian geocentric equatorial") == 0){
				Outecliptic = 0;
				Outgeocentric = 1;
			}
			else if(strcmp(format, "orbital heliocentric ecliptic") == 0){
				Outorbital = 1;
				Outecliptic = 1;
				Outheliocentric = 1;
			}
			else{
				printf("Error: Output coordinates is not valid!\n");
				return 0;

			}
			str = fgets(sp, 3, paramfile);
			++count[2];
		}
		else if(strcmp(sp, "Initial condition file =") == 0){
			er = fscanf (paramfile, "%s", inputFilename);
			if(er <= 0){
				printf("Error: Initial condition file is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
			++count[3];
		}
		else if(strcmp(sp, "Output file name =") == 0){
			er = fscanf (paramfile, "%s", name);
			if(er <= 0){
				printf("Error: Output file is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
			++count[4];
		}
		else if(strcmp(sp, "Reference Time =") == 0){
			er = fscanf (paramfile, "%lf", &time_reference);
			if(er <= 0){
				printf("Error: Reference Time value is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
			++count[5];
		}
		else if(strcmp(sp, "Start Time =") == 0){
			er = fscanf (paramfile, "%lf", &timeStart);
			if(er <= 0){
				printf("Error: Start Time value is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
			++count[6];
		}
		else if(strcmp(sp, "End Time =") == 0){
			er = fscanf (paramfile, "%lf", &timeEnd);
			if(er <= 0){
				printf("Error: End Time value is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
			++count[7];
		}
		else if(strcmp(sp, "Output Interval =") == 0){
			er = fscanf (paramfile, "%lf", &outputInterval);
			if(er <= 0){
				printf("Error: Output Interval value is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
			++count[8];
		}
		else if(strcmp(sp, "Integrator =") == 0){
			er = fscanf (paramfile, "%s", integratorName);
			if(er <= 0){
				printf("Error: Integrator Name is not valid!\n");
				return 0;
			}
			if(strcmp(integratorName, "LF") == 0){
				nStage = 1;
			}
			else if(strcmp(integratorName, "RK4") == 0){
				RKFn = 4;
				nStage = 4;
			}
			else if(strcmp(integratorName, "RK7") == 0){
				RKFn = 9;
				nStage = 9;
			}
			else if(strcmp(integratorName, "RKF45") == 0){
				RKFn = 6;
				nStage = 6;
			}
			else if(strcmp(integratorName, "DP54") == 0){
				RKFn = 7;
				nStage = 7;
			}
			else if(strcmp(integratorName, "RKF78") == 0){
				RKFn = 13;
				nStage = 13;
			}
			else if(strcmp(integratorName, "BS") == 0){
				BSn = 80;
				nStage = 80;
			}
			else if(strcmp(integratorName, "IMM") == 0){
				nStage = 1;
			}
			else{
				printf("Error, Integrator not valid\n");
				return 0;

			}
			str = fgets(sp, 3, paramfile);
			++count[9];
		}
		else if(strcmp(sp, "Integrator absolute tolerance =") == 0){
			er = fscanf (paramfile, "%lf", &atol);
			if(er <= 0){
				printf("Error: Integrator absolute tolerance value is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
			++count[10];
		}
		else if(strcmp(sp, "Integrator relative tolerance =") == 0){
			er = fscanf (paramfile, "%lf", &rtol);
			if(er <= 0){
				printf("Error: Integrator relative tolerance value is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
			++count[11];
		}
		else if(strcmp(sp, "Time Step =") == 0){
			er = fscanf (paramfile, "%lf", &dt);
			if(er <= 0){
				printf("Error: Time Step value is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
			++count[12];
		}
		else if(strcmp(sp, "Time Step Levels =") == 0){
			for(int i = 0; i < def_nLMax; ++i){
				//store file position
				long pos = ftell(paramfile); 

				er = fscanf (paramfile, "%lf", &dtlimit[i]);
				if(er <= 0){
					//set back read position
					fseek(paramfile, pos, SEEK_SET);
					nL = i + 1;
					break;
				}
				dtlimit[i] = abs(dtlimit[i]);
				if(i == def_nLMax - 1){
					printf("More levels than available %d, def_nLMax =  %d\n", i + 1, def_nLMax);
					return 0;
				}
			}
			str = fgets(sp, 3, paramfile);
			++count[13];
		}
		else if(strcmp(sp, "Path to perturbers file =") == 0){
			er = fscanf (paramfile, "%s", perturbersFilePath);
			if(er <= 0){
				printf("Error: Path to perturbers file is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
			++count[14];
		}
		else if(strcmp(sp, "Use GR correction =") == 0){
			er = fscanf (paramfile, "%d", &useGR);
			if(er <= 0){
				printf("Error: Use GR correction is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
			++count[15];
		}
		else if(strcmp(sp, "Use J2 force =") == 0){
			er = fscanf (paramfile, "%d", &useJ2);
			if(er <= 0){
				printf("Error: Use J2 force is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
			++count[16];
		}
		else if(strcmp(sp, "Use non-gravitational force =") == 0){
			er = fscanf (paramfile, "%d", &useNonGrav);
			if(er <= 0){
				printf("Error: Use non-gravitational force is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
			++count[17];
		}
		else if(strcmp(sp, "Use comets =") == 0){
			er = fscanf (paramfile, "%d", &cometFlag);
			if(er <= 0){
				printf("Error: Use comets is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
			++count[18];
		}
		else if(strcmp(sp, "comet alpha =") == 0){
			er = fscanf (paramfile, "%lf", &nonGrav_alpha);
			if(er <= 0){
				printf("Error: comet alpha is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
			++count[19];
		}
		else if(strcmp(sp, "comet nk =") == 0){
			er = fscanf (paramfile, "%lf", &nonGrav_nk);
			if(er <= 0){
				printf("Error: comet nk is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
			++count[20];
		}
		else if(strcmp(sp, "comet nm =") == 0){
			er = fscanf (paramfile, "%lf", &nonGrav_nm);
			if(er <= 0){
				printf("Error: comet nm is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
			++count[21];
		}
		else if(strcmp(sp, "comet nn =") == 0){
			er = fscanf (paramfile, "%lf", &nonGrav_nn);
			if(er <= 0){
				printf("Error: comet nn is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
			++count[22];
		}
		else if(strcmp(sp, "comet r0 =") == 0){
			er = fscanf (paramfile, "%lf", &nonGrav_r0);
			if(er <= 0){
				printf("Error: comet r0 is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
			++count[23];
		}
		else if(strcmp(sp, "comet tau =") == 0){
			er = fscanf (paramfile, "%lf", &nonGrav_tau);
			if(er <= 0){
				printf("Error: comet tau is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
			++count[24];
		}
		else if(strcmp(sp, "Use binary output format =") == 0){
			er = fscanf (paramfile, "%d", &outBinary);
			if(er <= 0){
				printf("Error: outBinary is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
			++count[25];
		}
		else if(strcmp(sp, "Print time steps =") == 0){
			er = fscanf (paramfile, "%d", &printdt);
			if(er <= 0){
				printf("Error: Print time steps is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
			++count[26];
		}
		else if(strcmp(sp, "Use individual time step mode =") == 0){
			er = fscanf (paramfile, "%d", &individualSteps);
			if(er <= 0){
				printf("Error: Time step mode is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
			++count[27];
		}
		else if(strcmp(sp, "GPU block mode =") == 0){
			er = fscanf (paramfile, "%d", &GPUMode);
			if(er <= 0){
				printf("Error: GPU block mode is not valid!\n");
				return 0;
			}
			str = fgets(sp, 3, paramfile);
			++count[28];
		}

		else{
			printf("Error: param.dat file is not valid! %s\n", sp);
			return 0;
		}
	}
	fclose(paramfile);

	for(int i = 0; i < 1000; ++i){
		if(count[i] > 1){
			printf("Error, param.dat file contains multiple entries\n");
			return 0;
		}
	}

	// ----------------------------------
	//read console arguments
	// ----------------------------------
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
		else if(strcmp(argv[i], "-GPUmode") == 0){
			GPUMode = atoi(argv[i + 1]);
		}
		else{
			printf("Error, console argument is not valid.\n");
		return 0;
		}
	}
	// ----------------------------------




	if(outBinary == 0){
		sprintf(outputFilename, "Out_%s.dat", name);
	}
	else{
		sprintf(outputFilename, "Out_%s.bin", name);
	}

	sprintf(infoFilename, "Info_%s.dat", name);
	sprintf(dtFilename, "dt_%s.dat", name);



	for(int i = 1; i < nL; ++i){
		if(dtlimit[i] >= dtlimit[i - 1]){
			printf("dtlimits not in order %g %g\n", dtlimit[i], dtlimit[i - 1]);
			return 0;
		}
	}



	if(strcmp(integratorName, "LF") == 0){
		individualSteps = 0;
	}
	else if(strcmp(integratorName, "RK4") == 0){
		individualSteps = 0;
	}
	else if(strcmp(integratorName, "RK7") == 0){
		individualSteps = 0;
	}
	else if(strcmp(integratorName, "IMM") == 0){
		individualSteps = 0;
	}


	if(individualSteps == 1){
		nL = 1;
	}


	return 1;

}

void asteroid::printInfo(){
	fprintf(infoFile, "Code version = %g\n", def_version);
	fprintf(infoFile, "Use GPU = %d\n", USEGPU);
	fprintf(infoFile, "GPU block mode = %d\n", GPUMode);
	fprintf(infoFile, "Use individual time steps = %d\n", individualSteps);

	fprintf(infoFile, "Initial condition file format = %d\n", ICformat);
	fprintf(infoFile, "Initial condition file coordinates: orbital = %d\n", ICorbital);
	fprintf(infoFile, "Initial condition file coordinates: ecliptic = %d\n", ICecliptic);
	fprintf(infoFile, "Initial condition file coordinates: heliocentric= %d\n", ICheliocentric);
	fprintf(infoFile, "Initial condition file = %s\n", inputFilename);
	fprintf(infoFile, "Output file name = %s\n", outputFilename);
	fprintf(infoFile, "Output Interval = %g\n", outputInterval);
	fprintf(infoFile, "Output coordinates: orbital = %d\n", Outorbital);
	fprintf(infoFile, "Output coordinates: ecliptic = %d\n", Outecliptic);
	fprintf(infoFile, "Out coordinates: heliocentric= %d\n", Outheliocentric);
	fprintf(infoFile, "Use binary output format = %d\n", outBinary);
	fprintf(infoFile, "Path to perturbers file = %s\n", perturbersFilePath);
	fprintf(infoFile, "Start Time = %.20g\n", timeStart);
	fprintf(infoFile, "End Time = %.20g\n", timeEnd);
	fprintf(infoFile, "Time Step = %.20g\n", dt);
	fprintf(infoFile, "Number of Time Step levels = %d\n", nL);
	for(int i = 0; i < nL - 1; ++i){
		fprintf(infoFile, "Time Step level %d = %.20g\n", i, dtlimit[i]);
	}
	fprintf(infoFile, "Integrator name = %s\n", integratorName);
	fprintf(infoFile, "Integrator absolute tolerance = %.20g\n", atol);
	fprintf(infoFile, "Integrator relative tolerance = %.20g\n", rtol);

	fprintf(infoFile, "Use GR correction = %d\n", useGR);
	fprintf(infoFile, "Use J2 force = %d\n", useJ2);
	fprintf(infoFile, "Use non-gravitational force = %d\n", useNonGrav);
	fprintf(infoFile, "Use comets = %d\n", cometFlag);
	fprintf(infoFile, "comet alpha = %g\n", nonGrav_alpha);
	fprintf(infoFile, "comet nk = %g\n", nonGrav_nk);
	fprintf(infoFile, "comet nm = %g\n", nonGrav_nm);
	fprintf(infoFile, "comet nn = %g\n", nonGrav_nn);
	fprintf(infoFile, "comet r0 = %g\n", nonGrav_r0);
	fprintf(infoFile, "comet tau = %g\n", nonGrav_tau);
	fprintf(infoFile, "Print time steps = %d\n", printdt);
	//fprintf(infoFile, "useFIFO = %d\n", useFIFO);
	//fprintf(infoFile, "InVersion = %g\n", InVersion);

	fprintf(infoFile, "Nperturbers = %d\n", Nperturbers);

	//fprintf(infoFile, "outStart: %.20g, time0: %.20g, time1: %.20g, outInterval: %lld\n", outStart, time0, time1, outInterval);
} 




int asteroid::allocate(){
	nCm = 0;
	dts = (dt > 0.0) ? 1.0 : -1.0;      //sign of time step
	dtsave = dt;
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
	startTime_h = (double*)malloc(Nperturbers * nStage * sizeof(double));
	endTime_h = (double*)malloc(Nperturbers * nStage * sizeof(double));
	offset0_h = (int*)malloc(Nperturbers * nStage * sizeof(int));
	offset1_h = (int*)malloc(Nperturbers * nStage * sizeof(int));
#endif
	GM_h = (double*)malloc(Nperturbers * sizeof(double));

	//read header
	int er;
	printf("  Start reading header of perturbers file\n");

	er = fread(&time0, sizeof(double), 1, perturbersFile);
	if(er <= 0){
		return 0;
	}

	er = fread(&time1, sizeof(double), 1, perturbersFile);
	if(er <= 0){
		return 0;
	}

	er = fread(&AUtokm, sizeof(double), 1, perturbersFile);
	if(er <= 0){
		return 0;
	}

	er = fread(&EM, sizeof(double), 1, perturbersFile);
	if(er <= 0){
		return 0;
	}

	er = fread(&CLIGHT, sizeof(double), 1, perturbersFile); 
	if(er <= 0){
		return 0;
	}

	er = fread(&RE, sizeof(double), 1, perturbersFile);
	if(er <= 0){
		return 0;
	}

	er = fread(&J2E, sizeof(double), 1, perturbersFile);
	if(er <= 0){
		return 0;
	}

	printf("  Reading header of perturbers file OK\n");


	//Check if perturbers file time span is OK
	if(timeStart > time1 || timeEnd > time1 || timeStart < time0 || timeEnd < time0){
		printf("Error, the integration time is longer than the data in the perturbers file\n");
		printf("Integration time: %.20g %.20g\n", timeStart, timeEnd);
		printf("Time in perturbers file: %.20g %.20g\n", time0, time1);
		return 0;

	}
	


#if USEGPU == 0
	for(int i = 0; i < Nperturbers; ++i){
		er = fread(&idp_h[i], sizeof(int), 1, perturbersFile);
		if(er <= 0){
			printf("Error, reading perturber %d failed \n", i);
			return 0;
		}

		er = fread(&nChebyshev_h[i], sizeof(int), 1, perturbersFile);
		if(er <= 0){
			printf("Error, reading perturber %d failed \n", i);
			return 0;
		}

		er = fread(&offset0_h[i], sizeof(int), 1, perturbersFile);
		if(er <= 0){
			printf("Error, reading perturber %d failed \n", i);
			return 0;
		}

		er = fread(&offset1_h[i], sizeof(int), 1, perturbersFile);
		if(er <= 0){
			printf("Error, reading perturber %d failed \n", i);
			return 0;
		}

		er = fread(&GM_h[i], sizeof(double), 1, perturbersFile);
		if(er <= 0){
			printf("Error, reading perturber %d failed \n", i);
			return 0;
		}

		nCm = (nCm > nChebyshev_h[i]) ? nCm : nChebyshev_h[i];

		offset0_h[i] += (3 * Nperturbers + 7);    //add size of header 7*double + Nperturbers * (4 int + double)
		offset1_h[i] += (3 * Nperturbers + 7);    //add size of header

		startTime_h[i] = 100000000.0;     //large number
		endTime_h[i] = 0.0; 
//printf(" %d %d %d %d %d %.20g\n", i, idp_h[i], nChebyshev_h[i], offset0_h[i], offset1_h[i], GM_h[i]);
	}
#else
	for(int i = 0; i < Nperturbers; ++i){
		er = fread(&idp_h[i], sizeof(int), 1, perturbersFile);
		if(er <= 0){
			printf("Error, reading perturber %d failed \n", i);
			return 0;
		}

		er = fread(&nChebyshev_h[i], sizeof(int), 1, perturbersFile);
		if(er <= 0){
			printf("Error, reading perturber %d failed \n", i);
			return 0;
		}

		er = fread(&offset0_h[i * nStage], sizeof(int), 1, perturbersFile);
		if(er <= 0){
			printf("Error, reading perturber %d failed \n", i);
			return 0;
		}

		er = fread(&offset1_h[i * nStage], sizeof(int), 1, perturbersFile);
		if(er <= 0){
			printf("Error, reading perturber %d failed \n", i);
			return 0;
		}

		er = fread(&GM_h[i], sizeof(double), 1, perturbersFile);
		if(er <= 0){
			printf("Error, reading perturber %d failed \n", i);
			return 0;
		}

		nCm = (nCm > nChebyshev_h[i]) ? nCm : nChebyshev_h[i];

		startTime_h[i * nStage] = 100000000.0;     //large number
		endTime_h[i * nStage] = 0.0; 

		for(int j = 0; j < nStage; ++j){
			//copy values from stage 0 to all other stages
			offset0_h[i * nStage + j] = offset0_h[i * nStage];
			offset1_h[i * nStage + j] = offset1_h[i * nStage];
			startTime_h[i * nStage + j] = startTime_h[i * nStage];
			endTime_h[i * nStage + j] = endTime_h[i * nStage];
		}

//printf("%d %d %d %d %d %.20g\n", i, idp_h[i], nChebyshev_h[i], offset0_h[i], offset1_h[i], GM_h[i]);
	}
#endif
	printf("nCm %d\n", nCm);
	if(nCm > def_nCMax){
		printf("Error: nCm larger than def_nCMax: %d %d\n", nCm, def_nCMax);
	}

	//Find size of entire data file  

	cdata_h = (double*)malloc(Nperturbers * nCm * 3 * sizeof(double));
#if USEGPU == 1
	datasize = offset1_h[(Nperturbers - 1) * nStage] - offset0_h[0 * nStage];
	printf("size of perturbers data table %d\n", datasize);
	data_h = (double*)malloc(datasize * sizeof(double));
#endif
	double c = (CLIGHT / AUtokm) * 86400.0;
	c2 = c * c;
	REAU = RE / AUtokm;   //Earth radius in AU

	xTable_h = (double*)malloc(Nperturbers * sizeof(double));
	yTable_h = (double*)malloc(Nperturbers * sizeof(double));
	zTable_h = (double*)malloc(Nperturbers * sizeof(double));

	vxTable_h = (double*)malloc(Nperturbers * sizeof(double));
	vyTable_h = (double*)malloc(Nperturbers * sizeof(double));
	vzTable_h = (double*)malloc(Nperturbers * sizeof(double));

	x_h = (double*)malloc(N * sizeof(double));
	y_h = (double*)malloc(N * sizeof(double));
	z_h = (double*)malloc(N * sizeof(double));

	vx_h = (double*)malloc(N * sizeof(double));
	vy_h = (double*)malloc(N * sizeof(double));
	vz_h = (double*)malloc(N * sizeof(double));

	xout_h = (double*)malloc(N * sizeof(double));
	yout_h = (double*)malloc(N * sizeof(double));
	zout_h = (double*)malloc(N * sizeof(double));

	vxout_h = (double*)malloc(N * sizeof(double));
	vyout_h = (double*)malloc(N * sizeof(double));
	vzout_h = (double*)malloc(N * sizeof(double));

	xt_h = (double*)malloc(N * sizeof(double));
	yt_h = (double*)malloc(N * sizeof(double));
	zt_h = (double*)malloc(N * sizeof(double));

	vxt_h = (double*)malloc(N * sizeof(double));
	vyt_h = (double*)malloc(N * sizeof(double));
	vzt_h = (double*)malloc(N * sizeof(double));

	xp_h = (double*)malloc(N * sizeof(double));
	yp_h = (double*)malloc(N * sizeof(double));
	zp_h = (double*)malloc(N * sizeof(double));

	vxp_h = (double*)malloc(N * sizeof(double));
	vyp_h = (double*)malloc(N * sizeof(double));
	vzp_h = (double*)malloc(N * sizeof(double));

	scalex_h = (double*)malloc(N * sizeof(double));
	scaley_h = (double*)malloc(N * sizeof(double));
	scalez_h = (double*)malloc(N * sizeof(double));

	scalevx_h = (double*)malloc(N * sizeof(double));
	scalevy_h = (double*)malloc(N * sizeof(double));
	scalevz_h = (double*)malloc(N * sizeof(double));


	if(strcmp(integratorName, "BS") == 0){
		dx_h = (double*)malloc(N * 8 * sizeof(double));
		dy_h = (double*)malloc(N * 8 * sizeof(double));
		dz_h = (double*)malloc(N * 8 * sizeof(double));

		dvx_h = (double*)malloc(N * 8 * sizeof(double));
		dvy_h = (double*)malloc(N * 8 * sizeof(double));
		dvz_h = (double*)malloc(N * 8 * sizeof(double));
	}
	else{
		dx_h = (double*)malloc(N * sizeof(double));
		dy_h = (double*)malloc(N * sizeof(double));
		dz_h = (double*)malloc(N * sizeof(double));

		dvx_h = (double*)malloc(N * sizeof(double));
		dvy_h = (double*)malloc(N * sizeof(double));
		dvz_h = (double*)malloc(N * sizeof(double));
	}


	A1_h = (double*)malloc(N * sizeof(double));
	A2_h = (double*)malloc(N * sizeof(double));
	A3_h = (double*)malloc(N * sizeof(double));

	dt_h = (double*)malloc(N * sizeof(double));
	dtsave_h = (double*)malloc(N * sizeof(double));
	dtmin_h = (double*)malloc(N * sizeof(double));
	dtminlevel_h = (double*)malloc(nL * sizeof(double));
	time_h = (double*)malloc(N * sizeof(double));
	timeStep_h = (long long int*)malloc(N * sizeof(long long int));
	snew_h = (double*)malloc(N * sizeof(double));
	snewlevel_h = (double*)malloc(nL * sizeof(double));
	jd_init_h = (double*)malloc(N * sizeof(double));
	Nlevel_h = (int*)malloc(nL * sizeof(int));
	stop_h = (int*)malloc(nL * sizeof(int));


	id_h = (long long int*)malloc(N * sizeof(long long int));
	index_h = (int*)malloc(N * nL * sizeof(int));

	if(RKFn > 0){

		kx_h = (double*)malloc(N * RKFn * sizeof(double));
		ky_h = (double*)malloc(N * RKFn * sizeof(double));
		kz_h = (double*)malloc(N * RKFn * sizeof(double));

		kvx_h = (double*)malloc(N * RKFn * sizeof(double));
		kvy_h = (double*)malloc(N * RKFn * sizeof(double));
		kvz_h = (double*)malloc(N * RKFn * sizeof(double));

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
	else{
		kx_h = NULL;
		ky_h = NULL;
		kz_h = NULL;

		kvx_h = NULL;
		kvy_h = NULL;
		kvz_h = NULL;

		RKFa_h = NULL;
		RKFb_h = NULL;
		RKFbb_h = NULL;
		RKFc_h = NULL;
	}


	if(BSn > 0){
		BSddt_h = (double*)malloc(8 * sizeof(double));
		BSt0_h = (double*)malloc(8 * 8 * sizeof(double));
		BSc_h = (double*)malloc(BSn * sizeof(double));


		for(int n = 1; n <= 8; ++n){
			BSddt_h[n-1] = 0.25 / (n*n);
		}

		for(int n = 1; n <= 8; ++n){
			for(int j = n-1; j >=1; --j){
				BSt0_h[(n-1) * 8 + (j - 1)] = 1.0 / (BSddt_h[j-1] - BSddt_h[n-1]);
			}
		}


		for(int i = 0; i < BSn; ++i){
			BSc_h[i] = 0.0;
		}

		int cc = 0;	
		for(int n = 1; n <= 8; ++n){

			double dt2 = 1.0 / (2.0 * n);
			double dt22 = dt2 * 2.0;

			BSc_h[cc] = 0.0;
//printf("%d %d %d %g\n", n, 0, cc, BSc_h[cc]);
			++cc;

			BSc_h[cc] = dt2;
//printf("%d %d %d %g\n", n, 0, cc, BSc_h[cc]);
			++cc;

			for(int m = 2; m <= n; ++m){
				BSc_h[cc] = (m-1) * dt22;
//printf("%d %d %d %g\n", n, m, cc, BSc_h[cc]);
				++cc;

				BSc_h[cc] = (m-1) * dt22 + dt2;
//printf("%d %d %d %g\n", n, m, cc, BSc_h[cc]);
				++cc;
			}

			BSc_h[cc] = 1.0;
//printf("%d %d %d %g\n", n, n, cc, BSc_h[cc]);
			++cc;
		}

	}
	else{
		BSddt_h = NULL;
		BSt0_h = NULL;
		BSc_h = NULL;
	}



	if(cometFlag > 0){
		Rsave_h = (double*)malloc(Rbuffersize * N * sizeof(double));
		Tsave_h = (double*)malloc(Rbuffersize * sizeof(double));
		if(Rsave_h == NULL){
			printf("Error, not enough memory for Rsave\n");
			return 0;
		}
	}
	else{
		Rsave_h = NULL;
		Tsave_h = NULL;
	}

	for(int i = 0; i < N; ++i){
		snew_h[i] = 0;
	}

	return 1;
}


void asteroid::printOutput(){
	printf("Reached time %.20g dtmin: ", time_reference + time);

	for(int i = 0; i < nL; ++i){
		printf("%.8g ", dtminlevel_h[i]);
	}
	printf("\n");

	for(int i = 0; i < N; ++i){
		if(outBinary == 0){
			fprintf(outputFile, "%.20g %lld %.20g %.20g %.20g %.20g %.20g %.20g %.20g\n", time_reference + time, id_h[i], xout_h[i], yout_h[i], zout_h[i], vxout_h[i], vyout_h[i], vzout_h[i], dtmin_h[i]);


//fprintf(outputFile, "%.20g %d %.20g %.20g %.20g %.20g %.20g %.20g %.20g\n", time_reference + time, 2, xTable_h[2], yTable_h[2], zTable_h[2], vxTable_h[2], vyTable_h[2], vzTable_h[2], dtmin_h[i]);
//fprintf(outputFile, "%.20g %d %.20g %.20g %.20g %.20g %.20g %.20g %.20g\n", time_reference + time, 9, xTable_h[9], yTable_h[9], zTable_h[9], vxTable_h[9], vyTable_h[9], vzTable_h[9], dtmin_h[i]);
		}
		else{
			long long int id_ = id_h[i];
			//unsigned long long int id = __builtin_bswap64 (id_h[i]);
			double t_ = time_reference + time;
			double x_ = xout_h[i];
			double y_ = yout_h[i];
			double z_ = zout_h[i];
			double vx_ = vxout_h[i];
			double vy_ = vyout_h[i];
			double vz_ = vzout_h[i];
			double dtmin_ = dtmin_h[i];

			fwrite(&t_, sizeof(double), 1, outputFile);
			fwrite(&id_, sizeof(long long int), 1, outputFile);
			fwrite(&x_, sizeof(double), 1, outputFile);
			fwrite(&y_, sizeof(double), 1, outputFile);
			fwrite(&z_, sizeof(double), 1, outputFile);
			fwrite(&vx_, sizeof(double), 1, outputFile);
			fwrite(&vy_, sizeof(double), 1, outputFile);
			fwrite(&vz_, sizeof(double), 1, outputFile);
			fwrite(&dtmin_, sizeof(double), 1, outputFile);


		}
	}
}

void asteroid::freeMemory(){

	free(idp_h);
	free(nChebyshev_h);
	free(startTime_h);
	free(endTime_h);
	free(offset0_h);
	free(offset1_h);
	free(GM_h);
	
	free(cdata_h);

#if USEGPU == 1

	free(data_h);
#endif

	free(xTable_h);
	free(yTable_h);
	free(zTable_h);

	free(vxTable_h);
	free(vyTable_h);
	free(vzTable_h);

	free(x_h);
	free(y_h);
	free(z_h);

	free(vx_h);
	free(vy_h);
	free(vz_h);

	free(xout_h);
	free(yout_h);
	free(zout_h);

	free(vxout_h);
	free(vyout_h);
	free(vzout_h);

	free(xt_h);
	free(yt_h);
	free(zt_h);

	free(vxt_h);
	free(vyt_h);
	free(vzt_h);

	free(xp_h);
	free(yp_h);
	free(zp_h);

	free(vxp_h);
	free(vyp_h);
	free(vzp_h);


	free(scalex_h);
	free(scaley_h);
	free(scalez_h);

	free(scalevx_h);
	free(scalevy_h);
	free(scalevz_h);

	free(dx_h);
	free(dy_h);
	free(dz_h);

	free(dvx_h);
	free(dvy_h);
	free(dvz_h);

	free(A1_h);
	free(A2_h);
	free(A3_h);

	free(dt_h);
	free(dtsave_h);
	free(dtmin_h);
	free(dtminlevel_h);
	free(time_h);
	free(timeStep_h);

	free(snew_h);
	free(snewlevel_h);
	free(jd_init_h);
	free(id_h);
	free(index_h);
	free(Nlevel_h);
	free(stop_h);

	if(RKFn > 0){
		free(kx_h);
		free(ky_h);
		free(kz_h);

		free(kvx_h);
		free(kvy_h);
		free(kvz_h);

		free(RKFa_h);
		free(RKFb_h);
		free(RKFbb_h);
		free(RKFc_h);
	}

	if(BSn > 0){
		free(BSddt_h);
		free(BSt0_h);
		free(BSc_h);
	}

	if(cometFlag > 0){
		free(Rsave_h);
		free(Tsave_h);

	}


}




#include "asteroid.h"

// *************************************************************
//Initial conditions in text format
// *************************************************************
//Read the size of the initial conditions file in text format
int asteroid::readICSize(){


	N = 0;

	long long int id;
	double x, y, z;
	double vx, vy, vz;
	double A1, A2, A3;

	for(int i = 0; i < 1024 * 1024; ++i){
		int er = 0;
		er = fscanf(inputFile, "%lld", &id);
		er = fscanf(inputFile, "%lf", &x);
		er = fscanf(inputFile, "%lf", &y);
		er = fscanf(inputFile, "%lf", &z);
		er = fscanf(inputFile, "%lf", &vx);
		er = fscanf(inputFile, "%lf", &vy);
		er = fscanf(inputFile, "%lf", &vz);
		er = fscanf(inputFile, "%lf", &A1);
		er = fscanf(inputFile, "%lf", &A2);
		er = fscanf(inputFile, "%lf", &A3);
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

	return 1;

}

//Read the initial conditions file in text format
int asteroid::readIC(){

	for(int i = 0; i < N; ++i){
		int er = 0;
		er = fscanf(inputFile, "%lld", &id_h[i]);
		er = fscanf(inputFile, "%lf", &x_h[i]);
		er = fscanf(inputFile, "%lf", &y_h[i]);
		er = fscanf(inputFile, "%lf", &z_h[i]);
		er = fscanf(inputFile, "%lf", &vx_h[i]);
		er = fscanf(inputFile, "%lf", &vy_h[i]);
		er = fscanf(inputFile, "%lf", &vz_h[i]);
		er = fscanf(inputFile, "%lf", &A1_h[i]);
		er = fscanf(inputFile, "%lf", &A2_h[i]);
		er = fscanf(inputFile, "%lf", &A3_h[i]);
		if(er < 0){
			printf("Error, reading initial conditions file failed.\n");
			return 0;
		}
	//printf("xyz %.40g %.40g %.40g %.40g %.40g %.40g %.40g %.40g %.40g\n", x_h[i], y_h[i], z_h[i], vx_h[i], vy_h[i], vz_h[i], A1_h[i], A2_h[i], A3_h[i]);
	}

	return 1;

}
// *************************************************************


// *************************************************************
//Initial conditions in binary format
// *************************************************************
//Read the header
//The binary file contains 10 8byte entries.
//The header contains:
//0x8100::0000,	uint8
//starting JD,	float8
//ending JD,	float8
//JD step(days)	float8
//comet flag	float8
//-		float8	
//-		float8	
//-		float8	
//-		float8	


//Header Record
//0-0 header byte constant 0x81 Identifies this record as a “header”
//1-7 bytes constant 0x00 0x00 0x00 0x00 0x00 0x00 0x00
//8-15 jd_start float8 Starting julian date for integrated output position, velocity vectors
//16-23 jd_end float8 Ending julian date for integrated output position, velocity vectors
//24-31 epoch days float8 Number of days between outputs (negative indicates backward in time)
//32-39 object type float8 0.0 = asteroid, 1.0 = comet
//40-47 number of objects float8 Number of detail records in batch
//48-55 jd_init float8 epoch of initial conditions
//56-63 Format Version float or (unused) float8 0.0
//64-71 (unused) float8 0.0
//72-79 (unused) float8 0.0
//80-87 (unused) float8 0.0 (in Version > 1.0)
int asteroid::readHeader(){

	double temp;
	long long int header;
	int er;

	//read header:
	er = fread(&header, sizeof(unsigned long long int), 1, inputFile);

	printf("header %lld %d\n", header, er);
	if(header != 129 || er <= 0.0){
		printf("Error in reading file header\n");
		return 0;
	}

	er = fread(&outStart, sizeof(double), 1, inputFile);
	er = fread(&timeEnd, sizeof(double), 1, inputFile);
	er = fread(&outputInterval, sizeof(double), 1, inputFile);


	er = fread(&temp, sizeof(double), 1, inputFile);
	cometFlag = int(temp);
	printf("comet flag %d\n", cometFlag);

	er = fread(&temp, sizeof(double), 1, inputFile);
	N = int(temp);
	printf("N %d\n", N);
	if(N >= 1024 * 1024 - 1){
		printf("Error, N is too large for scan kernels\n");
		return 0;
	}

	er = fread(&timeStart, sizeof(double), 1, inputFile);
	printf("epoch (start time), used in Version >= 1: %.20g\n", timeStart);


	er = fread(&inputFileVersion, sizeof(double), 1, inputFile);
	printf("File format version number: %g\n", inputFileVersion);

	er = fread(&temp, sizeof(double), 1, inputFile);
	er = fread(&temp, sizeof(double), 1, inputFile);


	printf("Start time: %.20g\n", timeStart);
	printf("End time: %.20g\n", timeEnd);
	printf("Output Interval: %.20g\n", outputInterval);
	printf("Output start: %.20g\n", outStart);

	if(inputFileVersion > 1.0){
		er = fread(&temp, sizeof(double), 1, inputFile);
	}

	return 1;

}

// read initial condition file 
//Detail Record
//0-0 detail byte constant 0x01 Identifies this record as a “detail”
//0-7 id 8 bytes or Big-endian unsigned long Object identifier
//8-15 x float8 Initial x position
//16-23 y float8 initial y position
//24-31 z float8 initial z position
//32-39 vx float8 initial x velocityinputFileVersion
//40-47 vy float8 initial y velocity
//48-55 vz float8 initial z velocity
//56-63 a1 float8 a1 non-grav term
//64-71 a2 float8 a2 non-grav term
//72-79 a3 float8 a3 non-grav term
//Trailer Record
//0-0 trailer byte constant 0x82 identifies this record as a “trailer”
//1-79 (padding) bytes constant 0x00 (all zero bytes)
int asteroid::readFile(){

	int er = 0;

	if(inputFileVersion > 1.0){
		//search for smallest starting time
		timeStart = 1.0e10;
	}
	for(int i = 0; i < N; ++i){
		er = fread(&id_h[i], sizeof(unsigned long long int), 1, inputFile);
		er = fread(&x_h[i], sizeof(double), 1, inputFile);
		er = fread(&y_h[i], sizeof(double), 1, inputFile);
		er = fread(&z_h[i], sizeof(double), 1, inputFile);
		er = fread(&vx_h[i], sizeof(double), 1, inputFile);
		er = fread(&vy_h[i], sizeof(double), 1, inputFile);
		er = fread(&vz_h[i], sizeof(double), 1, inputFile);
		er = fread(&A1_h[i], sizeof(double), 1, inputFile);
		er = fread(&A2_h[i], sizeof(double), 1, inputFile);
		er = fread(&A3_h[i], sizeof(double), 1, inputFile);
		if(inputFileVersion > 1.0){
			er = fread(&jd_init_h[i], sizeof(double), 1, inputFile);
			timeStart = (timeStart < jd_init_h[i]) ? timeStart : jd_init_h[i];
			//Check if initial times are different for some bodies
			//If this is the case then a pre-integration is needed to the starting point
			if(i > 0 && jd_init_h[i] != jd_init_h[i-1]){
				//DoPreIntegration = 1;
				printf("Error, pre integration required, this is not yet supported\n");
				return 0;
			}
		}

		if(er <= 0.0){
			printf("Error in reading initial conditions\n");
			return 0;
		}

		unsigned long long int j = __builtin_bswap64 (id_h[i]);
		//j is the correct index

		if(i < 10 || i > N - 10) printf("%d %llu %llu | %.20g %.20g %.20g %.20g %.20g %.20g %.20g %.20g %.20g %.20g\n", i, id_h[i], j, jd_init_h[i], x_h[i], y_h[i], z_h[i], vx_h[i], vy_h[i], vz_h[i], A1_h[i], A2_h[i], A3_h[i]);
		id_h[i] = j;
		//if(j == 72057594038644513) printf("%d %llu %llu | %.20g %.20g %.20g %.20g %.20g %.20g %.20g %.20g %.20g %.20g\n", i, id_h[i], j, jd_init_h[i], x_h[i], y_h[i], z_h[i], vx_h[i], vy_h[i], vz_h[i], A1_h[i], A2_h[i], A3_h[i]);

		//if(A1_h[i] != 0.0 || A2_h[i] != 0.0 || A3_h[i] != 0.0){
		//	printf("A %d %g %g %g\n", i, A1_h[i], A2_h[i], A3_h[i]);
		//}
	}
	printf("Start time: %.20g\n", timeStart);

	//read trailer
	double temp;
	long long int header;
	er = fread(&header, sizeof(unsigned long long int), 1, inputFile);
	printf("trailer %lld %d\n", header, er);
	if(header != 130 || er <= 0.0){
		printf("Error in reading file trailer\n");
		return 0;
	}
	for(int i = 0; i < 9; ++i){
		er = fread(&temp, sizeof(double), 1, inputFile);
		printf(" %g %d\n", temp, er);
		if(temp != 0.0 || er <= 0.0){
			printf("Error in reading file trailer\n");
			return 0;
		}
	}
/*
	//check for new header

	er = fread(&header, sizeof(unsigned long long int), 1, infile);

	printf("header %lld %d\n", header, er);
	if(header != 129 || er <= 0.0){
		printf("Error in reading file header\n");
		return 0;
	}
*/ 
	return 1; 

} 


// *************************************************************


//code to read binary input file
#include "Host.h"


__host__ int Host::readHeader(FILE *infile, int &N){

	double temp;
	long long int header;
	int er;

	//read header:
	er = fread(&header, sizeof(unsigned long long int), 1, infile);

	printf("header %lld %d\n", header, er);
	if(header != 129 || er <= 0.0){
		printf("Error in reading file header\n");
		return 0;
	}

	er = fread(&outStart, sizeof(double), 1, infile);
	er = fread(&time1, sizeof(double), 1, infile);
	er = fread(&temp, sizeof(double), 1, infile);
	outInterval = (long long int)(temp);

	printf("start: %.20g, end: %.20g, interval: %lld\n", outStart, time1, outInterval);

	er = fread(&temp, sizeof(double), 1, infile);
	comet = int(temp);
	printf("comet flag %d\n", comet);

	er = fread(&temp, sizeof(double), 1, infile);
	N = int(temp);
	printf("N %d\n", N);

	er = fread(&time0, sizeof(double), 1, infile);
	printf("epoch (start time), used in Version 1: %.20g\n", time0);


	er = fread(&InVersion, sizeof(double), 1, infile);
	printf("File format version number: %g\n", InVersion);

	er = fread(&temp, sizeof(double), 1, infile);
	er = fread(&temp, sizeof(double), 1, infile);

	if(InVersion > 1.0){
		er = fread(&temp, sizeof(double), 1, infile);
	}
	return 1;

}
 
__host__ int Host::readFile(FILE *infile){


	//printf("size %lu %lu\n", sizeof(unsigned long long int), sizeof(double));

	int er = 0;

	if(InVersion > 1.0){
		//search for smallest starting time
		time0 = 1.0e10;
	}
	for(int i = Nperturbers; i < N; ++i){
		er = fread(&id_h[i], sizeof(unsigned long long int), 1, infile);
		er = fread(&x_h[i], sizeof(double), 1, infile);
		er = fread(&y_h[i], sizeof(double), 1, infile);
		er = fread(&z_h[i], sizeof(double), 1, infile);
		er = fread(&vx_h[i], sizeof(double), 1, infile);
		er = fread(&vy_h[i], sizeof(double), 1, infile);
		er = fread(&vz_h[i], sizeof(double), 1, infile);
		er = fread(&A1_h[i], sizeof(double), 1, infile);
		er = fread(&A2_h[i], sizeof(double), 1, infile);
		er = fread(&A3_h[i], sizeof(double), 1, infile);
		if(InVersion > 1.0){
			er = fread(&jd_init_h[i], sizeof(double), 1, infile);
			time0 = fmin(time0, jd_init_h[i]);
			//Check if initial times are different for some bodies
			//If this is the case then a pre-integration is needed to the starting point
			if(i > Nperturbers && jd_init_h[i] != jd_init_h[i-1]){
				DoPreIntegration = 1;
			}
		}

		if(er <= 0.0){
			printf("Error in reading initial conditions\n");
			return 0;
		}

		unsigned long long int j = __builtin_bswap64 (id_h[i]);
		//j is the correct index

		//if(i < 10 + Nperturbers || i > N - 10) printf("%d %llu %llu | %.20g %.20g %.20g %.20g %.20g %.20g %.20g %.20g %.20g %.20g\n", i, id_h[i], j, jd_init_h[i], x_h[i], y_h[i], z_h[i], vx_h[i], vy_h[i], vz_h[i], A1_h[i], A2_h[i], A3_h[i]);
		id_h[i] = j;
		if(j == 72057594038644513) printf("%d %llu %llu | %.20g %.20g %.20g %.20g %.20g %.20g %.20g %.20g %.20g %.20g\n", i, id_h[i], j, jd_init_h[i], x_h[i], y_h[i], z_h[i], vx_h[i], vy_h[i], vz_h[i], A1_h[i], A2_h[i], A3_h[i]);

		if(A1_h[i] != 0.0 || A2_h[i] != 0.0 || A3_h[i] != 0.0){
			printf("A %d %g %g %g\n", i, A1_h[i], A2_h[i], A3_h[i]);
		}
	}


	if(NMax >= N){
		//read trailer
		double temp;
		long long int header;
		er = fread(&header, sizeof(unsigned long long int), 1, infile);
		printf("trailer %lld %d\n", header, er);
		if(header != 130 || er <= 0.0){
			printf("Error in reading file trailer\n");
			return 0;
		}
		for(int i = 0; i < 9; ++i){
			er = fread(&temp, sizeof(double), 1, infile);
			printf(" %g %d\n", temp, er);
			if(temp != 0.0 || er <= 0.0){
				printf("Error in reading file trailer\n");
				return 0;
			}
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


/*
int readFIFO(double *x_h, double *y_h, double *z_h, double *vx_h, double *vy_h, double *vz_h, double *A1_h, double *A2_h, double *A3_h, unsigned long long int *id_h){

	int fd = open(myfifo,O_RDONLY);
	const int size = 7 * sizeof(double) + sizeof(int);
	char buffer[size];
	double time;

	for(int i = Nperturbers; i < NN; ++i){
		read(fd, &buffer, size);
		time = *reinterpret_cast<double*>(&buffer);
		id_h[i] = *reinterpret_cast<int*>(&buffer[8]);
		x_h[i] = *reinterpret_cast<double*>(&buffer[8+4]);
		y_h[i] = *reinterpret_cast<double*>(&buffer[2*8+4]);
		z_h[i] = *reinterpret_cast<double*>(&buffer[3*8+4]);
		vx_h[i] = *reinterpret_cast<double*>(&buffer[4*8+4]);
		vy_h[i] = *reinterpret_cast<double*>(&buffer[5*8+4]);
		vz_h[i] = *reinterpret_cast<double*>(&buffer[6*8+4]);
		printf("er %d %d %d %.20g %.20g %.20g\n", i, id_h[i], N, x_h[i], y_h[i], z_h[i]);
		++N;
	}
	close(fd);

}
*/


//Read the size of the initial conditions file
__host__ int Host::readICSize(){
	FILE *infile;

	infile = fopen(infilename, "r");
	if(infile == NULL){
		printf("Error, initial conditions file not found %s\n", infilename);
		return 0;
	}
	double x, y, z, vx, vy, vz;
	double A1, A2, A3;
	unsigned long long int id;
	
	int N = 0;

	for(int i = 0; i < NMax; ++i){
		int er = 0;
		fscanf(infile, "%lf", &time0);
		fscanf(infile, "%llu", &id);
		fscanf(infile, "%lf", &x);
		fscanf(infile, "%lf", &y);
		fscanf(infile, "%lf", &z);
		fscanf(infile, "%lf", &vx);
		fscanf(infile, "%lf", &vy);
		er = fscanf(infile, "%lf", &vz);
		er = fscanf(infile, "%lf", &A1);
		er = fscanf(infile, "%lf", &A2);
		er = fscanf(infile, "%lf", &A3);
//printf("read %d %g\n", er, A3);
		//fscanf(infile, "%lf", &ALN);
		//fscanf(infile, "%lf", &NK);
		//fscanf(infile, "%lf", &NM);
		//fscanf(infile, "%lf", &Nn);
		//er = fscanf(infile, "%lf", &R0);
		if(er < 0) break;
		++N;
	}
	fclose(infile);
	return N;
}

//Read the initial conditions file
__host__ int Host::readIC(){
	FILE *infile;

	infile = fopen(infilename, "r");
	for(int i = Nperturbers; i < N; ++i){
		int er = 0;
		fscanf(infile, "%lf", &time0);
		fscanf(infile, "%llu", &id_h[i]);
		fscanf(infile, "%lf", &x_h[i]);
		fscanf(infile, "%lf", &y_h[i]);
		fscanf(infile, "%lf", &z_h[i]);
		fscanf(infile, "%lf", &vx_h[i]);
		fscanf(infile, "%lf", &vy_h[i]);
		er = fscanf(infile, "%lf", &vz_h[i]);
		er = fscanf(infile, "%lf", &A1_h[i]);
		er = fscanf(infile, "%lf", &A2_h[i]);
		er = fscanf(infile, "%lf", &A3_h[i]);
		//fscanf(infile, "%lf", &ALN[i]);
		//fscanf(infile, "%lf", &NK[i]);
		//fscanf(infile, "%lf", &NM[i]);
		//fscanf(infile, "%lf", &Nn[i]);
		//er = fscanf(infile, "%lf", &R0[i]);
		if(er < 0){
			printf("Error, reading initial conditions file failed.\n");
			return 0;
		}
printf("xyz %.40g %.40g %.40g %.40g %.40g %.40g %.40g %.40g %.40g\n", x_h[i], y_h[i], z_h[i], vx_h[i], vy_h[i], vz_h[i], A1_h[i], A2_h[i], A3_h[i]);
//printf("er %d %llu %d %d %.20g %.20g %.20g\n", i, id_h[i], er, N, x_h[i], y_h[i], z_h[i]);
	}
	fclose(infile);

	return 1;

}

//convert velocities to code units
__host__ void Host::convertV(){

	for(int i = Nperturbers; i < N; ++i){
		vx_h[i] /= dayUnit;
		vy_h[i] /= dayUnit;
		vz_h[i] /= dayUnit;

		A1_h[i] /= (dayUnit * dayUnit);
		A2_h[i] /= (dayUnit * dayUnit);
		A3_h[i] /= (dayUnit * dayUnit);
	}
}

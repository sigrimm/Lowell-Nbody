//code to read binary input file

#include <stdio.h> 
#include <stdlib.h> 
 


int readHeader(FILE *infile, double &time0, double &time1, long long int &outInterval, double &outStart, int &N, int &comet){

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

	printf("epoch (start time): %.20g\n", time0);


	er = fread(&temp, sizeof(double), 1, infile);
	er = fread(&temp, sizeof(double), 1, infile);
	er = fread(&temp, sizeof(double), 1, infile);

	return 0;

}
 
int readFile(FILE *infile, int Nperturbers, double *x_h, double *y_h, double *z_h, double *vx_h, double *vy_h, double *vz_h, double *A1_h, double *A2_h, double *A3_h, unsigned long long int *id_h, int N){


	//printf("size %lu %lu\n", sizeof(unsigned long long int), sizeof(double));

	int er = 0;

	for(int i = Nperturbers; i < N + Nperturbers; ++i){
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

		if(er <= 0.0){
			printf("Error in reading initial conditions\n");
			return 1;
		}

		unsigned long long int j = __builtin_bswap64 (id_h[i]);
		//j is the correct index
		id_h[i] = j;

		/*if(i < 10 + Nperturbers || i > N + Nperturbers - 10)*/ printf("%d %llu %llu | %.20g %.20g %.20g %.20g %.20g %.20g %.20g %.20g %.20g\n", i, id_h[i], j, x_h[i], y_h[i], z_h[i], vx_h[i], vy_h[i], vz_h[i], A1_h[i], A2_h[i], A3_h[i]);

		if(A1_h[i] != 0.0 || A2_h[i] != 0.0 || A3_h[i] != 0.0){
			printf("A %d %g %g %g\n", i, A1_h[i], A2_h[i], A3_h[i]);
		}
	}

/*
	//read trailer
	double temp;
	long long int header;
	er = fread(&header, sizeof(unsigned long long int), 1, infile);
	printf("header %lld %d\n", header, er);
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


	//check for new header

	er = fread(&header, sizeof(unsigned long long int), 1, infile);

	printf("header %lld %d\n", header, er);
	if(header != 129 || er <= 0.0){
		printf("Error in reading file header\n");
		return 0;
	}
*/ 
    return 0; 
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

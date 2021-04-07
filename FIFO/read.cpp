//code to read binary input file

#include <stdio.h> 
#include <stdlib.h> 
  
int main(){


	FILE *infile;
	infile = fopen("210327_1522_genga_req.bin", "rb");

	unsigned long long int id;
	double x, y, z, vx, vy, vz, A1, A2, A3;

	printf("size %lu %lu\n", sizeof(unsigned long long int), sizeof(double));

	int er = 0;

	double start, end, outInterval, epoch, temp;
	int comet, N;
	long long int header;

	//read header:
	er = fread(&header, sizeof(unsigned long long int), 1, infile);

	printf("header %lld %d\n", header, er);
	if(header != 129 || er <= 0.0){
		printf("Error in reading file header\n");
		return 0;
	}

	er = fread(&start, sizeof(double), 1, infile);
	er = fread(&end, sizeof(double), 1, infile);
	er = fread(&outInterval, sizeof(double), 1, infile);

	printf("start: %.20g, end: %.20g, interval: %.20g\n", start, end, outInterval);

	er = fread(&temp, sizeof(double), 1, infile);
	comet = int(temp);
	printf("comet flag %d\n", comet);

	er = fread(&temp, sizeof(double), 1, infile);
	N = int(temp);
	printf("N %d\n", N);

	er = fread(&epoch, sizeof(double), 1, infile);

	printf("epoch %.20g\n", epoch);


	er = fread(&temp, sizeof(double), 1, infile);
	er = fread(&temp, sizeof(double), 1, infile);
	er = fread(&temp, sizeof(double), 1, infile);


	for(int t = 0; t < N; ++t){
		er = fread(&id, sizeof(unsigned long long int), 1, infile);
		er = fread(&x, sizeof(double), 1, infile);
		er = fread(&y, sizeof(double), 1, infile);
		er = fread(&z, sizeof(double), 1, infile);
		er = fread(&vx, sizeof(double), 1, infile);
		er = fread(&vy, sizeof(double), 1, infile);
		er = fread(&vz, sizeof(double), 1, infile);
		er = fread(&A1, sizeof(double), 1, infile);
		er = fread(&A2, sizeof(double), 1, infile);
		er = fread(&A3, sizeof(double), 1, infile);

		unsigned long long int j = __builtin_bswap64 (id);
		//j is the correct index

		if(t < 10 || t > N - 10) printf("%d | %llu %llu %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g\n", t, id, j, x, y, z, vx, vy, vz, A1, A2, A3);
	}

	//read trailer
	er = fread(&header, sizeof(unsigned long long int), 1, infile);
	printf("header %lld %d\n", header, er);
	if(header != 130 || er <= 0.0){
		printf("Error in reading file header\n");
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

	fclose(infile);
 
    return 0; 
} 


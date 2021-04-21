#include <stdio.h>
#include <stdlib.h>

int main(){


	FILE *infile;
	infile = fopen("Outhelio10.bin", "rb");

	unsigned long long int id;
	double time, x, y, z, vx, vy, vz;

	int N = 1;

	int er;

	for(int t = 0; t < 1e6; ++t){
		for(int i = 0; i < N; ++i){
			er = fread(&id, sizeof(unsigned long long int), 1, infile);
			er = fread(&time, sizeof(double), 1, infile);
			er = fread(&x, sizeof(double), 1, infile);
			er = fread(&y, sizeof(double), 1, infile);
			er = fread(&z, sizeof(double), 1, infile);
			er = fread(&vx, sizeof(double), 1, infile);
			er = fread(&vy, sizeof(double), 1, infile);
			er = fread(&vz, sizeof(double), 1, infile);

			if(er <= 0) break;
			unsigned long long int j = __builtin_bswap64 (id);
			//j is the correct index

			printf("%.15g %llu %g %.15g %.15g %.15g %.15g %.15g %.15g\n", time, j, 0.0, x, y, z, vx, vy, vz);
		}
		if(er <= 0) break;
	}


	fclose(infile);
	return 0;
}

#include <stdio.h>
#include <stdlib.h>

int main(){


	FILE *infile;
	infile = fopen("Outhelio.bin", "rb");
	//infile = fopen("Out_220705_1512_genga_in_query_genga_input_100k.bin", "rb");

	unsigned long long int id;
	double time, x, y, z, vx, vy, vz, dtmin;

	int N = 1024;

	int er;

	for(int t = 0; t < 1e8; ++t){
	//for(int t = 0; t < 10; ++t){
		for(int i = 0; i < N; ++i){
			er = fread(&id, sizeof(unsigned long long int), 1, infile);
			er = fread(&time, sizeof(double), 1, infile);
			er = fread(&x, sizeof(double), 1, infile);
			er = fread(&y, sizeof(double), 1, infile);
			er = fread(&z, sizeof(double), 1, infile);
			er = fread(&vx, sizeof(double), 1, infile);
			er = fread(&vy, sizeof(double), 1, infile);
			er = fread(&vz, sizeof(double), 1, infile);
			er = fread(&dtmin, sizeof(double), 1, infile);

			if(er <= 0) break;
			unsigned long long int j = __builtin_bswap64 (id);
			//j is the correct index

			printf("%.15g %llu %g %.15g %.15g %.15g %.15g %.15g %.15g %.15g\n", time, j, 0.0, x, y, z, vx, vy, vz, dtmin);
		}
		if(er <= 0) break;
	}


	fclose(infile);
	return 0;
}

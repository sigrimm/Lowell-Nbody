#include <stdio.h>
#include <stdlib.h>

//code to convert plpos.dat binary file into perturbers file.

int main(){

	int useHelioCentric = 1;
	int useBinary = 1;


	FILE *infile;
	infile = fopen("plpos.dat", "rb");

	const int N = 27;

	FILE *outfile;


	if(useHelioCentric == 1){
		if(useBinary == 0){
			outfile = fopen("All3_h.dat", "w");
		}
		else{
			outfile = fopen("All3_h.bin", "wb");
		}
	}
	else{
		if(useBinary == 0){
			outfile = fopen("All3_b.dat", "w");
		}
		else{
			outfile = fopen("All3_b.bin", "wb");
		}
	}



	double t0, t1, dt;
	double skip;
	double x[N], y[N], z[N];


	int er;


	//read header
	er = fread(&t0, sizeof(double), 1, infile);
	er = fread(&t1, sizeof(double), 1, infile);
	er = fread(&dt, sizeof(double), 1, infile);
	for(int i = 0; i < N-1; ++i){
		er = fread(&skip, sizeof(double), 1, infile);
		er = fread(&skip, sizeof(double), 1, infile);
		er = fread(&skip, sizeof(double), 1, infile);
		//printf("%.20g %.20g %.20g\n", t0, t1, dt);
	}

	//printf("\n");

	for(int t = 0; t < 1e8; ++t){
	//for(int t = 0; t < 100; ++t){
		for(int i = 0; i < N; ++i){
			er = fread(&x[i], sizeof(double), 1, infile);
		}
		for(int i = 0; i < N; ++i){
			er = fread(&y[i], sizeof(double), 1, infile);
		}
		for(int i = 0; i < N; ++i){
			er = fread(&z[i], sizeof(double), 1, infile);
		}

		if(er <= 0) break;


		if(useBinary == 0){
			if(useHelioCentric == 0){
				fprintf(outfile, "%.20g %d %.30g %.30g %.30g\n", t0 + t, 27, x[N-1], y[N-1], z[N-1]);
			}
			else{
				fprintf(outfile,"%.20g %d %.30g %.30g %.30g\n", t0 + t, 27, 0.0, 0.0, 0.0);
			}

			for(int i = 0; i < N-1; ++i){
				if(useHelioCentric == 0){
					fprintf(outfile,"%.20g %d %.30g %.30g %.30g\n", t0 + t, i, x[i], y[i], z[i]);
				}
				else{
					fprintf(outfile,"%.20g %d %.30g %.30g %.30g\n", t0 + t, i, x[i]-x[N-1], y[i]-y[N-1], z[i]-z[N-1]);
				}
			}
		}
		else{
			double tt = t0 + t;
			if(useHelioCentric == 0){
				fwrite(&tt, sizeof(double), 1, outfile);
				fwrite(&x[N-1], sizeof(double), 1, outfile);
				fwrite(&y[N-1], sizeof(double), 1, outfile);
				fwrite(&z[N-1], sizeof(double), 1, outfile);
				for(int i = 0; i < N-1; ++i){
					fwrite(&tt, sizeof(double), 1, outfile);
					fwrite(&x[i], sizeof(double), 1, outfile);
					fwrite(&y[i], sizeof(double), 1, outfile);
					fwrite(&z[i], sizeof(double), 1, outfile);
				}
			}
			else{
				double xx = 0.0;
				double yy = 0.0;
				double zz = 0.0;
				fwrite(&tt, sizeof(double), 1, outfile);
				fwrite(&xx, sizeof(double), 1, outfile);
				fwrite(&yy, sizeof(double), 1, outfile);
				fwrite(&zz, sizeof(double), 1, outfile);
				for(int i = 0; i < N-1; ++i){
					xx = x[i]-x[N-1];
					yy = y[i]-y[N-1];
					zz = z[i]-z[N-1];
					fwrite(&tt, sizeof(double), 1, outfile);
					fwrite(&xx, sizeof(double), 1, outfile);
					fwrite(&yy, sizeof(double), 1, outfile);
					fwrite(&zz, sizeof(double), 1, outfile);
					//printf("%.20g %.30g %.30g %.30g\n", tt, xx, yy, zz);
				}
			}
		}
	}


	fclose(infile);
	return 0;
}

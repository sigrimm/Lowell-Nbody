#include <stdio.h>
#include <stdlib.h>


int main(){

	int useHelioCentric = 1;
	int useBinary = 0;

	const int N = 27;

	const char *names[N] = {"Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto", 
		 "Ceres", "Pallas", "Juno", "Vesta", "Hygiea", "Eunomia", "Euphrosyne", "Europa", "Davida", "Interamnia",
		 "Psyche", "Cybele", "Thisbe", "Camilla", "Iris", "Sylvia", "Moon"};


	FILE *infile[N];
	char infilename[N][160];

	FILE *outfile;

	double time, x, y, z, vx, vy, vz;

	if(useHelioCentric == 1){
		if(useBinary == 0){
			outfile = fopen("All2_h.dat", "w");
		}
		else{
			outfile = fopen("All2_h.bin", "wb");
		}
	}
	else{
		if(useBinary == 0){
			outfile = fopen("All2_b.dat", "w");
		}
		else{
			outfile = fopen("All2_b.bin", "wb");
		}
	}

	for(int i = 0; i < N; ++i){
		//printf("%s_h.dat\n", names[i]);
		sprintf(infilename[i], "%s_h.dat", names[i]);
		infile[i] = fopen(infilename[i], "r");

	}

	int er; 
	for(int t = 0; t < 1e6; ++t){
		for(int i = 0; i < N; ++i){
			fscanf(infile[i], "%lf", &time);
			fscanf(infile[i], "%lf", &x);
			fscanf(infile[i], "%lf", &y);
			fscanf(infile[i], "%lf", &z);
			fscanf(infile[i], "%lf", &vx);
			fscanf(infile[i], "%lf", &vy);
			er = fscanf(infile[i], "%lf", &vz);

			if(er <= 0){
				break;
			}
	
			if(useBinary == 0){
				fprintf(outfile, "%.20g %d %.30g %.30g %.30g %.30g %.30g %.30g\n", time, i, x, y, z, vx, vy, vz);
			}
			else{
				fwrite(&time, sizeof(double), 1, outfile);
				fwrite(&x, sizeof(double), 1, outfile);
				fwrite(&y, sizeof(double), 1, outfile);
				fwrite(&z, sizeof(double), 1, outfile);
			}
		}
		if(er <= 0){
			break;
		}
	}

	for(int i = 0; i < N; ++i){
		fclose(infile[i]);
	}	

	return 0;
}

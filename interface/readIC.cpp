#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>


//This code serves to test the interface betwen the database and the integrator
//The code reads the file initial.dat and prints its output.

//The goal is to pass an object with the initial conditions to the code
//without reading a file


//compile with:
//g++ -o readIC readIC.cpp

int main(int argc, char*argv[]){

	//Number of planets
	const int NN = 28;	//22
	const int Nperturbers = 27;	//21
	const int Ninterpolate = 10;	//number of interpolation points

	int Nsteps = 0;
	int outInterval = 0;
	double dt = 0.0;

	for(int i = 1; i < argc; i += 2){

		if(strcmp(argv[i], "-Nsteps") == 0){
			Nsteps = atoll(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-outInterval") == 0){
			outInterval = atoll(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-dt") == 0){
			dt = atof(argv[i + 1]);
		}
	}

	
	int id[NN];
	double m[NN];
	double x[NN];
	double y[NN];
	double z[NN];
	double vx[NN];
	double vy[NN];
	double vz[NN];

	double time = 0.0;

	int N = Nperturbers;

	FILE *infile;
	char infilename[160];

	sprintf(infilename, "initial.dat");
	infile = fopen(infilename, "r");
	//sun
	id[0] = 20;
	x[0] = 0.0;
	y[0] = 0.0;
	z[0] = 0.0;
	vx[0] = 0.0;
	vy[0] = 0.0;
	vz[0] = 0.0;

	for(int i = 1; i < NN; ++i){
		id[i] = -1;
		x[i] = 0.0;
		y[i] = 0.0;
		z[i] = 0.0;
		vx[i] = 0.0;
		vy[i] = 0.0;
		vz[i] = 0.0;
	}
	//read test particle
	for(int i = Nperturbers; i < NN; ++i){
		int er = 0;
		fscanf(infile, "%lf", &time);
		fscanf(infile, "%lf", &x[i]);
		fscanf(infile, "%lf", &y[i]);
		fscanf(infile, "%lf", &z[i]);
		fscanf(infile, "%lf", &vx[i]);
		fscanf(infile, "%lf", &vy[i]);
		er = fscanf(infile, "%lf", &vz[i]);
		if(er < 0) break;
		++N;
	}
	fclose(infile);

	for(int i = Nperturbers; i < N; ++i){
		printf("%d %.20g %.20g %.20g %.20g %.20g %.20g %.20g\n", i, time, x[i], y[i], z[i], vx[i], vy[i], vz[i]);
	}
}
	

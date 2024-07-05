#ifndef ASTEROID_H
#define ASTEROID_H

#include <stdio.h>
#include <stdlib.h>


class asteroid{


public:
        FILE *infile;

	int N;
        int Nperturbers;
        double time_reference;
        double timeStart;		//start time of integration
        double timeEnd;			//end time of integration
        double dt;
	int dts;			//sign of time step
	double dt1;
        double time0;			//start time of the data file
        double time1;			//end time of the data file
	double time;			//integration time
	int stop;			//used to refine last time step


	double AUtokm;          //AU to km
	double EM;              //Earth to moon mass ratio
	double CLIGHT;          //speed of light
	double RE;              //Radius of Earth
	double J2E;             //J2 of Earth
	double c2;		//speed of light squared
	double REAU;		//Earth radius in AU


	int nCm;		//Maximum number of Chebyshev coefficients


	//perturbers data
	double *startTime;
	double *endTime;
	int *id;
	int *nChebyshev;
	int *offset0;
	int *offset1;
	double *GM;

	double *cdata;

	double *x;
	double *y;
	double *z;

	double *vx;
	double *vy;
	double *vz;

	double *ax;
	double *ay;
	double *az;

	//Non-grav terms
	double *A1;
	double *A2;
	double *A3;


	void allocate();
	inline void update_Chebyshev(double);
	inline void update_perturbers(double);
	inline void NonGrav();
	inline void GR();
	inline void J2();
	inline void Gravity();
	inline void leapfrog_step();
	inline int loop();
};
#endif

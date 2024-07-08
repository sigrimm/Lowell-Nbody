#ifndef ASTEROID_H
#define ASTEROID_H

#include <stdio.h>
#include <stdlib.h>
#include <cstring>

class asteroid{


public:
        FILE *inputFile;
	FILE *perturbersFile;
	char inputFilename[160];
        FILE *outputFile;

	int N;
        int Nperturbers;
        double time_reference;
        double timeStart;		//start time of integration
        double timeEnd;			//end time of integration
        double dt;			//time step
	int dts;			//sign of time step
	double dt1;
	double outputInterval;		
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



	//RKF arrays
	int RKFn;
        double *a_h;
	double *b_h;
	double *bb_h;
	double *c_h;
        double RKF_ee;
	double RKF_atol;
	double RKF_rtol;
	double RKF_fac;
	double RKF_facmin;
	double RKF_facmax;

	double *xt;
	double *yt;
	double *zt;

	double *vxt;
	double *vyt;
	double *vzt;

	double *dx;
	double *dy;
	double *dz;

	double *dvx;
	double *dvy;
	double *dvz;

	double *kx;
	double *ky;
	double *kz;

	double *kvx;
	double *kvy;
	double *kvz;


	int readParam();
	int readIC();
	void allocate();
	inline void update_Chebyshev(double);
	inline void update_perturbers(double);
	inline void NonGrav(double *, double *, double *, double *, double *, double *);
	inline void GR(double *, double *, double *, double *, double *, double *);
	inline void J2(double *, double *, double *);
	inline void Gravity(double *, double *, double *);
	inline int loop();

	inline void leapfrog_step();
	inline void RK_step();
	inline void RKF_step();

	void setRK4();
	void setRKF45();
	void setDP54();
	void setRKF78();

};
#endif

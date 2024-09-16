#ifndef ASTEROID_H
#define ASTEROID_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstring>


#define def_NP 32		//used for shared memory, must be at least the number of perturbers

class asteroid{


public:
        FILE *inputFile;
	FILE *perturbersFile;
	char perturbersFilePath[160];
	char perturbersFileName[240]; // 160 + 80
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

	double AUtokm;			//AU to km
	double EM;			//Earth to moon mass ratio
	double CLIGHT;			//speed of light
	double RE;			//Radius of Earth
	double J2E;			//J2 of Earth
	double c2;			//speed of light squared
	double REAU;			//Earth radius in AU


	int nCm;			//Maximum number of Chebyshev coefficients
	int datasize;

	//perturbers data
	double *startTime_h, *startTime_d;	//Start time of perturbers data block
	double *endTime_h, *endTime_d;		//End time of perturbers data block
	int *id_h, *id_d;		
	int *nChebyshev_h, *nChebyshev_d;	//Number of Chebyshev coefficients of perturbers in current data block
	int *offset0_h, *offset0_d;
	int *offset1_h, *offset1_d;
	double *GM_h, *GM_d;

	double *cdata_h, *cdata_d;		//contains one record of data for each perturber
	double *data_h, *data_d;		//entire data file, only in GPU version


	double *xTable_h, *xTable_d;		//Perturbers table contains the position of the perturbers for all stages
	double *yTable_h, *yTable_d;
	double *zTable_h, *zTable_d;

	double *vxTable_h, *vxTable_d;
	double *vyTable_h, *vyTable_d;
	double *vzTable_h, *vzTable_d;


	double *x_h, *x_d;			//Contains only integration bodies, not perturbers
	double *y_h, *y_d;
	double *z_h, *z_d;

	double *vx_h, *vx_d;
	double *vy_h, *vy_d;
	double *vz_h, *vz_d;

	double *ax_h, *ax_d;
	double *ay_h, *ay_d;
	double *az_h, *az_d;

	//Non-grav terms
	double *A1_h, *A1_d;
	double *A2_h, *A2_d;
	double *A3_h, *A3_d;



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

	double *xt_h, *xt_d;
	double *yt_h, *yt_d;
	double *zt_h, *zt_d;

	double *vxt_h, *vxt_d;
	double *vyt_h, *vyt_d;
	double *vzt_h, *vzt_d;

	double *dx_h, *dx_d;
	double *dy_h, *dy_d;
	double *dz_h, *dz_d;

	double *dvx_h, *dvx_d;
	double *dvy_h, *dvy_d;
	double *dvz_h, *dvz_d;

	double *kx_h, *kx_d;
	double *ky_h, *ky_d;
	double *kz_h, *kz_d;

	double *kvx_h, *kvx_d;
	double *kvy_h, *kvy_d;
	double *kvz_h, *kvz_d;


	double *snew_h, *snew_d;


	int readParam();
	int readIC();
	int readICSize();
	int readData();
	int copyIC();
	void allocate();
	int allocateGPU();
	void copyOutput();
	int copyConst();

	inline void update_Chebyshev(double);
	inline void update_perturbers(double);
	inline void NonGrav(double, double, double, double, double, double, double, double, double, double, double &, double &, double &);
	inline void GR(double, double, double, double, double, double, double, double &, double &, double &, double);
	inline void J2(double, double, double, double &, double &, double &, double);
	inline void Gravity(double, double, double, double *, double *, double *, double &, double &, double &, int);
	int loop();

	inline void leapfrog_step();
	inline void RK_step();
	inline void RKF_step();

	void setRK4();
	void setRKF45();
	void setDP54();
	void setRKF78();

};
#endif

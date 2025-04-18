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
	FILE *outputFile;
	FILE *infoFile;

	char name[140];
	char inputFilename[160];
	char perturbersFilePath[160];
	char perturbersFileName[240]; // 160 + 80
	char outputFilename[160];
	char infoFilename[160];

	int ICformat = 0;			//Format of the initial conditions file, 0 = "text" or 1 = "binary"
	int ICorbital = 0;			//Coordinate system of initial conditions, 0 = cartesian, 1 = orbital
	int ICecliptic = 0;			//Coordinate system of initial conditions, 1 = ecliptic, 0 = equatorial
	int ICheliocentric = 0;			//Coordinate system of initial conditions, 1 = heliocentric, 0 = barycentric
	int Outorbital = 0;			//Coordinate system of output, 0 = cartesian, 1 = orbital
	int Outecliptic = 0;			//Coordinate system of output, 1 = ecliptic, 0 = equatorial
	int Outheliocentric = 0;		//Coordinate system of output, 1 = heliocentric, 0 = barycentric
	double Obliquity = 84381.448;		//Obliquity of ecliptic IAU76/J2000	in arcsec

	int N;					//Number of bodies to integrate, not including perturbers
	int Nperturbers;			//Number of perturbers
	double inputFileVersion = 1;
	double time_reference = 2451545.0;	//Reference time in JD 
	double timeStart = 2450800.5;		//Start time of the integration, in JD
	double timeEnd = 2461000.5;		//End time of the integration, in JD
	double dt = 1.0;			//time step, in days
	int dts;				//sign of time step
	double dt1;
	double outputInterval = 10;		//outut interval, in days		
	double outStart = 0.0;			//start time when outputs are written, in JD
	double time0;				//start time of the data file
	double time1;				//end time of the data file
	double time;				//integration time
	long long int timeStep = 0ll;
	int stop;				//used to refine last time step

	double AUtokm;				//AU to km
	double EM;				//Earth to moon mass ratio
	double CLIGHT;				//speed of light
	double RE;				//Radius of Earth
	double J2E;				//J2 of Earth
	double c2;				//speed of light squared
	double REAU;				//Earth radius in AU

	int useGR = 1;
	int useJ2 = 1;
	int useNonGrav = 1;
	int outBinary = 0;
	int cometFlag = 0;

	double nonGrav_alpha = 0.1112620426;	//nongrav normalizing factor ALN
	double nonGrav_nk = 4.6142;		//nongrav model constant NK
	double nonGrav_nm = 2.15;		//nongrav model constant NM
	double nonGrav_nn = 5.093;		//nongrav model constant NN
	double nonGrav_r0 = 2.808;		//nongrav normalizing distance in AU
	double nonGrav_tau = 0.0;		//nongrav time delay in days

	int Rbuffersize = 500000;		//array size of Rsave and Tsave
	double *Rsave_h, *Rsave_d;			//Used for time delay
	double *Tsave_h, *Tsave_d;			//Used for time delay


	int nCm;			//Maximum number of Chebyshev coefficients
	int datasize;

	int WarpSize = 0;
	int GPUMode = 0;

	//perturbers data
	double *startTime_h, *startTime_d;	//Start time of perturbers data block
	double *endTime_h, *endTime_d;		//End time of perturbers data block
	int *idp_h, *idp_d;		
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

	double *xout_h, *xout_d;			//Contains only integration bodies, not perturbers
	double *yout_h, *yout_d;
	double *zout_h, *zout_d;

	double *vxout_h, *vxout_d;
	double *vyout_h, *vyout_d;
	double *vzout_h, *vzout_d;


	double *ax_h, *ax_d;
	double *ay_h, *ay_d;
	double *az_h, *az_d;

	//Non-grav terms
	double *A1_h, *A1_d;
	double *A2_h, *A2_d;
	double *A3_h, *A3_d;
	double *R_h, *R_d;



	//RKF arrays
	int RKFn = 6;
        double *RKFa_h;
	double *RKFb_h;
	double *RKFbb_h;
	double *RKFc_h;

        double RKF_ee;
	double RKF_atol = 1.0e-16;
	double RKF_rtol = 1.0e-16;
	double RKF_fac = 0.84;
	double RKF_facmin = 0.8;
	double RKF_facmax = 1.5;

	double *xt_h;
	double *yt_h;
	double *zt_h;

	double *vxt_h;
	double *vyt_h;
	double *vzt_h;

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


	long long int *id_h, *id_d;

	double *snew_h, *snew_d;
	double *ssum_d;
	double *jd_init_h;

	int readParam(int , char*argv[]);
	int readIC();
	int readICkeplerian();
	int readICSize();
	int readData();
	int readHeader();
	int readFile();
	int copyIC();
	int allocate();
	int allocateGPU();
	void copyOutput();
	int copyConst();
	void printInfo();
	void printOutput(double);
	

	void KepToCart_M(double *, double *, double *, double *, double *, double *, int, double, double, double, double, double, double);
	void KepToCart_E(double *, double *, double *, double *, double *, double *, int, double, double, double, double, double, double);
	void CartToKep(double *, double *, double *, double *, double *, double *, int, double &, double &, double &, double &, double &, double &, double &, double &);
	void convertOutput();
	void EcpliptictoEquatorial(double *, double *, double *, double *, double *, double *);
	void EquatorialtoEcliptic(double *, double *, double *, double *, double *, double *);
	void HelioToBary(double *, double *, double *, double *, double *, double *);
	void BaryToHelio(double *, double *, double *, double *, double *, double *);

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

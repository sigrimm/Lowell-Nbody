#include "define.h"

#ifndef HOST_CLASS
#define HOST_CLASS


class Host{

public:
	int WarpSize = 32;

	const int dtimep = 1.0;		//interval between stored time steps
	const int NTable = 65000;	//length of perturbers table, number of days
	double timep0;			//start time from perturbers file


	int comet;			//flag for asteroids of comets

	int useGR;			//2: Sitarski 1982, heliocentric coordinates
	int useJ2;
	int useNonGrav;

	int useGPU;
	int useAdaptiveTimeSteps;
	int useIndividualTimeSteps;

	int useFIFO;
	double InVersion;		//version of Input file format
	int DoPreIntegration;		//If this is 1 then do a pre-Integration to synchronize all initial conditions


	int useHelio;
	int outHelio;			//1 print output in heliocentric coordinates
					//0 print output in barycentric coordinates


	int outBinary;			//0 print output files in ASCII format  
					//1 print output files in binary format 

	long long int Nsteps;
	long long int outInterval;


	double time0;			//start time from simulation
	double time1;			//end time from simulation
	double outStart;		//start time of output files
	double time;

	double dti;			//time step size in days
	double dt;			//time step size in code units
	double *dtiMin;			//Minimal time step size

	double dts;			//interval of interpolated points

	int N;				//Number of bodies
	int NMax;			//Maximum number of bodies

	//integration data
	unsigned long long int *id_h, *id_d;
	unsigned int *index_h, *index_d;
	double *m_h, *m_d;
	double *x_h, *x_d;
	double *y_h, *y_d;
	double *z_h, *z_d;
	double *vx_h, *vx_d;
	double *vy_h, *vy_d;
	double *vz_h, *vz_d;
	double *A1_h, *A1_d;
	double *A2_h, *A2_d;
	double *A3_h, *A3_d;
	double *jd_init_h, *jd_init_d;

	double *dx_h, *dx_d;
	double *dy_h, *dy_d;
	double *dz_h, *dz_d;
	double *dvx_h, *dvx_d;
	double *dvy_h, *dvy_d;
	double *dvz_h, *dvz_d;


	//coordinates from data table
	double *timep_h, *timep_d;
	double *xp_h, *xp_d;
	double *yp_h, *yp_d;
	double *zp_h, *zp_d;

	//Backup for forward/backward integration
	double *x0_h, *x0_d;
	double *y0_h, *y0_d;
	double *z0_h, *z0_d;
	double *vx0_h, *vx0_d;
	double *vy0_h, *vy0_d;
	double *vz0_h, *vz0_d;
	double *A10_h, *A10_d;
	double *A20_h, *A20_d;
	double *A30_h, *A30_d;
	double *m0_h, *m0_d;
	unsigned long long int *id0_h, *id0_d;
	unsigned int *index0_h, *index0_d;

	double *xb_h;
	double *yb_h;
	double *zb_h;
	double *vxb_h;
	double *vyb_h;
	double *vzb_h;
	double *A1b_h;
	double *A2b_h;
	double *A3b_h;
	double *mb_h;
	unsigned long long int *idb_h;


	//Files
	FILE *infile;
	char *infilename;  
	FILE *outfile;
	char *outfilename;  
	FILE *dtfile;
	char *dtfilename;
	FILE *XVfile; 


	//perturber buffers
	//the buffer contains time, x, y, ,z from all perturbers 
	//the buffer has 2 swaps
	double *readBufferA_h;
	double *readBufferB_h;
	double *XYdata_h, *XYdata_d;


	//additional arrays
	double *xt_h;
	double *yt_h;
	double *zt_h;
	double *vxt_h;
	double *vyt_h;
	double *vzt_h;

	double *kx_h, *kx_d;
	double *ky_h, *ky_d;
	double *kz_h, *kz_d;
	double *kvx_h, *kvx_d;
	double *kvy_h, *kvy_d;
	double *kvz_h, *kvz_d;

	double2 *snew_h, *snew_d;
	double *dtmin_h, *dtmin_d;

	int2 *scan_d;
	int *N_d;

	//interpolation table
	double *xTable_h, *xTable_d;
	double *yTable_h, *yTable_d;
	double *zTable_h, *zTable_d;

	//Runge Kutta Fehlberg variables
	double *a_h, *b_h, *bb_h, *c_h;
	double ee;

	int nRuns;
	int *runsN;
	double *runsdt;

	//Functions
	__host__ Host();
	__host__ int readparam(int, char*argv[]);

	__host__ int readHeader(FILE*, int &);
	__host__ int readFile(FILE *);
	__host__ int readICSize();
	__host__ int readIC();
	__host__ void Alloc1();
	__host__ void Alloc2();
	__host__ void initialize1();
	__host__ void initialize2();
	__host__ void initialize3();
	__host__ void restore3();
	__host__ void setSnew();
	__host__ void save(double);
	__host__ void save1();
	__host__ void perturbersMass();
	__host__ void perturbersIDs();
	__host__ void convertV();
	__host__ int readTable();
	__host__ void copy1();
	__host__ void copyConst();
	__host__ void copyPerturbersCPU();

	__host__ void setRKF78();
	__host__ void setDP54();
	__host__ void setRKF45();

	__host__ int preIntegration();
	__host__ void IntegrationLoop(int, unsigned long long int, double, double&);
	__host__ void interpolate(double, int);
	__host__ void interpolate2(double, int);
	__host__ void interpolateTable(double);
	__host__ void update(int);
	__host__ void stageStep2(int, double, double, double &);
	__host__ void stageStep(double, double, double &);
	__host__ void stageStep1(double, int, int, double &);

	__host__ void reduceCall(int);	
	__host__ void reduce(int);
	__host__ void output(unsigned long long int, double);

//privat:

};
#endif

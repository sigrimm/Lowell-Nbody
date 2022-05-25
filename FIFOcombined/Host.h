#include "define.h"

#ifndef HOST_CLASS
#define HOST_CLASS


class Host{

public:

	const int dtimep = 1.0;		//interval between stored time steps
	const int NTable = 65000;	//length of perturbers table, number of days
	double timep0;			//start time from perturbers file


	int comet;			//flag for asteroids of comets

	int useGR;			//2: Sitarski 1982, heliocentric coordinates
	int useJ2;
	int useNonGrav;

	int useGPU;

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
	double dtiMin;			//Minimal time step size

	double dts;			//interval of interpolated points
	unsigned long long int outI;	

	int N;				//Number of bodies


	//integration data
	unsigned long long int *id_h, *id_d;
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

	//Backup for repeated integrations
	double *x0_h;
	double *y0_h;
	double *z0_h;
	double *vx0_h;
	double *vy0_h;
	double *vz0_h;


	//Files
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
	double *dtmin_h;


	//interpolation table
	double *xTable_h, *xTable_d;
	double *yTable_h, *yTable_d;
	double *zTable_h, *zTable_d;

	//Runge Kutta Fehlberg variables
	double *a_h, *b_h, *bb_h, *c_h;
	double ee;

	int Sn[4];	

	//Functions
	__host__ Host();
	__host__ int readHeader(FILE*, int &);
	__host__ int readFile(FILE *);
	__host__ int readICSize();
	__host__ int readIC();
	__host__ void Alloc1();
	__host__ void Alloc2();
	__host__ void Initialize1();
	__host__ void Initialize2();
	__host__ void Initialize3();
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
	__host__ void IntegrationLoop(int);
	__host__ void interpolate(double, int);
	__host__ void interpolate2(double, int);
	__host__ void interpolateTable(double);
	__host__ void update(int);
	__host__ void stageStep2(int, double, double, double &);
	__host__ void stageStep(double, double, double &);
	__host__ void stageStep1(double, int, double &);

	__host__ void reduce(int);
	__host__ void output(long long int, double, int);

//privat:

};
#endif

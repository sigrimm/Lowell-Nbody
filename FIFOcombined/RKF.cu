//for FIFO
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

//#include "define.h"
#include "Host.h"

#include "force.h"
#include "integrator.h"



int main(int argc, char*argv[]){


	Host H;

	FILE *binfile;
	if(H.useFIFO == 2){	
		//binfile = fopen("210801_2342_genga_de440_perturbers.bin", "rb");
		//binfile = fopen("210921_2148_genga_in_yarkovsky_elements.bin", "rb");
		//binfile = fopen("211208_1916_genga_in_2021-12-08_specific_desig.bin", "rb");
		//binfile = fopen("210801_2104_genga_in_GA.bin", "rb");
		//binfile = fopen("210705_2315_genga_req.bin", "rb");
		//binfile = fopen("220301_2048_genga_in_new_last_14_days.bin", "rb");
		//binfile = fopen("220524_2256_genga_in_query_genga_input_40k.bin", "rb");
		binfile = fopen("220524_2258_genga_in_query_genga_input_10k.bin", "rb");

		if(binfile == NULL){
			printf("Error, input file not found.\n");
			return 0;
		}
	}	

	double timing[6];
	for(int i = 0; i < 6; ++i){
		timing[i] = 0.0;
	}
	float milliseconds = 0.0f;
	cudaError_t error;
	int er;

	cudaEvent_t tt1;		//start time for timing
	cudaEvent_t tt2;		//end time for timing
	cudaEventCreate(&tt1);
	cudaEventCreate(&tt2);

	cudaEventRecord(tt1);

	// **************************************************
	// Read header or size of the initial conditions file
	// **************************************************

	if(H.useFIFO == 1){
		//This is only for testing
		printf("Error, useFIFO = 1 not supported.\n");
		return 0;

		const char *myfifo = "myfifo";
		const char *fifoCheck = "fifoCheck";
		// ###############################
		//create FIFO
		// ###############################
		int nn = 0;
		int fd;
		mkfifo(myfifo, 0666); //path, permission mode
		mkfifo(fifoCheck, 0666); //path, permission mode

		// ###############################
		// read N
		// ###############################
		fd = open(myfifo,O_RDONLY);
		read(fd, &nn, sizeof(int));
		close(fd);
		printf("fifo n: %d\n", nn);
		// ###############################
		// send back N to check
		// ###############################
		int fd1;
		fd1 = open(fifoCheck, O_WRONLY);	
		write(fd1, &nn, sizeof(int));
		close(fd1);
		printf("sent back\n");
	}
	

	if(H.useFIFO == 2){
		int NTP = 0;
		printf("read file\n");

		er = H.readHeader(binfile, NTP);
		if(er == 0){
			return 0;
		}

		//H.NTP = 1;

		H.N += NTP;
	}

	if(H.useFIFO == 0){
		int n = H.readICSize();
		if(n == 0){
			return 0;
		}
		H.N += n;
		H.outStart = H.time0;
H.time1 = 2461000.5;
//H.time1 = H.time0 + H.dt * H.Nsteps;

	}
	H.time = H.time0;
	// **************************************************
	//H.time1 = 2451000.5;

	//move this to parameter file
	H.dti = 5.0;
	H.dts = 0.1;
	H.dtiMin = 2.0;
	//H.outInterval = 10.0;
	//H.outStart = 2450800.5;
	H.outI = (H.outInterval + 0.5 * H.dts) / H.dts;

	//	if(fabs(outStart - time0) > fabs(outInterval)){
	//		outI = (outStart - time0 + 0.5 * dts) / dts;
	//	}

	printf("outStart: %.20g, time0: %.20g, time1: %.20g, dts: %g, outI: %llu,  outInterval: %lld\n", H.outStart, H.time0, H.time1, H.dts, H.outI, H.outInterval);



	for(int i = 1; i < argc; i += 2){

		if(strcmp(argv[i], "-Nsteps") == 0){
			H.Nsteps = atoll(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-outInterval") == 0){
			H.outInterval = atoll(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-dt") == 0){
			H.dti = atof(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-endTime") == 0){
			H.time1 = atof(argv[i + 1]);
		}
	}
	H.dt = H.dti * dayUnit;
	printf("Nperturbers: %d N:%d\n", Nperturbers, H.N);


	// **************************************************
	// Allocate CPU and GPU memory
	// **************************************************
	H.Alloc1();

	cudaDeviceSynchronize();
	error = cudaGetLastError();
	printf("Alloc error = %d = %s\n",error, cudaGetErrorString(error));
	if(error != 0.0){
		return 0;
	}

	// **************************************************

	H.Initialize1();


	// **************************************************
	// Read initial conditions
	// **************************************************

	printf("Read initial conditions\n");
	if(H.useFIFO == 2){	
		//read test particles
		er = H.readFile(binfile);
		printf("read file OK\n");
		fclose(binfile);

		/*					
		// -----------------------------------
		// Use this to extract a single object
		int ii = 29;//166;//29; //84;
		//int ii = 83;//166;//29; //84;
		H.N = H.Nperturbers + 1;
		int n = H.Nperturbers;

		H.id_h[n] = H.id_h[ii];
		H.x_h[n] = H.x_h[ii];
		H.y_h[n] = H.y_h[ii];
		H.z_h[n] = H.z_h[ii];
		H.vx_h[n] = H.vx_h[ii];
		H.vy_h[n] = H.vy_h[ii];
		H.vz_h[n] = H.vz_h[ii];
		H.A1_h[n] = H.A1_h[ii];
		H.A2_h[n] = H.A2_h[ii];
		H.A3_h[n] = H.A3_h[ii];
		H.jd_init_h[n] = H.jd_init_h[ii];
		
printf("xyz %.40g %.40g %.40g %.40g %.40g %.40g %.40g %.40g %.40g %.20g %llu\n", H.x_h[n], H.y_h[n], H.z_h[n], H.vx_h[n], H.vy_h[n], H.vz_h[n], H.A1_h[n], H.A2_h[n], H.A3_h[n], H.jd_init_h[n], H.id_h[n]);
		
		// -----------------------------------
		*/
		if(er == 0){
			return 0;
		}
	}
	if(H.useFIFO == 0){

		int er = H.Host::readIC();
		if(er == 0){
			return 0;
		}


	}


	//convert velocities and nonGrav terms
	H.convertV();
	// **************************************************

	// **************************************************
	//copy the data to the device
	if(H.useGPU > 0){
		H.copy1();
	}	
	// **************************************************

	// *******************************************************************
	// Allocate and set parameters for the Runge-Kutta-Fehlberg integrator
	// *******************************************************************
	H.Alloc2();
	H.Initialize2();


	if(RKFn == 6){
		H.setRKF45();
	}
	else if(RKFn == 7){
		H.setDP54();
	}
	else if(RKFn == 13){
		H.setRKF78();
	}
	else{
		printf("RKFn values not valid %d\n", RKFn);
		return 0;
	}

	if(H.useGPU > 0){
		H.copyConst();
	}
	// *******************************************************************

	cudaDeviceSynchronize();
	error = cudaGetLastError();
	printf("copy error = %d = %s\n",error, cudaGetErrorString(error));
	if(error != 0.0){
		return 0;
	}

        cudaEventRecord(tt2);
        cudaEventSynchronize(tt2);
        cudaEventElapsedTime(&milliseconds, tt1, tt2);
	printf("Time for ic and allocation, %g seconds\n", milliseconds * 0.001);
        timing[0] += milliseconds * 0.001;


	cudaEventRecord(tt1);

	// **************************************************
	// Read perturbers masses from perturbers.h file
	// **************************************************

	H.perturbersMass();
	H.perturbersIDs();
//m[Nperturbers] = 1.e-11; //ca mass of Flora

	// **************************************************


	// **************************************************
	//perturbers table
	// **************************************************
	
	er = H.readTable();
	if(er == 0){
		return 0;
	}
	// **************************************************

	error = cudaGetLastError();
	printf("Perturbers error = %d = %s\n",error, cudaGetErrorString(error));
	if(error != 0.0){
		return 0;
	}

	cudaEventRecord(tt2);
	cudaEventSynchronize(tt2);
	
	cudaEventElapsedTime(&milliseconds, tt1, tt2);
	printf("Time for perturbers table, %g seconds\n", milliseconds * 0.001);
        timing[1] += milliseconds * 0.001;

	cudaEventRecord(tt1);


	


	H.dtfile = fopen(H.dtfilename, "w");

	//###########################################
	// Start pre-integration
	//###########################################

	if(H.DoPreIntegration == 1){
		double dtiOld = H.dti;
		printf("dtiOld %g\n", dtiOld);

		er = H.preIntegration();
		if(er == 0){
			return 0;
		}
		H.dti = dtiOld;
		H.dt = H.dti * dayUnit;
		H.time = H.outStart;
	}

	//###########################################
	// End pre-integration
	//###########################################


	cudaEventRecord(tt2);
	cudaEventSynchronize(tt2);
	
	cudaEventElapsedTime(&milliseconds, tt1, tt2);
	printf("Time for pre-integration, %g seconds\n", milliseconds * 0.001);
        timing[2] += milliseconds * 0.001;

	cudaEventRecord(tt1);


	//###########################################
	// First output
	//###########################################

	H.output(0, H.time, 0);


	//save coordinates for repeated integrations
	H.Initialize3();

	printf("dti %g\n", H.dti);

	//###########################################
	// Start time step loop
	//###########################################

	H.Sn[0] = H.N;
for(int S = 0; S < 3; ++S){

	//###########################################
	//Time step loop
	//###########################################
	H.IntegrationLoop(S);

	cudaDeviceSynchronize();

	cudaEventRecord(tt2);
	cudaEventSynchronize(tt2);
	
	cudaEventElapsedTime(&milliseconds, tt1, tt2);
        timing[3 + S] += milliseconds * 0.001;

	printf("Time for integration %d, %g seconds\n", S + 1, timing[3 + S]);

	cudaEventRecord(tt1);

	//###########################################
	// End time step loop
	//###########################################

	//reduce arrays for repeated integration with a smaller time step
	H.reduce(S);
}
	fclose(H.dtfile);
	
	printf("Time for ic and allocation, %g seconds\n", timing[0]);
	printf("Time for perturbers table, %g seconds\n", timing[1]);
	printf("Time for pre-integration, %g seconds\n", timing[2]);
	printf("Time for integration 1, %g seconds\n", timing[3]);
	printf("Time for integration 2, %g seconds\n", timing[4]);
	printf("Time for integration 3, %g seconds\n", timing[5]);

	for(int i = 0; i < 4; ++i){
		printf("N %d %d\n", i, H.Sn[i]-Nperturbers);
	}	


}
	

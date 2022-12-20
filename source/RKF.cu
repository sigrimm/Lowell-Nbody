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

	int er = H.readparam(argc, argv);
	if(er <= 0){
		printf("Error in reading param.dat file\n");
		return 0;
	}
	printf("Reading param.dat file OK\n");

	//Erase Outbinary file
	if(H.outBinary == 1){
		if(H.outHelio == 1){
			sprintf(H.outfilename, "Outhelio.bin");
		}
		else{
			sprintf(H.outfilename, "Outbary.bin");
		}
		H.outfile = fopen(H.outfilename, "wb");
		fclose(H.outfile);
	}


	if(H.useFIFO == 2){	
		//sprintf(H.infilename, "210921_2148_genga_in_yarkovsky_elements.bin");
		//sprintf(H.infilename, "211208_1916_genga_in_2021-12-08_specific_desig.bin");
		//sprintf(H.infilename, "210801_2104_genga_in_GA.bin");
		//sprintf(H.infilename, "220301_2048_genga_in_new_last_14_days.bin");
		//sprintf(H.infilename, "220524_2256_genga_in_query_genga_input_40k.bin");
		//sprintf(H.infilename, "220524_2258_genga_in_query_genga_input_10k.bin");
	}
	
	//read console arguments for the binary file name
	//other arguments are checked later to overwright the head data	
	for(int i = 1; i < argc; i += 2){
		if(strcmp(argv[i], "-in") == 0){
			sprintf(H.infilename, "%s", argv[i + 1]);
		}
	}

	double timing[6];
	for(int i = 0; i < 6; ++i){
		timing[i] = 0.0;
	}
	float milliseconds = 0.0f;
	cudaError_t error;

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
		H.infile = fopen(H.infilename, "rb");

		if(H.infile == NULL){
			printf("Error, input file not found %s\n", H.infilename);
			return 0;
		}

		int NTP = 0;
		printf("read file\n");

		er = H.readHeader(H.infile, NTP);
		if(er == 0){
			return 0;
		}

		//H.NTP = 1;

		H.N += NTP;
	}

	if(H.useFIFO == 0){
		int n = H.readICSize();
		if(n == 0){
			printf("Error, reading initial conditions file failed\n");
			return 0;
		}
		H.N += n;
	}

	if(H.N >= 1024 * 1024){
		printf("Error, N is too large for scan kernels\n");
		return 0;
	}


	H.time = H.time0;
	// **************************************************
//H.time1 = 2451000.5;
//H.outInterval = 10.0;
//H.outStart = 2450800.5;

	//move this to parameter file
	H.dti = 2.0;



	//read console arguments
	for(int i = 1; i < argc; i += 2){
		if(strcmp(argv[i], "-Nsteps") == 0){
			H.Nsteps = atoll(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-outInterval") == 0){
			H.outInterval = atoll(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-outStart") == 0){
			H.outStart = atof(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-dt") == 0){
			H.dti = atof(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-dtMin") == 0){
			H.dtiMin[0] = atof(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-endTime") == 0){
			H.time1 = atof(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-useGPU") == 0){
			H.useGPU = atoi(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-useFIFO") == 0){
			H.useFIFO = atoi(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-outBinary") == 0){
			H.outBinary = atoi(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-useAdaptive") == 0){
			H.useAdaptiveTimeSteps = atoi(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-useIndividual") == 0){
			H.useIndividualTimeSteps = atoi(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-NMax") == 0){
			H.NMax = atoi(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-in") == 0){
			//this is already done
		}
		else{
			printf("Error, console argument not valid.\n");
			return 0;
		}
	}
	H.dt = H.dti * dayUnit;


	if(H.N - Nperturbers > H.NMax){
		printf("Number of bodies larger than Nmax\n");
		H.N = Nperturbers + H.NMax;
	}

	printf("infile name %s\n", H.infilename);

	printf("outStart: %.20g, time0: %.20g, time1: %.20g, outInterval: %lld\n", H.outStart, H.time0, H.time1, H.outInterval);

	printf("Nperturbers: %d N: %d\n", Nperturbers, H.N);


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

	H.initialize1();


	// **************************************************
	// Read initial conditions
	// **************************************************

	printf("Read initial conditions\n");
	if(H.useFIFO == 2){	
		//read particles
		er = H.readFile(H.infile);
		printf("read file OK\n");
		fclose(H.infile);

		/*						
		// -----------------------------------
		// Use this to extract a single object
		int ii = 29;//166;//29; //84;
		//int ii = 83;//166;//29; //84;
		H.N = Nperturbers + 1;
		int n = Nperturbers;

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
	// Read perturbers masses from perturbers.h file
	// **************************************************

	H.perturbersMass();
	H.perturbersIDs();
//m[Nperturbers] = 1.e-11; //ca mass of Flora

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
	H.initialize2();


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
	else{
		H.time = H.time0;
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


	//save coordinates for repeated integrations
	H.initialize3();

	unsigned long long int outCount = 0;
	//###########################################
	// First output
	//###########################################

	H.output(outCount, H.time);
	++outCount;


	H.runsN[0] = H.N;

	double timeb = H.time;
	double time1b = H.time1;
	double outStartb = H.outStart;
	int N0 = H.N;


	//###########################################
	//loop for forward backward integration
	//###########################################
	//for(int b = 0; b < 1; ++b){	
	for(int b = 0; b < 2; ++b){	

		int DT0 = 2;
		int DT = DT0;
		unsigned long long int outI = H.outInterval;
		unsigned long long int cOut = 0llu;

		if(b == 0 && H.time1 <= H.time0){
			//skip forward integration
			continue;
		}


		if(b == 1){
			//backward integration
			H.restore3();
			H.time1 = outStartb;
			H.time = timeb;
			H.outStart = time1b;
			DT0 = -DT0;
			DT = -DT;
		}

		H.time0 = timeb;


		if(DT > 0 && H.outStart > H.time){
			outI = H.outStart - H.time;
		}
		if(DT < 0 && H.outStart < H.time){
			outI = H.time - H.outStart;
printf("outI %.20g %.20g\n", H.time, H.outStart);
		}


		//###########################################
		// Start time step loop
		//###########################################
		for(int tt = 0; tt < 100000; ++tt){


			if(cOut + fabs(DT) > outI && ((H.dt > 0 && H.time0 + H.dti >= H.outStart) || (H.dt < 0 && H.time + H.dti <= H.outStart))){
				DT = (outI - cOut);
				if(b == 1) DT = -DT;

printf("   correctT %.20g %.20g %d %llu %llu\n", H.time0, H.time0 + DT, DT, cOut, outI);
			}


printf("D %.20g %.20g %llu %llu %d\n", H.time0, H.time1, cOut, outI, DT);
			// loop for different time step ranges
			H.N = N0;
			H.save1();
			double time1 = H.time0 + DT;

			//###########################################
			// Loop over different time step levels
			//###########################################
			for(int S = 0; S < H.nRuns; ++S){
			//for(int S = 0; S < 1; ++S){

				H.dti = fabs(H.runsdt[S]);
				double dtmin = H.dti;

				if(fabs(H.dti) > fabs(DT)) H.dti = fabs(DT);

				//set dts, round dti down to power of 10
				{
					double l = log10(H.dti);
					double f = floor(l);
					double s = pow(10.0, f);
					H.dts = s * 0.1;
//printf("H.dti %g %g %g %g %g | %g %g\n", H.dti, l, f, s, H.dts, H.dti, H.dtiMin[S]);

					int dtt;
					if(H.dti >= 0.0){
					dtt = (fabs(H.dti) + (0.5 * H.dts)) / H.dts;
					}
					else{
					dtt = -(fabs(H.dti) - (0.5 * H.dts)) / H.dts;
					}
					H.dti = dtt * H.dts;
//printf("H.dtiB %.20g %d %g\n", H.dti, dtt, H.dts);


				}
				unsigned long long int Nci = (fabs(DT) + 0.5 * H.dts) / H.dts;

				if(b == 1){
					//backward integration
					H.dti = -H.dti;
				}


				H.time = H.time0;

				H.setSnew();


				//###########################################
				//Time step 
				//###########################################
				H.IntegrationLoop(S, Nci, time1, dtmin);

				cudaDeviceSynchronize();

				cudaEventRecord(tt2);
				cudaEventSynchronize(tt2);
				
				cudaEventElapsedTime(&milliseconds, tt1, tt2);
				timing[3 + S] += milliseconds * 0.001;

				printf("Time for integration %d, %g seconds\n", S + 1, timing[3 + S]);
				printf("With %d bodies\n", H.N);

				cudaEventRecord(tt1);

				//###########################################
				// End time step
				//###########################################

				H.save(dtmin);

				//reduce arrays for repeated integration with a smaller time step
				if(H.useGPU > 0){
					H.reduceCall(S);
				}
				else{
					H.reduce(S);
				}

				if(H.N == Nperturbers){
					break;
				}

			} //end of S loop
			H.N = N0;
			H.time0 = H.time;
			cOut += fabs(DT);


			if(cOut >= outI && ((DT > 0 && H.time >= H.outStart) || (DT < 0 && H.time <= H.outStart))){
			//if(t % 10 == 0){

				H.output(outCount, H.time);
				++outCount;
				if(outCount >= 1e12){
					printf("Error, too many output files, %llu \n", outCount);
					return 0;
				}

				DT = DT0;
				cOut = 0llu;
				outI = H.outInterval; //needed only at the first time

			}


			if((DT > 0 && H.time >= H.time1) || (DT < 0 && H.time <= H.time1)){
				printf("Reached the end\n");
				break;
			}


		} // end of tt loop
		//###########################################
		// End time step loop
		//###########################################


	} // end of b loop

	fclose(H.dtfile);
	
	printf("Time for ic and allocation, %g seconds\n", timing[0]);
	printf("Time for perturbers table, %g seconds\n", timing[1]);
	printf("Time for pre-integration, %g seconds\n", timing[2]);
	printf("Time for integration 1, %g seconds\n", timing[3]);
	printf("Time for integration 2, %g seconds\n", timing[4]);
	printf("Time for integration 3, %g seconds\n", timing[5]);

	for(int i = 0; i < 4; ++i){
		printf("N %d %d %g\n", i, H.runsN[i] - Nperturbers, H.dtiMin[i]);
	}	

	FILE *timefile;
	timefile = fopen("Timing.dat", "a");
	//for(int i = 0; i < 4; ++i){
	for(int i = 0; i < 1; ++i){
		fprintf(timefile, "%d %d %g %g\n", i, H.runsN[i] - Nperturbers, H.dtiMin[i], timing[3 + i]);
	}
	fclose(timefile);
	
	return 0;

}
	

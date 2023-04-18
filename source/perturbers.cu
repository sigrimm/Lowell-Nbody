#include "Host.h"
__host__ int Host::perturbersMass(){

	FILE *pFile;

	pFile = fopen("perturbers.dat", "r");
	if(pFile == NULL){
		printf("Error, perturbers.dat file not found\n");
		return 0;
	}

	char name[160];
	double im;
	long long int id1;
	long long int id2;


	//read header
	fscanf(pFile, "%[^\n]", name);
//printf("%s\n", name);

	//The order must agree with the perturbers coordinate file
	//Units are 1/(mass of object in solar masses)

	int er = 0;
	for(int i = 0; i < 1000; ++i){
		er = fscanf(pFile, "%lf", &im);
		er = fscanf(pFile, "%s", name);
		er = fscanf(pFile, "%lld", &id1);
		er = fscanf(pFile, "%lld", &id2);


		if(er <= 0){

			if(i != Nperturbers){
				printf("Error, perturbers.dat contains not 'Nperturbers' objects, %d %d\n", i, Nperturbers);
				return 0;

			}
			break;
		}
		m_h[i] = 1.0 / im;
		id_h[i] = id2;
printf("m %d %.20g %lld\n", i, m_h[i], id_h[i]);


	}

	fclose(pFile);

	return 1;
}



__global__ void bufferToX_kernel(double *XYdata_d, double *timep_d, double *xp_d, double *yp_d, double *zp_d, int N){

	int id = threadIdx.x + blockDim.x * blockIdx.x;

	if(id < N){
		timep_d[id] = XYdata_d[id * 4];
		xp_d[id] = XYdata_d[id * 4 + 1];
		yp_d[id] = XYdata_d[id * 4 + 2];
		zp_d[id] = XYdata_d[id * 4 + 3];
if(id < 30 || id > N - 10 ) printf("buffer %d %.20g %.20g %.20g %.20g\n", id, timep_d[id], xp_d[id], yp_d[id], zp_d[id]);
	}
}
__host__ void bufferToX(double *XYdata_h, double *timep_h, double *xp_h, double *yp_h, double *zp_h, int N){


	for(int id = 0; id < N; ++id){
		timep_h[id] = XYdata_h[id * 4];
		xp_h[id] = XYdata_h[id * 4 + 1];
		yp_h[id] = XYdata_h[id * 4 + 2];
		zp_h[id] = XYdata_h[id * 4 + 3];
if(id < 30 || id > N - 10) printf("buffer %d %.20g %.20g %.20g %.20g\n", id, timep_h[id], xp_h[id], yp_h[id], zp_h[id]);
	}
}



//read the perturbers table
__host__ int Host::readTable(){

	char *pfilename;
	pfilename = new char[320];
	if(useHelio == 1){
		sprintf(pfilename, "%s/All3_h.bin", perturbersPath);
		XVfile = fopen(pfilename, "rb");
printf("Read Heliocentric perturbers\n");
	}
	else{
		sprintf(pfilename, "%s/All3_b.bin", perturbersPath);
		XVfile = fopen(pfilename, "rb");
printf("Read Barycentric perturbers\n");
	}

	if(XVfile == NULL){
		printf("Error, perturber file not found\n");
		return 0;
	}

	int NTableC = NTable;
	//use two data buffers for interleaved data transfer
	for(int t = 0; t < NTable; ++t){

	//for(int t = 0; t < 10; ++t){
//printf("t %d\n", t);
		int er;

		if(t % 2 == 0){
//printf("start read A\n");
			er = fread(readBufferA_h, Nperturbers * 4 * sizeof(double), 1, XVfile);
//printf("end read A\n");
		}
		else{
//printf("start read B\n");
			er = fread(readBufferB_h, Nperturbers * 4 * sizeof(double), 1, XVfile);
//printf("end read B\n");
		}

		/*
		//only here for checking
		if(t < 4){
			for(int i = 0; i < Nperturbers; ++i){
				if(t % 2 == 0) printf("XYa %d %d %.20g %g %g %g\n", t, i, readBufferA_h[i * 4 + 0], readBufferA_h[i * 4 + 1], readBufferA_h[i * 4 + 2], readBufferA_h[i * 4 + 3]);
				if(t % 2 == 1) printf("XYb %d %d %.20g %g %g %g\n", t, i, readBufferB_h[i * 4 + 0], readBufferB_h[i * 4 + 1], readBufferB_h[i * 4 + 2], readBufferB_h[i * 4 + 3]);
			}
		}
		*/

		if(t == 0){
			//set start time of perturbers file
			timep0 = readBufferA_h[0];
		}

		if(useGPU == 0){
			if(t % 2 == 0){
				memcpy(XYdata_h + t * Nperturbers * 4, readBufferA_h, Nperturbers * 4 * sizeof(double));
			}
			else{
				memcpy(XYdata_h + t * Nperturbers * 4, readBufferB_h, Nperturbers * 4 * sizeof(double));
			}
		}
		else{
			cudaDeviceSynchronize(); //this must be here

			//both buffers A and B use the same stream, so copy can overlap with the next read, but not with the
			//next copy. 
			if(t % 2 == 0){
//printf("start copy A\n");
				cudaMemcpyAsync(XYdata_d + t * Nperturbers * 4, readBufferA_h, Nperturbers * 4 * sizeof(double), cudaMemcpyHostToDevice);
			}
			else{
//printf("start copy B\n");
				cudaMemcpyAsync(XYdata_d + t * Nperturbers * 4, readBufferB_h, Nperturbers * 4 * sizeof(double), cudaMemcpyHostToDevice);
			}
		}

		if(er <= 0){
//printf("readbuffer %d %d %d %.20g %g %g %g\n", er, t, bSwap, readBuffer_h[0], readBuffer_h[1], readBuffer_h[2], readBuffer_h[3]);
			NTableC = t;
			break;
		}

//if(t < 4) printf("readbuffer %d %d %d %.20g %g %g %g\n", er, t, bSwap, readBuffer_h[0], readBuffer_h[1], readBuffer_h[2], readBuffer_h[3]);
	}
	timep1 = fmax(readBufferA_h[0], readBufferB_h[0]);

	printf("NTableC: %d %.20g %.20g | %.20g %.20g\n", NTableC, timep0, timep1, time0, time1);

	if(time0 >= timep1 || time1 >= timep1){
		printf("Error, NTable is not large enough to store perturbers table for required time period.\n");
		return 0;
	}


	fclose(XVfile);
	cudaDeviceSynchronize();

	if(useGPU == 0){
		bufferToX (XYdata_h, timep_h, xp_h, yp_h, zp_h, NTableC * Nperturbers);
		free(XYdata_h);
		free(readBufferA_h);
		free(readBufferB_h);
	}
	else{
		bufferToX_kernel <<< (NTableC * Nperturbers + 127) / 128, 128 >>> (XYdata_d, timep_d, xp_d, yp_d, zp_d, NTableC * Nperturbers);
		cudaFree(XYdata_d);
		cudaFreeHost(readBufferA_h);
		cudaFreeHost(readBufferB_h);
	}

	return 1;
}

//needed with useGPU > 0 but still using the CPU for the preintegration
__host__ void Host::copyPerturbersCPU(){
	cudaMemcpy(timep_h, timep_d, Nperturbers * NTable * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(xp_h, xp_d, Nperturbers * NTable * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(yp_h, yp_d, Nperturbers * NTable * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(zp_h, zp_d, Nperturbers * NTable * sizeof(double), cudaMemcpyDeviceToHost);
}

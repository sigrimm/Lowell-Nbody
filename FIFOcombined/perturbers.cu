#include "Host.h"
__host__ void Host::perturbersMass(){

	//Units are 1/(mass of object in solar masses)
	double pmass[] = {
		//This values come from the Larry Wasserman aiset.f
		//From DE-440
		//The order must agree with the perturbers coordinate file
		1.000000000000000e0,      	// Sun        10

		6.02365794492933620e+06,	// Mercury    1
		4.08523718656268204e+05,	// Venus      2
		3.32946048773048213e+05,	// Earth      399
		3.09870354673725273e+06,	// Mars       4
		1.04734863124400454e+03,	// Jupiter    5
		3.49790180079320044e+03,	// Saturn     6
		2.29029507834766191e+04,	// Uranus     7
		1.94122597758754673e+04,	// Neptune    8
		1.36045556167380244e+08,	// Pluto      9

		2.11902913996742201e+09,	// Ceres      1
		9.71122664959812355e+09,	// Pallas     2
		6.91005231035131531e+10,	// Juno       3
		7.67646068680440331e+09,	// Vesta      4
		2.35927034111103439e+10,	// Hygiea     10
		6.56011187658897781e+10,	// Eunomia    15
		1.22953445817123718e+11,	// Euphrosyne 31
		4.94635345133462448e+10,	// Europa     52
		3.40770353836609001e+10,	// Davida     511
		4.68880681429191055e+10,	// Interamnia 704
		8.34848877284898376e+10,	// Psyche     16
		1.41468527548287018e+11,	// Cybele     65
		1.11541082696269424e+11,	// Thisbe     88
		9.19227747550098419e+10,	// Camilla    107
		1.16427460635738602e+11,	// Iris       7
		6.12076731319770355e+10,	// Sylvia     87

		2.70687029523511454e+07,	// Moon       301
	};
	for(int i = 0; i < Nperturbers; ++i){
		m_h[i] = 1.0/pmass[i];
printf("m %d %.20g\n", i, m_h[i]);
//printf("m %d %.20g %.20g\n", i, m_h[i], mass1[i] / mass1[0]);
	}

}


__host__ void Host::perturbersIDs(){
	//perturbers indices
	unsigned long long int id[] = {
		10,			//sun

		1,
		2,
		399,
		4,
		5,
		6,
		7,
		8,
		72057594038083239,      // 9,
		72057594037948899,	//Ceres      1
		72057594037948900,	//Pallas     2
		72057594037948901,	//Juno       3
		72057594037948902,	//Vesta      4
		72057594037948908,	//Hygiea     10
		72057594037948913,	//Eunomia    15
		72057594037948929,	//Euphrosyne 31
		72057594037948950,	//Europa     52
		72057594037949409,	//Davida     511
		72057594037949602,	//Interamnia 704
		72057594037948914,	//Psyche     16
		72057594037948963,	//Cybele     65
		72057594037948986,	//Thisbe     88
		72057594037949005,	//Camilla    107
		72057594037948905,	//Iris       7
		72057594037948985,	//Sylvia     87
		301,			//Moon
	};

	for(int i = 0; i < Nperturbers; ++i){
		id_h[i] = id[i];
	}
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

	if(useHelio == 1){
		XVfile = fopen("All3_h.bin", "rb");
	}
	else{
		XVfile = fopen("All3_b.bin", "rb");
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

#include "asteroid.h"


int asteroid::allocateGPU(){

	//Check warp size
	cudaDeviceProp devProp;
	int dev = 0;					//Set device number
	cudaGetDeviceProperties(&devProp, dev);
	WarpSize = devProp.warpSize;


	cudaMalloc((void **) &startTime_d, Nperturbers * RKFn * sizeof(double));
	cudaMalloc((void **) &endTime_d, Nperturbers * RKFn * sizeof(double));
	cudaMalloc((void **) &id_d, Nperturbers * sizeof(int));
	cudaMalloc((void **) &nChebyshev_d, Nperturbers * sizeof(int));
	cudaMalloc((void **) &offset0_d, Nperturbers * RKFn * sizeof(int));
	cudaMalloc((void **) &offset1_d, Nperturbers * RKFn * sizeof(int));
	cudaMalloc((void **) &GM_d, Nperturbers * sizeof(double));

	cudaMalloc((void **) &cdata_d, Nperturbers * RKFn * nCm * 3 * sizeof(double));
	cudaMalloc((void **) &data_d, datasize * sizeof(double));

	cudaMalloc((void **) &xTable_d, Nperturbers * RKFn * sizeof(double));
	cudaMalloc((void **) &yTable_d, Nperturbers * RKFn * sizeof(double));
	cudaMalloc((void **) &zTable_d, Nperturbers * RKFn * sizeof(double));

	cudaMalloc((void **) &vxTable_d, Nperturbers * RKFn * sizeof(double));
	cudaMalloc((void **) &vyTable_d, Nperturbers * RKFn * sizeof(double));
	cudaMalloc((void **) &vzTable_d, Nperturbers * RKFn * sizeof(double));

	cudaMalloc((void **) &x_d, N * sizeof(double));
	cudaMalloc((void **) &y_d, N * sizeof(double));
	cudaMalloc((void **) &z_d, N * sizeof(double));

	cudaMalloc((void **) &vx_d, N * sizeof(double));
	cudaMalloc((void **) &vy_d, N * sizeof(double));
	cudaMalloc((void **) &vz_d, N * sizeof(double));

	cudaMalloc((void **) &dx_d, N * sizeof(double));
	cudaMalloc((void **) &dy_d, N * sizeof(double));
	cudaMalloc((void **) &dz_d, N * sizeof(double));

	cudaMalloc((void **) &dvx_d, N * sizeof(double));
	cudaMalloc((void **) &dvy_d, N * sizeof(double));
	cudaMalloc((void **) &dvz_d, N * sizeof(double));

	cudaMalloc((void **) &kx_d, N * RKFn * sizeof(double));
	cudaMalloc((void **) &ky_d, N * RKFn * sizeof(double));
	cudaMalloc((void **) &kz_d, N * RKFn * sizeof(double));

	cudaMalloc((void **) &kvx_d, N * RKFn * sizeof(double));
	cudaMalloc((void **) &kvy_d, N * RKFn * sizeof(double));
	cudaMalloc((void **) &kvz_d, N * RKFn * sizeof(double));

	cudaMalloc((void **) &ax_d, N * sizeof(double));
	cudaMalloc((void **) &ay_d, N * sizeof(double));
	cudaMalloc((void **) &az_d, N * sizeof(double));

	cudaMalloc((void **) &A1_d, N * sizeof(double));
	cudaMalloc((void **) &A2_d, N * sizeof(double));
	cudaMalloc((void **) &A3_d, N * sizeof(double));


	cudaMalloc((void **) &snew_d, N * sizeof(double));
	cudaMalloc((void **) &ssum_d, N * sizeof(double));

	cudaDeviceSynchronize();
	cudaError_t error = cudaGetLastError();
	printf("allocate error = %d = %s\n",error, cudaGetErrorString(error));
	if(error != 0.0){
		return 0;
	}
	return 1;
}



int asteroid::readData(){

	int er;

	er = fread(data_h, sizeof(double), datasize, perturbersFile);

	if(er <= 0){
		return 0;
	}

	/*
	for(int i = 0; i < 20; ++i){
		printf("%.20g ", data_h[i]);
	}
	printf("\n");
	for(int i = datasize - 20; i < datasize; ++i){
		printf("%.20g ", data_h[i]);
	}
	printf("\n");
	*/
	
	return 1;
}


int asteroid::copyIC(){

	cudaMemcpy(startTime_d, startTime_h, Nperturbers * RKFn * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(endTime_d, endTime_h, Nperturbers * RKFn * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(id_d, id_h, Nperturbers * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(nChebyshev_d, nChebyshev_h, Nperturbers * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(offset0_d, offset0_h, Nperturbers * RKFn * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(offset1_d, offset1_h, Nperturbers * RKFn * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(GM_d, GM_h, Nperturbers * sizeof(double), cudaMemcpyHostToDevice);

	cudaMemcpy(data_d, data_h, datasize * sizeof(double), cudaMemcpyHostToDevice);

	cudaMemcpy(x_d, x_h, N * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(y_d, y_h, N * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(z_d, z_h, N * sizeof(double), cudaMemcpyHostToDevice);

	cudaMemcpy(vx_d, vx_h, N * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(vy_d, vy_h, N * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(vz_d, vz_h, N * sizeof(double), cudaMemcpyHostToDevice);

	cudaMemcpy(A1_d, A1_h, N * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(A2_d, A2_h, N * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(A3_d, A3_h, N * sizeof(double), cudaMemcpyHostToDevice);

	cudaDeviceSynchronize();
	cudaError_t error = cudaGetLastError();
	printf("copy error = %d = %s\n",error, cudaGetErrorString(error));
	if(error != 0.0){
		return 0;
	}

	return 1;
}

void asteroid::copyOutput(){

	cudaMemcpy(x_h, x_d, N * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(y_h, y_d, N * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(z_h, z_d, N * sizeof(double), cudaMemcpyDeviceToHost);

	cudaMemcpy(vx_h, vx_d, N * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(vy_h, vy_d, N * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(vz_h, vz_d, N * sizeof(double), cudaMemcpyDeviceToHost);

}


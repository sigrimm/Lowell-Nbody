#include "asteroid.h"



int asteroid::allocateGPU(){

	cudaMalloc((void **) &x_d, N * sizeof(double));
	cudaMalloc((void **) &y_d, N * sizeof(double));
	cudaMalloc((void **) &z_d, N * sizeof(double));

	cudaMalloc((void **) &vx_d, N * sizeof(double));
	cudaMalloc((void **) &vy_d, N * sizeof(double));
	cudaMalloc((void **) &vz_d, N * sizeof(double));

	cudaMalloc((void **) &xt_d, N * sizeof(double));
	cudaMalloc((void **) &yt_d, N * sizeof(double));
	cudaMalloc((void **) &zt_d, N * sizeof(double));

	cudaMalloc((void **) &vxt_d, N * sizeof(double));
	cudaMalloc((void **) &vyt_d, N * sizeof(double));
	cudaMalloc((void **) &vzt_d, N * sizeof(double));

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


	cudaMalloc((void **) &a_d, RKFn * RKFn * sizeof(double));
	cudaMalloc((void **) &b_d, RKFn * sizeof(double));
	cudaMalloc((void **) &bb_d, RKFn * sizeof(double));
	cudaMalloc((void **) &c_d, RKFn * sizeof(double));

	cudaDeviceSynchronize();
	cudaError_t error = cudaGetLastError();
	printf("allocate error = %d = %s\n",error, cudaGetErrorString(error));
	if(error != 0.0){
		return 0;
	}
	return 1;
}



int asteroid::copyIC(){

	cudaMemcpy(x_d, x_h, N * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(y_d, y_h, N * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(z_d, z_h, N * sizeof(double), cudaMemcpyHostToDevice);

	cudaMemcpy(vx_d, vx_h, N * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(vy_d, vy_h, N * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(vz_d, vz_h, N * sizeof(double), cudaMemcpyHostToDevice);

	cudaMemcpy(A1_d, A1_h, N * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(A2_d, A2_h, N * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(A3_d, A3_h, N * sizeof(double), cudaMemcpyHostToDevice);

	cudaMemcpy(a_d, a_h, RKFn * RKFn * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(b_d, b_h, RKFn * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(bb_d, bb_h, RKFn * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(c_d, c_h, RKFn * sizeof(double), cudaMemcpyHostToDevice);

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


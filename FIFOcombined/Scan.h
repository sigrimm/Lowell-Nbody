//**************************************
//This kernel performs a scan operation, used for  stream compactions
//
//It works for the case of multiple blocks
//must be followed by Scan32d2 and Scan32d3
//
//Uses shuffle instructions
//Authors: Simon Grimm
//June 2022
//  *****************************************
__global__ void Scan32d1_kernel(double2 *snew_d, int2 *scan_d, unsigned int *index_d, unsigned int *index0_d, const int N){

	int idx = threadIdx.x;
	int id = blockIdx.x * blockDim.x + idx;

	int t1 = 0;
	int t2 = 0;

	extern __shared__ int Scand1_s[];
	int *t_s = Scand1_s;

	int lane = threadIdx.x % warpSize;
	int warp = threadIdx.x / warpSize;

	if(warp == 0){
		t_s[threadIdx.x] = 0;
	}

	if(id < N){
		if(snew_d[id].y < 1.0) t1 = 1;
		index0_d[id] = index_d[id];
	}
	__syncthreads();
//if(id < 128) printf("Scan a %d %d %d\n", id, idx, t1);

	for(int i = 1; i < warpSize; i*=2){
#if def_OldShuffle == 0
		t2 = __shfl_up_sync(0xffffffff, t1, i, warpSize);
#else
		t2 = __shfl_up(t1, i);
#endif
		if(idx % warpSize >= i) t1 += t2;
	}
	__syncthreads();

	int t0 = t1;

	if(blockDim.x > warpSize){
	//reduce across warps

		if(lane == warpSize - 1){
			t_s[warp] = t1;
//if(blockIdx.x == 0) printf("w %d %d\n", warp, t1);
		}

		__syncthreads();
		//reduce previous warp results in the first warp
		if(warp == 0){
			t1 = t_s[threadIdx.x];
			for(int i = 1; i < warpSize; i*=2){
#if def_OldShuffle == 0
				t2 = __shfl_up_sync(0xffffffff, t1, i, warpSize);
#else
				t2 = __shfl_up(t1, i);
#endif
				if(lane >= i) t1 += t2;
			}
		}
		if(idx < blockDim.x / warpSize){
			t_s[idx] = t1;
		}

		__syncthreads();

		if(idx >= warpSize){
			t0 += t_s[warp - 1];
		}
	}
	__syncthreads();
//if(id == 1874) printf("Scan C %d %d %d | %d\n", id, idx, t0, t0 + Nperturbers - 1);

	if(id < N){
		scan_d[id].x = t0;
	}

	if(idx == blockDim.x - 1){
		scan_d[blockIdx.x].y = t0;
//printf("ScanD %d %d\n", blockIdx.x, t0);
	}

}

//**************************************
//This kernel reads the result from the multiple thread block kernel Scan32d1
//and performs the last summation step in
// --a single thread block --
//
//must be followed by Scan32d3
//
//Uses shuffle instructions
//Authors: Simon Grimm
//June 2022
//  *****************************************
__global__ void Scan32d2_kernel(int2 *scan_d, const int N){

	int idx = threadIdx.x;

	int t1 = 0;
	int t2 = 0;

	extern __shared__ int Scand2_s[];
	int *t_s = Scand2_s;

	int lane = threadIdx.x % warpSize;
	int warp = threadIdx.x / warpSize;

	if(warp == 0){
		t_s[threadIdx.x] = 0;
	}

	t1 = scan_d[idx].y;
	if(t1 < 0) t1 = 0;

	__syncthreads();
//if(idx < 32) printf("Scan a %d %d\n", idx, t1);

	for(int i = 1; i < warpSize; i*=2){
#if def_OldShuffle == 0
		t2 = __shfl_up_sync(0xffffffff, t1, i, warpSize);
#else
		t2 = __shfl_up(t1, i);
#endif
		if(idx % warpSize >= i) t1 += t2;
	}
	__syncthreads();
//if(idx < 32) printf("Scan b %d %d\n", idx, t1);

	int t0 = t1;

	if(blockDim.x > warpSize){
	//reduce across warps

		if(lane == warpSize - 1){
			t_s[warp] = t1;
		}

		__syncthreads();
		//reduce previous warp results in the first warp
		if(warp == 0){
			t1 = t_s[threadIdx.x];
			for(int i = 1; i < warpSize; i*=2){
	#if def_OldShuffle == 0
				t2 = __shfl_up_sync(0xffffffff, t1, i, warpSize);
	#else
				t2 = __shfl_up(t1, i);
	#endif
				if(lane >= i) t1 += t2;
			}
		}
		if(idx < blockDim.x / warpSize){
			t_s[idx] = t1;
		}

		__syncthreads();

		if(idx >= warpSize){
			t0 += t_s[warp - 1];
		}
	}
	__syncthreads();
//printf("Scan CC %d %d\n", idx, t0);
	if(idx < (N + 1023) / 1024){
//printf("Scan CC1 %d %d\n", idx, t0);
		scan_d[idx].y = t0;
	}
}

__global__ void Scan32d3_kernel(double2 *snew_d, int2 *scan_d, double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *A1_d, double *A2_d, double *A3_d, double *m_d, unsigned long long int *id_d, unsigned int *index_d, double *x0_d, double *y0_d, double *z0_d, double *vx0_d, double *vy0_d, double *vz0_d, double *A10_d, double *A20_d, double *A30_d, double *m0_d, unsigned long long int *id0_d, unsigned int *index0_d, int *N_d, const int N){

	int idx = threadIdx.x;
	int id = blockIdx.x * blockDim.x + idx;

	if(id < N){
		int ii = id / 1024;
		int t = scan_d[id].x;
		if(id >= 1024){
//if(id == 1874) printf("Scan E %d %d %d %d\n", id, ii, t, scan_d[ii - 1].y);
			t += scan_d[ii - 1].y;
		}
		scan_d[id].x = t;
		
		if(snew_d[id].y < 1.0){
			//reduce the arrays
			int k = scan_d[id].x - 1 + Nperturbers;
			unsigned int ii = index0_d[id];

			x_d[k] = x0_d[ii];
			y_d[k] = y0_d[ii];
			z_d[k] = z0_d[ii];
			vx_d[k] = vx0_d[ii];
			vy_d[k] = vy0_d[ii];
			vz_d[k] = vz0_d[ii];
			A1_d[k] = A10_d[ii];
			A2_d[k] = A20_d[ii];
			A3_d[k] = A30_d[ii];
			m_d[k] = m0_d[ii];
			id_d[k] = id0_d[ii];
			index_d[k] = index0_d[ii];

//if(id == 1874) printf("Scan3b %d %d %d\n", id, k, ii);
//printf("Scan3b %d %d %u %llu %.20g %.20g\n", id, k, ii, id_d[k], x_d[k], A1_d[k]);
		}

		if(id == N - 1){
			N_d[0] = t + Nperturbers;
printf("Scan3 F %d %d\n",  t, t + Nperturbers);
		}
	}
}


__host__ void Host::reduceCall(int S){
	int nct = 1024;
	int ncb = min((N + nct - 1) / nct, 1024);
	Scan32d1_kernel <<< ncb, nct, WarpSize * sizeof(int) >>> (snew_d, scan_d, index_d, index0_d, N);
	Scan32d2_kernel <<< 1, ((ncb + WarpSize - 1) / WarpSize) * WarpSize, WarpSize * sizeof(int)  >>> (scan_d, N);
	Scan32d3_kernel <<< ncb, nct >>> (snew_d, scan_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, A1_d, A2_d, A3_d, m_d, id_d, index_d, x0_d, y0_d, z0_d, vx0_d, vy0_d, vz0_d, A10_d, A20_d, A30_d, m0_d, id0_d, index0_d, N_d, N);
	cudaMemcpy(&N, N_d, sizeof(int), cudaMemcpyDeviceToHost);

printf("N %d\n", N);
	runsN[S + 1] = N;

}

#define dayUnit 0.01720209895

//check RKN12(10)17M

void setRKF45(double *a_h, double *b_h, double *bb_h, double *c_h){

	a_h[1 * 6 + 0] = 1.0/4.0;       //21

	a_h[2 * 6 + 0] = 3.0/32.0;      //31
	a_h[2 * 6 + 1] = 9.0/32.0;      //32

	a_h[3 * 6 + 0] = 1932.0/2197.0; //41
	a_h[3 * 6 + 1] = -7200.0/2197.0;//42
	a_h[3 * 6 + 2] = 7296.0/2197.0; //43

	a_h[4 * 6 + 0] = 439.0/216.0;   //51
	a_h[4 * 6 + 1] = -8.0;          //52
	a_h[4 * 6 + 2] = 3680.0/513.0;  //53
	a_h[4 * 6 + 3] = -845.0/4104.0; //54

	a_h[5 * 6 + 0] = -8.0/27.0;     //61
	a_h[5 * 6 + 1] = 2.0;           //62
	a_h[5 * 6 + 2] = -3544/2565.0;  //63
	a_h[5 * 6 + 3] = 1859.0/4104.0; //64
	a_h[5 * 6 + 4] = -11.0/40.0;    //65


	b_h[0] = 25.0/216.0;
	b_h[1] = 0.0;
	b_h[2] = 1408.0/2565.0;
	b_h[3] = 2197.0/4104.0;
	b_h[4] = -1.0/5.0;
	b_h[5] = 0.0;

	bb_h[0] = 16.0/135.0;
	bb_h[1] = 0.0;
	bb_h[2] = 6656.0/12825.0;
	bb_h[3] = 28561.0/56430.0;
	bb_h[4] = -9.0/50.0;
	bb_h[5] = 2.0/55.0;

	c_h[0] = 0.0;
	c_h[1] = 0.25;
	c_h[2] = 3.0 / 8.0;
	c_h[3] = 12.0 / 13.0;
	c_h[4] = 1.0;
	c_h[5] = 0.5;

}

void setDP54(double *a_h, double *b_h, double *bb_h, double *c_h){

	a_h[1 * 7 + 0] = 1.0/5.0;       //21

	a_h[2 * 7 + 0] = 3.0/40.0;      //31
	a_h[2 * 7 + 1] = 9.0/40.0;      //32

	a_h[3 * 7 + 0] = 44.0/45.0; //41
	a_h[3 * 7 + 1] = -56.0/15.0;//42
	a_h[3 * 7 + 2] = 32.0/9.0; //43

	a_h[4 * 7 + 0] = 19372.0/6561.0;   //51
	a_h[4 * 7 + 1] = -25360.0/2187.0;          //52
	a_h[4 * 7 + 2] = 64448.0/6561.0;  //53
	a_h[4 * 7 + 3] = -212.0/729.0; //54

	a_h[5 * 7 + 0] = 9017.0/3168.0;  //61
	a_h[5 * 7 + 1] = -355.0/33.0;           //62
	a_h[5 * 7 + 2] = 46732.0/5247.0;  //63
	a_h[5 * 7 + 3] = 49.0/176.0; //64
	a_h[5 * 7 + 4] = -5103.0/18656.0;    //65

	a_h[6 * 7 + 0] = 35.0/384.0;  //61
	a_h[6 * 7 + 1] = 0.0;           //62
	a_h[6 * 7 + 2] = 500.0/1113.0;  //63
	a_h[6 * 7 + 3] = 125.0/192.0; //64
	a_h[6 * 7 + 4] = -2187.0/6784.0;    //65
	a_h[6 * 7 + 5] = 11.0/84.0;    //65

	b_h[0] = 35.0/384.0;
	b_h[1] = 0.0;
	b_h[2] = 500.0/1113.0;
	b_h[3] = 125.0/192.0;
	b_h[4] = -2187.0/6784.0;
	b_h[5] = 11.0/84.0;
	b_h[6] = 0.0;

	bb_h[0] = 5179.0/57600.0;
	bb_h[1] = 0.0;
	bb_h[2] = 7571.0/16695.0;
	bb_h[3] = 393.0/640.0;
	bb_h[4] = -92097.0/339200.0;
	bb_h[5] = 187.0/2100.0;
	bb_h[6] = 1.0/40.0;

	c_h[0] = 0.0;
	c_h[1] = 1.0/5.0;
	c_h[2] = 3.0 / 10.0;
	c_h[3] = 4.0 / 5.0;
	c_h[4] = 8.0/9.0;
	c_h[5] = 1.0;
	c_h[6] = 1.0;

}


void setRKF78(double *a_h, double *b_h, double *bb_h, double *c_h){


	a_h[1 * 13 + 0] = 2.0/27.0;

	a_h[2 * 13 + 0] = 1.0/36.0;
	a_h[2 * 13 + 1] = 1.0/12.0;

	a_h[3 * 13 + 0] = 1.0/24.0;
	a_h[3 * 13 + 1] = 0.0;
	a_h[3 * 13 + 2] = 1.0/8.0;

	a_h[4 * 13 + 0] = 5.0/12.0;
	a_h[4 * 13 + 1] = 0.0;
	a_h[4 * 13 + 2] = -25.0/16.0;
	a_h[4 * 13 + 3] = 25.0/16.0;

	a_h[5 * 13 + 0] = 1.0/20.0;
	a_h[5 * 13 + 1] = 0.0;
	a_h[5 * 13 + 2] = 0.0;
	a_h[5 * 13 + 3] = 1.0/4.0;
	a_h[5 * 13 + 4] = 1.0/5.0;

	a_h[6 * 13 + 0] = -25.0/ 108.0;
	a_h[6 * 13 + 1] = 0.0;
	a_h[6 * 13 + 2] = 0.0;
	a_h[6 * 13 + 3] = 125.0/ 108.0;
	a_h[6 * 13 + 4] = -65.0/ 27.0;
	a_h[6 * 13 + 5] = 125.0/ 54.0;

	a_h[7 * 13 + 0] = 31.0/300.0;
	a_h[7 * 13 + 1] = 0.0;
	a_h[7 * 13 + 2] = 0.0;
	a_h[7 * 13 + 3] = 0.0;
	a_h[7 * 13 + 4] = 61.0/225.0;
	a_h[7 * 13 + 5] = -2.0/9.0;
	a_h[7 * 13 + 6] = 13.0/900.0;

	a_h[8 * 13 + 0] = 2.0;
	a_h[8 * 13 + 1] = 0.0;
	a_h[8 * 13 + 2] = 0.0;
	a_h[8 * 13 + 3] = -53.0/6.0;
	a_h[8 * 13 + 4] = 704.0/45.0;
	a_h[8 * 13 + 5] = -107.0/9.0;
	a_h[8 * 13 + 6] = 67.0/90.0;
	a_h[8 * 13 + 7] = 3.0;

	a_h[9 * 13 + 0] = -91.0/108.0;
	a_h[9 * 13 + 1] = 0.0;
	a_h[9 * 13 + 2] = 0.0;
	a_h[9 * 13 + 3] = 23.0/108.0;
	a_h[9 * 13 + 4] = -976.0/135.0;
	a_h[9 * 13 + 5] = 311.0/54.0;
	a_h[9 * 13 + 6] = -19.0/60.0;
	a_h[9 * 13 + 7] = 17.0/6.0;
	a_h[9 * 13 + 8] = -1.0/12.0;

	a_h[10 * 13 + 0] = 2383.0/4100.0;
	a_h[10 * 13 + 1] = 0.0;
	a_h[10 * 13 + 2] = 0.0;
	a_h[10 * 13 + 3] = -341.0/164.0;
	a_h[10 * 13 + 4] = 4496.0/1025.0;
	a_h[10 * 13 + 5] = -301.0/82.0;
	a_h[10 * 13 + 6] = 2133.0/4100.0;
	a_h[10 * 13 + 7] = 45.0/82.0;
	a_h[10 * 13 + 8] = 45.0/164.0;
	a_h[10 * 13 + 9] = 18.0/41.0;

	a_h[11 * 13 + 0] = 3.0/205.0;
	a_h[11 * 13 + 1] = 0.0;
	a_h[11 * 13 + 2] = 0.0;
	a_h[11 * 13 + 3] = 0.0;
	a_h[11 * 13 + 4] = 0.0;
	a_h[11 * 13 + 5] = -6.0/41.0;
	a_h[11 * 13 + 6] = -3.0/205.0;
	a_h[11 * 13 + 7] = -3.0/41.0;
	a_h[11 * 13 + 8] = 3.0/41.0;
	a_h[11 * 13 + 9] = 6.0/41.0;
	a_h[11 * 13 + 10] = 0.0;

	a_h[12 * 13 + 0] = -1777.0/4100.0;
	a_h[12 * 13 + 1] = 0.0;
	a_h[12 * 13 + 2] = 0.0;
	a_h[12 * 13 + 3] = -341.0/164.0;
	a_h[12 * 13 + 4] = 4496.0/1025.0;
	a_h[12 * 13 + 5] = -289.0/82.0;
	a_h[12 * 13 + 6] = 2193.0/4100.0;
	a_h[12 * 13 + 7] = 51.0/82.0;
	a_h[12 * 13 + 8] = 33.0/164.0;
	a_h[12 * 13 + 9] = 12.0/41.0;		//19 / 41 or 12/41 ?
	a_h[12 * 13 + 10] = 0.0;
	a_h[12 * 13 + 11] = 1.0;		//0 or 1 ?



	b_h[0] = 41.0/840.0;
	b_h[1] = 0.0;
	b_h[2] = 0.0;
	b_h[3] = 0.0;
	b_h[4] = 0.0;
	b_h[5] = 34.0/105.0;
	b_h[6] = 9.0/35.0;
	b_h[7] = 9.0/35.0;
	b_h[8] = 9.0/280.0;
	b_h[9] = 9.0/280.0;
	b_h[10] = 41.0/840.0;
	b_h[11] = 0.0;
	b_h[12] = 0.0;

	bb_h[0] = 0.0;
	bb_h[1] = 0.0;
	bb_h[2] = 0.0;
	bb_h[3] = 0.0;
	bb_h[4] = 0.0;
	bb_h[5] = 34.0/105.0;
	bb_h[6] = 9.0/35.0;
	bb_h[7] = 9.0/35.0;
	bb_h[8] = 9.0/280.0;
	bb_h[9] = 9.0/280.0;
	bb_h[10] = 0.0;
	bb_h[11] = 41.0/840.0;
	bb_h[12] = 41.0/840.0;


	c_h[0] = 0.0;
	c_h[1] = 2.0/27.0;
	c_h[2] = 1.0/9.0;
	c_h[3] = 1.0/6.0;
	c_h[4] = 5.0/12.0;
	c_h[5] = 1.0/2.0;
	c_h[6] = 5.0/6.0; 
	c_h[7] = 1.0/6.0;
	c_h[8] = 2.0/3.0;
	c_h[9] = 1.0/3.0;
	c_h[10] = 1.0;
	c_h[11] = 0.0;
	c_h[12] = 1.0;
/*
	
	for(int i = 0; i < 13; ++i){
		for(int j = 0; j < i + 1; ++j){

			printf("%d %d %.20g\n", i, j, a_h[i * 13 + j]);
		}
	}
	for(int i = 0; i < 13; ++i){
		printf("%d %.20g\n", i, b_h[i]);
	}
	for(int i = 0; i < 13; ++i){
		printf("%d %.20g\n", i, bb_h[i]);
	}
	for(int i = 0; i < 13; ++i){
		printf("%d %.20g\n", i, c_h[i]);
	}
*/	
}

//compute adaptive time step for only 1 particle
void computeError1(double2 *snew_h, double *kx_h, double *ky_h, double *kz_h, double *kvx_h, double *kvy_h, double *kvz_h, double *b_h, double *bb_h, int RKFn, int i, int N, double &snew, double dt, double ee){



	double sc = def_sc * RKFn;

	//error estimation
	double errorkx = 0.0;
	double errorky = 0.0;
	double errorkz = 0.0;
	double errorkvx = 0.0;
	double errorkvy = 0.0;
	double errorkvz = 0.0;

	for(int S = 0; S < RKFn; ++S){
		errorkx += (b_h[S] - bb_h[S]) * kx_h[S * N + i];
		errorky += (b_h[S] - bb_h[S]) * ky_h[S * N + i];
		errorkz += (b_h[S] - bb_h[S]) * kz_h[S * N + i];

		errorkvx += (b_h[S] - bb_h[S]) * kvx_h[S * N + i];
		errorkvy += (b_h[S] - bb_h[S]) * kvy_h[S * N + i];
		errorkvz += (b_h[S] - bb_h[S]) * kvz_h[S * N + i];
//printf("error %d %d %g %g\n", i, S, errorkx, kx_h[S * N + i]);
	}

	errorkx = sqrt(errorkx * errorkx) / sc;
	errorky = sqrt(errorky * errorky) / sc;
	errorkz = sqrt(errorkz * errorkz) / sc;
	errorkvx = sqrt(errorkvx * errorkvx) / sc;
	errorkvy = sqrt(errorkvy * errorkvy) / sc;
	errorkvz = sqrt(errorkvz * errorkvz) / sc;

	double errmax = errorkx;
	errmax = fmax(errmax, errorky);
	errmax = fmax(errmax, errorkz);
	errmax = fmax(errmax, errorkvx);
	errmax = fmax(errmax, errorkvy);
	errmax = fmax(errmax, errorkvz);

	double s = pow( 1.0  / errmax, ee);
//printf("error %d %g %g %g\n", i, errmax, s, errorkx);
	s = fmax(def_facmin, def_fac * s);
	s = fmin(def_facmax, s);

	snew_h[i].x = s;
	snew_h[i].y = fmin(snew_h[i].y, s);
	snew = s;

}


//compute adaptive time step for all particles
void computeError(double2 *snew_h, double *dx_h, double *dy_h, double *dz_h, double *dvx_h, double *dvy_h, double *dvz_h, double *kx_h, double *ky_h, double *kz_h, double *kvx_h, double *kvy_h, double *kvz_h, double *b_h, double *bb_h, int RKFn, int Nperturbers, int N, double &snew, double dt, double dti, double dtiMin, double ee){


	double sc = def_sc * RKFn;
	double s;	//large number

	//error estimation
	for(int i = Nperturbers; i < N; ++i){
		//update
		double dx = 0.0;
		double dy = 0.0;
		double dz = 0.0;
		double dvx = 0.0;
		double dvy = 0.0;
		double dvz = 0.0;

		for(int S = 0; S < RKFn; ++S){
			double dtb = dt * b_h[S];
			dx += dtb * kx_h[i + S * N];
			dy += dtb * ky_h[i + S * N];
			dz += dtb * kz_h[i + S * N];

			dvx += dtb * kvx_h[i + S * N];
			dvy += dtb * kvy_h[i + S * N];
			dvz += dtb * kvz_h[i + S * N];
		}

		dx_h[i] = dx;
		dy_h[i] = dy;
		dz_h[i] = dz;
		dvx_h[i] = dvx;
		dvy_h[i] = dvy;
		dvz_h[i] = dvz;



		double errorkx = 0.0;
		double errorky = 0.0;
		double errorkz = 0.0;
		double errorkvx = 0.0;
		double errorkvy = 0.0;
		double errorkvz = 0.0;

		for(int S = 0; S < RKFn; ++S){
			errorkx += (b_h[S] - bb_h[S]) * kx_h[S * N + i];
			errorky += (b_h[S] - bb_h[S]) * ky_h[S * N + i];
			errorkz += (b_h[S] - bb_h[S]) * kz_h[S * N + i];

			errorkvx += (b_h[S] - bb_h[S]) * kvx_h[S * N + i];
			errorkvy += (b_h[S] - bb_h[S]) * kvy_h[S * N + i];
			errorkvz += (b_h[S] - bb_h[S]) * kvz_h[S * N + i];
//printf("error %d %d %g %g\n", i, S, errorkx, kx_h[S * N + i]);
		}

		errorkx = sqrt(errorkx * errorkx) / sc;
		errorky = sqrt(errorky * errorky) / sc;
		errorkz = sqrt(errorkz * errorkz) / sc;
		errorkvx = sqrt(errorkvx * errorkvx) / sc;
		errorkvy = sqrt(errorkvy * errorkvy) / sc;
		errorkvz = sqrt(errorkvz * errorkvz) / sc;

		double errmax = errorkx;
		errmax = fmax(errmax, errorky);
		errmax = fmax(errmax, errorkz);
		errmax = fmax(errmax, errorkvx);
		errmax = fmax(errmax, errorkvy);
		errmax = fmax(errmax, errorkvz);
//printf("error, %d %d %g %g %g %d\n", id + i, blockIdx.x, errork, s, errorkx, i);

		s = pow( 1.0  / errmax, ee);
//printf("error %d %g %g %g %g %g %g\n", i, errmax, s, errorkx, s * dt, dt, sc);
		s = fmax(def_facmin, def_fac * s);
		s = fmin(def_facmax, s);


		// ****************************
		if(snew_h[i].y >= 1.0){
			snew_h[i].x = s;
		}
		else{
			snew_h[i].x = 1.5;
		}
		if(s * dti < dtiMin){
			snew_h[i].y = fmin(snew_h[i].y, s);
		}
//printf("snew %d %.20g %.20g %.20g %.20g\n", i, snew_h[i].y, snew_h[i].x, s, errmax);
		// ****************************
	
	}

// remove this
	snew = 1.0e6;	//large number
	for(int i = Nperturbers; i < N; ++i){
		snew = fmin(snew, snew_h[i].x);
	}

}


__global__ void computeError_d1_kernel(double2 *snew_d, int Nperturbers, int N){

	int id = blockIdx.x * blockDim.x + threadIdx.x;

	double s = 1.0e6;	//large number

	extern __shared__ double se_s[];
	double *s_s = se_s;

	int lane = threadIdx.x % warpSize;
	int warp = threadIdx.x / warpSize;

	if(warp == 0){
		s_s[threadIdx.x] = 1.0e6;
	}
	__syncthreads();

	if(id >= Nperturbers && id < N){	
		s = snew_d[id].x;
	}
//printf("snew %d %g\n", id, s);

	__syncthreads();

	double sr = s;

	for(int i = 1; i < warpSize; i*=2){
#if def_OldShuffle == 0
		sr = __shfl_xor_sync(0xffffffff, s, i, warpSize);
#else
                sr = __shfld_xor(s, i);
#endif
		s = fmin(s, sr);
//printf("s reduce  %d %d %d %.20g\n", blockIdx.x, i, threadIdx.x, s);
        }
//printf("s reduce  %d %d %d %d %.20g\n", blockIdx.x, 0, id, threadIdx.x, s);

	__syncthreads();
	if(blockDim.x > warpSize){
		//reduce across warps
		if(lane == 0){
			s_s[warp] = s;
		}
		__syncthreads();
		//reduce previous warp results in the first warp
		if(warp == 0){
			s = s_s[threadIdx.x];
//printf("r reduce 1  %d %d %d %.20g %d %d\n", blockIdx.x, 0, threadIdx.x, s, int(blockDim.x), warpSize);
			for(int i = 1; i < warpSize; i*=2){
#if def_OldShuffle == 0
				sr = __shfl_xor_sync(0xffffffff, s, i, warpSize);
#else
				sr = __shfld_xor(s, i);
#endif
				s = fmin(s, sr);
//printf("r reduce 2  %d %d %d %.20g\n", blockIdx.x, i, threadIdx.x, s);
			}
		}
	}
	__syncthreads();
	if(threadIdx.x == 0){

		snew_d[blockIdx.x].x = s;
//printf("r reduce 3  %d %d %.20g\n", blockIdx.x, threadIdx.x, s);
	}

}

__global__ void computeError_d2_kernel(double2 *snew_d, const int N){

	int idy = threadIdx.x;
	int idx = blockIdx.x;

	double s = 1.0e6;

	extern __shared__ double se2_s[];
	double *s_s = se2_s;

	int lane = threadIdx.x % warpSize;
	int warp = threadIdx.x / warpSize;

	if(warp == 0){
		s_s[threadIdx.x] = 1.0e6;
	}


	if(idy < N){
		s = snew_d[idy].x;
	}
//printf("s reduce d2  %d %.20g\n", idy, s);

	__syncthreads();

	double sr = s;

	for(int i = 1; i < warpSize; i*=2){
#if def_OldShuffle == 0
		sr = __shfl_xor_sync(0xffffffff, s, i, warpSize);
#else
		sr = __shfld_xor(s, i);
#endif
		s = fmin(s, sr);
//if(idx == 0) printf("HC2bx %d %d %.20g\n", i, idy, a);
	}
	__syncthreads();

	if(blockDim.x > warpSize){
		//reduce across warps
		if(lane == 0){
			s_s[warp] = s;
		}
		__syncthreads();
		//reduce previous warp results in the first warp
		if(warp == 0){
			s = s_s[threadIdx.x];
			//if(idx == 0) printf("HC2cx %d %d %.20g %d %d\n", 0, idy, a, int(blockDim.x), warpSize);
			for(int i = 1; i < warpSize; i*=2){
#if def_OldShuffle == 0
				sr = __shfl_xor_sync(0xffffffff, s, i, warpSize);
#else
				sr = __shfld_xor(s, i);
#endif
				s = fmin(s, sr);
//printf("s reduce d2b %d %d %.20g\n", i, idy, s);
			}
		}
	}
	__syncthreads();
	if(threadIdx.x == 0){
		if(idx == 0){
//printf("s reduce d2c %d %.20g\n", idy, s);
			snew_d[0].x = s;
		}
	}
}



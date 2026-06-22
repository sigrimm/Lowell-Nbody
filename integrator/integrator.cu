#include <stdio.h>

//Define Constant Memory
__constant__ double RKFa_c[20 * 20];	//20 is considered here to be large enough (>RKFn)
__constant__ double RKFb_c[20];
__constant__ double RKFbb_c[20];
__constant__ double RKFc_c[80];		//maximum of number of RK stages os Bulish-Stoer stages


__constant__ double BSddt_c[8];
__constant__ double BSt0_c[8 * 8];

__constant__ double atol_c;
__constant__ double rtol_c;
__constant__ double RKF_ee_c;
__constant__ double RKF_fac_c;
__constant__ double RKF_facmin_c;
__constant__ double RKF_facmax_c;


__constant__ double REAU_c;
__constant__ double J2E_c;
__constant__ double c2_c;

__constant__ int useGR_c;
__constant__ int useJ2_c;
__constant__ int useNonGrav_c;
__constant__ int cometFlag_c;
__constant__ int Rbuffersize_c;

__constant__ double nonGrav_alpha_c;
__constant__ double nonGrav_nk_c;
__constant__ double nonGrav_nm_c;
__constant__ double nonGrav_nn_c;
__constant__ double nonGrav_r0_c;
__constant__ double nonGrav_tau_c;

#include "asteroid.h"
#include "ChebyshevGPU.h"
#include "forceGPU.h"


__host__ int asteroid::copyConst(){
{
cudaDeviceSynchronize();
cudaError_t error = cudaGetLastError();
printf("copy const error 1= %d = %s\n",error, cudaGetErrorString(error));
if(error != 0.0){
	return 0;
}
}

	if(RKFn > 0){
		cudaMemcpyToSymbol(RKFa_c, RKFa_h, RKFn * RKFn * sizeof(double), 0, cudaMemcpyHostToDevice);
		cudaMemcpyToSymbol(RKFb_c, RKFb_h, RKFn * sizeof(double), 0, cudaMemcpyHostToDevice);
		cudaMemcpyToSymbol(RKFbb_c, RKFbb_h, RKFn * sizeof(double), 0, cudaMemcpyHostToDevice);
		cudaMemcpyToSymbol(RKFc_c, RKFc_h, RKFn * sizeof(double), 0, cudaMemcpyHostToDevice);

		cudaMemcpyToSymbol(RKF_ee_c, &RKF_ee, sizeof(double), 0, cudaMemcpyHostToDevice);
		cudaMemcpyToSymbol(RKF_fac_c, &RKF_fac, sizeof(double), 0, cudaMemcpyHostToDevice);
		cudaMemcpyToSymbol(RKF_facmin_c, &RKF_facmin, sizeof(double), 0, cudaMemcpyHostToDevice);
		cudaMemcpyToSymbol(RKF_facmax_c, &RKF_facmax, sizeof(double), 0, cudaMemcpyHostToDevice);
	}
	
	if(BSn > 0){
		cudaMemcpyToSymbol(BSddt_c, BSddt_h, 8 * sizeof(double), 0, cudaMemcpyHostToDevice);
		cudaMemcpyToSymbol(BSt0_c, BSt0_h, 8 * 8 * sizeof(double), 0, cudaMemcpyHostToDevice);
		cudaMemcpyToSymbol(RKFc_c, BSc_h, BSn * sizeof(double), 0, cudaMemcpyHostToDevice);
	}

	cudaMemcpyToSymbol(atol_c, &atol, sizeof(double), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(rtol_c, &rtol, sizeof(double), 0, cudaMemcpyHostToDevice);

	cudaMemcpyToSymbol(REAU_c, &REAU, sizeof(double), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(J2E_c, &J2E, sizeof(double), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c2_c, &c2, sizeof(double), 0, cudaMemcpyHostToDevice);

	cudaMemcpyToSymbol(useGR_c, &useGR, sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(useJ2_c, &useJ2, sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(useNonGrav_c, &useNonGrav, sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(cometFlag_c, &cometFlag, sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(Rbuffersize_c, &Rbuffersize, sizeof(int), 0, cudaMemcpyHostToDevice);

	cudaMemcpyToSymbol(nonGrav_alpha_c, &nonGrav_alpha, sizeof(double), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(nonGrav_nk_c, &nonGrav_nk, sizeof(double), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(nonGrav_nm_c, &nonGrav_nm, sizeof(double), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(nonGrav_nn_c, &nonGrav_nn, sizeof(double), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(nonGrav_r0_c, &nonGrav_r0, sizeof(double), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(nonGrav_tau_c, &nonGrav_tau, sizeof(double), 0, cudaMemcpyHostToDevice);

	cudaDeviceSynchronize();
	cudaError_t error = cudaGetLastError();
	printf("copy const error = %d = %s\n",error, cudaGetErrorString(error));
	if(error != 0.0){
		return 0;
	}

	return 1;

}




__global__ void HelioToBary_kernel(double *xTable_d, double *yTable_d, double *zTable_d, double *vxTable_d, double *vyTable_d, double *vzTable_d, double *xx_d, double *yy_d, double *zz_d, double *vxx_d, double *vyy_d, double *vzz_d, const int N, const int nStage){

	int id = blockIdx.x * blockDim.x + threadIdx.x;

	if(id < N){
		int ii = 10 * nStage;

		xx_d[id] += xTable_d[ii];
		yy_d[id] += yTable_d[ii];
		zz_d[id] += zTable_d[ii];

		vxx_d[id] += vxTable_d[ii];
		vyy_d[id] += vyTable_d[ii];
		vzz_d[id] += vzTable_d[ii];

//printf("xyz bary GPU %.40g %.40g %.40g %.40g %.40g %.40g\n", x_d[id], y_d[id], z_d[id], vx_d[id], vy_d[id], vz_d[id]);
	}
}


//Leapfrog step with fixed time step
//Every body runs on a thread
__global__ void leapfrog_stepA_kernel(double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, const int N, const double dt){

	int id = blockIdx.x * blockDim.x + threadIdx.x;

	//Drift
	if(id < N){
		x_d[id] += 0.5 * dt * vx_d[id];
		y_d[id] += 0.5 * dt * vy_d[id];
		z_d[id] += 0.5 * dt * vz_d[id];
	}

}

//Leapfrog step with fixed time step
//Kernel uses at least Nperturbers threads. The perturbers are loaded into shared memory
//Every body runs on a thread
__global__ void leapfrog_stepB_kernel(double *xTable_d, double *yTable_d, double *zTable_d, double *vxTable_d, double *vyTable_d, double *vzTable_d, double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *A1_d, double *A2_d, double *A3_d, double *Tsave_d, double *Rsave_d, double time, long long timeStep, double *GM_d, const int Nperturbers, const int N, const double dt){

	int itx = threadIdx.x;
	int id = blockIdx.x * blockDim.x + threadIdx.x;

	double ax = 0.0;
	double ay = 0.0;
	double az = 0.0;

	//compute force

	double xi;
	double yi;
	double zi;

	double vxi;
	double vyi;
	double vzi;

	__shared__ double xTable_s[def_NP];
	__shared__ double yTable_s[def_NP];
	__shared__ double zTable_s[def_NP];

	__shared__ double vxTable_s[def_NP];
	__shared__ double vyTable_s[def_NP];
	__shared__ double vzTable_s[def_NP];

	__shared__ double GM_s[def_NP];

	if(itx < Nperturbers){
		xTable_s[itx] = xTable_d[itx];
		yTable_s[itx] = yTable_d[itx];
		zTable_s[itx] = zTable_d[itx];

		vxTable_s[itx] = vxTable_d[itx];
		vyTable_s[itx] = vyTable_d[itx];
		vzTable_s[itx] = vzTable_d[itx];

		GM_s[itx] = GM_d[itx];
	}

	__syncthreads();


	if(id < N){

		xi = x_d[id];
		yi = y_d[id];
		zi = z_d[id];

		vxi = vx_d[id];
		vyi = vy_d[id];
		vzi = vz_d[id];

		//heliocentric coordinates
		double xih = xi - xTable_s[10];
		double yih = yi - yTable_s[10];
		double zih = zi - zTable_s[10];

		double vxih = vxi - vxTable_s[10];
		double vyih = vyi - vyTable_s[10];
		double vzih = vzi - vzTable_s[10];

		//r is used in multiple forces, so reuse it
		double rsq = __dmul_rn(xih, xih) + __dmul_rn(yih, yih) + __dmul_rn(zih, zih);
		double r = sqrt(rsq);

		if(cometFlag_c > 0){
			if(id == 0){
				Tsave_d[timeStep] = time;
			}
			Rsave_d[id * Rbuffersize_c + timeStep] = r;
		}

		//Earth centric coordinates
		double xiE = xi - xTable_s[2];
		double yiE = yi - yTable_s[2];
		double ziE = zi - zTable_s[2];

		if(useNonGrav_c == 1){
			NonGrav(xih, yih, zih, vxih, vyih, vzih, A1_d[id], A2_d[id], A3_d[id], r, ax, ay, az);
		}
		if(useGR_c == 1){
			GR(xih, yih, zih, vxih, vyih, vzih, r, ax, ay, az, GM_s[10]);
		}
		if(useJ2_c == 1){
			J2(xiE, yiE, ziE, ax, ay, az, GM_s[2]);
		}
		for(int p = 0; p < Nperturbers; ++p){
			Gravity(xi, yi, zi, xTable_s, yTable_s, zTable_s, ax, ay, az, GM_s, p);
		}
	}
	// ----------------------------------------------------------------------------

	__syncthreads();

	if(id < N){
		//Kick
		vx_d[id] += dt * ax;
		vy_d[id] += dt * ay;
		vz_d[id] += dt * az;

		//Drift
		x_d[id] += 0.5 * dt * vx_d[id];
		y_d[id] += 0.5 * dt * vy_d[id];
		z_d[id] += 0.5 * dt * vz_d[id];
	}

}
//Leapfrog step with fixed time step
// Every body runs on a sepparate thread block. Gravity calculation is spread along threads in the block and reuced in registers
// Kernel uses at least Nperturbers threads. The perturbers are loaded into shared memory
__global__ void leapfrog_stepB2_kernel(double *xTable_d, double *yTable_d, double *zTable_d, double *vxTable_d, double *vyTable_d, double *vzTable_d, double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *A1_d, double *A2_d, double *A3_d, double *Tsave_d, double *Rsave_d, double time, long long timeStep, double *GM_d, const int Nperturbers, const int N, const double dt){

	int itx = threadIdx.x;
	int id = blockIdx.x;    //particle index

	double ax = 0.0;
	double ay = 0.0;
	double az = 0.0;

	//compute force

	double xi;
	double yi;
	double zi;

	double vxi;
	double vyi;
	double vzi;

	__shared__ double xTable_s[def_NP];
	__shared__ double yTable_s[def_NP];
	__shared__ double zTable_s[def_NP];

	__shared__ double vxTable_s[def_NP];
	__shared__ double vyTable_s[def_NP];
	__shared__ double vzTable_s[def_NP];

	__shared__ double GM_s[def_NP];

	if(itx < Nperturbers){
		xTable_s[itx] = xTable_d[itx];
		yTable_s[itx] = yTable_d[itx];
		zTable_s[itx] = zTable_d[itx];

		vxTable_s[itx] = vxTable_d[itx];
		vyTable_s[itx] = vyTable_d[itx];
		vzTable_s[itx] = vzTable_d[itx];

		GM_s[itx] = GM_d[itx];
	}

	__syncthreads();


	if(id < N){

		xi = x_d[id];
		yi = y_d[id];
		zi = z_d[id];

		vxi = vx_d[id];
		vyi = vy_d[id];
		vzi = vz_d[id];

		if(itx == 0){

			//heliocentric coordinates
			double xih = xi - xTable_s[10];
			double yih = yi - yTable_s[10];
			double zih = zi - zTable_s[10];

			double vxih = vxi - vxTable_s[10];
			double vyih = vyi - vyTable_s[10];
			double vzih = vzi - vzTable_s[10];

			//r is used in multiple forces, so reuse it
			double rsq = __dmul_rn(xih, xih) + __dmul_rn(yih, yih) + __dmul_rn(zih, zih);
			double r = sqrt(rsq);

			if(cometFlag_c > 0){
				if(id == 0){
					Tsave_d[timeStep] = time;
				}
				Rsave_d[id * Rbuffersize_c + timeStep] = r;
			}

			//Earth centric coordinates
			double xiE = xi - xTable_s[2];
			double yiE = yi - yTable_s[2];
			double ziE = zi - zTable_s[2];

			if(useNonGrav_c == 1){
				NonGrav(xih, yih, zih, vxih, vyih, vzih, A1_d[id], A2_d[id], A3_d[id], r, ax, ay, az);
			}
			if(useGR_c == 1){
				GR(xih, yih, zih, vxih, vyih, vzih, r, ax, ay, az, GM_s[10]);
			}
			if(useJ2_c == 1){
				J2(xiE, yiE, ziE, ax, ay, az, GM_s[2]);
			}
		}
		__syncthreads();
		if(itx < Nperturbers){
			Gravity(xi, yi, zi, xTable_s, yTable_s, zTable_s, ax, ay, az, GM_s, itx);
		}
		__syncthreads();
		//reduce a
		// ----------------------------------------------------------------------------
		for(int i = 1; i < warpSize; i*=2){
			ax += __shfl_xor_sync(0xffffffff, ax, i, warpSize);
			ay += __shfl_xor_sync(0xffffffff, ay, i, warpSize);
			az += __shfl_xor_sync(0xffffffff, az, i, warpSize);
		}
		__syncthreads();

	}
	// ----------------------------------------------------------------------------

	if(id < N  && itx == 0){
		//Kick
		vx_d[id] += dt * ax;
		vy_d[id] += dt * ay;
		vz_d[id] += dt * az;

		//Drift
		x_d[id] += 0.5 * dt * vx_d[id];
		y_d[id] += 0.5 * dt * vy_d[id];
		z_d[id] += 0.5 * dt * vz_d[id];
	}
}

//Implicit Midpoint Method
//Kernel uses at least Nperturbers threads. The perturbers are loaded into shared memory
//Every body runs on a thread
__global__ void IMM_step_kernel(double *xTable_d, double *yTable_d, double *zTable_d, double *vxTable_d, double *vyTable_d, double *vzTable_d, double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *A1_d, double *A2_d, double *A3_d, double *Tsave_d, double *Rsave_d, double *snew_d, double time, long long timeStep, double *GM_d, const int Nperturbers, const int N, const double dt){

	int itx = threadIdx.x;
	int id = blockIdx.x * blockDim.x + threadIdx.x;

	double ax = 0.0;
	double ay = 0.0;
	double az = 0.0;

	//compute force

	double xi;
	double yi;
	double zi;

	double vxi;
	double vyi;
	double vzi;

	double xti;
	double yti;
	double zti;

	double vxti;
	double vyti;
	double vzti;

	double xpi;
	double ypi;
	double zpi;

	double vxpi;
	double vypi;
	double vzpi;

	double A1;
	double A2;
	double A3;

	__shared__ double xTable_s[def_NP];
	__shared__ double yTable_s[def_NP];
	__shared__ double zTable_s[def_NP];

	__shared__ double vxTable_s[def_NP];
	__shared__ double vyTable_s[def_NP];
	__shared__ double vzTable_s[def_NP];

	__shared__ double GM_s[def_NP];

	__shared__ int stop_s[1];


	if(itx < Nperturbers){
		xTable_s[itx] = xTable_d[itx];
		yTable_s[itx] = yTable_d[itx];
		zTable_s[itx] = zTable_d[itx];

		vxTable_s[itx] = vxTable_d[itx];
		vyTable_s[itx] = vyTable_d[itx];
		vzTable_s[itx] = vzTable_d[itx];

		GM_s[itx] = GM_d[itx];
	}


	if(id < N){
			xi = x_d[id];
			yi = y_d[id];
			zi = z_d[id];

			vxi = vx_d[id];
			vyi = vy_d[id];
			vzi = vz_d[id];

			xti = xi;
			yti = yi;
			zti = zi;

			vxti = vxi;
			vyti = vyi;
			vzti = vzi;

			A1 = A1_d[id];
			A2 = A2_d[id];
			A3 = A3_d[id];

	}

	__syncthreads();

	for(int k = 0; k < 30; ++k){


		// ----------------------------------------------------------------------------
		//compute forces
		if(id < N){

			ax = 0.0;
			ay = 0.0;
			az = 0.0;


			//heliocentric coordinates
			double xih = xti - xTable_s[10];
			double yih = yti - yTable_s[10];
			double zih = zti - zTable_s[10];

			double vxih = vxti - vxTable_s[10];
			double vyih = vyti - vyTable_s[10];
			double vzih = vzti - vzTable_s[10];

			//r is used in multiple forces, so reuse it
			double rsq = __dmul_rn(xih, xih) + __dmul_rn(yih, yih) + __dmul_rn(zih, zih);
			double r = sqrt(rsq);

			if(cometFlag_c > 0 && k == 0){
				if(id == 0){
					Tsave_d[timeStep] = time;
				}
				Rsave_d[id * Rbuffersize_c + timeStep] = r;
			}

			//Earth centric coordinates
			double xiE = xti - xTable_s[2];
			double yiE = yti - yTable_s[2];
			double ziE = zti - zTable_s[2];

			if(useNonGrav_c == 1){
				NonGrav(xih, yih, zih, vxih, vyih, vzih, A1, A2, A3, r, ax, ay, az);
			}
			if(useGR_c == 1){
				GR(xih, yih, zih, vxih, vyih, vzih, r, ax, ay, az, GM_s[10]);
			}
			if(useJ2_c == 1){
				J2(xiE, yiE, ziE, ax, ay, az, GM_s[2]);
			}
			for(int p = 0; p < Nperturbers; ++p){
				Gravity(xti, yti, zti, xTable_s, yTable_s, zTable_s, ax, ay, az, GM_s, p);
			}
		}
		// ----------------------------------------------------------------------------
		stop_s[0] = 1;

		__syncthreads();

		if(id < N){
//printf("k %d %.20g\n", k, ax);
			//store old values to check for convergence
			xpi = xti;
			ypi = yti;
			zpi = zti;

			vxpi = vxti;
			vypi = vyti;
			vzpi = vzti;

			xti = xi + 0.5 * dt * vxti;
			yti = yi + 0.5 * dt * vyti;
			zti = zi + 0.5 * dt * vzti;

			vxti = vxi + 0.5 * dt * ax;
			vyti = vyi + 0.5 * dt * ay;
			vzti = vzi + 0.5 * dt * az;

			if(fabs(xpi - xti) >= 1.0e-20) stop_s[0] = 0;
			if(fabs(ypi - yti) >= 1.0e-20) stop_s[0] = 0;
			if(fabs(zpi - zti) >= 1.0e-20) stop_s[0] = 0;

			if(fabs(vxpi - vxti) >= 1.0e-20) stop_s[0] = 0;
			if(fabs(vypi - vyti) >= 1.0e-20) stop_s[0] = 0;
			if(fabs(vzpi - vzti) >= 1.0e-20) stop_s[0] = 0;
		}

		__syncthreads();
		if(stop_s[0] == 1){
//if(k > 1) printf("k %d\n", k);	
			break;
		}
		if(k >= 29){
			snew_d[0] = -1;
		}



	}//end of k loop


	//update
	if(id < N){
		x_d[id] += dt * vxti;
		y_d[id] += dt * vyti;
		z_d[id] += dt * vzti;

		vx_d[id] += dt * ax;
		vy_d[id] += dt * ay;
		vz_d[id] += dt * az;
//printf("ax %.20g\n", ax);
	}
}

//Implicit Midpoint Method
// Every body runs on a sepparate thread block. Gravity calculation is spread along threads in the block and reuced in registers
// Kernel uses at least Nperturbers threads. The perturbers are loaded into shared memory
__global__ void IMM_step2_kernel(double *xTable_d, double *yTable_d, double *zTable_d, double *vxTable_d, double *vyTable_d, double *vzTable_d, double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *A1_d, double *A2_d, double *A3_d, double *Tsave_d, double *Rsave_d, double *snew_d, double time, long long timeStep, double *GM_d, const int Nperturbers, const int N, const double dt){

	int itx = threadIdx.x;
	int id = blockIdx.x;    //particle index

	double ax = 0.0;
	double ay = 0.0;
	double az = 0.0;

	//compute force

	double xi;
	double yi;
	double zi;

	double vxi;
	double vyi;
	double vzi;

	double xti;
	double yti;
	double zti;

	double vxti;
	double vyti;
	double vzti;

	double xpi;
	double ypi;
	double zpi;

	double vxpi;
	double vypi;
	double vzpi;

	double A1;
	double A2;
	double A3;

	__shared__ double xTable_s[def_NP];
	__shared__ double yTable_s[def_NP];
	__shared__ double zTable_s[def_NP];

	__shared__ double vxTable_s[def_NP];
	__shared__ double vyTable_s[def_NP];
	__shared__ double vzTable_s[def_NP];

	__shared__ double GM_s[def_NP];

	__shared__ int stop_s[1];


	if(itx < Nperturbers){
		xTable_s[itx] = xTable_d[itx];
		yTable_s[itx] = yTable_d[itx];
		zTable_s[itx] = zTable_d[itx];

		vxTable_s[itx] = vxTable_d[itx];
		vyTable_s[itx] = vyTable_d[itx];
		vzTable_s[itx] = vzTable_d[itx];

		GM_s[itx] = GM_d[itx];
	}


	if(id < N){
			xi = x_d[id];
			yi = y_d[id];
			zi = z_d[id];

			vxi = vx_d[id];
			vyi = vy_d[id];
			vzi = vz_d[id];

			xti = xi;
			yti = yi;
			zti = zi;

			vxti = vxi;
			vyti = vyi;
			vzti = vzi;

			A1 = A1_d[id];
			A2 = A2_d[id];
			A3 = A3_d[id];

	}

	__syncthreads();

	for(int k = 0; k < 30; ++k){

		ax = 0.0;
		ay = 0.0;
		az = 0.0;

		// ----------------------------------------------------------------------------
		//compute forces
		if(id < N){

			if(itx == 0){

				//heliocentric coordinates
				double xih = xti - xTable_s[10];
				double yih = yti - yTable_s[10];
				double zih = zti - zTable_s[10];

				double vxih = vxti - vxTable_s[10];
				double vyih = vyti - vyTable_s[10];
				double vzih = vzti - vzTable_s[10];

				//r is used in multiple forces, so reuse it
				double rsq = __dmul_rn(xih, xih) + __dmul_rn(yih, yih) + __dmul_rn(zih, zih);
				double r = sqrt(rsq);

				if(cometFlag_c > 0 && k == 0){
					if(id == 0){
						Tsave_d[timeStep] = time;
					}
					Rsave_d[id * Rbuffersize_c + timeStep] = r;
				}

				//Earth centric coordinates
				double xiE = xti - xTable_s[2];
				double yiE = yti - yTable_s[2];
				double ziE = zti - zTable_s[2];

				if(useNonGrav_c == 1){
					NonGrav(xih, yih, zih, vxih, vyih, vzih, A1, A2, A3, r, ax, ay, az);
				}
				if(useGR_c == 1){
					GR(xih, yih, zih, vxih, vyih, vzih, r, ax, ay, az, GM_s[10]);
				}
				if(useJ2_c == 1){
					J2(xiE, yiE, ziE, ax, ay, az, GM_s[2]);
				}
			}
			__syncthreads();
			if(itx < Nperturbers){
				Gravity(xti, yti, zti, xTable_s, yTable_s, zTable_s, ax, ay, az, GM_s, itx);
			}
			__syncthreads();
			//reduce a
			// ----------------------------------------------------------------------------
			for(int i = 1; i < warpSize; i*=2){
				ax += __shfl_xor_sync(0xffffffff, ax, i, warpSize);
				ay += __shfl_xor_sync(0xffffffff, ay, i, warpSize);
				az += __shfl_xor_sync(0xffffffff, az, i, warpSize);
			}
			__syncthreads();

			// ----------------------------------------------------------------------------

			//store old values to check for convergence
			xpi = xti;
			ypi = yti;
			zpi = zti;

			vxpi = vxti;
			vypi = vyti;
			vzpi = vzti;

			xti = xi + 0.5 * dt * vxti;
			yti = yi + 0.5 * dt * vyti;
			zti = zi + 0.5 * dt * vzti;

			vxti = vxi + 0.5 * dt * ax;
			vyti = vyi + 0.5 * dt * ay;
			vzti = vzi + 0.5 * dt * az;


			if(itx == 0){
//printf("k %d %.20g\n", k, ax);

				stop_s[0] = 1;

				if(fabs(xpi - xti) >= 1.0e-20) stop_s[0] = 0;
				if(fabs(ypi - yti) >= 1.0e-20) stop_s[0] = 0;
				if(fabs(zpi - zti) >= 1.0e-20) stop_s[0] = 0;

				if(fabs(vxpi - vxti) >= 1.0e-20) stop_s[0] = 0;
				if(fabs(vypi - vyti) >= 1.0e-20) stop_s[0] = 0;
				if(fabs(vzpi - vzti) >= 1.0e-20) stop_s[0] = 0;

				if(k >= 29){
					snew_d[0] = -1;
				}
			}
			__syncthreads();
			if(stop_s[0] == 1){
//if(k > 1) printf("k %d\n", k);	
				break;
			}
		}

	}//end of k loop


	//update
	if(id < N && itx == 0){
		x_d[id] += dt * vxti;
		y_d[id] += dt * vyti;
		z_d[id] += dt * vzti;

		vx_d[id] += dt * ax;
		vy_d[id] += dt * ay;
		vz_d[id] += dt * az;
//printf("ax %.20g\n", ax);
	}
}


// Runge Kutta step with fixed time step
// Every body runs on a thread
// Kernel uses at least Nperturbers threads. The perturbers are loaded into shared memory
__global__ void RK_step_kernel(double *xTable_d, double *yTable_d, double *zTable_d, double *vxTable_d, double *vyTable_d, double *vzTable_d, double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *kx_d, double *ky_d, double *kz_d, double *kvx_d, double *kvy_d, double *kvz_d, double *GM_d, double *A1_d, double *A2_d, double *A3_d, double *Tsave_d, double *Rsave_d, double time, long long int timeStep, const double dt, const int RKFn, const int Nperturbers, const int N){


	int itx = threadIdx.x;
	int id = blockIdx.x * blockDim.x + threadIdx.x;

	//shared memory contains only the perturbers
	__shared__ double xTable_s[def_NP];
	__shared__ double yTable_s[def_NP];
	__shared__ double zTable_s[def_NP];

	__shared__ double vxTable_s[def_NP];
	__shared__ double vyTable_s[def_NP];
	__shared__ double vzTable_s[def_NP];

	__shared__ double GM_s[def_NP];


	double xti;
	double yti;
	double zti;

	double vxti;
	double vyti;
	double vzti;

	double x0;
	double y0;
	double z0;

	double vx0;
	double vy0;
	double vz0;

	double A1;
	double A2;
	double A3;

	if(itx < Nperturbers){
		GM_s[itx] = GM_d[itx];
	}
	if(id < N){
		x0 = x_d[id];
		y0 = y_d[id];
		z0 = z_d[id];

		vx0 = vx_d[id];
		vy0 = vy_d[id];
		vz0 = vz_d[id];

		A1 = A1_d[id];
		A2 = A2_d[id];
		A3 = A3_d[id];
	}
	__syncthreads();


	for(int S = 0; S < RKFn; ++S){

		// ----------------------------------------------------------------------------
		//Read the perturbers position
		if(itx < Nperturbers){
			int ii = itx * RKFn + S;

			xTable_s[itx] = xTable_d[ii];
			yTable_s[itx] = yTable_d[ii];
			zTable_s[itx] = zTable_d[ii];

			vxTable_s[itx] = vxTable_d[ii];
			vyTable_s[itx] = vyTable_d[ii];
			vzTable_s[itx] = vzTable_d[ii];
//printf("%d %d %.20g %.20g %.20g\n", S, itx, xTable_s[itx], yTable_s[itx], zTable_s[itx]);
		}
		__syncthreads();

//if(id < Nperturbers){
//printf("p %d %.20g %.20g %.20g %.20g %.20g %.20g %.20g\n", id, time + RKFc_c[S] * dt, xt_s[id], yt_s[id], zt_s[id], vxt_s[id], vyt_s[id], vzt_s[id]);
//}

		if(id < N){

			xti = x0;
			yti = y0;
			zti = z0;

			vxti = vx0;
			vyti = vy0;
			vzti = vz0;

			for(int s = 0; s < S; ++s){
				double dtaa = dt * RKFa_c[S * RKFn + s];
				int ii = id + s * N;
				xti  += dtaa * kx_d[ii];
				yti  += dtaa * ky_d[ii];
				zti  += dtaa * kz_d[ii];
				vxti += dtaa * kvx_d[ii];
				vyti += dtaa * kvy_d[ii];
				vzti += dtaa * kvz_d[ii];
			}

			int ik = id + S * N;

			kx_d[ik] = vxti;
			ky_d[ik] = vyti;
			kz_d[ik] = vzti;


			// ----------------------------------------------------------------------------
			//compute forces
			double ax = 0.0;
			double ay = 0.0;
			double az = 0.0;

			//heliocentric coordinates
			double xih = xti - xTable_s[10];
			double yih = yti - yTable_s[10];
			double zih = zti - zTable_s[10];

			double vxih = vxti - vxTable_s[10];
			double vyih = vyti - vyTable_s[10];
			double vzih = vzti - vzTable_s[10];

			//r is used in multiple forces, so reuse it
			double rsq = __dmul_rn(xih, xih) + __dmul_rn(yih, yih) + __dmul_rn(zih, zih);
			double r = sqrt(rsq);

			if(cometFlag_c > 0 && S == 0){
				if(id == 0){
					Tsave_d[timeStep] = time;
				}
				Rsave_d[id * Rbuffersize_c + timeStep] = r;
			}

			//Earth centric coordinates
			double xiE = xti - xTable_s[2];
			double yiE = yti - yTable_s[2];
			double ziE = zti - zTable_s[2];

			if(useNonGrav_c == 1){
				NonGrav(xih, yih, zih, vxih, vyih, vzih, A1, A2, A3, r, ax, ay, az);
			}
			if(useGR_c == 1){
				GR(xih, yih, zih, vxih, vyih, vzih, r, ax, ay, az, GM_s[10]);
			}
			if(useJ2_c == 1){
				J2(xiE, yiE, ziE, ax, ay, az, GM_s[2]);
			}
			for(int p = 0; p < Nperturbers; ++p){
				Gravity(xti, yti, zti, xTable_s, yTable_s, zTable_s, ax, ay, az, GM_s, p);
			}
			// ----------------------------------------------------------------------------

			kvx_d[ik] = ax;
			kvy_d[ik] = ay;
			kvz_d[ik] = az;
		}
		__syncthreads();
	}

	//update
	if(id < N){

		double dx = 0.0;
		double dy = 0.0;
		double dz = 0.0;

		double dvx = 0.0;
		double dvy = 0.0;
		double dvz = 0.0;

		for(int S = 0; S < RKFn; ++S){
			double dtb = dt * RKFb_c[S];
			int ik = id + S * N;
			dx += dtb * kx_d[ik];
			dy += dtb * ky_d[ik];
			dz += dtb * kz_d[ik];

			dvx += dtb * kvx_d[ik];
			dvy += dtb * kvy_d[ik];
			dvz += dtb * kvz_d[ik];
		}

		x_d[id] += dx;
		y_d[id] += dy;
		z_d[id] += dz;

		vx_d[id] += dvx;
		vy_d[id] += dvy;
		vz_d[id] += dvz;

//printf("%d %.20g %g %g %g\n", id, time, x_d[id], y_d[id], z_d[id]);
	}
}

// Runge Kutta step with fixed time step
// Every body runs on a sepparate thread block. Gravity calculation is spread along threads in the block and reuced in registers
// Kernel uses at least Nperturbers threads. The perturbers are loaded into shared memory
__global__ void RK_step2_kernel(double *xTable_d, double *yTable_d, double *zTable_d, double *vxTable_d, double *vyTable_d, double *vzTable_d, double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *GM_d, double *A1_d, double *A2_d, double *A3_d, double *Tsave_d, double *Rsave_d, double time, long long int timeStep, const double dt, const int RKFn, const int Nperturbers, const int N){


	int itx = threadIdx.x;
	int id = blockIdx.x;	//particle index

	//shared memory contains only the perturbers
	__shared__ double xTable_s[def_NP];
	__shared__ double yTable_s[def_NP];
	__shared__ double zTable_s[def_NP];

	__shared__ double vxTable_s[def_NP];
	__shared__ double vyTable_s[def_NP];
	__shared__ double vzTable_s[def_NP];

	__shared__ double kx_s[def_NP];
	__shared__ double ky_s[def_NP];
	__shared__ double kz_s[def_NP];

	__shared__ double kvx_s[def_NP];
	__shared__ double kvy_s[def_NP];
	__shared__ double kvz_s[def_NP];

	__shared__ double GM_s[def_NP];


	double xti;
	double yti;
	double zti;

	double vxti;
	double vyti;
	double vzti;

	double x0;
	double y0;
	double z0;

	double vx0;
	double vy0;
	double vz0;

	double A1;
	double A2;
	double A3;

	if(itx < Nperturbers){
		GM_s[itx] = GM_d[itx];
	}
	if(id < N){
		x0 = x_d[id];
		y0 = y_d[id];
		z0 = z_d[id];

		vx0 = vx_d[id];
		vy0 = vy_d[id];
		vz0 = vz_d[id];

		A1 = A1_d[id];
		A2 = A2_d[id];
		A3 = A3_d[id];
	}
	__syncthreads();


	for(int S = 0; S < RKFn; ++S){

		// ----------------------------------------------------------------------------
		//Read the perturbers position
		if(itx < Nperturbers){
			int ii = itx * RKFn + S;

			xTable_s[itx] = xTable_d[ii];
			yTable_s[itx] = yTable_d[ii];
			zTable_s[itx] = zTable_d[ii];

			vxTable_s[itx] = vxTable_d[ii];
			vyTable_s[itx] = vyTable_d[ii];
			vzTable_s[itx] = vzTable_d[ii];
//printf("%d %d %.20g %.20g %.20g\n", S, itx, xTable_s[itx], yTable_s[itx], zTable_s[itx]);
		}
		__syncthreads();

//if(id < Nperturbers){
//printf("p %d %.20g %.20g %.20g %.20g %.20g %.20g %.20g\n", id, time + RKFc_c[S] * dt, xt_s[id], yt_s[id], zt_s[id], vxt_s[id], vyt_s[id], vzt_s[id]);
//}

		if(id < N){

			xti = x0;
			yti = y0;
			zti = z0;

			vxti = vx0;
			vyti = vy0;
			vzti = vz0;

			double ax = 0.0;
			double ay = 0.0;
			double az = 0.0;

			for(int s = 0; s < S; ++s){
				double dtaa = dt * RKFa_c[S * RKFn + s];
				xti  += dtaa * kx_s[s];
				yti  += dtaa * ky_s[s];
				zti  += dtaa * kz_s[s];
				vxti += dtaa * kvx_s[s];
				vyti += dtaa * kvy_s[s];
				vzti += dtaa * kvz_s[s];
//printf("update 2 %d %d %g %g %g %g %g %g\n", S, id, xti, yti, zti, RKFa_c[S * RKFn + s], kx_s[s], dt);

			}

			if(itx == 0){
				kx_s[S] = vxti;
				ky_s[S] = vyti;
				kz_s[S] = vzti;


				// ----------------------------------------------------------------------------
				//compute forces

				//heliocentric coordinates
				double xih = xti - xTable_s[10];
				double yih = yti - yTable_s[10];
				double zih = zti - zTable_s[10];

				double vxih = vxti - vxTable_s[10];
				double vyih = vyti - vyTable_s[10];
				double vzih = vzti - vzTable_s[10];

				//r is used in multiple forces, so reuse it
				double rsq = __dmul_rn(xih, xih) + __dmul_rn(yih, yih) + __dmul_rn(zih, zih);
				double r = sqrt(rsq);

				if(cometFlag_c > 0 && S == 0){
					if(id == 0){
						Tsave_d[timeStep] = time;
					}
					Rsave_d[id * Rbuffersize_c + timeStep] = r;
				}

				//Earth centric coordinates
				double xiE = xti - xTable_s[2];
				double yiE = yti - yTable_s[2];
				double ziE = zti - zTable_s[2];

				if(useNonGrav_c == 1){
					NonGrav(xih, yih, zih, vxih, vyih, vzih, A1, A2, A3, r, ax, ay, az);
				}
				if(useGR_c == 1){
					GR(xih, yih, zih, vxih, vyih, vzih, r, ax, ay, az, GM_s[10]);
				}
				if(useJ2_c == 1){
					J2(xiE, yiE, ziE, ax, ay, az, GM_s[2]);
				}
			}

			__syncthreads();
			if(itx < Nperturbers){
				Gravity(xti, yti, zti, xTable_s, yTable_s, zTable_s, ax, ay, az, GM_s, itx);
			}
			__syncthreads();
			//reduce a
			// ----------------------------------------------------------------------------
			for(int i = 1; i < warpSize; i*=2){
				ax += __shfl_xor_sync(0xffffffff, ax, i, warpSize);
				ay += __shfl_xor_sync(0xffffffff, ay, i, warpSize);
				az += __shfl_xor_sync(0xffffffff, az, i, warpSize);
			}
			__syncthreads();

			if(itx == 0){
				kvx_s[S] = ax;
				kvy_s[S] = ay;
				kvz_s[S] = az;
			}
			__syncthreads();
		}
	}

	//update
	if(id < N && itx == 0){

		double dx = 0.0;
		double dy = 0.0;
		double dz = 0.0;

		double dvx = 0.0;
		double dvy = 0.0;
		double dvz = 0.0;

		for(int S = 0; S < RKFn; ++S){
			double dtb = dt * RKFb_c[S];
			dx += dtb * kx_s[S];
			dy += dtb * ky_s[S];
			dz += dtb * kz_s[S];

			dvx += dtb * kvx_s[S];
			dvy += dtb * kvy_s[S];
			dvz += dtb * kvz_s[S];
		}

		x_d[id] += dx;
		y_d[id] += dy;
		z_d[id] += dz;

		vx_d[id] += dvx;
		vy_d[id] += dvy;
		vz_d[id] += dvz;

//printf("%d %.20g %g %g %g\n", id, time, x_d[id], y_d[id], z_d[id]);
	}
}


//Runge Kutta Fehlberg step with adaptive time step
// Every body runs on a thread
// Kernel uses at least Nperturbers threads. The perturbers are loaded into shared memory
__global__ void RKF_step_kernel(double *xTable_d, double *yTable_d, double *zTable_d, double *vxTable_d, double *vyTable_d, double *vzTable_d, double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *dx_d, double *dy_d, double *dz_d, double *dvx_d, double *dvy_d, double *dvz_d, double *kx_d, double *ky_d, double *kz_d, double *kvx_d, double *kvy_d, double *kvz_d, double *GM_d, double *A1_d, double *A2_d, double *A3_d, int *index_d, double *Tsave_d, double *Rsave_d, double *snew_d, double time, long long int timeStep, const double dt, const int dts, const int RKFn, const int Nperturbers, const int N, const int N0, const int level, const int nL, double dtlimit){


	int itx = threadIdx.x;
	int jj = blockIdx.x * blockDim.x + threadIdx.x;

	int id = jj;

	if(level > 0){
		if(jj < N){
			id = index_d[(level - 1) * N0 + jj];
		}
	}


	//shared memory contains only the perturbers
	__shared__ double xTable_s[def_NP];
	__shared__ double yTable_s[def_NP];
	__shared__ double zTable_s[def_NP];

	__shared__ double vxTable_s[def_NP];
	__shared__ double vyTable_s[def_NP];
	__shared__ double vzTable_s[def_NP];

	__shared__ double GM_s[def_NP];


	double xti;
	double yti;
	double zti;

	double vxti;
	double vyti;
	double vzti;

	double x0;
	double y0;
	double z0;

	double vx0;
	double vy0;
	double vz0;

	double A1;
	double A2;
	double A3;

	if(itx < Nperturbers){
		GM_s[itx] = GM_d[itx];
	}
	if(jj < N){
		x0 = x_d[id];
		y0 = y_d[id];
		z0 = z_d[id];

		vx0 = vx_d[id];
		vy0 = vy_d[id];
		vz0 = vz_d[id];

		A1 = A1_d[id];
		A2 = A2_d[id];
		A3 = A3_d[id];
	}

	__syncthreads();


	for(int S = 0; S < RKFn; ++S){

		// ----------------------------------------------------------------------------
		//Read the perturbers position
		if(itx < Nperturbers){
			int ii = itx * RKFn + S;

			xTable_s[itx] = xTable_d[ii];
			yTable_s[itx] = yTable_d[ii];
			zTable_s[itx] = zTable_d[ii];

			vxTable_s[itx] = vxTable_d[ii];
			vyTable_s[itx] = vyTable_d[ii];
			vzTable_s[itx] = vzTable_d[ii];
//printf("%d %d %.20g %.20g %.20g\n", S, itx, xTable_s[itx], yTable_s[itx], zTable_s[itx]);
		}
		__syncthreads();

//if(id < Nperturbers){
//printf("p %d %.20g %.20g %.20g %.20g %.20g %.20g %.20g\n", id, time + RKFc_c[S] * dt, xt_s[id], yt_s[id], zt_s[id], vxt_s[id], vyt_s[id], vzt_s[id]);
//}

		if(jj < N){

			xti = x0;
			yti = y0;
			zti = z0;

			vxti = vx0;
			vyti = vy0;
			vzti = vz0;

			for(int s = 0; s < S; ++s){
				double dtaa = dt * RKFa_c[S * RKFn + s];
				int ii = jj + s * N;
				xti  += dtaa * kx_d[ii];
				yti  += dtaa * ky_d[ii];
				zti  += dtaa * kz_d[ii];
				vxti += dtaa * kvx_d[ii];
				vyti += dtaa * kvy_d[ii];
				vzti += dtaa * kvz_d[ii];
//printf("update 2 %d %d %g %g %g %g %g %g\n", S, id, xti, yti, zti, RKFa_c[S * RKFn + s], kx_d[s], dt);

			}

			int ik = jj + S * N;
			kx_d[ik] = vxti;
			ky_d[ik] = vyti;
			kz_d[ik] = vzti;

			// ----------------------------------------------------------------------------
			//compute forces
			double ax = 0.0;
			double ay = 0.0;
			double az = 0.0;

			//heliocentric coordinates
			double xih = xti - xTable_s[10];
			double yih = yti - yTable_s[10];
			double zih = zti - zTable_s[10];

			double vxih = vxti - vxTable_s[10];
			double vyih = vyti - vyTable_s[10];
			double vzih = vzti - vzTable_s[10];

			//r is used in multiple forces, so reuse it
			double rsq = __dmul_rn(xih, xih) + __dmul_rn(yih, yih) + __dmul_rn(zih, zih);
			double r = sqrt(rsq);

			if(cometFlag_c > 0 && S == 0){
				if(id == 0){
					Tsave_d[timeStep] = time;
				}
				Rsave_d[id * Rbuffersize_c + timeStep] = r;
			}

			//Earth centric coordinates
			double xiE = xti - xTable_s[2];
			double yiE = yti - yTable_s[2];
			double ziE = zti - zTable_s[2];

			if(useNonGrav_c == 1){
				NonGrav(xih, yih, zih, vxih, vyih, vzih, A1, A2, A3, r, ax, ay, az);
			}
			if(useGR_c == 1){
				GR(xih, yih, zih, vxih, vyih, vzih, r, ax, ay, az, GM_s[10]);
			}
			if(useJ2_c == 1){
				J2(xiE, yiE, ziE, ax, ay, az, GM_s[2]);
			}
			for(int p = 0; p < Nperturbers; ++p){
				Gravity(xti, yti, zti, xTable_s, yTable_s, zTable_s, ax, ay, az, GM_s, p);
			}
			// ----------------------------------------------------------------------------

			kvx_d[ik] = ax;
			kvy_d[ik] = ay;
			kvz_d[ik] = az;
		}
		__syncthreads();
	}

	//update

	double dx = 0.0;
	double dy = 0.0;
	double dz = 0.0;

	double dvx = 0.0;
	double dvy = 0.0;
	double dvz = 0.0;

	if(jj < N){

		for(int S = 0; S < RKFn; ++S){
			double dtb = dt * RKFb_c[S];
			int ii = jj + S * N;
			dx += dtb * kx_d[ii];
			dy += dtb * ky_d[ii];
			dz += dtb * kz_d[ii];

			dvx += dtb * kvx_d[ii];
			dvy += dtb * kvy_d[ii];
			dvz += dtb * kvz_d[ii];
		}

		dx_d[id] = dx;
		dy_d[id] = dy;
		dz_d[id] = dz;

		dvx_d[id] = dvx;
		dvy_d[id] = dvy;
		dvz_d[id] = dvz;

		//compute integration error

		double scalex  = atol_c + fabs(x0) * rtol_c;
		double scaley  = atol_c + fabs(y0) * rtol_c;
		double scalez  = atol_c + fabs(z0) * rtol_c;

		double scalevx = atol_c + fabs(vx0) * rtol_c;
		double scalevy = atol_c + fabs(vy0) * rtol_c;
		double scalevz = atol_c + fabs(vz0) * rtol_c;

		//error estimation
		double errorkx = 0.0;
		double errorky = 0.0;
		double errorkz = 0.0;
		
		double errorkvx = 0.0;
		double errorkvy = 0.0;
		double errorkvz = 0.0;

		for(int S = 0; S < RKFn; ++S){
			double f = (RKFb_c[S] - RKFbb_c[S]) * dt;
			int ii = jj + S * N;

			errorkx += __dmul_rn(f, kx_d[ii]);
			errorky += __dmul_rn(f, ky_d[ii]);
			errorkz += __dmul_rn(f, kz_d[ii]);

			errorkvx += __dmul_rn(f, kvx_d[ii]);
			errorkvy += __dmul_rn(f, kvy_d[ii]);
			errorkvz += __dmul_rn(f, kvz_d[ii]);

//printf("error %d %d %.20g %.20g %.20g %.20g %.20g %.20g %.20g\n", id, S, f, kx_d[ii], ky_d[ii], kz_d[ii], kvx_d[ii], kvy_d[ii], kvz_d[ii]);
//printf("error %d %d %.20g %.20g %.20g %.20g %.20g %.20g %.20g\n", id, S, f, errorkx, errorky, errorkz, errorkvx, errorkvy, errorkvz);
		}

		double errork = 0.0;
		errork += errorkx * errorkx / (scalex * scalex);
		errork += errorky * errorky / (scaley * scaley);
		errork += errorkz * errorkz / (scalez * scalez);
		errork += errorkvx * errorkvx / (scalevx * scalevx);
		errork += errorkvy * errorkvy / (scalevy * scalevy);
		errork += errorkvz * errorkvz / (scalevz * scalevz);

		errork = sqrt(errork / 6.0);	//6 is the number of dimensions

		double s = pow( 1.0  / errork, RKF_ee_c);

		s = (RKF_fac_c * s > RKF_facmin_c) ? RKF_fac_c * s : RKF_facmin_c;
		s = (RKF_facmax_c < s) ? RKF_facmax_c : s;
//printf("snew %d %d %g %g %g \n", level, id, s, dt, s * dt);

		if(s * dt * dts >= dtlimit || level >= nL - 1){
			snew_d[id] = s;
		}
		else{
			snew_d[id] = 1.0e6;      //mark body for higher level integration
		}
	}
}

// Runge Kutta Fehlberg step with adaptive time step
// Every body runs on a thread block with the gravitation of the perturbers on an individual thread
// Kernel uses at least Nperturbers threads. The perturbers are loaded into shared memory
__global__ void RKF_step2_kernel(double *xTable_d, double *yTable_d, double *zTable_d, double *vxTable_d, double *vyTable_d, double *vzTable_d, double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *dx_d, double *dy_d, double *dz_d, double *dvx_d, double *dvy_d, double *dvz_d,double *GM_d, double *A1_d, double *A2_d, double *A3_d, int *index_d, double *Tsave_d, double *Rsave_d, double *snew_d, double time, long long int timeStep, const double dt, const int dts, const int RKFn, const int Nperturbers, const int N, const int N0, const int level, const int nL, double dtlimit){


	int itx = threadIdx.x;
	int jj = blockIdx.x;

	int id = jj;		//particle id

	if(level > 0){
		id = index_d[(level - 1) * N0 + jj];
	}

	//shared memory contains only the perturbers
	__shared__ double xTable_s[def_NP];
	__shared__ double yTable_s[def_NP];
	__shared__ double zTable_s[def_NP];

	__shared__ double vxTable_s[def_NP];
	__shared__ double vyTable_s[def_NP];
	__shared__ double vzTable_s[def_NP];

	__shared__ double kx_s[def_NP];
	__shared__ double ky_s[def_NP];
	__shared__ double kz_s[def_NP];

	__shared__ double kvx_s[def_NP];
	__shared__ double kvy_s[def_NP];
	__shared__ double kvz_s[def_NP];


	__shared__ double GM_s[def_NP];


	double xti;
	double yti;
	double zti;

	double vxti;
	double vyti;
	double vzti;

	double x0 = x_d[id];
	double y0 = y_d[id];
	double z0 = z_d[id];

	double vx0 = vx_d[id];
	double vy0 = vy_d[id];
	double vz0 = vz_d[id];

	double A1 = A1_d[id];
	double A2 = A2_d[id];
	double A3 = A3_d[id];


	if(itx < Nperturbers){
		GM_s[itx] = GM_d[itx];
	}
	__syncthreads();


	for(int S = 0; S < RKFn; ++S){

		// ----------------------------------------------------------------------------
		//Read the perturbers position
		if(itx < Nperturbers){
			int ii = itx * RKFn + S;

			xTable_s[itx] = xTable_d[ii];
			yTable_s[itx] = yTable_d[ii];
			zTable_s[itx] = zTable_d[ii];

			vxTable_s[itx] = vxTable_d[ii];
			vyTable_s[itx] = vyTable_d[ii];
			vzTable_s[itx] = vzTable_d[ii];
//printf("%d %d %.20g %.20g %.20g\n", S, itx, xTable_s[itx], yTable_s[itx], zTable_s[itx]);
		}
		__syncthreads();

//if(id < Nperturbers){
//printf("p %d %.20g %.20g %.20g %.20g %.20g %.20g %.20g\n", id, time + RKFc_c[S] * dt, xt_s[id], yt_s[id], zt_s[id], vxt_s[id], vyt_s[id], vzt_s[id]);
//}

		xti = x0;
		yti = y0;
		zti = z0;

		vxti = vx0;
		vyti = vy0;
		vzti = vz0;

		double ax = 0.0;
		double ay = 0.0;
		double az = 0.0;

		for(int s = 0; s < S; ++s){
			double dtaa = dt * RKFa_c[S * RKFn + s];
			xti  += dtaa * kx_s[s];
			yti  += dtaa * ky_s[s];
			zti  += dtaa * kz_s[s];
			vxti += dtaa * kvx_s[s];
			vyti += dtaa * kvy_s[s];
			vzti += dtaa * kvz_s[s];
//printf("update 2 %d %d %g %g %g %g %g %g\n", S, id, xti, yti, zti, RKFa_c[S * RKFn + s], kx_d[s], dt);

		}

		if(itx == 0){
			kx_s[S] = vxti;
			ky_s[S] = vyti;
			kz_s[S] = vzti;

			// ----------------------------------------------------------------------------
			//compute forces

			//heliocentric coordinates
			double xih = xti - xTable_s[10];
			double yih = yti - yTable_s[10];
			double zih = zti - zTable_s[10];

			double vxih = vxti - vxTable_s[10];
			double vyih = vyti - vyTable_s[10];
			double vzih = vzti - vzTable_s[10];

			//r is used in multiple forces, so reuse it
			double rsq = __dmul_rn(xih, xih) + __dmul_rn(yih, yih) + __dmul_rn(zih, zih);
			double r = sqrt(rsq);

			if(cometFlag_c > 0 && S == 0){
				if(id == 0){
					Tsave_d[timeStep] = time;
				}
				Rsave_d[id * Rbuffersize_c + timeStep] = r;
			}

			//Earth centric coordinates
			double xiE = xti - xTable_s[2];
			double yiE = yti - yTable_s[2];
			double ziE = zti - zTable_s[2];

			if(useNonGrav_c == 1){
				NonGrav(xih, yih, zih, vxih, vyih, vzih, A1, A2, A3, r, ax, ay, az);
			}
			if(useGR_c == 1){
				GR(xih, yih, zih, vxih, vyih, vzih, r, ax, ay, az, GM_s[10]);
			}
			if(useJ2_c == 1){
				J2(xiE, yiE, ziE, ax, ay, az, GM_s[2]);
			}
		}

		__syncthreads();
	
		if(itx < Nperturbers){
			Gravity(xti, yti, zti, xTable_s, yTable_s, zTable_s, ax, ay, az, GM_s, itx);
		}
		__syncthreads();

		//reduce a
		// ----------------------------------------------------------------------------
		for(int i = 1; i < warpSize; i*=2){
			ax += __shfl_xor_sync(0xffffffff, ax, i, warpSize);
			ay += __shfl_xor_sync(0xffffffff, ay, i, warpSize);
			az += __shfl_xor_sync(0xffffffff, az, i, warpSize);
		}
		__syncthreads();

		if(itx == 0){
			kvx_s[S] = ax;
			kvy_s[S] = ay;
			kvz_s[S] = az;
		}
		__syncthreads();
	}

	//update

	double dx = 0.0;
	double dy = 0.0;
	double dz = 0.0;

	double dvx = 0.0;
	double dvy = 0.0;
	double dvz = 0.0;


	if(itx == 0){

		for(int S = 0; S < RKFn; ++S){
			double dtb = dt * RKFb_c[S];
			dx += dtb * kx_s[S];
			dy += dtb * ky_s[S];
			dz += dtb * kz_s[S];

			dvx += dtb * kvx_s[S];
			dvy += dtb * kvy_s[S];
			dvz += dtb * kvz_s[S];
		}

		dx_d[id] = dx;
		dy_d[id] = dy;
		dz_d[id] = dz;

		dvx_d[id] = dvx;
		dvy_d[id] = dvy;
		dvz_d[id] = dvz;


		//compute integration error
		double scalex  = atol_c + fabs(x0) * rtol_c;
		double scaley  = atol_c + fabs(y0) * rtol_c;
		double scalez  = atol_c + fabs(z0) * rtol_c;

		double scalevx = atol_c + fabs(vx0) * rtol_c;
		double scalevy = atol_c + fabs(vy0) * rtol_c;
		double scalevz = atol_c + fabs(vz0) * rtol_c;

		//error estimation
		double errorkx = 0.0;
		double errorky = 0.0;
		double errorkz = 0.0;
		
		double errorkvx = 0.0;
		double errorkvy = 0.0;
		double errorkvz = 0.0;

		for(int S = 0; S < RKFn; ++S){
			double f = (RKFb_c[S] - RKFbb_c[S]) * dt;
			errorkx += __dmul_rn(f, kx_s[S]);
			errorky += __dmul_rn(f, ky_s[S]);
			errorkz += __dmul_rn(f, kz_s[S]);

			errorkvx += __dmul_rn(f, kvx_s[S]);
			errorkvy += __dmul_rn(f, kvy_s[S]);
			errorkvz += __dmul_rn(f, kvz_s[S]);
//printf("error %d %d %.20g %.20g %.20g %.20g %.20g %.20g\n", id, S, errorkx, errorky, errorkz, errorkvx, errorkvy, errorkvz);
		}

		double errork = 0.0;
		errork += errorkx * errorkx / (scalex * scalex);
		errork += errorky * errorky / (scaley * scaley);
		errork += errorkz * errorkz / (scalez * scalez);
		errork += errorkvx * errorkvx / (scalevx * scalevx);
		errork += errorkvy * errorkvy / (scalevy * scalevy);
		errork += errorkvz * errorkvz / (scalevz * scalevz);

		errork = sqrt(errork / 6.0);	//6 is the number of dimensions

		double s = pow( 1.0  / errork, RKF_ee_c);

		s = (RKF_fac_c * s > RKF_facmin_c) ? RKF_fac_c * s : RKF_facmin_c;
		s = (RKF_facmax_c < s) ? RKF_facmax_c : s;

		snew_d[id] = s;
//printf("snew %d %d %g %g %g \n", level, id, snew_d[id], dt, s * dt);

		if(s * dt * dts >= dtlimit || level >= nL - 1){
			snew_d[id] = s;
		}
		else{
			snew_d[id] = 1.0e6;      //mark body for higher level integration
		}
	}
}


// Bulirsh-Stoer step with adaptive time step
// Every body runs on a thread
// Kernel uses at least Nperturbers threads. The perturbers are loaded into shared memory
__global__ void BS_step_kernel(double *xTable_d, double *yTable_d, double *zTable_d, double *vxTable_d, double *vyTable_d, double *vzTable_d, double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *dx_d, double *dy_d, double *dz_d, double *dvx_d, double *dvy_d, double *dvz_d, double *kx_d, double *ky_d, double *kz_d, double *kvx_d, double *kvy_d, double *kvz_d, double *GM_d, double *A1_d, double *A2_d, double *A3_d, int *index_d, double *Tsave_d, double *Rsave_d, double *snew_d, double time, long long int timeStep, const double dt, const int dts, const int BSn, const int Nperturbers, const int N, const int N0, const int level, const int nL, double dtlimit){


	int itx = threadIdx.x;
	int jj = blockIdx.x * blockDim.x + threadIdx.x;

	int id = jj;

	if(level > 0){
		if(jj < N){
			id = index_d[(level - 1) * N0 + jj];
		}
	}


	//shared memory contains only the perturbers
	__shared__ double xTable_s[def_NP];
	__shared__ double yTable_s[def_NP];
	__shared__ double zTable_s[def_NP];

	__shared__ double vxTable_s[def_NP];
	__shared__ double vyTable_s[def_NP];
	__shared__ double vzTable_s[def_NP];

	__shared__ double GM_s[def_NP];


	__shared__ double dx_s[def_N][8];
	__shared__ double dy_s[def_N][8];
	__shared__ double dz_s[def_N][8];

	__shared__ double dvx_s[def_N][8];
	__shared__ double dvy_s[def_N][8];
	__shared__ double dvz_s[def_N][8];


	__shared__ int f_s[1];

	double x0;
	double y0;
	double z0;

	double vx0;
	double vy0;
	double vz0;

	double xp;
	double yp;
	double zp;

	double vxp;
	double vyp;
	double vzp;

	double xt;
	double yt;
	double zt;

	double vxt;
	double vyt;
	double vzt;

	double A1;
	double A2;
	double A3;

	double ax;
	double ay;
	double az;

	double scalex;
	double scaley;
	double scalez;

	double scalevx;
	double scalevy;
	double scalevz;


	double snew = -1000.0;

	if(itx < Nperturbers){
		GM_s[itx] = GM_d[itx];
	}
	if(jj < N){
		x0 = x_d[id];
		y0 = y_d[id];
		z0 = z_d[id];

		vx0 = vx_d[id];
		vy0 = vy_d[id];
		vz0 = vz_d[id];

		A1 = A1_d[id];
		A2 = A2_d[id];
		A3 = A3_d[id];

		scalex = atol_c + fabs(x0) * rtol_c;
		scaley = atol_c + fabs(y0) * rtol_c;
		scalez = atol_c + fabs(z0) * rtol_c;

		scalevx = atol_c + fabs(vx0) * rtol_c;
		scalevy = atol_c + fabs(vy0) * rtol_c;
		scalevz = atol_c + fabs(vz0) * rtol_c;


	}

	__syncthreads();

	

	int cc = 0;
	for(int n = 1; n <= 8; ++n){
		double dt2 = dt / (2.0 * n);
		double dt22 = dt2 * 2.0;

		f_s[0] = 0;

		// ----------------------------------------------------------------------------
		//Read the perturbers position
		if(itx < Nperturbers){
			int ii = itx * BSn + cc;

			xTable_s[itx] = xTable_d[ii];
			yTable_s[itx] = yTable_d[ii];
			zTable_s[itx] = zTable_d[ii];

			vxTable_s[itx] = vxTable_d[ii];
			vyTable_s[itx] = vyTable_d[ii];
			vzTable_s[itx] = vzTable_d[ii];
//printf("%d %d %.20g %.20g %.20g\n", S, itx, xTable_s[itx], yTable_s[itx], zTable_s[itx]);
		}
		++cc;
		// ----------------------------------------------------------------------------
		__syncthreads();

//if(id < Nperturbers){
//printf("p %d %.20g %.20g %.20g %.20g %.20g %.20g %.20g\n", id, time + RKFc_c[S] * dt, xt_s[id], yt_s[id], zt_s[id], vxt_s[id], vyt_s[id], vzt_s[id]);
//}

		if(jj < N && snew < 0.0){
			// ----------------------------------------------------------------------------
			//compute forces
			ax = 0.0;
			ay = 0.0;
			az = 0.0;

			//heliocentric coordinates
			double xih = x0 - xTable_s[10];
			double yih = y0 - yTable_s[10];
			double zih = z0 - zTable_s[10];

			double vxih = vx0 - vxTable_s[10];
			double vyih = vy0 - vyTable_s[10];
			double vzih = vz0 - vzTable_s[10];

			//r is used in multiple forces, so reuse it
			double rsq = __dmul_rn(xih, xih) + __dmul_rn(yih, yih) + __dmul_rn(zih, zih);
			double r = sqrt(rsq);

			if(cometFlag_c > 0 && n == 0){
				if(id == 0){
					Tsave_d[timeStep] = time;
				}
				Rsave_d[id * Rbuffersize_c + timeStep] = r;
			}

			//Earth centric coordinates
			double xiE = x0 - xTable_s[2];
			double yiE = y0 - yTable_s[2];
			double ziE = z0 - zTable_s[2];

			if(useNonGrav_c == 1){
				NonGrav(xih, yih, zih, vxih, vyih, vzih, A1, A2, A3, r, ax, ay, az);
			}
			if(useGR_c == 1){
				GR(xih, yih, zih, vxih, vyih, vzih, r, ax, ay, az, GM_s[10]);
			}
			if(useJ2_c == 1){
				J2(xiE, yiE, ziE, ax, ay, az, GM_s[2]);
			}
			for(int p = 0; p < Nperturbers; ++p){
				Gravity(x0, y0, z0, xTable_s, yTable_s, zTable_s, ax, ay, az, GM_s, p);
			}
			// ----------------------------------------------------------------------------

			xp = x0 + dt2 * vx0;
			yp = y0 + dt2 * vy0;
			zp = z0 + dt2 * vz0;

			vxp = vx0 + dt2 * ax;
			vyp = vy0 + dt2 * ay;
			vzp = vz0 + dt2 * az;

		}
		__syncthreads();


		// ----------------------------------------------------------------------------
		//Read the perturbers position
		if(itx < Nperturbers){
			int ii = itx * BSn + cc;

			xTable_s[itx] = xTable_d[ii];
			yTable_s[itx] = yTable_d[ii];
			zTable_s[itx] = zTable_d[ii];

			vxTable_s[itx] = vxTable_d[ii];
			vyTable_s[itx] = vyTable_d[ii];
			vzTable_s[itx] = vzTable_d[ii];
		}
		++cc;
		// ----------------------------------------------------------------------------
		__syncthreads();

		if(jj < N && snew < 0.0){

			// ----------------------------------------------------------------------------
			//compute forces
			ax = 0.0;
			ay = 0.0;
			az = 0.0;

			//heliocentric coordinates
			double xih = xp - xTable_s[10];
			double yih = yp - yTable_s[10];
			double zih = zp - zTable_s[10];

			double vxih = vxp - vxTable_s[10];
			double vyih = vyp - vyTable_s[10];
			double vzih = vzp - vzTable_s[10];

			//r is used in multiple forces, so reuse it
			double rsq = __dmul_rn(xih, xih) + __dmul_rn(yih, yih) + __dmul_rn(zih, zih);
			double r = sqrt(rsq);

			//Earth centric coordinates
			double xiE = xp - xTable_s[2];
			double yiE = yp - yTable_s[2];
			double ziE = zp - zTable_s[2];

			if(useNonGrav_c == 1){
				NonGrav(xih, yih, zih, vxih, vyih, vzih, A1, A2, A3, r, ax, ay, az);
			}
			if(useGR_c == 1){
				GR(xih, yih, zih, vxih, vyih, vzih, r, ax, ay, az, GM_s[10]);
			}
			if(useJ2_c == 1){
				J2(xiE, yiE, ziE, ax, ay, az, GM_s[2]);
			}
			for(int p = 0; p < Nperturbers; ++p){
				Gravity(xp, yp, zp, xTable_s, yTable_s, zTable_s, ax, ay, az, GM_s, p);
			}
			// ----------------------------------------------------------------------------

			xt = x0 + dt22 * vxp;
			yt = y0 + dt22 * vyp;
			zt = z0 + dt22 * vzp;

			vxt = vx0 + dt22 * ax;
			vyt = vy0 + dt22 * ay;
			vzt = vz0 + dt22 * az;

		}
		__syncthreads();

		for(int m = 2; m <= n; ++m){
			// ----------------------------------------------------------------------------
			//Read the perturbers position
			if(itx < Nperturbers){
				int ii = itx * BSn + cc;

				xTable_s[itx] = xTable_d[ii];
				yTable_s[itx] = yTable_d[ii];
				zTable_s[itx] = zTable_d[ii];

				vxTable_s[itx] = vxTable_d[ii];
				vyTable_s[itx] = vyTable_d[ii];
				vzTable_s[itx] = vzTable_d[ii];
			}
			++cc;
			// ----------------------------------------------------------------------------
			__syncthreads();


			if(jj < N && snew < 0.0){

				// ----------------------------------------------------------------------------
				//compute forces
				ax = 0.0;
				ay = 0.0;
				az = 0.0;

				//heliocentric coordinates
				double xih = xt - xTable_s[10];
				double yih = yt - yTable_s[10];
				double zih = zt - zTable_s[10];

				double vxih = vxt - vxTable_s[10];
				double vyih = vyt - vyTable_s[10];
				double vzih = vzt - vzTable_s[10];

				//r is used in multiple forces, so reuse it
				double rsq = __dmul_rn(xih, xih) + __dmul_rn(yih, yih) + __dmul_rn(zih, zih);
				double r = sqrt(rsq);


				//Earth centric coordinates
				double xiE = xt - xTable_s[2];
				double yiE = yt - yTable_s[2];
				double ziE = zt - zTable_s[2];

				if(useNonGrav_c == 1){
					NonGrav(xih, yih, zih, vxih, vyih, vzih, A1, A2, A3, r, ax, ay, az);
				}
				if(useGR_c == 1){
					GR(xih, yih, zih, vxih, vyih, vzih, r, ax, ay, az, GM_s[10]);
				}
				if(useJ2_c == 1){
					J2(xiE, yiE, ziE, ax, ay, az, GM_s[2]);
				}
				for(int p = 0; p < Nperturbers; ++p){
					Gravity(xt, yt, zt, xTable_s, yTable_s, zTable_s, ax, ay, az, GM_s, p);
				}
				// ----------------------------------------------------------------------------

				xp += dt22 * vxt;
				yp += dt22 * vyt;
				zp += dt22 * vzt;

				vxp += dt22 * ax;
				vyp += dt22 * ay;
				vzp += dt22 * az;

			}
			__syncthreads();

			// ----------------------------------------------------------------------------
			//Read the perturbers position
			if(itx < Nperturbers){
				int ii = itx * BSn + cc;

				xTable_s[itx] = xTable_d[ii];
				yTable_s[itx] = yTable_d[ii];
				zTable_s[itx] = zTable_d[ii];

				vxTable_s[itx] = vxTable_d[ii];
				vyTable_s[itx] = vyTable_d[ii];
				vzTable_s[itx] = vzTable_d[ii];
			}
			++cc;
			// ----------------------------------------------------------------------------
			__syncthreads();


			if(jj < N && snew < 0.0){

				// ----------------------------------------------------------------------------
				//compute forces
				ax = 0.0;
				ay = 0.0;
				az = 0.0;

				//heliocentric coordinates
				double xih = xp - xTable_s[10];
				double yih = yp - yTable_s[10];
				double zih = zp - zTable_s[10];

				double vxih = vxp - vxTable_s[10];
				double vyih = vyp - vyTable_s[10];
				double vzih = vzp - vzTable_s[10];

				//r is used in multiple forces, so reuse it
				double rsq = __dmul_rn(xih, xih) + __dmul_rn(yih, yih) + __dmul_rn(zih, zih);
				double r = sqrt(rsq);


				//Earth centric coordinates
				double xiE = xp - xTable_s[2];
				double yiE = yp - yTable_s[2];
				double ziE = zp - zTable_s[2];

				if(useNonGrav_c == 1){
					NonGrav(xih, yih, zih, vxih, vyih, vzih, A1, A2, A3, r, ax, ay, az);
				}
				if(useGR_c == 1){
					GR(xih, yih, zih, vxih, vyih, vzih, r, ax, ay, az, GM_s[10]);
				}
				if(useJ2_c == 1){
					J2(xiE, yiE, ziE, ax, ay, az, GM_s[2]);
				}
				for(int p = 0; p < Nperturbers; ++p){
					Gravity(xp, yp, zp, xTable_s, yTable_s, zTable_s, ax, ay, az, GM_s, p);
				}
				// ----------------------------------------------------------------------------

				xt += dt22 * vxp;
				yt += dt22 * vyp;
				zt += dt22 * vzp;

				vxt += dt22 * ax;
				vyt += dt22 * ay;
				vzt += dt22 * az;

			}
			__syncthreads();

		} // end of m loop

		// ----------------------------------------------------------------------------
		//Read the perturbers position
		if(itx < Nperturbers){
			int ii = itx * BSn + cc;

			xTable_s[itx] = xTable_d[ii];
			yTable_s[itx] = yTable_d[ii];
			zTable_s[itx] = zTable_d[ii];

			vxTable_s[itx] = vxTable_d[ii];
			vyTable_s[itx] = vyTable_d[ii];
			vzTable_s[itx] = vzTable_d[ii];
		}
		++cc;
		// ----------------------------------------------------------------------------
		__syncthreads();

		double error = 0.0;

		if(jj < N && snew < 0.0){

			// ----------------------------------------------------------------------------
			//compute forces
			ax = 0.0;
			ay = 0.0;
			az = 0.0;

			//heliocentric coordinates
			double xih = xt - xTable_s[10];
			double yih = yt - yTable_s[10];
			double zih = zt - zTable_s[10];

			double vxih = vxt - vxTable_s[10];
			double vyih = vyt - vyTable_s[10];
			double vzih = vzt - vzTable_s[10];

			//r is used in multiple forces, so reuse it
			double rsq = __dmul_rn(xih, xih) + __dmul_rn(yih, yih) + __dmul_rn(zih, zih);
			double r = sqrt(rsq);


			//Earth centric coordinates
			double xiE = xt - xTable_s[2];
			double yiE = yt - yTable_s[2];
			double ziE = zt - zTable_s[2];

			if(useNonGrav_c == 1){
				NonGrav(xih, yih, zih, vxih, vyih, vzih, A1, A2, A3, r, ax, ay, az);
			}
			if(useGR_c == 1){
				GR(xih, yih, zih, vxih, vyih, vzih, r, ax, ay, az, GM_s[10]);
			}
			if(useJ2_c == 1){
				J2(xiE, yiE, ziE, ax, ay, az, GM_s[2]);
			}
			for(int p = 0; p < Nperturbers; ++p){
				Gravity(xt, yt, zt, xTable_s, yTable_s, zTable_s, ax, ay, az, GM_s, p);
			}
			// ----------------------------------------------------------------------------

			xp += dt2 * vxt;
			yp += dt2 * vyt;
			zp += dt2 * vzt;

			vxp += dt2 * ax;
			vyp += dt2 * ay;
			vzp += dt2 * az;


			dx_s[itx][n-1] = 0.5 * (xt + xp);
			dy_s[itx][n-1] = 0.5 * (yt + yp);
			dz_s[itx][n-1] = 0.5 * (zt + zp);

			dvx_s[itx][n-1] = 0.5 * (vxt + vxp);
			dvy_s[itx][n-1] = 0.5 * (vyt + vyp);
			dvz_s[itx][n-1] = 0.5 * (vzt + vzp);

			for(int j = n - 1; j >= 1; --j){
				double t0 = BSt0_c[(n-1) * 8 + (j-1)];
				double t1 = t0 * BSddt_c[j];
				double t2 = t0 * BSddt_c[n-1];

				dx_s[itx][j-1] = (t1 * dx_s[itx][j]) - (t2 * dx_s[itx][j-1]);
				dy_s[itx][j-1] = (t1 * dy_s[itx][j]) - (t2 * dy_s[itx][j-1]);
				dz_s[itx][j-1] = (t1 * dz_s[itx][j]) - (t2 * dz_s[itx][j-1]);

				dvx_s[itx][j-1] = (t1 * dvx_s[itx][j]) - (t2 * dvx_s[itx][j-1]);
				dvy_s[itx][j-1] = (t1 * dvy_s[itx][j]) - (t2 * dvy_s[itx][j-1]);
				dvz_s[itx][j-1] = (t1 * dvz_s[itx][j]) - (t2 * dvz_s[itx][j-1]);

			}

			double error1 = 0.0;

			error1 = dx_s[itx][0] / scalex;
			error += error1 * error1;

			error1 = dy_s[itx][0] / scaley;
			error += error1 * error1;

			error1 = dz_s[itx][0] / scalez;
			error += error1 * error1;

			error1 = dvx_s[itx][0] / scalevx;
			error += error1 * error1;

			error1 = dvy_s[itx][0] / scalevy;
			error += error1 * error1;

			error1 = dvz_s[itx][0] / scalevz;
			error += error1 * error1;


			error = sqrt(error / 6.0);      //6 is the number of dimensions
//printf("error %d %d %d %.20g\n", level, id, n, error);


			if(error < 1.0){
				//update
				xt = dx_s[itx][0];
				yt = dy_s[itx][0];
				zt = dz_s[itx][0];

				vxt = dvx_s[itx][0];
				vyt = dvy_s[itx][0];
				vzt = dvz_s[itx][0];

				for(int j = 1; j < n; ++j){
					xt += dx_s[itx][j];
					yt += dy_s[itx][j];
					zt += dz_s[itx][j];

					vxt += dvx_s[itx][j];
					vyt += dvy_s[itx][j];
					vzt += dvz_s[itx][j];
				}
				dx_d[id] = xt;
				dy_d[id] = yt;
				dz_d[id] = zt;

				dvx_d[id] = vxt;
				dvy_d[id] = vyt;
				dvz_d[id] = vzt;


				snew = 1.0;

				if(n >= 8){
					snew = 0.55;
				}
				if(n < 7){
					snew = 1.3;
				}
//printf("Accecpt %d %d %d %g %g\n", level, id, n, dt, snew);

			}
			else{
				f_s[0] = 1;
				if(n == 8){
					snew = 0.5;
//printf("Repeat %d %d %d %g %g\n", level, id, n, dt, snew);
				}
			}


		}
		__syncthreads();

		if(f_s[0] == 0){
			break;
		}


	}//end of n loop

	if(jj < N){
		if(snew * dt * dts >= dtlimit || level >= nL - 1){
			snew_d[id] = snew;
		}
		else{
			snew_d[id] = 1.0e6;      //mark body for higher level integration
		}
//printf("snew %d %d %g\n", level, id, snew_d[id]); 
	}

}


// Bulirsh-Stoer step with adaptive time step
// Every body runs on a thread block with the gravitation of the perturbers on an individual thread
// Kernel uses at least Nperturbers threads. The perturbers are loaded into shared memory
__global__ void BS_step2_kernel(double *xTable_d, double *yTable_d, double *zTable_d, double *vxTable_d, double *vyTable_d, double *vzTable_d, double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *dx_d, double *dy_d, double *dz_d, double *dvx_d, double *dvy_d, double *dvz_d, double *kx_d, double *ky_d, double *kz_d, double *kvx_d, double *kvy_d, double *kvz_d, double *GM_d, double *A1_d, double *A2_d, double *A3_d, int *index_d, double *Tsave_d, double *Rsave_d, double *snew_d, double time, long long int timeStep, const double dt, const int dts, const int BSn, const int Nperturbers, const int N, const int N0, const int level, const int nL, double dtlimit){


	int itx = threadIdx.x;
	int jj = blockIdx.x;

	int id = jj;		//particle id

	if(level > 0){
		id = index_d[(level - 1) * N0 + jj];
	}


	//shared memory contains only the perturbers
	__shared__ double xTable_s[def_NP];
	__shared__ double yTable_s[def_NP];
	__shared__ double zTable_s[def_NP];

	__shared__ double vxTable_s[def_NP];
	__shared__ double vyTable_s[def_NP];
	__shared__ double vzTable_s[def_NP];

	__shared__ double GM_s[def_NP];


	__shared__ double dx_s[8];
	__shared__ double dy_s[8];
	__shared__ double dz_s[8];

	__shared__ double dvx_s[8];
	__shared__ double dvy_s[8];
	__shared__ double dvz_s[8];

	__shared__ int f_s[1];

	double xp;
	double yp;
	double zp;

	double vxp;
	double vyp;
	double vzp;

	double xt;
	double yt;
	double zt;

	double vxt;
	double vyt;
	double vzt;

	double ax;
	double ay;
	double az;

	double x0 = x_d[id];
	double y0 = y_d[id];
	double z0 = z_d[id];

	double vx0 = vx_d[id];
	double vy0 = vy_d[id];
	double vz0 = vz_d[id];

	double A1 = A1_d[id];
	double A2 = A2_d[id];
	double A3 = A3_d[id];

	double scalex = atol_c + fabs(x0) * rtol_c;
	double scaley = atol_c + fabs(y0) * rtol_c;
	double scalez = atol_c + fabs(z0) * rtol_c;

	double scalevx = atol_c + fabs(vx0) * rtol_c;
	double scalevy = atol_c + fabs(vy0) * rtol_c;
	double scalevz = atol_c + fabs(vz0) * rtol_c;


	double snew = -1000.0;

	if(itx < Nperturbers){
		GM_s[itx] = GM_d[itx];
	}



	__syncthreads();

	

	int cc = 0;
	for(int n = 1; n <= 8; ++n){
		double dt2 = dt / (2.0 * n);
		double dt22 = dt2 * 2.0;

		f_s[0] = 0;

		// ----------------------------------------------------------------------------
		//Read the perturbers position
		if(itx < Nperturbers){
			int ii = itx * BSn + cc;

			xTable_s[itx] = xTable_d[ii];
			yTable_s[itx] = yTable_d[ii];
			zTable_s[itx] = zTable_d[ii];

			vxTable_s[itx] = vxTable_d[ii];
			vyTable_s[itx] = vyTable_d[ii];
			vzTable_s[itx] = vzTable_d[ii];
//printf("%d %d %.20g %.20g %.20g\n", S, itx, xTable_s[itx], yTable_s[itx], zTable_s[itx]);
		}
		++cc;
		// ----------------------------------------------------------------------------
		__syncthreads();

//if(id < Nperturbers){
//printf("p %d %.20g %.20g %.20g %.20g %.20g %.20g %.20g\n", id, time + RKFc_c[S] * dt, xt_s[id], yt_s[id], zt_s[id], vxt_s[id], vyt_s[id], vzt_s[id]);
//}
		ax = 0.0;
		ay = 0.0;
		az = 0.0;

		if(itx == 0){

			// ----------------------------------------------------------------------------
			//compute forces

			//heliocentric coordinates
			double xih = x0 - xTable_s[10];
			double yih = y0 - yTable_s[10];
			double zih = z0 - zTable_s[10];

			double vxih = vx0 - vxTable_s[10];
			double vyih = vy0 - vyTable_s[10];
			double vzih = vz0 - vzTable_s[10];

			//r is used in multiple forces, so reuse it
			double rsq = __dmul_rn(xih, xih) + __dmul_rn(yih, yih) + __dmul_rn(zih, zih);
			double r = sqrt(rsq);

			if(cometFlag_c > 0 && n == 0){
				if(id == 0){
					Tsave_d[timeStep] = time;
				}
				Rsave_d[id * Rbuffersize_c + timeStep] = r;
			}

			//Earth centric coordinates
			double xiE = x0 - xTable_s[2];
			double yiE = y0 - yTable_s[2];
			double ziE = z0 - zTable_s[2];

			if(useNonGrav_c == 1){
				NonGrav(xih, yih, zih, vxih, vyih, vzih, A1, A2, A3, r, ax, ay, az);
			}
			if(useGR_c == 1){
				GR(xih, yih, zih, vxih, vyih, vzih, r, ax, ay, az, GM_s[10]);
			}
			if(useJ2_c == 1){
				J2(xiE, yiE, ziE, ax, ay, az, GM_s[2]);
			}

		}
		__syncthreads();

		if(itx < Nperturbers){
			Gravity(x0, y0, z0, xTable_s, yTable_s, zTable_s, ax, ay, az, GM_s, itx);
		}
		__syncthreads();

		//reduce a
		// ----------------------------------------------------------------------------
		for(int i = 1; i < warpSize; i*=2){
			ax += __shfl_xor_sync(0xffffffff, ax, i, warpSize);
			ay += __shfl_xor_sync(0xffffffff, ay, i, warpSize);
			az += __shfl_xor_sync(0xffffffff, az, i, warpSize);
		}
		__syncthreads();
		// inside a warp all threads have now the same value of ax ay and az

		if(itx < Nperturbers){
			xp = x0 + dt2 * vx0;
			yp = y0 + dt2 * vy0;
			zp = z0 + dt2 * vz0;

			vxp = vx0 + dt2 * ax;
			vyp = vy0 + dt2 * ay;
			vzp = vz0 + dt2 * az;
		}
		__syncthreads();


		// ----------------------------------------------------------------------------
		//Read the perturbers position
		if(itx < Nperturbers){
			int ii = itx * BSn + cc;

			xTable_s[itx] = xTable_d[ii];
			yTable_s[itx] = yTable_d[ii];
			zTable_s[itx] = zTable_d[ii];

			vxTable_s[itx] = vxTable_d[ii];
			vyTable_s[itx] = vyTable_d[ii];
			vzTable_s[itx] = vzTable_d[ii];
		}
		++cc;
		// ----------------------------------------------------------------------------
		__syncthreads();

		ax = 0.0;
		ay = 0.0;
		az = 0.0;

		if(itx == 0){

			// ----------------------------------------------------------------------------
			//compute forces

			//heliocentric coordinates
			double xih = xp - xTable_s[10];
			double yih = yp - yTable_s[10];
			double zih = zp - zTable_s[10];

			double vxih = vxp - vxTable_s[10];
			double vyih = vyp - vyTable_s[10];
			double vzih = vzp - vzTable_s[10];

			//r is used in multiple forces, so reuse it
			double rsq = __dmul_rn(xih, xih) + __dmul_rn(yih, yih) + __dmul_rn(zih, zih);
			double r = sqrt(rsq);

			//Earth centric coordinates
			double xiE = xp - xTable_s[2];
			double yiE = yp - yTable_s[2];
			double ziE = zp - zTable_s[2];

			if(useNonGrav_c == 1){
				NonGrav(xih, yih, zih, vxih, vyih, vzih, A1, A2, A3, r, ax, ay, az);
			}
			if(useGR_c == 1){
				GR(xih, yih, zih, vxih, vyih, vzih, r, ax, ay, az, GM_s[10]);
			}
			if(useJ2_c == 1){
				J2(xiE, yiE, ziE, ax, ay, az, GM_s[2]);
			}

		}
		__syncthreads();

		if(itx < Nperturbers){
			Gravity(xp, yp, zp, xTable_s, yTable_s, zTable_s, ax, ay, az, GM_s, itx);
		}
		__syncthreads();

		//reduce a
		// ----------------------------------------------------------------------------
		for(int i = 1; i < warpSize; i*=2){
			ax += __shfl_xor_sync(0xffffffff, ax, i, warpSize);
			ay += __shfl_xor_sync(0xffffffff, ay, i, warpSize);
			az += __shfl_xor_sync(0xffffffff, az, i, warpSize);
		}
		__syncthreads();


		if(itx < Nperturbers){
			xt = x0 + dt22 * vxp;
			yt = y0 + dt22 * vyp;
			zt = z0 + dt22 * vzp;

			vxt = vx0 + dt22 * ax;
			vyt = vy0 + dt22 * ay;
			vzt = vz0 + dt22 * az;
		}
		__syncthreads();

		for(int m = 2; m <= n; ++m){
			// ----------------------------------------------------------------------------
			//Read the perturbers position
			if(itx < Nperturbers){
				int ii = itx * BSn + cc;

				xTable_s[itx] = xTable_d[ii];
				yTable_s[itx] = yTable_d[ii];
				zTable_s[itx] = zTable_d[ii];

				vxTable_s[itx] = vxTable_d[ii];
				vyTable_s[itx] = vyTable_d[ii];
				vzTable_s[itx] = vzTable_d[ii];
			}
			++cc;
			// ----------------------------------------------------------------------------
			__syncthreads();

			ax = 0.0;
			ay = 0.0;
			az = 0.0;

			if(itx == 0){

				// ----------------------------------------------------------------------------
				//compute forces

				//heliocentric coordinates
				double xih = xt - xTable_s[10];
				double yih = yt - yTable_s[10];
				double zih = zt - zTable_s[10];

				double vxih = vxt - vxTable_s[10];
				double vyih = vyt - vyTable_s[10];
				double vzih = vzt - vzTable_s[10];

				//r is used in multiple forces, so reuse it
				double rsq = __dmul_rn(xih, xih) + __dmul_rn(yih, yih) + __dmul_rn(zih, zih);
				double r = sqrt(rsq);


				//Earth centric coordinates
				double xiE = xt - xTable_s[2];
				double yiE = yt - yTable_s[2];
				double ziE = zt - zTable_s[2];

				if(useNonGrav_c == 1){
					NonGrav(xih, yih, zih, vxih, vyih, vzih, A1, A2, A3, r, ax, ay, az);
				}
				if(useGR_c == 1){
					GR(xih, yih, zih, vxih, vyih, vzih, r, ax, ay, az, GM_s[10]);
				}
				if(useJ2_c == 1){
					J2(xiE, yiE, ziE, ax, ay, az, GM_s[2]);
				}
			}
			__syncthreads();

			if(itx < Nperturbers){
				Gravity(xt, yt, zt, xTable_s, yTable_s, zTable_s, ax, ay, az, GM_s, itx);
			}
			__syncthreads();

			//reduce a
			// ----------------------------------------------------------------------------
			for(int i = 1; i < warpSize; i*=2){
				ax += __shfl_xor_sync(0xffffffff, ax, i, warpSize);
				ay += __shfl_xor_sync(0xffffffff, ay, i, warpSize);
				az += __shfl_xor_sync(0xffffffff, az, i, warpSize);
			}
			__syncthreads();


			if(itx < Nperturbers){
				xp += dt22 * vxt;
				yp += dt22 * vyt;
				zp += dt22 * vzt;

				vxp += dt22 * ax;
				vyp += dt22 * ay;
				vzp += dt22 * az;
			}
			__syncthreads();

			// ----------------------------------------------------------------------------
			//Read the perturbers position
			if(itx < Nperturbers){
				int ii = itx * BSn + cc;

				xTable_s[itx] = xTable_d[ii];
				yTable_s[itx] = yTable_d[ii];
				zTable_s[itx] = zTable_d[ii];

				vxTable_s[itx] = vxTable_d[ii];
				vyTable_s[itx] = vyTable_d[ii];
				vzTable_s[itx] = vzTable_d[ii];
			}
			++cc;
			// ----------------------------------------------------------------------------
			__syncthreads();


			ax = 0.0;
			ay = 0.0;
			az = 0.0;

			if(itx == 0){

				// ----------------------------------------------------------------------------
				//compute forces

				//heliocentric coordinates
				double xih = xp - xTable_s[10];
				double yih = yp - yTable_s[10];
				double zih = zp - zTable_s[10];

				double vxih = vxp - vxTable_s[10];
				double vyih = vyp - vyTable_s[10];
				double vzih = vzp - vzTable_s[10];

				//r is used in multiple forces, so reuse it
				double rsq = __dmul_rn(xih, xih) + __dmul_rn(yih, yih) + __dmul_rn(zih, zih);
				double r = sqrt(rsq);


				//Earth centric coordinates
				double xiE = xp - xTable_s[2];
				double yiE = yp - yTable_s[2];
				double ziE = zp - zTable_s[2];

				if(useNonGrav_c == 1){
					NonGrav(xih, yih, zih, vxih, vyih, vzih, A1, A2, A3, r, ax, ay, az);
				}
				if(useGR_c == 1){
					GR(xih, yih, zih, vxih, vyih, vzih, r, ax, ay, az, GM_s[10]);
				}
				if(useJ2_c == 1){
					J2(xiE, yiE, ziE, ax, ay, az, GM_s[2]);
				}
			}
			__syncthreads();

			if(itx < Nperturbers){
				Gravity(xp, yp, zp, xTable_s, yTable_s, zTable_s, ax, ay, az, GM_s, itx);
			}
			__syncthreads();

			//reduce a
			// ----------------------------------------------------------------------------
			for(int i = 1; i < warpSize; i*=2){
				ax += __shfl_xor_sync(0xffffffff, ax, i, warpSize);
				ay += __shfl_xor_sync(0xffffffff, ay, i, warpSize);
				az += __shfl_xor_sync(0xffffffff, az, i, warpSize);
			}
			__syncthreads();


			if(itx < Nperturbers){
				xt += dt22 * vxp;
				yt += dt22 * vyp;
				zt += dt22 * vzp;

				vxt += dt22 * ax;
				vyt += dt22 * ay;
				vzt += dt22 * az;
			}
			__syncthreads();

		} // end of m loop

		// ----------------------------------------------------------------------------
		//Read the perturbers position
		if(itx < Nperturbers){
			int ii = itx * BSn + cc;

			xTable_s[itx] = xTable_d[ii];
			yTable_s[itx] = yTable_d[ii];
			zTable_s[itx] = zTable_d[ii];

			vxTable_s[itx] = vxTable_d[ii];
			vyTable_s[itx] = vyTable_d[ii];
			vzTable_s[itx] = vzTable_d[ii];
		}
		++cc;
		// ----------------------------------------------------------------------------
		__syncthreads();

		double error = 0.0;

		ax = 0.0;
		ay = 0.0;
		az = 0.0;

		if(itx == 0){

			// ----------------------------------------------------------------------------
			//compute forces

			//heliocentric coordinates
			double xih = xt - xTable_s[10];
			double yih = yt - yTable_s[10];
			double zih = zt - zTable_s[10];

			double vxih = vxt - vxTable_s[10];
			double vyih = vyt - vyTable_s[10];
			double vzih = vzt - vzTable_s[10];

			//r is used in multiple forces, so reuse it
			double rsq = __dmul_rn(xih, xih) + __dmul_rn(yih, yih) + __dmul_rn(zih, zih);
			double r = sqrt(rsq);


			//Earth centric coordinates
			double xiE = xt - xTable_s[2];
			double yiE = yt - yTable_s[2];
			double ziE = zt - zTable_s[2];

			if(useNonGrav_c == 1){
				NonGrav(xih, yih, zih, vxih, vyih, vzih, A1, A2, A3, r, ax, ay, az);
			}
			if(useGR_c == 1){
				GR(xih, yih, zih, vxih, vyih, vzih, r, ax, ay, az, GM_s[10]);
			}
			if(useJ2_c == 1){
				J2(xiE, yiE, ziE, ax, ay, az, GM_s[2]);
			}
		}
		__syncthreads();

		if(itx < Nperturbers){
			Gravity(xt, yt, zt, xTable_s, yTable_s, zTable_s, ax, ay, az, GM_s, itx);
		}
		__syncthreads();

		//reduce a
		// ----------------------------------------------------------------------------
		for(int i = 1; i < warpSize; i*=2){
			ax += __shfl_xor_sync(0xffffffff, ax, i, warpSize);
			ay += __shfl_xor_sync(0xffffffff, ay, i, warpSize);
			az += __shfl_xor_sync(0xffffffff, az, i, warpSize);
		}
		__syncthreads();


		if(itx < Nperturbers){
			xp += dt2 * vxt;
			yp += dt2 * vyt;
			zp += dt2 * vzt;

			vxp += dt2 * ax;
			vyp += dt2 * ay;
			vzp += dt2 * az;

		}

		__syncthreads();

		if(itx == 0){

			dx_s[n-1] = 0.5 * (xt + xp);
			dy_s[n-1] = 0.5 * (yt + yp);
			dz_s[n-1] = 0.5 * (zt + zp);

			dvx_s[n-1] = 0.5 * (vxt + vxp);
			dvy_s[n-1] = 0.5 * (vyt + vyp);
			dvz_s[n-1] = 0.5 * (vzt + vzp);

			for(int j = n - 1; j >= 1; --j){
				double t0 = BSt0_c[(n-1) * 8 + (j-1)];
				double t1 = t0 * BSddt_c[j];
				double t2 = t0 * BSddt_c[n-1];

				dx_s[j-1] = (t1 * dx_s[j]) - (t2 * dx_s[j-1]);
				dy_s[j-1] = (t1 * dy_s[j]) - (t2 * dy_s[j-1]);
				dz_s[j-1] = (t1 * dz_s[j]) - (t2 * dz_s[j-1]);

				dvx_s[j-1] = (t1 * dvx_s[j]) - (t2 * dvx_s[j-1]);
				dvy_s[j-1] = (t1 * dvy_s[j]) - (t2 * dvy_s[j-1]);
				dvz_s[j-1] = (t1 * dvz_s[j]) - (t2 * dvz_s[j-1]);

			}

			double error1 = 0.0;

			error1 = dx_s[0] / scalex;
			error += error1 * error1;

			error1 = dy_s[0] / scaley;
			error += error1 * error1;

			error1 = dz_s[0] / scalez;
			error += error1 * error1;

			error1 = dvx_s[0] / scalevx;
			error += error1 * error1;

			error1 = dvy_s[0] / scalevy;
			error += error1 * error1;

			error1 = dvz_s[0] / scalevz;
			error += error1 * error1;


			error = sqrt(error / 6.0);      //6 is the number of dimensions
//printf("error %d %d %d %.20g\n", level, id, n, error);


			if(error < 1.0){
				//update
				xt = dx_s[0];
				yt = dy_s[0];
				zt = dz_s[0];

				vxt = dvx_s[0];
				vyt = dvy_s[0];
				vzt = dvz_s[0];

				for(int j = 1; j < n; ++j){
					xt += dx_s[j];
					yt += dy_s[j];
					zt += dz_s[j];

					vxt += dvx_s[j];
					vyt += dvy_s[j];
					vzt += dvz_s[j];
				}
				dx_d[id] = xt;
				dy_d[id] = yt;
				dz_d[id] = zt;

				dvx_d[id] = vxt;
				dvy_d[id] = vyt;
				dvz_d[id] = vzt;


				snew = 1.0;

				if(n >= 8){
					snew = 0.55;
				}
				if(n < 7){
					snew = 1.3;
				}
//printf("Accecpt %d %d %d %g %g\n", level, id, n, dt, snew);
//				break;

			}
			else{
				f_s[0] = 1;
				if(n == 8){
					snew = 0.5;
//printf("Repeat %d %d %d %g %g\n", level, id, n, dt, snew);
				}
			}


		}
		__syncthreads();

		if(f_s[0] == 0){
			break;
		}


	}//end of n loop

	if(itx == 0){
		if(snew * dt * dts >= dtlimit || level >= nL - 1){
			snew_d[id] = snew;
		}
		else{
			snew_d[id] = 1.0e6;      //mark body for higher level integration
		}
//printf("snew %d %d %g\n", level, id, snew_d[id]); 
	}

}


__global__ void computeError_d1_kernel(double *snew_d, double *ssum_d, int *index_d, const int N, const int N0, const int level){

	int jj = blockIdx.x * blockDim.x + threadIdx.x;
	int id = jj;		//particle id

	if(level > 0){
		if(jj < N){
			id = index_d[(level - 1) * N0 + jj];
		}
	}

	double s = 1.0e6;	//large number

	extern __shared__ double se_s[];
	double *s_s = se_s;

	int lane = threadIdx.x % warpSize;
	int warp = threadIdx.x / warpSize;

	if(warp == 0){
		s_s[threadIdx.x] = 1.0e6;	//large number
	}
	__syncthreads();

	if(jj < N){	
		s = snew_d[id];
//printf("snew %d %g\n", id, s);
	}	

	__syncthreads();

	double sr = s;

	for(int i = 1; i < warpSize; i*=2){
		sr = __shfl_xor_sync(0xffffffff, s, i, warpSize);
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
				sr = __shfl_xor_sync(0xffffffff, s, i, warpSize);
				s = fmin(s, sr);
//printf("r reduce 2  %d %d %d %.20g\n", blockIdx.x, i, threadIdx.x, s);
			}
		}
	}
	__syncthreads();
	if(threadIdx.x == 0){

		ssum_d[blockIdx.x] = s;
//printf("r reduce 3  %d %d %.20g\n", blockIdx.x, threadIdx.x, s);
	}

}


//This is the second part of the parallel reduction scheme
//It runs on only one thread block
__global__ void computeError_d2_kernel(double *ssum_d, const int N){

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
		s = ssum_d[idy];
	}
//printf("s reduce d2  %d %.20g\n", idy, s);

	__syncthreads();

	double sr = s;

	for(int i = 1; i < warpSize; i*=2){
		sr = __shfl_xor_sync(0xffffffff, s, i, warpSize);
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
				sr = __shfl_xor_sync(0xffffffff, s, i, warpSize);
				s = fmin(s, sr);
//printf("s reduce d2b %d %d %.20g\n", i, idy, s);
			}
		}
	}
	__syncthreads();
	if(threadIdx.x == 0){
		if(idx == 0){
//printf("s reduce d2c %d %.20g\n", idy, s);
			ssum_d[0] = s;
		}
	}
}


__global__ void update_kernel(double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *dx_d, double *dy_d, double *dz_d, double *dvx_d, double *dvy_d, double *dvz_d, double *snew_d, double *dtmin_d, int *index_d, int *Nlevel_d, const int N, const int N0, const int level, const int stopFlag, const double dt){

	int jj = blockIdx.x * blockDim.x + threadIdx.x;
	int id = jj;		//particle id


	if(jj < N){

		if(level > 0){
			id = index_d[(level - 1) * N0 + jj];
		}

		//accept step
		if(snew_d[id] < 1.0e6){
			x_d[id] += dx_d[id];
			y_d[id] += dy_d[id];
			z_d[id] += dz_d[id];

			vx_d[id] += dvx_d[id];
			vy_d[id] += dvy_d[id];
			vz_d[id] += dvz_d[id];
//printf("update %d\n", id);

			if(stopFlag == 0){
				if(dt < 0){
					dtmin_d[id] = dt > dtmin_d[id] ? dt : dtmin_d[id];
				}
				else{
					dtmin_d[id] = dt < dtmin_d[id] ? dt : dtmin_d[id];
				}
			}
//printf("dtmin %d %d %d %g\n", id, level, stopFlag, dtmin_d[id]);

		}
		else{
			int j = atomicAdd(&Nlevel_d[level + 1], 1);
			index_d[level * N0 + j] = id;
			//add to list
//printf("add to list %d %d\n", id, j);
		}
	}
}

__global__ void update_BS_kernel(double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *dx_d, double *dy_d, double *dz_d, double *dvx_d, double *dvy_d, double *dvz_d, double *snew_d, double *dtmin_d, int *index_d, int *Nlevel_d, const int N, const int N0, const int level, const int stopFlag, const double dt){

	int jj = blockIdx.x * blockDim.x + threadIdx.x;
	int id = jj;		//particle id


	if(jj < N){

		if(level > 0){
			id = index_d[(level - 1) * N0 + jj];
		}

		//accept step
		if(snew_d[id] < 1.0e6){
			x_d[id] = dx_d[id];
			y_d[id] = dy_d[id];
			z_d[id] = dz_d[id];

			vx_d[id] = dvx_d[id];
			vy_d[id] = dvy_d[id];
			vz_d[id] = dvz_d[id];
//printf("update %d\n", id);

			if(stopFlag == 0){
				if(dt < 0){
					dtmin_d[id] = dt > dtmin_d[id] ? dt : dtmin_d[id];
				}
				else{
					dtmin_d[id] = dt < dtmin_d[id] ? dt : dtmin_d[id];
				}
			}
//printf("dtmin %d %d %d %g\n", id, level, stopFlag, dtmin_d[id]);

		}
		else{
			int j = atomicAdd(&Nlevel_d[level + 1], 1);
			index_d[level * N0 + j] = id;
			//add to list
//printf("add to list %d %d\n", id, j);
		}
	}
}


__device__ void CartToKepGPU(double *xx_d, double *yy_d, double *zz_d, double *vxx_d, double *vyy_d, double *vzz_d, int i, double &a, double &e, double &inc, double &Omega, double &w, double &Theta, double &M, double &E, double Msun){

	double mu = Msun;

	double x = xx_d[i];
	double y = yy_d[i];
	double z = zz_d[i];
	double vx = vxx_d[i];
	double vy = vyy_d[i];
	double vz = vzz_d[i];


	double rsq = x * x + y * y + z * z;
	double vsq = vx * vx + vy * vy + vz * vz;
	double u =  x * vx + y * vy + z * vz;
	double ir = 1.0 / sqrt(rsq);
	double ia = 2.0 * ir - vsq / mu;

	a = 1.0 / ia;

	//inclination
	double hx, hy, hz;

	hx = ( y * vz) - (z * vy);
	hy = (-x * vz) + (z * vx);
	hz = ( x * vy) - (y * vx);

	double h2 = hx * hx + hy * hy + hz * hz;
	double h = sqrt(h2);

	double t = hz / h;
	if(t < -1.0) t = -1.0;
	if(t > 1.0) t = 1.0;

	inc = acos(t);

	//longitude of ascending node
	double n = sqrt(hx * hx + hy * hy);
	Omega = acos(-hy / n);
	if(hx < 0.0){
		Omega = 2.0 * M_PI - Omega;
	}

	if(inc < 1.0e-10 || n == 0) Omega = 0.0;

	//argument of periapsis
	double ex, ey, ez;

	ex = ( vy * hz - vz * hy) / mu - x * ir;
	ey = (-vx * hz + vz * hx) / mu - y * ir;
	ez = ( vx * hy - vy * hx) / mu - z * ir;


	e = sqrt(ex * ex + ey * ey + ez * ez);

	t = (-hy * ex + hx * ey) / (n * e);
	if(t < -1.0) t = -1.0;
	if(t > 1.0) t = 1.0;
	w = acos(t);
	if(ez < 0.0) w = 2.0 * M_PI - w;
	if(n == 0) w = 0.0;

	//True Anomaly
	t = (ex * x + ey * y + ez * z) / e * ir;
	if(t < -1.0) t = -1.0;
	if(t > 1.0) t = 1.0;
	Theta = acos(t);

	if(u < 0.0){
		if(e < 1.0 - 1.0e-10){
			//elliptic
			Theta = 2.0 * M_PI - Theta;
		}
		else if(e > 1.0 + 1.0e-10){
			//hyperbolic
			Theta = -Theta;
		}
		else{
			//parabolic
			Theta = - Theta;
		}
	}

	//Non circular, equatorial orbit
	if(e > 1.0e-10 && inc < 1.0e-10){
		Omega = 0.0;
		w = acos(ex / e);
		if(ey < 0.0) w = 2.0 * M_PI - w;
	}

	//circular, inclinded orbit
	if(e <= 1.0e-10 && inc > 1.0e-11){
		w = 0.0;
	}

	//circular, equatorial orbit
	if(e <= 1.0e-10 && inc <= 1.0e-11){
		w = 0.0;
		Omega = 0.0;
	}

	if(w == 0 && Omega != 0.0){
		t = (-hy * x + hx * y) / n * ir;
		if(t < -1.0) t = -1.0;
		if(t > 1.0) t = 1.0;
		Theta = acos(t);
		if(z < 0.0){
			if(e < 1.0 - 1.0e-10){
				//elliptic
				Theta = 2.0 * M_PI - Theta;
			}
			else if(e > 1.0 + 1.0e-10){
				//hyperbolic
				Theta = -Theta;
			}
			else{
				//parabolic
				Theta = -Theta;
			}
		}
	}
	if(w == 0 && Omega == 0.0){
		Theta = acos(x * ir);
		if(y < 0.0){
			if(e < 1.0 - 1.0e-10){
				//elliptic
				Theta = 2.0 * M_PI - Theta;
			}
			else if(e > 1.0 + 1.0e-10){
				//hyperbolic
				Theta = -Theta;
			}
			else{
				//parabolic
				Theta = -Theta;
			}
		}
	}

	if(e < 1.0 - 1.0e-10){
		//Eccentric Anomaly
		E = acos((e + cos(Theta)) / (1.0 + e * cos(Theta)));
		if(M_PI < Theta && Theta < 2.0 * M_PI) E = 2.0 * M_PI - E;

		//Mean Anomaly
		M = E - e * sin(E);
//printf("%g %g %g %g\n", Theta, E, M, w);
	}
	else if(e > 1.0 + 1.0e-10){
		//Hyperbolic Anomaly
		//named still E instead of H or F
		E = acosh((e + t) / (1.0 + e * t));
		if(Theta < 0.0) E = - E;

		M = e * sinh(E) - E;
	}
	else{
		//Parabolic Anomaly
		E = tan(Theta * 0.5);
		if(E > M_PI) E = E - 2.0 * M_PI;

		M = E + E * E * E / 3.0;

		//use a to store q
		a = h * h / mu * 0.5;
	}

}


__global__ void convertOutput_kernel(double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *xout_d, double *yout_d, double *zout_d, double *vxout_d, double *vyout_d, double *vzout_d, double *xTable_d, double *yTable_d, double *zTable_d, double *vxTable_d, double *vyTable_d, double *vzTable_d, const int nStage, const int N, const int Outecliptic, const int Outheliocentric, const int Outgeocentric, const int Outorbital, const double Obliquity, double Msun){

	int id = blockIdx.x * blockDim.x + threadIdx.x;

	if(id < N){
		xout_d[id] = x_d[id];
		yout_d[id] = y_d[id];
		zout_d[id] = z_d[id];

		vxout_d[id] = vx_d[id];
		vyout_d[id] = vy_d[id];
		vzout_d[id] = vz_d[id];
	}

	//If needed, convert from barycentric equatorial coordinates to barycentric ecliptic coordinates
	//if(Outecliptic == 1){
	//	EquatorialtoEcliptic(xout_h, yout_h, zout_h, vxout_h, vyout_h, vzout_h);
	//}

	if(Outheliocentric == 1){
		//Convert Barycentric coordinates to HelioCentric coordinates
		if(id < N){
			int ii = 10 * nStage;
			xout_d[id] -= xTable_d[ii];
			yout_d[id] -= yTable_d[ii];
			zout_d[id] -= zTable_d[ii];

			vxout_d[id] -= vxTable_d[ii];
			vyout_d[id] -= vyTable_d[ii];
			vzout_d[id] -= vzTable_d[ii];
		}
	}
	if(Outgeocentric == 1){
		//Convert Barycentric coordinates to geoCentric coordinates
		if(id < N){
			int ii = 2 * nStage;
			xout_d[id] -= xTable_d[ii];
			yout_d[id] -= yTable_d[ii];
			zout_d[id] -= zTable_d[ii];

			vxout_d[id] -= vxTable_d[ii];
			vyout_d[id] -= vyTable_d[ii];
			vzout_d[id] -= vzTable_d[ii];
//printf("Earth %.20g %.20g %.20g\n", xTable_d[ii], yTable_d[ii], zTable_d[ii]);
		}
	}
	if(Outecliptic == 1){
		if(id < N){
			double eps = Obliquity / 3600.0 / 180.0 * M_PI; // convert arcseconds to radians

			double ceps = cos(eps);
			double seps = sin(eps);


			double x = xout_d[id];
			double y = yout_d[id];
			double z = zout_d[id];
			double vx = vxout_d[id];
			double vy = vyout_d[id];
			double vz = vzout_d[id];


			xout_d[id] = x;
			yout_d[id] =  ceps * y + seps * z;
			zout_d[id] = -seps * y + ceps * z;

			vxout_d[id] = vx;
			vyout_d[id] =  ceps * vy + seps * vz;
			vzout_d[id] = -seps * vy + ceps * vz;
		}
	}

	if(Outorbital == 1){
		if(id < N){
			double a, e, inc, Omega, w, Theta, M, E;
			CartToKepGPU(xout_d, yout_d, zout_d, vxout_d, vyout_d, vzout_d, id, a, e, inc, Omega, w, Theta, M, E, Msun);


			inc = inc * 180.0 / M_PI;	//convert rad to deg
			Omega = Omega * 180.0 / M_PI;	//convert rad to deg
			w = w * 180.0 / M_PI;		//convert rad to deg
			M = M * 180.0 / M_PI;		//convert rad to deg

			xout_d[id] = a;
			yout_d[id] = e;
			zout_d[id] = inc;
			vxout_d[id] = Omega;
			vyout_d[id] = w;
			vzout_d[id] = M;


		}
	}


}


__global__ void set_dtmin_kernel(double *dtmin_d, double dtinit, int N){

	int id = blockIdx.x * blockDim.x + threadIdx.x;

	if(id < N){
		dtmin_d[id] = dtinit;
	}

}

int asteroid::loop_individual(){


	return 1;
}
	
int asteroid::loop(){

	//If needed, convert from heliocentric coordinates to barycentric coordinates
	if(ICheliocentric == 1){
		update_perturbers_kernel <<< nStage, 32 >>>(xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, data_d, cdata_d, idp_d, startTime_d, endTime_d, nChebyshev_d, offset0_d, timeStart, time_reference, dt, nStage, nCm, EM, AUtokm, Nperturbers);
		HelioToBary_kernel <<< (N + 255) / 256 , 256 >>> (xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, N, nStage);

	}

	//At this point, the initial conditions coordinates are cartesian barycentric equatorial
	//The integration is also done in cartesian barycentric equatorial coordinates

	if(outBinary == 0){
		outputFile = fopen(outputFilename, "w");
	}
	else{
		outputFile = fopen(outputFilename, "wb");
	}
	if(printdt == 1){
		dtFile = fopen(dtFilename, "w");
	}
	printf("Start integration %.20g\n", timeStart + time_reference);



	dt_h[0] = dt;           

	if(nL > 1 && dt_h[0] * dts < dtlimit[0]){
		dt_h[0] = dtlimit[0] * dts;
	}

	for(int i = 1; i < nL; ++i){
		dt_h[i] = dtlimit[i - 1] * dts;
	}                                       
	for(int i = 0; i < nL; ++i){    
		dtsave_h[i] = dt;       
		time_h[i] = time;               
		stop_h[i] = 0;          
		dtminlevel_h[i] = dt_h[i];
	}

	//reduce loop to numer of levels
	for(int i = 0; i < N; ++i){
		dtmin_h[i] = dt;
		timeStep_h[i] = 0ll;
	}

	set_dtmin_kernel <<< (N + 255) / 256 , 256 >>> (dtmin_d, dt, N);

	if(time_reference + time >= outStart){
		if(Outheliocentric == 1 || Outgeocentric == 1){
			update_perturbers_kernel <<< nStage, 32 >>>(xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, data_d, cdata_d, idp_d, startTime_d, endTime_d, nChebyshev_d, offset0_d, time, time_reference, dt, nStage, nCm, EM, AUtokm, Nperturbers);

		}
		convertOutput_kernel <<< (N + 255) / 256 , 256 >>> (x_d, y_d, z_d, vx_d, vy_d, vz_d, xout_d, yout_d, zout_d, vxout_d, vyout_d, vzout_d, xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, nStage, N, Outecliptic, Outheliocentric, Outgeocentric, Outorbital, Obliquity, GM_h[10]);
		copyOutput();
		printOutput();
	}

	for(int tt = 0; tt < MaxTimeSteps1; ++tt){

	//dtmin is the minimum time step of an output intervall, only used for diagnostics
	for(int i = 0; i < nL; ++i){
			dtminlevel_h[i] = 1000.0 * dts;
		}

		set_dtmin_kernel <<< (N + 255) / 256 , 256 >>> (dtmin_d, outputInterval * dts, N);


		//next output time
		double timett1 = timeStart + dts * (tt + 1) * outputInterval;

		double snew = 10.0;
//printf("integrate %.20g %.20g %.20g\n", timeStart + dts * tt * 10.0, timett1, dt);

		//integrate until the next output interval
		for(int ttt = 0; ttt < MaxTimeSteps2; ++ttt){

			//refine last time step of interval to match output time
			if(dts < 0){
				if(dt_h[0] < -outputInterval){
					dt_h[0] = -outputInterval;
				}

				//refine last time step of interval to match output time
				if(time + dt_h[0] < timett1){
printf("refine %d %.20g | %.20g %.20g %.20g\n", 0, dt_h[0], time + dt_h[0], timett1, timett1 - time);
					dtsave_h[0] = dt_h[0];
					dt_h[0] = timett1 - time;
					stop_h[0] = 1;
				}
				else{
					dtminlevel_h[0] = dt_h[0] > dtminlevel_h[0] ? dt_h[0] : dtminlevel_h[0];
				}
			}
			else{
				if(dt_h[0] > outputInterval){
					dt_h[0] = outputInterval;
				}

				//refine last time step of interval to match output time
				if(time + dt_h[0] > timett1){
printf("refine %d %.20g | %.20g %.20g %.20g\n", 0, dt_h[0], time + dt_h[0], timett1, timett1 - time);
					dtsave_h[0] = dt_h[0];
					dt_h[0] = timett1 - time;
					stop_h[0] = 1;
				}
				else{
					dtminlevel_h[0] = dt_h[0] < dtminlevel_h[0] ? dt_h[0] : dtminlevel_h[0];
				}

			}

			Nlevel_h[0] = N;
			for(int i = 1; i < nL; ++i){
				Nlevel_h[i] = 0;
			}

			//do a time step of length dt

			if(strcmp(integratorName, "LF") == 0){
				leapfrog_stepA_kernel <<< (N + 255) / 256 , 256 >>> (x_d, y_d, z_d, vx_d, vy_d, vz_d, N, dt);
				time += dt * 0.5;
				//Needs at least Nperturbers threads per block
				update_perturbers_kernel <<< nStage, 32 >>>(xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, data_d, cdata_d, idp_d, startTime_d, endTime_d, nChebyshev_d, offset0_d, time, time_reference, dt, nStage, nCm, EM, AUtokm, Nperturbers);

				//Needs at least Nperturbers threads per block
				if(GPUMode == 0){
					leapfrog_stepB_kernel <<< (N + 255) / 256 , 256 >>> (xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, A1_d, A2_d, A3_d, Tsave_d, Rsave_d, time, timeStep, GM_d, Nperturbers, N, dt);
				}
				else{
					leapfrog_stepB2_kernel <<< N, def_NP >>> (xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, A1_d, A2_d, A3_d, Tsave_d, Rsave_d, time, timeStep, GM_d, Nperturbers, N, dt);
				}
				time += dt * 0.5;
				++timeStep;
			}


			if(strcmp(integratorName, "RK4") == 0 || strcmp(integratorName, "RK7") == 0 ){
				//Needs at least Nperturbers threads per block
				update_perturbers_kernel <<< nStage, 32 >>>(xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, data_d, cdata_d, idp_d, startTime_d, endTime_d, nChebyshev_d, offset0_d, time, time_reference, dt, nStage, nCm, EM, AUtokm, Nperturbers);
	
				//Needs at least Nperturbers threads per block
				if(GPUMode == 0){
					RK_step_kernel <<< (N + 63) / 64 , dim3(64, 1, 1) >>> (xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, kx_d, ky_d, kz_d, kvx_d, kvy_d, kvz_d, GM_d, A1_d, A2_d, A3_d, Tsave_d, Rsave_d, time, timeStep, dt, RKFn, Nperturbers, N);
				}
				else{
					RK_step2_kernel <<< N, def_NP >>> (xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, GM_d, A1_d, A2_d, A3_d, Tsave_d, Rsave_d, time, timeStep, dt, RKFn, Nperturbers, N);
				}


				time += dt;
				++timeStep;
			}


			if(strcmp(integratorName, "RKF45") == 0 || strcmp(integratorName, "DP54") == 0 || strcmp(integratorName, "RKF78") == 0){
				int N0 = Nlevel_h[0];
				int level = 0;
				
				//Needs at least Nperturbers threads per block
				update_perturbers_kernel <<< nStage, 32 >>>(xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, data_d, cdata_d, idp_d, startTime_d, endTime_d, nChebyshev_d, offset0_d, time, time_reference, dt_h[level], nStage, nCm, EM, AUtokm, Nperturbers);


		
				int stopFlag = 0;
				for(int l = 0; l <= level; ++l){
					if(stop_h[l] == 1){
						stopFlag = 1;
					}
				}

				double dtlimit_ = dtlimit[level];
				if(stopFlag == 1){
					dtlimit_ = dt_h[level] * dts;
				}

				if(GPUMode == 0){
					RKF_step_kernel <<< (Nlevel_h[level] + 63) / 64 , 64 >>> (xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, dx_d, dy_d, dz_d, dvx_d, dvy_d, dvz_d, kx_d, ky_d, kz_d, kvx_d, kvy_d, kvz_d, GM_d, A1_d, A2_d, A3_d, index_d, Tsave_d, Rsave_d, snew_d, time, timeStep, dt_h[level], dts, RKFn, Nperturbers, Nlevel_h[level], N0, level, nL, dtlimit_);
				}
				else{
					RKF_step2_kernel <<< Nlevel_h[level], def_NP >>> (xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, dx_d, dy_d, dz_d, dvx_d, dvy_d, dvz_d, GM_d, A1_d, A2_d, A3_d, index_d, Tsave_d, Rsave_d, snew_d, time, timeStep, dt_h[level], dts, RKFn, Nperturbers, Nlevel_h[level], N0, level, nL, dtlimit_);
				}


				//Calculate the minimal time step value
				//Using a parallel reduction sum
				int nct = 512;
				int ncb = min((Nlevel_h[level] + nct - 1) / nct, 1024);
				computeError_d1_kernel <<< ncb, nct, WarpSize * sizeof(double)  >>> (snew_d, ssum_d, index_d, Nlevel_h[level], N0, level);
				if(ncb > 1){
					computeError_d2_kernel <<< 1, ((ncb + WarpSize - 1) / WarpSize) * WarpSize, WarpSize * sizeof(double)  >>> (ssum_d, ncb);
				}

				cudaDeviceSynchronize();
				cudaMemcpy(&snew, ssum_d, sizeof(double), cudaMemcpyDeviceToHost);

				if(snew == 1.0e6){
					snew = 1.0;
				}

printf("snewMin %g %g\n", snew, dtlimit[level]);
				snewlevel_h[level] = snew;

				if(snew >= RKF_fac){
					//accept step
					update_kernel <<< (Nlevel_h[level] + 63) / 64 , 64 >>> (x_d, y_d, z_d, vx_d, vy_d, vz_d, dx_d, dy_d, dz_d, dvx_d, dvy_d, dvz_d, snew_d, dtmin_d, index_d, Nlevel_d, Nlevel_h[level], N0, level, stopFlag, dt_h[level]);

					cudaDeviceSynchronize();
					cudaMemcpy(&Nlevel_h[1], &Nlevel_d[1], sizeof(int), cudaMemcpyDeviceToHost);
printf("Nlevel %d %d\n", level + 1, Nlevel_h[1]);

					if(Nlevel_h[1] > 0){ 
						time_h[level + 1] = time_h[level];
printf("set time %d %d %g %g\n", level, level + 1, time_h[level], time_h[level + 1]);						
					}

					time_h[level] += dt_h[level];
					++timeStep_h[level];

					if(Nlevel_h[1] > 0){ 
						loop_recursive(1);
					}
				}
				else{
printf("repeat level %d %g\n", level, snew);
				}


				timeStep = timeStep_h[level];
				time = time_h[level];

			}


			if(strcmp(integratorName, "BS") == 0){
				int N0 = Nlevel_h[0];
				int level = 0;

				//Needs at least Nperturbers threads per block
				update_perturbers_kernel <<< nStage, 32 >>>(xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, data_d, cdata_d, idp_d, startTime_d, endTime_d, nChebyshev_d, offset0_d, time, time_reference, dt_h[level], nStage, nCm, EM, AUtokm, Nperturbers);

				int stopFlag = 0;
				for(int l = 0; l <= level; ++l){
					if(stop_h[l] == 1){
						stopFlag = 1;
					}
				}

				double dtlimit_ = dtlimit[level];
				if(stopFlag == 1){
					dtlimit_ = dt_h[level] * dts;
				}


				if(GPUMode == 0){
					BS_step_kernel <<< (Nlevel_h[level] + def_N - 1) / def_N , def_N >>> (xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, dx_d, dy_d, dz_d, dvx_d, dvy_d, dvz_d, kx_d, ky_d, kz_d, kvx_d, kvy_d, kvz_d, GM_d, A1_d, A2_d, A3_d, index_d, Tsave_d, Rsave_d, snew_d, time, timeStep, dt_h[level], dts, BSn, Nperturbers, Nlevel_h[level], N0, level, nL, dtlimit_);
				}
				else{
					BS_step2_kernel <<< Nlevel_h[level], def_NP >>> (xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, dx_d, dy_d, dz_d, dvx_d, dvy_d, dvz_d, kx_d, ky_d, kz_d, kvx_d, kvy_d, kvz_d, GM_d, A1_d, A2_d, A3_d, index_d, Tsave_d, Rsave_d, snew_d, time, timeStep, dt_h[level], dts, BSn, Nperturbers, Nlevel_h[level], N0, level, nL, dtlimit_);
				}


				//Calculate the minimal time step value
				//Using a parallel reduction sum
				int nct = 512;
				int ncb = min((Nlevel_h[level] + nct - 1) / nct, 1024);
				computeError_d1_kernel <<< ncb, nct, WarpSize * sizeof(double)  >>> (snew_d, ssum_d, index_d, Nlevel_h[level], N0, level);
				if(ncb > 1){
					computeError_d2_kernel <<< 1, ((ncb + WarpSize - 1) / WarpSize) * WarpSize, WarpSize * sizeof(double)  >>> (ssum_d, ncb);
				}

				cudaDeviceSynchronize();
				cudaMemcpy(&snew, ssum_d, sizeof(double), cudaMemcpyDeviceToHost);

				if(snew == 1.0e6){
					snew = 1.0;
				}

printf("snewMin %g %g\n", snew, dtlimit[level]);

				snewlevel_h[level] = snew;

				if(snew > 0.5){
					//accept step
					update_BS_kernel <<< (Nlevel_h[level] + 63) / 64 , 64 >>> (x_d, y_d, z_d, vx_d, vy_d, vz_d, dx_d, dy_d, dz_d, dvx_d, dvy_d, dvz_d, snew_d, dtmin_d, index_d, Nlevel_d, Nlevel_h[level], N0, level, stopFlag, dt_h[level]);

					cudaDeviceSynchronize();
					cudaMemcpy(&Nlevel_h[1], &Nlevel_d[1], sizeof(int), cudaMemcpyDeviceToHost);
printf("Nlevel %d %d\n", level + 1, Nlevel_h[1]);

					if(Nlevel_h[1] > 0){ 
						time_h[level + 1] = time_h[level];
printf("set time %d %d %g %g\n", level, level + 1, time_h[level], time_h[level + 1]);						
					}

					time_h[level] += dt_h[level];
					++timeStep_h[level];

					if(Nlevel_h[1] > 0){ 
						loop_recursive(1);
					}

				}
				else{
printf("repeat level %d %g\n", level, snew);
				}

				timeStep = timeStep_h[level];
				time = time_h[level];
			}
			if(strcmp(integratorName, "IMM") == 0 ){
				//Needs at least Nperturbers threads per block

				update_perturbers_kernel <<< nStage, 32 >>>(xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, data_d, cdata_d, idp_d, startTime_d, endTime_d, nChebyshev_d, offset0_d, time + dt * 0.5, time_reference, dt, nStage, nCm, EM, AUtokm, Nperturbers);
	
				//Needs at least Nperturbers threads per block
				if(GPUMode == 0){
					IMM_step_kernel <<< (N + 63) / 64 , 64 >>> (xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, A1_d, A2_d, A3_d, Tsave_d, Rsave_d, snew_d, time, timeStep, GM_d, Nperturbers, N, dt);
				}
				else{
					IMM_step2_kernel <<< N, def_NP >>> (xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, A1_d, A2_d, A3_d, Tsave_d, Rsave_d, snew_d, time, timeStep, GM_d, Nperturbers, N, dt);

				}

				cudaDeviceSynchronize();
				cudaMemcpy(snew_h, snew_d, sizeof(double), cudaMemcpyDeviceToHost);
				if(snew_h[0] == -1){
					printf("Error, implicit midpoint method did not converge\n");
					return 0;
				}


				time += dt;
				++timeStep;
			}


			if(printdt == 1){
				fprintf(dtFile, "%-25.20g %lld %-25.20g\n", time + time_reference, timeStep, dt_h[0] / snew);
			}
printf("dt %d %lld %g %g %g\n", 0, timeStep, snew, dt_h[0], dtsave_h[0]);



			dt_h[0] *= snew;


			if(dts < 0 && time <= timett1){
				//set time step equal to the last accepted full time step

				if(snew >= 1.0 && stop_h[0] == 1){

					for(int l = 0; l < nL; ++l){
						dt_h[l] = dtsave_h[l];
printf("reset %d  %.20g\n", l, dt_h[l]);

					}
				}

				stop_h[0] = 0;
				break;
			}
			if(dts > 0 && time >= timett1){
				//set time step equal to the last accepted full time step


				if(snew >= 1.0 && stop_h[0] == 1){

					for(int l = 0; l < nL; ++l){
						dt_h[l] = dtsave_h[l];
printf("reset %d  %.20g\n", l, dt_h[l]);

					}

				}
				stop_h[0] = 0;
				break;
			}


			if(time + time_reference > time1 || time + time_reference < time0){
				cudaDeviceSynchronize();
				printf("Reached the end of the Chebyshev data file\n");
				return 0;
			}



			if(ttt >= MaxTimeSteps2 - 1){

				printf("Error, time step loop2 did not finish\n");
				return 0;
			}

		}//end of ttt loop
		cudaDeviceSynchronize();
		if(time_reference + time >= outStart){
			if(Outheliocentric == 1 || Outgeocentric == 1){
				update_perturbers_kernel <<< nStage, 32 >>>(xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, data_d, cdata_d, idp_d, startTime_d, endTime_d, nChebyshev_d, offset0_d, time, time_reference, dt, nStage, nCm, EM, AUtokm, Nperturbers);
			}
			convertOutput_kernel <<< (N + 255) / 256 , 256 >>> (x_d, y_d, z_d, vx_d, vy_d, vz_d, xout_d, yout_d, zout_d, vxout_d, vyout_d, vzout_d, xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, nStage, N, Outecliptic, Outheliocentric, Outgeocentric, Outorbital, Obliquity, GM_h[10]);
			copyOutput();
			printOutput();
			fflush(outputFile);
		}

		if(dts < 0 && time < timeEnd){
			cudaDeviceSynchronize();
			printf("Reached the end of the integration\n");
			return 0;
		}
		if(dts > 0 && time > timeEnd){
			cudaDeviceSynchronize();
			printf("Reached the end of the integration\n");
			return 0;
		}

		if(tt >= MaxTimeSteps1 - 1){

			printf("Error, time step loop2 did not finish\n");
			return 0;
		}

	}//end of tt loop
	return 1;
}

int asteroid::loop_recursive(int level){
	double snew = 10.0;	
printf("\nStart level %d %g %g\n", level, dt_h[level - 1], dt_h[level]);
	for(int t = 0; t < MaxTimeSteps2; ++t){


		if(dts < 0){
			if(dt_h[level] < dt_h[level - 1]){
				dt_h[level] = dt_h[level - 1];
			}


			//refine last time step of interval to match output time
			if(time_h[level] + dt_h[level] < time_h[level - 1]){
printf("refine %d %.20g | %.20g %.20g %.20g\n", level, dt_h[level - 1], time_h[level] + dt_h[level], time_h[level - 1], time_h[level - 1] - time_h[level]);
				dtsave_h[level] = dt_h[level];
				dt_h[level] = time_h[level - 1] - time_h[level];
				stop_h[level] = 1;
			}
			else{
				dtminlevel_h[level] = dt_h[level] > dtminlevel_h[level] ? dt_h[level] : dtminlevel_h[level];
			}
		}
		else{
			if(dt_h[level] > dt_h[level - 1]){
				dt_h[level] = dt_h[level - 1];
			}

			//refine last time step of interval to match output time
			if(time_h[level] + dt_h[level] > time_h[level - 1]){
printf("refine %d %.20g | %.20g %.20g %.20g\n", 1, dt_h[level - 1], time_h[level] + dt_h[level], time_h[level - 1], time_h[level - 1] - time_h[level]);
				dtsave_h[level] = dt_h[level];
				dt_h[level] = time_h[level - 1] - time_h[level];
				stop_h[level] = 1;
			}
			else{
				dtminlevel_h[level] = dt_h[level] < dtminlevel_h[level] ? dt_h[level] : dtminlevel_h[level];
			}

		}

		if(strcmp(integratorName, "RKF45") == 0 || strcmp(integratorName, "DP54") == 0 || strcmp(integratorName, "RKF78") == 0){
			int N0 = Nlevel_h[0];

			//Needs at least Nperturbers threads per block
			update_perturbers_kernel <<< nStage, 32 >>>(xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, data_d, cdata_d, idp_d, startTime_d, endTime_d, nChebyshev_d, offset0_d, time_h[level], time_reference, dt_h[level], nStage, nCm, EM, AUtokm, Nperturbers);

			int stopFlag = 0;
			for(int l = 0; l <= level; ++l){
				if(stop_h[l] == 1){
					stopFlag = 1;
				}
			}

			double dtlimit_ = dtlimit[level];
			if(stopFlag == 1){
				dtlimit_ = dt_h[level] * dts;
			}

			if(GPUMode == 0){
				RKF_step_kernel <<< (Nlevel_h[level] + 63) / 64 , 64 >>> (xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, dx_d, dy_d, dz_d, dvx_d, dvy_d, dvz_d, kx_d, ky_d, kz_d, kvx_d, kvy_d, kvz_d, GM_d, A1_d, A2_d, A3_d, index_d, Tsave_d, Rsave_d, snew_d, time_h[level], timeStep_h[level], dt_h[level], dts, RKFn, Nperturbers, Nlevel_h[level], N0, level, nL, dtlimit_);
			}
			else{
				RKF_step2_kernel <<< Nlevel_h[level], def_NP >>> (xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, dx_d, dy_d, dz_d, dvx_d, dvy_d, dvz_d, GM_d, A1_d, A2_d, A3_d, index_d, Tsave_d, Rsave_d, snew_d, time_h[level], timeStep_h[level], dt_h[level], dts, RKFn, Nperturbers, Nlevel_h[level], N0, level, nL, dtlimit_);
			}


			//Calculate the minimal time step value
			//Using a parallel reduction sum
			int nct = 512;
			int ncb = min((Nlevel_h[level] + nct - 1) / nct, 1024);
			computeError_d1_kernel <<< ncb, nct, WarpSize * sizeof(double)  >>> (snew_d, ssum_d, index_d, Nlevel_h[level], N0, level);
			if(ncb > 1){
				computeError_d2_kernel <<< 1, ((ncb + WarpSize - 1) / WarpSize) * WarpSize, WarpSize * sizeof(double)  >>> (ssum_d, ncb);
			}

			cudaDeviceSynchronize();
			cudaMemcpy(&snew, ssum_d, sizeof(double), cudaMemcpyDeviceToHost);

			if(snew == 1.0e6){
				snew = 1.0;
			}

printf("snewMin %g %g\n", snew, dtlimit[level]);
			snewlevel_h[level] = snew;

			if(snew >= RKF_fac){
				//accept step
				update_kernel <<< (Nlevel_h[level] + 63) / 64 , 64 >>> (x_d, y_d, z_d, vx_d, vy_d, vz_d, dx_d, dy_d, dz_d, dvx_d, dvy_d, dvz_d, snew_d, dtmin_d, index_d, Nlevel_d, Nlevel_h[level], N0, level, stopFlag, dt_h[level]);

				cudaDeviceSynchronize();
				cudaMemcpy(&Nlevel_h[level + 1], &Nlevel_d[level + 1], sizeof(int), cudaMemcpyDeviceToHost);

				if(Nlevel_h[level + 1] > 0){
					time_h[level + 1] = time_h[level];


				}

				time_h[level] += dt_h[level];
				++timeStep_h[level];

				if(Nlevel_h[level + 1] > 0){

					loop_recursive(level + 1);

				}


			}
			else{
printf("repeat level %d %g\n", level, snew);
			}


		}

		if(strcmp(integratorName, "BS") == 0){
			int N0 = Nlevel_h[0];

			//Needs at least Nperturbers threads per block
			update_perturbers_kernel <<< nStage, 32 >>>(xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, data_d, cdata_d, idp_d, startTime_d, endTime_d, nChebyshev_d, offset0_d, time_h[level], time_reference, dt_h[level], nStage, nCm, EM, AUtokm, Nperturbers);

			int stopFlag = 0;
			for(int l = 0; l <= level; ++l){
				if(stop_h[l] == 1){
					stopFlag = 1;
				}
			}

			double dtlimit_ = dtlimit[level];
			if(stopFlag == 1){
				dtlimit_ = dt_h[level] * dts;
			}

			if(GPUMode == 0){
				BS_step_kernel <<< (Nlevel_h[level] + def_N - 1) / def_N , def_N >>> (xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, dx_d, dy_d, dz_d, dvx_d, dvy_d, dvz_d, kx_d, ky_d, kz_d, kvx_d, kvy_d, kvz_d, GM_d, A1_d, A2_d, A3_d, index_d, Tsave_d, Rsave_d, snew_d, time, timeStep, dt_h[level], dts, BSn, Nperturbers, Nlevel_h[level], N0, level, nL, dtlimit_);
			}
			else{
				BS_step2_kernel <<< Nlevel_h[level], def_NP >>> (xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, dx_d, dy_d, dz_d, dvx_d, dvy_d, dvz_d, kx_d, ky_d, kz_d, kvx_d, kvy_d, kvz_d, GM_d, A1_d, A2_d, A3_d, index_d, Tsave_d, Rsave_d, snew_d, time, timeStep, dt_h[level], dts, BSn, Nperturbers, Nlevel_h[level], N0, level, nL, dtlimit_);
			}


			//Calculate the minimal time step value
			//Using a parallel reduction sum
			int nct = 512;
			int ncb = min((Nlevel_h[level] + nct - 1) / nct, 1024);
			computeError_d1_kernel <<< ncb, nct, WarpSize * sizeof(double)  >>> (snew_d, ssum_d, index_d, Nlevel_h[level], N0, level);
			if(ncb > 1){
				computeError_d2_kernel <<< 1, ((ncb + WarpSize - 1) / WarpSize) * WarpSize, WarpSize * sizeof(double)  >>> (ssum_d, ncb);
			}

			cudaDeviceSynchronize();
			cudaMemcpy(&snew, ssum_d, sizeof(double), cudaMemcpyDeviceToHost);

			if(snew == 1.0e6){
				snew = 1.0;
			}

printf("snewMin %g %g\n", snew, dtlimit[level]);
			snewlevel_h[level] = snew;

			if(snew > 0.5){
				//accept step
				update_BS_kernel <<< (Nlevel_h[level] + 63) / 64 , 64 >>> (x_d, y_d, z_d, vx_d, vy_d, vz_d, dx_d, dy_d, dz_d, dvx_d, dvy_d, dvz_d, snew_d, dtmin_d, index_d, Nlevel_d, Nlevel_h[level], N0, level, stopFlag, dt_h[level]);

				cudaDeviceSynchronize();
				cudaMemcpy(&Nlevel_h[level + 1], &Nlevel_d[level + 1], sizeof(int), cudaMemcpyDeviceToHost);

				if(Nlevel_h[level + 1] > 0){
					time_h[level + 1] = time_h[level];


				}

				time_h[level] += dt_h[level];
				++timeStep_h[level];

				if(Nlevel_h[level + 1] > 0){

					loop_recursive(level + 1);

				}


			}
			else{
printf("repeat level %d %g\n", level, snew);
			}

		}

		if(Nlevel_h[level + 1] > 0){

			loop_recursive(level + 1);

		}


		if(printdt == 1){
			fprintf(dtFile, "%-25.20g %lld %-25.20g %d\n", time_h[level] + time_reference, timeStep, dt_h[level] / snewlevel_h[level], level);
		}
printf("dt %d %lld %g %g %g\n", level, timeStep, snewlevel_h[level], dt_h[level], dtsave_h[level]);

		dt_h[level] *= snewlevel_h[level];



		if(dts < 0 && time_h[level] <= time_h[level - 1]){
			//set time step equal to the last accepted full time step

			if(snewlevel_h[level] >= 1.0 && stop_h[level] == 1){

				for(int l = level; l < nL; ++l){
					dt_h[l] = dtsave_h[l];
printf("reset %d  %.20g\n", l, dt_h[l]);

				}
			}

			stop_h[level] = 0;
			break;
		}
		if(dts > 0 && time_h[level] >= time_h[level - 1]){
			//set time step equal to the last accepted full time step


			if(snewlevel_h[level] >= 1.0 && stop_h[level] == 1){

				for(int l = level; l < nL; ++l){
					dt_h[l] = dtsave_h[l];
printf("reset %d  %.20g\n", l, dt_h[l]);

				}
			}

			stop_h[level] = 0;
			break;
		}
		if(dts > 0 && time_h[level] >= time_h[level - 1]){
			//set time step equal to the last accepted full time step


			if(snewlevel_h[level] >= 1.0 && stop_h[level] == 1){

				for(int l = level; l < nL; ++l){
					dt_h[l] = dtsave_h[l];
printf("reset %d  %.20g\n", l, dt_h[l]);

				}


			}
			stop_h[level] = 0;
			break;
		}
		if(t >= MaxTimeSteps2 - 1){

			printf("Error, time step loop2 in level %d did not finish\n", level);
			return 0;
		}


	} // end of t loop
	Nlevel_h[level] = 0;
	cudaMemset(&Nlevel_d[level], 0, sizeof(int));

	return 1;
}


#include <stdio.h>

//Define Constant Memory
__constant__ double RKFa_c[20 * 20];       //20 is considered here to be large enough (>RKFn)
__constant__ double RKFb_c[20];
__constant__ double RKFbb_c[20];
__constant__ double RKFc_c[20];


__constant__ double RKF_atol_c;
__constant__ double RKF_rtol_c;
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
	cudaMemcpyToSymbol(RKFa_c, RKFa_h, RKFn * RKFn * sizeof(double), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(RKFb_c, RKFb_h, RKFn * sizeof(double), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(RKFbb_c, RKFbb_h, RKFn * sizeof(double), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(RKFc_c, RKFc_h, RKFn * sizeof(double), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(RKF_atol_c, &RKF_atol, sizeof(double), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(RKF_rtol_c, &RKF_rtol, sizeof(double), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(RKF_ee_c, &RKF_ee, sizeof(double), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(RKF_fac_c, &RKF_fac, sizeof(double), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(RKF_facmin_c, &RKF_facmin, sizeof(double), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(RKF_facmax_c, &RKF_facmax, sizeof(double), 0, cudaMemcpyHostToDevice);
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




__global__ void HelioToBary_kernel(double *xTable_d, double *yTable_d, double *zTable_d, double *vxTable_d, double *vyTable_d, double *vzTable_d, double *xx_d, double *yy_d, double *zz_d, double *vxx_d, double *vyy_d, double *vzz_d, const int N, const int RKFn){

	int id = blockIdx.x * blockDim.x + threadIdx.x;

	if(id < N){
		int ii = 10 * RKFn;

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
__global__ void leapfrog_stepB_kernel(double *xTable_d, double *yTable_d, double *zTable_d, double *vxTable_d, double *vyTable_d, double *vzTable_d, double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *A1_d, double *A2_d, double *A3_d, double *Tsave_d, double *Rsave_d, double time, long long timeStep, double *GM_d, const int Nperturbers, const int N, const int RKFn, const double dt){

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
		int ii = itx * RKFn;
		xTable_s[itx] = xTable_d[ii];
		yTable_s[itx] = yTable_d[ii];
		zTable_s[itx] = zTable_d[ii];

		vxTable_s[itx] = vxTable_d[ii];
		vyTable_s[itx] = vyTable_d[ii];
		vzTable_s[itx] = vzTable_d[ii];

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

// Runge Kutta step with fixed time step
// Every body runs on a sepparate thread block. Gravity calculation is spread along threads in the block and reuced in registers
// Kernel uses at least Nperturbers threads. The perturbers are loaded into shared memory
__global__ void RK_stage_kernel(double *xTable_d, double *yTable_d, double *zTable_d, double *vxTable_d, double *vyTable_d, double *vzTable_d, double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *ax_d, double *ay_d, double *az_d, double *kx_d, double *ky_d, double *kz_d, double *kvx_d, double *kvy_d, double *kvz_d, double *GM_d, double *A1_d, double *A2_d, double *A3_d, double *Tsave_d, double *Rsave_d, double time, long long int timeStep, const double dt, const int RKFn, const int Nperturbers, const int N, const int S){

	int id = blockIdx.x;	//particle index

//if(id < Nperturbers){
//printf("p %d %.20g %.20g %.20g %.20g %.20g %.20g %.20g\n", id, time + RKFc_c[S] * dt, xt_s[id], yt_s[id], zt_s[id], vxt_s[id], vyt_s[id], vzt_s[id]);
//}

	if(id < N){

		double xti = x_d[id];
		double yti = y_d[id];
		double zti = z_d[id];

		double vxti = vx_d[id];
		double vyti = vy_d[id];
		double vzti = vz_d[id];

		double ax = 0.0;
		double ay = 0.0;
		double az = 0.0;

		for(int s = 0; s < S; ++s){
			double dtaa = dt * RKFa_c[S * RKFn + s];
			int ii = id + s * N;
			xti  += dtaa * kx_d[ii];
			yti  += dtaa * ky_d[ii];
			zti  += dtaa * kz_d[ii];
			vxti += dtaa * kvx_d[ii];
			vyti += dtaa * kvy_d[ii];
			vzti += dtaa * kvz_d[ii];
//printf("update 2 %d %d %g %g %g %g %g %g\n", S, id, xti, yti, zti, RKFa_c[S * RKFn + s], kx_d[s], dt);

		}

		// ----------------------------------------------------------------------------
		//compute forces

		//heliocentric coordinates
		int iS = 10 * RKFn + S;
		double xih = xti - xTable_d[iS];
		double yih = yti - yTable_d[iS];
		double zih = zti - zTable_d[iS];

		double vxih = vxti - vxTable_d[iS];
		double vyih = vyti - vyTable_d[iS];
		double vzih = vzti - vzTable_d[iS];

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
		int iE = 2 * RKFn + S;
		double xiE = xti - xTable_d[iE];
		double yiE = yti - yTable_d[iE];
		double ziE = zti - zTable_d[iE];

		if(useNonGrav_c == 1){
			NonGrav(xih, yih, zih, vxih, vyih, vzih, A1_d[id], A2_d[id], A3_d[id], r, ax, ay, az);
		}
		if(useGR_c == 1){
			GR(xih, yih, zih, vxih, vyih, vzih, r, ax, ay, az, GM_d[10]);
		}
		if(useJ2_c == 1){
			J2(xiE, yiE, ziE, ax, ay, az, GM_d[2]);
		}
		ax_d[id] = ax;
		ay_d[id] = ay;
		az_d[id] = az;

	}

}


//Runge Kutta Fehlberg step with adaptive time step
// Every body runs on a thread
// Kernel uses at least Nperturbers threads. The perturbers are loaded into shared memory
__global__ void RKF_step_kernel(double *xTable_d, double *yTable_d, double *zTable_d, double *vxTable_d, double *vyTable_d, double *vzTable_d, double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *dx_d, double *dy_d, double *dz_d, double *dvx_d, double *dvy_d, double *dvz_d, double *kx_d, double *ky_d, double *kz_d, double *kvx_d, double *kvy_d, double *kvz_d, double *GM_d, double *A1_d, double *A2_d, double *A3_d, double *Tsave_d, double *Rsave_d, double *snew_d, double time, long long int timeStep, const double dt, const int dts, const int RKFn, const int Nperturbers, const int N, const int stop){


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
		//Update the Chebyshev coefficients if necessary
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
//printf("update 2 %d %d %g %g %g %g %g %g\n", S, id, xti, yti, zti, RKFa_c[S * RKFn + s], kx_d[s], dt);

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

	double dx = 0.0;
	double dy = 0.0;
	double dz = 0.0;

	double dvx = 0.0;
	double dvy = 0.0;
	double dvz = 0.0;

	double snew = 10.0;

	if(id < N){

		for(int S = 0; S < RKFn; ++S){
			double dtb = dt * RKFb_c[S];
			int ii = id + S * N;
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

		double ym = 0.0;
		ym = (fabs(x0) > ym) ? fabs(x0) : ym;
		ym = (fabs(y0) > ym) ? fabs(y0) : ym;
		ym = (fabs(z0) > ym) ? fabs(z0) : ym;

		ym = (fabs(vx0) > ym) ? fabs(vx0) : ym;
		ym = (fabs(vy0) > ym) ? fabs(vy0) : ym;
		ym = (fabs(vz0) > ym) ? fabs(vz0) : ym;

		double isc = 1.0 / (RKF_atol_c + ym * RKF_rtol_c);
		isc *= isc;

		//error estimation
		double errorkx = 0.0;
		double errorky = 0.0;
		double errorkz = 0.0;
		
		double errorkvx = 0.0;
		double errorkvy = 0.0;
		double errorkvz = 0.0;

		for(int S = 0; S < RKFn; ++S){
			double f = (RKFb_c[S] - RKFbb_c[S]) * dt;
			int ii = id + S * N;

			errorkx += __dmul_rn(f, kx_d[ii]);
			errorky += __dmul_rn(f, ky_d[ii]);
			errorkz += __dmul_rn(f, kz_d[ii]);

			errorkvx += __dmul_rn(f, kvx_d[ii]);
			errorkvy += __dmul_rn(f, kvy_d[ii]);
			errorkvz += __dmul_rn(f, kvz_d[ii]);

//printf("error %d %d %.20g %.20g %.20g %.20g %.20g %.20g %.20g\n", id, S, f, errorkx, errorky, errorkz, errorkvx, errorkvy, errorkvz);
		}

		double errork = 0.0;

		errork += __dmul_rn(errorkx, errorkx) * isc;
		errork += __dmul_rn(errorky, errorky) * isc;
		errork += __dmul_rn(errorkz, errorkz) * isc;
		errork += __dmul_rn(errorkvx, errorkvx) * isc;
		errork += __dmul_rn(errorkvy, errorkvy) * isc;
		errork += __dmul_rn(errorkvz, errorkvz) * isc;

		errork = sqrt(errork / 6.0);    //6 is the number of dimensions

		double s = pow( 1.0  / errork, RKF_ee_c);
//printf("%.20g %.20g\n", errork, s);

		s = (RKF_fac_c * s > RKF_facmin_c) ? RKF_fac_c * s : RKF_facmin_c;
		s = (RKF_facmax_c < s) ? RKF_facmax_c : s;

		snew = (snew < s) ? snew : s;

		snew_d[id] = snew;
//printf("id %d %g %g\n", id, s, snew);
	}
/*
	__syncthreads();

	if(snew >= 1.0){
		//accept step
		if(id < N){
			x_d[id] += dx;
			y_d[id] += dy;
			z_d[id] += dz;

			vx_d[id] += dvx;
			vy_d[id] += dvy;
			vz_d[id] += dvz;
		}
	}
*/
}

// Runge Kutta Fehlberg step with adaptive time step
// Every body runs on a thread
// Kernel uses at least Nperturbers threads. The perturbers are loaded into shared memory
__global__ void RKF_step2_kernel(double *xTable_d, double *yTable_d, double *zTable_d, double *vxTable_d, double *vyTable_d, double *vzTable_d, double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *dx_d, double *dy_d, double *dz_d, double *dvx_d, double *dvy_d, double *dvz_d,double *GM_d, double *A1_d, double *A2_d, double *A3_d, double *Tsave_d, double *Rsave_d, double *snew_d, double time, long long int timeStep, const double dt, const int dts, const int RKFn, const int Nperturbers, const int N, const int stop){


	int itx = threadIdx.x;
	int id = blockIdx.x; 	//particle id

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
		//Update the Chebyshev coefficients if necessary
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

	double snew = 10.0;

	if(id < N && itx == 0){

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

		double ym = 0.0;
		ym = (fabs(x0) > ym) ? fabs(x0) : ym;
		ym = (fabs(y0) > ym) ? fabs(y0) : ym;
		ym = (fabs(z0) > ym) ? fabs(z0) : ym;

		ym = (fabs(vx0) > ym) ? fabs(vx0) : ym;
		ym = (fabs(vy0) > ym) ? fabs(vy0) : ym;
		ym = (fabs(vz0) > ym) ? fabs(vz0) : ym;

		double isc = 1.0 / (RKF_atol_c + ym * RKF_rtol_c);
		isc *= isc;

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
		errork += __dmul_rn(errorkx, errorkx) * isc;
		errork += __dmul_rn(errorky, errorky) * isc;
		errork += __dmul_rn(errorkz, errorkz) * isc;
		errork += __dmul_rn(errorkvx, errorkvx) * isc;
		errork += __dmul_rn(errorkvy, errorkvy) * isc;
		errork += __dmul_rn(errorkvz, errorkvz) * isc;

		errork = sqrt(errork / 6.0);    //6 is the number of dimensions

		double s = pow( 1.0  / errork, RKF_ee_c);

		s = (RKF_fac_c * s > RKF_facmin_c) ? RKF_fac_c * s : RKF_facmin_c;
		s = (RKF_facmax_c < s) ? RKF_facmax_c : s;

		snew = (snew < s) ? snew : s;

		snew_d[id] = snew;
//printf("id %d %g %g\n", id, s, snew);
	}
/*
	__syncthreads();

	if(snew >= 1.0){
		//accept step
		if(id < N){
			x_d[id] += dx;
			y_d[id] += dy;
			z_d[id] += dz;

			vx_d[id] += dvx;
			vy_d[id] += dvy;
			vz_d[id] += dvz;
		}
	}
*/
}

__global__ void computeError_d1_kernel(double *snew_d, double *ssum_d, const int N){

	int id = blockIdx.x * blockDim.x + threadIdx.x;

	double s = 1.0e6;	//large number

	extern __shared__ double se_s[];
	double *s_s = se_s;

	int lane = threadIdx.x % warpSize;
	int warp = threadIdx.x / warpSize;

	if(warp == 0){
		s_s[threadIdx.x] = 1.0e6;	//large number
	}
	__syncthreads();

	if(id < N){	
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


__global__ void update_kernel(double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *dx_d, double *dy_d, double *dz_d, double *dvx_d, double *dvy_d, double *dvz_d, const int N){

	int id = blockIdx.x * blockDim.x + threadIdx.x;

//	if(snew_d[id] >= 1.0){
//Add a flag for the different time step classes
		//accept step
		if(id < N){
			x_d[id] += dx_d[id];
			y_d[id] += dy_d[id];
			z_d[id] += dz_d[id];

			vx_d[id] += dvx_d[id];
			vy_d[id] += dvy_d[id];
			vz_d[id] += dvz_d[id];
		}
//	}
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


__global__ void convertOutput_kernel(double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *xout_d, double *yout_d, double *zout_d, double *vxout_d, double *vyout_d, double *vzout_d, double *xTable_d, double *yTable_d, double *zTable_d, double *vxTable_d, double *vyTable_d, double *vzTable_d, const int RKFn, const int N, const int Outecliptic, const int Outheliocentric, const int Outorbital, const double Obliquity, double Msun){

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
			int ii = 10 * RKFn;
			xout_d[id] -= xTable_d[ii];
			yout_d[id] -= yTable_d[ii];
			zout_d[id] -= zTable_d[ii];

			vxout_d[id] -= vxTable_d[ii];
			vyout_d[id] -= vyTable_d[ii];
			vzout_d[id] -= vzTable_d[ii];
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


			inc = inc * 180.0 / M_PI;       //convert rad to deg
			Omega = Omega * 180.0 / M_PI;   //convert rad to deg
			w = w * 180.0 / M_PI;           //convert rad to deg
			M = M * 180.0 / M_PI;           //convert rad to deg

			xout_d[id] = a;
			yout_d[id] = e;
			zout_d[id] = inc;
			vxout_d[id] = Omega;
			vyout_d[id] = w;
			vzout_d[id] = M;


		}
	}


}


	
int asteroid::loop(){

	//If needed, convert from heliocentric coordinates to barycentric coordinates
	if(ICheliocentric == 1){
		update_perturbers_kernel <<< RKFn, 32 >>>(xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, data_d, cdata_d, idp_d, startTime_d, endTime_d, nChebyshev_d, offset0_d, timeStart, time_reference, dt, RKFn, nCm, EM, AUtokm, Nperturbers);
		HelioToBary_kernel <<< (N + 255) / 256 , 256 >>> (xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, N, RKFn);

	}

	if(outBinary == 0){
		outputFile = fopen(outputFilename, "w");
	}
	else{
		outputFile = fopen(outputFilename, "wb");
	}

	printf("Start integration\n");

	if(time_reference + time >= outStart){
		if(Outheliocentric == 1){
			update_perturbers_kernel <<< RKFn, 32 >>>(xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, data_d, cdata_d, idp_d, startTime_d, endTime_d, nChebyshev_d, offset0_d, time, time_reference, dt, RKFn, nCm, EM, AUtokm, Nperturbers);

		}
		convertOutput_kernel <<< (N + 255) / 256 , 256 >>> (x_d, y_d, z_d, vx_d, vy_d, vz_d, xout_d, yout_d, zout_d, vxout_d, vyout_d, vzout_d, xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, RKFn, N, Outecliptic, Outheliocentric, Outorbital, Obliquity, GM_h[10]);
		copyOutput();
		printOutput(dt);
	}

	//for(int tt = 0; tt < 2; ++tt){
	for(int tt = 0; tt < 1000000; ++tt){
		double dtmin = dt;

		double timett1 = timeStart + dts * (tt + 1) * outputInterval;

		double snew = 10.0;

//printf("integrate %.20g %.20g\n", timeStart + dts * tt * 10.0, timett1);

		//integrate until the next output interval
		for(int ttt = 0; ttt < 1000000; ++ttt){

			//refine last time step of interval to match output time
			if(dts < 0){
				if(time + dt < timett1){
					dt1 = dt;
					dt = (timett1 - time);
					stop = 1;
//printf("refine %.20g\n", timett1 - time);
				}
			}
			else{
				if(time + dt > timett1){
					dt1 = dt;
					dt = (timett1 - time);
					stop = 1;
//printf("refine %.20g\n", timett1 - time);
				}

			}

			if(RKFn == 1){
				leapfrog_stepA_kernel <<< (N + 255) / 256 , 256 >>> (x_d, y_d, z_d, vx_d, vy_d, vz_d, N, dt);
				time += dt * 0.5;
				//Needs at least Nperturbers threads per block
				update_perturbers_kernel <<< RKFn, 32 >>>(xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, data_d, cdata_d, idp_d, startTime_d, endTime_d, nChebyshev_d, offset0_d, time, time_reference, dt, RKFn, nCm, EM, AUtokm, Nperturbers);

				//Needs at least Nperturbers threads per block
				leapfrog_stepB_kernel <<< (N + 255) / 256 , 256 >>> (xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, A1_d, A2_d, A3_d, Tsave_d, Rsave_d, time, timeStep, GM_d, Nperturbers, N, RKFn, dt);
				time += dt * 0.5;
				++timeStep;
			}
			if(RKFn == 4){
				//Needs at least Nperturbers threads per block
				update_perturbers_kernel <<< RKFn, 32 >>>(xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, data_d, cdata_d, idp_d, startTime_d, endTime_d, nChebyshev_d, offset0_d, time, time_reference, dt, RKFn, nCm, EM, AUtokm, Nperturbers);
	
				//Needs at least Nperturbers threads per block
				if(GPUMode == 0){
					RK_step_kernel <<< (N + 63) / 64 , dim3(64, 1, 1) >>> (xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, kx_d, ky_d, kz_d, kvx_d, kvy_d, kvz_d, GM_d, A1_d, A2_d, A3_d, Tsave_d, Rsave_d, time, timeStep, dt, RKFn, Nperturbers, N);
				}
				else{
					RK_step2_kernel <<< N, def_NP >>> (xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, GM_d, A1_d, A2_d, A3_d, Tsave_d, Rsave_d, time, timeStep, dt, RKFn, Nperturbers, N);
				}

//				for(int S = 0; S < RKFn; ++S){

//					RK_stage_kernel <<< (N + 255) / 256 , 256 >>> (xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d,vzTable_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, ax_d, ay_d, az_d, kx_d, ky_d, kz_d, kvx_d, kvy_d, kvz_d, GM_d, A1_d, A2_d, A3_d, Tsave_d, Rsave_d, time, dt, RKFn, Nperturbers, N, S);
//				}


				time += dt;
				++timeStep;
			}
			if(RKFn == 6 || RKFn == 7 || RKFn == 13){
				//Needs at least Nperturbers threads per block
				update_perturbers_kernel <<< RKFn, 32 >>>(xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, data_d, cdata_d, idp_d, startTime_d, endTime_d, nChebyshev_d, offset0_d, time, time_reference, dt, RKFn, nCm, EM, AUtokm, Nperturbers);

				if(GPUMode == 0){
					RKF_step_kernel <<< (N + 63) / 64 , 64 >>> (xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, dx_d, dy_d, dz_d, dvx_d, dvy_d, dvz_d, kx_d, ky_d, kz_d, kvx_d, kvy_d, kvz_d, GM_d, A1_d, A2_d, A3_d, Tsave_d, Rsave_d, snew_d, time, timeStep, dt, dts, RKFn, Nperturbers, N, stop);
				}
				else{
					RKF_step2_kernel <<< N, def_NP >>> (xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, dx_d, dy_d, dz_d, dvx_d, dvy_d, dvz_d, GM_d, A1_d, A2_d, A3_d, Tsave_d, Rsave_d, snew_d, time, timeStep, dt, dts, RKFn, Nperturbers, N, stop);
				}


				//Calculate the minimal time step value
				//Using a parallel reduction sum
				int nct = 512;
				int ncb = min((N + nct - 1) / nct, 1024);
				computeError_d1_kernel <<< ncb, nct, WarpSize * sizeof(double)  >>> (snew_d, ssum_d, N);
				if(ncb > 1){
					computeError_d2_kernel <<< 1, ((ncb + WarpSize - 1) / WarpSize) * WarpSize, WarpSize * sizeof(double)  >>> (ssum_d, ncb);
				}

				cudaDeviceSynchronize();
				cudaMemcpy(snew_h, ssum_d, sizeof(double), cudaMemcpyDeviceToHost);
				snew = snew_h[0];

				if(snew >= 1.0){
					//accept step
					update_kernel <<< (N + 63) / 64 , 64 >>> (x_d, y_d, z_d, vx_d, vy_d, vz_d, dx_d, dy_d, dz_d, dvx_d, dvy_d, dvz_d, N);
					time += dt;
					++timeStep;
					if(stop != 1){
						 //only increase time step when stop == 0
						dt *= snew;
					}
				}
				else{
					//redo step
					dt *= snew;
				}

			}

			if(stop == 0){
				dtmin = (abs(dt) < abs(dtmin)) ? dt : dtmin;
			}

			if(time + time_reference > time1 || time + time_reference < time0){
				cudaDeviceSynchronize();
				printf("Reached the end of the Chebyshev data file\n");
				return 0;
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

			if(stop == 1){
				stop = 0;
				if(snew >= 1.0){
					//set time step equal to the last accepted full time step
					dt = dt1;
					break;
				}
			}

			if(ttt >= 1000000 - 1){

				printf("Error, time step loop did not finish\n");
				return 0;
			}

		}//end of ttt loop
		cudaDeviceSynchronize();
		if(time_reference + time >= outStart){
			if(Outheliocentric == 1){
				update_perturbers_kernel <<< RKFn, 32 >>>(xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, data_d, cdata_d, idp_d, startTime_d, endTime_d, nChebyshev_d, offset0_d, time, time_reference, dt, RKFn, nCm, EM, AUtokm, Nperturbers);
			}
			convertOutput_kernel <<< (N + 255) / 256 , 256 >>> (x_d, y_d, z_d, vx_d, vy_d, vz_d, xout_d, yout_d, zout_d, vxout_d, vyout_d, vzout_d, xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, RKFn, N, Outecliptic, Outheliocentric, Outorbital, Obliquity, GM_h[10]);
			copyOutput();
			printOutput(dtmin);
		}

	}//end of tt loop
	fclose(outputFile);
	return 1;
}

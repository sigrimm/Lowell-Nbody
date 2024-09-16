//Define Constant Memory
__constant__ double a_c[20 * 20];       //20 is considered here to be large enough (>RKFn)
__constant__ double b_c[20];
__constant__ double bb_c[20];
__constant__ double c_c[20];

#include "asteroid.h"
#include "ChebyshevGPU.h"
#include "forceGPU.h"


__host__ int asteroid::copyConst(){
        cudaMemcpyToSymbol(a_c, a_h, RKFn * RKFn * sizeof(double), 0, cudaMemcpyHostToDevice);
        cudaMemcpyToSymbol(b_c, b_h, RKFn * sizeof(double), 0, cudaMemcpyHostToDevice);
        cudaMemcpyToSymbol(bb_c, bb_h, RKFn * sizeof(double), 0, cudaMemcpyHostToDevice);
        cudaMemcpyToSymbol(c_c, c_h, RKFn * sizeof(double), 0, cudaMemcpyHostToDevice);

	cudaDeviceSynchronize();
	cudaError_t error = cudaGetLastError();
	printf("copy error = %d = %s\n",error, cudaGetErrorString(error));
	if(error != 0.0){
		return 0;
	}

	return 1;

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
__global__ void leapfrog_stepB_kernel(double *xTable_d, double *yTable_d, double *zTable_d, double *vxTable_d, double *vyTable_d, double *vzTable_d, double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *A1_d, double *A2_d, double *A3_d, double *GM_d, const int Nperturbers, const int N, const int RKFn, const double dt, const double c2, const double J2E, const double REAU){

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
		xTable_s[itx] = xTable_d[itx * RKFn];
		yTable_s[itx] = yTable_d[itx * RKFn];
		zTable_s[itx] = zTable_d[itx * RKFn];

		vxTable_s[itx] = vxTable_d[itx * RKFn];
		vyTable_s[itx] = vyTable_d[itx * RKFn];
		vzTable_s[itx] = vzTable_d[itx * RKFn];

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
		double rsq = xih * xih + yih * yih + zih * zih;
		double r = sqrt(rsq);

		//Earth centric coordinates
		double xiE = xi - xTable_s[2];
		double yiE = yi - yTable_s[2];
		double ziE = zi - zTable_s[2];

		NonGrav(xih, yih, zih, vxih, vyih, vzih, A1_d[id], A2_d[id], A3_d[id], r, ax, ay, az);
		GR(xih, yih, zih, vxih, vyih, vzih, r, ax, ay, az, GM_s[10], c2);
		J2(xiE, yiE, ziE, ax, ay, az, REAU, J2E, GM_s[2]);
		Gravity(xi, yi, zi, xTable_s, yTable_s, zTable_s, ax, ay, az, GM_s, Nperturbers);
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


//Runge Kutta step with fixed time step
//Kernel uses at least Nperturbers threads. The perturbers are loaded into shared memory
__global__ void RK_step_kernel(double *xTable_d, double *yTable_d, double *zTable_d, double *vxTable_d, double *vyTable_d, double *vzTable_d, double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *kx_d, double *ky_d, double *kz_d, double *kvx_d, double *kvy_d, double *kvz_d, double *GM_d, double *A1_d, double *A2_d, double *A3_d, double time, double time_reference, const double dt, const double REAU, const double J2E, const double c2, const int RKFn, const int Nperturbers, const int N){


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
//printf("p %d %.20g %.20g %.20g %.20g %.20g %.20g %.20g\n", id, time + c_c[S] * dt, xt_s[id], yt_s[id], zt_s[id], vxt_s[id], vyt_s[id], vzt_s[id]);
//}

		if(id < N){

			xti = x0;
			yti = y0;
			zti = z0;

			vxti = vx0;
			vyti = vy0;
			vzti = vz0;

			for(int s = 0; s < S; ++s){
				double dtaa = dt * a_c[S * RKFn + s];
				xti  += dtaa * kx_d[id + s * N];
				yti  += dtaa * ky_d[id + s * N];
				zti  += dtaa * kz_d[id + s * N];
				vxti += dtaa * kvx_d[id + s * N];
				vyti += dtaa * kvy_d[id + s * N];
				vzti += dtaa * kvz_d[id + s * N];
//printf("update 2 %d %d %g %g %g %g %g %g\n", S, id, xti, yti, zti, a_c[S * RKFn + s], kx_d[s], dt);

			}

			kx_d[id + S * N] = vxti;
			ky_d[id + S * N] = vyti;
			kz_d[id + S * N] = vzti;


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
			double rsq = xih * xih + yih * yih + zih * zih;
			double r = sqrt(rsq);

			//Earth centric coordinates
			double xiE = xti - xTable_s[2];
			double yiE = yti - yTable_s[2];
			double ziE = zti - zTable_s[2];

			NonGrav(xih, yih, zih, vxih, vyih, vzih, A1, A2, A3, r, ax, ay, az);
			GR(xih, yih, zih, vxih, vyih, vzih, r, ax, ay, az, GM_s[10], c2);
			J2(xiE, yiE, ziE, ax, ay, az, REAU, J2E, GM_s[2]);
			Gravity(xti, yti, zti, xTable_s, yTable_s, zTable_s, ax, ay, az, GM_s, Nperturbers);
			// ----------------------------------------------------------------------------

			kvx_d[id + S * N] = ax;
			kvy_d[id + S * N] = ay;
			kvz_d[id + S * N] = az;
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
			double dtb = dt * b_c[S];
			dx += dtb * kx_d[id + S * N];
			dy += dtb * ky_d[id + S * N];
			dz += dtb * kz_d[id + S * N];

			dvx += dtb * kvx_d[id + S * N];
			dvy += dtb * kvy_d[id + S * N];
			dvz += dtb * kvz_d[id + S * N];
		}

		x_d[id] += dx;
		y_d[id] += dy;
		z_d[id] += dz;

		vx_d[id] += dvx;
		vy_d[id] += dvy;
		vz_d[id] += dvz;

//printf("%d %.20g %.20g %g %g %g\n", id, time, time + time_reference, x_d[id], y_d[id], z_d[id]);
	}
}



//Runge Kutta Fehlberg step with adaptive time step
__global__ void RKF_step_kernel(double *xTable_d, double *yTable_d, double *zTable_d, double *vxTable_d, double *vyTable_d, double *vzTable_d, double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *kx_d, double *ky_d, double *kz_d, double *kvx_d, double *kvy_d, double *kvz_d, double *GM_d, double *A1_d, double *A2_d, double *A3_d, double *snew_d, double time, double time_reference, const double dt, const int dts, const double REAU, const double J2E, const double c2, const int RKFn, const int Nperturbers, const int N, const double RKF_atol, const double RKF_rtol, const double RKF_ee, const double RKF_fac, const double RKF_facmin, const double RKF_facmax, const int stop){


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

	if(id < Nperturbers){
		GM_s[id] = GM_d[id];
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
//printf("p %d %.20g %.20g %.20g %.20g %.20g %.20g %.20g\n", id, time + c_c[S] * dt, xt_s[id], yt_s[id], zt_s[id], vxt_s[id], vyt_s[id], vzt_s[id]);
//}

		if(id < N){

			xti = x0;
			yti = y0;
			zti = z0;

			vxti = vx0;
			vyti = vy0;
			vzti = vz0;

			for(int s = 0; s < S; ++s){
				double dtaa = dt * a_c[S * RKFn + s];
				xti  += dtaa * kx_d[id + s * N];
				yti  += dtaa * ky_d[id + s * N];
				zti  += dtaa * kz_d[id + s * N];
				vxti += dtaa * kvx_d[id + s * N];
				vyti += dtaa * kvy_d[id + s * N];
				vzti += dtaa * kvz_d[id + s * N];
//printf("update 2 %d %d %g %g %g %g %g %g\n", S, id, xti, yti, zti, a_c[S * RKFn + s], kx_d[s], dt);

			}

			kx_d[id + S * N] = vxti;
			ky_d[id + S * N] = vyti;
			kz_d[id + S * N] = vzti;

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
			double rsq = xih * xih + yih * yih + zih * zih;
			double r = sqrt(rsq);

			//Earth centric coordinates
			double xiE = xti - xTable_s[2];
			double yiE = yti - yTable_s[2];
			double ziE = zti - zTable_s[2];

			NonGrav(xih, yih, zih, vxih, vyih, vzih, A1, A2, A3, r, ax, ay, az);
			GR(xih, yih, zih, vxih, vyih, vzih, r, ax, ay, az, GM_s[10], c2);
			J2(xiE, yiE, ziE, ax, ay, az, REAU, J2E, GM_s[2]);
			Gravity(xti, yti, zti, xTable_s, yTable_s, zTable_s, ax, ay, az, GM_s, Nperturbers);
			// ----------------------------------------------------------------------------

			kvx_d[id + S * N] = ax;
			kvy_d[id + S * N] = ay;
			kvz_d[id + S * N] = az;
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
			double dtb = dt * b_c[S];
			dx += dtb * kx_d[id + S * N];
			dy += dtb * ky_d[id + S * N];
			dz += dtb * kz_d[id + S * N];

			dvx += dtb * kvx_d[id + S * N];
			dvy += dtb * kvy_d[id + S * N];
			dvz += dtb * kvz_d[id + S * N];
		}


		//compute integration error

		double ym = 0.0;
		ym = (fabs(x_d[id]) > ym) ? fabs(x_d[id]) : ym;
		ym = (fabs(y_d[id]) > ym) ? fabs(y_d[id]) : ym;
		ym = (fabs(z_d[id]) > ym) ? fabs(z_d[id]) : ym;

		ym = (fabs(vx_d[id]) > ym) ? fabs(vx_d[id]) : ym;
		ym = (fabs(vy_d[id]) > ym) ? fabs(vy_d[id]) : ym;
		ym = (fabs(vz_d[id]) > ym) ? fabs(vz_d[id]) : ym;

		double isc = 1.0 / (RKF_atol + ym * RKF_rtol);

		//error estimation
		double errorkx = 0.0;
		double errorky = 0.0;
		double errorkz = 0.0;
		
		double errorkvx = 0.0;
		double errorkvy = 0.0;
		double errorkvz = 0.0;

		for(int S = 0; S < RKFn; ++S){
			double f = (b_c[S] - bb_c[S]) * dt;
			errorkx += f * kx_d[id + S * N];
			errorky += f * ky_d[id + S * N];
			errorkz += f * kz_d[id + S * N];

			errorkvx += f * kvx_d[id + S * N];
			errorkvy += f * kvy_d[id + S * N];
			errorkvz += f * kvz_d[id + S * N];
//printf("error %d %d %g %g\n", id, S, errorkx, kx[S]);
		}

		double errork = 0.0;
		errork += errorkx * errorkx * isc * isc;
		errork += errorky * errorky * isc * isc;
		errork += errorkz * errorkz * isc * isc;
		errork += errorkvx * errorkvx * isc * isc;
		errork += errorkvy * errorkvy * isc * isc;
		errork += errorkvz * errorkvz * isc * isc;

		errork = sqrt(errork / 6.0);    //6 is the number of dimensions

		double s = pow( 1.0  / errork, RKF_ee);
		s = fmax(RKF_facmin, RKF_fac * s);
		s = fmin(RKF_facmax, s);

		snew = (snew < s) ? snew : s;

		snew_d[0] = snew;
//printf("id %d %g %g\n", id, s, snew);
	}

	__syncthreads();

	if(stop == 1){
		//accept step
		if(id < N){
			x_d[id] += dx;
			y_d[id] += dy;
			z_d[id] += dz;

			vx_d[id] += dvx;
			vy_d[id] += dvy;
			vz_d[id] += dvz;

			snew_d[0] = 1.0;
		}
	}
	else if(snew >= 1.0){
		//accept step
		if(id < N){
			x_d[id] += dx;
			y_d[id] += dy;
			z_d[id] += dz;

			vx_d[id] += dvx;
			vy_d[id] += dvy;
			vz_d[id] += dvz;

		}
		//dt *= snew;

		//set maximum time step to 1
		//if(abs(dt) > 1.0){
		//	dt = dts * 1.0;
		//}

	}
	else{
		//redo step
		//dt *= snew;
	}
}


	
int asteroid::loop(){


	outputFile = fopen("Out.dat", "w");
	copyOutput();

	for(int p = 0; p < N; ++p){
		printf("Start integration\n");
		printf("Reached time %.20g dtmin %.8g\n", time_reference + time, dt);
		fprintf(outputFile, "%.20g %d %.20g %.20g %.20g %.20g %.20g %.20g %.20g\n", time_reference + time, p, x_h[p], y_h[p], z_h[p], vx_h[p], vy_h[p], vz_h[p], dt);
	}
	//for(int tt = 0; tt < 2; ++tt){
	for(int tt = 0; tt < 1000000; ++tt){
		double dtmin = dt;

		double timett1 = timeStart + dts * (tt + 1) * outputInterval;

//printf("integrate %.20g %.20g\n", timeStart + dts * tt * 10.0, timett1);


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
//printf("refine %.20g\n", timett1 - A.time);
				}

			}

			if(RKFn == 1){
				leapfrog_stepA_kernel <<< (N + 255) / 256 , 256 >>> (x_d, y_d, z_d, vx_d, vy_d, vz_d, N, dt);
				time += dt * 0.5;
				//Needs at least Nperturbers threads per block
				update_perturbers_kernel <<< RKFn, 32 >>>(xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, data_d, cdata_d, id_d, startTime_d, endTime_d, nChebyshev_d, offset0_d, time, time_reference, dt, RKFn, nCm, EM, AUtokm, Nperturbers);

				//Needs at least Nperturbers threads per block
				leapfrog_stepB_kernel <<< (N + 255) / 256 , 256 >>> (xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, A1_d, A2_d, A3_d, GM_d, Nperturbers, N, RKFn, dt, c2, J2E, REAU);
				time += dt * 0.5;
			}
			if(RKFn == 4){
				//Needs at least Nperturbers threads per block
				update_perturbers_kernel <<< RKFn, 32 >>>(xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, data_d, cdata_d, id_d, startTime_d, endTime_d, nChebyshev_d, offset0_d, time, time_reference, dt, RKFn, nCm, EM, AUtokm, Nperturbers);
				//Needs at least Nperturbers threads per block
				RK_step_kernel <<< (N + 255) / 256 , 256 >>> (xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, kx_d, ky_d, kz_d, kvx_d, kvy_d, kvz_d, GM_d, A1_d, A2_d, A3_d, time, time_reference, dt, REAU, J2E, c2, RKFn, Nperturbers, N);
				time += dt;
			}
			if(RKFn == 6 || RKFn == 7 || RKFn == 13){
				//Needs at least Nperturbers threads per block
				update_perturbers_kernel <<< RKFn, 32 >>>(xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, data_d, cdata_d, id_d, startTime_d, endTime_d, nChebyshev_d, offset0_d, time, time_reference, dt, RKFn, nCm, EM, AUtokm, Nperturbers);
				RKF_step_kernel <<< (N + 255) / 256 , 256 >>> (xTable_d, yTable_d, zTable_d, vxTable_d, vyTable_d, vzTable_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, kx_d, ky_d, kz_d, kvx_d, kvy_d, kvz_d, GM_d, A1_d, A2_d, A3_d, snew_d, time, time_reference, dt, dts, REAU, J2E, c2, RKFn, Nperturbers, N, RKF_atol, RKF_rtol, RKF_ee, RKF_fac, RKF_facmin, RKF_facmax, stop);
				cudaMemcpy(snew_h, snew_d, sizeof(double), cudaMemcpyDeviceToHost);
				cudaDeviceSynchronize();
				double snew = snew_h[0];

				if(snew >= 1.0){
					time += dt;
				}
				dt *= snew;


			}

			dtmin = (abs(dt) < abs(dtmin)) ? dt : dtmin;

			if(time + time_reference > time1 || time + time_reference < time0){
				printf("Reached the end of the Chebyshev data file\n");
				return 0;
			}

			if(dts < 0 && time < timeEnd){
				printf("Reached the end of the integration\n");
				return 0;
			}
			if(dts > 0 && time > timeEnd){
				printf("Reached the end of the integration\n");
				return 0;
			}

			if(stop == 1){
				stop = 0;
				dt = dt1;
				break;
			}

			if(ttt >= 1000000 - 1){

				printf("Error, time step loop did not finish\n");
				return 0;
			}

		}//end of ttt loop
		cudaDeviceSynchronize();
		copyOutput();
		for(int p = 0; p < N; ++p){
			printf("Reached time %.20g dtmin %.8g\n", time_reference + time, dt);
			fprintf(outputFile, "%.20g %d %.20g %.20g %.20g %.20g %.20g %.20g %.20g\n", time_reference + time, p, x_h[p], y_h[p], z_h[p], vx_h[p], vy_h[p], vz_h[p], dtmin);
		}

	}//end of tt loop
	fclose(outputFile);
	return 1;
}

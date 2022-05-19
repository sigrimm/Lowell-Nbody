#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//for FIFO
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "define.h"

#include "read.h"
#include "perturbers.h"
#include "force.h"
#include "integrator.h"
#include "interpolate.h"
#include "output.h"


__host__ void stageStep(unsigned long long int *id_h, double *m_h, double *x_h, double *y_h, double *z_h, double *vx_h, double *vy_h, double *vz_h, double *dx_h, double *dy_h, double *dz_h, double *dvx_h, double *dvy_h, double *dvz_h, double *xTable_h, double *yTable_h, double *zTable_h, double *kx_h, double *ky_h, double *kz_h, double *kvx_h, double *kvy_h, double *kvz_h, double *A1_h, double *A2_h, double *A3_h, double2 *snew_h, double *a_h, double *b_h, double *bb_h, double dt, double dti, double dtiMin, int RKFn, int Nperturbers, int N, int useHelio, int useGR, int useJ2, int useNonGrav, double ee, double &snew){


	double xt[Nperturbers];
	double yt[Nperturbers];
	double zt[Nperturbers];

	for(int S = 0; S < RKFn; ++S){

		for(int i = 0; i < Nperturbers; ++i){
			int ii = i * RKFn + S;
			xt[i] = xTable_h[ii];
			yt[i] = yTable_h[ii];
			zt[i] = zTable_h[ii];
//printf("%d %d %.20g\n", i, S, xt[i]);
		}

		for(int i = Nperturbers; i < N; ++i){

			//update
			double xi = x_h[i];
			double yi = y_h[i];
			double zi = z_h[i];
			double vxi = vx_h[i];
			double vyi = vy_h[i];
			double vzi = vz_h[i];

			for(int s = 0; s < S; ++s){
				double dtaa = dt * a_h[S * RKFn + s];
				xi  += dtaa * kx_h[i + s * N];
				yi  += dtaa * ky_h[i + s * N];
				zi  += dtaa * kz_h[i + s * N];
				vxi += dtaa * kvx_h[i + s * N];
				vyi += dtaa * kvy_h[i + s * N];
				vzi += dtaa * kvz_h[i + s * N];
//printf("x %d %d %.20g %.20g\n", S, s, dtaa, xi);
			}



			kx_h[i + S * N] = vxi;
			ky_h[i + S * N] = vyi;
			kz_h[i + S * N] = vzi;

			double ax = 0.0;
			double ay = 0.0;
			double az = 0.0;

			if(useHelio == 0){
				for(int j = Nperturbers-1; j >= 0; --j){
					if(id_h[i] != id_h[j]){
						accP(m_h[j], xt[j], yt[j], zt[j], xi, yi, zi, ax, ay, az);
					}
				}
			}
			else{
				for(int j = Nperturbers-1; j >= 1; --j){
					if(id_h[i] != id_h[j]){
						accP(m_h[j], xt[j], yt[j], zt[j], xi, yi, zi, ax, ay, az);
//printf("Nij %d %d %llu %llu %.20g %.20g %.20g\n", i, j, id_h[i], id_h[j], ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit);
					}
				}
				accS(m_h[0] + m_h[i], xi, yi, zi, ax, ay, az);
//if(i == 27) printf("Nij %d %d %.20g %.20g %.20g\n", i, 0, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit);
//printf("N0 %d %.20g %.20g %.20g\n", i, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit);
				for(int j = Nperturbers-1; j >= 1; --j){
					if(id_h[i] != id_h[j]){
						accP2(m_h[j], xt[j], yt[j], zt[j], ax, ay, az);
					}
				}
//if(id_h[i] == 72057594037948950) printf("Np %d %.20g %.20g %.20g %d\n", i, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit, S);
			}


			if(useGR == 2){
				acchGR2(xi, yi, zi, vxi, vyi, vzi, ax, ay, az);
			}

			if(useNonGrav == 1){
				NonGrav(xi, yi, zi, vxi, vyi, vzi, ax, ay, az, A1_h[i], A2_h[i], A3_h[i]);
			}

			if(useJ2 == 1){
				J2(m_h[3], xt[3], yt[3], zt[3], xi, yi, zi, ax, ay, az);
			}

			kvx_h[i + S * N] = ax;
			kvy_h[i + S * N] = ay;
			kvz_h[i + S * N] = az;
		}
	}


	snew = 1.0e6;	//large number
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


		//compute error
		double ym = fabs(x_h[i]);
		ym = fmax(ym, fabs(y_h[i]));
		ym = fmax(ym, fabs(z_h[i]));
		ym = fmax(ym, fabs(vx_h[i]));
		ym = fmax(ym, fabs(vy_h[i]));
		ym = fmax(ym, fabs(vz_h[i]));

		ym = fmax(ym, fabs(x_h[i] + dx));
		ym = fmax(ym, fabs(y_h[i] + dy));
		ym = fmax(ym, fabs(z_h[i] + dz));
		ym = fmax(ym, fabs(vx_h[i] + dvx));
		ym = fmax(ym, fabs(vy_h[i] + dvy));
		ym = fmax(ym, fabs(vz_h[i] + dvz));


		double isc = 1.0 / (def_atol + ym * def_rtol);

		double s = 1.0e6;       //large number

		double errorkx = 0.0;
		double errorky = 0.0;
		double errorkz = 0.0;
		double errorkvx = 0.0;
		double errorkvy = 0.0;
		double errorkvz = 0.0;

		for(int S = 0; S < RKFn; ++S){
			//this is y1i - y^1i
			double f = (b_h[S] - bb_h[S]) * dt;
			errorkx += f * kx_h[S * N + i];
			errorky += f * ky_h[S * N + i];
			errorkz += f * kz_h[S * N + i];

			errorkvx += f * kvx_h[S * N + i];
			errorkvy += f * kvy_h[S * N + i];
			errorkvz += f * kvz_h[S * N + i];
		}

		double errork = 0.0;
		errork += errorkx * errorkx * isc * isc;
		errork += errorky * errorky * isc * isc;
		errork += errorkz * errorkz * isc * isc;
		errork += errorkvx * errorkvx * isc * isc;
		errork += errorkvy * errorkvy * isc * isc;
		errork += errorkvz * errorkvz * isc * isc;

		errork = sqrt(errork / 6.0);    //6 is the number of dimensions

		s = pow( 1.0  / errork, ee);
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

		snew = fmin(snew, snew_h[i].x);
//printf("snew %d %.20g %.20g %.20g %.20g\n", i, snew_h[i].y, snew_h[i].x, errork, snew);
		// ****************************

	}

}


template < int Nperturbers >
__global__ void stageStep1_kernel(unsigned long long int *id_d, double *m_d, double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *dx_d, double *dy_d, double *dz_d, double *dvx_d, double *dvy_d, double *dvz_d, double *xTable_d, double *yTable_d, double *zTable_d, double *kx_d, double *ky_d, double *kz_d, double *kvx_d, double *kvy_d, double *kvz_d, double *A1_d, double *A2_d, double *A3_d, double2 *snew_d, double dt, double dti, double dtiMin, int RKFn, int N, int useHelio, int useGR, int useJ2, int useNonGrav, double ee){

	int itx = threadIdx.x;
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	//shared memory contains only the perturbers
	//the particle idx is stored in registers
	__shared__ double x_s[Nperturbers];
	__shared__ double y_s[Nperturbers];
	__shared__ double z_s[Nperturbers];
	__shared__ double m_s[Nperturbers];
	__shared__ unsigned long long int id_s[Nperturbers];

	double xi0 = x_d[idx];
	double yi0 = y_d[idx];
	double zi0 = z_d[idx];
	double vxi0 = vx_d[idx];
	double vyi0 = vy_d[idx];
	double vzi0 = vz_d[idx];

	double mi = m_d[idx];
	double A1i = A1_d[idx];
	double A2i = A2_d[idx];
	double A3i = A3_d[idx];
	unsigned long long int idi = id_d[idx];

	for(int S = 0; S < RKFn; ++S){


		if(itx < Nperturbers){
			int ii = itx * RKFn + S;
			x_s[itx] = xTable_d[ii];
			y_s[itx] = yTable_d[ii];
			z_s[itx] = zTable_d[ii];
			m_s[itx] = m_d[itx];
			id_s[itx] = id_d[itx];
		}
		__syncthreads();

		if(idx >= Nperturbers && idx < N){


			// ***********************
			// update
			double xi = xi0;
			double yi = yi0;
			double zi = zi0;
			double vxi = vxi0;
			double vyi = vyi0;
			double vzi = vzi0;

			for(int s = 0; s < S; ++s){
				double dtaa = dt * a_c[S * RKFn + s];
				xi  += dtaa * kx_d[idx + s * N];
				yi  += dtaa * ky_d[idx + s * N];
				zi  += dtaa * kz_d[idx + s * N];
				vxi += dtaa * kvx_d[idx + s * N];
				vyi += dtaa * kvy_d[idx + s * N];
				vzi += dtaa * kvz_d[idx + s * N];
//printf("x %d %d %.20g %.20g\n", S, s, dtaa, xi);
			}



			// end update
			// *****************************

			kx_d[idx + S * N] = vxi;
			ky_d[idx + S * N] = vyi;
			kz_d[idx + S * N] = vzi;

			double ax = 0.0;
			double ay = 0.0;
			double az = 0.0;

			if(useHelio == 0){
				for(int j = Nperturbers-1; j >= 0; --j){
					if(idi != id_s[j]){
						accP(m_s[j], x_s[j], y_s[j], z_s[j], xi, yi, zi, ax, ay, az);
					}
				}
			}
			else{
				for(int j = Nperturbers-1; j >= 1; --j){
					if(idi != id_s[j]){
						accP(m_s[j], x_s[j], y_s[j], z_s[j], xi, yi, zi, ax, ay, az);
//printf("Nij %d %d %.20g %.20g %.20g\n", idx, j, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit);
					}
				}
				accS(m_s[0] + mi, xi, yi, zi, ax, ay, az);

//printf("N0 %d %.20g %.20g %.20g\n", idx, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit);
				for(int j = Nperturbers-1; j >= 1; --j){
					if(idi != id_s[j]){
						accP2(m_s[j], x_s[j], y_s[j], z_s[j], ax, ay, az);
//printf("Npi %d %d %.20g %.20g %.20g %d\n", idx, j, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit, S);
					}
				}
//printf("Np %d %.20g %.20g %.20g %d\n", idx, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit, S);
			}


			if(useGR == 2){
				acchGR2(xi, yi, zi, vxi, vyi, vzi, ax, ay, az);
			}

			if(useNonGrav == 1){
				NonGrav(xi, yi, zi, vxi, vyi, vzi, ax, ay, az, A1i, A2i, A3i);
			}

			if(useJ2 == 1){
				J2(m_s[3], x_s[3], y_s[3], z_s[3], xi, yi, zi, ax, ay, az);
			}

			kvx_d[idx + S * N] = ax;
			kvy_d[idx + S * N] = ay;
			kvz_d[idx + S * N] = az;
		}
	}

	if(idx >= Nperturbers && idx < N){
		//update
		double dx = 0.0;
		double dy = 0.0;
		double dz = 0.0;
		double dvx = 0.0;
		double dvy = 0.0;
		double dvz = 0.0;
		
		for(int S = 0; S < RKFn; ++S){
			double dtb = dt * b_c[S];
			dx += dtb * kx_d[idx + S * N];
			dy += dtb * ky_d[idx + S * N];
			dz += dtb * kz_d[idx + S * N];

			dvx += dtb * kvx_d[idx + S * N];
			dvy += dtb * kvy_d[idx + S * N];
			dvz += dtb * kvz_d[idx + S * N];
		}
		
		dx_d[idx] = dx;
		dy_d[idx] = dy;
		dz_d[idx] = dz;
		dvx_d[idx] = dvx;
		dvy_d[idx] = dvy;
		dvz_d[idx] = dvz;

		//compute error
		double ym = fabs(xi0);
		ym = fmax(ym, fabs(yi0));
		ym = fmax(ym, fabs(zi0));
		ym = fmax(ym, fabs(vxi0));
		ym = fmax(ym, fabs(vyi0));
		ym = fmax(ym, fabs(vzi0));

		ym = fmax(ym, fabs(xi0 + dx));
		ym = fmax(ym, fabs(yi0 + dy));
		ym = fmax(ym, fabs(zi0 + dz));
		ym = fmax(ym, fabs(vxi0 + dvx));
		ym = fmax(ym, fabs(vyi0 + dvy));
		ym = fmax(ym, fabs(vzi0 + dvz));


		double isc = 1.0 / (def_atol + ym * def_rtol);

		double s = 1.0e6;       //large number

		double errorkx = 0.0;
		double errorky = 0.0;
		double errorkz = 0.0;
		double errorkvx = 0.0;
		double errorkvy = 0.0;
		double errorkvz = 0.0;

		for(int S = 0; S < RKFn; ++S){
			//this is y1i - y^1i
			double f = (b_c[S] - bb_c[S]) * dt;
			errorkx += f * kx_d[S * N + idx];
			errorky += f * ky_d[S * N + idx];
			errorkz += f * kz_d[S * N + idx];

			errorkvx += f * kvx_d[S * N + idx];
			errorkvy += f * kvy_d[S * N + idx];
			errorkvz += f * kvz_d[S * N + idx];
		}

		double errork = 0.0;
		errork += errorkx * errorkx * isc * isc;
		errork += errorky * errorky * isc * isc;
		errork += errorkz * errorkz * isc * isc;
		errork += errorkvx * errorkvx * isc * isc;
		errork += errorkvy * errorkvy * isc * isc;
		errork += errorkvz * errorkvz * isc * isc;

		errork = sqrt(errork / 6.0);    //6 is the number of dimensions

		s = pow( 1.0  / errork, ee);
		s = fmax(def_facmin, def_fac * s);
		s = fmin(def_facmax, s);

		// ****************************
		if(snew_d[idx].y >= 1.0){
			snew_d[idx].x = s;
		}
		else{
			snew_d[idx].x = 1.5;
		}
		if(s * dti < dtiMin){
			snew_d[idx].y = fmin(snew_d[idx].y, s);
		}
//printf("snew %d %.20g %.20g %.20g\n", idx, snew_d[idx].y, snew_d[idx].x, errork);
		// ****************************
	}


}

//every particle runs on an own thread block
//the force calculation is parallelized along the threads and reduced in registers
template < int Nperturbers >
__global__ void stageStep2_kernel(unsigned long long int *id_d, double *m_d, double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *dx_d, double *dy_d, double *dz_d, double *dvx_d, double *dvy_d, double *dvz_d, double *xTable_d, double *yTable_d, double *zTable_d, double *kx_d, double *ky_d, double *kz_d, double *kvx_d, double *kvy_d, double *kvz_d, double *A1_d, double *A2_d, double *A3_d, double2 *snew_d, double dt, double dti, double dtiMin, int RKFn, int N, int useHelio, int useGR, int useJ2, int useNonGrav, double ee){

	int itx = threadIdx.x;
	int idx = blockIdx.x + Nperturbers;

	//shared memory contains only the perturbers
	__shared__ double x_s[Nperturbers];
	__shared__ double y_s[Nperturbers];
	__shared__ double z_s[Nperturbers];
	__shared__ double m_s[Nperturbers];
	__shared__ unsigned long long int id_s[Nperturbers];


	__shared__ double xi_s[1];
	__shared__ double yi_s[1];
	__shared__ double zi_s[1];
	__shared__ double vxi_s[1];
	__shared__ double vyi_s[1];
	__shared__ double vzi_s[1];

	if(idx < N){
		double xi0 = x_d[idx];
		double yi0 = y_d[idx];
		double zi0 = z_d[idx];
		double vxi0 = vx_d[idx];
		double vyi0 = vy_d[idx];
		double vzi0 = vz_d[idx];

		double mi = m_d[idx];
		double A1i = A1_d[idx];
		double A2i = A2_d[idx];
		double A3i = A3_d[idx];
		unsigned long long int idi = id_d[idx];

		for(int S = 0; S < RKFn; ++S){

			if(threadIdx.x < Nperturbers){
					int ii = itx * RKFn + S;
					x_s[threadIdx.x] = xTable_d[ii];
					y_s[threadIdx.x] = yTable_d[ii];
					z_s[threadIdx.x] = zTable_d[ii];
					m_s[threadIdx.x] = m_d[itx];
					id_s[threadIdx.x] = id_d[itx];
//printf("%d %d %.20g\n", itx, S, x_s[itx]);
			}
			__syncthreads();

			if(itx == 0){

				// ***********************
				// update
				xi_s[0] = xi0;
				yi_s[0] = yi0;
				zi_s[0] = zi0;
				vxi_s[0] = vxi0;
				vyi_s[0] = vyi0;
				vzi_s[0] = vzi0;

				for(int s = 0; s < S; ++s){
					double dtaa = dt * a_c[S * RKFn + s];
					xi_s[0]  += dtaa * kx_d[idx + s * N];
					yi_s[0]  += dtaa * ky_d[idx + s * N];
					zi_s[0]  += dtaa * kz_d[idx + s * N];
					vxi_s[0] += dtaa * kvx_d[idx + s * N];
					vyi_s[0] += dtaa * kvy_d[idx + s * N];
					vzi_s[0] += dtaa * kvz_d[idx + s * N];
	//printf("x %d %d %.20g %.20g\n", S, s, aa, xi);
				}


				// end update
				// *****************************

				kx_d[idx + S * N] = vxi_s[0];
				ky_d[idx + S * N] = vyi_s[0];
				kz_d[idx + S * N] = vzi_s[0];

				if(useHelio > 0){
					m_s[0] += mi;
				}
			}
			__syncthreads();
		
			double ax = 0.0;
			double ay = 0.0;
			double az = 0.0;

			if(useHelio == 0){
				if(itx < Nperturbers){
					if(idi != id_s[itx]){
						accP(m_s[itx], x_s[itx], y_s[itx], z_s[itx], xi_s[0], yi_s[0], zi_s[0], ax, ay, az);
					}
				}
				__syncthreads();
			}
			else{
				if(itx < Nperturbers){
					if(idi != id_s[itx]){
						accP(m_s[itx], x_s[itx], y_s[itx], z_s[itx], xi_s[0], yi_s[0], zi_s[0], ax, ay, az);
	//printf("Nij %d %d %.20g %.20g %.20g\n", idx, j, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit);
					}
				}

	//printf("N0 %d %.20g %.20g %.20g\n", idx, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit);
				if(itx > 0 && itx < Nperturbers){
					if(idi != id_s[itx]){
						accP2(m_s[itx], x_s[itx], y_s[itx], z_s[itx], ax, ay, az);
	//printf("Npi %d %d %.20g %.20g %.20g %d\n", idx, j, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit, S);
					}
				}
			}
			__syncthreads();
			//reduce

			for(int i = 1; i < warpSize; i*=2){
			#if def_OldShuffle == 0
				ax += __shfl_xor_sync(0xffffffff, ax, i, warpSize);
				ay += __shfl_xor_sync(0xffffffff, ay, i, warpSize);
				az += __shfl_xor_sync(0xffffffff, az, i, warpSize);
			#else
				ax += __shfld_xor(ax, i);
				ay += __shfld_xor(ay, i);
				az += __shfld_xor(az, i);
			#endif
			//printf("s reduce  %d %d %d %.20g\n", blockIdx.x, i, threadIdx.x, s);
			}

			__syncthreads();
	//if(itx == 0){
	//printf("Np %d %.20g %.20g %.20g %d\n", idx, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit, S);
	//}

			if(itx == 0){

				if(useGR == 2){
					acchGR2(xi_s[0], yi_s[0], zi_s[0], vxi_s[0], vyi_s[0], vzi_s[0], ax, ay, az);
				}

				if(useNonGrav == 1){
					NonGrav(xi_s[0], yi_s[0], zi_s[0], vxi_s[0], vyi_s[0], vzi_s[0], ax, ay, az, A1i, A2i, A3i);
				}

				if(useJ2 == 1){
					J2(m_s[3], x_s[3], y_s[3], z_s[3], xi_s[0], yi_s[0], zi_s[0], ax, ay, az);
				}

				kvx_d[idx + S * N] = ax;
				kvy_d[idx + S * N] = ay;
				kvz_d[idx + S * N] = az;

			}
			__syncthreads();
		}

		if(itx == 0){
			//update

			double dx = 0.0;
			double dy = 0.0;
			double dz = 0.0;
			double dvx = 0.0;
			double dvy = 0.0;
			double dvz = 0.0;
			
			for(int S = 0; S < RKFn; ++S){
				double dtb = dt * b_c[S];
				dx += dtb * kx_d[idx + S * N];
				dy += dtb * ky_d[idx + S * N];
				dz += dtb * kz_d[idx + S * N];

				dvx += dtb * kvx_d[idx + S * N];
				dvy += dtb * kvy_d[idx + S * N];
				dvz += dtb * kvz_d[idx + S * N];
			}
			
			dx_d[idx] = dx;
			dy_d[idx] = dy;
			dz_d[idx] = dz;
			dvx_d[idx] = dvx;
			dvy_d[idx] = dvy;
			dvz_d[idx] = dvz;


			//compute error

			double ym = fabs(xi0);
			ym = fmax(ym, fabs(yi0));
			ym = fmax(ym, fabs(zi0));
			ym = fmax(ym, fabs(vxi0));
			ym = fmax(ym, fabs(vyi0));
			ym = fmax(ym, fabs(vzi0));

			ym = fmax(ym, fabs(xi0 + dx));
			ym = fmax(ym, fabs(yi0 + dy));
			ym = fmax(ym, fabs(zi0 + dz));
			ym = fmax(ym, fabs(vxi0 + dvx));
			ym = fmax(ym, fabs(vyi0 + dvy));
			ym = fmax(ym, fabs(vzi0 + dvz));


			double isc = 1.0 / (def_atol + ym * def_rtol);

			double s = 1.0e6;       //large number

			double errorkx = 0.0;
			double errorky = 0.0;
			double errorkz = 0.0;
			double errorkvx = 0.0;
			double errorkvy = 0.0;
			double errorkvz = 0.0;

			for(int S = 0; S < RKFn; ++S){
				//this is y1i - y^1i
				double f = (b_c[S] - bb_c[S]) * dt;
				errorkx += f * kx_d[S * N + idx];
				errorky += f * ky_d[S * N + idx];
				errorkz += f * kz_d[S * N + idx];

				errorkvx += f * kvx_d[S * N + idx];
				errorkvy += f * kvy_d[S * N + idx];
				errorkvz += f * kvz_d[S * N + idx];
			}
			
			double errork = 0.0;
			errork += errorkx * errorkx * isc * isc;
			errork += errorky * errorky * isc * isc;
			errork += errorkz * errorkz * isc * isc;
			errork += errorkvx * errorkvx * isc * isc;
			errork += errorkvy * errorkvy * isc * isc;
			errork += errorkvz * errorkvz * isc * isc;

			errork = sqrt(errork / 6.0);	//6 is the number of dimensions

			s = pow( 1.0  / errork, ee);
			s = fmax(def_facmin, def_fac * s);
			s = fmin(def_facmax, s);

			// ****************************
			if(snew_d[idx].y >= 1.0){
				snew_d[idx].x = s;
			}
			else{
				snew_d[idx].x = 1.5;
			}
			if(s * dti < dtiMin){
				snew_d[idx].y = fmin(snew_d[idx].y, s);
			}
//printf("snew %d %.20g %.20g %.20g\n", idx, snew_d[idx].y, snew_d[idx].x, errork);
			// ****************************

		}
	}
}


void update(double *x_h, double *y_h, double *z_h, double *vx_h, double *vy_h, double *vz_h, double *dx_h, double *dy_h, double *dz_h, double *dvx_h, double *dvy_h, double *dvz_h, int i){

	x_h[i] += dx_h[i];
	y_h[i] += dy_h[i];
	z_h[i] += dz_h[i];

	vx_h[i] += dvx_h[i];
	vy_h[i] += dvy_h[i];
	vz_h[i] += dvz_h[i];

//printf("dx %d %.20g %.20g %.20g %.20g %.20g %.20g\n", i, dx_h[i], dy_h[i], dz_h[i], dvx_h[i], dvy_h[i], dvz_h[i]);
//printf("dx %d %.20g %.20g %.20g %.20g %.20g %.20g\n", i, x_h[i], y_h[i], z_h[i], vx_h[i], vy_h[i], vz_h[i]);
}

__global__ void update_kernel(double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *dx_d, double *dy_d, double *dz_d, double *dvx_d, double *dvy_d, double *dvz_d, int N, int Nperturbers){

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if(idx >= Nperturbers && idx < N){
		x_d[idx] += dx_d[idx];
		y_d[idx] += dy_d[idx];
		z_d[idx] += dz_d[idx];

		vx_d[idx] += dvx_d[idx];
		vy_d[idx] += dvy_d[idx];
		vz_d[idx] += dvz_d[idx];
	}
}




__global__ void bufferToX_kernel(double *XYdata_d, double *timep_d, double *xp_d, double *yp_d, double *zp_d, int N){

	int id = threadIdx.x + blockDim.x * blockIdx.x;

	if(id < N){
		timep_d[id] = XYdata_d[id * 4];
		xp_d[id] = XYdata_d[id * 4 + 1];
		yp_d[id] = XYdata_d[id * 4 + 2];
		zp_d[id] = XYdata_d[id * 4 + 3];
if(id < 30 || id > N - 10 ) printf("buffer %d %d %.20g %.20g %.20g %.20g\n", id / 27, id, timep_d[id], xp_d[id], yp_d[id], zp_d[id]);
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

int main(int argc, char*argv[]){

	//Number of planets
	int NTP = 1;			//number of small particles
	const int Nperturbers = 27;
	const int Ninterpolate = 16;	//number of interpolation points
	const double dtimep = 1.0;	//interval between stored time steps
	double time0 = 0.0;		//start time from simulation
	double time1 = 0.0;		//end time from simulation
	double outStart = 0;		//start time of output files

	int comet = 0;

	int useGR = 2;		//2
	//2: Sitarski 1982, heliocentric coordinates

	int useJ2 = 1;		//1
	int useNonGrav = 1;	//1

	int useFIFO = 2;	//use 0 or 2
	int useGPU = 1;		// 0 or 1

	FILE *binfile;
	if(useFIFO == 2){	
		//binfile = fopen("210801_2342_genga_de440_perturbers.bin", "rb");
		//binfile = fopen("210921_2148_genga_in_yarkovsky_elements.bin", "rb");
		binfile = fopen("211208_1916_genga_in_2021-12-08_specific_desig.bin", "rb");
		//binfile = fopen("210801_2104_genga_in_GA.bin", "rb");
		//binfile = fopen("210705_2315_genga_req.bin", "rb");
		//binfile = fopen("220301_2048_genga_in_new_last_14_days.bin", "rb");
	}	

	int useHelio = 1;
	int outHelio = 1;
	//1 print output in heliocentric coordinates
	//0 print output in barycentric coordinates

	int outBinary = 0;
	//0 print output files in ASCII format	
	//1 print output files in binary format	

	//long long int Nsteps = 40000;
	long long int Nsteps = 1750000;
	long long int outInterval = 1;

	double dti = 0.05;
	double dt = dti * dayUnit;
	double dtiMin = 1.0e-5;	//Minimal time step

	double dts = 0.01;	//interval of interpolated points

	double InVersion = 0.0;	//version of Input file format
	int DoPreIntegration = 0; //If this is 1 then do a pre-Integration to synchronize all initial conditions

	//Integrator
	//const int RKFn = 6; //RKF45
	//const int RKFn = 7; //DP54
	const int RKFn = 13; //RKF78

	double timing[6];
	for(int i = 0; i < 6; ++i){
		timing[i] = 0.0;
	}
	cudaEvent_t tt1;                        //start time
	cudaEvent_t tt2;                        //end time
	cudaEventCreate(&tt1);
	cudaEventCreate(&tt2);

	cudaEventRecord(tt1);

	time1 = time0 + dt * Nsteps;
	unsigned long long int outI = 1llu;

	if(useFIFO == 2){
		printf("read file\n");
		readHeader(binfile, time0, time1, outInterval, outStart, NTP, comet, InVersion);

		//change this later to be more general
		Nsteps = 1e9;
		dti = 10.0;
		dts = 0.1;
		dtiMin = 1.0e-5;
		//dti = 5.0;
		//dts = 0.1;
		//dtiMin = 5.0;
		//NTP = 1;

		outInterval = 10.0;
		//outInterval = 1e4;
		//outStart = 2450800.5;


		outI = (outInterval + 0.5 * dts) / dts;

	//	if(fabs(outStart - time0) > fabs(outInterval)){
	//		outI = (outStart - time0 + 0.5 * dts) / dts;
	//	}


printf("outStart: %.20g, time0: %.20g, dts: %g, outI: %llu,  outInterval: %lld\n", outStart, time0, dts, outI, outInterval);

	}
	dt = dti * dayUnit;

	for(int i = 1; i < argc; i += 2){

		if(strcmp(argv[i], "-Nsteps") == 0){
			Nsteps = atoll(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-outInterval") == 0){
			outInterval = atoll(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-dt") == 0){
			dti = atof(argv[i + 1]);
			dt = dti * dayUnit;
		}
		else if(strcmp(argv[i], "-N") == 0){
			NTP = atoi(argv[i + 1]);
		}
	}
	int NN = Nperturbers + NTP;//+ 8192 * 4; //28
	printf("NN %d %d %d\n", Nperturbers, NTP, NN);

	float milliseconds = 0.0f;

	const char *myfifo = "myfifo";
	const char *fifoCheck = "fifoCheck";
	if(useFIFO == 1){
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

	cudaError_t error;

	unsigned long long int *id_h, *id_d;
	double *m_h, *m_d;
	double *x_h, *x_d;
	double *y_h, *y_d;
	double *z_h, *z_d;
	double *vx_h, *vx_d;
	double *vy_h, *vy_d;
	double *vz_h, *vz_d;
	double *A1_h, *A1_d;
	double *A2_h, *A2_d;
	double *A3_h, *A3_d;
	double *jd_init_h, *jd_init_d;

	double *dx_h, *dx_d;
	double *dy_h, *dy_d;
	double *dz_h, *dz_d;
	double *dvx_h, *dvx_d;
	double *dvy_h, *dvy_d;
	double *dvz_h, *dvz_d;


	//coordinates from data table
	double *timep_h, *timep_d;
	double *xp_h, *xp_d;
	double *yp_h, *yp_d;
	double *zp_h, *zp_d;

	//Backup for repeated integrations
	double *x0_h;
	double *y0_h;
	double *z0_h;
	double *vx0_h;
	double *vy0_h;
	double *vz0_h;

	const int NTable = 65000;	//length of perturbers table, number of days
	//store the entire perturbers file

	//allocate data on host
	id_h = (unsigned long long int*)malloc(NN * sizeof(unsigned long long int));
	m_h = (double*)malloc(NN * sizeof(double));
	x_h = (double*)malloc(NN * sizeof(double));
	y_h = (double*)malloc(NN * sizeof(double));
	z_h = (double*)malloc(NN * sizeof(double));
	vx_h = (double*)malloc(NN * sizeof(double));
	vy_h = (double*)malloc(NN * sizeof(double));
	vz_h = (double*)malloc(NN * sizeof(double));
	A1_h = (double*)malloc(NN * sizeof(double));
	A2_h = (double*)malloc(NN * sizeof(double));
	A3_h = (double*)malloc(NN * sizeof(double));
	jd_init_h = (double*)malloc(NN * sizeof(double));

	timep_h = (double*)malloc(Nperturbers * NTable * sizeof(double));
	xp_h = (double*)malloc(Nperturbers * NTable * sizeof(double));
	yp_h = (double*)malloc(Nperturbers * NTable * sizeof(double));
	zp_h = (double*)malloc(Nperturbers * NTable * sizeof(double));

	x0_h = (double*)malloc(NN * sizeof(double));
	y0_h = (double*)malloc(NN * sizeof(double));
	z0_h = (double*)malloc(NN * sizeof(double));
	vx0_h = (double*)malloc(NN * sizeof(double));
	vy0_h = (double*)malloc(NN * sizeof(double));
	vz0_h = (double*)malloc(NN * sizeof(double));

	dx_h = (double*)malloc(NN * sizeof(double));
	dy_h = (double*)malloc(NN * sizeof(double));
	dz_h = (double*)malloc(NN * sizeof(double));
	dvx_h = (double*)malloc(NN * sizeof(double));
	dvy_h = (double*)malloc(NN * sizeof(double));
	dvz_h = (double*)malloc(NN * sizeof(double));

	//allocate data on the device
	cudaMalloc((void **) &id_d, NN * sizeof(unsigned long long int));
	cudaMalloc((void **) &m_d, NN * sizeof(double));
	cudaMalloc((void **) &x_d, NN * sizeof(double));
	cudaMalloc((void **) &y_d, NN * sizeof(double));
	cudaMalloc((void **) &z_d, NN * sizeof(double));
	cudaMalloc((void **) &vx_d, NN * sizeof(double));
	cudaMalloc((void **) &vy_d, NN * sizeof(double));
	cudaMalloc((void **) &vz_d, NN * sizeof(double));
	cudaMalloc((void **) &A1_d, NN * sizeof(double));
	cudaMalloc((void **) &A2_d, NN * sizeof(double));
	cudaMalloc((void **) &A3_d, NN * sizeof(double));
	cudaMalloc((void **) &jd_init_d, NN * sizeof(double));

	cudaMalloc((void **) &timep_d, Nperturbers * NTable * sizeof(double));
	cudaMalloc((void **) &xp_d, Nperturbers * NTable * sizeof(double));
	cudaMalloc((void **) &yp_d, Nperturbers * NTable * sizeof(double));
	cudaMalloc((void **) &zp_d, Nperturbers * NTable * sizeof(double));

	cudaMalloc((void **) &dx_d, NN * sizeof(double));
	cudaMalloc((void **) &dy_d, NN * sizeof(double));
	cudaMalloc((void **) &dz_d, NN * sizeof(double));
	cudaMalloc((void **) &dvx_d, NN * sizeof(double));
	cudaMalloc((void **) &dvy_d, NN * sizeof(double));
	cudaMalloc((void **) &dvz_d, NN * sizeof(double));
	

	for(int i = 0; i < NN; ++i){
		A1_h[i] = 0.0;
		A2_h[i] = 0.0;
		A3_h[i] = 0.0;
		jd_init_h[i] = 0.0;
	}
	
	//Sun
	id_h[0] = 10;
	m_h[0] = 1.0;
	x_h[0] = 0.0;
	y_h[0] = 0.0;
	z_h[0] = 0.0;
	vx_h[0] = 0.0;
	vy_h[0] = 0.0;
	vz_h[0] = 0.0;

	for(int i = 1; i < NN; ++i){
		id_h[i] = Nperturbers + i;
		m_h[i] = 0.0;
		x_h[i] = 0.0;
		y_h[i] = 0.0;
		z_h[i] = 0.0;
		vx_h[i] = 0.0;
		vy_h[i] = 0.0;
		vz_h[i] = 0.0;
	}

	FILE *outfile;
	char outfilename[160];	
	FILE *dtfile;
	char dtfilename[160];	
	sprintf(dtfilename, "timesteps.dat");


	double time = 0.0;


	// ********************************
	// Read perturbers masses from perturbers.h file
	perturbersMass(m_h, Nperturbers);
	perturbersIDs(id_h, Nperturbers);
	// ********************************

	//Erase Outbinary file
	if(outBinary == 1){
		if(outHelio == 1){
			sprintf(outfilename, "Outhelio10.bin");
		}
		else{
			sprintf(outfilename, "Outbary10.bin");
		}
		outfile = fopen(outfilename, "wb");
		fclose(outfile);
	}


//m[Nperturbers] = 1.e-11; //ca mass of Flora

	int N = Nperturbers;
	printf("Read initial conditions\n");

	if(useFIFO == 2){	
		//read test particles
		int er = 0;
		er = readFile(binfile, Nperturbers, x_h, y_h, z_h, vx_h, vy_h, vz_h, A1_h, A2_h, A3_h, id_h, jd_init_h, NTP, time0, InVersion, DoPreIntegration);
		printf("read file OK\n");
		fclose(binfile);

							
		// -----------------------------------
		// Use this to extract a single object
		int ii = 29;//166;//29; //84;
		//int ii = 83;//166;//29; //84;
		id_h[N] = id_h[ii];
		x_h[N] = x_h[ii];
		y_h[N] = y_h[ii];
		z_h[N] = z_h[ii];
		vx_h[N] = vx_h[ii];
		vy_h[N] = vy_h[ii];
		vz_h[N] = vz_h[ii];
		A1_h[N] = A1_h[ii];
		A2_h[N] = A2_h[ii];
		A3_h[N] = A3_h[ii];
		jd_init_h[N] = jd_init_h[ii];
		NTP = 1;
		
printf("xyz %.40g %.40g %.40g %.40g %.40g %.40g %.40g %.40g %.40g %.20g %llu\n", x_h[N], y_h[N], z_h[N], vx_h[N], vy_h[N], vz_h[N], A1_h[N], A2_h[N], A3_h[N], jd_init_h[N], id_h[N]);
//x_h[N] += 1.0e-6;
		
		// -----------------------------------
		

		N += NTP;
		if(er == 1) return 0;
	}
	else{

		FILE *infile;
		char infilename[160];

		sprintf(infilename, "initial.dat");
		infile = fopen(infilename, "r");
		for(int i = Nperturbers; i < NN; ++i){
			int er = 0;
			fscanf(infile, "%lf", &time);
			fscanf(infile, "%llu", &id_h[i]);
			fscanf(infile, "%lf", &x_h[i]);
			fscanf(infile, "%lf", &y_h[i]);
			fscanf(infile, "%lf", &z_h[i]);
			fscanf(infile, "%lf", &vx_h[i]);
			fscanf(infile, "%lf", &vy_h[i]);
			er = fscanf(infile, "%lf", &vz_h[i]);
			er = fscanf(infile, "%lf", &A1_h[i]);
			er = fscanf(infile, "%lf", &A2_h[i]);
			er = fscanf(infile, "%lf", &A3_h[i]);
			//fscanf(infile, "%lf", &ALN[i]);
			//fscanf(infile, "%lf", &NK[i]);
			//fscanf(infile, "%lf", &NM[i]);
			//fscanf(infile, "%lf", &Nn[i]);
			//er = fscanf(infile, "%lf", &R0[i]);
			if(er < 0) break;
			++N;
//printf("xyz %.40g %.40g %.40g %.40g %.40g %.40g %.40g %.40g %.40g\n", x_h[i], y_h[i], z_h[i], vx_h[i], vy_h[i], vz_h[i], A1_h[i], A2_h[i], A3_h[i]);
//printf("er %d %llu %d %d %.20g %.20g %.20g\n", i, id_h[i], er, N, x_h[i], y_h[i], z_h[i]);
		}
		fclose(infile);
		time0 = time;	//start time from simulation

time1 = 2461000.5;
	}

	double timep0 = 0.0;	//start time from perturbers file

	for(int i = Nperturbers; i < N; ++i){
		vx_h[i] /= dayUnit;
		vy_h[i] /= dayUnit;
		vz_h[i] /= dayUnit;

		A1_h[i] /= (dayUnit * dayUnit);
		A2_h[i] /= (dayUnit * dayUnit);
		A3_h[i] /= (dayUnit * dayUnit);
	}

        cudaEventRecord(tt2);
        cudaEventSynchronize(tt2);
        cudaEventElapsedTime(&milliseconds, tt1, tt2);
	printf("Time for ic and allocation, %g seconds\n", milliseconds * 0.001);
        timing[0] += milliseconds * 0.001;

	//###########################################
	//perturbers table
	//###########################################
	cudaEventRecord(tt1);

	FILE *XVfile;
	if(useHelio == 1){
		XVfile = fopen("All_h.bin", "rb");
	}
	else{
		XVfile = fopen("All_b.bin", "rb");
	}

	if(XVfile == NULL){
		printf("Error, perturber file not found\n");
		return 0;
	}

	// -----------------------------------------
	//Read table

	double *readBufferA_h;
	double *readBufferB_h;
	//the buffer contains time, x, y, ,z from all perturbers 
	//the buffer has 2 swaps

	double *XYdata_h, *XYdata_d;

	if(useGPU == 0){
		XYdata_h = (double*)malloc(Nperturbers * NTable * 4 * sizeof(double));
		readBufferA_h = (double*)malloc(Nperturbers * 4 * sizeof(double));
		readBufferB_h = (double*)malloc(Nperturbers * 4 * sizeof(double));
	}
	else{
		cudaHostAlloc((void **) &readBufferA_h, Nperturbers * 4 * sizeof(double), cudaHostAllocDefault);
		cudaHostAlloc((void **) &readBufferB_h, Nperturbers * 4 * sizeof(double), cudaHostAllocDefault);
		cudaMalloc((void **) &XYdata_d, Nperturbers * NTable * 4 * sizeof(double));
	}

	int NTableC = 0;
	for(int t = 0; t < 1000000; ++t){
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
				if(t % 2 == 0) printf("XYa %d %.20g %g %g %g\n", i, readBufferA_h[i * 4 + 0], readBufferA_h[i * 4 + 1], readBufferA_h[i * 4 + 2], readBufferA_h[i * 4 + 3]);
				if(t % 2 == 1) printf("XYb %d %.20g %g %g %g\n", i, readBufferB_h[i * 4 + 0], readBufferB_h[i * 4 + 1], readBufferB_h[i * 4 + 2], readBufferB_h[i * 4 + 3]);
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

	printf("NTableC: %d\n", NTableC);

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

	
	//###########################################
	// end perturbers table
	//###########################################
	

	//copy the data to the device
	if(useGPU > 0){
		cudaMemcpy(m_d, m_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(id_d, id_h, N * sizeof(unsigned long long int), cudaMemcpyHostToDevice);
		cudaMemcpy(x_d, x_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(y_d, y_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(z_d, z_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(vx_d, vx_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(vy_d, vy_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(vz_d, vz_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(A1_d, A1_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(A2_d, A2_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(A3_d, A3_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(jd_init_d, jd_init_h, N * sizeof(double), cudaMemcpyHostToDevice);
	}
	

	//remove this later to shared memory
	double *kx_h, *kx_d;
	double *ky_h, *ky_d;
	double *kz_h, *kz_d;
	double *kvx_h, *kvx_d;
	double *kvy_h, *kvy_d;
	double *kvz_h, *kvz_d;

	double2 *snew_h, *snew_d;
	double *dtmin_h;

	kx_h = (double*)malloc(N * RKFn * sizeof(double));
	ky_h = (double*)malloc(N * RKFn * sizeof(double));
	kz_h = (double*)malloc(N * RKFn * sizeof(double));
	kvx_h = (double*)malloc(N * RKFn * sizeof(double));
	kvy_h = (double*)malloc(N * RKFn * sizeof(double));
	kvz_h = (double*)malloc(N * RKFn * sizeof(double));

	snew_h = (double2*)malloc(N * sizeof(double2));
	dtmin_h = (double*)malloc(N * sizeof(double));

	for(int i = 0; i < N; ++i){
		dtmin_h[i] = 1.0e6;
		snew_h[i].x = 1.5;
		snew_h[i].y = 1.5;
	}


	//interpolation table
	double *xTable_h, *xTable_d;
	double *yTable_h, *yTable_d;
	double *zTable_h, *zTable_d;

	xTable_h = (double*)malloc(Nperturbers * RKFn * sizeof(double));
	yTable_h = (double*)malloc(Nperturbers * RKFn * sizeof(double));
	zTable_h = (double*)malloc(Nperturbers * RKFn * sizeof(double));

	if(useGPU > 0){
		cudaMalloc((void **) &kx_d, N * RKFn * sizeof(double));
		cudaMalloc((void **) &ky_d, N * RKFn * sizeof(double));
		cudaMalloc((void **) &kz_d, N * RKFn * sizeof(double));
		cudaMalloc((void **) &kvx_d, N * RKFn * sizeof(double));
		cudaMalloc((void **) &kvy_d, N * RKFn * sizeof(double));
		cudaMalloc((void **) &kvz_d, N * RKFn * sizeof(double));

		cudaMalloc((void **) &snew_d, N * sizeof(double2));

		cudaMalloc((void **) &xTable_d, Nperturbers * RKFn * sizeof(double));
		cudaMalloc((void **) &yTable_d, Nperturbers * RKFn * sizeof(double));
		cudaMalloc((void **) &zTable_d, Nperturbers * RKFn * sizeof(double));
	}
	


	// *******************************************************************
	// Allocate and set parameters for the Runge-Kutta-Fehlberg integrator
	// *******************************************************************
	double *a_h, *b_h, *bb_h, *c_h;

	a_h = (double*)malloc(RKFn * RKFn * sizeof(double));
	b_h = (double*)malloc(RKFn * sizeof(double));
	bb_h = (double*)malloc(RKFn * sizeof(double));
	c_h = (double*)malloc(RKFn * sizeof(double));

	for(int i = 0; i < RKFn; ++i){
		for(int j = 0; j < RKFn; ++j){
			a_h[i * RKFn + j] = 0.0;
		}
		b_h[i] = 0.0;
		bb_h[i] = 0.0;
		c_h[i] = 0.0;
	}

	double ee;

	if(RKFn == 6){
		setRKF45(a_h, b_h, bb_h, c_h);
		ee = 1.0 / 5.0;
	}
	else if(RKFn == 7){
		setDP54(a_h, b_h, bb_h, c_h);
		ee = 1.0 / 5.0;
	}
	else if(RKFn == 13){
		setRKF78(a_h, b_h, bb_h, c_h);
		ee = 1.0 / 8.0;
	}
	else{
		printf("RKFn values not valid %d\n", RKFn);
		return 0;
	}

	if(useGPU > 0){
		cudaMemcpyToSymbol(a_c, a_h, RKFn * RKFn * sizeof(double), 0, cudaMemcpyHostToDevice);
		cudaMemcpyToSymbol(b_c, b_h, RKFn * sizeof(double), 0, cudaMemcpyHostToDevice);
		cudaMemcpyToSymbol(bb_c, bb_h, RKFn * sizeof(double), 0, cudaMemcpyHostToDevice);
		cudaMemcpyToSymbol(c_c, c_h, RKFn * sizeof(double), 0, cudaMemcpyHostToDevice);
	}

	cudaDeviceSynchronize();
	error = cudaGetLastError();
	printf("copy error = %d = %s\n",error, cudaGetErrorString(error));
	if(error != 0.0){
		return 0;
	}

	cudaEventRecord(tt2);
	cudaEventSynchronize(tt2);
	
	cudaEventElapsedTime(&milliseconds, tt1, tt2);
	printf("Time for perturbers table, %g seconds\n", milliseconds * 0.001);
        timing[1] += milliseconds * 0.001;

	cudaEventRecord(tt1);

	error = cudaGetLastError();
	printf("Perturbers error = %d = %s\n",error, cudaGetErrorString(error));
	if(error != 0.0){
		return 0;
	}
	

	double dtiOld = dti;
	time = time0;
	printf("dtiOld %g\n", dtiOld);

	dtfile = fopen(dtfilename, "w");

	//###########################################
	// Start pre-integration
	//###########################################

	if(DoPreIntegration == 1){
		if(useGPU > 0){
			//copy perturbers data to CPU, because its not yet here
			cudaMemcpy(timep_h, timep_d, Nperturbers * NTable * sizeof(double), cudaMemcpyDeviceToHost);
			cudaMemcpy(xp_h, xp_d, Nperturbers * NTable * sizeof(double), cudaMemcpyDeviceToHost);
			cudaMemcpy(yp_h, yp_d, Nperturbers * NTable * sizeof(double), cudaMemcpyDeviceToHost);
			cudaMemcpy(zp_h, zp_d, Nperturbers * NTable * sizeof(double), cudaMemcpyDeviceToHost);
		}


		//Do pre-Integration to synchronize all initial conditions
		printf("Start pre-Integration up to time %.20g\n", outStart);
		for(int i = Nperturbers; i < N; ++i){
			time = jd_init_h[i];
			double snew = 1.0;
			dti = 10.0;
			dt = dti * dayUnit;
			int stop = 0;

printf("preIntegration %d %.20g %.20g\n", i, time, outStart);

			if(time == outStart){
				//These particles are already at the correct epoch
				continue;
			}
			if(time > outStart){
				dti = -dti;
				dt = -dt;
			}
			
			for(long long int t = 1; t <= Nsteps; ++t){
//printf("A %d %lld %.20g %.20g %.20g %.20g\n", i, id_h[i], time, x_h[i], y_h[i], z_h[i]);

				interpolateTable(Ninterpolate, Nperturbers, RKFn, xp_h, yp_h, zp_h, timep_h, timep0, dtimep, time, dti, xTable_h, yTable_h, zTable_h, c_h);
				//for(int S = 0; S < RKFn; ++S){
					//for(int p = 0; p < Nperturbers; ++p){
						//interpolate(Ninterpolate, Nperturbers, xp_h, yp_h, zp_h, timep_h, timep0, dtimep, time + c_h[S] * dti, xt_h, yt_h, zt_h, p);
						//interpolate2(Ninterpolate, Nperturbers, xp_h, yp_h, zp_h, timep_h, timep0, dtimep, time + c_h[S] * dti, xt_h, yt_h, zt_h, p);
					//}
					//stageStep(id_h, m_h, xt_h, yt_h, zt_h, vxt_h, vyt_h, vzt_h, kx_h, ky_h, kz_h, kvx_h, kvy_h, kvz_h, A1_h, A2_h, A3_h, S, i, Nperturbers, N, useHelio, useGR, useJ2, useNonGrav);

				//}
			//	computeError1(snew_h, kx_h, ky_h, kz_h, kvx_h, kvy_h, kvz_h, b_h, bb_h, RKFn, i, N, snew, dt, ee);

				if(snew >= 1.0){		
					update(x_h, y_h, z_h, vx_h, vy_h, vz_h, dx_h, dy_h, dz_h, dvx_h, dvy_h, dvz_h, i);	
		
					time += dti;

					dti *= snew;
					dt = dti * dayUnit;

					if(stop == 1){
						break;
					}
					else{
						dtmin_h[i] = fmin(fabs(dt / dayUnit), dtmin_h[i]);
					}
				}
				else{
					dti *= snew;
					dt = dti * dayUnit;
					stop = 0;
				}

//printf("B %d %lld %.20g %.20g %.20g %.20g %.20g %.20g\n", i, id_h[i], time, dti, snew, x_h[i], y_h[i], z_h[i]);


				if(dti > 0.0 && time + dti > outStart){
					dti = (outStart - time);
					dt = dti * dayUnit;
					stop = 1;
//printf("Final time step %.20g\n", dti);
				}
				if(dti < 0.0 && time + dti < outStart){
					dti = (outStart - time);
					dt = dti * dayUnit;
					stop = 1;
//printf("Final time step %.20g\n", dti);
				}

			}
printf("C %d %lld %.20g %.20g %.20g %.20g %.20g %.20g\n", i, id_h[i], time, dti, snew, x_h[i], y_h[i], z_h[i]);
			if(fabs(time - outStart) > 1.0e-10){
				printf("Error in pre-integration of particle %d, start time not reached\n", i);
				return 0;
			}
		}
		dti = dtiOld;
		dt = dti * dayUnit;
		time = outStart;
		if(useGPU > 0){
			cudaMemcpy(x_d, x_h, N * sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(y_d, y_h, N * sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(z_d, z_h, N * sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(vx_d, vx_h, N * sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(vy_d, vy_h, N * sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(vz_d, vz_h, N * sizeof(double), cudaMemcpyHostToDevice);
		}
	}
	//###########################################
	// End pre-integration
	//###########################################

	//save coordinates for repeated integrations
	for(int i = 0; i < N; ++i){
		x0_h[i] = x_h[i];
		y0_h[i] = y_h[i];
		z0_h[i] = z_h[i];
		vx0_h[i] = vx_h[i];
		vy0_h[i] = vy_h[i];
		vz0_h[i] = vz_h[i];
	}

	cudaEventRecord(tt2);
	cudaEventSynchronize(tt2);
	
	cudaEventElapsedTime(&milliseconds, tt1, tt2);
	printf("Time for pre-integration, %g seconds\n", milliseconds * 0.001);
        timing[2] += milliseconds * 0.001;

	cudaEventRecord(tt1);

	//###########################################
	// First output
	//###########################################

	output(snew_h, dtmin_h, x_h, y_h, z_h, vx_h, vy_h, vz_h, xTable_h, yTable_h, zTable_h, m_h, id_h, 0, time, N, Nperturbers, useGPU, useHelio, outHelio, outBinary, 0);

	for(int i = 0; i < N; ++i){
		dtmin_h[i] = 1.0e6;
		snew_h[i].x = 1.5;
		snew_h[i].y = 1.5;
		if(useGPU > 0){
			cudaMemcpy(snew_d, snew_h, N * sizeof(double2), cudaMemcpyHostToDevice);
		}
	}
	printf("dti %g\n", dti);

	//###########################################
	// Start time step loop
	//###########################################


	unsigned long long int cOut = 0;		//counter for output
	int ci = 0;

	int Sn[4];
	Sn[0] = N;
for(int S = 0; S < 1; ++S){

	for(long long int t = 1; t <= Nsteps; ++t){
	//for(long long int t = 1; t < 5000; ++t){
	//for(long long int t = 1; t < 10; ++t){

//cudaDeviceSynchronize();
//printf("%lld %d | %d %d\n", t, 0, NTP, Nperturbers);	

		double snew = 1.0;

		if(useGPU == 0){

			interpolateTable(Ninterpolate, Nperturbers, RKFn, xp_h, yp_h, zp_h, timep_h, timep0, dtimep, time, dti, xTable_h, yTable_h, zTable_h, c_h);
			stageStep(id_h, m_h, x_h, y_h, z_h, vx_h, vy_h, vz_h, dx_h, dy_h, dz_h, dvx_h, dvy_h, dvz_h, xTable_h, yTable_h, zTable_h, kx_h, ky_h, kz_h, kvx_h, kvy_h, kvz_h, A1_h, A2_h, A3_h, snew_h, a_h, b_h, bb_h, dt, dti, dtiMin, RKFn, Nperturbers, N, useHelio, useGR, useJ2, useNonGrav, ee, snew);

			//for(int S = 0; S < RKFn; ++S){
				//for(int p = 0; p < Nperturbers; ++p){
					//interpolate(Ninterpolate, Nperturbers, xp_h, yp_h, zp_h, timep_h, timep0, dtimep, time + c_h[S] * dti, xt_h, yt_h, zt_h, p);
					//interpolate2(Ninterpolate, Nperturbers, xp_h, yp_h, zp_h, timep_h, timep0, dtimep, time + c_h[S] * dti, xt_h, yt_h, zt_h, p);
				//}
				//for(int i = Nperturbers; i < N; ++i){
				//	stageStep(id_h, m_h, xt_h, yt_h, zt_h, vxt_h, vyt_h, vzt_h, kx_h, ky_h, kz_h, kvx_h, kvy_h, kvz_h, A1_h, A2_h, A3_h, S, i, Nperturbers, N, useHelio, useGR, useJ2, useNonGrav);
				//}


			//}
		// -----------------------------------
		}
		else{
			if(useGPU == 1){				
				interpolateTable_kernel < Ninterpolate > <<< dim3(Nperturbers, RKFn, 1), dim3(Ninterpolate) >>> (Nperturbers, RKFn, xp_d, yp_d, zp_d, timep_d, timep0, dtimep, time, dti, xTable_d, yTable_d, zTable_d);
				if(N > 300){
					stageStep1_kernel < Nperturbers > <<< (N + 127) / 128, 128 >>> (id_d, m_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, dx_d, dy_d, dz_d, dvx_d, dvy_d, dvz_d, xTable_d, yTable_d, zTable_d, kx_d, ky_d, kz_d, kvx_d, kvy_d, kvz_d, A1_d, A2_d, A3_d, snew_d, dt, dti, dtiMin, RKFn, N, useHelio, useGR, useJ2, useNonGrav, ee);
				}
				else{
					stageStep2_kernel < Nperturbers > <<< (N - Nperturbers), 32 >>> (id_d, m_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, dx_d, dy_d, dz_d, dvx_d, dvy_d, dvz_d, xTable_d, yTable_d, zTable_d, kx_d, ky_d, kz_d, kvx_d, kvy_d, kvz_d, A1_d, A2_d, A3_d, snew_d, dt, dti, dtiMin, RKFn, N, useHelio, useGR, useJ2, useNonGrav, ee);
				}

			}

			int nct = 512;
			int ncb = min((N + nct - 1) / nct, 1024);
			int WarpSize = 32;
			computeError_d1_kernel <<< ncb, nct, WarpSize * sizeof(double)  >>> (snew_d, Nperturbers, N);
			if(ncb > 1){
				computeError_d2_kernel <<< 1, ((ncb + WarpSize - 1) / WarpSize) * WarpSize, WarpSize * sizeof(double)  >>> (snew_d, ncb);
			}
			cudaMemcpy(snew_h, snew_d, sizeof(double2), cudaMemcpyDeviceToHost);
			cudaDeviceSynchronize();
			snew = snew_h[0].x;
		}
//if(useGPU > 0){
//cudaMemcpy(x_h, x_d, N * sizeof(double), cudaMemcpyDeviceToHost);
//cudaMemcpy(y_h, y_d, N * sizeof(double), cudaMemcpyDeviceToHost);
//cudaMemcpy(z_h, z_d, N * sizeof(double), cudaMemcpyDeviceToHost);
//}
printf("%.20g %llu dt: %.20g %.g %g %d\n", time, cOut, dti, dts, snew, S);
fprintf(dtfile, "%.20g %llu dt: %.20g %.g %g %d\n", time, cOut, dti, dts, snew, S);


//snew = 1.0;

		if(snew >= 1.0){
printf("---- accept step %llu %llu %.20g %g----\n", cOut, outI, time + dti, dti);		
			if(useGPU == 0){
				for(int i = Nperturbers; i < N; ++i){
					update(x_h, y_h, z_h, vx_h, vy_h, vz_h, dx_h, dy_h, dz_h, dvx_h, dvy_h, dvz_h, i);	
				}
			}
			else{
				update_kernel <<< (N + 127) / 128, 128 >>> (x_d, y_d, z_d, vx_d, vy_d, vz_d, dx_d, dy_d, dz_d, dvx_d, dvy_d, dvz_d, N, Nperturbers);	
			}
			ci = (dti + 0.5 * dts) / dts;
			cOut += ci;		

			time += dti;

			if(cOut >= outI){
			//if(t % 10 == 0){
				if(useGPU > 0){
					cudaMemcpy(snew_h, snew_d, N * sizeof(double2), cudaMemcpyDeviceToHost);
					cudaMemcpy(x_h, x_d, N * sizeof(double), cudaMemcpyDeviceToHost);
					cudaMemcpy(y_h, y_d, N * sizeof(double), cudaMemcpyDeviceToHost);
					cudaMemcpy(z_h, z_d, N * sizeof(double), cudaMemcpyDeviceToHost);
					cudaMemcpy(vx_h, vx_d, N * sizeof(double), cudaMemcpyDeviceToHost);
					cudaMemcpy(vy_h, vy_d, N * sizeof(double), cudaMemcpyDeviceToHost);
					cudaMemcpy(vz_h, vz_d, N * sizeof(double), cudaMemcpyDeviceToHost);
				}
				output(snew_h, dtmin_h, x_h, y_h, z_h, vx_h, vy_h, vz_h, xTable_h, yTable_h, zTable_h, m_h, id_h, t, time, N, Nperturbers, useGPU, useHelio, outHelio, outBinary, S);
				dti = dtiOld;
				dt = dti * dayUnit;
				snew = 1.0;
				cOut = 0;
				outI = (outInterval + 0.5 * dts) / dts; //needed only at the first time
			}
		}
		else{
printf(" ---- repeat step %llu %llu %.20g ----\n", cOut, outI, time);		
		
		}


		dti *= snew;

		if(useGPU < 2 && dtiOld >= 65.0 * dts && cOut % 10 == 0 && snew >= 1.0){
			//synchronize with coarser time steps
			dts *= 10;
			outI /= 10;
			cOut /= 10;
			dti = 5.0 * dts;
printf("increase time step C %g %g\n", dti, dts);
		}

		//round dti to dts intervals
		int dtt = dti / dts;
printf("dta %.20g %d %g %llu\n", dti, dtt, dts, cOut);
		dti = dtt * dts;
		if(dti < dts) dti = dts;


printf("dtb %.20g %d %g %llu\n", dti, dtt, dts, cOut);
		ci = (dti + 0.5 * dts) / dts;

if(dti < dtiMin) dti = dtiMin;


		dtiOld = dti;

		dt = dti * dayUnit;
//printf("%llu %llu %.20g, %.20g %.20g\n", cOut + ci, outI, time, time + dti, outStart);

		if(cOut + ci > outI && time + dti >= outStart){
			dti = (outI - cOut) * dts;

			dt = dti * dayUnit;
printf("   correct %.20g %.20g %.20g %.20g %llu %llu\n", time, time + dti, dti, dtiOld, cOut, outI);
printf("dtc %.20g %g\n", dti, dts);
		}


		else{

			if(useGPU < 2 && dti >= 65.0 * dts){
				//synchronize with coarser time steps
printf("increase time step A %g %g\n", dti, dts);
				dti = (((cOut + ci) / 10) * 10 - cOut) * dts;
				dt = dti * dayUnit;
				//dts *= 10;
				//outI /= 10;
				//cOut /= 10;
printf("increase time step B %g %g\n", dti, dts);
			}
			if(useGPU < 2 && dti <= 2.0 * dts){
printf("refine time steps %g %g\n", dti, dts);
				//refine interpolation points
				dts *= 0.1;
				outI *= 10;
				cOut *= 10;
			}
		}


//add condition for backward integration
		if(time >= time1){
			printf("Reached the end\n");
			break;
		}
		
	}	// end of time step loop

	cudaDeviceSynchronize();

	cudaEventRecord(tt2);
	cudaEventSynchronize(tt2);
	
	cudaEventElapsedTime(&milliseconds, tt1, tt2);
        timing[3 + S] += milliseconds * 0.001;

	cudaEventRecord(tt1);

	//###########################################
	// End time step loop
	//###########################################

	//reduce arrays for repeated integration with a smaller time step
	int k = Nperturbers;
	for(int i = Nperturbers; i < N; ++i){
		if(snew_h[i].y < 1.0){
			m_h[k] = m_h[i];
			id_h[k] = id_h[i];
			x_h[k] = x0_h[i];
			y_h[k] = y0_h[i];
			z_h[k] = z0_h[i];
			vx_h[k] = vx0_h[i];
			vy_h[k] = vy0_h[i];
			vz_h[k] = vz0_h[i];
			A1_h[k] = A1_h[i];
			A2_h[k] = A2_h[i];
			A3_h[k] = A3_h[i];

			x0_h[k] = x0_h[i];
			y0_h[k] = y0_h[i];
			z0_h[k] = z0_h[i];
			vx0_h[k] = vx0_h[i];
			vy0_h[k] = vy0_h[i];
			vz0_h[k] = vz0_h[i];
			snew_h[k].x = 1.5;
			snew_h[k].y = 1.5;
printf("%d %d %llu\n", i, k, id_h[k]);
			++k; 
		}
	}
	N = k;
	Sn[S + 1] = N;
	if(useGPU > 0){
		cudaMemcpy(m_d, m_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(id_d, id_h, N * sizeof(unsigned long long int), cudaMemcpyHostToDevice);
		cudaMemcpy(x_d, x_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(y_d, y_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(z_d, z_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(vx_d, vx_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(vy_d, vy_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(vz_d, vz_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(A1_d, A1_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(A2_d, A2_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(A3_d, A3_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(snew_d, snew_h, N * sizeof(double2), cudaMemcpyHostToDevice);
	}
	if(S == 0){
		dti = 0.1;
		dts = 0.001;
		dtiMin = 0.1;
		dt = dti * dayUnit;
		outI = (outInterval + 0.5 * dts) / dts;
		dtiOld = dti;
		time = outStart;
		cOut = 0llu;		//counter for output
		ci = 0;
	}
	if(S == 1){
		dti = 0.01;
		dts = 1.0e-3;
		dtiMin = 1.0e-4;
		dt = dti * dayUnit;
		outI = (outInterval + 0.5 * dts) / dts;
		dtiOld = dti;
		time = outStart;
		cOut = 0llu;		//counter for output
		ci = 0;
	}

}
	fclose(dtfile);
	
	printf("Time for ic and allocation, %g seconds\n", timing[0]);
	printf("Time for perturbers table, %g seconds\n", timing[1]);
	printf("Time for pre-integration, %g seconds\n", timing[2]);
	printf("Time for integration 1, %g seconds\n", timing[3]);
	printf("Time for integration 2, %g seconds\n", timing[4]);
	printf("Time for integration 3, %g seconds\n", timing[5]);

	for(int i = 0; i < 4; ++i){
		printf("N %d %d\n", i, Sn[i]-Nperturbers);
	}	


}
	

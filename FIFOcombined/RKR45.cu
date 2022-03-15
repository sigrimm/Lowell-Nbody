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
#include "integrator.h"
#include "interpolate.h"


// --------------------------------
//barycentric coordinates
__device__ __host__ void acc(double *m, double *x, double *y, double *z, double &ax, double &ay, double &az, int i, int j){

	double rx = x[j] - x[i];
	double ry = y[j] - y[i];
	double rz = z[j] - z[i];
	double rsq = rx * rx + ry * ry + rz * rz;
	double r = sqrt(rsq);

	double s = m[j] / (r * rsq);
	
	ax += s * rx;
	ay += s * ry;
	az += s * rz;
}
// ---------------------------------

// --------------------------------
//heliocentric coordinates
//sun part
__host__ void accS(double *m, double *x, double *y, double *z, double &ax, double &ay, double &az, int i){

	if(i > 0){
		double rx = -x[i];
		double ry = -y[i];
		double rz = -z[i];
		double rsq = rx * rx + ry * ry + rz * rz;
		double r = sqrt(rsq);

		double s = (m[0] + m[i]) / (r * rsq);

		ax += s * rx;
		ay += s * ry;
		az += s * rz;
	}
}
__device__ void accS_device(double mu, double xi, double yi, double zi, double &ax, double &ay, double &az){

	double rx = -xi;
	double ry = -yi;
	double rz = -zi;
	double rsq = rx * rx + ry * ry + rz * rz;
	double r = sqrt(rsq);

	double s = mu / (r * rsq);

	ax += s * rx;
	ay += s * ry;
	az += s * rz;
}
//planet part
__host__ void accP(double *m, double *x, double *y, double *z, double &ax, double &ay, double &az, int i, int j){

	if(i != j){
		double rx = x[j] - x[i];
		double ry = y[j] - y[i];
		double rz = z[j] - z[i];
		double rsq = rx * rx + ry * ry + rz * rz;
		double r = sqrt(rsq);

		double s = m[j] / (r * rsq);

		ax += s * rx;
		ay += s * ry;
		az += s * rz;
	}
}
__device__ void accP_device(double mj, double xj, double yj, double zj, double xi, double yi, double zi, double &ax, double &ay, double &az){

	double rx = xj - xi;
	double ry = yj - yi;
	double rz = zj - zi;
	double rsq = rx * rx + ry * ry + rz * rz;
	double r = sqrt(rsq);

	double s = mj / (r * rsq);

	ax += s * rx;
	ay += s * ry;
	az += s * rz;
}
//planet part 2
__host__ void accP2(double *m, double *x, double *y, double *z, double &ax, double &ay, double &az, int i, int j){

	if(i != j){
		double rx = -x[j];
		double ry = -y[j];
		double rz = -z[j];
		double rsq = rx * rx + ry * ry + rz * rz;
		double r = sqrt(rsq);

		double s = m[j] / (r * rsq);
	
		ax += s * rx;
		ay += s * ry;
		az += s * rz;
//printf("%d %d %g %g %g %g %g %g\n", i, j, rx, ry, rz, rsq, s, m[j]);
	}
}
__device__ void accP2_device(double mj, double xj, double yj, double zj, double &ax, double &ay, double &az){

	double rx = -xj;
	double ry = -yj;
	double rz = -zj;
	double rsq = rx * rx + ry * ry + rz * rz;
	double r = sqrt(rsq);

	double s = mj / (r * rsq);
	ax += s * rx;
	ay += s * ry;
	az += s * rz;
}
// --------------------------------------


//Sitarski 1982, Isotropic equation 5, heliocentric
//modified k2 to dayUnit
//should be equivalent to the Quinn et all function, assuming m[0] = 1.0
//heliocentric
__device__ __host__ void acchGR2(double xi, double yi, double zi, double vxi, double vyi, double vzi, double &ax, double &ay, double &az){
	
	double c2 = def_c * def_c;

	double rsq = xi * xi + yi * yi + zi * zi;
	double r = sqrt(rsq);
	double vsq = vxi * vxi + vyi * vyi + vzi * vzi;

	double rv = xi * vxi + yi * vyi + zi * vzi;

	double f1 = 1.0 / (r * rsq * c2);
	double t1 = 4.0 / r;
	double t2 = -vsq;
	double t3 = 4.0 * rv;
//printf("a %d %.20g %.20g %.20g\n", i, ax, ay, az);

//printf("A %d %.20g %.20g %.20g %.20g %.20g %.20g\n", i, xi, yi, zi, vxi, vyi, vzi);
//printf("B %d %.20g %.20g %.20g %.20g\n", i, f1, t1, t2, t3);
//printf("C %d %.20g %.20g %.20g %.20g\n", i, t1 + t2, (t1 + t2) * xi, ((t1 + t2) * xi + t3 * vxi), f1 * ((t1 + t2) * xi + t3 * vxi));

	ax += f1 * ((t1 + t2) * xi + t3 * vxi);
	ay += f1 * ((t1 + t2) * yi + t3 * vyi);
	az += f1 * ((t1 + t2) * zi + t3 * vzi);
//printf("D %d %.20g %.20g %.20g\n", i, ax, ay, az);

}





//A1, A2 and A3 terms for asteroids on heliocentric coordinates
__host__ __device__ void NonGrav(double xi, double yi, double zi, double vxi, double vyi, double vzi, double &ax, double &ay, double &az, double A1i, double A2i, double A3i){

	double rsq = xi * xi + yi * yi + zi * zi;
	double r = sqrt(rsq);

	//angular momenrum h = r x v
	double hx = yi * vzi - zi * vyi;
	double hy =-xi * vzi + zi * vxi;
	double hz = xi * vyi - yi * vxi;

	double hsq = hx * hx + hy * hy + hz * hz;
	double h = sqrt(hsq);

	//Transverse velocity t = h x r
	double tx = hy * zi - hz * yi;
	double ty =-hx * zi + hz * xi;
	double tz = hx * yi - hy * xi;

	double tsq = tx * tx + ty * ty + tz * tz;
	double t = sqrt(tsq);

	double gr = 1.0 / rsq;	//only valid for asteroids, not for comets 
/*
	double rr = r / R0[i];
	double g1 = pow(rr, -NM[i]);
	double g2 = pow(rr, Nn[i]);
	double g3 = pow(1.0 + g2, -NK[i]);
	double gr = ALN[i] * g1 * g3;

printf("gr %.20g %.20g\n", gr1, gr);
*/

	double f1 = A1i * gr / r;
	double f2 = A2i * gr / t;
	double f3 = A3i * gr / h;
	
	
	ax += f1 * xi + f2 * tx + f3 * hx;
	ay += f1 * yi + f2 * ty + f3 * hy;
	az += f1 * zi + f2 * tz + f3 * hz;
//printf("NonGrav %d %.20g %.20g %.20g\n", i, (f1 * x[i] + f2 * tx + f3 * hx) * dayUnit * dayUnit, (f1 * y[i] + f2 * ty + f3 * hy) * dayUnit * dayUnit, (f1 * z[i] + f2 * tz + f3 * hz) * dayUnit * dayUnit);

}

//J2 perturbation from Earth
//Walter 2018, 12.2.10
__host__ __device__ void J2(double *m, double *x, double *y, double *z, double &ax, double &ay, double &az, int i){

	double J2E = 0.00108262545; // J2 Earth from DE 430
	double RE = 6378136.3; // Earth radius in m from DE 430

	//double J2E = 1.08263e-3; //1.08262668e-3;
	//double RE = 6371.009; // Earth radius in km
	//double muE = 398600.44 // in km^3 /s^2	G * mEarth

	RE /= def_AU;	//Earth radius in AU

	int iE = 3; 	//index of Earth

	double muE = m[iE];

	//muE = 3.986004415e5 km^3s-2

	double xE = x[i] - x[iE];
	double yE = y[i] - y[iE];
	double zE = z[i] - z[iE];

	double rsq = xE * xE + yE * yE + zE * zE;
	double r = sqrt(rsq);
	double r5 = rsq * rsq * r;

	double t1 = 3.0 * J2E * muE * RE * RE / (2.0 * r5);
	double t2 = 5.0 * zE * zE / rsq;

//printf("rm %.20g %.20g %.20g\n", RE, muE, t1);

	double tx = t1 * (t2 - 1.0) * xE;
	double ty = t1 * (t2 - 1.0) * yE;
	double tz = t1 * (t2 - 3.0) * zE;
	
	ax += tx;
	ay += ty;
	az += tz;

//printf("J2 %d %.20g %.20g %.20g %.20g | %.20g %.20g %.20g\n", i, r, tx * dayUnit * dayUnit, ty * dayUnit * dayUnit, tz * dayUnit * dayUnit, xE, yE, zE); 


}

__host__ void stageStep(unsigned long long int *id, double *m, double *xt, double *yt, double *zt, double *vxt, double *vyt, double *vzt, double *kx, double *ky, double *kz, double *kvx, double *kvy, double *kvz, double *A1, double *A2, double *A3, int S, int i, int Nperturbers, int N, int useHelio, int useGR, int useJ2, int useNonGrav){

	kx[i + S * N] = vxt[i];
	ky[i + S * N] = vyt[i];
	kz[i + S * N] = vzt[i];

	double ax = 0.0;
	double ay = 0.0;
	double az = 0.0;

	if(useHelio == 0){
		for(int j = Nperturbers-1; j >= 0; --j){
			if(id[i] != id[j]){
				accP(m, xt, yt, zt, ax, ay, az, i, j);
			}
		}
	}
	else{
		for(int j = Nperturbers-1; j >= 1; --j){
			if(id[i] != id[j]){
				accP(m, xt, yt, zt, ax, ay, az, i, j);
//if(i == 27) printf("Nij %d %d %llu %llu %.20g %.20g %.20g\n", i, j, id[i], id[j], ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit);
			}
		}
		accS(m, xt, yt, zt, ax, ay, az, i);
//if(i == 27) printf("Nij %d %d %.20g %.20g %.20g\n", i, 0, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit);
//printf("N0 %d %.20g %.20g %.20g\n", i, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit);
		for(int j = Nperturbers-1; j >= 1; --j){
			if(id[i] != id[j]){
				accP2(m, xt, yt, zt, ax, ay, az, i, j);
			}
		}
//if(id[i] == 72057594037948950) printf("Np %d %.20g %.20g %.20g %d\n", i, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit, S);
	}


	if(useGR == 2){
		acchGR2(xt[i], yt[i], zt[i], vxt[i], vyt[i], vzt[i], ax, ay, az);
	}

	if(useNonGrav == 1){
		NonGrav(xt[i], yt[i], zt[i], vxt[i], vyt[i], vzt[i], ax, ay, az, A1[i], A2[i], A3[i]);
	}

	if(useJ2 == 1){
		J2(m, xt, yt, zt, ax, ay, az, i);
	}

	kvx[i + S * N] = ax;
	kvy[i + S * N] = ay;
	kvz[i + S * N] = az;
}

template < int NN >
__global__ void stageStep_kernel(unsigned long long int *id_d, double *m_d, double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *xt_d, double *yt_d, double *zt_d, double *kx_d, double *ky_d, double *kz_d, double *kvx_d, double *kvy_d, double *kvz_d, double *A1_d, double *A2_d, double *A3_d, double dt, int S, int RKFn, int Nperturbers, int N, int useHelio, int useGR, int useJ2, int useNonGrav){

	int itx = threadIdx.x;
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	//shared memory contains only the perturbers
	//the particle idx is stored in registers
	__shared__ double x_s[NN];
	__shared__ double y_s[NN];
	__shared__ double z_s[NN];
	__shared__ double m_s[NN];
	__shared__ unsigned long long int id_s[NN];


	if(threadIdx.x < Nperturbers){
			x_s[threadIdx.x] = xt_d[itx];
			y_s[threadIdx.x] = yt_d[itx];
			z_s[threadIdx.x] = zt_d[itx];
			m_s[threadIdx.x] = m_d[itx];
			id_s[threadIdx.x] = id_d[itx];
	}
	__syncthreads();

	if(idx >= Nperturbers && idx < N){


		// ***********************
		// update
		double xi = x_d[idx];
		double yi = y_d[idx];
		double zi = z_d[idx];
		double vxi = vx_d[idx];
		double vyi = vy_d[idx];
		double vzi = vz_d[idx];
		double mi = m_d[idx];
		double A1i = A1_d[idx];
		double A2i = A2_d[idx];
		double A3i = A3_d[idx];
		unsigned long long int idi = id_d[idx];

		for(int s = 0; s < S; ++s){
			double aa = a_c[S * RKFn + s];
			xi  += dt * aa * kx_d[idx + s * N];
			yi  += dt * aa * ky_d[idx + s * N];
			zi  += dt * aa * kz_d[idx + s * N];
			vxi += dt * aa * kvx_d[idx + s * N];
			vyi += dt * aa * kvy_d[idx + s * N];
			vzi += dt * aa * kvz_d[idx + s * N];
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
					accP_device(m_s[j], x_s[j], y_s[j], z_s[j], xi, yi, zi, ax, ay, az);
				}
			}
		}
		else{
			for(int j = Nperturbers-1; j >= 1; --j){
//if(idi == 72057594037948950) printf("%llu %llu\n", idi, id_s[j]);
				if(idi != id_s[j]){
					accP_device(m_s[j], x_s[j], y_s[j], z_s[j], xi, yi, zi, ax, ay, az);
///*if(idi == 72057594037948950) */printf("Nij %d %d %.20g %.20g %.20g\n", idx, j, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit);
				}
			}
			accS_device(m_s[0] + mi, xi, yi, zi, ax, ay, az);

//printf("N0 %d %.20g %.20g %.20g\n", idx, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit);
			for(int j = Nperturbers-1; j >= 1; --j){
				if(idi != id_s[j]){
					accP2_device(m_s[j], x_s[j], y_s[j], z_s[j], ax, ay, az);
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
			J2(m_d, x_d, y_d, z_d, ax, ay, az, idx);
		}

		kvx_d[idx + S * N] = ax;
		kvy_d[idx + S * N] = ay;
		kvz_d[idx + S * N] = az;
	}

}


template < int NN >
__global__ void stageStep1_kernel(unsigned long long int *id_d, double *m_d, double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *xTable_d, double *yTable_d, double *zTable_d, double *kx_d, double *ky_d, double *kz_d, double *kvx_d, double *kvy_d, double *kvz_d, double *A1_d, double *A2_d, double *A3_d, double dt, long long int sid, int RKFn, int Nperturbers, int NTable, int N, int useHelio, int useGR, int useJ2, int useNonGrav){

	int itx = threadIdx.x;
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	//shared memory contains only the perturbers
	//the particle idx is stored in registers
	__shared__ double x_s[NN];
	__shared__ double y_s[NN];
	__shared__ double z_s[NN];
	__shared__ double m_s[NN];
	__shared__ unsigned long long int id_s[NN];


	for(int S = 0; S < RKFn; ++S){


		if(threadIdx.x < Nperturbers){
				int ii = itx * NTable * RKFn + sid * RKFn + S;
				x_s[threadIdx.x] = xTable_d[ii];
				y_s[threadIdx.x] = yTable_d[ii];
				z_s[threadIdx.x] = zTable_d[ii];
				m_s[threadIdx.x] = m_d[itx];
				id_s[threadIdx.x] = id_d[itx];
	//if(itx == 1 && sid < 10) printf("%lld %d %.20g\n", sid, S, x_s[itx]);
		}
		__syncthreads();

		if(idx >= Nperturbers && idx < N){


			// ***********************
			// update
			double xi = x_d[idx];
			double yi = y_d[idx];
			double zi = z_d[idx];
			double vxi = vx_d[idx];
			double vyi = vy_d[idx];
			double vzi = vz_d[idx];
			double mi = m_d[idx];
			double A1i = A1_d[idx];
			double A2i = A2_d[idx];
			double A3i = A3_d[idx];
			unsigned long long int idi = id_d[idx];

			for(int s = 0; s < S; ++s){
				double aa = a_c[S * RKFn + s];
				xi  += dt * aa * kx_d[idx + s * N];
				yi  += dt * aa * ky_d[idx + s * N];
				zi  += dt * aa * kz_d[idx + s * N];
				vxi += dt * aa * kvx_d[idx + s * N];
				vyi += dt * aa * kvy_d[idx + s * N];
				vzi += dt * aa * kvz_d[idx + s * N];
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
						accP_device(m_s[j], x_s[j], y_s[j], z_s[j], xi, yi, zi, ax, ay, az);
					}
				}
			}
			else{
				for(int j = Nperturbers-1; j >= 1; --j){
	//if(idi == 72057594037948950) printf("%llu %llu\n", idi, id_s[j]);
					if(idi != id_s[j]){
						accP_device(m_s[j], x_s[j], y_s[j], z_s[j], xi, yi, zi, ax, ay, az);
	///*if(idi == 72057594037948950)*/ printf("Nij %d %d %.20g %.20g %.20g\n", idx, j, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit);
					}
				}
				accS_device(m_s[0] + mi, xi, yi, zi, ax, ay, az);

	//printf("N0 %d %.20g %.20g %.20g\n", idx, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit);
				for(int j = Nperturbers-1; j >= 1; --j){
					if(idi != id_s[j]){
						accP2_device(m_s[j], x_s[j], y_s[j], z_s[j], ax, ay, az);
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
				J2(m_d, x_d, y_d, z_d, ax, ay, az, idx);
			}

			kvx_d[idx + S * N] = ax;
			kvy_d[idx + S * N] = ay;
			kvz_d[idx + S * N] = az;
		}
	}
}

__host__ void update1(double *xt, double *yt, double *zt, double *vxt, double *vyt, double *vzt, double *x, double *y, double *z, double *vx, double *vy, double *vz, int i){

	xt[i] = x[i];
	yt[i] = y[i];
	zt[i] = z[i];
	vxt[i] = vx[i];
	vyt[i] = vy[i];
	vzt[i] = vz[i];
//printf("update 1 %d %g %g %g %g %g %g\n", i, xt[i], yt[i], zt[i], vxt[i], vyt[i], vzt[i]);
}


__host__ void update2(double *xt, double *yt, double *zt, double *vxt, double *vyt, double *vzt, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *kx, double *ky, double *kz, double *kvx, double *kvy, double *kvz, int i, int N, double dt, int S, int RKFn, double *a){

	xt[i]  = x[i];
	yt[i]  = y[i];
	zt[i]  = z[i];
	vxt[i] = vx[i];
	vyt[i] = vy[i];
	vzt[i] = vz[i];

	for(int s = 0; s < S; ++s){
		xt[i]  += dt * a[S * RKFn + s] * kx[i + s * N];
		yt[i]  += dt * a[S * RKFn + s] * ky[i + s * N];
		zt[i]  += dt * a[S * RKFn + s] * kz[i + s * N];
		vxt[i] += dt * a[S * RKFn + s] * kvx[i + s * N];
		vyt[i] += dt * a[S * RKFn + s] * kvy[i + s * N];
		vzt[i] += dt * a[S * RKFn + s] * kvz[i + s * N];
	}
//printf("update 2 %d %g %g %g %g %g %g\n", i, xt[i], yt[i], zt[i], vxt[i], vyt[i], vzt[i]);
}


__host__ void update(double *x, double *y, double *z, double *vx, double *vy, double *vz, double *kx, double *ky, double *kz, double *kvx, double *kvy, double *kvz, int i, int N, double dt, int RKFn, double *b){


	for(int s = 0; s < RKFn; ++s){
		x[i] += dt * b[s] * kx[i + s * N];
		y[i] += dt * b[s] * ky[i + s * N];
		z[i] += dt * b[s] * kz[i + s * N];

		vx[i] += dt * b[s] * kvx[i + s * N];
		vy[i] += dt * b[s] * kvy[i + s * N];
		vz[i] += dt * b[s] * kvz[i + s * N];
	}
}

__global__ void update_kernel(double *x, double *y, double *z, double *vx, double *vy, double *vz, double *kx, double *ky, double *kz, double *kvx, double *kvy, double *kvz, int N, int Nperturbers, double dt, int RKFn){
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if(idx >= Nperturbers && idx < N){
		for(int s = 0; s < RKFn; ++s){
			x[idx] += dt * b_c[s] * kx[idx + s * N];
			y[idx] += dt * b_c[s] * ky[idx + s * N];
			z[idx] += dt * b_c[s] * kz[idx + s * N];

			vx[idx] += dt * b_c[s] * kvx[idx + s * N];
			vy[idx] += dt * b_c[s] * kvy[idx + s * N];
			vz[idx] += dt * b_c[s] * kvz[idx + s * N];
		}
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
	int useGPU = 0;		// 0 or 1

	FILE *binfile;
	if(useFIFO == 2){	
		//binfile = fopen("210801_2342_genga_de440_perturbers.bin", "rb");
		//binfile = fopen("210921_2148_genga_in_yarkovsky_elements.bin", "rb");
		//binfile = fopen("211208_1916_genga_in_2021-12-08_specific_desig.bin", "rb");
		//binfile = fopen("210801_2104_genga_in_GA.bin", "rb");
		//binfile = fopen("210705_2315_genga_req.bin", "rb");
		binfile = fopen("220301_2048_genga_in_new_last_14_days.bin", "rb");
	}	

	int useHelio = 1;
	int outHelio = 1;
	//1 print output in heliocentric coordinates
	//0 print output in barycentric coordinates

	int OutBinary = 1;
	//0 print output files in ASCII format	
	//1 print output files in binary format	

	//long long int Nsteps = 40000;
	long long int Nsteps = 1750000;
	long long int outInterval = 1;
	double dt = 0.1 * dayUnit;

	double InVersion = 0.0;	//version of Input file format
	int DoPreIntegration = 0; //If this is 1 then do a pre-Integration to synchronize all initial conditions

	// short comparison
	//long long int Nsteps = 40000;
	//long long int Nsteps = 400;
	//long long int outInterval = 100;
	//double dt = 0.01 * dayUnit;
	
	//for timing
	//long long int Nsteps = 1000;
	//long long int outInterval = 1000;
	//double dt = 0.01 * dayUnit;


	//Integrator
	//const int RKFn = 6; //RKF45
	//const int RKFn = 7; //DP54
	const int RKFn = 13; //RKF78



	time1 = time0 + dt * Nsteps;


	if(useFIFO == 2){
		printf("read file\n");
		readHeader(binfile, time0, time1, outInterval, outStart, NTP, comet, InVersion);

		//change this later to be more general
		Nsteps = 1e9;
		dt = 0.1 * dayUnit;
		//NTP = 1;

		outInterval = 1.0;
		//outInterval = 1e8;
		//outStart = 2450800.5;
	}

	for(int i = 1; i < argc; i += 2){

		if(strcmp(argv[i], "-Nsteps") == 0){
			Nsteps = atoll(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-outInterval") == 0){
			outInterval = atoll(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-dt") == 0){
			dt = atof(argv[i + 1]);
			dt *= dayUnit;
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

	//coordinates from data table
	double *timep_h, *timep_d;
	double *xp_h, *xp_d;
	double *yp_h, *yp_d;
	double *zp_h, *zp_d;

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


	double time = 0.0;


	// ********************************
	// Read perturbers masses from perturbers.h file
	perturbersMass(m_h, Nperturbers);
	perturbersIDs(id_h, Nperturbers);
	// ********************************

	//Erase Outbinary file
	if(OutBinary == 1){
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

		/*				
		// -----------------------------------
		// Use this to extract a single object
		int ii = 29;//166;//29; //84;
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

printf("xyz %.40g %.40g %.40g %.40g %.40g %.40g %.40g %.40g %.40g %.20g\n", x_h[N], y_h[N], z_h[N], vx_h[N], vy_h[N], vz_h[N], A1_h[N], A2_h[N], A3_h[N], jd_init_h[N]);
//x_h[N] += 1.0e-6;

		// -----------------------------------
		*/

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

	//###########################################
	//perturbers table
	//###########################################
	cudaEvent_t ReadStart, ReadStop;
	cudaEventCreate(&ReadStart);
	cudaEventCreate(&ReadStop);

	cudaEventRecord(ReadStart);

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

	cudaEventRecord(ReadStop);
	cudaEventSynchronize(ReadStop);
	
	cudaEventElapsedTime(&milliseconds, ReadStart, ReadStop);
	printf("Time for perturbers table, %g seconds\n", milliseconds * 0.001);

	error = cudaGetLastError();
	printf("Perturbers error = %d = %s\n",error, cudaGetErrorString(error));
	if(error != 0.0){
		return 0;
	}
	
	//###########################################
	// end perturbers table
	//###########################################
	


	//first output
	double comx = 0.0;
	double comy = 0.0;
	double comz = 0.0;
	double vcomx = 0.0;
	double vcomy = 0.0;
	double vcomz = 0.0;


	if(useHelio == 0 && outHelio == 1){
		//convert to heliocentric output

		comx = -x_h[0];
		comy = -y_h[0];
		comz = -z_h[0];
		vcomx = -vx_h[0];
		vcomy = -vy_h[0];
		vcomz = -vz_h[0];
	}
	
	if(outHelio == 1){
		if(OutBinary == 0){
			sprintf(outfilename, "Outhelio10_%.12d.dat", 0);
		}
		else{
			sprintf(outfilename, "Outhelio10.bin");
		}
	}
	else{
		if(OutBinary == 0){
			sprintf(outfilename, "Outbary10_%.12d.dat", 0);
		}
		else{
			sprintf(outfilename, "Outbary10.bin");
		}
	}
//add condition for backward integration

	if(time0 >= outStart){
		if(OutBinary == 0){
			outfile = fopen(outfilename, "w");
		}
		else{
			outfile = fopen(outfilename, "wb");
		}

		printf("%s\n", outfilename);
		if(OutBinary == 0){
			for(int i = Nperturbers; i < N; ++i){
				fprintf(outfile, "%.10g %llu %.40g %.40g %.40g %.40g %.40g %.40g %.40g\n", time0, id_h[i], m_h[i], comx + x_h[i], comy + y_h[i], comz + z_h[i], (vcomx + vx_h[i]) * dayUnit, (vcomy + vy_h[i]) * dayUnit, (vcomz + vz_h[i]) * dayUnit);
			}
		}
		else{
			for(int i = Nperturbers; i < N; ++i){
				//unsigned long long int id = id_h[i];
				unsigned long long int id = __builtin_bswap64 (id_h[i]);
				double xx = comx + x_h[i];
				double yy = comy + y_h[i];
				double zz = comz + z_h[i];
				double vxx = (vcomx + vx_h[i]) * dayUnit;
				double vyy = (vcomy + vy_h[i]) * dayUnit;
				double vzz = (vcomz + vz_h[i]) * dayUnit;

				fwrite(&id, 1, sizeof(unsigned long long int), outfile);
				fwrite(&time0, 1, sizeof(double), outfile);
				fwrite(&xx, 1, sizeof(double), outfile);
				fwrite(&yy, 1, sizeof(double), outfile);
				fwrite(&zz, 1, sizeof(double), outfile);
				fwrite(&vxx, 1, sizeof(double), outfile);
				fwrite(&vyy, 1, sizeof(double), outfile);
				fwrite(&vzz, 1, sizeof(double), outfile);
			}

		}
		fclose(outfile);
	}

	//copy the data to the device
	if(useGPU == 1){
		cudaMemcpy(m_d, m_h, NN * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(id_d, id_h, NN * sizeof(unsigned long long int), cudaMemcpyHostToDevice);
		cudaMemcpy(x_d, x_h, NN * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(y_d, y_h, NN * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(z_d, z_h, NN * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(vx_d, vx_h, NN * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(vy_d, vy_h, NN * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(vz_d, vz_h, NN * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(A1_d, A1_h, NN * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(A2_d, A2_h, NN * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(A3_d, A3_h, NN * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(jd_init_d, jd_init_h, NN * sizeof(double), cudaMemcpyHostToDevice);
	}
	

	//remove this later to shared memory
	double *xt_h, *xt_d;
	double *yt_h, *yt_d;
	double *zt_h, *zt_d;
	double *vxt_h;
	double *vyt_h;
	double *vzt_h;
	
	double *kx_h, *kx_d;
	double *ky_h, *ky_d;
	double *kz_h, *kz_d;
	double *kvx_h, *kvx_d;
	double *kvy_h, *kvy_d;
	double *kvz_h, *kvz_d;

	double *snew_h, *snew_d;

	xt_h = (double*)malloc(NN * sizeof(double));
	yt_h = (double*)malloc(NN * sizeof(double));
	zt_h = (double*)malloc(NN * sizeof(double));
	vxt_h = (double*)malloc(NN * sizeof(double));
	vyt_h = (double*)malloc(NN * sizeof(double));
	vzt_h = (double*)malloc(NN * sizeof(double));

	kx_h = (double*)malloc(NN * RKFn * sizeof(double));
	ky_h = (double*)malloc(NN * RKFn * sizeof(double));
	kz_h = (double*)malloc(NN * RKFn * sizeof(double));
	kvx_h = (double*)malloc(NN * RKFn * sizeof(double));
	kvy_h = (double*)malloc(NN * RKFn * sizeof(double));
	kvz_h = (double*)malloc(NN * RKFn * sizeof(double));

	snew_h = (double*)malloc(NN * sizeof(double));


	//interpolation table
	const int NTable1 = 10000;
	double *xTable_d;
	double *yTable_d;
	double *zTable_d;


	if(useGPU == 1){
		cudaMalloc((void **) &xt_d, NN * sizeof(double));
		cudaMalloc((void **) &yt_d, NN * sizeof(double));
		cudaMalloc((void **) &zt_d, NN * sizeof(double));

		cudaMalloc((void **) &kx_d, NN * RKFn * sizeof(double));
		cudaMalloc((void **) &ky_d, NN * RKFn * sizeof(double));
		cudaMalloc((void **) &kz_d, NN * RKFn * sizeof(double));
		cudaMalloc((void **) &kvx_d, NN * RKFn * sizeof(double));
		cudaMalloc((void **) &kvy_d, NN * RKFn * sizeof(double));
		cudaMalloc((void **) &kvz_d, NN * RKFn * sizeof(double));

		cudaMalloc((void **) &snew_d, NN * sizeof(double));

		cudaMalloc((void **) &xTable_d, Nperturbers * NTable1 * RKFn * sizeof(double));
		cudaMalloc((void **) &yTable_d, Nperturbers * NTable1 * RKFn * sizeof(double));
		cudaMalloc((void **) &zTable_d, Nperturbers * NTable1 * RKFn * sizeof(double));
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
		ee = 1.0 / 4.0;
	}
	else if(RKFn == 7){
		setDP54(a_h, b_h, bb_h, c_h);
		ee = 1.0 / 5.0;
	}
	else if(RKFn == 13){
		setRKF78(a_h, b_h, bb_h, c_h);
		ee = 1.0 / 7.0;
	}
	else{
		printf("RKFn values not valid %d\n", RKFn);
		return 0;
	}

	if(useGPU == 1){
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
	
	//###########################################
	// Start time step loop
	//###########################################

	cudaEvent_t LoopStart, LoopStop;
	cudaEventCreate(&LoopStart);
	cudaEventCreate(&LoopStop);

	cudaEventRecord(LoopStart);

	double dtOld = dt;

	if(DoPreIntegration == 1){
		//Do pre-Integration to syncrhonize all initial conditions
		printf("Start pre-Integration up to time %.20g\n", outStart);
		for(int i = Nperturbers; i < N; ++i){
			time = jd_init_h[i];
			double snew = 1.0;
			dt = 1.0 * dayUnit;
			int stop = 0;

printf("preIntegration %d %.20g %.20g\n", i, time, outStart);

			if(time == outStart){
				continue;
			}
			if(time > outStart){
				dt = -dt;
			}
			
			for(long long int t = 1; t <= Nsteps; ++t){
//printf("A %d %lld %.20g %.20g %.20g %.20g\n", i, id_h[i], time, x_h[i], y_h[i], z_h[i]);

				for(int S = 0; S < RKFn; ++S){
					for(int p = 0; p < Nperturbers; ++p){
						interpolate(Ninterpolate, Nperturbers, xp_h, yp_h, zp_h, timep_h, timep0, dtimep, time + c_h[S] * dt / dayUnit, xt_h, yt_h, zt_h, p);
						//interpolate2(Ninterpolate, Nperturbers, xp_h, yp_h, zp_h, timep_h, timep0, dtimep, time + c_h[S] * dt / dayUnit, xt_h, yt_h, zt_h, p);
					}
					if(S == 0) update1(xt_h, yt_h, zt_h, vxt_h, vyt_h, vzt_h, x_h, y_h, z_h, vx_h, vy_h, vz_h, i);
					else update2(xt_h, yt_h, zt_h, vxt_h, vyt_h, vzt_h, x_h, y_h, z_h, vx_h, vy_h, vz_h, kx_h, ky_h, kz_h, kvx_h, kvy_h, kvz_h, i, N, dt, S, RKFn, a_h);

					stageStep(id_h, m_h, xt_h, yt_h, zt_h, vxt_h, vyt_h, vzt_h, kx_h, ky_h, kz_h, kvx_h, kvy_h, kvz_h, A1_h, A2_h, A3_h, S, i, Nperturbers, N, useHelio, useGR, useJ2, useNonGrav);

				}
				computeError1(snew_h, kx_h, ky_h, kz_h, kvx_h, kvy_h, kvz_h, b_h, bb_h, RKFn, i, N, snew, dt, ee);

				if(snew >= 1.0){		
					update(x_h, y_h, z_h, vx_h, vy_h, vz_h, kx_h, ky_h, kz_h, kvx_h, kvy_h, kvz_h, i, N, dt, RKFn, b_h);	
		
					time += dt / dayUnit;

					dt *= snew;

					if(stop == 1){
						break;
					}
				}
				else{
					dt *= snew;
					stop = 0;
				}


//printf("B %d %lld %.20g %.20g %.20g %.20g %.20g %.20g\n", i, id_h[i], time, dt / dayUnit, snew, x_h[i], y_h[i], z_h[i]);


				if(dt > 0.0 && time + dt / dayUnit > outStart){
					dt = (outStart - time) * dayUnit;
					stop = 1;
//printf("Final time step %.20g\n", dt/dayUnit);
				}
				if(dt < 0.0 && time + dt / dayUnit < outStart){
					dt = (outStart - time) * dayUnit;
					stop = 1;
//printf("Final time step %.20g\n", dt/dayUnit);
				}

			}
printf("C %d %lld %.20g %.20g %.20g %.20g %.20g %.20g\n", i, id_h[i], time, dt / dayUnit, snew, x_h[i], y_h[i], z_h[i]);
			if(fabs(time - outStart) > 1.0e-10){
				printf("Error in pre-integration of particle %d, start time not reached\n", i);
				return 0;
			}
		}
		dt = dtOld;
	}


return 0;



	time = time0;
	int printOutput = 0;
	double nextOuttime = outInterval + time0;

	int cTable = NTable1;	//counter for interpolation table

	for(long long int t = 1; t <= Nsteps; ++t){
	//for(long long int t = 1; t < 3; ++t){

cudaDeviceSynchronize();
printf("%lld %d | %d %d\n", t, 0, NTP, Nperturbers);	

		if(useGPU == 1 && cTable >= NTable1){
			//Precalculate interpolation points on 0.1 time steps for NTable1 steps and store them in global memory

			//interpolateTable_kernel < Ninterpolate > <<< dim3(Nperturbers, NTable1, RKFn), dim3(Ninterpolate) >>> (Nperturbers, NTable1, RKFn, xp_d, yp_d, zp_d, timep_d, timep0, dtimep, time, xTable_d, yTable_d, zTable_d);
			interpolate2bTable_kernel < Ninterpolate > <<< dim3((NTable1 + 255) / 256, Nperturbers, RKFn), dim3(256) >>> (Nperturbers, NTable1, RKFn, xp_d, yp_d, zp_d, timep_d, timep0, dtimep, time, xTable_d, yTable_d, zTable_d);

			cTable = 0;

			cudaDeviceSynchronize();
			error = cudaGetLastError();
			printf("Interpolate table  error = %d = %s\n",error, cudaGetErrorString(error));
			if(error != 0.0){
				return 0;
			}

			//return 0;
		}
		++cTable;

		double newTime = time + dt / dayUnit;
		double snew = 1.0;

		if(useGPU == 0){
			for(int S = 0; S < RKFn; ++S){
				for(int p = 0; p < Nperturbers; ++p){
					interpolate(Ninterpolate, Nperturbers, xp_h, yp_h, zp_h, timep_h, timep0, dtimep, time + c_h[S] * dt / dayUnit, xt_h, yt_h, zt_h, p);
					//interpolate2(Ninterpolate, Nperturbers, xp_h, yp_h, zp_h, timep_h, timep0, dtimep, time + c_h[S] * dt / dayUnit, xt_h, yt_h, zt_h, p);
				}
				for(int i = Nperturbers; i < N; ++i){
					if(S == 0) update1(xt_h, yt_h, zt_h, vxt_h, vyt_h, vzt_h, x_h, y_h, z_h, vx_h, vy_h, vz_h, i);
					else update2(xt_h, yt_h, zt_h, vxt_h, vyt_h, vzt_h, x_h, y_h, z_h, vx_h, vy_h, vz_h, kx_h, ky_h, kz_h, kvx_h, kvy_h, kvz_h, i, N, dt, S, RKFn, a_h);
				}
				for(int i = Nperturbers; i < N; ++i){
					stageStep(id_h, m_h, xt_h, yt_h, zt_h, vxt_h, vyt_h, vzt_h, kx_h, ky_h, kz_h, kvx_h, kvy_h, kvz_h, A1_h, A2_h, A3_h, S, i, Nperturbers, N, useHelio, useGR, useJ2, useNonGrav);
				}


			}
			//computeError(snew_h, kx_h, ky_h, kz_h, kvx_h, kvy_h, kvz_h, b_h, bb_h, RKFn, Nperturbers, N, snew, dt, ee);
		// -----------------------------------
		}
		else{
			//for(int S = 0; S < RKFn; ++S){
//cudaDeviceSynchronize();
//printf("%lld %d\n", t, S);	
				//interpolate_kernel < Ninterpolate > <<< dim3(Nperturbers, NTP,1), Ninterpolate >>> (Nperturbers, xp_d, yp_d, zp_d, timep_d, timep0, dtimep, time + c_h[S] * dt / dayUnit, xt_d, yt_d, zt_d);
				//interpolate2b_kernel < Ninterpolate > <<< NTP, Nperturbers >>> (Nperturbers, xp_d, yp_d, zp_d, timep_d, timep0, dtimep, time + c_h[S] * dt / dayUnit, xt_d, yt_d, zt_d);
				//stageStep_kernel < Nperturbers > <<< (NN + 127) / 128, 128 >>> (id_d, m_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, xt_d, yt_d, zt_d, kx_d, ky_d, kz_d, kvx_d, kvy_d, kvz_d, A1_d, A2_d, A3_d, dt, S, RKFn, Nperturbers, N, useHelio, useGR, useJ2, useNonGrav);

			//}
			stageStep1_kernel < Nperturbers > <<< (NN + 127) / 128, 128 >>> (id_d, m_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, xTable_d, yTable_d, zTable_d, kx_d, ky_d, kz_d, kvx_d, kvy_d, kvz_d, A1_d, A2_d, A3_d, dt, cTable - 1, RKFn, Nperturbers, NTable1, N, useHelio, useGR, useJ2, useNonGrav);
cudaDeviceSynchronize();
		}


//printf("%.20g %.20g dt: %.20g %.20g\n", time, newTime, dt / dayUnit, snew);
printf("%.20g %.20g dt: %.20g %.20g | %.20g %.20g %.20g\n", time, newTime, dt / dayUnit, snew, x_h[Nperturbers], y_h[Nperturbers], z_h[Nperturbers]);
//printf("%.20g %.20g\n", nextOuttime, outStart);


/*
		if(printOutput == 2){
			printOutput = 1;
		}

		if(newTime > nextOuttime){
			dtOld = dt;
			dt = (nextOuttime - time) * dayUnit;
//printf("   correct %.20g %.20g %.20g %.20g\n", newTime, nextOuttime, dt / dayUnit, dtOld / dayUnit);
			snew = 1.0;
			nextOuttime += outInterval;
			printOutput = 2;
		}
*/
		if(snew >= 1.0 && printOutput != 2){		
			if(useGPU == 0){
				for(int i = Nperturbers; i < N; ++i){
					update(x_h, y_h, z_h, vx_h, vy_h, vz_h, kx_h, ky_h, kz_h, kvx_h, kvy_h, kvz_h, i, N, dt, RKFn, b_h);	
				}
			}
			else{
				update_kernel <<< (NN + 127) / 128, 128 >>> (x_d, y_d, z_d, vx_d, vy_d, vz_d, kx_d, ky_d, kz_d, kvx_d, kvy_d, kvz_d, N, Nperturbers, dt, RKFn);	
			}
		
			time += dt / dayUnit;

			dt *= snew;
		}
		else{
		
			dt *= snew;
		}


//add condition for backward integration
		//if(printOutput == 1 && time >= outStart){
		if(t % 10 == 0){

			printf("Output %.20g %lld\n", time, t);

			//printOutput = 0;
			//dt = dtOld;


			if(useGPU == 1){
				cudaMemcpy(x_h, x_d, NN * sizeof(double), cudaMemcpyDeviceToHost);
				cudaMemcpy(y_h, y_d, NN * sizeof(double), cudaMemcpyDeviceToHost);
				cudaMemcpy(z_h, z_d, NN * sizeof(double), cudaMemcpyDeviceToHost);
				cudaMemcpy(vx_h, vx_d, NN * sizeof(double), cudaMemcpyDeviceToHost);
				cudaMemcpy(vy_h, vy_d, NN * sizeof(double), cudaMemcpyDeviceToHost);
				cudaMemcpy(vz_h, vz_d, NN * sizeof(double), cudaMemcpyDeviceToHost);
			}
			
			if(outHelio == 1){
				if(OutBinary == 0){	
					sprintf(outfilename, "Outhelio10_%.12lld.dat", t);
				}
				else{
					sprintf(outfilename, "Outhelio10.bin");
				}
			}
			else{
				if(OutBinary == 0){
					sprintf(outfilename, "Outbary10_%.12lld.dat", t);
				}
				else{
					sprintf(outfilename, "Outbary10.bin");
				}
			}
			if(OutBinary == 0){
				outfile = fopen(outfilename, "w");
			}
			else{
				outfile = fopen(outfilename, "ab");
			}
//			printf("%s\n", outfilename);

			comx = 0.0;
			comy = 0.0;
			comz = 0.0;
			vcomx = 0.0;
			vcomy = 0.0;
			vcomz = 0.0;

			if(useHelio == 0 && outHelio == 1){
				//convert to heliocentric output
				comx = -x_h[0];
				comy = -y_h[0];
				comz = -z_h[0];
				vcomx = -vx_h[0];
				vcomy = -vy_h[0];
				vcomz = -vz_h[0];
			}
			
			if(OutBinary == 0){
				for(int i = Nperturbers; i < N; ++i){
					fprintf(outfile, "%.10g %llu %.40g %.40g %.40g %.40g %.40g %.40g %.40g\n", time, id_h[i], m_h[i], comx + x_h[i], comy + y_h[i], comz + z_h[i], (vcomx + vx_h[i]) * dayUnit, (vcomy + vy_h[i]) * dayUnit, (vcomz + vz_h[i]) * dayUnit);
				}
			}
			else{
				for(int i = Nperturbers; i < N; ++i){

					//unsigned long long int id = id_h[i];
					unsigned long long int id = __builtin_bswap64 (id_h[i]);
					double xx = comx + x_h[i];
					double yy = comy + y_h[i];
					double zz = comz + z_h[i];
					double vxx = (vcomx + vx_h[i]) * dayUnit;
					double vyy = (vcomy + vy_h[i]) * dayUnit;
					double vzz = (vcomz + vz_h[i]) * dayUnit;

					fwrite(&id, 1, sizeof(unsigned long long int), outfile);
					fwrite(&time, 1, sizeof(double), outfile);
					fwrite(&xx, 1, sizeof(double), outfile);
					fwrite(&yy, 1, sizeof(double), outfile);
					fwrite(&zz, 1, sizeof(double), outfile);
					fwrite(&vxx, 1, sizeof(double), outfile);
					fwrite(&vyy, 1, sizeof(double), outfile);
					fwrite(&vzz, 1, sizeof(double), outfile);

					//printf("%llu %g %g %g %g %g %g %g\n", id, time, xx, yy, zz, vxx, vyy, vzz);
				}
				
			}
			fclose(outfile);


		}
		if(time > time1){
			printf("Reached the end\n");
			break;
		}
		
	}	// end of time step loop

	//###########################################
	// End time step loop
	//###########################################

	cudaDeviceSynchronize();

	cudaEventRecord(LoopStop);
	cudaEventSynchronize(LoopStop);
	
	cudaEventElapsedTime(&milliseconds, LoopStart, LoopStop);
	printf("Time for integration, %g seconds\n", milliseconds * 0.001);
	
}
	

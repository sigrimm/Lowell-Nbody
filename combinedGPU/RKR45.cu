#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define dayUnit 0.01720209895
//#define dayUnit 0.01720209894846
//#define dayUnit 1.0

//#define def_c 10065.3201686
#define def_c 10065.320121
#define def_AU 149597870700.0           //AU in m

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
__device__ __host__ void accS(double *m, double *x, double *y, double *z, double &ax, double &ay, double &az, int i){

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
//planet part
__device__ __host__ void accP(double *m, double *x, double *y, double *z, double &ax, double &ay, double &az, int i, int j){

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
//planet part 2
__device__ __host__ void accP2(double *m, double *x, double *y, double *z, double &ax, double &ay, double &az, int i, int j){

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
	}
}
// --------------------------------------


//Sitarski 1982, Isotropic equation 5, heliocentric
//modified k2 to dayUnit
//should be equivalent to the Quinn et all function, assuming m[0] = 1.0
//heliocentric
__device__ __host__ void acchGR2(double *m, double *x, double *y, double *z, double *vx, double *vy, double *vz, double &ax, double &ay, double &az, int i){
	
	double c2 = def_c * def_c;

	double rsq = x[i] * x[i] + y[i] * y[i] + z[i] * z[i];
	double r = sqrt(rsq);
	double vsq = vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];

	double rv = x[i] * vx[i] + y[i] * vy[i] + z[i] * vz[i];

	double f1 = 1.0 / (r * rsq * c2);
	double t1 = 4.0 / r;
	double t2 = -vsq;
	double t3 = 4.0 * rv;
//printf("a %d %.20g %.20g %.20g\n", i, ax, ay, az);

//printf("A %d %.20g %.20g %.20g %.20g %.20g %.20g\n", i, x[i], y[i], z[i], vx[i], vy[i], vz[i]);
//printf("B %d %.20g %.20g %.20g %.20g\n", i, f1, t1, t2, t3);
//printf("C %d %.20g %.20g %.20g %.20g\n", i, t1 + t2, (t1 + t2) * x[i], ((t1 + t2) * x[i] + t3 * vx[i]), f1 * ((t1 + t2) * x[i] + t3 * vx[i]));


	ax += f1 * ((t1 + t2) * x[i] + t3 * vx[i]);
	ay += f1 * ((t1 + t2) * y[i] + t3 * vy[i]);
	az += f1 * ((t1 + t2) * z[i] + t3 * vz[i]);
//printf("D %d %.20g %.20g %.20g\n", i, ax, ay, az);

}
//barycentric
__device__ __host__ void accbGR2(double *m, double *x, double *y, double *z, double *vx, double *vy, double *vz, double &ax, double &ay, double &az, int i){
	
	double c2 = def_c * def_c;

	double rhx = x[i] - x[0];
	double rhy = y[i] - y[0];
	double rhz = z[i] - z[0];

	double vhx = vx[i] - vx[0];
	double vhy = vy[i] - vy[0];
	double vhz = vz[i] - vz[0];

	double rsq = rhx * rhx + rhy * rhy + rhz * rhz;
	double r = sqrt(rsq);
	double vsq = vhx * vhx + vhy * vhy + vhz * vhz;

	double rv = rhx * vhx + rhy * vhy + rhz * vhz;

	double f1 = 1.0 / (r * rsq * c2);
	double t1 = 4.0 / r;
	double t2 = -vsq;
	double t3 = 4.0 * rv;
//printf("a %d %.20g %.20g %.20g\n", i, ax, ay, az);

//printf("A %d %.20g %.20g %.20g %.20g %.20g %.20g\n", i, rhx, rhy, rhz, vhx, vhy, vhz);
//printf("B %d %.20g %.20g %.20g %.20g\n", i, f1, t1, t2, t3);
//printf("C %d %.20g %.20g %.20g %.20g\n", i, t1 + t2, (t1 + t2) * rhx, ((t1 + t2) * rhx + t3 * vhx), f1 * ((t1 + t2) * rhx + t3 * vhx));

	ax += f1 * ((t1 + t2) * rhx + t3 * vhx);
	ay += f1 * ((t1 + t2) * rhy + t3 * vhy);
	az += f1 * ((t1 + t2) * rhz + t3 * vhz);
//printf("D%d %.20g %.20g %.20g\n", i, ax, ay, az);

}


//Neville-Aitken interpolation
template <int NN>
__host__ __device__ void interpolate(int N, double *xp, double *yp, double *zp, double *timep, double time, double *xt, double *yt, double *zt, int p){


	double Px[NN][NN];
	double Py[NN][NN];
	double Pz[NN][NN];
	double tn[NN];

	for(int i = 0; i < N; ++i){
		Px[0][i] = xp[p * N + i];
		Py[0][i] = yp[p * N + i];
		Pz[0][i] = zp[p * N + i];
		tn[i] = timep[p * N + i];


//printf("interpolate %d %d %.20g %.20g %.20g\n", p, i, time, tn[i], P[0][i]);
	}

	for(int j = 1; j < N; ++j){
//printf("****\n");
		for(int i = 0; i < N - j; ++i){
			Px[j][i] = ((time - tn[i+j]) * Px[j-1][i] + (tn[i] - time) * Px[j-1][i+1]) / (tn[i] - tn[i+j]);
			Py[j][i] = ((time - tn[i+j]) * Py[j-1][i] + (tn[i] - time) * Py[j-1][i+1]) / (tn[i] - tn[i+j]);
			Pz[j][i] = ((time - tn[i+j]) * Pz[j-1][i] + (tn[i] - time) * Pz[j-1][i+1]) / (tn[i] - tn[i+j]);
//printf("%d %d %g %g %g %g %.20g\n", i, i+j, tn[i], tn[i+j], P[j-1][i], P[j-1][i+1], P[j][i]);

		}
	}
	xt[p] = Px[N-1][0];
	yt[p] = Py[N-1][0];
	zt[p] = Pz[N-1][0];
//printf("%.20g %d %.20g %.20g %.20g\n", time, p, xt[p], yt[p], zt[p]);

}

template <int NN>
__global__ void interpolate_kernel(int N, int Ninterpolate, int Nperturbers, double *xp, double *yp, double *zp, double *timep, double time, double *xt, double *yt, double *zt){

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if(idx < Nperturbers){
		interpolate < NN> (Ninterpolate, xp, yp, zp, timep, time, xt, yt, zt, idx);
	}

}


//A1, A2 and A3 terms for asteroids on heliocentric coordinates
void NonGrav(double *x, double *y, double *z, double *vx, double *vy, double *vz, double &ax, double &ay, double &az, double *A1, double *A2, double *A3, double *ALN, double *NK, double *NM, double *Nn, double *R0, int i){

	double rsq = x[i] * x[i] + y[i] * y[i] + z[i] * z[i];
	double r = sqrt(rsq);

	//angular momenrum h = r x v
	double hx = y[i] * vz[i] - z[i] * vy[i];
	double hy =-x[i] * vz[i] + z[i] * vx[i];
	double hz = x[i] * vy[i] - y[i] * vx[i];

	double hsq = hx * hx + hy * hy + hz * hz;
	double h = sqrt(hsq);

	//Transverse velocity t = h x r
	double tx = hy * z[i] - hz * y[i];
	double ty =-hx * z[i] + hz * x[i];
	double tz = hx * y[i] - hy * x[i];

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

	double f1 = A1[i] * gr / r;
	double f2 = A2[i] * gr / t;
	double f3 = A3[i] * gr / h;
	
	
	ax += f1 * x[i] + f2 * tx + f3 * hx;
	ay += f1 * y[i] + f2 * ty + f3 * hy;
	az += f1 * z[i] + f2 * tz + f3 * hz;
//printf("NonGrav %d %.20g %.20g %.20g\n", i, (f1 * x[i] + f2 * tx + f3 * hx) * dayUnit * dayUnit, (f1 * y[i] + f2 * ty + f3 * hy) * dayUnit * dayUnit, (f1 * z[i] + f2 * tz + f3 * hz) * dayUnit * dayUnit);

}

//J2 perturbation from Earth
void J2(double *m, double *x, double *y, double *z, double &ax, double &ay, double &az, int i){

	double J2E = 0.00108262545; // J2 Earth from DE 430
	double RE = 6378136.3; // Earth radius in m from DE 430

	//double J2E =  1.08263e-3; //  1.08262668e-3;
	//double RE = 6371.009; // Earth radius in km
	//double muE = 398600.44 // in km^3 /s^2	G * mEarth

	RE /= def_AU;	//Earth radius in AU

	int iE = 3; 	//index of Earth

	double muE = m[iE];

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
 
//printf("J2 %d %.20g %.20g %.20g | %.20g %.20g %.20g\n", i, tx * dayUnit * dayUnit, ty * dayUnit * dayUnit, tz * dayUnit * dayUnit, xE, yE, zE); 


}

__host__ __device__ void stageStep(double *m, double *xt, double *yt, double *zt, double *vxt, double *vyt, double *vzt, double *kx, double *ky, double *kz, double *kvx, double *kvy, double *kvz, int S, int i, int Nperturbers, int N, int useHelio, int GR){

	kx[i + S * N] = vxt[i];
	ky[i + S * N] = vyt[i];
	kz[i + S * N] = vzt[i];

	double ax = 0.0;
	double ay = 0.0;
	double az = 0.0;

	if(useHelio == 0){
		for(int j = Nperturbers-1; j >= 0; --j){
			accP(m, xt, yt, zt, ax, ay, az, i, j);
		}
	}
	else{
		for(int j = Nperturbers-1; j >= 1; --j){
			accP(m, xt, yt, zt, ax, ay, az, i, j);
//if(i == 27) printf("Nij %d %d %.20g %.20g %.20g\n", i, j, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit);
		}
		accS(m, xt, yt, zt, ax, ay, az, i);
//if(i == 27) printf("Nij %d %d %.20g %.20g %.20g\n", i, 0, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit);
//printf("N0 %d %.20g %.20g %.20g\n", i, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit);
		for(int j = Nperturbers-1; j >= 1; --j){
			accP2(m, xt, yt, zt, ax, ay, az, i, j);
		}
//printf("Np %d %.20g %.20g %.20g %d\n", i, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit, S);
	}


	if(GR == 2){
		acchGR2(m, xt, yt, zt, vxt, vyt, vzt, ax, ay, az, i);
	}

	//NonGrav(xt, yt, zt, vxt, vyt, vzt, ax, ay, az, A1, A2, A3, ALN, NK, NM, Nn, R0, i);
	//J2(m, xt, yt, zt, ax, ay, az, i);

	kvx[i + S * N] = ax;
	kvy[i + S * N] = ay;
	kvz[i + S * N] = az;
}

__global__ void stageStep_kernel(double *m, double *xt, double *yt, double *zt, double *vxt, double *vyt, double *vzt, double *kx, double *ky, double *kz, double *kvx, double *kvy, double *kvz, int S, int Nperturbers, int N, int useHelio, int GR){

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if(idx >= Nperturbers && idx < N){
                        stageStep(m, xt, yt, zt, vxt, vyt, vzt, kx, ky, kz, kvx, kvy, kvz, S, idx, Nperturbers, N, useHelio, GR);
	}

}

__host__ void update1(double *xt, double *yt, double *zt, double *vxt, double *vyt, double *vzt, double *x, double *y, double *z, double *vx, double *vy, double *vz, int i){

	xt[i] = x[i];
	yt[i] = y[i];
	zt[i] = z[i];
	vxt[i] = vx[i];
	vyt[i] = vy[i];
	vzt[i] = vz[i];
}

__global__ void update1_kernel(double *xt, double *yt, double *zt, double *vxt, double *vyt, double *vzt, double *x, double *y, double *z, double *vx, double *vy, double *vz, int N){

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if(idx < N){

		xt[idx] = x[idx];
		yt[idx] = y[idx];
		zt[idx] = z[idx];
		vxt[idx] = vx[idx];
		vyt[idx] = vy[idx];
		vzt[idx] = vz[idx];
	}
}

__host__ void update2(double *xt, double *yt, double *zt, double *vxt, double *vyt, double *vzt, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *kx, double *ky, double *kz, double *kvx, double *kvy, double *kvz, int i, int N, double dt, double a21){

	xt[i]  = x[i]  + dt * a21 * kx[i + 0 * N];
	yt[i]  = y[i]  + dt * a21 * ky[i + 0 * N];
	zt[i]  = z[i]  + dt * a21 * kz[i + 0 * N];
	vxt[i] = vx[i] + dt * a21 * kvx[i + 0 * N];
	vyt[i] = vy[i] + dt * a21 * kvy[i + 0 * N];
	vzt[i] = vz[i] + dt * a21 * kvz[i + 0 * N];
}

__global__ void update2_kernel(double *xt, double *yt, double *zt, double *vxt, double *vyt, double *vzt, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *kx, double *ky, double *kz, double *kvx, double *kvy, double *kvz, int N, int Nperturbers, double dt, double a21){
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if(idx >= Nperturbers && idx < N){

		xt[idx]  = x[idx]  + dt * a21 * kx[idx + 0 * N];
		yt[idx]  = y[idx]  + dt * a21 * ky[idx + 0 * N];
		zt[idx]  = z[idx]  + dt * a21 * kz[idx + 0 * N];
		vxt[idx] = vx[idx] + dt * a21 * kvx[idx + 0 * N];
		vyt[idx] = vy[idx] + dt * a21 * kvy[idx + 0 * N];
		vzt[idx] = vz[idx] + dt * a21 * kvz[idx + 0 * N];
	}
}

__host__ void update3(double *xt, double *yt, double *zt, double *vxt, double *vyt, double *vzt, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *kx, double *ky, double *kz, double *kvx, double *kvy, double *kvz, int i, int N, double dt, double a31, double a32){

	xt[i]  = x[i]  + dt * (a31 * kx[i + 0 * N]  + a32 * kx[i + 1 * N]);
	yt[i]  = y[i]  + dt * (a31 * ky[i + 0 * N]  + a32 * ky[i + 1 * N]);
	zt[i]  = z[i]  + dt * (a31 * kz[i + 0 * N]  + a32 * kz[i + 1 * N]);
	vxt[i] = vx[i] + dt * (a31 * kvx[i + 0 * N] + a32 * kvx[i + 1 * N]);
	vyt[i] = vy[i] + dt * (a31 * kvy[i + 0 * N] + a32 * kvy[i + 1 * N]);
	vzt[i] = vz[i] + dt * (a31 * kvz[i + 0 * N] + a32 * kvz[i + 1 * N]);
}
__global__ void update3_kernel(double *xt, double *yt, double *zt, double *vxt, double *vyt, double *vzt, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *kx, double *ky, double *kz, double *kvx, double *kvy, double *kvz, int N, int Nperturbers, double dt, double a31, double a32){
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if(idx >= Nperturbers && idx < N){
		xt[idx]  = x[idx]  + dt * (a31 * kx[idx + 0 * N]  + a32 * kx[idx + 1 * N]);
		yt[idx]  = y[idx]  + dt * (a31 * ky[idx + 0 * N]  + a32 * ky[idx + 1 * N]);
		zt[idx]  = z[idx]  + dt * (a31 * kz[idx + 0 * N]  + a32 * kz[idx + 1 * N]);
		vxt[idx] = vx[idx] + dt * (a31 * kvx[idx + 0 * N] + a32 * kvx[idx + 1 * N]);
		vyt[idx] = vy[idx] + dt * (a31 * kvy[idx + 0 * N] + a32 * kvy[idx + 1 * N]);
		vzt[idx] = vz[idx] + dt * (a31 * kvz[idx + 0 * N] + a32 * kvz[idx + 1 * N]);
	}
}

__host__ void update4(double *xt, double *yt, double *zt, double *vxt, double *vyt, double *vzt, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *kx, double *ky, double *kz, double *kvx, double *kvy, double *kvz, int i, int N, double dt, double a41, double a42, double a43){

	xt[i]  = x[i]  + dt * (a41 * kx[i + 0 * N]  + a42 * kx[i + 1 * N]  + a43 * kx[i + 2 * N]);
	yt[i]  = y[i]  + dt * (a41 * ky[i + 0 * N]  + a42 * ky[i + 1 * N]  + a43 * ky[i + 2 * N]);
	zt[i]  = z[i]  + dt * (a41 * kz[i + 0 * N]  + a42 * kz[i + 1 * N]  + a43 * kz[i + 2 * N]);
	vxt[i] = vx[i] + dt * (a41 * kvx[i + 0 * N] + a42 * kvx[i + 1 * N] + a43 * kvx[i + 2 * N]);
	vyt[i] = vy[i] + dt * (a41 * kvy[i + 0 * N] + a42 * kvy[i + 1 * N] + a43 * kvy[i + 2 * N]);
	vzt[i] = vz[i] + dt * (a41 * kvz[i + 0 * N] + a42 * kvz[i + 1 * N] + a43 * kvz[i + 2 * N]);
}
__global__ void update4_kernel(double *xt, double *yt, double *zt, double *vxt, double *vyt, double *vzt, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *kx, double *ky, double *kz, double *kvx, double *kvy, double *kvz, int N, int Nperturbers, double dt, double a41, double a42, double a43){
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if(idx >= Nperturbers && idx < N){

		xt[idx]  = x[idx]  + dt * (a41 * kx[idx + 0 * N]  + a42 * kx[idx + 1 * N]  + a43 * kx[idx + 2 * N]);
		yt[idx]  = y[idx]  + dt * (a41 * ky[idx + 0 * N]  + a42 * ky[idx + 1 * N]  + a43 * ky[idx + 2 * N]);
		zt[idx]  = z[idx]  + dt * (a41 * kz[idx + 0 * N]  + a42 * kz[idx + 1 * N]  + a43 * kz[idx + 2 * N]);
		vxt[idx] = vx[idx] + dt * (a41 * kvx[idx + 0 * N] + a42 * kvx[idx + 1 * N] + a43 * kvx[idx + 2 * N]);
		vyt[idx] = vy[idx] + dt * (a41 * kvy[idx + 0 * N] + a42 * kvy[idx + 1 * N] + a43 * kvy[idx + 2 * N]);
		vzt[idx] = vz[idx] + dt * (a41 * kvz[idx + 0 * N] + a42 * kvz[idx + 1 * N] + a43 * kvz[idx + 2 * N]);
	}
}

__host__ void update5(double *xt, double *yt, double *zt, double *vxt, double *vyt, double *vzt, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *kx, double *ky, double *kz, double *kvx, double *kvy, double *kvz, int i, int N, double dt, double a51, double a52, double a53, double a54){

	xt[i]  = x[i]  + dt * (a51 * kx[i + 0 * N]  + a52 * kx[i + 1 * N]  + a53 * kx[i + 2 * N]  + a54 * kx[i + 3 * N]);
	yt[i]  = y[i]  + dt * (a51 * ky[i + 0 * N]  + a52 * ky[i + 1 * N]  + a53 * ky[i + 2 * N]  + a54 * ky[i + 3 * N]);
	zt[i]  = z[i]  + dt * (a51 * kz[i + 0 * N]  + a52 * kz[i + 1 * N]  + a53 * kz[i + 2 * N]  + a54 * kz[i + 3 * N]);
	vxt[i] = vx[i] + dt * (a51 * kvx[i + 0 * N] + a52 * kvx[i + 1 * N] + a53 * kvx[i + 2 * N] + a54 * kvx[i + 3 * N]);
	vyt[i] = vy[i] + dt * (a51 * kvy[i + 0 * N] + a52 * kvy[i + 1 * N] + a53 * kvy[i + 2 * N] + a54 * kvy[i + 3 * N]);
	vzt[i] = vz[i] + dt * (a51 * kvz[i + 0 * N] + a52 * kvz[i + 1 * N] + a53 * kvz[i + 2 * N] + a54 * kvz[i + 3 * N]);
}

__global__ void update5_kernel(double *xt, double *yt, double *zt, double *vxt, double *vyt, double *vzt, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *kx, double *ky, double *kz, double *kvx, double *kvy, double *kvz, int N, int Nperturbers, double dt, double a51, double a52, double a53, double a54){
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if(idx >= Nperturbers && idx < N){

		xt[idx]  = x[idx]  + dt * (a51 * kx[idx + 0 * N]  + a52 * kx[idx + 1 * N]  + a53 * kx[idx + 2 * N]  + a54 * kx[idx + 3 * N]);
		yt[idx]  = y[idx]  + dt * (a51 * ky[idx + 0 * N]  + a52 * ky[idx + 1 * N]  + a53 * ky[idx + 2 * N]  + a54 * ky[idx + 3 * N]);
		zt[idx]  = z[idx]  + dt * (a51 * kz[idx + 0 * N]  + a52 * kz[idx + 1 * N]  + a53 * kz[idx + 2 * N]  + a54 * kz[idx + 3 * N]);
		vxt[idx] = vx[idx] + dt * (a51 * kvx[idx + 0 * N] + a52 * kvx[idx + 1 * N] + a53 * kvx[idx + 2 * N] + a54 * kvx[idx + 3 * N]);
		vyt[idx] = vy[idx] + dt * (a51 * kvy[idx + 0 * N] + a52 * kvy[idx + 1 * N] + a53 * kvy[idx + 2 * N] + a54 * kvy[idx + 3 * N]);
		vzt[idx] = vz[idx] + dt * (a51 * kvz[idx + 0 * N] + a52 * kvz[idx + 1 * N] + a53 * kvz[idx + 2 * N] + a54 * kvz[idx + 3 * N]);
	}
}
__host__ void update6(double *xt, double *yt, double *zt, double *vxt, double *vyt, double *vzt, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *kx, double *ky, double *kz, double *kvx, double *kvy, double *kvz, int i, int N, double dt, double a61, double a62, double a63, double a64, double a65){
	
	xt[i]  = x[i]  + dt * (a61 * kx[i + 0 * N]  + a62 * kx[i + 1 * N]  + a63 * kx[i + 2 * N]  + a64 * kx[i + 3 * N]  + a65 * kx[i + 4 * N]);
	yt[i]  = y[i]  + dt * (a61 * ky[i + 0 * N]  + a62 * ky[i + 1 * N]  + a63 * ky[i + 2 * N]  + a64 * ky[i + 3 * N]  + a65 * ky[i + 4 * N]);
	zt[i]  = z[i]  + dt * (a61 * kz[i + 0 * N]  + a62 * kz[i + 1 * N]  + a63 * kz[i + 2 * N]  + a64 * kz[i + 3 * N]  + a65 * kz[i + 4 * N]);
	vxt[i] = vx[i] + dt * (a61 * kvx[i + 0 * N] + a62 * kvx[i + 1 * N] + a63 * kvx[i + 2 * N] + a64 * kvx[i + 3 * N] + a65 * kvx[i + 4 * N]);
	vyt[i] = vy[i] + dt * (a61 * kvy[i + 0 * N] + a62 * kvy[i + 1 * N] + a63 * kvy[i + 2 * N] + a64 * kvy[i + 3 * N] + a65 * kvy[i + 4 * N]);
	vzt[i] = vz[i] + dt * (a61 * kvz[i + 0 * N] + a62 * kvz[i + 1 * N] + a63 * kvz[i + 2 * N] + a64 * kvz[i + 3 * N] + a65 * kvz[i + 4 * N]);
}
__global__ void update6_kernel(double *xt, double *yt, double *zt, double *vxt, double *vyt, double *vzt, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *kx, double *ky, double *kz, double *kvx, double *kvy, double *kvz, int N, int Nperturbers, double dt, double a61, double a62, double a63, double a64, double a65){
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if(idx >= Nperturbers && idx < N){
		
		xt[idx]  = x[idx]  + dt * (a61 * kx[idx + 0 * N]  + a62 * kx[idx + 1 * N]  + a63 * kx[idx + 2 * N]  + a64 * kx[idx + 3 * N]  + a65 * kx[idx + 4 * N]);
		yt[idx]  = y[idx]  + dt * (a61 * ky[idx + 0 * N]  + a62 * ky[idx + 1 * N]  + a63 * ky[idx + 2 * N]  + a64 * ky[idx + 3 * N]  + a65 * ky[idx + 4 * N]);
		zt[idx]  = z[idx]  + dt * (a61 * kz[idx + 0 * N]  + a62 * kz[idx + 1 * N]  + a63 * kz[idx + 2 * N]  + a64 * kz[idx + 3 * N]  + a65 * kz[idx + 4 * N]);
		vxt[idx] = vx[idx] + dt * (a61 * kvx[idx + 0 * N] + a62 * kvx[idx + 1 * N] + a63 * kvx[idx + 2 * N] + a64 * kvx[idx + 3 * N] + a65 * kvx[idx + 4 * N]);
		vyt[idx] = vy[idx] + dt * (a61 * kvy[idx + 0 * N] + a62 * kvy[idx + 1 * N] + a63 * kvy[idx + 2 * N] + a64 * kvy[idx + 3 * N] + a65 * kvy[idx + 4 * N]);
		vzt[idx] = vz[idx] + dt * (a61 * kvz[idx + 0 * N] + a62 * kvz[idx + 1 * N] + a63 * kvz[idx + 2 * N] + a64 * kvz[idx + 3 * N] + a65 * kvz[idx + 4 * N]);
	}
}

__host__ void update(double *x, double *y, double *z, double *vx, double *vy, double *vz, double *kx, double *ky, double *kz, double *kvx, double *kvy, double *kvz, int i, int N, double dt, double *b){
	//RKF45
	x[i] += dt * (b[0] * kx[i + 0 * N] + b[2] * kx[i + 2 * N] + b[3] * kx[i + 3 * N] + b[4] * kx[i + 4 * N] + b[5] * kx[i + 5 * N]);
	y[i] += dt * (b[0] * ky[i + 0 * N] + b[2] * ky[i + 2 * N] + b[3] * ky[i + 3 * N] + b[4] * ky[i + 4 * N] + b[5] * ky[i + 5 * N]);
	z[i] += dt * (b[0] * kz[i + 0 * N] + b[2] * kz[i + 2 * N] + b[3] * kz[i + 3 * N] + b[4] * kz[i + 4 * N] + b[5] * kz[i + 5 * N]);

	vx[i] += dt * (b[0] * kvx[i + 0 * N] + b[2] * kvx[i + 2 * N] + b[3] * kvx[i + 3 * N] + b[4] * kvx[i + 4 * N] + b[5] * kvx[i + 5 * N]);
	vy[i] += dt * (b[0] * kvy[i + 0 * N] + b[2] * kvy[i + 2 * N] + b[3] * kvy[i + 3 * N] + b[4] * kvy[i + 4 * N] + b[5] * kvy[i + 5 * N]);
	vz[i] += dt * (b[0] * kvz[i + 0 * N] + b[2] * kvz[i + 2 * N] + b[3] * kvz[i + 3 * N] + b[4] * kvz[i + 4 * N] + b[5] * kvz[i + 5 * N]);
}

__global__ void update_kernel(double *x, double *y, double *z, double *vx, double *vy, double *vz, double *kx, double *ky, double *kz, double *kvx, double *kvy, double *kvz, int N, int Nperturbers, double dt, double *b){
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if(idx >= Nperturbers && idx < N){
	//RKF45
		x[idx] += dt * (b[0] * kx[idx + 0 * N] + b[2] * kx[idx + 2 * N] + b[3] * kx[idx + 3 * N] + b[4] * kx[idx + 4 * N] + b[5] * kx[idx + 5 * N]);
		y[idx] += dt * (b[0] * ky[idx + 0 * N] + b[2] * ky[idx + 2 * N] + b[3] * ky[idx + 3 * N] + b[4] * ky[idx + 4 * N] + b[5] * ky[idx + 5 * N]);
		z[idx] += dt * (b[0] * kz[idx + 0 * N] + b[2] * kz[idx + 2 * N] + b[3] * kz[idx + 3 * N] + b[4] * kz[idx + 4 * N] + b[5] * kz[idx + 5 * N]);

		vx[idx] += dt * (b[0] * kvx[idx + 0 * N] + b[2] * kvx[idx + 2 * N] + b[3] * kvx[idx + 3 * N] + b[4] * kvx[idx + 4 * N] + b[5] * kvx[idx + 5 * N]);
		vy[idx] += dt * (b[0] * kvy[idx + 0 * N] + b[2] * kvy[idx + 2 * N] + b[3] * kvy[idx + 3 * N] + b[4] * kvy[idx + 4 * N] + b[5] * kvy[idx + 5 * N]);
		vz[idx] += dt * (b[0] * kvz[idx + 0 * N] + b[2] * kvz[idx + 2 * N] + b[3] * kvz[idx + 3 * N] + b[4] * kvz[idx + 4 * N] + b[5] * kvz[idx + 5 * N]);
	}
}

int main(int argc, char*argv[]){

	//Number of planets
	const int NN = 28;	//22
	const int Nperturbers = 27;	//21
	const int Ninterpolate = 10;	//number of interpolation points
	const double dtime = 1.0;     //interval between stored time steps

	int GR = 2;
	//2 Sitarski 1982, heliocentric coordinates

	int useHelio = 1;
	int outHelio = 1;
	//1 print output in heliocentric coordinates
	//0 print output in barycentric  coordinates

	//long long int Nsteps = 40000;	
	//long long int outInterval = 10;
	//double dt = 0.1 * dayUnit;

	//long long int Nsteps = 400000;	
	long long int Nsteps = 40000;
	long long int outInterval = 100;
	double dt = 0.01 * dayUnit;

	for(int i = 1; i < argc; i += 2){

		if(strcmp(argv[i], "-Nsteps") == 0){
			Nsteps = atoll(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-outInterval") == 0){
			outInterval = atoll(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-dt") == 0){
			dt = atof(argv[i + 1]);
		}
	}

	
	int *id_h, *id_d;
	double *m_h, *m_d;
	double *x_h, *x_d;
	double *y_h, *y_d;
	double *z_h, *z_d;
	double *vx_h, *vx_d;
	double *vy_h, *vy_d;
	double *vz_h, *vz_d;

	//coordinates from data table
	double *timep_h, *timep_d;
	double *xp_h, *xp_d;
	double *yp_h, *yp_d;
	double *zp_h, *zp_d;


	//allocate data on host
	id_h = (int*)malloc(NN * sizeof(int));
	m_h = (double*)malloc(NN * sizeof(double));
	x_h = (double*)malloc(NN * sizeof(double));
	y_h = (double*)malloc(NN * sizeof(double));
	z_h = (double*)malloc(NN * sizeof(double));
	vx_h = (double*)malloc(NN * sizeof(double));
	vy_h = (double*)malloc(NN * sizeof(double));
	vz_h = (double*)malloc(NN * sizeof(double));

	timep_h = (double*)malloc(Nperturbers * Ninterpolate * sizeof(double));
	xp_h = (double*)malloc(Nperturbers * Ninterpolate * sizeof(double));
	yp_h = (double*)malloc(Nperturbers * Ninterpolate * sizeof(double));
	zp_h = (double*)malloc(Nperturbers * Ninterpolate * sizeof(double));

	//allocate data on the device
	cudaMalloc((void **) &id_d, NN * sizeof(int));
	cudaMalloc((void **) &m_d, NN * sizeof(double));
	cudaMalloc((void **) &x_d, NN * sizeof(double));
	cudaMalloc((void **) &y_d, NN * sizeof(double));
	cudaMalloc((void **) &z_d, NN * sizeof(double));
	cudaMalloc((void **) &vx_d, NN * sizeof(double));
	cudaMalloc((void **) &vy_d, NN * sizeof(double));
	cudaMalloc((void **) &vz_d, NN * sizeof(double));

	cudaMalloc((void **) &timep_d, Nperturbers * Ninterpolate * sizeof(double));
	cudaMalloc((void **) &xp_d, Nperturbers * Ninterpolate * sizeof(double));
	cudaMalloc((void **) &yp_d, Nperturbers * Ninterpolate * sizeof(double));
	cudaMalloc((void **) &zp_d, Nperturbers * Ninterpolate * sizeof(double));

	/*
	//non Gravitational constants
	double A1[NN];
	double A2[NN];
	double A3[NN];
	double ALN[NN];
	double NK[NN];
	double NM[NN];
	double Nn[NN];
	double R0[NN];

	for(int i = 0; i < NN; ++i){
		A1[i] = 0.0;
		A2[i] = 0.0;
		A3[i] = 0.0;
		ALN[i] = 0.0;
		NK[i] = 0.0;
		NM[i] = 0.0;
		Nn[i] = 0.0;
		R0[i] = 1.0;

	}
	*/

	FILE *outfile;
	char outfilename[160];	


	double time = 0.0;

	//Units are 1/(mass of object in solar masses)
	double pmass[] = {
		1.000000000000000e0,      // Sun        (0)
		6023682.155592479e0,      // Mercury    (1)
		408523.7186582996e0,      // Venus      (2)
		332946.0488339480e0,      // Earth      (3)
		3098703.590290707e0,      // Mars       (4)
		1047.348625463337e0,      // Jupiter    (5)
		3497.901767786633e0,      // Saturn     (6)
		22902.98161308703e0,      // Uranus     (7)
		19412.25977597307e0,      // Neptune    (8)
		135836683.7686175e0,      // Pluto      (9)
		2112939391.8192508,                  // Ceres      (10)
		9531877787.0654011,                  // Pallas     (11)
		81799329362.428986,                  // Juno       (12)
		7676559929.1351004,                  // Vesta      (13)
		23944976514.662392,                  // Hygiea     (14)
		63251980219.354561,                  // Eunomia    (15)
		46649712166.264168,                  // Euphrosyne (16)
		119474172269.94408,                  // Europa     (17)
		56926698684.931702,                  // Davida     (18)
		56298080671.641434,                  // Interamnia (19)
		27068703.24120323e0,      // Moon       (20)

		86737410876.841156,		  // Psyche
		93034865412.812271,		  // Cybele
		114823090351.20033,		  // Thisbe
		116910898662.48077,		  // Doris
		128906361339.41116,		  // Patientia
		134548655333.38321,		  // Sylvia

		0.0			  // test particle

	};
	for(int i = 0; i < Nperturbers; ++i){
		m_h[i] = 1.0/pmass[i];
printf("m %d %.20g\n", i, m_h[i]);
	}

//m[Nperturbers] = 1.e-11; //ca mass of Flora

	int N = Nperturbers;

	//Sun
	FILE *infile;
	char infilename[160];

	sprintf(infilename, "initial.dat");
	infile = fopen(infilename, "r");
	//sun
	id_h[0] = 20;
	x_h[0] = 0.0;
	y_h[0] = 0.0;
	z_h[0] = 0.0;
	vx_h[0] = 0.0;
	vy_h[0] = 0.0;
	vz_h[0] = 0.0;

	for(int i = 1; i < NN; ++i){
		id_h[i] = -1;
		x_h[i] = 0.0;
		y_h[i] = 0.0;
		z_h[i] = 0.0;
		vx_h[i] = 0.0;
		vy_h[i] = 0.0;
		vz_h[i] = 0.0;
	}
	//read test particle
	for(int i = Nperturbers; i < NN; ++i){
		int er = 0;
		fscanf(infile, "%lf", &time);
		fscanf(infile, "%lf", &x_h[i]);
		fscanf(infile, "%lf", &y_h[i]);
		fscanf(infile, "%lf", &z_h[i]);
		fscanf(infile, "%lf", &vx_h[i]);
		fscanf(infile, "%lf", &vy_h[i]);
		er= fscanf(infile, "%lf", &vz_h[i]);
		//fscanf(infile, "%lf", &A1[i]);
		//fscanf(infile, "%lf", &A2[i]);
		//fscanf(infile, "%lf", &A3[i]);
		//fscanf(infile, "%lf", &ALN[i]);
		//fscanf(infile, "%lf", &NK[i]);
		//fscanf(infile, "%lf", &NM[i]);
		//fscanf(infile, "%lf", &Nn[i]);
		//er = fscanf(infile, "%lf", &R0[i]);
		if(er < 0) break;
		++N;
printf("er %d %d %d %d %.20g %.20g %.20g\n", i, id_h[i], er, N, x_h[i], y_h[i], z_h[i]);
	}
	fclose(infile);
	double time0 = time;	//start time from simulation
	double time1 = time;	//time from table position

	for(int i = Nperturbers; i < N; ++i){
		vx_h[i] /= dayUnit;
		vy_h[i] /= dayUnit;
		vz_h[i] /= dayUnit;

		//A1[i] /= (dayUnit * dayUnit);
		//A2[i] /= (dayUnit * dayUnit);
		//A3[i] /= (dayUnit * dayUnit);
	}



	FILE *XVfile;
	if(useHelio == 1){
		XVfile = fopen("All_h.dat", "r");
	}
	else{
		XVfile = fopen("All_b.dat", "r");
	}


	// -----------------------------------------
	//Read table
	int countNodes = 0;
	for(int t = 0; t < 1000000; ++t){
		int er;
printf("CountNodes %d\n", countNodes);
		for(int i = 0; i < Nperturbers; ++i){
			double skip;
			double timepp;
			int id;
			er = fscanf(XVfile, "%lf %d", &timepp, &id);
			fscanf(XVfile, "%lf %lf %lf", &xp_h[id * Ninterpolate + countNodes], &yp_h[id * Ninterpolate + countNodes], &zp_h[id * Ninterpolate + countNodes]);
			//remove velocities read later
			fscanf(XVfile, "%lf %lf %lf", &skip, &skip, &skip);

			if(er < 0) break;
			timep_h[id * Ninterpolate + countNodes] = timepp;

printf("read %.20g %d %.20g %.20g %.20g%d\n", timep_h[id * Ninterpolate + countNodes], id, xp_h[id * Ninterpolate + countNodes], yp_h[id * Ninterpolate + countNodes], zp_h[id * Ninterpolate + countNodes], id * Ninterpolate + countNodes);

			//vxp[id * Ninterpolate + countNodes] /= dayUnit;
			//vyp[id * Ninterpolate + countNodes] /= dayUnit;
			//vzp[id * Ninterpolate + countNodes] /= dayUnit;


			if(i == 0 && t == 0 && timep_h[id * Ninterpolate + countNodes] > time - (Ninterpolate/2 - 1) * dtime){
				printf("Error, time too small, not enough data before time\n");
				return 0;
			}
			if(i == Nperturbers - 1 && timep_h[id * Ninterpolate + countNodes] > time - Ninterpolate/2 * dtime){
				++countNodes;
			}
		}
		if(er < 0) break;
		if(countNodes >= Ninterpolate){
			break;
		}
	}
	if(countNodes < Ninterpolate){
		printf("Error, time too large, not enough data after time\n");
		return 0;
	}
	// ---------------------------------------

	

	// ---------------------------------------
	//interpolate on host
	for(int p = 0; p < Nperturbers; ++p){
		interpolate <NN> (Ninterpolate, xp_h, yp_h, zp_h, timep_h, time, x_h, y_h, z_h, p);
	}
	// ---------------------------------------

	

	//first output
	double comx = 0.0;
	double comy = 0.0;
	double comz = 0.0;
	double vcomx = 0.0;
	double vcomy = 0.0;
	double vcomz = 0.0;
	double mtot = 0.0;	


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
		sprintf(outfilename, "Outhelio10_%.12d.dat", 0);
	}
	else{
		sprintf(outfilename, "Outbary10_%.12d.dat", 0);
	}
	outfile = fopen(outfilename, "w");
	printf("%s\n", outfilename);
	for(int i = 0; i < N; ++i){
		fprintf(outfile, "%.10g %d %.40g %.40g %.40g %.40g %.40g %.40g %.40g\n", time, i, m_h[i], comx + x_h[i], comy + y_h[i], comz + z_h[i], vcomx + vx_h[i], vcomy + vy_h[i], vcomz + vz_h[i]);

	}
	fclose(outfile);


	//copy the data to the device
	cudaMemcpy(m_d, m_h, NN * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(x_d, x_h, NN * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(y_d, y_h, NN * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(z_d, z_h, NN * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(vx_d, vx_h, NN * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(vy_d, vy_h, NN * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(vz_d, vz_h, NN * sizeof(double), cudaMemcpyHostToDevice);

	cudaMemcpy(timep_d, timep_h, Nperturbers * Ninterpolate * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(xp_d, xp_h, Nperturbers * Ninterpolate * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(yp_d, yp_h, Nperturbers * Ninterpolate * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(zp_d, zp_h, Nperturbers * Ninterpolate * sizeof(double), cudaMemcpyHostToDevice);

	

	//remove this later to shared memory
	double *xt_h, *xt_d;
	double *yt_h, *yt_d;
	double *zt_h, *zt_d;
	double *vxt_h, *vxt_d;
	double *vyt_h, *vyt_d;
	double *vzt_h, *vzt_d;
	
	double *kx_h, *kx_d;
	double *ky_h, *ky_d;
	double *kz_h, *kz_d;
	double *kvx_h, *kvx_d;
	double *kvy_h, *kvy_d;
	double *kvz_h, *kvz_d;

	xt_h = (double*)malloc(NN * sizeof(double));
	yt_h = (double*)malloc(NN * sizeof(double));
	zt_h = (double*)malloc(NN * sizeof(double));
	vxt_h = (double*)malloc(NN * sizeof(double));
	vyt_h = (double*)malloc(NN * sizeof(double));
	vzt_h = (double*)malloc(NN * sizeof(double));

	kx_h = (double*)malloc(NN * 6 * sizeof(double));
	ky_h = (double*)malloc(NN * 6 * sizeof(double));
	kz_h = (double*)malloc(NN * 6 * sizeof(double));
	kvx_h = (double*)malloc(NN * 6 * sizeof(double));
	kvy_h = (double*)malloc(NN * 6 * sizeof(double));
	kvz_h = (double*)malloc(NN * 6 * sizeof(double));

	cudaMalloc((void **) &xt_d, NN * sizeof(double));
	cudaMalloc((void **) &yt_d, NN * sizeof(double));
	cudaMalloc((void **) &zt_d, NN * sizeof(double));
	cudaMalloc((void **) &vxt_d, NN * sizeof(double));
	cudaMalloc((void **) &vyt_d, NN * sizeof(double));
	cudaMalloc((void **) &vzt_d, NN * sizeof(double));

	cudaMalloc((void **) &kx_d, NN * 6 * sizeof(double));
	cudaMalloc((void **) &ky_d, NN * 6 * sizeof(double));
	cudaMalloc((void **) &kz_d, NN * 6 * sizeof(double));
	cudaMalloc((void **) &kvx_d, NN * 6 * sizeof(double));
	cudaMalloc((void **) &kvy_d, NN * 6 * sizeof(double));
	cudaMalloc((void **) &kvz_d, NN * 6 * sizeof(double));

	//double errorkx[N];
	//double errorky[N];
	//double errorkz[N];
	//double errorkvx[N];
	//double errorkvy[N];
	//double errorkvz[N];

	

	double a21 = 1.0/4.0;

	double a31 = 3.0/32.0;
	double a32 = 9.0/32.0;

	double a41 = 1932.0/2197.0;
	double a42 = -7200.0/2197.0;
	double a43 = 7296.0/2197.0;

	double a51 = 439.0/216.0;
	double a52 = -8.0;
	double a53 = 3680.0/513.0;
	double a54 = -845.0/4104.0;

	double a61 = -8.0/27.0;
	double a62 = 2.0;
	double a63 = -3544/2565.0;
	double a64 = 1859.0/4104.0;
	double a65 = -11.0/40.0;

	double b[6];
	b[0] = 25.0/216.0;
	b[1] = 0.0;
	b[2] = 1408.0/2565.0;
	b[3] = 2197.0/4104.0;
	b[4] = -1.0/5.0;
	b[5] = 0.0;

	double bb[6];
	bb[0] = 16.0/135.0;
	bb[1] = 0.0;
	bb[2] = 6656.0/12825.0;
	bb[3] = 28561.0/56430.0;
	bb[4] = -9.0/50.0;
	bb[5] = 2.0/55.0;


	double cc[6];
	cc[0] = 0.0;
	cc[1] = 0.25;
	cc[2] = 3.0 / 8.0;
	cc[3] = 12.0 / 13.0;
	cc[4] = 1.0;
	cc[5] = 0.5;

	int S;
	for(long long int t = 1; t <= Nsteps; ++t){

			
		//stage 1
		S = 0;
		for(int i = 0; i < N; ++i){
			update1(xt_h, yt_h, zt_h, vxt_h, vyt_h, vzt_h, x_h, y_h, z_h, vx_h, vy_h, vz_h, i);
		}
		update1_kernel <<< 1, NN >>> (xt_d, yt_d, zt_d, vxt_d, vyt_d, vzt_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, N);

		
		for(int i = Nperturbers; i < N; ++i){
			stageStep(m_h, xt_h, yt_h, zt_h, vxt_h, vyt_h, vzt_h, kx_h, ky_h, kz_h, kvx_h, kvy_h, kvz_h, S, i, Nperturbers, N, useHelio, GR);
		}
		stageStep_kernel <<< 1, NN >>> (m_d, xt_d, yt_d, zt_d, vxt_d, vyt_d, vzt_d, kx_d, ky_d, kz_d, kvx_d, kvy_d, kvz_d, S, Nperturbers, N, useHelio, GR);
	
				
		//stage 2
		S = 1;
		for(int p = 0; p < Nperturbers; ++p){
			interpolate < NN > (Ninterpolate, xp_h, yp_h, zp_h, timep_h, time + cc[S] * dt / dayUnit, xt_h, yt_h, zt_h, p);
		}
		interpolate_kernel < NN > <<< 1, NN >>> (N, Ninterpolate, Nperturbers, xp_d, yp_d, zp_d, timep_d, time + cc[S] * dt / dayUnit, xt_d, yt_d, zt_d);

		for(int i = Nperturbers; i < N; ++i){
			update2(xt_h, yt_h, zt_h, vxt_h, vyt_h, vzt_h, x_h, y_h, z_h, vx_h, vy_h, vz_h, kx_h, ky_h, kz_h, kvx_h, kvy_h, kvz_h, i, N, dt, a21);
		}
		update2_kernel <<< 1, NN >>> (xt_d, yt_d, zt_d, vxt_d, vyt_d, vzt_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, kx_d, ky_d, kz_d, kvx_d, kvy_d, kvz_h, N, Nperturbers, dt, a21);

		for(int i = Nperturbers; i < N; ++i){
			stageStep(m_h, xt_h, yt_h, zt_h, vxt_h, vyt_h, vzt_h, kx_h, ky_h, kz_h, kvx_h, kvy_h, kvz_h, S, i, Nperturbers, N, useHelio, GR);
		}
		stageStep_kernel <<< 1, NN >>> (m_d, xt_d, yt_d, zt_d, vxt_d, vyt_d, vzt_d, kx_d, ky_d, kz_d, kvx_d, kvy_d, kvz_d, S, Nperturbers, N, useHelio, GR);

		
		//stage 3
		S = 2;
		for(int p = 0; p < Nperturbers; ++p){
			interpolate < NN > (Ninterpolate, xp_h, yp_h, zp_h, timep_h, time + cc[S] * dt / dayUnit, xt_h, yt_h, zt_h, p);
		}
		interpolate_kernel < NN > <<< 1, NN >>> (N, Ninterpolate, Nperturbers, xp_d, yp_d, zp_d, timep_d, time + cc[S] * dt / dayUnit, xt_d, yt_d, zt_d);

		for(int i = Nperturbers; i < N; ++i){
			update3(xt_h, yt_h, zt_h, vxt_h, vyt_h, vzt_h, x_h, y_h, z_h, vx_h, vy_h, vz_h, kx_h, ky_h, kz_h, kvx_h, kvy_h, kvz_h, i, N, dt, a31, a32);
		}
		update3_kernel <<< 1, NN >>> (xt_d, yt_d, zt_d, vxt_d, vyt_d, vzt_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, kx_d, ky_d, kz_d, kvx_d, kvy_d, kvz_h, N, Nperturbers, dt, a31, a32);

		for(int i = Nperturbers; i < N; ++i){
			stageStep(m_h, xt_h, yt_h, zt_h, vxt_h, vyt_h, vzt_h, kx_h, ky_h, kz_h, kvx_h, kvy_h, kvz_h, S, i, Nperturbers, N, useHelio, GR);
		}
		stageStep_kernel <<< 1, NN >>> (m_d, xt_d, yt_d, zt_d, vxt_d, vyt_d, vzt_d, kx_d, ky_d, kz_d, kvx_d, kvy_d, kvz_d, S, Nperturbers, N, useHelio, GR);

		//stage 4
		S = 3;
		for(int p = 0; p < Nperturbers; ++p){
			interpolate < NN > (Ninterpolate, xp_h, yp_h, zp_h, timep_h, time + cc[S] * dt / dayUnit, xt_h, yt_h, zt_h, p);
		}
		interpolate_kernel < NN > <<< 1, NN >>> (N, Ninterpolate, Nperturbers, xp_d, yp_d, zp_d, timep_d, time + cc[S] * dt / dayUnit, xt_d, yt_d, zt_d);

		for(int i = Nperturbers; i < N; ++i){
			update4(xt_h, yt_h, zt_h, vxt_h, vyt_h, vzt_h, x_h, y_h, z_h, vx_h, vy_h, vz_h, kx_h, ky_h, kz_h, kvx_h, kvy_h, kvz_h, i, N, dt, a41, a42, a43);
		}
		update4_kernel <<< 1, NN >>> (xt_d, yt_d, zt_d, vxt_d, vyt_d, vzt_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, kx_d, ky_d, kz_d, kvx_d, kvy_d, kvz_h, N, Nperturbers, dt, a41, a42, a43);

		for(int i = Nperturbers; i < N; ++i){
			stageStep(m_h, xt_h, yt_h, zt_h, vxt_h, vyt_h, vzt_h, kx_h, ky_h, kz_h, kvx_h, kvy_h, kvz_h, S, i, Nperturbers, N, useHelio, GR);
		}
		stageStep_kernel <<< 1, NN >>> (m_d, xt_d, yt_d, zt_d, vxt_d, vyt_d, vzt_d, kx_d, ky_d, kz_d, kvx_d, kvy_d, kvz_d, S, Nperturbers, N, useHelio, GR);

		//stage 5
		S = 4;
		for(int p = 0; p < Nperturbers; ++p){
			interpolate < NN > (Ninterpolate, xp_h, yp_h, zp_h, timep_h, time + cc[S] * dt / dayUnit, xt_h, yt_h, zt_h, p);
		}
		interpolate_kernel < NN > <<< 1, NN >>> (N, Ninterpolate, Nperturbers, xp_d, yp_d, zp_d, timep_d, time + cc[S] * dt / dayUnit, xt_d, yt_d, zt_d);

		for(int i = Nperturbers; i < N; ++i){
			update5(xt_h, yt_h, zt_h, vxt_h, vyt_h, vzt_h, x_h, y_h, z_h, vx_h, vy_h, vz_h, kx_h, ky_h, kz_h, kvx_h, kvy_h, kvz_h, i, N, dt, a51, a52, a53, a54);
		}
		update5_kernel <<< 1, NN >>> (xt_d, yt_d, zt_d, vxt_d, vyt_d, vzt_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, kx_d, ky_d, kz_d, kvx_d, kvy_d, kvz_h, N, Nperturbers, dt, a51, a52, a53, a54);

		for(int i = Nperturbers; i < N; ++i){
			stageStep(m_h, xt_h, yt_h, zt_h, vxt_h, vyt_h, vzt_h, kx_h, ky_h, kz_h, kvx_h, kvy_h, kvz_h, S, i, Nperturbers, N, useHelio, GR);
		}
		stageStep_kernel <<< 1, NN >>> (m_d, xt_d, yt_d, zt_d, vxt_d, vyt_d, vzt_d, kx_d, ky_d, kz_d, kvx_d, kvy_d, kvz_d, S, Nperturbers, N, useHelio, GR);

		//stage 6
		S = 5;
		for(int p = 0; p < Nperturbers; ++p){
			interpolate < NN > (Ninterpolate, xp_h, yp_h, zp_h, timep_h, time + cc[S] * dt / dayUnit, xt_h, yt_h, zt_h, p);
		}
		interpolate_kernel < NN > <<< 1, NN >>> (N, Ninterpolate, Nperturbers, xp_d, yp_d, zp_d, timep_d, time + cc[S] * dt / dayUnit, xt_d, yt_d, zt_d);
		for(int i = Nperturbers; i < N; ++i){
			update6(xt_h, yt_h, zt_h, vxt_h, vyt_h, vzt_h, x_h, y_h, z_h, vx_h, vy_h, vz_h, kx_h, ky_h, kz_h, kvx_h, kvy_h, kvz_h, i, N, dt, a61, a62, a63, a64, a65);
		}
		update6_kernel <<< 1, NN >>> (xt_d, yt_d, zt_d, vxt_d, vyt_d, vzt_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, kx_d, ky_d, kz_d, kvx_d, kvy_d, kvz_h, N, Nperturbers, dt, a61, a62, a63, a64, a65);

		for(int i = Nperturbers; i < N; ++i){
			stageStep(m_h, xt_h, yt_h, zt_h, vxt_h, vyt_h, vzt_h, kx_h, ky_h, kz_h, kvx_h, kvy_h, kvz_h, S, i, Nperturbers, N, useHelio, GR);
		}
		stageStep_kernel <<< 1, NN >>> (m_d, xt_d, yt_d, zt_d, vxt_d, vyt_d, vzt_d, kx_d, ky_d, kz_d, kvx_d, kvy_d, kvz_d, S, Nperturbers, N, useHelio, GR);
		
		double sc = 1.0e-15;

		/*	
		//error estimation
		for(int i = 0; i < N; ++i){
			errorkx[i] = ((b1 - bb1) * kx[i][0] + (b3 - bb3) * kx[i][2] + (b4 - bb4) * kx[i][3] + (b5 - bb5) * kx[i][4] + (b6 - bb6) * kx[i][5]) / sc;
			errorky[i] = ((b1 - bb1) * ky[i][0] + (b3 - bb3) * ky[i][2] + (b4 - bb4) * ky[i][3] + (b5 - bb5) * ky[i][4] + (b6 - bb6) * ky[i][5]) / sc;
			errorkz[i] = ((b1 - bb1) * kz[i][0] + (b3 - bb3) * kz[i][2] + (b4 - bb4) * kz[i][3] + (b5 - bb5) * kz[i][4] + (b6 - bb6) * kz[i][5]) / sc;
			errorkvx[i] = ((b1 - bb1) * kvx[i][0] + (b3 - bb3) * kvx[i][2] + (b4 - bb4) * kvx[i][3] + (b5 - bb5) * kvx[i][4] + (b6 - bb6) * kvx[i][5]) / sc;
			errorkvy[i] = ((b1 - bb1) * kvy[i][0] + (b3 - bb3) * kvy[i][2] + (b4 - bb4) * kvy[i][3] + (b5 - bb5) * kvy[i][4] + (b6 - bb6) * kvy[i][5]) / sc;
			errorkvz[i] = ((b1 - bb1) * kvz[i][0] + (b3 - bb3) * kvz[i][2] + (b4 - bb4) * kvz[i][3] + (b5 - bb5) * kvz[i][4] + (b6 - bb6) * kvz[i][5]) / sc;
		}


		double error = 0.0;
		for(int i = 0; i < N; ++i){
			error += errorkx[i] * errorkx[i];
			error += errorky[i] * errorky[i];
			error += errorkz[i] * errorkz[i];
			error += errorkvx[i] * errorkvx[i];
			error += errorkvy[i] * errorkvy[i];
			error += errorkvz[i] * errorkvz[i];
		}
		

		
		double errmax = 0.0;
		for(int i = 0; i < N; ++i){
			ermax = fmax(ermax, fabs(errorkx[i]));
			ermax = fmax(ermax, fabs(errorky[i]));
			ermax = fmax(ermax, fabs(errorkz[i]));
			ermax = fmax(ermax, fabs(errorkvx[i]));
			ermax = fmax(ermax, fabs(errorkvy[i]));
			ermax = fmax(ermax, fabs(errorkvz[i]));
		}
		

		
		double ee = 1.0/5.0;	

		double s = pow( 1.0  / error, ee);
		*/

		//printf("%g %g\n", dt, s);			

		
		for(int i = Nperturbers; i < N; ++i){
			update(x_h, y_h, z_h, vx_h, vy_h, vz_h, kx_h, ky_h, kz_h, kvx_h, kvy_h, kvz_h, i, N, dt, b);	
		}
		update_kernel <<< 1, NN >>> (x_d, y_d, z_d, vx_d, vy_d, vz_d, kx_d, ky_d, kz_d, kvx_d, kvy_d, kvz_d, N, Nperturbers, dt, b);	

		
		time = time0 + t * dt / dayUnit;

		//update table
		if(time - time1 >= dtime){
			int countNodes = Ninterpolate - 1;
			int er;
			for(int j = 0; j < Ninterpolate - 1; ++j){
				for(int i = 0; i < Nperturbers; ++i){
					xp_h[i * Ninterpolate + j] = xp_h[i * Ninterpolate + j + 1];
					yp_h[i * Ninterpolate + j] = yp_h[i * Ninterpolate + j + 1];
					zp_h[i * Ninterpolate + j] = zp_h[i * Ninterpolate + j + 1];
					timep_h[i * Ninterpolate + j] = timep_h[i * Ninterpolate + j + 1];
				}
			}

//printf("CountNodes %d\n", countNodes);
			for(int i = 0; i < Nperturbers; ++i){
				double skip;
				double timepp;
				int id;
				er = fscanf(XVfile, "%lf %d", &timepp, &id);
				fscanf(XVfile, "%lf %lf %lf", &xp_h[id * Ninterpolate + countNodes], &yp_h[id * Ninterpolate + countNodes], &zp_h[id * Ninterpolate + countNodes]);
				//remove velocities read later
				fscanf(XVfile, "%lf %lf %lf", &skip, &skip, &skip);

				//vxp[id * Ninterpolate + countNodes] /= dayUnit;
				//vyp[id * Ninterpolate + countNodes] /= dayUnit;
				//vzp[id * Ninterpolate + countNodes] /= dayUnit;

				if(er < 0) break;
				timep_h[id * Ninterpolate + countNodes] = timepp;
	
//printf("%.20g %d %.20g %.20g %.20g %d\n", timep_h[id * Ninterpolate + countNodes], id, xp_h[id * Ninterpolate + countNodes], yp_h[id * Ninterpolate + countNodes], zp_h[id * Ninterpolate + countNodes], id * Ninterpolate + countNodes);

			}
			if(er < 0){
				printf("Error, time too large, not enough data after time\n");
				return 0;
			}
			time1 = time;
			cudaMemcpy(timep_d, timep_h, Nperturbers * Ninterpolate * sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(xp_d, xp_h, Nperturbers * Ninterpolate * sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(yp_d, yp_h, Nperturbers * Ninterpolate * sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(zp_d, zp_h, Nperturbers * Ninterpolate * sizeof(double), cudaMemcpyHostToDevice);
		}
		
		// ---------------------------------------
		//interpolate on host
		for(int p = 0; p < Nperturbers; ++p){
			interpolate <NN> (Ninterpolate, xp_h, yp_h, zp_h, timep_h, time, x_h, y_h, z_h, p);
		}
		interpolate_kernel < NN > <<< 1, NN >>> (N, Ninterpolate, Nperturbers, xp_d, yp_d, zp_d, timep_d, time, x_d, y_d, z_d);
		// ---------------------------------------

		
		
		if(t % outInterval == 0){

			cudaMemcpy(x_h, x_d, NN * sizeof(double), cudaMemcpyDeviceToHost);
			cudaMemcpy(y_h, y_d, NN * sizeof(double), cudaMemcpyDeviceToHost);
			cudaMemcpy(z_h, z_d, NN * sizeof(double), cudaMemcpyDeviceToHost);
			cudaMemcpy(vx_h, vx_d, NN * sizeof(double), cudaMemcpyDeviceToHost);
			cudaMemcpy(vy_h, vy_d, NN * sizeof(double), cudaMemcpyDeviceToHost);
			cudaMemcpy(vz_h, vz_d, NN * sizeof(double), cudaMemcpyDeviceToHost);
			
			if(outHelio == 1){
				sprintf(outfilename, "Outhelio10_%.12lld.dat", t);
			}
			else{
				sprintf(outfilename, "Outbary10_%.12lld.dat", t);
			}
			outfile = fopen(outfilename, "w");
//			printf("%s\n", outfilename);

			comx = 0.0;
			comy = 0.0;
			comz = 0.0;
			vcomx = 0.0;
			vcomy = 0.0;
			vcomz = 0.0;
			mtot = 0.0;

			if(useHelio == 0 && outHelio == 1){
				//convert to heliocentric output
				comx = -x_h[0];
				comy = -y_h[0];
				comz = -z_h[0];
				vcomx = -vx_h[0];
				vcomy = -vy_h[0];
				vcomz = -vz_h[0];
			}
			
			for(int i = 0; i < N; ++i){
				fprintf(outfile, "%.10g %d %.40g %.40g %.40g %.40g %.40g %.40g %.40g\n", time, i, m_h[i], comx + x_h[i], comy + y_h[i], comz + z_h[i], vcomx + vx_h[i], vcomy + vy_h[i], vcomz + vz_h[i]);
				
			}
			fclose(outfile);
		}
		
	
		/*
		if( s < 1.0){
			dt *= 0.5;
		}
		if( s > 3.0){
			dt *= 2.0;
		}
		*/
	}	// end of time step loop
	fclose(XVfile);
	
}
	

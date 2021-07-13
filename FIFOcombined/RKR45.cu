#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//for FIFO
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "read.h"

#define dayUnit 0.01720209895
//#define dayUnit 0.01720209894846
//#define dayUnit 1.0

//#define def_c 10065.3201686
#define def_c 10065.320121
#define def_AU 149597870700.0           //AU in m

//constant memory
__constant__ double  a_c[6 * 6];
__constant__ double  b_c[6];
__constant__ double  c_c[6];


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

//Neville-Aitken interpolation
// p is perturber index
__host__ void interpolate(int Ninterpolate, double *xp, double *yp, double *zp, double *timep, double time, double *xt, double *yt, double *zt, int p){


	double Px[Ninterpolate][Ninterpolate];
	double Py[Ninterpolate][Ninterpolate];
	double Pz[Ninterpolate][Ninterpolate];
	double tn[Ninterpolate];

	for(int i = 0; i < Ninterpolate; ++i){
		Px[0][i] = xp[p * Ninterpolate + i];
		Py[0][i] = yp[p * Ninterpolate + i];
		Pz[0][i] = zp[p * Ninterpolate + i];
		tn[i] = timep[p * Ninterpolate + i];

//printf("interpolate %d %d %.20g %.20g %.20g\n", p, i, time, tn[i], P[0][i]);
	}

	for(int j = 1; j < Ninterpolate; ++j){
//printf("****\n");
		for(int i = 0; i < Ninterpolate - j; ++i){
			Px[j][i] = ((time - tn[i+j]) * Px[j-1][i] + (tn[i] - time) * Px[j-1][i+1]) / (tn[i] - tn[i+j]);
			Py[j][i] = ((time - tn[i+j]) * Py[j-1][i] + (tn[i] - time) * Py[j-1][i+1]) / (tn[i] - tn[i+j]);
			Pz[j][i] = ((time - tn[i+j]) * Pz[j-1][i] + (tn[i] - time) * Pz[j-1][i+1]) / (tn[i] - tn[i+j]);
//printf("%d %d %g %g %g %g %.20g\n", i, i+j, tn[i], tn[i+j], P[j-1][i], P[j-1][i+1], P[j][i]);

		}
	}
	xt[p] = Px[Ninterpolate-1][0];
	yt[p] = Py[Ninterpolate-1][0];
	zt[p] = Pz[Ninterpolate-1][0];
//printf("interpolate %.20g %d %.20g %.20g %.20g\n", time, p, xt[p], yt[p], zt[p]);

}
__host__ void interpolate2(int Ninterpolate, double *xp, double *yp, double *zp, double *timep, double timep0, double dtimep, double time, double *xt, double *yt, double *zt, int p){


	//p is the perturber index

	double Cx[Ninterpolate];
	double Cy[Ninterpolate];
	double Cz[Ninterpolate];

	double Dx[Ninterpolate];
	double Dy[Ninterpolate];
	double Dz[Ninterpolate];

	double tn[Ninterpolate];

	for(int i = 0; i < Ninterpolate; ++i){
		Cx[i] = xp[p * Ninterpolate + i];
		Cy[i] = yp[p * Ninterpolate + i];
		Cz[i] = zp[p * Ninterpolate + i];

		Dx[i] = Cx[i];		
		Dy[i] = Cy[i];		
		Dz[i] = Cz[i];		

		tn[i] = timep[p * Ninterpolate + i];

//printf("interpolate %d %d %.20g %.20g %.20g\n", p, i, time, tn[i], Cx[i]);
	}

	//initialize with closest solution
	//Assume that the closest solution is in the middle

	int ii = Ninterpolate / 2 - 1;
	xt[p] = Cx[ii];
	yt[p] = Cy[ii];
	zt[p] = Cz[ii];

	--ii;

	for(int j = 1; j < Ninterpolate; ++j){
//printf("**** %d %d %g\n", j, ii, xt[p]);
		for(int i = 0; i < Ninterpolate - j; ++i){

			double dtn0 = tn[i] - time;
			double dtn1 = tn[i + j] - time;
			double dtn = tn[i] - tn[i + j];

			double dPx = (Cx[i + 1] - Dx[i]) / dtn;
			double dPy = (Cy[i + 1] - Dy[i]) / dtn;
			double dPz = (Cz[i + 1] - Dz[i]) / dtn;

			Dx[i] = dtn1 * dPx;
			Dy[i] = dtn1 * dPy;
			Dz[i] = dtn1 * dPz;

			Cx[i] = dtn0 * dPx;
			Cy[i] = dtn0 * dPy;
			Cz[i] = dtn0 * dPz;

	
//printf("%d %d %g %g %g %g %.20g\n", i, i+j, tn[i], tn[i+j], P[j-1][i], P[j-1][i+1], P[j][i]);

		}

		if(2 * ii < Ninterpolate - j){
			xt[p] += Cx[ii + 1];
			yt[p] += Cy[ii + 1];
			zt[p] += Cz[ii + 1];
		}
		else{
			xt[p] += Dx[ii];
			yt[p] += Dy[ii];
			zt[p] += Dz[ii];
			--ii;
		}

	}
//printf("interpolate %.20g %d %.20g %.20g %.20g\n", time, p, xt[p], yt[p], zt[p]);

}


template <int Ninterpolate>
__global__ void interpolate_kernel(int Nperturbers, double *xp_d, double *yp_d, double *zp_d, double *timep_d, double timep0, double dtimep, double time, double *xt_d, double *yt_d, double *zt_d){

	int pid = blockIdx.x;	//perturber index, Nperturbers
	int idx = threadIdx.x;

	if(pid < Nperturbers){

		__shared__ double Px_s[Ninterpolate][Ninterpolate];
		__shared__ double Py_s[Ninterpolate][Ninterpolate];
		__shared__ double Pz_s[Ninterpolate][Ninterpolate];
		__shared__ double tn_s[Ninterpolate];

		int it = floor((time - timep0) / dtimep);
		it -= (Ninterpolate / 2 + 1);

//printf("interpolateA %d %.20g %.20g %d\n", pid, timep0, time, it);

		if(idx < Ninterpolate){
			Px_s[0][idx] = xp_d[Nperturbers * (idx + it) + pid];
			Py_s[0][idx] = yp_d[Nperturbers * (idx + it) + pid];
			Pz_s[0][idx] = zp_d[Nperturbers * (idx + it) + pid];
			tn_s[idx] = timep_d[Nperturbers * (idx + it) + pid];

//if(pid == 1) printf("interpolateB %d %d %.20g %.20g %.20g %.20g\n", pid, idx, timep0, time, tn_s[idx], Px_s[0][idx]);
		}
		__syncthreads();

		for(int j = 1; j < Ninterpolate; ++j){
//if(pid == 1) printf("****\n");
			if(idx < Ninterpolate - j){
				Px_s[j][idx] = ((time - tn_s[idx + j]) * Px_s[j-1][idx] + (tn_s[idx] - time) * Px_s[j-1][idx + 1]) / (tn_s[idx] - tn_s[idx + j]);
				Py_s[j][idx] = ((time - tn_s[idx + j]) * Py_s[j-1][idx] + (tn_s[idx] - time) * Py_s[j-1][idx + 1]) / (tn_s[idx] - tn_s[idx + j]);
				Pz_s[j][idx] = ((time - tn_s[idx + j]) * Pz_s[j-1][idx] + (tn_s[idx] - time) * Pz_s[j-1][idx + 1]) / (tn_s[idx] - tn_s[idx + j]);
//if(pid == 1) printf("%d %d %.20g %.20g %g %g %.20g\n", idx, idx+j, tn_s[idx], tn_s[idx + j], Px_s[j-1][idx], Px_s[j-1][idx + 1], Px_s[j][idx]);
			}
			__syncthreads();
		}

		if(idx == 0){
			xt_d[pid] = Px_s[Ninterpolate-1][0];
			yt_d[pid] = Py_s[Ninterpolate-1][0];
			zt_d[pid] = Pz_s[Ninterpolate-1][0];
//printf("interpolateB %.20g %d %.20g %.20g %.20g\n", time, pid, xt_d[pid], yt_d[pid], zt_d[pid]);
		}

	}
}


//using shared memory
template <int Ninterpolate, int Nperturbers>
__global__ void interpolate2_kernel(double *xp_d, double *yp_d, double *zp_d, double *timep_d, double timep0, double dtimep, double time, double *xt_d, double *yt_d, double *zt_d){

	//every perturber is a thread
	//every particle is a block

	int pid = threadIdx.x;	//perturber index, Nperturbers
	int id = blockIdx.x;

	__shared__ double Cx_s[Ninterpolate][Nperturbers];
	__shared__ double Cy_s[Ninterpolate][Nperturbers];
	__shared__ double Cz_s[Ninterpolate][Nperturbers];

	__shared__ double Dx_s[Ninterpolate][Nperturbers];
	__shared__ double Dy_s[Ninterpolate][Nperturbers];
	__shared__ double Dz_s[Ninterpolate][Nperturbers];

	__shared__ double tn_s[Ninterpolate][Nperturbers];

	if(pid < Nperturbers){

		double x, y, z;

		int it = floor((time - timep0) / dtimep);
		it -= (Ninterpolate / 2 + 1);

//printf("interpolateA %d %.20g %.20g %d\n", pid, timep0, time, it);


		for(int i = 0; i < Ninterpolate; ++i){
			Cx_s[i][pid] = xp_d[Nperturbers * (i + it) + pid];
			Cy_s[i][pid] = yp_d[Nperturbers * (i + it) + pid];
			Cz_s[i][pid] = zp_d[Nperturbers * (i + it) + pid];

			Dx_s[i][pid] = Cx_s[i][pid];
			Dy_s[i][pid] = Cy_s[i][pid];
			Dz_s[i][pid] = Cz_s[i][pid];

			tn_s[i][pid] = timep_d[Nperturbers * (i + it) + pid];

		}

//if(pid == 1) printf("interpolateB %d %d %.20g %.20g %.20g %.20g\n", pid, id, timep0, time, tn_s[0][pid], Cx_s[0][pid]);

		//initialize with closest solution
		//Assume that the closest solution is in the middle

		int ii = Ninterpolate / 2 - 1;
		x = Cx_s[ii][pid];
		y = Cy_s[ii][pid];
		z = Cz_s[ii][pid];

		--ii;

//not yet up to date from here

		for(int j = 1; j < Ninterpolate; ++j){
//if(pid == 1) printf("****\n");
			for(int i = 0; i < Ninterpolate - j; ++i){

				double dtn0 = tn_s[i][pid] - time;
				double dtn1 = tn_s[i + j][pid] - time;
				double dtn = tn_s[i][pid] - tn_s[i + j][pid];

				double dPx = (Cx_s[i + 1][pid] - Dx_s[i][pid]) / dtn;
				double dPy = (Cy_s[i + 1][pid] - Dy_s[i][pid]) / dtn;
				double dPz = (Cz_s[i + 1][pid] - Dz_s[i][pid]) / dtn;

				Dx_s[i][pid] = dtn1 * dPx;
				Dy_s[i][pid] = dtn1 * dPy;
				Dz_s[i][pid] = dtn1 * dPz;

				Cx_s[i][pid] = dtn0 * dPx;
				Cy_s[i][pid] = dtn0 * dPy;
				Cz_s[i][pid] = dtn0 * dPz;
			}

			if(2 * ii < Ninterpolate - j){
				x += Cx_s[ii + 1][pid];
				y += Cy_s[ii + 1][pid];
				z += Cz_s[ii + 1][pid];
			}
			else{
				x += Dx_s[ii][pid];
				y += Dy_s[ii][pid];
				z += Dz_s[ii][pid];
				--ii;
			}


		}

		xt_d[pid] = x;
		yt_d[pid] = y;
		zt_d[pid] = z;
//if(pid == 1) printf("interpolateC %.20g %d %.20g %.20g %.20g\n", time, pid, xt_d[pid], yt_d[pid], zt_d[pid]);

	}
}

//using local memory arrays
template <int Ninterpolate>
__global__ void interpolate2b_kernel(int Nperturbers, double *xp_d, double *yp_d, double *zp_d, double *timep_d, double timep0, double dtimep, double time, double *xt_d, double *yt_d, double *zt_d){

	//every perturber is a thread
	//every particle is a block

	int pid = threadIdx.x;	//perturber index, Nperturbers
	int idx = blockIdx.x;

	if(pid < Nperturbers){

		double x, y, z;

		double Cx[Ninterpolate];
		double Cy[Ninterpolate];
		double Cz[Ninterpolate];

		double Dx[Ninterpolate];
		double Dy[Ninterpolate];
		double Dz[Ninterpolate];

		double tn[Ninterpolate];

		int it = floor((time - timep0) / dtimep);
		it -= (Ninterpolate / 2 + 1);

//printf("interpolateA %d %.20g %.20g %d\n", pid, timep0, time, it);


		for(int i = 0; i < Ninterpolate; ++i){
			Cx[i] = xp_d[Nperturbers * (i + it) + pid];
			Cy[i] = yp_d[Nperturbers * (i + it) + pid];
			Cz[i] = zp_d[Nperturbers * (i + it) + pid];

			Dx[i] = Cx[i];
			Dy[i] = Cy[i];
			Dz[i] = Cz[i];

			tn[i] = timep_d[Nperturbers * (i + it) + pid];

		}

//if(pid == 1) printf("interpolateB %d %d %.20g %.20g %.20g %.20g\n", pid, idx, timep0, time, tn_s[idx], Px_s[0][idx]);

		//initialize with closest solution
		//Assume that the closest solution is in the middle

		int ii = Ninterpolate / 2 - 1;
		x = Cx[ii];
		y = Cy[ii];
		z = Cz[ii];

		--ii;

//not yet up to date from here

		for(int j = 1; j < Ninterpolate; ++j){
//if(pid == 1) printf("****\n");
			for(int i = 0; i < Ninterpolate - j; ++i){

				double dtn0 = tn[i] - time;
				double dtn1 = tn[i + j] - time;
				double dtn = tn[i] - tn[i + j];

				double dPx = (Cx[i + 1] - Dx[i]) / dtn;
				double dPy = (Cy[i + 1] - Dy[i]) / dtn;
				double dPz = (Cz[i + 1] - Dz[i]) / dtn;

				Dx[i] = dtn1 * dPx;
				Dy[i] = dtn1 * dPy;
				Dz[i] = dtn1 * dPz;

				Cx[i] = dtn0 * dPx;
				Cy[i] = dtn0 * dPy;
				Cz[i] = dtn0 * dPz;
			}

			if(2 * ii < Ninterpolate - j){
				x += Cx[ii + 1];
				y += Cy[ii + 1];
				z += Cz[ii + 1];
			}
			else{
				x += Dx[ii];
				y += Dy[ii];
				z += Dz[ii];
				--ii;
			}


		}

		xt_d[pid] = x;
		yt_d[pid] = y;
		zt_d[pid] = z;
//printf("interpolateB %.20g %d %.20g %.20g %.20g\n", time, pid, xt_d[pid], yt_d[pid], zt_d[pid]);

	}
}

//using local memory arrays
template <int Ninterpolate>
__global__ void interpolate3_kernel(int Nperturbers, double *xp_d, double *yp_d, double *zp_d, double *timep_d, double timep0, double dtimep, double time, double *xt_d, double *yt_d, double *zt_d){

	//every perturber is a thread
	//every particle is a block


	for(int pid = 0; pid < Nperturbers; ++pid){

		double x, y, z;

		double Cx[Ninterpolate];
		double Cy[Ninterpolate];
		double Cz[Ninterpolate];

		double Dx[Ninterpolate];
		double Dy[Ninterpolate];
		double Dz[Ninterpolate];

		double tn[Ninterpolate];

		int it = floor((time - timep0) / dtimep);
		it -= (Ninterpolate / 2 + 1);

//printf("interpolateA %d %.20g %.20g %d\n", pid, timep0, time, it);


		for(int i = 0; i < Ninterpolate; ++i){
			Cx[i] = xp_d[Nperturbers * (i + it) + pid];
			Cy[i] = yp_d[Nperturbers * (i + it) + pid];
			Cz[i] = zp_d[Nperturbers * (i + it) + pid];

			Dx[i] = Cx[i];
			Dy[i] = Cy[i];
			Dz[i] = Cz[i];

			tn[i] = timep_d[Nperturbers * (i + it) + pid];

		}

//if(pid == 1) printf("interpolateB %d %.20g %.20g %.20g\n", pid, timep0, time, tn[i]);

		//initialize with closest solution
		//Assume that the closest solution is in the middle

		int ii = Ninterpolate / 2 - 1;
		x = Cx[ii];
		y = Cy[ii];
		z = Cz[ii];

		--ii;

//not yet up to date from here

		for(int j = 1; j < Ninterpolate; ++j){
//if(pid == 1) printf("****\n");
			for(int i = 0; i < Ninterpolate - j; ++i){

				double dtn0 = tn[i] - time;
				double dtn1 = tn[i + j] - time;
				double dtn = tn[i] - tn[i + j];

				double dPx = (Cx[i + 1] - Dx[i]) / dtn;
				double dPy = (Cy[i + 1] - Dy[i]) / dtn;
				double dPz = (Cz[i + 1] - Dz[i]) / dtn;

				Dx[i] = dtn1 * dPx;
				Dy[i] = dtn1 * dPy;
				Dz[i] = dtn1 * dPz;

				Cx[i] = dtn0 * dPx;
				Cy[i] = dtn0 * dPy;
				Cz[i] = dtn0 * dPz;
			}

			if(2 * ii < Ninterpolate - j){
				x += Cx[ii + 1];
				y += Cy[ii + 1];
				z += Cz[ii + 1];
			}
			else{
				x += Dx[ii];
				y += Dy[ii];
				z += Dz[ii];
				--ii;
			}


		}

		xt_d[pid] = x;
		yt_d[pid] = y;
		zt_d[pid] = z;
//printf("interpolateB %.20g %d %.20g %.20g %.20g\n", time, pid, xt_d[pid], yt_d[pid], zt_d[pid]);

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

__host__ void stageStep(double *m, double *xt, double *yt, double *zt, double *vxt, double *vyt, double *vzt, double *kx, double *ky, double *kz, double *kvx, double *kvy, double *kvz, int S, int i, int Nperturbers, int N, int useHelio, int GR){

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
//if(i == 27) printf("Np %d %.20g %.20g %.20g %d\n", i, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit, S);
	}


	if(GR == 2){
		acchGR2(xt[i], yt[i], zt[i], vxt[i], vyt[i], vzt[i], ax, ay, az);
	}

	//NonGrav(xt, yt, zt, vxt, vyt, vzt, ax, ay, az, A1, A2, A3, ALN, NK, NM, Nn, R0, i);
	//J2(m, xt, yt, zt, ax, ay, az, i);

	kvx[i + S * N] = ax;
	kvy[i + S * N] = ay;
	kvz[i + S * N] = az;
}

template < int NN >
__global__ void stageStep_kernel(double *m_d, double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *xt_d, double *yt_d, double *zt_d, double *vxt_d, double *vyt_d, double *vzt_d, double *kx_d, double *ky_d, double *kz_d, double *kvx_d, double *kvy_d, double *kvz_d, double dt, int S, int Nperturbers, int N, int useHelio, int GR){

	int itx = threadIdx.x;
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	//shared memory contains only the perturbers
	//the particle idx is stored in registers
	__shared__ double x_s[NN];
	__shared__ double y_s[NN];
	__shared__ double z_s[NN];
	__shared__ double m_s[NN];


	if(threadIdx.x < Nperturbers){
			x_s[threadIdx.x] = xt_d[itx];
			y_s[threadIdx.x] = yt_d[itx];
			z_s[threadIdx.x] = zt_d[itx];
			m_s[threadIdx.x] = m_d[itx];
	}
	__syncthreads();

	if(idx >= Nperturbers && idx < N){


// ***********************
		//update
		double xi = x_d[idx];
		double yi = y_d[idx];
		double zi = z_d[idx];
		double vxi = vx_d[idx];
		double vyi = vy_d[idx];
		double vzi = vz_d[idx];
		double mi = m_d[idx];

		for(int s = 0; s < S; ++s){
			double aa = a_c[S * 6 + s];
			xi  += dt * aa * kx_d[idx + s * N];
			yi  += dt * aa * ky_d[idx + s * N];
			zi  += dt * aa * kz_d[idx + s * N];
			vxi += dt * aa * kvx_d[idx + s * N];
			vyi += dt * aa * kvy_d[idx + s * N];
			vzi += dt * aa * kvz_d[idx + s * N];
		}
// *****************************

		kx_d[idx + S * N] = vxi;
		ky_d[idx + S * N] = vyi;
		kz_d[idx + S * N] = vzi;

		double ax = 0.0;
		double ay = 0.0;
		double az = 0.0;

		if(useHelio == 0){
			for(int j = Nperturbers-1; j >= 0; --j){
				accP_device(m_s[j], x_s[j], y_s[j], z_s[j], xi, yi, zi, ax, ay, az);
			}
		}
		else{
			for(int j = Nperturbers-1; j >= 1; --j){
				accP_device(m_s[j], x_s[j], y_s[j], z_s[j], xi, yi, zi, ax, ay, az);
//printf("Nij %d %d %.20g %.20g %.20g\n", idx, j, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit);
			}
			accS_device(m_s[0] + mi, xi, yi, zi, ax, ay, az);

//printf("N0 %d %.20g %.20g %.20g\n", idx, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit);
			for(int j = Nperturbers-1; j >= 1; --j){
				accP2_device(m_s[j], x_s[j], y_s[j], z_s[j], ax, ay, az);
//printf("Npi %d %d %.20g %.20g %.20g %d\n", idx, j, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit, S);
			}
//printf("Np %d %.20g %.20g %.20g %d\n", idx, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit, S);
		}


		if(GR == 2){
			acchGR2(xi, yi, zi, vxi, vyi, vzi, ax, ay, az);
		}

		//NonGrav(xt, yt, zt, vxt, vyt, vzt, ax, ay, az, A1, A2, A3, ALN, NK, NM, Nn, R0, idx);
		//J2(m, xt, yt, zt, ax, ay, az, idx);

		kvx_d[idx + S * N] = ax;
		kvy_d[idx + S * N] = ay;
		kvz_d[idx + S * N] = az;
	}

}

//combine stage step and interpolation kernel
template < const int NN, const int Ninterpolate >
__global__ void stageStep1_kernel(double *m_d, double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *xp_d, double *yp_d, double *zp_d, double *timep_d, double timep0, double dtimep, double time, double *kx_d, double *ky_d, double *kz_d, double *kvx_d, double *kvy_d, double *kvz_d, double dt, int S, int Nperturbers, int N, int useHelio, int GR){

	int idx = threadIdx.x;
	int id = blockIdx.x * blockDim.x + idx;

	//shared memory contains only the perturbers
	//the particle idx is stored in registers
	__shared__ double x_s[NN];
	__shared__ double y_s[NN];
	__shared__ double z_s[NN];
	__shared__ double m_s[NN];

	
	//needs at least Nperturbers threads per block
	// *******************************************************
	//interpolate

	//do interpolation in every thread block
	//p is the particle index

	double Cx[Ninterpolate];
	double Cy[Ninterpolate];
	double Cz[Ninterpolate];

	double Dx[Ninterpolate];
	double Dy[Ninterpolate];
	double Dz[Ninterpolate];

	double tn[Ninterpolate];

	int it = floor((time - timep0) / dtimep);
	it -= (Ninterpolate / 2 + 1);

	if(idx < Nperturbers){

		for(int i = 0; i < Ninterpolate; ++i){
			Cx[i] = xp_d[Nperturbers * (i + it) + idx];
			Cy[i] = yp_d[Nperturbers * (i + it) + idx];
			Cz[i] = zp_d[Nperturbers * (i + it) + idx];

			Dx[i] = Cx[i];		
			Dy[i] = Cy[i];		
			Dz[i] = Cz[i];		

			tn[i] = timep_d[Nperturbers * (i + it) + idx];

	//printf("interpolate %d %d %.20g %.20g %.20g\n", idx, i, time, tn[i], Cx[i]);
		}

		//initialize with closest solution
		//Assume that the closest solution is in the middle

		int ii = Ninterpolate / 2 - 1;
		x_s[idx] = Cx[ii];
		y_s[idx] = Cy[ii];
		z_s[idx] = Cz[ii];
		m_s[idx] = m_d[idx]; 

		--ii;

		for(int j = 1; j < Ninterpolate; ++j){
	//printf("**** %d %d %g\n", j, ii, x_s[idx]);
			for(int i = 0; i < Ninterpolate - j; ++i){

				double dtn0 = tn[i] - time;
				double dtn1 = tn[i + j] - time;
				double dtn = tn[i] - tn[i + j];

				double dPx = (Cx[i + 1] - Dx[i]) / dtn;
				double dPy = (Cy[i + 1] - Dy[i]) / dtn;
				double dPz = (Cz[i + 1] - Dz[i]) / dtn;

				Dx[i] = dtn1 * dPx;
				Dy[i] = dtn1 * dPy;
				Dz[i] = dtn1 * dPz;

				Cx[i] = dtn0 * dPx;
				Cy[i] = dtn0 * dPy;
				Cz[i] = dtn0 * dPz;

		
	//printf("%d %d %g %g %g %g %.20g\n", i, i+j, tn[i], tn[i+j], P[j-1][i], P[j-1][i+1], P[j][i]);

			}

			if(2 * ii < Ninterpolate - j){
				x_s[idx] += Cx[ii + 1];
				y_s[idx] += Cy[ii + 1];
				z_s[idx] += Cz[ii + 1];
			}
			else{
				x_s[idx] += Dx[ii];
				y_s[idx] += Dy[ii];
				z_s[idx] += Dz[ii];
				--ii;
			}

		}

	}

	if(id >= Nperturbers && id < N){


// ***********************
		//update
		double xi = x_d[id];
		double yi = y_d[id];
		double zi = z_d[id];
		double vxi = vx_d[id];
		double vyi = vy_d[id];
		double vzi = vz_d[id];
		double mi = m_d[id];

		for(int s = 0; s < S; ++s){
			double aa = a_c[S * 6 + s];
			xi  += dt * aa * kx_d[id + s * N];
			yi  += dt * aa * ky_d[id + s * N];
			zi  += dt * aa * kz_d[id + s * N];
			vxi += dt * aa * kvx_d[id + s * N];
			vyi += dt * aa * kvy_d[id + s * N];
			vzi += dt * aa * kvz_d[id + s * N];
		}
// *****************************

		kx_d[id + S * N] = vxi;
		ky_d[id + S * N] = vyi;
		kz_d[id + S * N] = vzi;

		double ax = 0.0;
		double ay = 0.0;
		double az = 0.0;

		if(useHelio == 0){
			for(int j = Nperturbers-1; j >= 0; --j){
				accP_device(m_s[j], x_s[j], y_s[j], z_s[j], xi, yi, zi, ax, ay, az);
			}
		}
		else{
			for(int j = Nperturbers-1; j >= 1; --j){
				accP_device(m_s[j], x_s[j], y_s[j], z_s[j], xi, yi, zi, ax, ay, az);
//printf("Nij %d %d %.20g %.20g %.20g\n", id, j, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit);
			}
			accS_device(m_s[0] + mi, xi, yi, zi, ax, ay, az);

//printf("N0 %d %.20g %.20g %.20g\n", id, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit);
			for(int j = Nperturbers-1; j >= 1; --j){
				accP2_device(m_s[j], x_s[j], y_s[j], z_s[j], ax, ay, az);
//printf("Npi %d %d %.20g %.20g %.20g %d\n", id, j, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit, S);
			}
//printf("Np %d %.20g %.20g %.20g %d\n", id, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit, S);
		}


		if(GR == 2){
			acchGR2(xi, yi, zi, vxi, vyi, vzi, ax, ay, az);
		}

		//NonGrav(xt, yt, zt, vxt, vyt, vzt, ax, ay, az, A1, A2, A3, ALN, NK, NM, Nn, R0, id);
		//J2(m, xt, yt, zt, ax, ay, az, id);

		kvx_d[id + S * N] = ax;
		kvy_d[id + S * N] = ay;
		kvz_d[id + S * N] = az;
	}

}

//combine stage step and interpolation kernel
//all stages in the same kernel
template < const int Nperturbers, const int Ninterpolate, const int nn >
__global__ void stageStepAll_kernel(double *m_d, double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *xp_d, double *yp_d, double *zp_d, double *timep_d, double time0, double dt, int N, int useHelio, int GR){

	int idx = threadIdx.x;
	int id = blockIdx.x * blockDim.x + idx;

	//shared memory contains only the perturbers
	//the particle idx is stored in registers
	__shared__ double x_s[Nperturbers];
	__shared__ double y_s[Nperturbers];
	__shared__ double z_s[Nperturbers];
	__shared__ double m_s[Nperturbers];


	__shared__ double kx_s[nn * 6];
	__shared__ double ky_s[nn * 6];
	__shared__ double kz_s[nn * 6];

	__shared__ double kvx_s[nn * 6];
	__shared__ double kvy_s[nn * 6];
	__shared__ double kvz_s[nn * 6];
	
	//needs at least Nperturbers threads per block
	// *******************************************************
	//interpolate

	//do interpolation in every thread block
	//p is the particle index

	double Cx[Ninterpolate];
	double Cy[Ninterpolate];
	double Cz[Ninterpolate];

	double Dx[Ninterpolate];
	double Dy[Ninterpolate];
	double Dz[Ninterpolate];

	double tn[Ninterpolate];


	for(int S = 0; S < 6; ++S){

		double time = time0 + c_c[S] * dt / dayUnit;
		// **********************************3
		//interpolation
		if(idx < Nperturbers){

			for(int i = 0; i < Ninterpolate; ++i){
				Cx[i] = xp_d[idx * Ninterpolate + i];
				Cy[i] = yp_d[idx * Ninterpolate + i];
				Cz[i] = zp_d[idx * Ninterpolate + i];

				Dx[i] = Cx[i];		
				Dy[i] = Cy[i];		
				Dz[i] = Cz[i];		

				tn[i] = timep_d[idx * Ninterpolate + i];

		//printf("interpolate %d %d %.20g %.20g %.20g\n", idx, i, time, tn[i], Cx[i]);
			}

			//initialize with closest solution
			//Assume that the closest solution is in the middle

			int ii = Ninterpolate / 2 - 1;
			x_s[idx] = Cx[ii];
			y_s[idx] = Cy[ii];
			z_s[idx] = Cz[ii];
			m_s[idx] = m_d[idx]; 

			--ii;

			for(int j = 1; j < Ninterpolate; ++j){
		//printf("**** %d %d %g\n", j, ii, x_s[idx]);
				for(int i = 0; i < Ninterpolate - j; ++i){

					double dtn0 = tn[i] - time;
					double dtn1 = tn[i + j] - time;
					double dtn = tn[i] - tn[i + j];

					double dPx = (Cx[i + 1] - Dx[i]) / dtn;
					double dPy = (Cy[i + 1] - Dy[i]) / dtn;
					double dPz = (Cz[i + 1] - Dz[i]) / dtn;

					Dx[i] = dtn1 * dPx;
					Dy[i] = dtn1 * dPy;
					Dz[i] = dtn1 * dPz;

					Cx[i] = dtn0 * dPx;
					Cy[i] = dtn0 * dPy;
					Cz[i] = dtn0 * dPz;

			
		//printf("%d %d %g %g %g %g %.20g\n", i, i+j, tn[i], tn[i+j], P[j-1][i], P[j-1][i+1], P[j][i]);

				}

				if(2 * ii < Ninterpolate - j){
					x_s[idx] += Cx[ii + 1];
					y_s[idx] += Cy[ii + 1];
					z_s[idx] += Cz[ii + 1];
				}
				else{
					x_s[idx] += Dx[ii];
					y_s[idx] += Dy[ii];
					z_s[idx] += Dz[ii];
					--ii;
				}

			}

		}
		// end interpolation **********************************************

		__syncthreads();

		if(id >= Nperturbers && id < N){

			// ***********************
			//update
			double xi = x_d[id];
			double yi = y_d[id];
			double zi = z_d[id];
			double vxi = vx_d[id];
			double vyi = vy_d[id];
			double vzi = vz_d[id];
			double mi = m_d[id];

			for(int s = 0; s < S; ++s){
				double aa = a_c[S * 6 + s];
				xi  += dt * aa * kx_s[idx + s * nn];
				yi  += dt * aa * ky_s[idx + s * nn];
				zi  += dt * aa * kz_s[idx + s * nn];
				vxi += dt * aa * kvx_s[idx + s * nn];
				vyi += dt * aa * kvy_s[idx + s * nn];
				vzi += dt * aa * kvz_s[idx + s * nn];
			}
			// *****************************

			kx_s[idx + S * nn] = vxi;
			ky_s[idx + S * nn] = vyi;
			kz_s[idx + S * nn] = vzi;

			double ax = 0.0;
			double ay = 0.0;
			double az = 0.0;

			if(useHelio == 0){
				for(int j = Nperturbers-1; j >= 0; --j){
					accP_device(m_s[j], x_s[j], y_s[j], z_s[j], xi, yi, zi, ax, ay, az);
				}
			}
			else{
				for(int j = Nperturbers-1; j >= 1; --j){
					accP_device(m_s[j], x_s[j], y_s[j], z_s[j], xi, yi, zi, ax, ay, az);
	//printf("Nij %d %d %.20g %.20g %.20g\n", id, j, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit);
				}
				accS_device(m_s[0] + mi, xi, yi, zi, ax, ay, az);

	//printf("N0 %d %.20g %.20g %.20g\n", id, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit);
				for(int j = Nperturbers-1; j >= 1; --j){
					accP2_device(m_s[j], x_s[j], y_s[j], z_s[j], ax, ay, az);
	//printf("Npi %d %d %.20g %.20g %.20g %d\n", id, j, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit, S);
				}
	//printf("Np %d %.20g %.20g %.20g %d\n", id, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit, S);
			}


			if(GR == 2){
				acchGR2(xi, yi, zi, vxi, vyi, vzi, ax, ay, az);
			}

			//NonGrav(xt, yt, zt, vxt, vyt, vzt, ax, ay, az, A1, A2, A3, ALN, NK, NM, Nn, R0, id);
			//J2(m, xt, yt, zt, ax, ay, az, id);

			kvx_s[idx + S * nn] = ax;
			kvy_s[idx + S * nn] = ay;
			kvz_s[idx + S * nn] = az;
		}
		// end of stage
		__syncthreads();
	}

	if(id >= Nperturbers && id < N){
	//RKF45
		x_d[id] += dt * (b_c[0] * kx_s[idx + 0 * nn] + b_c[2] * kx_s[idx + 2 * nn] + b_c[3] * kx_s[idx + 3 * nn] + b_c[4] * kx_s[idx + 4 * nn] + b_c[5] * kx_s[idx + 5 * nn]);
		y_d[id] += dt * (b_c[0] * ky_s[idx + 0 * nn] + b_c[2] * ky_s[idx + 2 * nn] + b_c[3] * ky_s[idx + 3 * nn] + b_c[4] * ky_s[idx + 4 * nn] + b_c[5] * ky_s[idx + 5 * nn]);
		z_d[id] += dt * (b_c[0] * kz_s[idx + 0 * nn] + b_c[2] * kz_s[idx + 2 * nn] + b_c[3] * kz_s[idx + 3 * nn] + b_c[4] * kz_s[idx + 4 * nn] + b_c[5] * kz_s[idx + 5 * nn]);

		vx_d[id] += dt * (b_c[0] * kvx_s[idx + 0 * nn] + b_c[2] * kvx_s[idx + 2 * nn] + b_c[3] * kvx_s[idx + 3 * nn] + b_c[4] * kvx_s[idx + 4 * nn] + b_c[5] * kvx_s[idx + 5 * nn]);
		vy_d[id] += dt * (b_c[0] * kvy_s[idx + 0 * nn] + b_c[2] * kvy_s[idx + 2 * nn] + b_c[3] * kvy_s[idx + 3 * nn] + b_c[4] * kvy_s[idx + 4 * nn] + b_c[5] * kvy_s[idx + 5 * nn]);
		vz_d[id] += dt * (b_c[0] * kvz_s[idx + 0 * nn] + b_c[2] * kvz_s[idx + 2 * nn] + b_c[3] * kvz_s[idx + 3 * nn] + b_c[4] * kvz_s[idx + 4 * nn] + b_c[5] * kvz_s[idx + 5 * nn]);

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


__host__ void update2(double *xt, double *yt, double *zt, double *vxt, double *vyt, double *vzt, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *kx, double *ky, double *kz, double *kvx, double *kvy, double *kvz, int i, int N, double dt, int S, double *a){

	xt[i]  = x[i];
	yt[i]  = y[i];
	zt[i]  = z[i];
	vxt[i] = vx[i];
	vyt[i] = vy[i];
	vzt[i] = vz[i];

	for(int s = 0; s < S; ++s){
		xt[i]  += dt * a[S * 6 + s] * kx[i + s * N];
		yt[i]  += dt * a[S * 6 + s] * ky[i + s * N];
		zt[i]  += dt * a[S * 6 + s] * kz[i + s * N];
		vxt[i] += dt * a[S * 6 + s] * kvx[i + s * N];
		vyt[i] += dt * a[S * 6 + s] * kvy[i + s * N];
		vzt[i] += dt * a[S * 6 + s] * kvz[i + s * N];
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

__global__ void update_kernel(double *x, double *y, double *z, double *vx, double *vy, double *vz, double *kx, double *ky, double *kz, double *kvx, double *kvy, double *kvz, int N, int Nperturbers, double dt){
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if(idx >= Nperturbers && idx < N){
	//RKF45
		x[idx] += dt * (b_c[0] * kx[idx + 0 * N] + b_c[2] * kx[idx + 2 * N] + b_c[3] * kx[idx + 3 * N] + b_c[4] * kx[idx + 4 * N] + b_c[5] * kx[idx + 5 * N]);
		y[idx] += dt * (b_c[0] * ky[idx + 0 * N] + b_c[2] * ky[idx + 2 * N] + b_c[3] * ky[idx + 3 * N] + b_c[4] * ky[idx + 4 * N] + b_c[5] * ky[idx + 5 * N]);
		z[idx] += dt * (b_c[0] * kz[idx + 0 * N] + b_c[2] * kz[idx + 2 * N] + b_c[3] * kz[idx + 3 * N] + b_c[4] * kz[idx + 4 * N] + b_c[5] * kz[idx + 5 * N]);

		vx[idx] += dt * (b_c[0] * kvx[idx + 0 * N] + b_c[2] * kvx[idx + 2 * N] + b_c[3] * kvx[idx + 3 * N] + b_c[4] * kvx[idx + 4 * N] + b_c[5] * kvx[idx + 5 * N]);
		vy[idx] += dt * (b_c[0] * kvy[idx + 0 * N] + b_c[2] * kvy[idx + 2 * N] + b_c[3] * kvy[idx + 3 * N] + b_c[4] * kvy[idx + 4 * N] + b_c[5] * kvy[idx + 5 * N]);
		vz[idx] += dt * (b_c[0] * kvz[idx + 0 * N] + b_c[2] * kvz[idx + 2 * N] + b_c[3] * kvz[idx + 3 * N] + b_c[4] * kvz[idx + 4 * N] + b_c[5] * kvz[idx + 5 * N]);
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
	const int Ninterpolate = 10;	//number of interpolation points
	const double dtimep = 1.0;     //interval between stored time steps
	double time0 = 0.0;		//start time from simulation
	double time1 = 0.0;		//end time from simulation
	double outStart = 0;		//start time of output files

	int comet = 0;

	int GR = 2;
	//2 Sitarski 1982, heliocentric coordinates

	int useFIFO = 2;	//use 0 or 2
	int useGPU = 1;

	FILE *binfile;
	if(useFIFO == 2){	
		binfile = fopen("210705_2315_genga_req.bin", "rb");
	}	

	int useHelio = 1;
	int outHelio = 1;
	//1 print output in heliocentric coordinates
	//0 print output in barycentric  coordinates

	int OutBinary = 1;

	//long long int Nsteps = 40000;	
	long long int Nsteps = 18000;
	long long int outInterval = 10;
	double dt = 0.1 * dayUnit;

	// short comparison
	//long long int Nsteps = 40000;
	//long long int Nsteps = 400;
	//long long int outInterval = 100;
	//double dt = 0.01 * dayUnit;
	
	//for timing
	//long long int Nsteps = 1000;
	//long long int outInterval = 1000;
	//double dt = 0.01 * dayUnit;

	time1 = time0 + dt * Nsteps;


	if(useFIFO == 2){
		printf("read file\n");
		readHeader(binfile, time0, time1, outInterval, outStart, NTP, comet);

		//change this later to be more general
		double ddt = 0.1;
		outInterval /= ddt;
		dt = ddt * dayUnit;
		Nsteps = (time1 - time0) / ddt;

		printf("Nsteps: %lld %.20g %.20g\n", Nsteps, time1, time0);

		//NTP = 1;
		//Nsteps = 18000;
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
		}
		else if(strcmp(argv[i], "-N") == 0){
			NTP = atoi(argv[i + 1]);
		}
	}
	int NN = Nperturbers + NTP;//+ 8192 * 4; //28

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

	cudaMalloc((void **) &timep_d, Nperturbers * NTable * sizeof(double));
	cudaMalloc((void **) &xp_d, Nperturbers * NTable * sizeof(double));
	cudaMalloc((void **) &yp_d, Nperturbers * NTable * sizeof(double));
	cudaMalloc((void **) &zp_d, Nperturbers * NTable * sizeof(double));

	
	//non Gravitational constants
	double A1_h[NN];
	double A2_h[NN];
	double A3_h[NN];

	for(int i = 0; i < NN; ++i){
		A1_h[i] = 0.0;
		A2_h[i] = 0.0;
		A3_h[i] = 0.0;
	}
	

	FILE *outfile;
	char outfilename[160];	


	double time = 0.0;

	/*
	//Units are 1/(mass of object in solar masses)
	double pmass[] = {
		//The sun has to be at the first position
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
	*/

	double mass1[] = {
		//masses form DE440
		//in units of GM
		1.3271244004127942E+11,	//Sun		10
		2.2031868551400003E+04,	//Mercury	1
		3.2485859200000000E+05,	//Venus		2
		3.9860043550702266E+05,	//Earth		399
		4.2828375815756102E+04,	//Mars		4
		1.2671276409999998E+08,	//Jupiter	5
		3.7940584841799997E+07,	//Saturn	6
		5.7945563999999985E+06,	//Uranus	7
		6.8365271005803989E+06,	//Neptune	8
		9.7550000000000000E+02,	//Pluto		9
		4.9028001184575496E+03,	//Moon		301

		6.2628888644409933E+01,	//Ceres		1
		1.3665878145967422E+01,	//Pallas	2
		1.9205707002025889E+00,	//Juno		3
		1.7288232879171513E+01,	//Vesta		4
		1.1398723232184107E+00,	//Iris		7
		5.6251476453852289E+00,	//Hygiea	10
		2.0230209871098284E+00,	//Eunomia	15
		1.5896582441709424E+00,	//Psyche	16
		1.0793714577033560E+00,	//Euphrosyne	31
		2.6830359242821795E+00,	//Europa	52
		9.3810575639151328E-01,	//Cybele	65
		2.1682320736996910E+00,	//Sylvia	87
		1.1898077088121908E+00,	//Thisbe	88
		1.4437384031866001E+00,	//Camilla	107
		3.8944831481705644E+00,	//Davida	511
		2.8304096393299849E+00,	//Interamnia	704
	};
	for(int i = 0; i < Nperturbers; ++i){
		m_h[i] = mass1[i] / mass1[0];
printf("m %d %.20g\n", i, m_h[i]);
	}


//m[Nperturbers] = 1.e-11; //ca mass of Flora

	int N = Nperturbers;

	//Sun
	id_h[0] = 20;
	x_h[0] = 0.0;
	y_h[0] = 0.0;
	z_h[0] = 0.0;
	vx_h[0] = 0.0;
	vy_h[0] = 0.0;
	vz_h[0] = 0.0;

	for(int i = 1; i < NN; ++i){
		id_h[i] = Nperturbers + i;
		x_h[i] = 0.0;
		y_h[i] = 0.0;
		z_h[i] = 0.0;
		vx_h[i] = 0.0;
		vy_h[i] = 0.0;
		vz_h[i] = 0.0;
	}
	printf("Read initial conditions\n");

	if(useFIFO == 2){	
		//read test particles
		int er = 0;
		er = readFile(binfile, Nperturbers, x_h, y_h, z_h, vx_h, vy_h, vz_h, A1_h, A2_h, A3_h, id_h, NTP);
		fclose(binfile);
		N += NTP;
		if(er == 1) return 0;
	}
	else{

		FILE *infile;
		char infilename[160];

		//sprintf(infilename, "initial1.dat");
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
	printf("er %d %llu %d %d %.20g %.20g %.20g\n", i, id_h[i], er, N, x_h[i], y_h[i], z_h[i]);
		}
		fclose(infile);
		time0 = time;	//start time from simulation
	}

	double timep0 = 0.0;	//start time from perturbers file

	for(int i = Nperturbers; i < N; ++i){
		vx_h[i] /= dayUnit;
		vy_h[i] /= dayUnit;
		vz_h[i] /= dayUnit;

		//A1[i] /= (dayUnit * dayUnit);
		//A2[i] /= (dayUnit * dayUnit);
		//A3[i] /= (dayUnit * dayUnit);
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
				unsigned long long int id = id_h[i];
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
		cudaMemcpy(x_d, x_h, NN * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(y_d, y_h, NN * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(z_d, z_h, NN * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(vx_d, vx_h, NN * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(vy_d, vy_h, NN * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(vz_d, vz_h, NN * sizeof(double), cudaMemcpyHostToDevice);
	}
	

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

	
	double *a_h;
	a_h = (double*)malloc(6 * 6 * sizeof(double));

	a_h[1 * 6 + 0] = 1.0/4.0;	//21

	a_h[2 * 6 + 0] = 3.0/32.0;	//31
	a_h[2 * 6 + 1] = 9.0/32.0;	//32

	a_h[3 * 6 + 0] = 1932.0/2197.0;	//41
	a_h[3 * 6 + 1] = -7200.0/2197.0;//42
	a_h[3 * 6 + 2] = 7296.0/2197.0;	//43

	a_h[4 * 6 + 0] = 439.0/216.0;	//51
	a_h[4 * 6 + 1] = -8.0;		//52
	a_h[4 * 6 + 2] = 3680.0/513.0;	//53
	a_h[4 * 6 + 3] = -845.0/4104.0;	//54

	a_h[5 * 6 + 0] = -8.0/27.0;	//61
	a_h[5 * 6 + 1] = 2.0;		//62
	a_h[5 * 6 + 2] = -3544/2565.0;	//63
	a_h[5 * 6 + 3] = 1859.0/4104.0;	//64
	a_h[5 * 6 + 4] = -11.0/40.0;	//65

	double *b_h, *bb_h;
	b_h = (double*)malloc(6 * sizeof(double));
	bb_h = (double*)malloc(6 * sizeof(double));

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

	double *c_h;
	c_h = (double*)malloc(6 * sizeof(double));

	c_h[0] = 0.0;
	c_h[1] = 0.25;
	c_h[2] = 3.0 / 8.0;
	c_h[3] = 12.0 / 13.0;
	c_h[4] = 1.0;
	c_h[5] = 0.5;

	if(useGPU == 1){
		cudaMemcpyToSymbol(a_c, a_h, 6 * 6 * sizeof(double), 0, cudaMemcpyHostToDevice);
		cudaMemcpyToSymbol(b_c, b_h, 6 * sizeof(double), 0, cudaMemcpyHostToDevice);
		cudaMemcpyToSymbol(c_c, c_h, 6 * sizeof(double), 0, cudaMemcpyHostToDevice);
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

	int S;

	cudaEvent_t LoopStart, LoopStop;
	cudaEventCreate(&LoopStart);
	cudaEventCreate(&LoopStop);

	cudaEventRecord(LoopStart);

	time = time0;
	for(long long int t = 1; t <= Nsteps; ++t){

//cudaDeviceSynchronize();
//printf("%lld %d | %d %d\n", t, 0, NTP, Nperturbers);	
		//stage 1
		S = 0;
		if(useGPU == 0){
			for(int p = 0; p < Nperturbers; ++p){
				//interpolate(Ninterpolate, xp_h, yp_h, zp_h, timep_h, time + c_h[0] * dt / dayUnit, xt_h, yt_h, zt_h, p);
				interpolate2(Ninterpolate, xp_h, yp_h, zp_h, timep_h, timep0, dtimep, time + c_h[0] * dt / dayUnit, xt_h, yt_h, zt_h, p);
			}
			for(int i = 0; i < N; ++i){
				update1(xt_h, yt_h, zt_h, vxt_h, vyt_h, vzt_h, x_h, y_h, z_h, vx_h, vy_h, vz_h, i);
			}
			for(int i = Nperturbers; i < N; ++i){
				stageStep(m_h, xt_h, yt_h, zt_h, vxt_h, vyt_h, vzt_h, kx_h, ky_h, kz_h, kvx_h, kvy_h, kvz_h, S, i, Nperturbers, N, useHelio, GR);
			}
		}
		else{
			interpolate_kernel < Ninterpolate > <<< dim3(Nperturbers, NTP, 1), Ninterpolate >>> (Nperturbers, xp_d, yp_d, zp_d, timep_d, timep0, dtimep, time + c_h[0] * dt / dayUnit, xt_d, yt_d, zt_d);
			//interpolate2_kernel < Ninterpolate, Nperturbers > <<< NTP, Nperturbers >>> (xp_d, yp_d, zp_d, timep_d, timep0, dtimep, time + c_h[0] * dt / dayUnit, xt_d, yt_d, zt_d);
			//interpolate2b_kernel < Ninterpolate > <<< NTP, Nperturbers >>> (Nperturbers, xp_d, yp_d, zp_d, timep_d, timep0, dtimep, time + c_h[0] * dt / dayUnit, xt_d, yt_d, zt_d);
			//interpolate3_kernel < Ninterpolate > <<< (NTP + 127) / 128, 128 >>> (Nperturbers, xp_d, yp_d, zp_d, timep_d, timep0, dtimep, time + c_h[0] * dt / dayUnit, xt_d, yt_d, zt_d);
			stageStep_kernel < Nperturbers > <<< (NN + 127) / 128, 128 >>> (m_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, xt_d, yt_d, zt_d, vxt_d, vyt_d, vzt_d, kx_d, ky_d, kz_d, kvx_d, kvy_d, kvz_d, dt, S, Nperturbers, N, useHelio, GR);
			//stageStep1_kernel < Nperturbers, Ninterpolate > <<< (NN + 127) / 128, 128 >>> (m_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, xp_d, yp_d, zp_d, timep_d, timep0, dtimep, time + c_h[S] * dt / dayUnit, kx_d, ky_d, kz_d, kvx_d, kvy_d, kvz_d, dt, S, Nperturbers, N, useHelio, GR);
		}	
	
		//stage 2 - 6
		for(int S = 1; S < 6; ++S){
//cudaDeviceSynchronize();
//printf("%lld %d\n", t, S);	
			if(useGPU == 0){
				for(int p = 0; p < Nperturbers; ++p){
					//interpolate(Ninterpolate, xp_h, yp_h, zp_h, timep_h, time + c_h[S] * dt / dayUnit, xt_h, yt_h, zt_h, p);
					interpolate2(Ninterpolate, xp_h, yp_h, zp_h, timep_h, timep0, dtimep, time + c_h[S] * dt / dayUnit, xt_h, yt_h, zt_h, p);
				}
				for(int i = Nperturbers; i < N; ++i){
					update2(xt_h, yt_h, zt_h, vxt_h, vyt_h, vzt_h, x_h, y_h, z_h, vx_h, vy_h, vz_h, kx_h, ky_h, kz_h, kvx_h, kvy_h, kvz_h, i, N, dt, S, a_h);	//a21
				}
				for(int i = Nperturbers; i < N; ++i){
					stageStep(m_h, xt_h, yt_h, zt_h, vxt_h, vyt_h, vzt_h, kx_h, ky_h, kz_h, kvx_h, kvy_h, kvz_h, S, i, Nperturbers, N, useHelio, GR);
				}
			}
			else{
				interpolate_kernel < Ninterpolate > <<< dim3(Nperturbers, NTP,1), Ninterpolate >>> (Nperturbers, xp_d, yp_d, zp_d, timep_d, timep0, dtimep, time + c_h[S] * dt / dayUnit, xt_d, yt_d, zt_d);
				//interpolate2_kernel < Ninterpolate, Nperturbers > <<< NTP, Nperturbers >>> (xp_d, yp_d, zp_d, timep_d, timep0, dtimep, time + c_h[S] * dt / dayUnit, xt_d, yt_d, zt_d);
				//interpolate2b_kernel < Ninterpolate > <<< NTP, Nperturbers >>> (Nperturbers, xp_d, yp_d, zp_d, timep_d, timep0, dtimep, time + c_h[S] * dt / dayUnit, xt_d, yt_d, zt_d);
				//interpolate3_kernel < Ninterpolate > <<< (NTP + 127) / 128, 128 >>> (Nperturbers, xp_d, yp_d, zp_d, timep_d, timep0, dtimep, time + c_h[S] * dt / dayUnit, xt_d, yt_d, zt_d);
				stageStep_kernel < Nperturbers > <<< (NN + 127) / 128, 128 >>> (m_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, xt_d, yt_d, zt_d, vxt_d, vyt_d, vzt_d, kx_d, ky_d, kz_d, kvx_d, kvy_d, kvz_d, dt, S, Nperturbers, N, useHelio, GR);
				//stageStep1_kernel < Nperturbers, Ninterpolate > <<< (NN + 127) / 128, 128 >>> (m_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, xp_d, yp_d, zp_d, timep_d, timep0, dtimep, time + c_h[S] * dt / dayUnit, kx_d, ky_d, kz_d, kvx_d, kvy_d, kvz_d, dt, S, Nperturbers, N, useHelio, GR);

			}
		}
cudaDeviceSynchronize();
//		stageStepAll_kernel < Nperturbers, Ninterpolate, 128 > <<< (NN + 127) / 128, 128 >>> (m_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, xp_d, yp_d, zp_d, timep_d, time, dt, N, useHelio, GR);
	
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

		
		if(useGPU == 0){
			for(int i = Nperturbers; i < N; ++i){
				update(x_h, y_h, z_h, vx_h, vy_h, vz_h, kx_h, ky_h, kz_h, kvx_h, kvy_h, kvz_h, i, N, dt, b_h);	
			}
		}
		else{
			update_kernel <<< (NN + 127) / 128, 128 >>> (x_d, y_d, z_d, vx_d, vy_d, vz_d, kx_d, ky_d, kz_d, kvx_d, kvy_d, kvz_d, N, Nperturbers, dt);	
		}
		
		time = time0 + t * dt / dayUnit;


//add condition for backward integration
		if(t % outInterval == 0 && time >= outStart){
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

					unsigned long long int id = id_h[i];
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
		
	
		/*
		if( s < 1.0){
			dt *= 0.5;
		}
		if( s > 3.0){
			dt *= 2.0;
		}
		*/
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
	

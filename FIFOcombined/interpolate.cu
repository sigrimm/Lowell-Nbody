#include "Host.h"
//Neville-Aitken interpolation
// p is perturber index
__host__ void Host::interpolate(double time, int p){

	double Px[Ninterpolate][Ninterpolate];
	double Py[Ninterpolate][Ninterpolate];
	double Pz[Ninterpolate][Ninterpolate];
	double tn[Ninterpolate];

	int it = floor((time - timep0) / dtimep);
	it -= (Ninterpolate / 2);

	for(int i = 0; i < Ninterpolate; ++i){
		int ii = Nperturbers * (i + it) + p;
		Px[0][i] = xp_h[ii];
		Py[0][i] = yp_h[ii];
		Pz[0][i] = zp_h[ii];
		tn[i] = timep_h[ii];

//printf("interpolateP %d %d %d %.20g %.20g %.20g\n", p, i, p + Nperturbers * (i + it), time, tn[i], Px[0][i]);
	}

	for(int j = 1; j < Ninterpolate; ++j){
//printf("****\n");
		for(int i = 0; i < Ninterpolate - j; ++i){
			double t1 = time - tn[i + j];
			double t2 = tn[i] - time;
			double t3 = 1.0 / (tn[i] - tn[i + j]);
			Px[j][i] = (t1 * Px[j-1][i] + t2 * Px[j-1][i+1]) * t3;
			Py[j][i] = (t1 * Py[j-1][i] + t2 * Py[j-1][i+1]) * t3;
			Pz[j][i] = (t1 * Pz[j-1][i] + t2 * Pz[j-1][i+1]) * t3;
//printf("%d %d %g %g %g %g %.20g\n", i, i+j, tn[i], tn[i+j], P[j-1][i], P[j-1][i+1], P[j][i]);

		}
	}
	xt_h[p] = Px[Ninterpolate-1][0];
	yt_h[p] = Py[Ninterpolate-1][0];
	zt_h[p] = Pz[Ninterpolate-1][0];
//printf("interpolate %.20g %d %.20g %.20g %.20g\n", time, p, xt_h[p], yt_h[p], zt_h[p]);

}


__host__ void Host::interpolateTable(double time0){

	double Px[Ninterpolate][Ninterpolate];
	double Py[Ninterpolate][Ninterpolate];
	double Pz[Ninterpolate][Ninterpolate];
	double tn[Ninterpolate];
	double time;

	for(int S = 0; S < RKFn; ++S){

		time = time0 + c_h[S] * dti;

		int it = floor((time - timep0) / dtimep);
		it -= (Ninterpolate / 2);

//printf("interpolateA %d %.20g %.20g %.20g %.20g %.20g %d\n", S, timep0, c_h[S], time0, time, dtt, it);

		for(int p = 0; p < Nperturbers; ++p){

			for(int i = 0; i < Ninterpolate; ++i){
				int ii = Nperturbers * (i + it) + p;
				Px[0][i] = xp_h[ii];
				Py[0][i] = yp_h[ii];
				Pz[0][i] = zp_h[ii];
				tn[i] = timep_h[ii];

//if(p == 1 && S == 0)  printf("interpolateP %d %d %d %.20g %.20g %.20g %.20g\n", p, S, i, timep0, time, tn[i], Px[0][i]);
			}

			for(int j = 1; j < Ninterpolate; ++j){
//if(p == 1 && S == 0)  printf("****\n");
				for(int i = 0; i < Ninterpolate; ++i){
					double t1 = time - tn[i + j];
					double t2 = tn[i] - time;
					double t3 = 1.0 / (tn[i] - tn[i + j]);
					Px[j][i] = (t1 * Px[j-1][i] + t2 * Px[j-1][i + 1]) * t3;
					Py[j][i] = (t1 * Py[j-1][i] + t2 * Py[j-1][i + 1]) * t3;
					Pz[j][i] = (t1 * Pz[j-1][i] + t2 * Pz[j-1][i + 1]) * t3;

//if(p == 1 && S == 0) printf("%d %d %.20g %.20g %.20g %g %g %.20g\n", i, i+j, time, tn[i], tn[i + j], Px[j-1][i], Px[j-1][i + 1], Px[j][i]);
				}
			}

			int ii = p * RKFn + S;
			xTable_h[ii] = Px[Ninterpolate-1][0];
			yTable_h[ii] = Py[Ninterpolate-1][0];
			zTable_h[ii] = Pz[Ninterpolate-1][0];
//printf("interpolateB %.20g %d %d %.20g %.20g %.20g\n", time, S, p, xTable_h[ii], yTable_h[ii], zTable_h[ii]);
		}
	}
}

__host__ void Host::interpolate2(double time, int p){


	//p is the perturber index

	double Cx[Ninterpolate];
	double Cy[Ninterpolate];
	double Cz[Ninterpolate];

	double Dx[Ninterpolate];
	double Dy[Ninterpolate];
	double Dz[Ninterpolate];

	double tn[Ninterpolate];

	int it = floor((time - timep0) / dtimep);
	it -= (Ninterpolate / 2);

	for(int i = 0; i < Ninterpolate; ++i){
		Cx[i] = xp_h[Nperturbers * (i + it) + p];
		Cy[i] = yp_h[Nperturbers * (i + it) + p];
		Cz[i] = zp_h[Nperturbers * (i + it) + p];

		Dx[i] = Cx[i];		
		Dy[i] = Cy[i];		
		Dz[i] = Cz[i];		

		tn[i] = timep_h[Nperturbers * (i + it) + p];

//printf("interpolateC %d %d %.20g %.20g %.20g\n", p, i, time, tn[i], Cx[i]);
	}

	//initialize with closest solution
	//Assume that the closest solution is in the middle

	int ii = Ninterpolate / 2;
	xt_h[p] = Cx[ii];
	yt_h[p] = Cy[ii];
	zt_h[p] = Cz[ii];
//printf("%d %d %g %g %g\n", p, ii, xt_h[p], yt_h[p], zt_h[p]);
	--ii;

	for(int j = 1; j < Ninterpolate; ++j){
//printf("**** %d %d %g\n", j, ii, xt_h[p]);
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

	
//printf("%d %d %.15g %.15g %g %g %g %g\n", i, i+j, tn[i], tn[i+j], dPx, dtn0, dtn1, dtn);

		}

		if(2 * ii < Ninterpolate - j){
			xt_h[p] += Cx[ii + 1];
			yt_h[p] += Cy[ii + 1];
			zt_h[p] += Cz[ii + 1];
		}
		else{
			xt_h[p] += Dx[ii];
			yt_h[p] += Dy[ii];
			zt_h[p] += Dz[ii];
			--ii;
		}

	}
//printf("interpolate %.20g %d %.20g %.20g %.20g\n", time, p, xt[p], yt[p], zt[p]);

}


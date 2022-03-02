
//Neville-Aitken interpolation
// p is perturber index
__host__ void interpolate(int Ninterpolate, int Nperturbers, double *xp, double *yp, double *zp, double *timep, double timep0, double dtimep, double time, double *xt, double *yt, double *zt, int p){

	double Px[Ninterpolate][Ninterpolate];
	double Py[Ninterpolate][Ninterpolate];
	double Pz[Ninterpolate][Ninterpolate];
	double tn[Ninterpolate];

	int it = floor((time - timep0) / dtimep);
	it -= (Ninterpolate / 2);

	for(int i = 0; i < Ninterpolate; ++i){
		Px[0][i] = xp[Nperturbers * (i + it) + p];
		Py[0][i] = yp[Nperturbers * (i + it) + p];
		Pz[0][i] = zp[Nperturbers * (i + it) + p];
		tn[i] = timep[Nperturbers * (i + it) + p];

//printf("interpolateP %d %d %d %.20g %.20g %.20g\n", p, i, p + Nperturbers * (i + it), time, tn[i], Px[0][i]);
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

__host__ void interpolate2(int Ninterpolate, int Nperturbers, double *xp, double *yp, double *zp, double *timep, double timep0, double dtimep, double time, double *xt, double *yt, double *zt, int p){


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
		Cx[i] = xp[Nperturbers * (i + it) + p];
		Cy[i] = yp[Nperturbers * (i + it) + p];
		Cz[i] = zp[Nperturbers * (i + it) + p];

		Dx[i] = Cx[i];		
		Dy[i] = Cy[i];		
		Dz[i] = Cz[i];		

		tn[i] = timep[Nperturbers * (i + it) + p];

//printf("interpolateC %d %d %.20g %.20g %.20g\n", p, i, time, tn[i], Cx[i]);
	}

	//initialize with closest solution
	//Assume that the closest solution is in the middle

	int ii = Ninterpolate / 2;
	xt[p] = Cx[ii];
	yt[p] = Cy[ii];
	zt[p] = Cz[ii];
//printf("%d %d %g %g %g\n", p, ii, xt[p], yt[p], zt[p]);
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

	
//printf("%d %d %.15g %.15g %g %g %g %g\n", i, i+j, tn[i], tn[i+j], dPx, dtn0, dtn1, dtn);

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
		it -= (Ninterpolate / 2);

//printf("interpolateA %d %.20g %.20g %d\n", pid, timep0, time, it);

		if(idx < Ninterpolate){
			Px_s[0][idx] = xp_d[Nperturbers * (idx + it) + pid];
			Py_s[0][idx] = yp_d[Nperturbers * (idx + it) + pid];
			Pz_s[0][idx] = zp_d[Nperturbers * (idx + it) + pid];
			tn_s[idx] = timep_d[Nperturbers * (idx + it) + pid];

///*if(pid == 1)*/ printf("interpolateP %d %d %.20g %.20g %.20g %.20g\n", pid, idx, timep0, time, tn_s[idx], Px_s[0][idx]);
		}
		__syncthreads();

		for(int j = 1; j < Ninterpolate; ++j){
////if(pid == 1) printf("****\n");
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
		it -= (Ninterpolate / 2);

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

		int ii = Ninterpolate / 2;
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
		it -= (Ninterpolate / 2);

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

		int ii = Ninterpolate / 2;
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

template <int Ninterpolate, int Ntable, int RKFn>
__global__ void interpolateTable_kernel(int Nperturbers, double *xp_d, double *yp_d, double *zp_d, double *timep_d, double timep0, double dtimep, double time, double dt, double *xtable_d, double *ytable_d, double *ztable_d){

	int pid = blockIdx.x;	//perturber index, Nperturbers
	int sid = blockIdx.y;	//parallel time step index
	int S = blockIdx.z;	//Stage index
	int idx = threadIdx.x;	//Ninterpolate	, must be a tread index

	if(pid < Nperturbers && sid < Ntable){

		__shared__ double Px_s[Ninterpolate][Ninterpolate];
		__shared__ double Py_s[Ninterpolate][Ninterpolate];
		__shared__ double Pz_s[Ninterpolate][Ninterpolate];
		__shared__ double tn_s[Ninterpolate];

		time += (sid + c_c[S]) * 0.1;

		int it = floor((time - timep0) / dtimep);
		it -= (Ninterpolate / 2);

//if(idx == 0 && pid == 0) printf("interpolateA %d %d %d %.20g %.20g %.20g %d\n", pid, sid, S, timep0, c_c[S], time, it);

		if(idx < Ninterpolate){
			Px_s[0][idx] = xp_d[Nperturbers * (idx + it) + pid];
			Py_s[0][idx] = yp_d[Nperturbers * (idx + it) + pid];
			Pz_s[0][idx] = zp_d[Nperturbers * (idx + it) + pid];
			tn_s[idx] = timep_d[Nperturbers * (idx + it) + pid];

//if(pid == 1 && sid < 2 && S == 0)  printf("interpolateP %d %d %d %d %.20g %.20g %.20g %.20g\n", pid, sid, S, idx, timep0, time, tn_s[idx], Px_s[0][idx]);
		}
		__syncthreads();

		for(int j = 1; j < Ninterpolate; ++j){
//if(pid == 1 && sid < 2 && S == 0)  printf("****\n");
			if(idx < Ninterpolate - j){
				Px_s[j][idx] = ((time - tn_s[idx + j]) * Px_s[j-1][idx] + (tn_s[idx] - time) * Px_s[j-1][idx + 1]) / (tn_s[idx] - tn_s[idx + j]);
				Py_s[j][idx] = ((time - tn_s[idx + j]) * Py_s[j-1][idx] + (tn_s[idx] - time) * Py_s[j-1][idx + 1]) / (tn_s[idx] - tn_s[idx + j]);
				Pz_s[j][idx] = ((time - tn_s[idx + j]) * Pz_s[j-1][idx] + (tn_s[idx] - time) * Pz_s[j-1][idx + 1]) / (tn_s[idx] - tn_s[idx + j]);

//if(pid == 1 && sid < 2 && S == 0) printf("%d %d %d %.20g %.20g %.20g %g %g %.20g\n", sid, idx, idx+j, time, tn_s[idx], tn_s[idx + j], Px_s[j-1][idx], Px_s[j-1][idx + 1], Px_s[j][idx]);
			}
			__syncthreads();
		}

		if(idx == 0){
			xtable_d[pid * Ntable + sid] = Px_s[Ninterpolate-1][0];
			ytable_d[pid * Ntable + sid] = Py_s[Ninterpolate-1][0];
			ztable_d[pid * Ntable + sid] = Pz_s[Ninterpolate-1][0];
// /*if(pid == 1) */printf("interpolateB %.20g %d %d %d %.20g %.20g %.20g\n", time, sid, S, pid, xtable_d[pid * Ntable + sid], ytable_d[pid * Ntable + sid], ztable_d[pid * Ntable + sid]);
		}

	}
}


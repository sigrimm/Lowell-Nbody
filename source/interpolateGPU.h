
__global__ void interpolate_kernel(double *xp_d, double *yp_d, double *zp_d, double *timep_d, double timep0, double dtimep, double time, double *xt_d, double *yt_d, double *zt_d){
	
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
			int ii = Nperturbers * (idx + it) + pid;
			Px_s[0][idx] = xp_d[ii];
			Py_s[0][idx] = yp_d[ii];
			Pz_s[0][idx] = zp_d[ii];
			tn_s[idx] = timep_d[ii];
			
			///*if(pid == 1)*/ printf("interpolateP %d %d %.20g %.20g %.20g %.20g\n", pid, idx, timep0, time, tn_s[idx], Px_s[0][idx]);
		}
		__syncthreads();
		
		for(int j = 1; j < Ninterpolate; ++j){
			////if(pid == 1) printf("****\n");
			if(idx < Ninterpolate - j){
				double t1 = time - tn_s[idx + j];
				double t2 = tn_s[idx] - time;
				double t3 = 1.0 / (tn_s[idx] - tn_s[idx + j]);
				Px_s[j][idx] = (t1 * Px_s[j-1][idx] + t2 * Px_s[j-1][idx + 1]) * t3;
				Py_s[j][idx] = (t1 * Py_s[j-1][idx] + t2 * Py_s[j-1][idx + 1]) * t3;
				Pz_s[j][idx] = (t1 * Pz_s[j-1][idx] + t2 * Pz_s[j-1][idx + 1]) * t3;
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


//using local memory arrays
__global__ void interpolate2b_kernel(double *xp_d, double *yp_d, double *zp_d, double *timep_d, double timep0, double dtimep, double time, double *xt_d, double *yt_d, double *zt_d){
	
	//every perturber is a thread
	
	int pid = threadIdx.x;	//perturber index, Nperturbers
	
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
		
		//if(pid == 1) printf("interpolateB %d %.20g %.20g %.20g\n", pid, timep0, time, tn_s[0]);
		
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

__global__ void interpolateTable_kernel(double *xp_d, double *yp_d, double *zp_d, double *timep_d, double timep0, double dtimep, double time, double dtt, double *xTable_d, double *yTable_d, double *zTable_d){
	
	int pid = blockIdx.x;	//perturber index, Nperturbers
	int S = blockIdx.y;	//Stage index
	int idx = threadIdx.x;	//Ninterpolate	, must be a tread index
	
	if(pid < Nperturbers){
		
		__shared__ double Px_s[Ninterpolate][Ninterpolate];
		__shared__ double Py_s[Ninterpolate][Ninterpolate];
		__shared__ double Pz_s[Ninterpolate][Ninterpolate];
		__shared__ double tn_s[Ninterpolate];
		
		time += c_c[S] * dtt;
		
		int it = floor((time - timep0) / dtimep);
		it -= (Ninterpolate / 2);
		
		//if(idx == 0 && pid == 0) printf("interpolateA %d %d %.20g %.20g %.20g %d\n", pid, S, timep0, c_c[S], time, it);
		
		if(idx < Ninterpolate){
			int ii = Nperturbers * (idx + it) + pid;
			Px_s[0][idx] = xp_d[ii];
			Py_s[0][idx] = yp_d[ii];
			Pz_s[0][idx] = zp_d[ii];
			tn_s[idx] = timep_d[ii];
			
			//if(pid == 1 && S == 0)  printf("interpolateP %d %d %d %.20g %.20g %.20g %.20g\n", pid, S, idx, timep0, time, tn_s[idx], Px_s[0][idx]);
		}
		__syncthreads();
		
		for(int j = 1; j < Ninterpolate; ++j){
			//if(pid == 1 && S == 0)  printf("****\n");
			if(idx < Ninterpolate - j){
				double t1 = time - tn_s[idx + j];
				double t2 = tn_s[idx] - time;
				double t3 = 1.0 / (tn_s[idx] - tn_s[idx + j]);
				Px_s[j][idx] = (t1 * Px_s[j-1][idx] + t2 * Px_s[j-1][idx + 1]) * t3;
				Py_s[j][idx] = (t1 * Py_s[j-1][idx] + t2 * Py_s[j-1][idx + 1]) * t3;
				Pz_s[j][idx] = (t1 * Pz_s[j-1][idx] + t2 * Pz_s[j-1][idx + 1]) * t3;
				
				//if(pid == 1 && S == 0) printf("%d %d %.20g %.20g %.20g %g %g %.20g\n", idx, idx+j, time, tn_s[idx], tn_s[idx + j], Px_s[j-1][idx], Px_s[j-1][idx + 1], Px_s[j][idx]);
			}
			__syncthreads();
		}
		
		if(idx == 0){
			int ii = pid * RKFn + S;
			xTable_d[ii] = Px_s[Ninterpolate-1][0];
			yTable_d[ii] = Py_s[Ninterpolate-1][0];
			zTable_d[ii] = Pz_s[Ninterpolate-1][0];
			//printf("interpolateB %.20g %d %d %.20g %.20g %.20g\n", time, S, pid, xTable_d[ii], yTable_d[ii], zTable_d[ii]);
		}
		
	}
}

//using local memory arrays
__global__ void interpolate2bTable_kernel(double *xp_d, double *yp_d, double *zp_d, double *timep_d, double timep0, double dtimep, double time, double dts, double *xTable_d, double *yTable_d, double *zTable_d){
	
	//every perturber is a thread
	//every particle is a block
	
	int pid = blockIdx.y;	//perturber index, Nperturbers
	int S = blockIdx.z;
	
	if(pid < Nperturbers){
		
		double x, y, z;
		
		double Cx[Ninterpolate];
		double Cy[Ninterpolate];
		double Cz[Ninterpolate];
		
		double Dx[Ninterpolate];
		double Dy[Ninterpolate];
		double Dz[Ninterpolate];
		
		double tn[Ninterpolate];
		
		time += c_c[S] * dts;
		
		int it = floor((time - timep0) / dtimep);
		it -= (Ninterpolate / 2);
		
		//printf("interpolateA %d %d %g %.20g %.20g %d\n", pid, S, c_c[S], timep0, time, it);
		
		
		for(int i = 0; i < Ninterpolate; ++i){
			Cx[i] = xp_d[Nperturbers * (i + it) + pid];
			Cy[i] = yp_d[Nperturbers * (i + it) + pid];
			Cz[i] = zp_d[Nperturbers * (i + it) + pid];
			
			Dx[i] = Cx[i];
			Dy[i] = Cy[i];
			Dz[i] = Cz[i];
			
			tn[i] = timep_d[Nperturbers * (i + it) + pid];
			
			//if(pid == 1) printf("interpolateB %d %.20g %.20g %.20g %.20g\n", pid, timep0, time, tn_s[i], Cx[i]);
			
		}
		
		
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
		
		xTable_d[pid * RKFn + S] = x;
		yTable_d[pid * RKFn + S] = y;
		zTable_d[pid * RKFn + S] = z;
		//if(pid == 1) printf("interpolateB %.20g %d %.20g %.20g %.20g\n", time, pid, xTable_d[pid *  * RKFn + S], yTable_d[pid * RKFn + S], zTable_d[pid * RKFn + S]);
		
	}
}

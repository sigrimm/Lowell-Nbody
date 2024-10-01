__device__ void  update_ChebyshevGPU(double *data_d, double *cdata_d, int *id_d, double &startTime, double &endTime, const int nChebyshev, int &offset0, double time, double time_reference){

	//Update the Chebyshev coefficients if necessary
	//int p = threadIdx.x;	//Perturber index

	if(time + time_reference > endTime){
		for(int k = 0; k < 1000000; ++k){

			startTime = data_d[offset0];
			endTime = data_d[offset0 + 1];

//printf(" ++ %d %d %d %d %.20g %.20g\n", p, id_d[p], offset0, nChebyshev, startTime, endTime);

			for(int i = 0; i < nChebyshev * 3; ++i){
				cdata_d[i] = data_d[offset0 + 2 + i]; 
			}
			

			if(time + time_reference <= endTime){
				break;
			}

			offset0 += nChebyshev * 3 + 2;
		}
	}
	if(time + time_reference < startTime){
		for(int k = 0; k < 1000000; ++k){

			startTime = data_d[offset0];
			endTime = data_d[offset0 + 1];

//printf(" -- %d %d %d %d %.20g %.20g\n", p, id_d[p], offset0, nChebyshev, startTime, endTime);

			for(int i = 0; i < nChebyshev * 3; ++i){
				cdata_d[i] = data_d[offset0 + 2 + i]; 
			}

			if(time + time_reference >= startTime){
				break;
			}

			offset0 -= nChebyshev * 3 + 2;
		}
	}
}

__device__ void update_perturbersGPU(double *xTable_s, double *yTable_s, double *zTable_s, double *vxTable_s, double *vyTable_s, double *vzTable_s, double *cdata_d, double &startTime, double &endTime, const int nChebyshev, double time_reference, double time, const int nCm, const double EM, const double AUtokm, const int Nperturbers, int p){

	if(p < Nperturbers){

		//double Tx[nCm];
		//double Ty[nCm];
		//double Tz[nCm];
		//double Tvx[nCm];
		//double Tvy[nCm];
		//double Tvz[nCm];

		double Tx[14];
		double Ty[14];
		double Tz[14];
		double Tvx[14];
		double Tvy[14];
		double Tvz[14];


		//Calculate positions of perturbers
		double sizeSubInterval = endTime - startTime;
//printf("time %d %d %.20g %.20g %.20g %.20g\n", blockIdx.x, p, time + time_reference, sizeSubInterval, startTime, endTime);
		double subTime = (time_reference - startTime + time) / sizeSubInterval;   //normalized time in  0 - 1
		double t = 2.0 * subTime - 1.0;                         //mormalized time in -1 - 1

		//double ct = 2.0 / sizeSubInterval;                    //correction factor for time units

		//This is done so in Assist, remove the time factor later
		double ct = 2.0 / sizeSubInterval / 86400.0;            //correction factor for time units

//if(p > 10) printf("Chebyshev time %d %.20g %.20g %.20g\n", p, subTime, time, t);


		Tx[0] = 1.0;
		Tx[1] = t;
		Ty[0] = 1.0;
		Ty[1] = t;
		Tz[0] = 1.0;
		Tz[1] = t;

		Tvx[0] = 0.0;
		Tvx[1] = 1.0;
		Tvy[0] = 0.0;
		Tvy[1] = 1.0;
		Tvz[0] = 0.0;
		Tvz[1] = 1.0;

		double t2 = 2.0 * t;

		for(int j = 2; j < nChebyshev; ++j){
			Tx[j] = t2 * Tx[j - 1] - Tx[j - 2];
			Ty[j] = t2 * Ty[j - 1] - Ty[j - 2];
			Tz[j] = t2 * Tz[j - 1] - Tz[j - 2];
			Tvx[j] = t2 * Tvx[j - 1] + 2.0 * Tx[j - 1] - Tvx[j - 2];
			Tvy[j] = t2 * Tvy[j - 1] + 2.0 * Ty[j - 1] - Tvy[j - 2];
			Tvz[j] = t2 * Tvz[j - 1] + 2.0 * Tz[j - 1] - Tvz[j - 2];
		}

		double xp = 0.0;
		double yp = 0.0;
		double zp = 0.0;

		double vxp = 0.0;
		double vyp = 0.0;
		double vzp = 0.0;

		for(int j = 0; j < nChebyshev; ++j){
		//for(int j = nChebyshev - 1; j >= 0; --j){    //reduce floating point errors by revert order
			double cx = cdata_d[j];
			double cy = cdata_d[nChebyshev + j];
			double cz = cdata_d[2 * nChebyshev + j];

			xp += Tx[j] * cx;
			yp += Ty[j] * cy;
			zp += Tz[j] * cz;
//printf("Chebyshev %d %.20g %.20g %.20g %.20g %.20g\n", j, Tx[j], cx, xp, yp, zp);                    

			vxp += Tvx[j] * cx * ct;
			vyp += Tvy[j] * cy * ct;
			vzp += Tvz[j] * cz * ct;
		}

		xTable_s[p] = xp;
		yTable_s[p] = yp;
		zTable_s[p] = zp;

		vxTable_s[p] = vxp;
		vyTable_s[p] = vyp;
		vzTable_s[p] = vzp;

//printf("positionA %d %.20g %.20g %.20g %.20g %.20g\n", p, time, xTable_s[p], yTable_s[p], zTable_s[p], t);
//printf("positionvA %d %.20g %.20g %.20g %.20g %.20g\n", p, time, vxTable_s[p], vyTable_s[p], vzTable_s[p], t);


//printf("positionB %d %.20g %.20g %.20g %.20g %.20g\n", p, time, xTable_s[p] / AUtokm, yTable_s[p] / AUtokm, zTable_s[p] / AUtokm, t);

//printf("%d %.20g %.20g %.20g %.20g %.20g\n", p, time, xTable_s[p], yTable_s[p], zTable_s[p], t);


		//Calculate Earth and Moon positions, id 2 and 9
		//Up to here id = 2 is the Earth-Moon barycentrum, id = 9 is the geocentric position of the Moon        
	}

	__syncthreads();

	if(p == 0){

		double xB = xTable_s[2];
		double yB = yTable_s[2];
		double zB = zTable_s[2];

		double vxB = vxTable_s[2];
		double vyB = vyTable_s[2];
		double vzB = vzTable_s[2];

		double xM = xTable_s[9];
		double yM = yTable_s[9];
		double zM = zTable_s[9];

		double vxM = vxTable_s[9];
		double vyM = vyTable_s[9];
		double vzM = vzTable_s[9];

		double f = 1.0/(1.0 + EM);

		xTable_s[2] = xB - xM * f;
		yTable_s[2] = yB - yM * f;
		zTable_s[2] = zB - zM * f;

		vxTable_s[2] = vxB - vxM * f;
		vyTable_s[2] = vyB - vyM * f;
		vzTable_s[2] = vzB - vzM * f;

		xTable_s[9] = xB + xM * EM * f;
		yTable_s[9] = yB + yM * EM * f;
		zTable_s[9] = zB + zM * EM * f;

		vxTable_s[9] = vxB + vxM * EM * f;
		vyTable_s[9] = vyB + vyM * EM * f;
		vzTable_s[9] = vzB + vzM * EM * f;

	}

	__syncthreads();

	if(p < Nperturbers){

		//positions are in km
		//velocities are in km/day
		xTable_s[p] /= AUtokm;
		yTable_s[p] /= AUtokm;
		zTable_s[p] /= AUtokm;

		//vxTable_s[p] /= AUtokm;
		//vyTable_s[p] /= AUtokm;
		//vzTable_s[p] /= AUtokm;
		//remove time factor again
		vxTable_s[p] /= AUtokm / 86400.0;
		vyTable_s[p] /= AUtokm / 86400.0;
		vzTable_s[p] /= AUtokm / 86400.0;

//printf("positionB %d %d %.20g %.20g %.20g %.20g\n", blockIdx.x, p, time, xTable_s[p], yTable_s[p], zTable_s[p]);


		/*
		//print in the order of asssit
		int pp = p + 10;
		if(p == 0) pp = 11;     
		if(p == 1) pp = 10;     //Sun
		if(p == 17) pp = 8;     //Pluto
		if(p == 18) pp = 9;     //Moon
		if(p == 19) pp = 3;     //Mars
		if(p == 20) pp = 0;     //Mercurs
		if(p == 21) pp = 7;     //Neptune
		if(p == 22) pp = 6;     //Uranus
		if(p == 23) pp = 2;     //Earth
		if(p == 24) pp = 1;     //Venus
		if(p == 25) pp = 5;     //Saturn
		if(p == 26) pp = 4;     //Jupiter


		//printf("positionB %d %.20g %.20g %.20g %.20g\n", pp, time, xTable_s[pp], yTable_s[pp], zTable_s[pp]);
		//printf("positionvB %d %.20g %.20g %.20g %.20g\n", pp, time, vxTable_s[pp], vyTable_s[pp], vzTable_s[pp]);
		//printf("positionB %d %.20g %.20g %.20g %.20g\n", p, time, xTable_s[p], yTable_s[p], zTable_s[p]);
		*/

	}
	__syncthreads();
	//Translate asteroid orbits from Heliocentric to Barycentric coordinates
	//This is done so in Assist, but probably this is not correct, check!
	if(p >= 11 && p < Nperturbers){
		xTable_s[p] += xTable_s[10];
		yTable_s[p] += yTable_s[10];
		zTable_s[p] += zTable_s[10];

		vxTable_s[p] += vxTable_s[10];
		vyTable_s[p] += vyTable_s[10];
		vzTable_s[p] += vzTable_s[10];
	}
	__syncthreads();
}



// The calculation of the perturbers position can be parallelized along the perturbers and the stage index of the Runke Kutta integrator.
// Every perturber runs on a thread.
// Every stage runs on a different thread block.
// This way of scheduling is needed because perturbers must be in the same thread block because of the calcualtion of the Earth and Moon.
__global__ void update_perturbers_kernel(double *xTable_d, double *yTable_d, double *zTable_d, double *vxTable_d, double *vyTable_d, double *vzTable_d, double *data_d, double *cdata_d, int *id_d, double *startTime_d, double *endTime_d, int *nChebyshev_d, int *offset0_d, double time, double time_reference, const double dt, const int RKFn, const int nCm, const double EM, const double AUtokm, const int Nperturbers){

	int id = threadIdx.x;	//Perturber index
	int idx = blockIdx.x;	//Stage index

	__shared__ double xTable_s[def_NP];
	__shared__ double yTable_s[def_NP];
	__shared__ double zTable_s[def_NP];

	__shared__ double vxTable_s[def_NP];
	__shared__ double vyTable_s[def_NP];
	__shared__ double vzTable_s[def_NP];

	int offset0;
	double startTime;
	double endTime;
	int nChebyshev;

	const int pp = (id * RKFn + idx) * nCm * 3;

	if(id < Nperturbers && idx < RKFn){

		offset0 = offset0_d[id * RKFn + idx];
		startTime = startTime_d[id * RKFn + idx];
		endTime = endTime_d[id * RKFn + idx];
		nChebyshev = nChebyshev_d[id];


		update_ChebyshevGPU(data_d, cdata_d + pp, id_d, startTime, endTime, nChebyshev, offset0, time + RKFc_c[idx] * dt, time_reference);
	}
	__syncthreads();
	update_perturbersGPU(xTable_s, yTable_s, zTable_s, vxTable_s, vyTable_s, vzTable_s, cdata_d + pp, startTime, endTime, nChebyshev, time_reference, time + RKFc_c[idx] * dt, nCm, EM, AUtokm, Nperturbers, id);

	__syncthreads();
	if(id < Nperturbers && idx < RKFn){

		offset0_d[id * RKFn + idx] = offset0;
		startTime_d[id * RKFn + idx] = startTime;
		endTime_d[id * RKFn + idx] = endTime;

		xTable_d[id * RKFn + idx] = xTable_s[id];
		yTable_d[id * RKFn + idx] = yTable_s[id];
		zTable_d[id * RKFn + idx] = zTable_s[id];

		vxTable_d[id * RKFn + idx] = vxTable_s[id];
		vyTable_d[id * RKFn + idx] = vyTable_s[id];
		vzTable_d[id * RKFn + idx] = vzTable_s[id];
	}
}

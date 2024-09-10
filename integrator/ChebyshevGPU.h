__device__ void  update_ChebyshevGPU(double *data_d, double *cdata_d, int *id_d, double &startTime, double &endTime, const int nChebyshev, int &offset0, double time, double time_reference, const int nCm, int p){

	//Update the Chebyshev coefficients if necessary

	const int pp = p * nCm * 3;

	if(time + time_reference > endTime){
		for(int k = 0; k < 1000000; ++k){

			startTime = data_d[offset0];
			endTime = data_d[offset0 + 1];


//printf(" ++ %d %d %d %d %.20g %.20g\n", p, id_d[p], offset0, nChebyshev, startTime, endTime);


			for(int i = 0; i < nChebyshev * 3; ++i){
				cdata_d[pp + i] = data_d[offset0 + 2 + i]; 
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
				cdata_d[pp + i] = data_d[offset0 + 2 + i]; 
			}

			if(time + time_reference >= startTime){
				break;
			}

			offset0 -= nChebyshev * 3 + 2;
		}
	}
}

__device__ void update_perturbersGPU(double *x_s, double *y_s, double *z_s, double *vx_s, double *vy_s, double *vz_s, double *cdata_d, double &startTime, double &endTime, const int nChebyshev, double time_reference, double time, const int nCm, const double EM, const double AUtokm, const int Nperturbers, int p){

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
		const int pp = p * nCm * 3;

		double sizeSubInterval = endTime - startTime;
//printf("time %d %.20g %.20g %.20g %.20g\n", p, time + time_reference, sizeSubInterval, startTime, endTime);
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

		for(int j = 2; j < nChebyshev; ++j){
			Tx[j] = 2.0 * t * Tx[j - 1] - Tx[j - 2];
			Ty[j] = 2.0 * t * Ty[j - 1] - Ty[j - 2];
			Tz[j] = 2.0 * t * Tz[j - 1] - Tz[j - 2];
			Tvx[j] = 2.0 * t * Tvx[j - 1] + 2.0 * Tx[j - 1] - Tvx[j - 2];
			Tvy[j] = 2.0 * t * Tvy[j - 1] + 2.0 * Ty[j - 1] - Tvy[j - 2];
			Tvz[j] = 2.0 * t * Tvz[j - 1] + 2.0 * Tz[j - 1] - Tvz[j - 2];
		}

		double xp = 0.0;
		double yp = 0.0;
		double zp = 0.0;

		double vxp = 0.0;
		double vyp = 0.0;
		double vzp = 0.0;

		for(int j = 0; j < nChebyshev; ++j){
		//for(int j = nChebyshev - 1; j >= 0; --j){    //reduce floating point errors by revert order
			xp += Tx[j] * cdata_d[pp + j];
			yp += Ty[j] * cdata_d[pp + nChebyshev + j];
			zp += Tz[j] * cdata_d[pp + 2 * nChebyshev + j];
//printf("Chebyshev %d %.20g %.20g %.20g %.20g %.20g\n", j, Tx[j], cdata_d[pp + j], xp, yp, zp);                    

			vxp += Tvx[j] * cdata_d[pp + j] * ct;
			vyp += Tvy[j] * cdata_d[pp + nChebyshev + j] * ct;
			vzp += Tvz[j] * cdata_d[pp + 2 * nChebyshev + j] * ct;
		}

		x_s[p] = xp;
		y_s[p] = yp;
		z_s[p] = zp;

		vx_s[p] = vxp;
		vy_s[p] = vyp;
		vz_s[p] = vzp;

//printf("positionA %d %.20g %.20g %.20g %.20g %.20g\n", p, time, x_s[p], y_s[p], z_s[p], t);
//printf("positionvA %d %.20g %.20g %.20g %.20g %.20g\n", p, time, vx_s[p], vy_s[p], vz_s[p], t);


//printf("positionB %d %.20g %.20g %.20g %.20g %.20g\n", p, time, x_s[p] / AUtokm, y_s[p] / AUtokm, z_s[p] / AUtokm, t);

//printf("%d %.20g %.20g %.20g %.20g %.20g\n", p, time, x_s[p], y_s[p], z_s[p], t);


		//Calculate Earth and Moon positions, id 2 and 9
		//Up to here id = 2 is the Earth-Moon barycentrum, id = 9 is the geocentric position of the Moon        
	}

	__syncthreads();

	if(p == 0){

		double xB = x_s[2];
		double yB = y_s[2];
		double zB = z_s[2];

		double vxB = vx_s[2];
		double vyB = vy_s[2];
		double vzB = vz_s[2];

		double xM = x_s[9];
		double yM = y_s[9];
		double zM = z_s[9];

		double vxM = vx_s[9];
		double vyM = vy_s[9];
		double vzM = vz_s[9];

		double f = 1.0/(1.0 + EM);

		x_s[2] = xB - xM * f;
		y_s[2] = yB - yM * f;
		z_s[2] = zB - zM * f;

		vx_s[2] = vxB - vxM * f;
		vy_s[2] = vyB - vyM * f;
		vz_s[2] = vzB - vzM * f;

		x_s[9] = xB + xM * EM * f;
		y_s[9] = yB + yM * EM * f;
		z_s[9] = zB + zM * EM * f;

		vx_s[9] = vxB + vxM * EM * f;
		vy_s[9] = vyB + vyM * EM * f;
		vz_s[9] = vzB + vzM * EM * f;

	}

	__syncthreads();

	if(p < Nperturbers){

		//positions are in km
		//velocities are in km/day
		x_s[p] /= AUtokm;
		y_s[p] /= AUtokm;
		z_s[p] /= AUtokm;

		//vx_s[p] /= AUtokm;
		//vy_s[p] /= AUtokm;
		//vz_s[p] /= AUtokm;
		//remove time factor again
		vx_s[p] /= AUtokm / 86400.0;
		vy_s[p] /= AUtokm / 86400.0;
		vz_s[p] /= AUtokm / 86400.0;

//printf("positionB %d %.20g %.20g %.20g %.20g\n", p, time, x_s[p], y_s[p], z_s[p]);


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


		//printf("positionB %d %.20g %.20g %.20g %.20g\n", pp, time, x_s[pp], y_s[pp], z_s[pp]);
		//printf("positionvB %d %.20g %.20g %.20g %.20g\n", pp, time, vx_s[pp], vy_s[pp], vz_s[pp]);
		//printf("positionB %d %.20g %.20g %.20g %.20g\n", p, time, x_s[p], y_s[p], z_s[p]);
		*/

	//Translate asteroid orbits from Heliocentric to Barycentric coordinates
	//This is done so in Assist, but probably this is not correct, check!
	}
	__syncthreads();
	if(p >= 11 && p < Nperturbers){
		x_s[p] += x_s[10];
		y_s[p] += y_s[10];
		z_s[p] += z_s[10];
	}

}




__global__ void update_perturbers_kernel(double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *data_d, double *cdata_d, int *id_d, double *startTime_d, double *endTime_d, int *nChebyshev_d, int *offset0_d, double time, double time_reference, const int nCm, const double EM, const double AUtokm, const int Nperturbers){

	int id = blockIdx.x * blockDim.x + threadIdx.x;

	__shared__ double x_s[def_NP];
	__shared__ double y_s[def_NP];
	__shared__ double z_s[def_NP];

	__shared__ double vx_s[def_NP];
	__shared__ double vy_s[def_NP];
	__shared__ double vz_s[def_NP];

	int offset0;
	double startTime;
	double endTime;
	int nChebyshev;

	if(id < Nperturbers){

		x_s[id] = x_d[id];
		y_s[id] = y_d[id];
		z_s[id] = z_d[id];

		vx_s[id] = vx_d[id];
		vy_s[id] = vy_d[id];
		vz_s[id] = vz_d[id];

		offset0 = offset0_d[id];
		startTime = startTime_d[id];
		endTime = endTime_d[id];
		nChebyshev = nChebyshev_d[id];


		update_ChebyshevGPU(data_d, cdata_d, id_d, startTime, endTime, nChebyshev, offset0, time, time_reference, nCm, id);
	}
	__syncthreads();
	update_perturbersGPU(x_s, y_s, z_s, vx_s, vy_s, vz_s, cdata_d, startTime, endTime, nChebyshev, time_reference, time, nCm, EM, AUtokm, Nperturbers, id);

	__syncthreads();
	if(id < Nperturbers){

		offset0_d[id] = offset0;
		startTime_d[id] = startTime;
		endTime_d[id] = endTime;

		x_d[id] = x_s[id];
		y_d[id] = y_s[id];
		z_d[id] = z_s[id];

		vx_d[id] = vx_s[id];
		vy_d[id] = vy_s[id];
		vz_d[id] = vz_s[id];
	}
}

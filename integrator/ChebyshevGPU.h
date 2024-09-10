__device__ void  update_ChebyshevGPU(double *data_d, double *cdata_d, int *id_d, double *startTime_d, double *endTime_d, int *nChebyshev_d, int *offset0_d, double time, double time_reference, const int nCm, int p){

	//Update the Chebyshev coefficients if necessary

	const int pp = p * nCm * 3;
	if(time + time_reference > endTime_d[p]){
		for(int k = 0; k < 1000000; ++k){

			startTime_d[p] = data_d[offset0_d[p]];
			endTime_d[p] = data_d[offset0_d[p] + 1];


printf(" ++ %d %d %d %d %.20g %.20g\n", p, id_d[p], offset0_d[p], nChebyshev_d[p], startTime_d[p], endTime_d[p]);


			for(int i = 0; i < nChebyshev_d[p] * 3; ++i){
				cdata_d[pp + i] = data_d[offset0_d[p] + 2 + i]; 
			}
			

			if(time + time_reference <= endTime_d[p]){
				break;
			}

			offset0_d[p] += nChebyshev_d[p] * 3 + 2;
		}
	}
	if(time + time_reference < startTime_d[p]){
		for(int k = 0; k < 1000000; ++k){

			startTime_d[p] = data_d[offset0_d[p]];
			endTime_d[p] = data_d[offset0_d[p] + 1];


printf(" -- %d %d %d %d %.20g %.20g\n", p, id_d[p], offset0_d[p], nChebyshev_d[p], startTime_d[p], endTime_d[p]);

			for(int i = 0; i < nChebyshev_d[p] * 3; ++i){
				cdata_d[pp + i] = data_d[offset0_d[p] + 2 + i]; 
			}

			if(time + time_reference >= startTime_d[p]){
				break;
			}

			offset0_d[p] -= nChebyshev_d[p] * 3 + 2;
		}
	}
}

__device__ void update_perturbersGPU(double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *cdata_d, double *startTime_d, double *endTime_d, int *nChebyshev_d, double time_reference, double time, const int nCm, const double EM, const double AUtokm, const int Nperturbers, int p){

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
	const int nC = nChebyshev_d[p];

	const int pp = p * nCm * 3;

	double sizeSubInterval = endTime_d[p] - startTime_d[p];
//if(p > 10) printf("time %d %.20g %.20g %.20g %.20g\n", p, time + time_reference, sizeSubInterval, startTime[p], endTime[p]);
	double subTime = (time_reference - startTime_d[p] + time) / sizeSubInterval;   //normalized time in  0 - 1
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

	for(int j = 2; j < nC; ++j){
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

	for(int j = 0; j < nC; ++j){
	//for(int j = nChebyshev_h[p] - 1; j >= 0; --j){    //reduce floating point errors by revert order
		xp += Tx[j] * cdata_d[pp + j];
		yp += Ty[j] * cdata_d[pp + nC + j];
		zp += Tz[j] * cdata_d[pp + 2 * nC + j];
//printf("Chebyshev %d %.20g %.20g %.20g\n", j, Tx[j], cdata_d[pp + j], x[p]);                    

		vxp += Tvx[j] * cdata_d[pp + j] * ct;
		vyp += Tvy[j] * cdata_d[pp + nC + j] * ct;
		vzp += Tvz[j] * cdata_d[pp + 2 * nC + j] * ct;
	}

	x_d[p] = xp;
	y_d[p] = yp;
	z_d[p] = zp;

	vx_d[p] = vxp;
	vy_d[p] = vyp;
	vz_d[p] = vzp;

//printf("positionA %d %.20g %.20g %.20g %.20g %.20g\n", p, time, x_d[p], y_d[p], z_d[p], t);
//printf("positionvA %d %.20g %.20g %.20g %.20g %.20g\n", p, time, vx_d[p], vy_d[p], vz_d[p], t);


//printf("positionB %d %.20g %.20g %.20g %.20g %.20g\n", p, time, x_d[p] / AUtokm, y_d[p] / AUtokm, z_d[p] / AUtokm, t);

//printf("%d %.20g %.20g %.20g %.20g %.20g\n", p, time, x_d[p], y_d[p], z_d[p], t);


	//Calculate Earth and Moon positions, id 2 and 9
	//Up to here id = 2 is the Earth-Moon barycentrum, id = 9 is the geocentric position of the Moon        

	__syncthreads();

	if(p == 0){

		double xB = x_d[2];
		double yB = y_d[2];
		double zB = z_d[2];

		double vxB = vx_d[2];
		double vyB = vy_d[2];
		double vzB = vz_d[2];

		double xM = x_d[9];
		double yM = y_d[9];
		double zM = z_d[9];

		double vxM = vx_d[9];
		double vyM = vy_d[9];
		double vzM = vz_d[9];

		double f = 1.0/(1.0 + EM);

		x_d[2] = xB - xM * f;
		y_d[2] = yB - yM * f;
		z_d[2] = zB - zM * f;

		vx_d[2] = vxB - vxM * f;
		vy_d[2] = vyB - vyM * f;
		vz_d[2] = vzB - vzM * f;

		x_d[9] = xB + xM * EM * f;
		y_d[9] = yB + yM * EM * f;
		z_d[9] = zB + zM * EM * f;

		vx_d[9] = vxB + vxM * EM * f;
		vy_d[9] = vyB + vyM * EM * f;
		vz_d[9] = vzB + vzM * EM * f;

	}

	__syncthreads();

	//positions are in km
	//velocities are in km/day
	x_d[p] /= AUtokm;
	y_d[p] /= AUtokm;
	z_d[p] /= AUtokm;

	//vx_h[p] /= AUtokm;
	//vy_h[p] /= AUtokm;
	//vz_h[p] /= AUtokm;
	//remove time factor again
	vx_d[p] /= AUtokm / 86400.0;
	vy_d[p] /= AUtokm / 86400.0;
	vz_d[p] /= AUtokm / 86400.0;
//printf("positionB %d %.20g %.20g %.20g %.20g\n", p, time, x_h[p], y_h[p], z_h[p]);


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


	//printf("positionB %d %.20g %.20g %.20g %.20g\n", pp, time, x_d[pp], y_d[pp], z_d[pp]);
	//printf("positionvB %d %.20g %.20g %.20g %.20g\n", pp, time, vx_d[pp], vy_d[pp], vz_d[pp]);
	//printf("positionB %d %.20g %.20g %.20g %.20g\n", p, time, x_d[p], y_d[p], z_d[p]);
	*/

	//Translate asteroid orbits from Heliocentric to Barycentric coordinates
	//This is done so in Assist, but probably this is not correct, check!

	__syncthreads();
	if(p >= 11){
		x_d[p] += x_d[10];
		y_d[p] += y_d[10];
		z_d[p] += z_d[10];
	}

}




__global__ void update_perturbers_kernel(double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *data_d, double *cdata_d, int *id_d, double *startTime_d, double *endTime_d, int *nChebyshev_d, int *offset0_d, double time, double time_reference, const int nCm, const double EM, const double AUtokm, const int Nperturbers){

	int id = blockIdx.x * blockDim.x + threadIdx.x;

	if(id < Nperturbers){
		update_ChebyshevGPU(data_d, cdata_d, id_d, startTime_d, endTime_d, nChebyshev_d, offset0_d, time, time_reference, nCm, id);
		update_perturbersGPU(x_d, y_d, z_d, vx_d, vy_d, vz_d, cdata_d, startTime_d, endTime_d, nChebyshev_d, time_reference, time, nCm, EM, AUtokm, Nperturbers, id);

	}
}

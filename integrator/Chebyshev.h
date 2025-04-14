inline void asteroid::update_Chebyshev(double time){

	int er;
	//Update the Chebyshev coefficients if necessary
	for(int p = 0; p < Nperturbers; ++p){

		const int pp = p * nCm * 3;
		if(time + time_reference > endTime_h[p]){
			for(int k = 0; k < 1000000; ++k){

				fseek(perturbersFile, offset0_h[p] * sizeof(double), SEEK_SET);

				er = fread(&startTime_h[p], sizeof(double), 1, perturbersFile);
				er = fread(&endTime_h[p], sizeof(double), 1, perturbersFile);


//printf(" ++ %d %d %d %d %.20g %.20g | %.20g %d\n", p, idp_h[p], offset0_h[p], nChebyshev_h[p], startTime_h[p], endTime_h[p], time + time_reference, er);

				er = fread(cdata_h + pp, sizeof(double), nChebyshev_h[p] * 3, perturbersFile);

				if(time + time_reference <= endTime_h[p]){
					break;
				}

				offset0_h[p] += nChebyshev_h[p] * 3 + 2;
			}
		}
		if(time + time_reference < startTime_h[p]){
			for(int k = 0; k < 1000000; ++k){

				fseek(perturbersFile, offset0_h[p] * sizeof(double), SEEK_SET);

				er = fread(&startTime_h[p], sizeof(double), 1, perturbersFile);
				er = fread(&endTime_h[p], sizeof(double), 1, perturbersFile);


//printf(" -- %d %d %d %d %.20g %.20g | %.20g %d\n", p, idp_h[p], offset0_h[p], nChebyshev_h[p], startTime_h[p], endTime_h[p], time + time_reference, er);

				er = fread(cdata_h + pp, sizeof(double), nChebyshev_h[p] * 3, perturbersFile);

				if(time + time_reference >= startTime_h[p]){
					break;
				}

				offset0_h[p] -= nChebyshev_h[p] * 3 + 2;
			}
		}
	}
}

inline void asteroid::update_perturbers(double time){

	double Tx[nCm];
	double Ty[nCm];
	double Tz[nCm];
	double Tvx[nCm];
	double Tvy[nCm];
	double Tvz[nCm];


	//Calculate positions of perturbers
	for(int p = 0; p < Nperturbers; ++p){
		const int nC = nChebyshev_h[p];

		const int pp = p * nCm * 3;

		double sizeSubInterval = endTime_h[p] - startTime_h[p];
//printf("time %d %.20g %.20g %.20g %.20g\n", p, time + time_reference, sizeSubInterval, startTime_h[p], endTime_h[p]);
		double subTime = (time_reference - startTime_h[p] + time) / sizeSubInterval;   //normalized time in  0 - 1
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
			xp += Tx[j] * cdata_h[pp + j];
			yp += Ty[j] * cdata_h[pp + nC + j];
			zp += Tz[j] * cdata_h[pp + 2 * nC + j];
//printf("Chebyshev %d %.20g %.20g %.20g %.20g %.20g\n", j, Tx[j], cdata_h[pp + j], xp, yp, zp);                    

			vxp += Tvx[j] * cdata_h[pp + j] * ct;
			vyp += Tvy[j] * cdata_h[pp + nC + j] * ct;
			vzp += Tvz[j] * cdata_h[pp + 2 * nC + j] * ct;
		}

		xTable_h[p] = xp;
		yTable_h[p] = yp;
		zTable_h[p] = zp;

		vxTable_h[p] = vxp;
		vyTable_h[p] = vyp;
		vzTable_h[p] = vzp;

//printf("positionA %d %.20g %.20g %.20g %.20g %.20g\n", p, time, xTable_h[p], yTable_h[p], zTable_h[p], t);
//printf("positionvA %d %.20g %.20g %.20g %.20g %.20g\n", p, time, vxTable_h[p], vyTable_h[p], vzTable_h[p], t);


//printf("positionB %d %.20g %.20g %.20g %.20g %.20g\n", p, time, xTable_h[p] / AUtokm, yTable_h[p] / AUtokm, zTable_h[p] / AUtokm, t);

//printf("%d %.20g %.20g %.20g %.20g %.20g\n", p, time, xTable_h[p], yTable_h[p], zTable_h[p], t);


	}

	//Calculate Earth and Moon positions, id 2 and 9
	//Up to here id = 2 is the Earth-Moon barycentrum, id = 9 is the geocentric position of the Moon        

	double xB = xTable_h[2];
	double yB = yTable_h[2];
	double zB = zTable_h[2];

	double vxB = vxTable_h[2];
	double vyB = vyTable_h[2];
	double vzB = vzTable_h[2];

	double xM = xTable_h[9];
	double yM = yTable_h[9];
	double zM = zTable_h[9];

	double vxM = vxTable_h[9];
	double vyM = vyTable_h[9];
	double vzM = vzTable_h[9];

	double f = 1.0/(1.0 + EM);

	xTable_h[2] = xB - xM * f;
	yTable_h[2] = yB - yM * f;
	zTable_h[2] = zB - zM * f;

	vxTable_h[2] = vxB - vxM * f;
	vyTable_h[2] = vyB - vyM * f;
	vzTable_h[2] = vzB - vzM * f;

	xTable_h[9] = xB + xM * EM * f;
	yTable_h[9] = yB + yM * EM * f;
	zTable_h[9] = zB + zM * EM * f;

	vxTable_h[9] = vxB + vxM * EM * f;
	vyTable_h[9] = vyB + vyM * EM * f;
	vzTable_h[9] = vzB + vzM * EM * f;

	for(int p = 0; p < Nperturbers; ++p){
		//positions are in km
		//velocities are in km/day
		xTable_h[p] /= AUtokm;
		yTable_h[p] /= AUtokm;
		zTable_h[p] /= AUtokm;

		//vxTable_h[p] /= AUtokm;
		//vyTable_h[p] /= AUtokm;
		//vzTable_h[p] /= AUtokm;
		//remove time factor again
		vxTable_h[p] /= AUtokm / 86400.0;
		vyTable_h[p] /= AUtokm / 86400.0;
		vzTable_h[p] /= AUtokm / 86400.0;
//printf("positionB %d %.20g %.20g %.20g %.20g\n", p, time, xTable_h[p], yTable_h[p], zTable_h[p]);

	}

//Planets Coordinates are now Barycentric equatorial
//Asteroids Coordinates are now Heliocentric equatorial
//Jupiter = 4
//Vesta = 26
//printf("%.20g %.20g %.20g %.20g\n", time + time_reference, xTable_h[26], yTable_h[26], zTable_h[26]);


	/*
	//print in the order of asssit
	for(int pp = 0; pp < Nperturbers; ++pp){
		int p = pp + 10;
		if(pp == 0) p = 11;     
		if(pp == 1) p = 10;     //Sun
		if(pp == 17) p = 8;     //Pluto
		if(pp == 18) p = 9;     //Moon
		if(pp == 19) p = 3;     //Mars
		if(pp == 20) p = 0;     //Mercurs
		if(pp == 21) p = 7;     //Neptune
		if(pp == 22) p = 6;     //Uranus
		if(pp == 23) p = 2;     //Earth
		if(pp == 24) p = 1;     //Venus
		if(pp == 25) p = 5;     //Saturn
		if(pp == 26) p = 4;     //Jupiter


		//printf("positionB %d %.20g %.20g %.20g %.20g\n", p, time, xTable_h[p], yTable_h[p], zTable_h[p]);
		//printf("positionvB %d %.20g %.20g %.20g %.20g\n", p, time, vxTable_h[p], vyTable_h[p], vzTable_h[p]);
		//printf("positionB %d %.20g %.20g %.20g %.20g\n", pp, time, xTable_h[pp], yTable_h[pp], zTable_h[pp]);
	}
	*/

	//Translate asteroid orbits from Heliocentric to Barycentric coordinates
	for(int p = 11; p < Nperturbers; ++p){
		xTable_h[p] += xTable_h[10];
		yTable_h[p] += yTable_h[10];
		zTable_h[p] += zTable_h[10];

		vxTable_h[p] += vxTable_h[10];
		vyTable_h[p] += vyTable_h[10];
		vzTable_h[p] += vzTable_h[10];
	}

//All Coordinates are now Barycentric equatorial

//printf("%.20g %.20g %.20g %.20g\n", time + time_reference, xTable_h[26], yTable_h[26], zTable_h[26]);

}

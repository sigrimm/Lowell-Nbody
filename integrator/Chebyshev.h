inline void asteroid::update_Chebyshev(double time){

	int er;
	//Update the Chebyshev coefficients if necessary
	for(int p = 0; p < Nperturbers; ++p){

		const int pp = p * nCm * 3;
		if(time + time_reference > endTime[p]){
			for(int k = 0; k < 1000000; ++k){
				//fscanf(perturbersFile, "%d", &id);
				//fscanf(perturbersFile, "%d", &nChebyshev);
				//fscanf(perturbersFile, "%lf", &startTime[p]);
				//fscanf(perturbersFile, "%lf", &endTime[p]);

				fseek(perturbersFile, offset0[p] * sizeof(double), SEEK_SET);

				er = fread(&startTime[p], sizeof(double), 1, perturbersFile);
				er = fread(&endTime[p], sizeof(double), 1, perturbersFile);


//printf(" ++ %d %d %d %d %.20g %.20g\n", p, id[p], offset0[p], nChebyshev[p], startTime[p], endTime[p]);

				er = fread(cdata + pp, sizeof(double), nChebyshev[p] * 3, perturbersFile);
				//for(int i = 0; i < nChebyshev[p] * 3; ++i){
				//      fscanf(perturbersFile, "%lf", &cdata[pp + i]);
				//      printf("%.20g ", cdata[pp + i]);
				//}
				//printf("\n");

				if(time + time_reference <= endTime[p]){
					break;
				}

				offset0[p] += nChebyshev[p] * 3 + 2;
			}
		}
		if(time + time_reference < startTime[p]){
			for(int k = 0; k < 1000000; ++k){
				//fscanf(perturbersFile, "%d", &id);
				//fscanf(perturbersFile, "%d", &nChebyshev);
				//fscanf(perturbersFile, "%lf", &startTime[p]);
				//fscanf(perturbersFile, "%lf", &endTime[p]);

				fseek(perturbersFile, offset0[p] * sizeof(double), SEEK_SET);

				er = fread(&startTime[p], sizeof(double), 1, perturbersFile);
				er = fread(&endTime[p], sizeof(double), 1, perturbersFile);


//printf(" -- %d %d %d %d %.20g %.20g\n", p, id[p], offset0[p], nChebyshev[p], startTime[p], endTime[p]);

				er = fread(cdata + pp, sizeof(double), nChebyshev[p] * 3, perturbersFile);
				//for(int i = 0; i < nChebyshev[p] * 3; ++i){
				//      fscanf(perturbersFile, "%lf", &cdata[pp + i]);
				//      printf("%.20g ", cdata[pp + i]);
				//}
				//printf("\n");

				if(time + time_reference >= startTime[p]){
					break;
				}

				offset0[p] -= nChebyshev[p] * 3 + 2;
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
		const int nC = nChebyshev[p];

		const int pp = p * nCm * 3;

		double sizeSubInterval = endTime[p] - startTime[p];
//if(p > 10) printf("time %d %.20g %.20g %.20g %.20g\n", p, time + time_reference, sizeSubInterval, startTime[p], endTime[p]);
		double subTime = (time_reference - startTime[p] + time) / sizeSubInterval;   //normalized time in  0 - 1
		double t = 2.0 * subTime - 1.0;                         //mormalized time in -1 - 1

		//double ct = 2.0 / sizeSubInterval;                    //correction factor for time units

		//This is sone in Assist, remove the time factor later
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
		//for(int j = nChebyshev[p] - 1; j >= 0; --j){    //reduce floating point errors by revert order
			xp += Tx[j] * cdata[pp + j];
			yp += Ty[j] * cdata[pp + nC + j];
			zp += Tz[j] * cdata[pp + 2 * nC + j];
//printf("Chebyshev %d %.20g %.20g %.20g\n", j, Tx[j], cdata[pp + j], x[p]);                    

			vxp += Tvx[j] * cdata[pp + j] * ct;
			vyp += Tvy[j] * cdata[pp + nC + j] * ct;
			vzp += Tvz[j] * cdata[pp + 2 * nC + j] * ct;
		}

		x[p] = xp;
		y[p] = yp;
		z[p] = zp;

		vx[p] = vxp;
		vy[p] = vyp;
		vz[p] = vzp;

//printf("positionA %d %.20g %.20g %.20g %.20g %.20g\n", p, time, x[p], y[p], z[p], t);
//printf("positionvA %d %.20g %.20g %.20g %.20g %.20g\n", p, time, vx[p], vy[p], vz[p], t);


//printf("positionB %d %.20g %.20g %.20g %.20g %.20g\n", p, time, x[p] / AUtokm, y[p] / AUtokm, z[p] / AUtokm, t);

//printf("%d %.20g %.20g %.20g %.20g %.20g\n", p, time, x[p], y[p], z[p], t);


	}

	//Calculate Earth and Moon positions, id 2 and 9
	//Up to here id = 2 is the Earth-Moon barycentrum, id = 9 is the geocentric position of the Moon        

	double xB = x[2];
	double yB = y[2];
	double zB = z[2];

	double vxB = vx[2];
	double vyB = vy[2];
	double vzB = vz[2];

	double xM = x[9];
	double yM = y[9];
	double zM = z[9];

	double vxM = vx[9];
	double vyM = vy[9];
	double vzM = vz[9];

	double f = 1.0/(1.0 + EM);

	x[2] = xB - xM * f;
	y[2] = yB - yM * f;
	z[2] = zB - zM * f;

	vx[2] = vxB - vxM * f;
	vy[2] = vyB - vyM * f;
	vz[2] = vzB - vzM * f;

	x[9] = xB + xM * EM * f;
	y[9] = yB + yM * EM * f;
	z[9] = zB + zM * EM * f;

	vx[9] = vxB + vxM * EM * f;
	vy[9] = vyB + vyM * EM * f;
	vz[9] = vzB + vzM * EM * f;

	for(int p = 0; p < Nperturbers; ++p){
		//positions are in km
		//velocities are in km/day
		x[p] /= AUtokm;
		y[p] /= AUtokm;
		z[p] /= AUtokm;

		//vx[p] /= AUtokm;
		//vy[p] /= AUtokm;
		//vz[p] /= AUtokm;
		//remove time factor again
		vx[p] /= AUtokm / 86400.0;
		vy[p] /= AUtokm / 86400.0;
		vz[p] /= AUtokm / 86400.0;
//printf("positionB %d %.20g %.20g %.20g %.20g\n", p, time, x[p], y[p], z[p]);

	}

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


		//printf("positionB %d %.20g %.20g %.20g %.20g\n", p, time, x[p], y[p], z[p]);
		//printf("positionvB %d %.20g %.20g %.20g %.20g\n", p, time, vx[p], vy[p], vz[p]);
		//printf("positionB %d %.20g %.20g %.20g %.20g\n", pp, time, x[pp], y[pp], z[pp]);
	}
	*/

	//Translate asteroid orbits from Heliocentric to Barycentric coordinates
	//This is done so in Assist, but probably this is not correct, check!
	for(int p = 11; p < Nperturbers; ++p){
		x[p] += x[10];
		y[p] += y[10];
		z[p] += z[10];
	}

}

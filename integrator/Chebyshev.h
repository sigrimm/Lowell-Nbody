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

		x_h[p] = xp;
		y_h[p] = yp;
		z_h[p] = zp;

		vx_h[p] = vxp;
		vy_h[p] = vyp;
		vz_h[p] = vzp;

//printf("positionA %d %.20g %.20g %.20g %.20g %.20g\n", p, time, x_h[p], y_h[p], z_h[p], t);
//printf("positionvA %d %.20g %.20g %.20g %.20g %.20g\n", p, time, v_hx[p], vy_h[p], vz_h[p], t);


//printf("positionB %d %.20g %.20g %.20g %.20g %.20g\n", p, time, x_h[p] / AUtokm, y_h[p] / AUtokm, z_h[p] / AUtokm, t);

//printf("%d %.20g %.20g %.20g %.20g %.20g\n", p, time, x_h[p], y_h[p], z_h[p], t);


	}

	//Calculate Earth and Moon positions, id 2 and 9
	//Up to here id = 2 is the Earth-Moon barycentrum, id = 9 is the geocentric position of the Moon        

	double xB = x_h[2];
	double yB = y_h[2];
	double zB = z_h[2];

	double vxB = vx_h[2];
	double vyB = vy_h[2];
	double vzB = vz_h[2];

	double xM = x_h[9];
	double yM = y_h[9];
	double zM = z_h[9];

	double vxM = vx_h[9];
	double vyM = vy_h[9];
	double vzM = vz_h[9];

	double f = 1.0/(1.0 + EM);

	x_h[2] = xB - xM * f;
	y_h[2] = yB - yM * f;
	z_h[2] = zB - zM * f;

	vx_h[2] = vxB - vxM * f;
	vy_h[2] = vyB - vyM * f;
	vz_h[2] = vzB - vzM * f;

	x_h[9] = xB + xM * EM * f;
	y_h[9] = yB + yM * EM * f;
	z_h[9] = zB + zM * EM * f;

	vx_h[9] = vxB + vxM * EM * f;
	vy_h[9] = vyB + vyM * EM * f;
	vz_h[9] = vzB + vzM * EM * f;

	for(int p = 0; p < Nperturbers; ++p){
		//positions are in km
		//velocities are in km/day
		x_h[p] /= AUtokm;
		y_h[p] /= AUtokm;
		z_h[p] /= AUtokm;

		//vx_h[p] /= AUtokm;
		//vy_h[p] /= AUtokm;
		//vz_h[p] /= AUtokm;
		//remove time factor again
		vx_h[p] /= AUtokm / 86400.0;
		vy_h[p] /= AUtokm / 86400.0;
		vz_h[p] /= AUtokm / 86400.0;
//printf("positionB %d %.20g %.20g %.20g %.20g\n", p, time, x_h[p], y_h[p], z_h[p]);

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


		//printf("positionB %d %.20g %.20g %.20g %.20g\n", p, time, x_h[p], y_h[p], z_h[p]);
		//printf("positionvB %d %.20g %.20g %.20g %.20g\n", p, time, vx_h[p], vy_h[p], vz_h[p]);
		//printf("positionB %d %.20g %.20g %.20g %.20g\n", pp, time, x_h[pp], y_h[pp], z_h[pp]);
	}
	*/

	//Translate asteroid orbits from Heliocentric to Barycentric coordinates
	//This is done so in Assist, but probably this is not correct, check!
	for(int p = 11; p < Nperturbers; ++p){
		x_h[p] += x_h[10];
		y_h[p] += y_h[10];
		z_h[p] += z_h[10];
	}

}

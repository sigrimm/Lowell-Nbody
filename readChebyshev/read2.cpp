#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(){


	FILE *infile;
	//infile = fopen("PerturbersChebyshev.dat", "r");
	infile = fopen("PerturbersChebyshev.bin", "rb");


	int Nperturbers = 27;
	//double time = 2450800.5;
	double time_reference = 2451545.0;
	double time = 8255.495;
	double timeEnd = 6820.37;
	double dt = -0.01;
	double time0;	//start time of the data file
	double time1;	//end time of the data file

	double startTime[Nperturbers];
	double endTime[Nperturbers];
	int id[Nperturbers];
	int nChebyshev[Nperturbers];
	int nCm = 0;		//Maximum number of Chebyshev coefficients
	int offset0[Nperturbers];
	int offset1[Nperturbers];
	double GM[Nperturbers];

	double AUtokm;		//AU to km
	double EM;		//Earth to moon mass ratio
	double CLIGHT;		//speed of light
	double RE;		//Radius of Earth
	double J2E;		//J2 of Earth

	int N = Nperturbers + 1;

	double x[N];
	double y[N];
	double z[N];

	double vx[N];
	double vy[N];
	double vz[N];

	double ax[N];
	double ay[N];
	double az[N];

	double A1[N];
	double A2[N];
	double A3[N];


	// *********************************
	//Set initial conditions

        //Barycentric coordinates of 15TC25   
        x[Nperturbers + 0] = -0.621994840672447701912517459277;   // x in AU
        y[Nperturbers + 0] = -0.828032662228602056586623803014;  // y
        z[Nperturbers + 0] = -0.406943193813317449780697643291;  // z
        vx[Nperturbers + 0] = 0.0122692503194383149833779356186;   // vx in AU/day
        vy[Nperturbers + 0] = -0.00859370367531481910150503722434;   // vy
        vz[Nperturbers + 0] = -0.0046654983615674223973446288482;  // vz

	//Non-grav terms
	A1[Nperturbers + 0] = 1.599229127169e-10;
	A2[Nperturbers + 0] = -5.273644346744e-12;
	A3[Nperturbers + 0] = 0.0;



	// *********************************


	fread(&time0, sizeof(double), 1, infile);
	fread(&time1, sizeof(double), 1, infile);
	fread(&AUtokm, sizeof(double), 1, infile);
	fread(&EM, sizeof(double), 1, infile);
	fread(&CLIGHT, sizeof(double), 1, infile);
	fread(&RE, sizeof(double), 1, infile);
	fread(&J2E, sizeof(double), 1, infile);

	for(int i = 0; i < Nperturbers; ++i){
		fread(&id[i], sizeof(int), 1, infile);
		fread(&nChebyshev[i], sizeof(int), 1, infile);
		fread(&offset0[i], sizeof(int), 1, infile);
		fread(&offset1[i], sizeof(int), 1, infile);
		fread(&GM[i], sizeof(double), 1, infile);

		nCm = (nCm > nChebyshev[i]) ? nCm : nChebyshev[i];

		offset0[i] += (3 * Nperturbers + 7);	//add size of header 7*double + Nperturbers * (4 int + double)
		offset1[i] += (3 * Nperturbers + 7);	//add size of header
	
		startTime[i] = 100000000.0;	//large number
		endTime[i] = 0;
		printf("%d %d %d %d %.20g\n", id[i], nChebyshev[i], offset0[i], offset1[i], GM[i]);
	}


	double cdata[Nperturbers * nCm * 3];

	double Tx[nCm];
	double Ty[nCm];
	double Tz[nCm];
	double Tvx[nCm];
	double Tvy[nCm];
	double Tvz[nCm];

	
	for(int tt = 0; tt < 200000; ++tt){


		for(int p = 0; p < Nperturbers; ++p){

			int pp = p * nCm * 3;
			if(time + time_reference > endTime[p]){
				for(int k = 0; k < 1000000; ++k){
					//fscanf(infile, "%d", &id);
					//fscanf(infile, "%d", &nChebyshev);
					//fscanf(infile, "%lf", &startTime[p]);
					//fscanf(infile, "%lf", &endTime[p]);
		
					fseek(infile, offset0[p] * sizeof(double), SEEK_SET);

					fread(&startTime[p], sizeof(double), 1, infile);
					fread(&endTime[p], sizeof(double), 1, infile);

					
//printf(" ++ %d %d %d %d %.20g %.20g\n", p, id[p], offset0[p], nChebyshev[p], startTime[p], endTime[p]);

					fread(cdata + pp, sizeof(double), nChebyshev[p] * 3, infile);
					for(int i = 0; i < nChebyshev[p] * 3; ++i){
						//fscanf(infile, "%lf", &cdata[pp + i]);
						//printf("%.20g ", cdata[pp + i]);
					}
					//printf("\n");
				

					if(time + time_reference <= endTime[p]){
						break;
					}

					offset0[p] += nChebyshev[p] * 3 + 2;
				}
			}
			if(time + time_reference < startTime[p]){
				for(int k = 0; k < 1000000; ++k){
					//fscanf(infile, "%d", &id);
					//fscanf(infile, "%d", &nChebyshev);
					//fscanf(infile, "%lf", &startTime[p]);
					//fscanf(infile, "%lf", &endTime[p]);
		
					fseek(infile, offset0[p] * sizeof(double), SEEK_SET);

					fread(&startTime[p], sizeof(double), 1, infile);
					fread(&endTime[p], sizeof(double), 1, infile);

					
//printf(" -- %d %d %d %d %.20g %.20g\n", p, id[p], offset0[p], nChebyshev[p], startTime[p], endTime[p]);

					fread(cdata + pp, sizeof(double), nChebyshev[p] * 3, infile);
					for(int i = 0; i < nChebyshev[p] * 3; ++i){
						//fscanf(infile, "%lf", &cdata[pp + i]);
						//printf("%.20g ", cdata[pp + i]);
					}
					//printf("\n");
				

					if(time + time_reference >= startTime[p]){
						break;
					}

					offset0[p] -= nChebyshev[p] * 3 + 2;
				}
			}


			double sizeSubInterval = endTime[p] - startTime[p];
//if(p > 10) printf("time %d %.20g %.20g %.20g %.20g\n", p, time + time_reference, sizeSubInterval, startTime[p], endTime[p]);
			double subTime = (time_reference - startTime[p] + time) / sizeSubInterval;   //normalized time in  0 - 1
			double t = 2.0 * subTime - 1.0;                         //mormalized time in -1 - 1

			//double ct = 2.0 / sizeSubInterval;			//correction factor for time units

			//This is sone in Assist, remove the time factor later
			double ct = 2.0 / sizeSubInterval / 86400.0;		//correction factor for time units

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

			for(int j = 2; j < nChebyshev[p]; ++j){
				Tx[j] = 2.0 * t * Tx[j - 1] - Tx[j - 2];
				Ty[j] = 2.0 * t * Ty[j - 1] - Ty[j - 2];
				Tz[j] = 2.0 * t * Tz[j - 1] - Tz[j - 2];
				Tvx[j] = 2.0 * t * Tvx[j - 1] + 2.0 * Tx[j - 1] - Tvx[j - 2];
				Tvy[j] = 2.0 * t * Tvy[j - 1] + 2.0 * Ty[j - 1] - Tvy[j - 2];
				Tvz[j] = 2.0 * t * Tvz[j - 1] + 2.0 * Tz[j - 1] - Tvz[j - 2];
			}

			x[p] = 0.0;
			y[p] = 0.0;
			z[p] = 0.0;

			vx[p] = 0.0;
			vy[p] = 0.0;
			vz[p] = 0.0;

			for(int j = 0; j < nChebyshev[p]; ++j){    //reduce floating point errors by revert order
			//for(int j = nChebyshev[p] - 1; j >= 0; --j){    //reduce floating point errors by revert order
				x[p] += Tx[j] * cdata[pp + j];
				y[p] += Ty[j] * cdata[pp + nChebyshev[p] + j];
				z[p] += Tz[j] * cdata[pp + 2 * nChebyshev[p] + j];
//printf("Chebyshev %d %.20g %.20g %.20g\n", j, Tx[j], cdata[pp + j], x[p]);			

				vx[p] += Tvx[j] * cdata[pp + j] * ct;
				vy[p] += Tvy[j] * cdata[pp + nChebyshev[p] + j] * ct;
				vz[p] += Tvz[j] * cdata[pp + 2 * nChebyshev[p] + j] * ct;
			}

//printf("positionA %d %.20g %.20g %.20g %.20g %.20g\n", p, time, x[p], y[p], z[p], t);
//printf("positionvA %d %.20g %.20g %.20g %.20g %.20g\n", p, time, vx[p], vy[p], vz[p], t);

			
//printf("positionB %d %.20g %.20g %.20g %.20g %.20g\n", p, time, x[p] / AUtokm, y[p] / AUtokm, z[p] / AUtokm, t);

			//printf("%d %.20g %.20g %.20g %.20g %.20g\n", p, time, x[p], y[p], z[p], t);
			

		}

		//Calculate Earth and Moon positions, id 2 and 9
		//Up to hear id = 2 is the Earth-Moon barycentrum, id = 9 is the geocentric position of the Moon	
		
		double xB = x[2];
		double yB = y[2];
		double zB = z[2];

		double vxB = x[2];
		double vyB = y[2];
		double vzB = z[2];

		double xM = x[9];
		double yM = y[9];
		double zM = z[9];

		double vxM = x[9];
		double vyM = y[9];
		double vzM = z[9];

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

//			printf("%d %.20g %.20g %.20g %.20g\n", p, time, x[p], y[p], z[p]);

		}
		//print in the order of asssit
		for(int pp = 0; pp < Nperturbers; ++pp){
			int p = pp + 10;
			if(pp == 0) p = 11;	
			if(pp == 1) p = 10;	//Sun
			if(pp == 17) p = 8;	//Pluto
			if(pp == 18) p = 9;	//Moon
			if(pp == 19) p = 3;	//Mars
			if(pp == 20) p = 0;	//Mercurs
			if(pp == 21) p = 7;	//Neptune
			if(pp == 22) p = 6;	//Uranus
			if(pp == 23) p = 2;	//Earth
			if(pp == 24) p = 1;	//Venus
			if(pp == 25) p = 5;	//Saturn
			if(pp == 26) p = 4;	//Jupiter


//printf("positionB %d %.20g %.20g %.20g %.20g\n", p, time, x[p], y[p], z[p]);
printf("positionvB %d %.20g %.20g %.20g %.20g\n", p, time, vx[p], vy[p], vz[p]);
//printf("positionB %d %.20g %.20g %.20g %.20g\n", pp, time, x[pp], y[pp], z[pp]);
		}

		//Translate asteroid orbits from Heliocentric to Barycentric coordinates
		//This is done so in Assist, but probably this is not correct, check!
		for(int p = 11; p < Nperturbers; ++p){
			x[p] += x[10];
			y[p] += y[10];
			z[p] += z[10];
		}

		//Drift
		for(int p = Nperturbers; p < N; ++p){
printf("Drift %d %.20g %.20g %.20g %.20g\n", p, x[p], vx[p], dt, 0.5* dt * vx[p]);
			x[p] += 0.5* dt * vx[p];
			y[p] += 0.5* dt * vy[p];
			z[p] += 0.5* dt * vz[p];
		}
		time += dt / 2.0;	//Split in two steps because of the leapfrog integrator
	

		
		//Kick
		// ----------------------------------------------------------------------------
		for(int i = Nperturbers; i < N; ++i){
			ax[i] = 0.0;
			ay[i] = 0.0;
			az[i] = 0.0;
		}
		// ----------------------------------------------------------------------------
		//Non Grav
		for(int i = Nperturbers; i < N; ++i){

			//Translate to heliocentric coordinates
			double xi = x[i] - x[10];
			double yi = y[i] - y[10];
			double zi = z[i] - z[10];

			double vxi = vx[i] - vx[10];
			double vyi = vy[i] - vy[10];
			double vzi = vz[i] - vz[10];

			double rsq = xi * xi + yi * yi + zi * zi;
			double r = sqrt(rsq);
//printf("r %.20g %.20g %.20g %.20g\n", xi, yi ,zi, r);
//printf("v %.20g %.20g %.20g %.20g %.20g\n", vxi, vyi, vzi, vx[i], vx[10]);

			//angular momenrum h = r x v
			double hx = yi * vzi - zi * vyi;
			double hy =-xi * vzi + zi * vxi;
			double hz = xi * vyi - yi * vxi;

			double hsq = hx * hx + hy * hy + hz * hz;
			double h = sqrt(hsq);

			//Transverse velocity t = h x r
			double tx = hy * zi - hz * yi;
			double ty =-hx * zi + hz * xi;
			double tz = hx * yi - hy * xi;

			double tsq = tx * tx + ty * ty + tz * tz;
			double t = sqrt(tsq);

			//double gr = 1.0 / rsq;  //only valid for asteroids, not for comets 
			double alpha = 1.0;
			double nk = 0.0;
			double nm = 2.0;
			double nn = 5.093;
			double r0 = 1.0;

			const double gr = alpha*pow(r/r0, -nm) * pow(1.0 + pow(r / r0, nn), -nk);


			/*
			double rr = r / R0[i];
			double g1 = pow(rr, -NM[i]);
			double g2 = pow(rr, Nn[i]);
			double g3 = pow(1.0 + g2, -NK[i]);
			double gr = ALN[i] * g1 * g3;
			//Larry uses constant values for these parameters, for all comets the same
			//add time delay for comets

			printf("gr %.20g %.20g\n", gr1, gr);
			*/

			double f1 = A1[i] * gr / r;
			double f2 = A2[i] * gr / t;
			double f3 = A3[i] * gr / h;

			//printf("NonGrav  %.20g %.20g %.20g %.20g |%.20g %.20g %.20g\n", gr, r, t, h, f1, f2, f3);
//printf("NonGrav a %.20g %.20g %.20g\n", ax[i], ay[i], az[i]);

			//ax[i] += f1 * xi + f2 * tx + f3 * hx;
			//ay[i] += f1 * yi + f2 * ty + f3 * hy;
			//az[i] += f1 * zi + f2 * tz + f3 * hz;
			ax[i] += A1[i] * gr * xi / r + A2[i] * gr * tx / t + A3[i] * gr * hx / h;
			ay[i] += A1[i] * gr * yi / r + A2[i] * gr * ty / t + A3[i] * gr * hy / h;
			az[i] += A1[i] * gr * zi / r + A2[i] * gr * tz / t + A3[i] * gr * hz / h;

//printf("NonGrav %.20g %.20g %.20g %.20g | %.20g %.20g %.20g\n", gr, r, t, h, f1, f2, f3);
printf("NonGrav a %.20g %.20g %.20g\n", ax[i], ay[i], az[i]);

		}
		// ----------------------------------------------------------------------------
		// ----------------------------------------------------------------------------
		//GR
		for(int i = Nperturbers; i < N; ++i){
			double c = (CLIGHT / AUtokm) * 86400.0 ;
//printf("c %.30g %.30g\n", c, CLIGHT);
			double c2 = c * c;

			//Translate to heliocentric coordinates
			double xi = x[i] - x[10];
			double yi = y[i] - y[10];
			double zi = z[i] - z[10];

			double vxi = vx[i] - vx[10];
			double vyi = vy[i] - vy[10];
			double vzi = vz[i] - vz[10];

			double rsq = xi * xi + yi * yi + zi * zi;
			double r = sqrt(rsq);
			double vsq = vxi * vxi + vyi * vyi + vzi * vzi;

			double rv = xi * vxi + yi * vyi + zi * vzi;

			double f1 = GM[10] / (r * r * r * c2);
			double t1 = 4.0 * GM[10] / r - vsq;
			double t3 = 4.0 * rv;
//printf("GRa %.20g %.20g %.20g %.20g %.20g\n", vsq, r, t1, t3, f1);
			//printf("a %d %.20g %.20g %.20g\n", i, ax, ay, az);

			//printf("A %d %.20g %.20g %.20g %.20g %.20g %.20g\n", i, xi, yi, zi, vxi, vyi, vzi);
			//printf("B %d %.20g %.20g %.20g\n", i, f1, t1, t3);

			double aax = f1 * (t1 * xi + t3 * vxi);
			double aay = f1 * (t1 * yi + t3 * vyi);
			double aaz = f1 * (t1 * zi + t3 * vzi);

//printf("GR a %.20g %.20g %.20g\n", aax, aay, aaz);
			ax[i] += f1 * (t1 * xi + t3 * vxi);
			ay[i] += f1 * (t1 * yi + t3 * vyi);
			az[i] += f1 * (t1 * zi + t3 * vzi);
printf("GR a %d %.20g %.20g %.20g\n", i, ax[i], ay[i], az[i]);



		}
		// ----------------------------------------------------------------------------
		// ----------------------------------------------------------------------------
		//J2
		for(int i = Nperturbers; i < N; ++i){


			RE /= AUtokm;   //Earth radius in AU

printf("J2 %.20g %.20g %.20g\n", J2E, RE, RE * AUtokm);

			double xE = x[i] - x[2];
			double yE = y[i] - y[2];
			double zE = z[i] - z[2];
printf("dz %.20g\n", zE);

			double rsq = xE * xE + yE * yE + zE * zE;
			double r = sqrt(rsq);
			double r5 = rsq * rsq * r;

			//double t1 = GM[2] * 3.0 * J2E * RE * RE / (2.0 * r5);
			double t1 = 3.0 * J2E * RE * RE / rsq / rsq / r / 2.0;
			double t2 = 5.0 * (zE * zE / rsq) - 1.0;

printf("%.20g %.20g %.20g %.20g\n", GM[2], t1, t2, 3.0 * J2E * RE * RE/rsq/rsq/r/2.0);

			double tx = GM[2] * t1 * t2 * xE;
			double ty = GM[2] * t1 * t2 * yE;
			double tz = GM[2] * t1 * (t2 - 2.0) * zE;

			ax[i] += tx;
			ay[i] += ty;
			az[i] += tz;

printf("J2 a %.20g %.20g %.20g\n", tx, ty, tz);
		}
		// ----------------------------------------------------------------------------

		// ----------------------------------------------------------------------------
		//Gravity
		for(int i = Nperturbers; i < N; ++i){
			for(int pp = 0; pp < Nperturbers; ++pp){
				int p = pp + 11;
				if(pp == 16) p = 8;	//Pluto
				if(pp == 17) p = 9;	//Moon
				if(pp == 18) p = 3;	//Mars
				if(pp == 19) p = 0;	//Mercurs
				if(pp == 20) p = 7;	//Neptune
				if(pp == 21) p = 6;	//Uranus
				if(pp == 22) p = 2;	//Earth
				if(pp == 23) p = 1;	//Venus
				if(pp == 24) p = 5;	//Saturn
				if(pp == 25) p = 4;	//Jupiter
				if(pp == 26) p = 10;	//Sun

				double dx = x[i] - x[p];
				double dy = y[i] - y[p];
				double dz = z[i] - z[p];
				double rsq = dx*dx + dy*dy + dz*dz;
				double r = sqrt(rsq);
				double s = GM[p] / (r * r * r);

				ax[i] -= s*dx;
				ay[i] -= s*dy;
				az[i] -= s*dz;


printf("position %d %d %.20g %.20g %.20g %.20g %.20g\n", pp, p, time, x[p], y[p], z[p], GM[p]);
			}
		}
		// ----------------------------------------------------------------------------

		for(int i = Nperturbers; i < N; ++i){
printf("Kick %d %.20g %.20g %.20g %.20g\n", i, vx[i], ax[i], dt, dt * ax[i]);
			vx[i] += dt * ax[i];
			vy[i] += dt * ay[i];
			vz[i] += dt * az[i];
		}
		//Drift
		for(int p = Nperturbers; p < N; ++p){
printf("Drift %d %.20g %.20g %.20g %.20g\n", p, x[p], vx[p], dt, 0.5* dt * vx[p]);
			x[p]  += 0.5* dt * vx[p];
			y[p]  += 0.5* dt * vy[p];
			z[p]  += 0.5* dt * vz[p];
		}
		time += dt / 2.0;

		if(time + time_reference > time1 || time + time_reference < time0){
			printf("Reached the end of the Chebyshev data file\n");
			return 0;
		}

		if(time < timeEnd){
			printf("Reached the end of the integration\n");
			return 0;
		}

	}

	fclose(infile);
	
}


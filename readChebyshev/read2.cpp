#include <stdio.h>
#include <stdlib.h>


int main(){


	FILE *infile;
	//infile = fopen("PerturbersChebyshev.dat", "r");
	infile = fopen("PerturbersChebyshev.bin", "rb");


	int N = 27;
	//double time = 2450800.5;
	double time_reference = 2451545.0;
	double time = 8255.495;
	double timeEnd = 6820.37;
	double dt = -0.01;
	double time0;	//start time of the data file
	double time1;	//end time of the data file

	double startTime[N];
	double endTime[N];
	int id[N];
	int nChebyshev[N];
	int nCm = 0;		//Maximum number of Chebyshev coefficients
	int offset0[N];
	int offset1[N];
	double AUtokm;		//AU to km
	double EM;		//Earth to moon mass ratio

	double x[N];
	double y[N];
	double z[N];

	double vx[N];
	double vy[N];
	double vz[N];

	fread(&time0, sizeof(double), 1, infile);
	fread(&time1, sizeof(double), 1, infile);
	fread(&AUtokm, sizeof(double), 1, infile);
	fread(&EM, sizeof(double), 1, infile);

	for(int i = 0; i < N; ++i){
		fread(&id[i], sizeof(int), 1, infile);
		fread(&nChebyshev[i], sizeof(int), 1, infile);
		fread(&offset0[i], sizeof(int), 1, infile);
		fread(&offset1[i], sizeof(int), 1, infile);

		nCm = (nCm > nChebyshev[i]) ? nCm : nChebyshev[i];

		offset0[i] += (2*N + 4);	//add size of header
		offset1[i] += (2*N + 4);	//add size of header
	
		startTime[i] = 100000000.0;	//large number
		endTime[i] = 0;
		//printf("%d %d %d %d\n", id[i], nChebyshev[i], offset0[i], offset1[i]);
	}


	double cdata[N * nCm * 3];

	double Tx[nCm];
	double Ty[nCm];
	double Tz[nCm];
	double Tvx[nCm];
	double Tvy[nCm];
	double Tvz[nCm];

	
	for(int tt = 0; tt < 200000; ++tt){


		for(int p = 0; p < N; ++p){

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

				vx[p] += Tvx[j] * cdata[pp + j];       // * c
				vy[p] += Tvy[j] * cdata[pp + nChebyshev[p] + j];       // * c
				vz[p] += Tvz[j] * cdata[pp + 2 * nChebyshev[p] + j];   // * c
			}

//printf("positionA %d %.20g %.20g %.20g %.20g %.20g\n", p, time, x[p], y[p], z[p], t);

			
//printf("positionB %d %.20g %.20g %.20g %.20g %.20g\n", p, time, x[p] / AUtokm, y[p] / AUtokm, z[p] / AUtokm, t);

			//printf("%d %.20g %.20g %.20g %.20g %.20g\n", p, time, x[p], y[p], z[p], t);
			

		}

		//Calculate Earth and Moon positions, id 2 and 9
		//Up to hear id = 2 is the Earth-Moon barycentrum, id = 9 is the geocentric position of the Moon	
		
		double xB = x[2];
		double yB = y[2];
		double zB = z[2];

		double xM = x[9];
		double yM = y[9];
		double zM = z[9];

		double f = 1.0/(1.0 + EM);

		x[2] = xB - xM * f;
		y[2] = yB - yM * f;
		z[2] = zB - zM * f;

		x[9] = xB + xM * EM * f;
		y[9] = yB + yM * EM * f;
		z[9] = zB + zM * EM * f;

		for(int p = 0; p < N; ++p){
			x[p] /= AUtokm;
			y[p] /= AUtokm;
			z[p] /= AUtokm;
//printf("positionB %d %.20g %.20g %.20g %.20g\n", p, time, x[p], y[p], z[p]);

//			printf("%d %.20g %.20g %.20g %.20g\n", p, time, x[p], y[p], z[p]);

		}
		//print in the order of asssit
		for(int pp = 0; pp < N; ++pp){
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


printf("positionB %d %.20g %.20g %.20g %.20g\n", p, time, x[p], y[p], z[p]);
//printf("positionB %d %.20g %.20g %.20g %.20g\n", pp, time, x[pp], y[pp], z[pp]);
		}

		time += dt / 2.0;	//Split in to steps because of the leapfrog integrator
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


#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define dayUnit 0.01720209895
//#define dayUnit 0.01720209894846
//#define dayUnit 1.0
//Erics value:
//#define dayUnit 0.01720412578565452474399499749324604636058


void acc(double *m, double *x, double *y, double *z, double &ax, double &ay, double &az, int i, int j){
	double rx = x[j] - x[i];
	double ry = y[j] - y[i];
	double rz = z[j] - z[i];
	double rsq = rx * rx + ry * ry + rz * rz;
	double r = sqrt(rsq);
	
	double k = 1.0;
	//use tis if dayUnit = 1.0
	//double k = 0.01720209895 * 0.01720209895;
	//double k = 0.000295981944048623;

	ax += k * m[j] / (r * rsq) * rx;
	ay += k * m[j] / (r * rsq) * ry;
	az += k * m[j] / (r * rsq) * rz;
}


int main(){

	//Number of planets
	
	const int NN = 21;

	int outHelio = 1;
	//1 print output in heliocentric coordinates
	//0 print output in barycentric  coordinates

	long long int Nsteps = 40000;	
	long long int outInterval = 10;
	double dt = 0.1 * dayUnit;
	
	int id[NN];
	double m[NN];
	double rad[NN];
	double x[NN];
	double y[NN];
	double z[NN];
	double vx[NN];
	double vy[NN];
	double vz[NN];

	FILE *outfile;
	char outfilename[160];	

	FILE *infile;
	char infilename[160];

	sprintf(infilename, "initial.dat");
	infile = fopen(infilename, "r");

	double time = 0.0;

	//Units are 1/(mass of object in solar masses)
	double pmass[] = {
	    6023682.155592479e0,      // Mercury    (1)
	    408523.7186582996e0,      // Venus      (2)
	    332946.0488339480e0,      // Earth      (3)
	    3098703.590290707e0,      // Mars       (4)
	    1047.348625463337e0,      // Jupiter    (5)
	    3497.901767786633e0,      // Saturn     (6)
	    22902.98161308703e0,      // Uranus     (7)
	    19412.25977597307e0,      // Neptune    (8)
	    135836683.7686175e0,      // Pluto      (9)
	    2.11e09,                  // Ceres      (10)
	    9.83e09,                  // Pallas     (11)
	    7.41e10,                  // Juno       (12)
	    7.66e09,                  // Vesta      (13)
	    1.89e10,                  // Hygiea     (14)
	    7.63e10,                  // Eunomia    (15)
	    9.09e10,                  // Euphrosyne (16)
	    5.98e10,                  // Europa     (17)
	    6.74e10,                  // Davida     (18)
	    5.50e10,                  // Interamnia (19)
	    27068703.24120323e0,      // Moon       (20)
	    1.000000000000000e0       // Sun        (21)
	};

	printf("%.30g\n", 1.0/pmass[2]);

	int N = 1;

	//Sun
	id[0] = 20;
	m[0] = 1.0;
	rad[0] = 1.0e-3; //?
	x[0] = 0.0;
	y[0] = 0.0;
	z[0] = 0.0;
	vx[0] = 0.0;
	vy[0] = 0.0;
	vz[0] = 0.0;

	for(int i = 1; i < NN; ++i){
		id[i] = -1;
		m[i] = 0.0;
		rad[i] = 0.0;
		x[i] = 0.0;
		y[i] = 0.0;
		z[i] = 0.0;
		vx[i] = 0.0;
		vy[i] = 0.0;
		vz[i] = 0.0;
	}
	for(int i = 1; i < NN; ++i){
		int er = 0;
		fscanf(infile, "%lf", &time);
		fscanf(infile, "%d",  &id[i]);
		fscanf(infile, "%lf", &x[i]);
		fscanf(infile, "%lf", &y[i]);
		fscanf(infile, "%lf", &z[i]);
		fscanf(infile, "%lf", &vx[i]);
		fscanf(infile, "%lf", &vy[i]);
		er = fscanf(infile, "%lf", &vz[i]);
		if(er < 0) break;
		++N;
printf("er %d %d %d %d\n", i, id[i], er, N);


		/*
		//add perturbation to check convergence
		double eps = 1.0e-12;
		x[i] += eps;
		y[i] += eps;
		z[i] += eps;
		vx[i] += eps;
		vy[i] += eps;
		vz[i] += eps;
		*/

	}
	fclose(infile);

	for(int i = 1; i < N; ++i){
		m[i] = 1.0/pmass[id[i] -1];
		//m[i] += m[i] * 1.0e-6;

		//if(id[i] >= 10 && id[i] <= 19)
		//m[i] = 0.0;
printf("m %d %d %.20g\n", i, id[i], m[i]);
		vx[i] /= dayUnit;
		vy[i] /= dayUnit;
		vz[i] /= dayUnit;
	}

	//When Moon is not used
	//set Earth mass as Earth plus Moon,
	if(N < 11){
		for(int i = 1; i < N; ++i){
			if(id[i] == 3)
			m[3] = 1.0/pmass[2] + 1.0/pmass[19];
printf("%.20g %.20g %.20g\n", 1.0/pmass[2], 1.0/pmass[19], m[3]);


		}
	}

	
	//barycentric, all particles
	
	//convert to barycentric coordinates
	double comx = 0.0;
	double comy = 0.0;
	double comz = 0.0;
	double vcomx = 0.0;
	double vcomy = 0.0;
	double vcomz = 0.0;
	double mtot = 0.0;
	
	for(int i = 0; i < N; ++i){
		comx += m[i] * x[i];
		comy += m[i] * y[i];
		comz += m[i] * z[i];
		vcomx += m[i] * vx[i];
		vcomy += m[i] * vy[i];
		vcomz += m[i] * vz[i];
		mtot += m[i];
	}

	for(int i = 0; i < N; ++i){
		x[i] -= comx / mtot;
		y[i] -= comy / mtot;
		z[i] -= comz / mtot;
		vx[i] -= vcomx / mtot;
		vy[i] -= vcomy / mtot;
		vz[i] -= vcomz / mtot;
	}

	//first output
	
	//convert to heliocentric output
	comx = 0.0;
	comy = 0.0;
	comz = 0.0;
	vcomx = 0.0;
	vcomy = 0.0;
	vcomz = 0.0;
	
	if(outHelio == 1){
		for(int i = 1; i < N; ++i){
			comx += m[i] * x[i];
			comy += m[i] * y[i];
			comz += m[i] * z[i];
			vcomx += m[i] * vx[i];
			vcomy += m[i] * vy[i];
			vcomz += m[i] * vz[i];
		}
		comx /= m[0];
		comy /= m[0];
		comz /= m[0];
		vcomx /= m[0];
		vcomy /= m[0];
		vcomz /= m[0];
	}
	
	if(outHelio == 1){
		sprintf(outfilename, "Outhelio10_%.12d.dat", 0);
	}
	else{
		sprintf(outfilename, "Outbary10_%.12d.dat", 0);
	}
	outfile = fopen(outfilename, "w");
	printf("%s\n", outfilename);
	for(int i = 0; i < N; ++i){
		fprintf(outfile, "%.10g %d %.40g %.40g %.40g %.40g %.40g %.40g %.40g %.40g 0 0 0 0 0 0 0 0 0 0 0\n", time, i, m[i], rad[i], comx + x[i], comy + y[i], comz + z[i], vcomx + vx[i], vcomy + vy[i], vcomz + vz[i]);

	}
	fclose(outfile);

	double kx[N][4];
	double ky[N][4];
	double kz[N][4];
	double kvx[N][4];
	double kvy[N][4];
	double kvz[N][4];

	double xt[N];
	double yt[N];
	double zt[N];
	double vxt[N];
	double vyt[N];
	double vzt[N];

	//double aa[4] = {1.0, 0.0, 0.0, 0.0};	//Euler
	//double aa[4] = {0.0, 0.5, 0.0, 0.0};	//RK2
	double aa[4] = {0.0, 0.5, 0.5, 1.0};	//RK4

	for(long long int t = 1; t < Nsteps; ++t){
		double ax = 0.0;
		double ay = 0.0;
		double az = 0.0;
		
		//stage 1	
		for(int i = 0; i < N; ++i){
			xt[i] = x[i];
			yt[i] = y[i];
			zt[i] = z[i];
			vxt[i] = vx[i];
			vyt[i] = vy[i];
			vzt[i] = vz[i];
		}
		for(int i = 0; i < N; ++i){
			kx[i][0] = vxt[i];
			ky[i][0] = vyt[i];
			kz[i][0] = vzt[i];

			ax = 0.0;
			ay = 0.0;
			az = 0.0;

			for(int j = 0; j < N; ++j){
				if(i != j){
					acc(m, xt, yt, zt, ax, ay, az, i, j);
				}
			}

			kvx[i][0] = ax;
			kvy[i][0] = ay;
			kvz[i][0] = az;
		}

		
		//stage 2
		for(int i = 0; i < N; ++i){
			xt[i] = x[i] + dt * aa[1] * kx[i][0];
			yt[i] = y[i] + dt * aa[1] * ky[i][0];
			zt[i] = z[i] + dt * aa[1] * kz[i][0];
			vxt[i] = vx[i] + dt * aa[1] * kvx[i][0];
			vyt[i] = vy[i] + dt * aa[1] * kvy[i][0];
			vzt[i] = vz[i] + dt * aa[1] * kvz[i][0];
		}

		for(int i = 0; i < N; ++i){

			kx[i][1] = vxt[i];
			ky[i][1] = vyt[i];
			kz[i][1] = vzt[i];

			ax = 0.0;
			ay = 0.0;
			az = 0.0;

			for(int j = 0; j < N; ++j){
				if(i != j){
					acc(m, xt, yt, zt, ax, ay, az, i, j);
				}
			}

			kvx[i][1] = ax;
			kvy[i][1] = ay;
			kvz[i][1] = az;

		}
		//stage 3
		for(int i = 0; i < N; ++i){
			xt[i] = x[i] + dt * aa[2] * kx[i][1];
			yt[i] = y[i] + dt * aa[2] * ky[i][1];
			zt[i] = z[i] + dt * aa[2] * kz[i][1];
			vxt[i] = vx[i] + dt * aa[2] * kvx[i][1];
			vyt[i] = vy[i] + dt * aa[2] * kvy[i][1];
			vzt[i] = vz[i] + dt * aa[2] * kvz[i][1];
		}

		for(int i = 0; i < N; ++i){

			kx[i][2] = vxt[i];
			ky[i][2] = vyt[i];
			kz[i][2] = vzt[i];

			ax = 0.0;
			ay = 0.0;
			az = 0.0;

			for(int j = 0; j < N; ++j){
				if(i != j){
					acc(m, xt, yt, zt, ax, ay, az, i, j);
				}
			}

			kvx[i][2] = ax;
			kvy[i][2] = ay;
			kvz[i][2] = az;

		}
		//stage 4
		for(int i = 0; i < N; ++i){
			xt[i] = x[i] + dt * aa[3] * kx[i][2];
			yt[i] = y[i] + dt * aa[3] * ky[i][2];
			zt[i] = z[i] + dt * aa[3] * kz[i][2];
			vxt[i] = vx[i] + dt * aa[3] * kvx[i][2];
			vyt[i] = vy[i] + dt * aa[3] * kvy[i][2];
			vzt[i] = vz[i] + dt * aa[3] * kvz[i][2];
		}

		for(int i = 0; i < N; ++i){

			kx[i][3] = vxt[i];
			ky[i][3] = vyt[i];
			kz[i][3] = vzt[i];

			ax = 0.0;
			ay = 0.0;
			az = 0.0;

			for(int j = 0; j < N; ++j){
				if(i != j){
					acc(m, xt, yt, zt, ax, ay, az, i, j);
				}
			}

			kvx[i][3] = ax;
			kvy[i][3] = ay;
			kvz[i][3] = az;

		}
		
		for(int i = 0; i < N; ++i){
		
			//RK4	
			x[i] += dt * 1.0 / 6.0 * (kx[i][0] + 2.0 * kx[i][1] + 2.0 * kx[i][2] + kx[i][3]);
			y[i] += dt * 1.0 / 6.0 * (ky[i][0] + 2.0 * ky[i][1] + 2.0 * ky[i][2] + ky[i][3]);
			z[i] += dt * 1.0 / 6.0 * (kz[i][0] + 2.0 * kz[i][1] + 2.0 * kz[i][2] + kz[i][3]);

			vx[i] += dt * 1.0 / 6.0 * (kvx[i][0] + 2.0 * kvx[i][1] + 2.0 * kvx[i][2] + kvx[i][3]);
			vy[i] += dt * 1.0 / 6.0 * (kvy[i][0] + 2.0 * kvy[i][1] + 2.0 * kvy[i][2] + kvy[i][3]);
			vz[i] += dt * 1.0 / 6.0 * (kvz[i][0] + 2.0 * kvz[i][1] + 2.0 * kvz[i][2] + kvz[i][3]);
			


			/*
			//Euler 
			x[i] += dt * kx[i][0];
			y[i] += dt * ky[i][0];
			z[i] += dt * kz[i][0];

			vx[i] += dt * kvx[i][0];
			vy[i] += dt * kvy[i][0];
			vz[i] += dt * kvz[i][0];
			*/

			/*
			//RK2
			x[i] += dt * kx[i][1];
			y[i] += dt * ky[i][1];
			z[i] += dt * kz[i][1];

			vx[i] += dt * kvx[i][1];
			vy[i] += dt * kvy[i][1];
			vz[i] += dt * kvz[i][1];
			*/
		}




		if(t % outInterval == 0){
			
			if(outHelio == 1){
				sprintf(outfilename, "Outhelio10_%.12lld.dat", t);
			}
			else{
				sprintf(outfilename, "Outbary10_%.12lld.dat", t);
			}
			outfile = fopen(outfilename, "w");
			printf("%s\n", outfilename);
			
			//convert to heliocentric output
			comx = 0.0;
			comy = 0.0;
			comz = 0.0;
			vcomx = 0.0;
			vcomy = 0.0;
			vcomz = 0.0;
			
			if(outHelio == 1){
				for(int i = 1; i < N; ++i){
					comx += m[i] * x[i];
					comy += m[i] * y[i];
					comz += m[i] * z[i];
					vcomx += m[i] * vx[i];
					vcomy += m[i] * vy[i];
					vcomz += m[i] * vz[i];
				}
				comx /= m[0];
				comy /= m[0];
				comz /= m[0];
				vcomx /= m[0];
				vcomy /= m[0];
				vcomz /= m[0];
			}
			for(int i = 0; i < N; ++i){
				fprintf(outfile, "%.10g %d %.40g %.40g %.40g %.40g %.40g %.40g %.40g %.40g 0 0 0 0 0 0 0 0 0 0 0\n", time + t * dt / dayUnit, i, m[i], rad[i], comx + x[i], comy + y[i], comz + z[i], vcomx + vx[i], vcomy + vy[i], vcomz + vz[i]);
				
			}
			fclose(outfile);
		}
	}
	
}
	

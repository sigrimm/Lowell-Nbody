#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "fg.h"

#define dayUnit 0.01720209895
//#define dayUnit 1.0
//Erics value:
//#define dayUnit 0.01720412578565452474399499749324604636058

//g++ -o integrator integrator.cpp

void acc(double *m, double *x, double *y, double *z, double &ax, double &ay, double &az, int i, int j){
	double rx = x[j] - x[i];
	double ry = y[j] - y[i];
	double rz = z[j] - z[i];
	double rsq = rx * rx + ry * ry + rz * rz;
	double r = sqrt(rsq);
	
	ax += m[j] / (r * rsq) * rx;
	ay += m[j] / (r * rsq) * ry;
	az += m[j] / (r * rsq) * rz;
}
void force(double *m, double *x, double *y, double *z, double &fx, double &fy, double &fz, int i, int j){
	double rx = x[j] - x[i];
	double ry = y[j] - y[i];
	double rz = z[j] - z[i];
	double rsq = rx * rx + ry * ry + rz * rz;
	double r = sqrt(rsq);
	
	fx += m[i] * m[j] / (r * rsq) * rx;
	fy += m[i] * m[j] / (r * rsq) * ry;
	fz += m[i] * m[j] / (r * rsq) * rz;
}


void Drift(double *m, double *x, double *y, double *z, double *vx, double *vy, double *vz, double dt, int i, double mu, int GR){
	
	if(GR == 1){// GR time rescale (Saha & Tremaine 1994)
		double c = 10065.3201686;//c in AU / day * dayUnit
		double c2 = c * c;
		
		double rsq = x[i]*x[i] + y[i]*y[i] + z[i]*z[i];
		double vsq = vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
		double ir = 1.0 / sqrt(rsq);
		double ia = 2.0 * ir  -vsq  / mu;
		dt *= 1.0 - 1.5 * mu * ia / c2;
        }
	x[i] += dt * vx[i];
	y[i] += dt * vy[i];
	z[i] += dt * vz[i];
}

//in heliocentric coordinates
void GRh(double Msun, double m, double rx, double ry, double rz, double vx, double vy, double vz, double &ax, double &ay, double &az){
	
	double c = 10065.3201686;//c in AU / day * dayUnit
	double csq = c * c;
	
	double rsq = rx * rx + ry * ry + rz * rz;
	double r = sqrt(rsq);
	
        double vsq = vx * vx + vy * vy + vz * vz;
        double rd = (rx * vx + ry * vy + rz * vz) / r;

        double A = (Msun + m) / r;
        double B = A / (r * csq);

        double eta = Msun * m / ((Msun + m) * (Msun + m));

        double C = 2.0 * (2.0 - eta) * rd;
        double D = (1.0 + 3.0 * eta) * vsq - 1.5 * eta * rd * rd - 2.0 * (2.0 + eta) * A;

        ax += B * (C * vx - D * rx / r);
        ay += B * (C * vy - D * ry / r);
        az += B * (C * vz - D * rz / r);
}
void GRhP(double Msun, double m, double x, double y, double z, double &ax, double &ay, double &az){
	// GR symplectic
	// GR part depending on position only (see Saha & Tremaine 1994)
	double c = 10065.3201686;//c in AU / day * dayUnit
	double rsq = x * x + y * y + z * z;

	double mu = (Msun + m);
	double A = mu / (rsq * c);
	double B = 2.0 * A * A;
	ax -= B * x;
	ay -= B * y;
	az -= B * z;
}
void GRdP(double Msun, double m, double x, double y, double z, double &ax, double &ay, double &az){
	// GR symplectic
	// GR part depending on position only (see Saha & Tremaine 1994)
	double c = 10065.3201686;//c in AU / day * dayUnit
	double rsq = x * x + y * y + z * z;

	double mu = 1.0;//(Msun + m);
	double A = mu / (rsq * c);
	double B = 2.0 * A * A;
	ax -= B * x;
	ay -= B * y;
	az -= B * z;
}
void GRV(double *x, double *y, double *z, double *vx, double *vy, double *vz, double dt, int i){
	// GR part depending on velocity only (see Saha & Tremaine 1994)
	double c = 10065.3201686;//c in AU / day * dayUnit
	double c2 = c * c;
	double vsq = vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
	double vcdt = 2.0 * vsq/c2 * dt;
	x[i] -= vx[i] * vcdt;
	y[i] -= vy[i] * vcdt;
	z[i] -= vz[i] * vcdt;
}

//implicit midpoint
void KickhGRI(double *m, double *x, double *y, double *z, double *vx, double *vy, double *vz, double dt, int i){
	double ax = 0.0;
	double ay = 0.0;
	double az = 0.0;
	
	double c = 10065.3201686;//c in AU / day * dayUnit
	double csq = c * c;
	
	
	double rx = -x[i];
	double ry = -y[i];
	double rz = -z[i];
	double rsq = rx * rx + ry * ry + rz * rz;
	double ir = 1.0 / sqrt(rsq);

	double A = (m[0] + m[i]) * ir;
	double B = A * ir / csq;
	double eta = m[0] * m[i] / ((m[0] + m[i]) * (m[0] + m[i]));

	double vtx = vx[i];
	double vty = vy[i];
	double vtz = vz[i];

	//printf("%.40g %.40g %.40g\n", ax, ay, az);
	for(int k = 0; k < 30; ++k){
		ax = 0.0;
		ay = 0.0;
		az = 0.0;
		double vsq = (vtx * vtx + vty * vty + vtz * vtz);
		double rd = (x[i] * vtx + y[i] * vty + z[i] * vtz) * ir;

		double C = 2.0 * (2.0 - eta) * rd;
		double D = (1.0 + 3.0 * eta) * vsq - 1.5 * eta * rd * rd - 2.0 * (2.0 + eta) * A;

		ax += B * (C * vtx - D * x[i] * ir);
		ay += B * (C * vty - D * y[i] * ir);
		az += B * (C * vtz - D * z[i] * ir);
		
		vtx = vx[i] + 0.5 * dt * ax;
		vty = vy[i] + 0.5 * dt * ay;
		vtz = vz[i] + 0.5 * dt * az;


	//printf("%.40g %.40g %.40g\n", ax, ay, az);
	}

	vx[i] += dt * ax;
	vy[i] += dt * ay;
	vz[i] += dt * az;
}


void Kick(double *m, double *x, double *y, double *z, double *vx, double *vy, double *vz, double dt, int i, int N){
	double ax = 0.0;
	double ay = 0.0;
	double az = 0.0;

	for(int j = 1; j < N; ++j){
		if(i != j){
			acc(m, x, y, z, ax, ay, az, i, j);
		}
	}
	
	vx[i] += dt * ax;
	vy[i] += dt * ay;
	vz[i] += dt * az;	
}

//including the Sun
void KickAll(double *m, double *x, double *y, double *z, double *vx, double *vy, double *vz, double dt, int i, int N){
	double ax = 0.0;
	double ay = 0.0;
	double az = 0.0;

	for(int j = 0; j < N; ++j){
		if(i != j){
			acc(m, x, y, z, ax, ay, az, i, j);
		}
	}
	
	vx[i] += dt * ax;
	vy[i] += dt * ay;
	vz[i] += dt * az;	
}

//in heliocentric coordinates
//aj = Fj / mj + Sum,i=1toN Fi/mstar 
void KickAllh(double *m, double *x, double *y, double *z, double *vx, double *vy, double *vz, double dt, int N){
	double fx[N];
	double fy[N];
	double fz[N];

	for(int i = 0; i < N; ++i){
		fx[i] = 0.0;
		fy[i] = 0.0;
		fz[i] = 0.0;
	}
	
	for(int i = 1; i < N; ++i){
		for(int j = 0; j < N; ++j){
			if(i != j){
				force(m, x, y, z, fx[i], fy[i], fz[i], i, j);
			}
		}
	}

	for(int i = 1; i < N; ++i){
		double ax = fx[i] / m[i];
		double ay = fy[i] / m[i];
		double az = fz[i] / m[i];
		
		for(int j = 0; j < N; ++j){
			ax += fx[j] / m[0];
			ay += fy[j] / m[0];
			az += fz[j] / m[0];
		
		}
		
		vx[i] += dt * ax;
		vy[i] += dt * ay;
		vz[i] += dt * az;	

	}
}

//in heliocentric coordinates
//aj = Fj / mj + Sum,i=1toN Fi/mstar 
void Kickh(double *m, double *x, double *y, double *z, double *vx, double *vy, double *vz, double dt, int N){
	double fx[N];
	double fy[N];
	double fz[N];

	for(int i = 0; i < N; ++i){
		fx[i] = 0.0;
		fy[i] = 0.0;
		fz[i] = 0.0;
	}

	for(int i = 1; i < N; ++i){
		for(int j = 1; j < N; ++j){
			if(i != j){
				force(m, x, y, z, fx[i], fy[i], fz[i], i, j);
			}
		}
	}

	for(int i = 1; i < N; ++i){
		double ax = fx[i] / m[i];
		double ay = fy[i] / m[i];
		double az = fz[i] / m[i];

		for(int j = 1; j < N; ++j){
			ax += fx[j] / m[0];
			ay += fy[j] / m[0];
			az += fz[j] / m[0];

		}

		vx[i] += dt * ax;
		vy[i] += dt * ay;
		vz[i] += dt * az;

	}
}



void KickGRh(double *m, double *x, double *y, double *z, double *vx, double *vy, double *vz, double dt, int i, int N){
	double ax = 0.0;
	double ay = 0.0;
	double az = 0.0;

	GRhP(m[0], m[i], x[i], y[i], z[i], ax, ay, az);
	
	vx[i] += dt * ax;
	vy[i] += dt * ay;
	vz[i] += dt * az;	
}

void KickGRd(double *m, double *x, double *y, double *z, double *vx, double *vy, double *vz, double dt, int i, int N){
	double ax = 0.0;
	double ay = 0.0;
	double az = 0.0;

	GRdP(m[0], m[i], x[i], y[i], z[i], ax, ay, az);
	
	vx[i] += dt * ax;
	vy[i] += dt * ay;
	vz[i] += dt * az;	
}

//democratic
void SunKick(double *m, double *x, double *y, double *z, double *vx, double *vy, double *vz, double dt, int N){
	
	double ax = 0.0;
	double ay = 0.0;
	double az = 0.0;
	for(int i = 1; i < N; ++i){
		ax += m[i] * vx[i];
		ay += m[i] * vy[i];
		az += m[i] * vz[i];
	}
	for(int i = 1; i < N; ++i){
		x[i] += ax * dt / m[0];
		y[i] += ay * dt / m[0];
		z[i] += az * dt / m[0];
	}
}

//f = 1  convert bary to helio
//f = -1 convert helio to bary
void convert(double *m, double *vx, double *vy, double *vz, int N, int f){
	
	double px = 0.0;
	double py = 0.0;
	double pz = 0.0;
	double Mtot = 0.0;
	for(int i = 1; i < N; ++i){
		px += m[i] * vx[i];
		py += m[i] * vy[i];
		pz += m[i] * vz[i];
		Mtot += m[i];
	}
	
	for(int i = 1; i < N; ++i){
		if(f == 1){
			double iMsun = 1.0 / m[0];
			vx[i] += px * iMsun;
			vy[i] += py * iMsun;
			vz[i] += pz * iMsun;
		}
		if(f == -1){
			double iMsunp = 1.0 / (m[0] + Mtot);
			vx[i] -= px * iMsunp;
			vy[i] -= py * iMsunp;
			vz[i] -= pz * iMsunp;
		}
	}
}

int main(){

	//Number of planets
	
	int integrator = 30;
	// 10 barycentric integrating the sun
	// 20 heliocentric, not integrating the sun
	// 30 democratic
	
	const int NN = 21;
	int GR = 1;
	// 0 no GR
	// 2 implicit
	// 1 symplectic
	
	int outHelio = 1;
	//1 print output in heliocentric coordinates
	//0 print output in barycentric  coordinates
	int useFG = 1;

	//integrator 20, useFG = 0, GR 0, Nsteps 4000000, dt 0.001 * dayUnit, OK
	//integrator 20, useFG = 0, GR 2, Nsteps 4000000, dt 0.001 * dayUnit, OK
		//integrator 20, useFG = 0, GR 1, Nsteps 4000000, dt 0.001 * dayUnit, fail

	//integrator 30, useFG = 0, GR 0, Nsteps 4000000, dt 0.001 * dayUnit, OK
	//integrator 30, useFG = 0, GR 2, Nsteps 4000000, dt 0.001 * dayUnit, OK
		//integrator 30, useFG = 0, GR 1, Nsteps 4000000, dt 0.001 * dayUnit, fail

	//integrator 30, useFG = 1, GR 0, Nsteps 400000, dt 0.01 * dayUnit, OK
	//integrator 30, useFG = 1, GR 2, Nsteps 400000, dt 0.01 * dayUnit, OK

	long long int Nsteps = 400000;	
	long long int outInterval = 100;
	double dt = 0.01 * dayUnit;
	
	int id[NN];
	double m[NN];
	double rad[NN];
	double x[NN];
	double y[NN];
	double z[NN];
	double vx[NN];
	double vy[NN];
	double vz[NN];

	double a[NN];
	double e[NN];
	double inc[NN];
	double Omega[NN];
	double w[NN];
	double M[NN];
	double l[NN];

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
printf("Earth Moon %.20g %.20g %.20g\n", 1.0/pmass[2], 1.0/pmass[19], m[3]);


		}
	}

	
	if(integrator == 10){
	//barycentric, all particles
		
		//convert to barycentric positions
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
			mtot += m[i];
		}
		
		for(int i = 0; i < N; ++i){
			x[i] -= comx / mtot;
			y[i] -= comy / mtot;
			z[i] -= comz / mtot;
		}
		//convert to barycentric velocities
		
		for(int i = 0; i < N; ++i){
			vcomx += m[i] * vx[i];
			vcomy += m[i] * vy[i];
			vcomz += m[i] * vz[i];
		}
		
		for(int i = 0; i < N; ++i){
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
			}
			comx /= m[0];
			comy /= m[0];
			comz /= m[0];

			
			for(int i = 1; i < N; ++i){
				vcomx += m[i] * vx[i];
				vcomy += m[i] * vy[i];
				vcomz += m[i] * vz[i];
			}
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
		//printf("%s\n", outfilename);
		for(int i = 0; i < N; ++i){
			fprintf(outfile, "%.10g %d %.40g %.40g %.40g %.40g %.40g %.40g %.40g %.40g 0 0 0 0 0 0 0 0 0 0 0\n", time, i-1, m[i], rad[i], comx + x[i], comy + y[i], comz + z[i], vcomx + vx[i], vcomy + vy[i], vcomz + vz[i]);

		}
		fclose(outfile);
		for(long long int t = 1; t < Nsteps; ++t){
			
			if(useFG == 0){
				for(int i = 0; i < N; ++i){
			//			KickhGR(m, x, y, z, vx, vy, vz, dt, i, N);
					KickAll(m, x, y, z, vx, vy, vz, 0.5 * dt, i, N);
				}
			}
			else{
				for(int i = 1; i < N; ++i){
					Kick(m, x, y, z, vx, vy, vz, 0.5 * dt, i, N);
				}
			}
			
			for(int i = 0; i < N; ++i){
				if(useFG == 0){
					Drift(m, x, y, z, vx, vy, vz, dt, i, m[0], GR);
				}
				else{
					//this is not working correctly
					fgfull(x[i], y[i], z[i], vx[i], vy[i], vz[i], dt, m[0], GR);
				}
			}
			
			if(useFG == 0){
				for(int i = 0; i < N; ++i){
			//			KickhGR(m, x, y, z, vx, vy, vz, dt, i, N);
					KickAll(m, x, y, z, vx, vy, vz, 0.5 * dt, i, N);
				}
			}
			else{
				for(int i = 1; i < N; ++i){
					Kick(m, x, y, z, vx, vy, vz, 0.5 * dt, i, N);
				}
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
					}
					comx /= m[0];
					comy /= m[0];
					comz /= m[0];

					
					for(int i = 1; i < N; ++i){
						vcomx += m[i] * vx[i];
						vcomy += m[i] * vy[i];
						vcomz += m[i] * vz[i];
					}
					vcomx /= m[0];
					vcomy /= m[0];
					vcomz /= m[0];
				}
				for(int i = 0; i < N; ++i){
					fprintf(outfile, "%.10g %d %.40g %.40g %.40g %.40g %.40g %.40g %.40g %.40g 0 0 0 0 0 0 0 0 0 0 0\n", time + t * dt / dayUnit, i-1, m[i], rad[i], comx + x[i], comy + y[i], comz + z[i], vcomx + vx[i], vcomy + vy[i], vcomz + vz[i]);
					
				}
				fclose(outfile);
			}
		}
		
	}
	

	if(integrator == 20){
	//heliocentric, not integrating the sun
		//first output
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

		if(outHelio == 1){
			sprintf(outfilename, "Outhelio20_%.12d.dat", 0);
		}
		else{
			sprintf(outfilename, "Outbary20_%.12d.dat", 0);
		}
		outfile = fopen(outfilename, "w");
		printf("%s\n", outfilename);
		for(int i = 1; i < N; ++i){
			double xp, yp, zp, vxp, vyp, vzp;
			if(outHelio == 1){
				xp = x[i];
				yp = y[i];
				zp = z[i];
				vxp = x[i];
				vyp = y[i];
				vzp = z[i];
			}
			else{
				xp = x[i] - comx / mtot;
				yp = y[i] - comy / mtot;
				zp = z[i] - comz / mtot;
				vxp = x[i] - comx / mtot;
				vyp = y[i] - comy / mtot;
				vzp = z[i] - comz / mtot;
			}
			fprintf(outfile, "%.10g %d %.40g %.40g %.40g %.40g %.40g %.40g %.40g %.40g 0 0 0 0 0 0 0 0 0 0 0\n", time, i-1, m[i], rad[i], xp, yp, zp, vxp, vyp, vzp);

		}
		fclose(outfile);
	
		for(long long int t = 1; t < Nsteps; ++t){
			//time += dt / dayUnit;
			
			if(GR == 1){
				for(int i = 1; i < N; ++i){
					KickGRh(m, x, y, z, vx, vy, vz, 0.5 * dt, i, N);
				}
			}
			
			if(GR == 2){
				for(int i = 1; i < N; ++i){
					KickhGRI(m, x, y, z, vx, vy, vz, 0.5 * dt, i);
				}
			}
			if(useFG == 1){
				Kickh(m, x, y, z, vx, vy, vz, 0.5 * dt, N);
			}
			else{
				KickAllh(m, x, y, z, vx, vy, vz, 0.5 * dt, N);
				//Gaskickh(m[i], rad[i], x[i], y[i], z[i], vx[i], vy[i], vz[i], 0.5 * dt, m[0], time, i, N);
			}
			if(GR == 1){
				for(int i = 1; i < N; ++i){
					GRV(x, y, z, vx, vy, vz, 0.5 * dt, i);
				}
			}
			for(int i = 1; i < N; ++i){
				//Solve the Kepler equation
				if(useFG == 1){
					//this is not yet working correctly
					fgfull(x[i], y[i], z[i], vx[i], vy[i], vz[i], dt, m[0] + m[i], GR);
				}
				else{
				//no Kepler equation
					Drift(m, x, y, z, vx, vy, vz, dt, i, m[0] + m[i], GR);
				}
			}
			if(GR == 1){
				for(int i = 1; i < N; ++i){
					GRV(x, y, z, vx, vy, vz, 0.5 * dt, i);
				}
			}
			if(useFG == 1){
				Kickh(m, x, y, z, vx, vy, vz, 0.5 * dt, N);
			}
			else{
				KickAllh(m, x, y, z, vx, vy, vz, 0.5 * dt, N);
				//Gaskickh(m[i], rad[i], x[i], y[i], z[i], vx[i], vy[i], vz[i], 0.5 * dt, m[0], time, i, N);
			}
			if(GR == 1){
				for(int i = 1; i < N; ++i){
					KickGRh(m, x, y, z, vx, vy, vz, 0.5 * dt, i, N);
				}
			}
			if(GR == 2){
				for(int i = 1; i < N; ++i){
					KickhGRI(m, x, y, z, vx, vy, vz, 0.5 * dt, i);
				}
			}
			if(t % outInterval == 0){
				double comx = 0.0;
				double comy = 0.0;
				double comz = 0.0;
				double vcomx = 0.0;
				double vcomy = 0.0;
				double vcomz = 0.0;
				double mtot = 0.0;
				//convert to barycentric velocities

				for(int i = 0; i < N; ++i){
					comx += m[i] * x[i];
					comy += m[i] * y[i];
					comz += m[i] * z[i];
					vcomx += m[i] * vx[i];
					vcomy += m[i] * vy[i];
					vcomz += m[i] * vz[i];
					mtot += m[i];
				}
				if(outHelio == 1){
					sprintf(outfilename, "Outhelio20_%.12lld.dat", t);
				}
				else{
					sprintf(outfilename, "Outbary20_%.12lld.dat", t);
				}
				outfile = fopen(outfilename, "w");
				//printf("%s\n", outfilename);

				for(int i = 1; i < N; ++i){
					double xp, yp, zp, vxp, vyp, vzp;
					if(outHelio == 1){
						xp = x[i];
						yp = y[i];
						zp = z[i];
						vxp = x[i];
						vyp = y[i];
						vzp = z[i];
					}
					else{
						xp = x[i] - comx / mtot;
						yp = y[i] - comy / mtot;
						zp = z[i] - comz / mtot;
						vxp = x[i] - comx / mtot;
						vyp = y[i] - comy / mtot;
						vzp = z[i] - comz / mtot;
					}

					fprintf(outfile, "%.10g %d %.20g %.20g %.20g %.20g %.20g %.20g %.20g %.20g 0 0 0 0 0 0 0 0 0 0 0\n", time + t * dt / dayUnit, i, m[i], rad[i], xp, yp, zp, vxp, vyp, vzp);
				}
				fclose(outfile);		
                        }
		
		}
	}



	if(integrator == 30){
		//democratic
		double comx = 0.0;
		double comy = 0.0;
		double comz = 0.0;
		double vcomx = 0.0;
		double vcomy = 0.0;
		double vcomz = 0.0;
		double mtot = 0.0;
		//convert to barycentric velocities
		
		for(int i = 0; i < N; ++i){
			vcomx += m[i] * vx[i];
			vcomy += m[i] * vy[i];
			vcomz += m[i] * vz[i];
			mtot += m[i];
		}
		
		for(int i = 0; i < N; ++i){
			vx[i] -= vcomx / mtot;
			vy[i] -= vcomy / mtot;
			vz[i] -= vcomz / mtot;
		}

		//first output


		//convert to heliocentric or barycentric coordinates for output
                comx = 0.0;
                comy = 0.0;
                comz = 0.0;
		vcomx = 0.0;
		vcomy = 0.0;
		vcomz = 0.0;

		for(int i = 1; i < N; ++i){
			comx += m[i] * x[i];
			comy += m[i] * y[i];
			comz += m[i] * z[i];
			vcomx += m[i] * vx[i];
			vcomy += m[i] * vy[i];
			vcomz += m[i] * vz[i];
		}

		if(outHelio == 1){
			sprintf(outfilename, "Outhelio30_%.12d.dat", 0);
		}
		else{
			sprintf(outfilename, "Outbary30_%.12d.dat", 0);
		}
		outfile = fopen(outfilename, "w");
		printf("%s\n", outfilename);
		for(int i = 1; i < N; ++i){
			double xp, yp, zp, vxp, vyp, vzp;
			if(outHelio == 1){
				xp = x[i];
				yp = y[i];
				zp = z[i];
				vxp = vx[i] + vcomx / m[0];
				vyp = vy[i] + vcomy / m[0];
				vzp = vz[i] + vcomz / m[0];
			}
			else{
				xp = x[i] - comx / mtot;
				yp = y[i] - comy / mtot;
				zp = z[i] - comz / mtot;
				vxp = vx[i];
				vyp = vy[i];
				vzp = vz[i];
			}
			fprintf(outfile, "%.10g %d %.40g %.40g %.40g %.40g %.40g %.40g %.40g %.40g 0 0 0 0 0 0 0 0 0 0 0\n", time, i-1, m[i], rad[i], xp, yp, zp, vxp, vyp, vzp);

		}
		fclose(outfile);

		for(long long int t = 1; t < Nsteps; ++t){

			if(GR == 1){
				for(int i = 1; i < N; ++i){
					//KickGRd(m, x, y, z, vx, vy, vz, 0.5 * dt, i, N);
				}
			}
			if(GR == 2){
				convert(m, vx, vy, vz, N, 1);
				for(int i = 1; i < N; ++i){
					KickhGRI(m, x, y, z, vx, vy, vz, 0.5 * dt, i);
				}
				convert(m, vx, vy, vz, N, -1);
			}
			
			
			for(int i = 1; i < N; ++i){
				if(useFG == 1){
					Kick(m, x, y, z, vx, vy, vz, 0.5 * dt, i, N);
				}
				else{
					KickAll(m, x, y, z, vx, vy, vz, 0.5 * dt, i, N);
				}
			}
			
			SunKick(m, x, y, z, vx, vy, vz, 0.5 * dt, N);
			if(GR == 1){
				convert(m, vx, vy, vz, N, 1);
				for(int i = 1; i < N; ++i){
					//GRV(x, y, z, vx, vy, vz, 0.5 * dt, i);
				}
				convert(m, vx, vy, vz, N, -1);
			}		
			
			for(int i = 1; i < N; ++i){
				if(useFG == 0){
					Drift(m, x, y, z, vx, vy, vz, dt, i, m[0], GR);
				}
				else{
					fgfull(x[i], y[i], z[i], vx[i], vy[i], vz[i], dt, m[0], GR);
				}
			}
			SunKick(m, x, y, z, vx, vy, vz, 0.5 * dt, N);
			if(GR == 1){
				convert(m, vx, vy, vz, N, 1);
				for(int i = 1; i < N; ++i){
					//GRV(x, y, z, vx, vy, vz, 0.5 * dt, i);
				}
				convert(m, vx, vy, vz, N, -1);
			}
			for(int i = 1; i < N; ++i){
				if(useFG == 1){
					Kick(m, x, y, z, vx, vy, vz, 0.5 * dt, i, N);
				}
				else{
					KickAll(m, x, y, z, vx, vy, vz, 0.5 * dt, i, N);
				}
			}
			if(GR == 1){
				for(int i = 1; i < N; ++i){
					//KickGRd(m, x, y, z, vx, vy, vz, 0.5 * dt, i, N);
				}
			}
			if(GR == 2){
				convert(m, vx, vy, vz, N, 1);
				for(int i = 1; i < N; ++i){
					KickhGRI(m, x, y, z, vx, vy, vz, 0.5 * dt, i);
				}
				convert(m, vx, vy, vz, N, -1);
			}
			if(t % outInterval == 0){
				double comx = 0.0;
				double comy = 0.0;
				double comz = 0.0;
				double vcomx = 0.0;
				double vcomy = 0.0;
				double vcomz = 0.0;

				for(int i = 1; i < N; ++i){
					comx += m[i] * x[i];
					comy += m[i] * y[i];
					comz += m[i] * z[i];
					vcomx += m[i] * vx[i];
					vcomy += m[i] * vy[i];
					vcomz += m[i] * vz[i];
				}
				if(outHelio == 1){
					sprintf(outfilename, "Outhelio30_%.12lld.dat", t);
				}
				else{
					sprintf(outfilename, "Outbary30_%.12lld.dat", t);
				}
				outfile = fopen(outfilename, "w");
				printf("%s\n", outfilename);

				for(int i = 1; i < N; ++i){
					double xp, yp, zp, vxp, vyp, vzp;
					if(outHelio == 1){
						xp = x[i];
						yp = y[i];
						zp = z[i];
						vxp = vx[i] + vcomx / m[0];
						vyp = vy[i] + vcomy / m[0];
						vzp = vz[i] + vcomz / m[0];
					}
					else{
						xp = x[i] - comx / mtot;
						yp = y[i] - comy / mtot;
						zp = z[i] - comz / mtot;
						vxp = vx[i];
						vyp = vy[i];
						vzp = vz[i];
					}
					fprintf(outfile, "%.10g %d %.40g %.40g %.40g %.40g %.40g %.40g %.40g %.40g 0 0 0 0 0 0 0 0 0 0 0\n", time + t * dt / dayUnit, i-1, m[i], rad[i], xp, yp, zp, vxp, vyp, vzp);

				}
				fclose(outfile);
			}
		}
	}


}

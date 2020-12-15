#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define dayUnit 0.01720209895
//#define dayUnit 0.01720209894846
//#define dayUnit 1.0
//Erics value:
//#define dayUnit 0.01720412578565452474399499749324604636058

// --------------------------------
//barycentric coordinates
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
// ---------------------------------

// --------------------------------
//heliocentric coordinates
//sun part
void accS(double *m, double *x, double *y, double *z, double &ax, double &ay, double &az, int i){

	double rx = -x[i];
	double ry = -y[i];
	double rz = -z[i];
	double rsq = rx * rx + ry * ry + rz * rz;
	double r = sqrt(rsq);
	
	ax += (m[0] + m[i]) / (r * rsq) * rx;
	ay += (m[0] + m[i]) / (r * rsq) * ry;
	az += (m[0] + m[i]) / (r * rsq) * rz;
}
//planet part
void accP(double *m, double *x, double *y, double *z, double &ax, double &ay, double &az, int i, int j){

	double rx = x[j] - x[i];
	double ry = y[j] - y[i];
	double rz = z[j] - z[i];
	double rsq = rx * rx + ry * ry + rz * rz;
	double r = sqrt(rsq);
	
	ax += m[j] / (r * rsq) * rx;
	ay += m[j] / (r * rsq) * ry;
	az += m[j] / (r * rsq) * rz;
}
//planet part 2
void accP2(double *m, double *x, double *y, double *z, double &ax, double &ay, double &az, int j){

	double rx = -x[j];
	double ry = -y[j];
	double rz = -z[j];
	double rsq = rx * rx + ry * ry + rz * rz;
	double r = sqrt(rsq);
	
	ax += m[j] / (r * rsq) * rx;
	ay += m[j] / (r * rsq) * ry;
	az += m[j] / (r * rsq) * rz;
}

//GR Quinn, Tremaine, Duncan 1991, heliocentric
void accGR(double *m, double *x, double *y, double *z, double *vx, double *vy, double *vz, double &ax, double &ay, double &az, int i){
	
	double c2 = 10065.3201686 * 10065.3201686;

	double rsq = x[i] * x[i] + y[i] * y[i] + z[i] * z[i];
	double r = sqrt(rsq);
	double vsq = vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];

	double rv = x[i] * vx[i] + y[i] * vy[i] + z[i] * vz[i];

	double t1 = 2.0 * (1.0 + 1.0) * m[0] / (c2 * r) - 1.0 * vsq / c2;
	double t2 = (2.0 + 2.0 * 1.0) * m[0] / (c2 * r * rsq) * rv;

	ax += m[0] / (r * rsq) * x[i] * t1 + t2 * vx[i];
	ay += m[0] / (r * rsq) * y[i] * t1 + t2 * vy[i];
	az += m[0] / (r * rsq) * z[i] * t1 + t2 * vz[i];

}
//Sitarski 1982, Isotropic equation 5, heliocentric
//modified k2 to dayUnit
//should be equivalent to the Quin et all function, assuming m[0] = 1.0
void accGR2(double *m, double *x, double *y, double *z, double *vx, double *vy, double *vz, double &ax, double &ay, double &az, int i){
	
	double c2 = 10065.3201686 * 10065.3201686;

	double rsq = x[i] * x[i] + y[i] * y[i] + z[i] * z[i];
	double r = sqrt(rsq);
	double vsq = vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];

	double rv = x[i] * vx[i] + y[i] * vy[i] + z[i] * vz[i];

	double f1 = 1.0 / (r * rsq * c2);
	double t1 = 4.0 / r;
	double t2 = -vsq;
	double t3 = 4.0 * rv;


	ax += f1 * ((t1 + t2) * x[i] + t3 * vx[i]);
	ay += f1 * ((t1 + t2) * y[i] + t3 * vy[i]);
	az += f1 * ((t1 + t2) * z[i] + t3 * vz[i]);

}
//Fabricky 2010, heliocentric
void accGR3(double *m, double *x, double *y, double *z, double *vx, double *vy, double *vz, double &ax, double &ay, double &az, int i){
	
	double c2 = 10065.3201686 * 10065.3201686;

	double rsq = x[i] * x[i] + y[i] * y[i] + z[i] * z[i];
	double r = sqrt(rsq);

	double A = (m[0] + m[i]) / r;
	double B = A / ( r * c2);
	double eta = m[0] * m[i] / ((m[0] + m[i]) * (m[0] + m[i]));
	double vsq = (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
	double rd = (x[i] * vx[i] + y[i] * vy[i] + z[i] * vz[i]) / r; 

	double C = 2.0 * (2.0 - eta) * rd;
	double D = (1.0 + 3.0 * eta) * vsq - 1.5 * eta * rd * rd - 2.0 * (2.0 + eta) * A;
	ax += B * (C * vx[i] - D * x[i] / r);         
	ay += B * (C * vy[i] - D * y[i] / r);
	az += B * (C * vz[i] - D * z[i] / r);

}
// --------------------------------


int main(){

	//Number of planets
	const int NN = 22;

	int GR = 2;
	//1 Quin Tremaine Duncan 1991, heliocentric coordinates
	//2 Sitarski 1982, heliocentric coordinates
	//3 Fabricky 2010, heliocentric coordinates

	int useHelio = 1;
	int outHelio = 1;
	//1 print output in heliocentric coordinates
	//0 print output in barycentric  coordinates

	//long long int Nsteps = 40000;	
	//long long int outInterval = 10;
	//double dt = 0.1 * dayUnit;

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
	    1.000000000000000e0,      // Sun        (21)
	    0.0			      // Test particle        (22)
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

		if(i >=21){
			m[i] = 0.0;
		}
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
	
	double comx = 0.0;
	double comy = 0.0;
	double comz = 0.0;
	double vcomx = 0.0;
	double vcomy = 0.0;
	double vcomz = 0.0;
	double mtot = 0.0;

	if(useHelio == 0){	
		//convert to barycentric coordinates
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
	}

	//first output
	
	comx = 0.0;
	comy = 0.0;
	comz = 0.0;
	vcomx = 0.0;
	vcomy = 0.0;
	vcomz = 0.0;
	if(useHelio == 0 && outHelio == 1){	
		//convert to heliocentric output
	
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
	if(useHelio == 1 && outHelio == 0){	
		//convert to barycentric output
		// ************** //
printf("Error, funtion not supported\n");
		return 0;
		// add function here
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

	double kx[N][6];
	double ky[N][6];
	double kz[N][6];
	double kvx[N][6];
	double kvy[N][6];
	double kvz[N][6];

	double xt[N];
	double yt[N];
	double zt[N];
	double vxt[N];
	double vyt[N];
	double vzt[N];

	double errorkx[N];
	double errorky[N];
	double errorkz[N];
	double errorkvx[N];
	double errorkvy[N];
	double errorkvz[N];

	double a21 = 1.0/4.0;

	double a31 = 3.0/32.0;
	double a32 = 9.0/32.0;

	double a41 = 1932.0/2197.0;
	double a42 = -7200.0/2197.0;
	double a43 = 7296.0/2197.0;

	double a51 = 439.0/216.0;
	double a52 = -8.0;
	double a53 = 3680.0/513.0;
	double a54 = -845.0/4104.0;

	double a61 = -8.0/27.0;
	double a62 = 2.0;
	double a63 = -3544/2565.0;
	double a64 = 1859.0/4104.0;
	double a65 = -11.0/40.0;

	double b1 = 25.0/216.0;
	double b3 = 1408.0/2565.0;
	double b4 = 2197.0/4104.0;
	double b5 = -1.0/5.0;
	double b6 = 0.0;

	double bb1 = 16.0/135.0;
	double bb3 = 6656.0/12825.0;
	double bb4 = 28561.0/56430.0;
	double bb5 = -9.0/50.0;
	double bb6 = 2.0/55.0;

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

			if(useHelio == 0){
				for(int j = 0; j < N; ++j){
					if(i != j){
						acc(m, xt, yt, zt, ax, ay, az, i, j);
					}
				}
			}
			else{
				
				if(i > 0){
					accS(m, xt, yt, zt, ax, ay, az, i);
				}
				for(int j = 1; j < N; ++j){
					if(i != j){
						accP(m, xt, yt, zt, ax, ay, az, i, j);
					}
				}
				for(int j = 1; j < N; ++j){
					if(i != j){
						accP2(m, xt, yt, zt, ax, ay, az, j);
					}
				}
				if(GR == 1 && i > 0){
					accGR(m, xt, yt, zt, vxt, vyt, vzt, ax, ay, az, i);
				}
				if(GR == 2 && i > 0){
					accGR2(m, xt, yt, zt, vxt, vyt, vzt, ax, ay, az, i);
				}
				if(GR == 3 && i > 0){
					accGR3(m, xt, yt, zt, vxt, vyt, vzt, ax, ay, az, i);
				}
			}
			kvx[i][0] = ax;
			kvy[i][0] = ay;
			kvz[i][0] = az;
		}

		
		//stage 2
		for(int i = 0; i < N; ++i){
			xt[i]  = x[i]  + dt * a21 * kx[i][0];
			yt[i]  = y[i]  + dt * a21 * ky[i][0];
			zt[i]  = z[i]  + dt * a21 * kz[i][0];
			vxt[i] = vx[i] + dt * a21 * kvx[i][0];
			vyt[i] = vy[i] + dt * a21 * kvy[i][0];
			vzt[i] = vz[i] + dt * a21 * kvz[i][0];
		}

		for(int i = 0; i < N; ++i){

			kx[i][1] = vxt[i];
			ky[i][1] = vyt[i];
			kz[i][1] = vzt[i];

			ax = 0.0;
			ay = 0.0;
			az = 0.0;

			if(useHelio == 0){
				for(int j = 0; j < N; ++j){
					if(i != j){
						acc(m, xt, yt, zt, ax, ay, az, i, j);
					}
				}
			}
			else{
				
				if(i > 0){
					accS(m, xt, yt, zt, ax, ay, az, i);
				}
				for(int j = 1; j < N; ++j){
					if(i != j){
						accP(m, xt, yt, zt, ax, ay, az, i, j);
					}
				}
				for(int j = 1; j < N; ++j){
					if(i != j){
						accP2(m, xt, yt, zt, ax, ay, az, j);
					}
				}
				if(GR == 1 && i > 0){
					accGR(m, xt, yt, zt, vxt, vyt, vzt, ax, ay, az, i);
				}
				if(GR == 2 && i > 0){
					accGR2(m, xt, yt, zt, vxt, vyt, vzt, ax, ay, az, i);
				}
				if(GR == 3 && i > 0){
					accGR3(m, xt, yt, zt, vxt, vyt, vzt, ax, ay, az, i);
				}
			}

			kvx[i][1] = ax;
			kvy[i][1] = ay;
			kvz[i][1] = az;

		}
		//stage 3
		for(int i = 0; i < N; ++i){
			xt[i]  = x[i]  + dt * (a31 * kx[i][0]  + a32 * kx[i][1]);
			yt[i]  = y[i]  + dt * (a31 * ky[i][0]  + a32 * ky[i][1]);
			zt[i]  = z[i]  + dt * (a31 * kz[i][0]  + a32 * kz[i][1]);
			vxt[i] = vx[i] + dt * (a31 * kvx[i][0] + a32 * kvx[i][1]);
			vyt[i] = vy[i] + dt * (a31 * kvy[i][0] + a32 * kvy[i][1]);
			vzt[i] = vz[i] + dt * (a31 * kvz[i][0] + a32 * kvz[i][1]);
		}

		for(int i = 0; i < N; ++i){

			kx[i][2] = vxt[i];
			ky[i][2] = vyt[i];
			kz[i][2] = vzt[i];

			ax = 0.0;
			ay = 0.0;
			az = 0.0;

			if(useHelio == 0){
				for(int j = 0; j < N; ++j){
					if(i != j){
						acc(m, xt, yt, zt, ax, ay, az, i, j);
					}
				}
			}
			else{
				
				if(i > 0){
					accS(m, xt, yt, zt, ax, ay, az, i);
				}
				for(int j = 1; j < N; ++j){
					if(i != j){
						accP(m, xt, yt, zt, ax, ay, az, i, j);
					}
				}
				for(int j = 1; j < N; ++j){
					if(i != j){
						accP2(m, xt, yt, zt, ax, ay, az, j);
					}
				}
				if(GR == 1 && i > 0){
					accGR(m, xt, yt, zt, vxt, vyt, vzt, ax, ay, az, i);
				}
				if(GR == 2 && i > 0){
					accGR2(m, xt, yt, zt, vxt, vyt, vzt, ax, ay, az, i);
				}
				if(GR == 3 && i > 0){
					accGR3(m, xt, yt, zt, vxt, vyt, vzt, ax, ay, az, i);
				}
			}

			kvx[i][2] = ax;
			kvy[i][2] = ay;
			kvz[i][2] = az;

		}
		//stage 4
		for(int i = 0; i < N; ++i){
			xt[i]  = x[i]  + dt * (a41 * kx[i][0]  + a42 * kx[i][1]  + a43 * kx[i][2]);
			yt[i]  = y[i]  + dt * (a41 * ky[i][0]  + a42 * ky[i][1]  + a43 * ky[i][2]);
			zt[i]  = z[i]  + dt * (a41 * kz[i][0]  + a42 * kz[i][1]  + a43 * kz[i][2]);
			vxt[i] = vx[i] + dt * (a41 * kvx[i][0] + a42 * kvx[i][1] + a43 * kvx[i][2]);
			vyt[i] = vy[i] + dt * (a41 * kvy[i][0] + a42 * kvy[i][1] + a43 * kvy[i][2]);
			vzt[i] = vz[i] + dt * (a41 * kvz[i][0] + a42 * kvz[i][1] + a43 * kvz[i][2]);
		}

		for(int i = 0; i < N; ++i){

			kx[i][3] = vxt[i];
			ky[i][3] = vyt[i];
			kz[i][3] = vzt[i];

			ax = 0.0;
			ay = 0.0;
			az = 0.0;

			if(useHelio == 0){
				for(int j = 0; j < N; ++j){
					if(i != j){
						acc(m, xt, yt, zt, ax, ay, az, i, j);
					}
				}
			}
			else{
				
				if(i > 0){
					accS(m, xt, yt, zt, ax, ay, az, i);
				}
				for(int j = 1; j < N; ++j){
					if(i != j){
						accP(m, xt, yt, zt, ax, ay, az, i, j);
					}
				}
				for(int j = 1; j < N; ++j){
					if(i != j){
						accP2(m, xt, yt, zt, ax, ay, az, j);
					}
				}
				if(GR == 1 && i > 0){
					accGR(m, xt, yt, zt, vxt, vyt, vzt, ax, ay, az, i);
				}
				if(GR == 2 && i > 0){
					accGR2(m, xt, yt, zt, vxt, vyt, vzt, ax, ay, az, i);
				}
				if(GR == 3 && i > 0){
					accGR3(m, xt, yt, zt, vxt, vyt, vzt, ax, ay, az, i);
				}
			}

			kvx[i][3] = ax;
			kvy[i][3] = ay;
			kvz[i][3] = az;

		}
		//stage 5
		for(int i = 0; i < N; ++i){
			xt[i]  = x[i]  + dt * (a51 * kx[i][0]  + a52 * kx[i][1]  + a53 * kx[i][2]  + a54 * kx[i][3]);
			yt[i]  = y[i]  + dt * (a51 * ky[i][0]  + a52 * ky[i][1]  + a53 * ky[i][2]  + a54 * ky[i][3]);
			zt[i]  = z[i]  + dt * (a51 * kz[i][0]  + a52 * kz[i][1]  + a53 * kz[i][2]  + a54 * kz[i][3]);
			vxt[i] = vx[i] + dt * (a51 * kvx[i][0] + a52 * kvx[i][1] + a53 * kvx[i][2] + a54 * kvx[i][3]);
			vyt[i] = vy[i] + dt * (a51 * kvy[i][0] + a52 * kvy[i][1] + a53 * kvy[i][2] + a54 * kvy[i][3]);
			vzt[i] = vz[i] + dt * (a51 * kvz[i][0] + a52 * kvz[i][1] + a53 * kvz[i][2] + a54 * kvz[i][3]);
		}

		for(int i = 0; i < N; ++i){

			kx[i][4] = vxt[i];
			ky[i][4] = vyt[i];
			kz[i][4] = vzt[i];

			ax = 0.0;
			ay = 0.0;
			az = 0.0;

			if(useHelio == 0){
				for(int j = 0; j < N; ++j){
					if(i != j){
						acc(m, xt, yt, zt, ax, ay, az, i, j);
					}
				}
			}
			else{
				
				if(i > 0){
					accS(m, xt, yt, zt, ax, ay, az, i);
				}
				for(int j = 1; j < N; ++j){
					if(i != j){
						accP(m, xt, yt, zt, ax, ay, az, i, j);
					}
				}
				for(int j = 1; j < N; ++j){
					if(i != j){
						accP2(m, xt, yt, zt, ax, ay, az, j);
					}
				}
				if(GR == 1 && i > 0){
					accGR(m, xt, yt, zt, vxt, vyt, vzt, ax, ay, az, i);
				}
				if(GR == 2 && i > 0){
					accGR2(m, xt, yt, zt, vxt, vyt, vzt, ax, ay, az, i);
				}
				if(GR == 3 && i > 0){
					accGR3(m, xt, yt, zt, vxt, vyt, vzt, ax, ay, az, i);
				}
			}

			kvx[i][4] = ax;
			kvy[i][4] = ay;
			kvz[i][4] = az;

		}

		//stage 6
		for(int i = 0; i < N; ++i){
			xt[i]  = x[i]  + dt * (a61 * kx[i][0]  + a62 * kx[i][1]  + a63 * kx[i][2]  + a64 * kx[i][3]  + a65 * kx[i][4]);
			yt[i]  = y[i]  + dt * (a61 * ky[i][0]  + a62 * ky[i][1]  + a63 * ky[i][2]  + a64 * ky[i][3]  + a65 * ky[i][4]);
			zt[i]  = z[i]  + dt * (a61 * kz[i][0]  + a62 * kz[i][1]  + a63 * kz[i][2]  + a64 * kz[i][3]  + a65 * kz[i][4]);
			vxt[i] = vx[i] + dt * (a61 * kvx[i][0] + a62 * kvx[i][1] + a63 * kvx[i][2] + a64 * kvx[i][3] + a65 * kvx[i][4]);
			vyt[i] = vy[i] + dt * (a61 * kvy[i][0] + a62 * kvy[i][1] + a63 * kvy[i][2] + a64 * kvy[i][3] + a65 * kvy[i][4]);
			vzt[i] = vz[i] + dt * (a61 * kvz[i][0] + a62 * kvz[i][1] + a63 * kvz[i][2] + a64 * kvz[i][3] + a65 * kvz[i][4]);
		}

		for(int i = 0; i < N; ++i){

			kx[i][5] = vxt[i];
			ky[i][5] = vyt[i];
			kz[i][5] = vzt[i];

			ax = 0.0;
			ay = 0.0;
			az = 0.0;

			if(useHelio == 0){
				for(int j = 0; j < N; ++j){
					if(i != j){
						acc(m, xt, yt, zt, ax, ay, az, i, j);
					}
				}
			}
			else{
				
				if(i > 0){
					accS(m, xt, yt, zt, ax, ay, az, i);
				}
				for(int j = 1; j < N; ++j){
					if(i != j){
						accP(m, xt, yt, zt, ax, ay, az, i, j);
					}
				}
				for(int j = 1; j < N; ++j){
					if(i != j){
						accP2(m, xt, yt, zt, ax, ay, az, j);
					}
				}
				if(GR == 1 && i > 0){
					accGR(m, xt, yt, zt, vxt, vyt, vzt, ax, ay, az, i);
				}
				if(GR == 2 && i > 0){
					accGR2(m, xt, yt, zt, vxt, vyt, vzt, ax, ay, az, i);
				}
				if(GR == 3 && i > 0){
					accGR3(m, xt, yt, zt, vxt, vyt, vzt, ax, ay, az, i);
				}
			}

			kvx[i][5] = ax;
			kvy[i][5] = ay;
			kvz[i][5] = az;

		}
		double sc = 1.0e-15;

		
		//error estimation
		for(int i = 0; i < N; ++i){
			errorkx[i] = ((b1 - bb1) * kx[i][0] + (b3 - bb3) * kx[i][2] + (b4 - bb4) * kx[i][3] + (b5 - bb5) * kx[i][4] + (b6 - bb6) * kx[i][5]) / sc;
			errorky[i] = ((b1 - bb1) * ky[i][0] + (b3 - bb3) * ky[i][2] + (b4 - bb4) * ky[i][3] + (b5 - bb5) * ky[i][4] + (b6 - bb6) * ky[i][5]) / sc;
			errorkz[i] = ((b1 - bb1) * kz[i][0] + (b3 - bb3) * kz[i][2] + (b4 - bb4) * kz[i][3] + (b5 - bb5) * kz[i][4] + (b6 - bb6) * kz[i][5]) / sc;
			errorkvx[i] = ((b1 - bb1) * kvx[i][0] + (b3 - bb3) * kvx[i][2] + (b4 - bb4) * kvx[i][3] + (b5 - bb5) * kvx[i][4] + (b6 - bb6) * kvx[i][5]) / sc;
			errorkvy[i] = ((b1 - bb1) * kvy[i][0] + (b3 - bb3) * kvy[i][2] + (b4 - bb4) * kvy[i][3] + (b5 - bb5) * kvy[i][4] + (b6 - bb6) * kvy[i][5]) / sc;
			errorkvz[i] = ((b1 - bb1) * kvz[i][0] + (b3 - bb3) * kvz[i][2] + (b4 - bb4) * kvz[i][3] + (b5 - bb5) * kvz[i][4] + (b6 - bb6) * kvz[i][5]) / sc;
		}


		double error = 0.0;
		for(int i = 0; i < N; ++i){
			error += errorkx[i] * errorkx[i];
			error += errorky[i] * errorky[i];
			error += errorkz[i] * errorkz[i];
			error += errorkvx[i] * errorkvx[i];
			error += errorkvy[i] * errorkvy[i];
			error += errorkvz[i] * errorkvz[i];
		}

		/*
		double errmax = 0.0;
		for(int i = 0; i < N; ++i){
			ermax = fmax(ermax, fabs(errorkx[i]));
			ermax = fmax(ermax, fabs(errorky[i]));
			ermax = fmax(ermax, fabs(errorkz[i]));
			ermax = fmax(ermax, fabs(errorkvx[i]));
			ermax = fmax(ermax, fabs(errorkvy[i]));
			ermax = fmax(ermax, fabs(errorkvz[i]));
		}
		*/


		double ee = 1.0/5.0;	

		double s = pow( 1.0  / error, ee);

		//printf("%g %g\n", dt, s);			

		for(int i = 0; i < N; ++i){
		
			//RKF45
			x[i] += dt * (b1 * kx[i][0] + b3 * kx[i][2] + b4 * kx[i][3] + b5 * kx[i][4] + b6 * kx[i][5]);
			y[i] += dt * (b1 * ky[i][0] + b3 * ky[i][2] + b4 * ky[i][3] + b5 * ky[i][4] + b6 * ky[i][5]);
			z[i] += dt * (b1 * kz[i][0] + b3 * kz[i][2] + b4 * kz[i][3] + b5 * kz[i][4] + b6 * kz[i][5]);

			vx[i] += dt * (b1 * kvx[i][0] + b3 * kvx[i][2] + b4 * kvx[i][3] + b5 * kvx[i][4] + b6 * kvx[i][5]);
			vy[i] += dt * (b1 * kvy[i][0] + b3 * kvy[i][2] + b4 * kvy[i][3] + b5 * kvy[i][4] + b6 * kvy[i][5]);
			vz[i] += dt * (b1 * kvz[i][0] + b3 * kvz[i][2] + b4 * kvz[i][3] + b5 * kvz[i][4] + b6 * kvz[i][5]);

		}

		time += dt / dayUnit;



		if(t % outInterval == 0){
			
			if(outHelio == 1){
				sprintf(outfilename, "Outhelio10_%.12lld.dat", t);
			}
			else{
				sprintf(outfilename, "Outbary10_%.12lld.dat", t);
			}
			outfile = fopen(outfilename, "w");
//			printf("%s\n", outfilename);
			
			comx = 0.0;
			comy = 0.0;
			comz = 0.0;
			vcomx = 0.0;
			vcomy = 0.0;
			vcomz = 0.0;
			
			if(useHelio == 0 && outHelio == 1){
				//convert to heliocentric output
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
				fprintf(outfile, "%.10g %d %.40g %.40g %.40g %.40g %.40g %.40g %.40g %.40g 0 0 0 0 0 0 0 0 0 0 0\n", time, i, m[i], rad[i], comx + x[i], comy + y[i], comz + z[i], vcomx + vx[i], vcomy + vy[i], vcomz + vz[i]);
				
			}
			fclose(outfile);
		}
		/*
		if( s < 1.0){
			dt *= 0.5;
		}
		if( s > 3.0){
			dt *= 2.0;
		}
	*/
	}
}
	

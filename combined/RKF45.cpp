#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define dayUnit 0.01720209895
//#define dayUnit 0.01720209894846
//#define dayUnit 1.0

//#define def_c 10065.3201686
#define def_c 10065.320121
#define def_AU 149597870700.0           //AU in m

// --------------------------------
//barycentric coordinates
void acc(double *m, double *x, double *y, double *z, double &ax, double &ay, double &az, int i, int j){

	double rx = x[j] - x[i];
	double ry = y[j] - y[i];
	double rz = z[j] - z[i];
	double rsq = rx * rx + ry * ry + rz * rz;
	double r = sqrt(rsq);

	double s = m[j] / (r * rsq);
	
	ax += s * rx;
	ay += s * ry;
	az += s * rz;
}
// ---------------------------------

// --------------------------------
//heliocentric coordinates
//sun part
void accS(double *m, double *x, double *y, double *z, double &ax, double &ay, double &az, int i){

	if(i > 0){
		double rx = -x[i];
		double ry = -y[i];
		double rz = -z[i];
		double rsq = rx * rx + ry * ry + rz * rz;
		double r = sqrt(rsq);

		double s = (m[0] + m[i]) / (r * rsq);

		ax += s * rx;
		ay += s * ry;
		az += s * rz;
	}
}
//planet part
void accP(double *m, double *x, double *y, double *z, double &ax, double &ay, double &az, int i, int j){

	if(i != j){
		double rx = x[j] - x[i];
		double ry = y[j] - y[i];
		double rz = z[j] - z[i];
		double rsq = rx * rx + ry * ry + rz * rz;
		double r = sqrt(rsq);

		double s = m[j] / (r * rsq);

		ax += s * rx;
		ay += s * ry;
		az += s * rz;
	}
}
//planet part 2
void accP2(double *m, double *x, double *y, double *z, double &ax, double &ay, double &az, int i, int j){

	if(i != j){
		double rx = -x[j];
		double ry = -y[j];
		double rz = -z[j];
		double rsq = rx * rx + ry * ry + rz * rz;
		double r = sqrt(rsq);

		double s = m[j] / (r * rsq);
	
		ax += s * rx;
		ay += s * ry;
		az += s * rz;
	}
}
// --------------------------------------

//GR Quinn, Tremaine, Duncan 1991, heliocentric
//heliocentric
void acchGR(double *m, double *x, double *y, double *z, double *vx, double *vy, double *vz, double &ax, double &ay, double &az, int i){
	
	double c2 = def_c * def_c;

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
//barycentric
void accbGR(double *m, double *x, double *y, double *z, double *vx, double *vy, double *vz, double &ax, double &ay, double &az, int i){
	
	double c2 = def_c * def_c;

	double rhx = x[i] - x[0];
	double rhy = y[i] - y[0];
	double rhz = z[i] - z[0];

	double vhx = vx[i] - vx[0];
	double vhy = vy[i] - vy[0];
	double vhz = vz[i] - vz[0];

	double rsq = rhx * rhx + rhy * rhy + rhz * rhz;
	double r = sqrt(rsq);
	double vsq = vhx * vhx + vhy * vhy + vhz * vhz;

	double rv = rhx * vhx + rhy * vhy + rhz * vhz;

	double t1 = 2.0 * (1.0 + 1.0) * m[0] / (c2 * r) - 1.0 * vsq / c2;
	double t2 = (2.0 + 2.0 * 1.0) * m[0] / (c2 * r * rsq) * rv;

	ax += m[0] / (r * rsq) * rhx * t1 + t2 * vhx;
	ay += m[0] / (r * rsq) * rhy * t1 + t2 * vhy;
	az += m[0] / (r * rsq) * rhz * t1 + t2 * vhz;

}

//Sitarski 1982, Isotropic equation 5, heliocentric
//modified k2 to dayUnit
//should be equivalent to the Quinn et all function, assuming m[0] = 1.0
//heliocentric
void acchGR2(double *m, double *x, double *y, double *z, double *vx, double *vy, double *vz, double &ax, double &ay, double &az, int i){
	
	double c2 = def_c * def_c;

	double rsq = x[i] * x[i] + y[i] * y[i] + z[i] * z[i];
	double r = sqrt(rsq);
	double vsq = vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];

	double rv = x[i] * vx[i] + y[i] * vy[i] + z[i] * vz[i];

	double f1 = 1.0 / (r * rsq * c2);
	double t1 = 4.0 / r;
	double t2 = -vsq;
	double t3 = 4.0 * rv;
//printf("a %d %.20g %.20g %.20g\n", i, ax, ay, az);

//printf("A %d %.20g %.20g %.20g %.20g %.20g %.20g\n", i, x[i], y[i], z[i], vx[i], vy[i], vz[i]);
//printf("B %d %.20g %.20g %.20g %.20g\n", i, f1, t1, t2, t3);
//printf("C %d %.20g %.20g %.20g %.20g\n", i, t1 + t2, (t1 + t2) * x[i], ((t1 + t2) * x[i] + t3 * vx[i]), f1 * ((t1 + t2) * x[i] + t3 * vx[i]));


	ax += f1 * ((t1 + t2) * x[i] + t3 * vx[i]);
	ay += f1 * ((t1 + t2) * y[i] + t3 * vy[i]);
	az += f1 * ((t1 + t2) * z[i] + t3 * vz[i]);
//printf("D %d %.20g %.20g %.20g\n", i, ax, ay, az);

}
//barycentric
void accbGR2(double *m, double *x, double *y, double *z, double *vx, double *vy, double *vz, double &ax, double &ay, double &az, int i){
	
	double c2 = def_c * def_c;

	double rhx = x[i] - x[0];
	double rhy = y[i] - y[0];
	double rhz = z[i] - z[0];

	double vhx = vx[i] - vx[0];
	double vhy = vy[i] - vy[0];
	double vhz = vz[i] - vz[0];

	double rsq = rhx * rhx + rhy * rhy + rhz * rhz;
	double r = sqrt(rsq);
	double vsq = vhx * vhx + vhy * vhy + vhz * vhz;

	double rv = rhx * vhx + rhy * vhy + rhz * vhz;

	double f1 = 1.0 / (r * rsq * c2);
	double t1 = 4.0 / r;
	double t2 = -vsq;
	double t3 = 4.0 * rv;
//printf("a %d %.20g %.20g %.20g\n", i, ax, ay, az);

//printf("A %d %.20g %.20g %.20g %.20g %.20g %.20g\n", i, rhx, rhy, rhz, vhx, vhy, vhz);
//printf("B %d %.20g %.20g %.20g %.20g\n", i, f1, t1, t2, t3);
//printf("C %d %.20g %.20g %.20g %.20g\n", i, t1 + t2, (t1 + t2) * rhx, ((t1 + t2) * rhx + t3 * vhx), f1 * ((t1 + t2) * rhx + t3 * vhx));

	ax += f1 * ((t1 + t2) * rhx + t3 * vhx);
	ay += f1 * ((t1 + t2) * rhy + t3 * vhy);
	az += f1 * ((t1 + t2) * rhz + t3 * vhz);
//printf("D%d %.20g %.20g %.20g\n", i, ax, ay, az);

}

//Fabricky 2010, heliocentric
//heliocentric
void acchGR3(double *m, double *x, double *y, double *z, double *vx, double *vy, double *vz, double &ax, double &ay, double &az, int i){
	
	double c2 = def_c * def_c;

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
//barycentric
void accbGR3(double *m, double *x, double *y, double *z, double *vx, double *vy, double *vz, double &ax, double &ay, double &az, int i){
	
	double c2 = def_c * def_c;

	double rhx = x[i] - x[0];
	double rhy = y[i] - y[0];
	double rhz = z[i] - z[0];

	double vhx = vx[i] - vx[0];
	double vhy = vy[i] - vy[0];
	double vhz = vz[i] - vz[0];

	double rsq = rhx * rhx + rhy * rhy + rhz * rhz;
	double r = sqrt(rsq);

	double A = (m[0] + m[i]) / r;
	double B = A / ( r * c2);
	double eta = m[0] * m[i] / ((m[0] + m[i]) * (m[0] + m[i]));
	double vsq = (vhx * vhx + vhy * vhy + vhz * vhz);
	double rd = (rhx * vhx + rhy * vhy + rhz * vhz) / r; 

	double C = 2.0 * (2.0 - eta) * rd;
	double D = (1.0 + 3.0 * eta) * vsq - 1.5 * eta * rd * rd - 2.0 * (2.0 + eta) * A;
	ax += B * (C * vhx - D * rhx / r);         
	ay += B * (C * vhy - D * rhy / r);
	az += B * (C * vhz - D * rhz / r);

}

//Moyer, barycentric
void accbGR4(double *m, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *aNx, double *aNy, double *aNz, double &ax, double &ay, double &az, int N, int i, int j){
	
	double c2 = def_c * def_c;
	double beta = 1.0;
	double gamma = 1.0;

	//if(j != i){
	if(j < 1 && j != i){
		double dxij = x[j] - x[i];
		double dyij = y[j] - y[i];
		double dzij = z[j] - z[i];
		double rijsq = dxij * dxij + dyij * dyij + dzij * dzij;
		double rij = sqrt(rijsq);

		double sdisq = vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
		double sdjsq = vx[j] * vx[j] + vy[j] * vy[j] + vz[j] * vz[j];
//printf("v %d %d %g %g %g %g\n", i, j, vx[i], vx[j], vy[i], vy[j]);
		double rdirdj = vx[i] * vx[j] + vy[i] * vy[j] + vz[i] * vz[j];

		double drijdt = -(dxij * vx[j] + dyij * vy[j] + dzij * vz[j]) / rij;
//printf("aN %d %d %g %g %g\n", i, j, aNx[j], aNy[j], aNz[j]);
		double rijrjdd = (dxij * aNx[j] + dyij * aNy[j] + dzij * aNz[j]);


		double Sl = 0.0;	//Sum over muL/ril
		for(int l = 0; l < N; ++l){
		//for(int l = 0; l < 1; ++l){
			if(l != i){
				double dxil = x[l] - x[i];
				double dyil = y[l] - y[i];
				double dzil = z[l] - z[i];
				double rilsq = dxil * dxil + dyil * dyil + dzil * dzil;
				double ril = sqrt(rilsq);
				
				Sl += m[l] / ril;
			}
		}
		double Sk = 0.0;	//Sum over muk/rjk
		for(int k = 0; k < N; ++k){
		//for(int k = 0; k < 1; ++k){
			if(k != j){
				double dxjk = x[k] - x[j];
				double dyjk = y[k] - y[j];
				double dzjk = z[k] - z[j];
				double rjksq = dxjk * dxjk + dyjk * dyjk + dzjk * dzjk;
				double rjk = sqrt(rjksq);
				
				Sk += m[k] / rjk;
			}
		}

		double t1 = -2.0 * (beta + gamma) / c2 * Sl;
		double t2 = -(2.0 * beta - 1.0) / c2 * Sk;
		double t3 = gamma * sdisq / c2;
		double t4 = (1.0 + gamma) * sdjsq / c2;
		double t5 = -2.0 * (1.0 + gamma) / c2 * rdirdj;
		double t6 = -3.0 / (2.0 * c2) * drijdt * drijdt;
		double t7 = 1.0 / (2.0 * c2) * rijrjdd;

printf("A %d %d %g %g %g %g %g %g %g\n", i,j, t1, t2, t3, t4, t5, t6, t7);

		double t = t1 + t2 + t3 + t4 + t5 + t6 + t7;
		double tx = t * dxij / (rijsq * rij) * m[j];
		double ty = t * dyij / (rijsq * rij) * m[j];
		double tz = t * dzij / (rijsq * rij) * m[j];


		double u1x = (2.0 + 2.0 * gamma) * vx[i] - (1.0 + 2.0 * gamma) * vx[j];
		double u1y = (2.0 + 2.0 * gamma) * vy[i] - (1.0 + 2.0 * gamma) * vy[j];
		double u1z = (2.0 + 2.0 * gamma) * vz[i] - (1.0 + 2.0 * gamma) * vz[j];

		double u1 = -dxij * u1x - dyij * u1y - dzij * u1z;
		double u2 = m[j] / (c2 * rijsq * rij);		

		double ux = u2 * u1 * (vx[i] - vx[j]);
		double uy = u2 * u1 * (vy[i] - vy[j]);
		double uz = u2 * u1 * (vz[i] - vz[j]);
printf("B %d %d %g %g %g %g %g %g %g\n", i,j, t, tx, ty, tz, ux, uy, uz);

		
		double w1 = (3.0 + 4.0 * gamma) / (2.0 * c2) * m[j] / rij;

		double wx = w1 * aNx[j];
		double wy = w1 * aNy[j];
		double wz = w1 * aNz[j];

		ax += tx + ux + wx;
		ay += ty + uy + wy;
		az += tz + uz + wz;
printf("C %d %d %g %g %g\n", i,j, wx, wy, wz);

	}

}
// --------------------------------


//Neville-Aitken interpolation
void interpolate(int Nperturbers, int N, double *xp, double *yp, double *zp, double *vxp, double *vyp, double *vzp, double *timep, double time, double *xt, double *yt, double *zt, double *vxt, double *vyt, double *vzt){

	for(int p = 0; p < Nperturbers; ++p){
		double r3[6];

		for(int k = 0; k < 6; ++k){

			double P[N][N];
			double tn[N];

			for(int i = 0; i < N; ++i){
				if(k == 0) P[0][i] = xp[p * N + i];
				if(k == 1) P[0][i] = yp[p * N + i];
				if(k == 2) P[0][i] = zp[p * N + i];
				if(k == 3) P[0][i] = vxp[p * N + i];
				if(k == 4) P[0][i] = vyp[p * N + i];
				if(k == 5) P[0][i] = vzp[p * N + i];
				tn[i] = timep[p * N + i];


//printf("interpolate %d %d %d %.20g %.20g %.20g\n", p, i, k, time, tn[i], P[0][i]);
			}

			for(int j = 1; j < N; ++j){
//printf("****\n");
				for(int i = 0; i < N - j; ++i){
					P[j][i] = ((time - tn[i+j]) * P[j-1][i] + (tn[i] - time) * P[j-1][i+1]) / (tn[i] - tn[i+j]);
//printf("%d %d %g %g %g %g %.20g\n", i, i+j, tn[i], tn[i+j], P[j-1][i], P[j-1][i+1], P[j][i]);

				}
			}
			r3[k] = P[N-1][0];

		}
		xt[p] = r3[0];
		yt[p] = r3[1];
		zt[p] = r3[2];
		vxt[p] = r3[3];
		vyt[p] = r3[4];
		vzt[p] = r3[5];
//printf("%.20g %d %.20g %.20g %.20g %.20g %.20g %.20g\n", time, p, r3[0], r3[1], r3[2], r3[3], r3[4], r3[5]);

	}

}


void Newtonian(double *m, double *xt, double *yt, double *zt, double *aNx, double *aNy, double *aNz, double *aNxb, double *aNyb, double *aNzb, int N, int Nperturbers, int useHelio){

	for(int i = 0; i < N; ++i){
		aNx[i] = 0.0;
		aNy[i] = 0.0;
		aNz[i] = 0.0;

		if(useHelio == 0){
			//for(int j = 0; j < Nperturbers; ++j){
			for(int j = Nperturbers-1; j >= 0; --j){
				accP(m, xt, yt, zt, aNx[i], aNy[i], aNz[i], i, j);
<<<<<<< HEAD
//if(i == 27) printf("Nij %d %d %.20g %.20g %.20g\n", i, j, aNx[i] * dayUnit * dayUnit, aNy[i] * dayUnit * dayUnit, aNz[i] * dayUnit * dayUnit);
			}
//printf("N0 %d %.20g %.20g %.20g\n", i, aNx[i] * dayUnit * dayUnit, aNy[i] * dayUnit * dayUnit, aNz[i] * dayUnit * dayUnit);
=======
if(i == 27) printf("Nij %d %d %.20g %.20g %.20g\n", i, j, aNx[i] * dayUnit * dayUnit, aNy[i] * dayUnit * dayUnit, aNz[i] * dayUnit * dayUnit);
			}
printf("N0 %d %.20g %.20g %.20g\n", i, aNx[i] * dayUnit * dayUnit, aNy[i] * dayUnit * dayUnit, aNz[i] * dayUnit * dayUnit);
>>>>>>> f28efb018735bbd62b2845e662c43faea3666aee
		}
		else{
			//for(int j = 1; j < Nperturbers; ++j){
			for(int j = Nperturbers-1; j >= 1; --j){
				accP(m, xt, yt, zt, aNx[i], aNy[i], aNz[i], i, j);
<<<<<<< HEAD
//if(i == 27) printf("Nij %d %d %.20g %.20g %.20g\n", i, j, aNx[i] * dayUnit * dayUnit, aNy[i] * dayUnit * dayUnit, aNz[i] * dayUnit * dayUnit);
			}
			accS(m, xt, yt, zt, aNx[i], aNy[i], aNz[i], i);
//if(i == 27) printf("Nij %d %d %.20g %.20g %.20g\n", i, 0, aNx[i] * dayUnit * dayUnit, aNy[i] * dayUnit * dayUnit, aNz[i] * dayUnit * dayUnit);
//if(i == 27) printf("N0 %d %.20g %.20g %.20g\n", i, aNx[i] * dayUnit * dayUnit, aNy[i] * dayUnit * dayUnit, aNz[i] * dayUnit * dayUnit);
=======
if(i == 27) printf("Nij %d %d %.20g %.20g %.20g\n", i, j, aNx[i] * dayUnit * dayUnit, aNy[i] * dayUnit * dayUnit, aNz[i] * dayUnit * dayUnit);
			}
			accS(m, xt, yt, zt, aNx[i], aNy[i], aNz[i], i);
if(i == 27) printf("Nij %d %d %.20g %.20g %.20g\n", i, 0, aNx[i] * dayUnit * dayUnit, aNy[i] * dayUnit * dayUnit, aNz[i] * dayUnit * dayUnit);
printf("N0 %d %.20g %.20g %.20g\n", i, aNx[i] * dayUnit * dayUnit, aNy[i] * dayUnit * dayUnit, aNz[i] * dayUnit * dayUnit);
>>>>>>> f28efb018735bbd62b2845e662c43faea3666aee

			aNxb[i] = aNx[i];
			aNyb[i] = aNy[i];
			aNzb[i] = aNz[i];

			//for(int j = 1; j < Nperturbers; ++j){
			for(int j = Nperturbers-1; j >= 1; --j){
				accP2(m, xt, yt, zt, aNx[i], aNy[i], aNz[i], i, j);
			}
<<<<<<< HEAD
//if(i == 27) printf("Np %d %.20g %.20g %.20g\n", i, aNx[i] * dayUnit * dayUnit, aNy[i] * dayUnit * dayUnit, aNz[i] * dayUnit * dayUnit);
=======
//printf("Np %d %.20g %.20g %.20g\n", i, aNx[i] * dayUnit * dayUnit, aNy[i] * dayUnit * dayUnit, aNz[i] * dayUnit * dayUnit);
>>>>>>> f28efb018735bbd62b2845e662c43faea3666aee
		}

	}
}

void GRCall(double *m, double *xt, double *yt, double *zt, double *vxt, double *vyt, double *vzt, double &ax, double &ay, double &az, double *aNx, double *aNy, double *aNz, int N, int Nperturbers, int useHelio, int GR, int i){

	if(useHelio == 0){
		if(GR == 1){
			accbGR(m, xt, yt, zt, vxt, vyt, vzt, ax, ay, az, i);
		}
		if(GR == 2){
			accbGR2(m, xt, yt, zt, vxt, vyt, vzt, ax, ay, az, i);
		}
		if(GR == 3){
			accbGR3(m, xt, yt, zt, vxt, vyt, vzt, ax, ay, az, i);
		}
		if(GR == 4){
			for(int j = 0; j < Nperturbers; ++j){
				accbGR4(m, xt, yt, zt, vxt, vyt, vzt, aNx, aNy, aNz, ax, ay, az, Nperturbers, i, j);
			}
		}
	}
	else{
		if(GR == 1){
			acchGR(m, xt, yt, zt, vxt, vyt, vzt, ax, ay, az, i);
		}
		if(GR == 2){
			acchGR2(m, xt, yt, zt, vxt, vyt, vzt, ax, ay, az, i);
		}
		if(GR == 3){
			acchGR3(m, xt, yt, zt, vxt, vyt, vzt, ax, ay, az, i);
		}
		if(GR == 4){
			for(int j = 0; j < Nperturbers; ++j){
				accbGR4(m, xt, yt, zt, vxt, vyt, vzt, aNx, aNy, aNz, ax, ay, az, Nperturbers, i, j);
			}
		}
	}
}

//A1, A2 and A3 terms for asteroids on heliocentric coordinates
void NonGrav(double *x, double *y, double *z, double *vx, double *vy, double *vz, double &ax, double &ay, double &az, double *A1, double *A2, double *A3, double *ALN, double *NK, double *NM, double *Nn, double *R0, int i){

	double rsq = x[i] * x[i] + y[i] * y[i] + z[i] * z[i];
	double r = sqrt(rsq);

	//angular momenrum h = r x v
	double hx = y[i] * vz[i] - z[i] * vy[i];
	double hy =-x[i] * vz[i] + z[i] * vx[i];
	double hz = x[i] * vy[i] - y[i] * vx[i];

	double hsq = hx * hx + hy * hy + hz * hz;
	double h = sqrt(hsq);

	//Transverse velocity t = h x r
	double tx = hy * z[i] - hz * y[i];
	double ty =-hx * z[i] + hz * x[i];
	double tz = hx * y[i] - hy * x[i];

	double tsq = tx * tx + ty * ty + tz * tz;
	double t = sqrt(tsq);

	double gr = 1.0 / rsq;	//only valid for asteroids, not for comets 
/*
	double rr = r / R0[i];
	double g1 = pow(rr, -NM[i]);
	double g2 = pow(rr, Nn[i]);
	double g3 = pow(1.0 + g2, -NK[i]);
	double gr = ALN[i] * g1 * g3;

printf("gr %.20g %.20g\n", gr1, gr);
*/

	double f1 = A1[i] * gr / r;
	double f2 = A2[i] * gr / t;
	double f3 = A3[i] * gr / h;
	
	
	ax += f1 * x[i] + f2 * tx + f3 * hx;
	ay += f1 * y[i] + f2 * ty + f3 * hy;
	az += f1 * z[i] + f2 * tz + f3 * hz;
printf("NonGrav %d %.20g %.20g %.20g\n", i, (f1 * x[i] + f2 * tx + f3 * hx) * dayUnit * dayUnit, (f1 * y[i] + f2 * ty + f3 * hy) * dayUnit * dayUnit, (f1 * z[i] + f2 * tz + f3 * hz) * dayUnit * dayUnit);

}

//J2 perturbation from Earth
void J2(double *m, double *x, double *y, double *z, double &ax, double &ay, double &az, int i){

	double J2E = 0.00108262545; // J2 Earth from DE 430
	double RE = 6378136.3; // Earth radius in m from DE 430

	//double J2E =  1.08263e-3; //  1.08262668e-3;
	//double RE = 6371.009; // Earth radius in km
	//double muE = 398600.44 // in km^3 /s^2	G * mEarth

	RE /= def_AU;	//Earth radius in AU

	int iE = 3; 	//index of Earth

	double muE = m[iE];

	double xE = x[i] - x[iE];
	double yE = y[i] - y[iE];
	double zE = z[i] - z[iE];

	double rsq = xE * xE + yE * yE + zE * zE;
	double r = sqrt(rsq);
	double r5 = rsq * rsq * r;

	double t1 = 3.0 * J2E * muE * RE * RE / (2.0 * r5);
	double t2 = 5.0 * zE * zE / rsq;

printf("rm %.20g %.20g %.20g\n", RE, muE, t1);

<<<<<<< HEAD
	double tx = t1 * (t2 - 1.0) * xE;
	double ty = t1 * (t2 - 1.0) * yE;
	double tz = t1 * (t2 - 3.0) * zE;
=======
	double tx = t1 * (t2 - 1.0) * xE;// * xE / r;
	double ty = t1 * (t2 - 1.0) * yE;// * yE / r;
	double tz = t1 * (t2 - 3.0) * zE;// * zE / r;
>>>>>>> f28efb018735bbd62b2845e662c43faea3666aee
	
	ax += tx;
	ay += ty;
	az += tz;
 
printf("J2 %d %.20g %.20g %.20g | %.20g %.20g %.20g\n", i, tx * dayUnit * dayUnit, ty * dayUnit * dayUnit, tz * dayUnit * dayUnit, xE, yE, zE); 


}


int main(int argc, char*argv[]){

	//Number of planets
	const int NN = 28;	//22
	const int Nperturbers = 27;	//21
	const int Ninterpolate = 10;	//number of interpolation points
	const double dtime = 1.0;     //interval between stored time steps

	int GR = 4;
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

	//long long int Nsteps = 400000;	
<<<<<<< HEAD
	long long int Nsteps = 40000;
=======
	long long int Nsteps = 1;//4400000;	
>>>>>>> f28efb018735bbd62b2845e662c43faea3666aee
	long long int outInterval = 100;
	double dt = 0.01 * dayUnit;

	for(int i = 1; i < argc; i += 2){

		if(strcmp(argv[i], "-Nsteps") == 0){
			Nsteps = atoll(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-outInterval") == 0){
			outInterval = atoll(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-dt") == 0){
			dt = atof(argv[i + 1]);
		}
	}

	
	int id[NN];
	double m[NN];
	double x[NN];
	double y[NN];
	double z[NN];
	double vx[NN];
	double vy[NN];
	double vz[NN];

	double aNx[NN];	//Newtonian acceleration, used for GR=4
	double aNy[NN];
	double aNz[NN];

	double aNxb[NN];//Newtonian acceleration, used for GR=4, barycentric
	double aNyb[NN];
	double aNzb[NN];


	//non Gravitational constants
	double A1[NN];
	double A2[NN];
	double A3[NN];
	double ALN[NN];
	double NK[NN];
	double NM[NN];
	double Nn[NN];
	double R0[NN];

	for(int i = 0; i < NN; ++i){
		A1[i] = 0.0;
		A2[i] = 0.0;
		A3[i] = 0.0;
		ALN[i] = 0.0;
		NK[i] = 0.0;
		NM[i] = 0.0;
		Nn[i] = 0.0;
		R0[i] = 1.0;

	}

	FILE *outfile;
	char outfilename[160];	


	double time = 0.0;

	//Units are 1/(mass of object in solar masses)
	double pmass[] = {
		1.000000000000000e0,      // Sun        (0)
		6023682.155592479e0,      // Mercury    (1)
		408523.7186582996e0,      // Venus      (2)
		332946.0488339480e0,      // Earth      (3)
		3098703.590290707e0,      // Mars       (4)
		1047.348625463337e0,      // Jupiter    (5)
		3497.901767786633e0,      // Saturn     (6)
		22902.98161308703e0,      // Uranus     (7)
		19412.25977597307e0,      // Neptune    (8)
		135836683.7686175e0,      // Pluto      (9)
		2112939391.8192508,                  // Ceres      (10)
		9531877787.0654011,                  // Pallas     (11)
		81799329362.428986,                  // Juno       (12)
		7676559929.1351004,                  // Vesta      (13)
		23944976514.662392,                  // Hygiea     (14)
		63251980219.354561,                  // Eunomia    (15)
		46649712166.264168,                  // Euphrosyne (16)
		119474172269.94408,                  // Europa     (17)
		56926698684.931702,                  // Davida     (18)
		56298080671.641434,                  // Interamnia (19)
		27068703.24120323e0,      // Moon       (20)

		86737410876.841156,		  // Psyche
		93034865412.812271,		  // Cybele
		114823090351.20033,		  // Thisbe
		116910898662.48077,		  // Doris
		128906361339.41116,		  // Patientia
		134548655333.38321,		  // Sylvia

		0.0			  // test particle

	};
	for(int i = 0; i < Nperturbers; ++i){
		m[i] = 1.0/pmass[i];
printf("m %d %.20g\n", i, m[i]);
	}

//m[Nperturbers] = 1.e-11; //ca mass of Flora

	int N = Nperturbers;

	//Sun
	FILE *infile;
	char infilename[160];

	sprintf(infilename, "initial.dat");
	infile = fopen(infilename, "r");
	//sun
	id[0] = 20;
	x[0] = 0.0;
	y[0] = 0.0;
	z[0] = 0.0;
	vx[0] = 0.0;
	vy[0] = 0.0;
	vz[0] = 0.0;

	for(int i = 1; i < NN; ++i){
		id[i] = -1;
		x[i] = 0.0;
		y[i] = 0.0;
		z[i] = 0.0;
		vx[i] = 0.0;
		vy[i] = 0.0;
		vz[i] = 0.0;
	}
	//read test particle
	for(int i = Nperturbers; i < NN; ++i){
		int er = 0;
		fscanf(infile, "%lf", &time);
		fscanf(infile, "%lf", &x[i]);
		fscanf(infile, "%lf", &y[i]);
		fscanf(infile, "%lf", &z[i]);
		fscanf(infile, "%lf", &vx[i]);
		fscanf(infile, "%lf", &vy[i]);
<<<<<<< HEAD
		er = fscanf(infile, "%lf", &vz[i]);
		//fscanf(infile, "%lf", &A1[i]);
		//fscanf(infile, "%lf", &A2[i]);
		//fscanf(infile, "%lf", &A3[i]);
		//fscanf(infile, "%lf", &ALN[i]);
		//fscanf(infile, "%lf", &NK[i]);
		//fscanf(infile, "%lf", &NM[i]);
		//fscanf(infile, "%lf", &Nn[i]);
		//er = fscanf(infile, "%lf", &R0[i]);
=======
		fscanf(infile, "%lf", &vz[i]);
		fscanf(infile, "%lf", &A1[i]);
		fscanf(infile, "%lf", &A2[i]);
		fscanf(infile, "%lf", &A3[i]);
		fscanf(infile, "%lf", &ALN[i]);
		fscanf(infile, "%lf", &NK[i]);
		fscanf(infile, "%lf", &NM[i]);
		fscanf(infile, "%lf", &Nn[i]);
		er = fscanf(infile, "%lf", &R0[i]);
>>>>>>> f28efb018735bbd62b2845e662c43faea3666aee
		if(er < 0) break;
		++N;
printf("er %d %d %d %d %.20g %.20g %.20g\n", i, id[i], er, N, x[i], y[i], z[i]);
	}
	fclose(infile);
	double time0 = time;	//start time from simulation
	double time1 = time;	//time from table position

	for(int i = Nperturbers; i < N; ++i){
		vx[i] /= dayUnit;
		vy[i] /= dayUnit;
		vz[i] /= dayUnit;

		A1[i] /= (dayUnit * dayUnit);
		A2[i] /= (dayUnit * dayUnit);
		A3[i] /= (dayUnit * dayUnit);
	}


	//interpolate initial values

	double timep[Nperturbers * Ninterpolate];
	double xp[Nperturbers * Ninterpolate];
	double yp[Nperturbers * Ninterpolate];
	double zp[Nperturbers * Ninterpolate];
	double vxp[Nperturbers * Ninterpolate];
	double vyp[Nperturbers * Ninterpolate];
	double vzp[Nperturbers * Ninterpolate];

	//Barycenter
	double xpbc[Ninterpolate];
	double ypbc[Ninterpolate];
	double zpbc[Ninterpolate];
	double vxpbc[Ninterpolate];
	double vypbc[Ninterpolate];
	double vzpbc[Ninterpolate];

	double xbc[1];
	double ybc[1];
	double zbc[1];
	double vxbc[1];
	double vybc[1];
	double vzbc[1];



	FILE *XVfile;
	FILE *BCfile;	//Barycenter
	if(useHelio == 1){
		XVfile = fopen("All_h.dat", "r");
		BCfile = fopen("Barycenter_h.dat", "r");
	}
	else{
		XVfile = fopen("All_b.dat", "r");
	}


	// -----------------------------------------
	//Read table
	int countNodes = 0;
	for(int t = 0; t < 1000000; ++t){
		int er;
printf("CountNodes %d\n", countNodes);
		for(int i = 0; i < Nperturbers; ++i){
			double skip;
			double timepp;
			int id;
			er = fscanf(XVfile, "%lf %d", &timepp, &id);
			fscanf(XVfile, "%lf %lf %lf", &xp[id * Ninterpolate + countNodes], &yp[id * Ninterpolate + countNodes], &zp[id * Ninterpolate + countNodes]);
			fscanf(XVfile, "%lf %lf %lf", &vxp[id * Ninterpolate + countNodes], &vyp[id * Ninterpolate + countNodes], &vzp[id * Ninterpolate + countNodes]);

			if(useHelio == 1 && i == 0){
				double bcTime;
				fscanf(BCfile, "%lf", &bcTime);
				fscanf(BCfile, "%lf %lf %lf", &xpbc[countNodes], &ypbc[countNodes], &zpbc[countNodes]);
				fscanf(BCfile, "%lf %lf %lf", &vxpbc[countNodes], &vypbc[countNodes], &vzpbc[countNodes]);
			}

			if(er < 0) break;
			timep[id * Ninterpolate + countNodes] = timepp;
printf("read %.20g %d %.20g %.20g %.20g %.20g %.20g %.20g %d\n", timep[id * Ninterpolate + countNodes], id, xp[id * Ninterpolate + countNodes], yp[id * Ninterpolate + countNodes], zp[id * Ninterpolate + countNodes], vxp[id * Ninterpolate + countNodes], vyp[id * Ninterpolate + countNodes], vzp[id * Ninterpolate + countNodes], id * Ninterpolate + countNodes);

			vxp[id * Ninterpolate + countNodes] /= dayUnit;
			vyp[id * Ninterpolate + countNodes] /= dayUnit;
			vzp[id * Ninterpolate + countNodes] /= dayUnit;

			if(useHelio == 1 && i == 0){
				vxpbc[countNodes] /= dayUnit;
				vypbc[countNodes] /= dayUnit;
				vzpbc[countNodes] /= dayUnit;
			}

			if(i == 0 && t == 0 && timep[id * Ninterpolate + countNodes] > time - (Ninterpolate/2 - 1) * dtime){
				printf("Error, time too small, not enough data before time\n");
				return 0;
			}
			if(i == Nperturbers - 1 && timep[id * Ninterpolate + countNodes] > time - Ninterpolate/2 * dtime){
				++countNodes;
			}
		}
		if(er < 0) break;
		if(countNodes >= Ninterpolate){
			break;
		}
	}
	if(countNodes < Ninterpolate){
		printf("Error, time too large, not enough data after time\n");
		return 0;
	}
	// ---------------------------------------



	// ---------------------------------------
	//interpolate
	interpolate(Nperturbers, Ninterpolate, xp, yp, zp, vxp, vyp, vzp, timep, time, x, y, z, vx, vy, vz);
	if(useHelio == 1){
		interpolate(1, Ninterpolate, xpbc, ypbc, zpbc, vxpbc, vypbc, vzpbc, timep, time, xbc, ybc, zbc, vxbc, vybc, vzbc);
	}
	// ---------------------------------------




	//first output
	double comx = 0.0;
	double comy = 0.0;
	double comz = 0.0;
	double vcomx = 0.0;
	double vcomy = 0.0;
	double vcomz = 0.0;
	double mtot = 0.0;	


	if(useHelio == 0 && outHelio == 1){
		//convert to heliocentric output

		comx = -x[0];
		comy = -y[0];
		comz = -z[0];
		vcomx = -vx[0];
		vcomy = -vy[0];
		vcomz = -vz[0];
	}

	if(useHelio == 1 && outHelio == 0){	
		//convert to barycentric output
/*
		for(int i = 0; i < N; ++i){
			comx += m[i] * x[i];
			comy += m[i] * y[i];
			comz += m[i] * z[i];
			vcomx += m[i] * vx[i];
			vcomy += m[i] * vy[i];
			vcomz += m[i] * vz[i];
			mtot += m[i];
		}
		comx /= mtot;
		comy /= mtot;
		comz /= mtot;
		vcomx /= mtot;
		vcomy /= mtot;
		vcomz /= mtot;

		comx = -comx;
		comy = -comy;
		comz = -comz;
		vcomx = -vcomx;
		vcomy = -vcomy;
		vcomz = -vcomz;
*/
		comx = -xbc[0];
		comy = -ybc[0];
		comz = -zbc[0];
		vcomx = -vxbc[0];
		vcomy = -vybc[0];
		vcomz = -vzbc[0];
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
		fprintf(outfile, "%.10g %d %.40g %.40g %.40g %.40g %.40g %.40g %.40g\n", time, i, m[i], comx + x[i], comy + y[i], comz + z[i], vcomx + vx[i], vcomy + vy[i], vcomz + vz[i]);

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

	//Barycentric coordinates
	double xb[N];
	double yb[N];
	double zb[N];
	double vxb[N];
	double vyb[N];
	double vzb[N];


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

	double c2 = 0.25;
	double c3 = 3.0 / 8.0;
	double c4 = 12.0 / 13.0;
	double c5 = 1.0;
	double c6 = 0.5;

	for(long long int t = 1; t <= Nsteps; ++t){
		
		//stage 1	
		for(int i = 0; i < N; ++i){
			xt[i] = x[i];
			yt[i] = y[i];
			zt[i] = z[i];
			vxt[i] = vx[i];
			vyt[i] = vy[i];
			vzt[i] = vz[i];
		}

		//Calculate Newtonian accelerations
		Newtonian(m, xt, yt, zt, aNx, aNy, aNz, aNxb, aNyb, aNzb, N, Nperturbers, useHelio);

		//convert to barycentric coordinates
		if(useHelio == 1 && GR > 0){	
/*
			comx = 0.0;
			comy = 0.0;
			comz = 0.0;
			vcomx = 0.0;
			vcomy = 0.0;
			vcomz = 0.0;
			mtot = 0.0;
			for(int i = 0; i < N; ++i){
				comx += m[i] * xt[i];
				comy += m[i] * yt[i];
				comz += m[i] * zt[i];
				vcomx += m[i] * vxt[i];
				vcomy += m[i] * vyt[i];
				vcomz += m[i] * vzt[i];
				mtot += m[i];
			}
			comx /= mtot;
			comy /= mtot;
			comz /= mtot;
			vcomx /= mtot;
			vcomy /= mtot;
			vcomz /= mtot;
		
			for(int i = 0; i < N; ++i){
				xb[i] = xt[i] - comx;
				yb[i] = yt[i] - comy;
				zb[i] = zt[i] - comz;
				vxb[i] = vxt[i] - vcomx;
				vyb[i] = vyt[i] - vcomy;
				vzb[i] = vzt[i] - vcomz;
			}
*/
			for(int i = 0; i < N; ++i){
				xb[i] = xt[i] - xbc[0];
				yb[i] = yt[i] - ybc[0];
				zb[i] = zt[i] - zbc[0];
				vxb[i] = vxt[i] - vxbc[0];
				vyb[i] = vyt[i] - vybc[0];
				vzb[i] = vzt[i] - vzbc[0];
			}
		}

		for(int i = Nperturbers; i < N; ++i){
			kx[i][0] = vxt[i];
			ky[i][0] = vyt[i];
			kz[i][0] = vzt[i];

			double ax = aNx[i];
			double ay = aNy[i];
			double az = aNz[i];

			if(useHelio == 1 && GR == 4){
				GRCall(m, xb, yb, zb, vxb, vyb, vzb, ax, ay, az, aNxb, aNyb, aNzb, N, Nperturbers, useHelio, GR, i);
			}
			else{
				GRCall(m, xt, yt, zt, vxt, vyt, vzt, ax, ay, az, aNx, aNy, aNz, N, Nperturbers, useHelio, GR, i);
			}
<<<<<<< HEAD
			//NonGrav(xt, yt, zt, vxt, vyt, vzt, ax, ay, az, A1, A2, A3, ALN, NK, NM, Nn, R0, i);
			//J2(m, xt, yt, zt, ax, ay, az, i);
=======
			NonGrav(xt, yt, zt, vxt, vyt, vzt, ax, ay, az, A1, A2, A3, ALN, NK, NM, Nn, R0, i);
			J2(m, xt, yt, zt, ax, ay, az, i);
>>>>>>> f28efb018735bbd62b2845e662c43faea3666aee

			kvx[i][0] = ax;
			kvy[i][0] = ay;
			kvz[i][0] = az;
		}
		
		//stage 2
		interpolate(Nperturbers, Ninterpolate, xp, yp, zp, vxp, vyp, vzp, timep, time + c2 * dt / dayUnit, xt, yt, zt, vxt, vyt, vzt);
		if(useHelio == 1){
			interpolate(1, Ninterpolate, xpbc, ypbc, zpbc, vxpbc, vypbc, vzpbc, timep, time + c2 * dt / dayUnit, xbc, ybc, zbc, vxbc, vybc, vzbc);
		}
		for(int i = Nperturbers; i < N; ++i){
			xt[i]  = x[i]  + dt * a21 * kx[i][0];
			yt[i]  = y[i]  + dt * a21 * ky[i][0];
			zt[i]  = z[i]  + dt * a21 * kz[i][0];
			vxt[i] = vx[i] + dt * a21 * kvx[i][0];
			vyt[i] = vy[i] + dt * a21 * kvy[i][0];
			vzt[i] = vz[i] + dt * a21 * kvz[i][0];
		}
		//Calculate Newtonian accelerations
		Newtonian(m, xt, yt, zt, aNx, aNy, aNz, aNxb, aNyb, aNzb, N, Nperturbers, useHelio);

		//convert to barycentric coordinates
		if(useHelio == 1 && GR > 0){	
/*
			comx = 0.0;
			comy = 0.0;
			comz = 0.0;
			vcomx = 0.0;
			vcomy = 0.0;
			vcomz = 0.0;
			mtot = 0.0;
			for(int i = 0; i < N; ++i){
				comx += m[i] * xt[i];
				comy += m[i] * yt[i];
				comz += m[i] * zt[i];
				vcomx += m[i] * vxt[i];
				vcomy += m[i] * vyt[i];
				vcomz += m[i] * vzt[i];
				mtot += m[i];
			}
			comx /= mtot;
			comy /= mtot;
			comz /= mtot;
			vcomx /= mtot;
			vcomy /= mtot;
			vcomz /= mtot;
		
			for(int i = 0; i < N; ++i){
				xb[i] = xt[i] - comx;
				yb[i] = yt[i] - comy;
				zb[i] = zt[i] - comz;
				vxb[i] = vxt[i] - vcomx;
				vyb[i] = vyt[i] - vcomy;
				vzb[i] = vzt[i] - vcomz;
			}
*/
			for(int i = 0; i < N; ++i){
				xb[i] = xt[i] - xbc[0];
				yb[i] = yt[i] - ybc[0];
				zb[i] = zt[i] - zbc[0];
				vxb[i] = vxt[i] - vxbc[0];
				vyb[i] = vyt[i] - vybc[0];
				vzb[i] = vzt[i] - vzbc[0];
			}
		}

		for(int i = Nperturbers; i < N; ++i){

			kx[i][1] = vxt[i];
			ky[i][1] = vyt[i];
			kz[i][1] = vzt[i];

			double ax = aNx[i];
			double ay = aNy[i];
			double az = aNz[i];

			if(useHelio == 1 && GR == 4){
				GRCall(m, xb, yb, zb, vxb, vyb, vzb, ax, ay, az, aNxb, aNyb, aNzb, N, Nperturbers, useHelio, GR, i);
			}
			else{
				GRCall(m, xt, yt, zt, vxt, vyt, vzt, ax, ay, az, aNx, aNy, aNz, N, Nperturbers, useHelio, GR, i);
			}
<<<<<<< HEAD
			//NonGrav(xt, yt, zt, vxt, vyt, vzt, ax, ay, az, A1, A2, A3, ALN, NK, NM, Nn, R0, i);
			//J2(m, xt, yt, zt, ax, ay, az, i);
=======
			NonGrav(xt, yt, zt, vxt, vyt, vzt, ax, ay, az, A1, A2, A3, ALN, NK, NM, Nn, R0, i);
			J2(m, xt, yt, zt, ax, ay, az, i);
>>>>>>> f28efb018735bbd62b2845e662c43faea3666aee

			kvx[i][1] = ax;
			kvy[i][1] = ay;
			kvz[i][1] = az;

		}
		//stage 3
		interpolate(Nperturbers, Ninterpolate, xp, yp, zp, vxp, vyp, vzp, timep, time + c3 * dt / dayUnit, xt, yt, zt, vxt, vyt, vzt);
		if(useHelio == 1){
			interpolate(1, Ninterpolate, xpbc, ypbc, zpbc, vxpbc, vypbc, vzpbc, timep, time + c3 * dt / dayUnit, xbc, ybc, zbc, vxbc, vybc, vzbc);
		}
		for(int i = Nperturbers; i < N; ++i){
			xt[i]  = x[i]  + dt * (a31 * kx[i][0]  + a32 * kx[i][1]);
			yt[i]  = y[i]  + dt * (a31 * ky[i][0]  + a32 * ky[i][1]);
			zt[i]  = z[i]  + dt * (a31 * kz[i][0]  + a32 * kz[i][1]);
			vxt[i] = vx[i] + dt * (a31 * kvx[i][0] + a32 * kvx[i][1]);
			vyt[i] = vy[i] + dt * (a31 * kvy[i][0] + a32 * kvy[i][1]);
			vzt[i] = vz[i] + dt * (a31 * kvz[i][0] + a32 * kvz[i][1]);
		}
		//Calculate Newtonian accelerations
		Newtonian(m, xt, yt, zt, aNx, aNy, aNz, aNxb, aNyb, aNzb, N, Nperturbers, useHelio);

		//convert to barycentric coordinates
		if(useHelio == 1 && GR > 0){	
/*
			comx = 0.0;
			comy = 0.0;
			comz = 0.0;
			vcomx = 0.0;
			vcomy = 0.0;
			vcomz = 0.0;
			mtot = 0.0;
			for(int i = 0; i < N; ++i){
				comx += m[i] * xt[i];
				comy += m[i] * yt[i];
				comz += m[i] * zt[i];
				vcomx += m[i] * vxt[i];
				vcomy += m[i] * vyt[i];
				vcomz += m[i] * vzt[i];
				mtot += m[i];
			}
			comx /= mtot;
			comy /= mtot;
			comz /= mtot;
			vcomx /= mtot;
			vcomy /= mtot;
			vcomz /= mtot;
		
			for(int i = 0; i < N; ++i){
				xb[i] = xt[i] - comx;
				yb[i] = yt[i] - comy;
				zb[i] = zt[i] - comz;
				vxb[i] = vxt[i] - vcomx;
				vyb[i] = vyt[i] - vcomy;
				vzb[i] = vzt[i] - vcomz;
			}
*/
			for(int i = 0; i < N; ++i){
				xb[i] = xt[i] - xbc[0];
				yb[i] = yt[i] - ybc[0];
				zb[i] = zt[i] - zbc[0];
				vxb[i] = vxt[i] - vxbc[0];
				vyb[i] = vyt[i] - vybc[0];
				vzb[i] = vzt[i] - vzbc[0];
			}
		}

		for(int i = Nperturbers; i < N; ++i){

			kx[i][2] = vxt[i];
			ky[i][2] = vyt[i];
			kz[i][2] = vzt[i];

			double ax = aNx[i];
			double ay = aNy[i];
			double az = aNz[i];

			if(useHelio == 1 && GR == 4){
				GRCall(m, xb, yb, zb, vxb, vyb, vzb, ax, ay, az, aNxb, aNyb, aNzb, N, Nperturbers, useHelio, GR, i);
			}
			else{
				GRCall(m, xt, yt, zt, vxt, vyt, vzt, ax, ay, az, aNx, aNy, aNz, N, Nperturbers, useHelio, GR, i);
			}
<<<<<<< HEAD
			//NonGrav(xt, yt, zt, vxt, vyt, vzt, ax, ay, az, A1, A2, A3, ALN, NK, NM, Nn, R0, i);
			//J2(m, xt, yt, zt, ax, ay, az, i);
=======
			NonGrav(xt, yt, zt, vxt, vyt, vzt, ax, ay, az, A1, A2, A3, ALN, NK, NM, Nn, R0, i);
			J2(m, xt, yt, zt, ax, ay, az, i);
>>>>>>> f28efb018735bbd62b2845e662c43faea3666aee

			kvx[i][2] = ax;
			kvy[i][2] = ay;
			kvz[i][2] = az;

		}
		//stage 4
		interpolate(Nperturbers, Ninterpolate, xp, yp, zp, vxp, vyp, vzp, timep, time + c4 * dt / dayUnit, xt, yt, zt, vxt, vyt, vzt);
		if(useHelio == 1){
			interpolate(1, Ninterpolate, xpbc, ypbc, zpbc, vxpbc, vypbc, vzpbc, timep, time + c4 * dt / dayUnit, xbc, ybc, zbc, vxbc, vybc, vzbc);
		}
		for(int i = Nperturbers; i < N; ++i){
			xt[i]  = x[i]  + dt * (a41 * kx[i][0]  + a42 * kx[i][1]  + a43 * kx[i][2]);
			yt[i]  = y[i]  + dt * (a41 * ky[i][0]  + a42 * ky[i][1]  + a43 * ky[i][2]);
			zt[i]  = z[i]  + dt * (a41 * kz[i][0]  + a42 * kz[i][1]  + a43 * kz[i][2]);
			vxt[i] = vx[i] + dt * (a41 * kvx[i][0] + a42 * kvx[i][1] + a43 * kvx[i][2]);
			vyt[i] = vy[i] + dt * (a41 * kvy[i][0] + a42 * kvy[i][1] + a43 * kvy[i][2]);
			vzt[i] = vz[i] + dt * (a41 * kvz[i][0] + a42 * kvz[i][1] + a43 * kvz[i][2]);
		}
		//Calculate Newtonian accelerations
		Newtonian(m, xt, yt, zt, aNx, aNy, aNz, aNxb, aNyb, aNzb, N, Nperturbers, useHelio);

		//convert to barycentric coordinates
		if(useHelio == 1 && GR > 0){	
/*
			comx = 0.0;
			comy = 0.0;
			comz = 0.0;
			vcomx = 0.0;
			vcomy = 0.0;
			vcomz = 0.0;
			mtot = 0.0;
			for(int i = 0; i < N; ++i){
				comx += m[i] * xt[i];
				comy += m[i] * yt[i];
				comz += m[i] * zt[i];
				vcomx += m[i] * vxt[i];
				vcomy += m[i] * vyt[i];
				vcomz += m[i] * vzt[i];
				mtot += m[i];
			}
			comx /= mtot;
			comy /= mtot;
			comz /= mtot;
			vcomx /= mtot;
			vcomy /= mtot;
			vcomz /= mtot;
		
			for(int i = 0; i < N; ++i){
				xb[i] = xt[i] - comx;
				yb[i] = yt[i] - comy;
				zb[i] = zt[i] - comz;
				vxb[i] = vxt[i] - vcomx;
				vyb[i] = vyt[i] - vcomy;
				vzb[i] = vzt[i] - vcomz;
			}
*/
			for(int i = 0; i < N; ++i){
				xb[i] = xt[i] - xbc[0];
				yb[i] = yt[i] - ybc[0];
				zb[i] = zt[i] - zbc[0];
				vxb[i] = vxt[i] - vxbc[0];
				vyb[i] = vyt[i] - vybc[0];
				vzb[i] = vzt[i] - vzbc[0];
			}
		}

		for(int i = Nperturbers; i < N; ++i){

			kx[i][3] = vxt[i];
			ky[i][3] = vyt[i];
			kz[i][3] = vzt[i];

			double ax = aNx[i];
			double ay = aNy[i];
			double az = aNz[i];

			if(useHelio == 1 && GR == 4){
				GRCall(m, xb, yb, zb, vxb, vyb, vzb, ax, ay, az, aNxb, aNyb, aNzb, N, Nperturbers, useHelio, GR, i);
			}
			else{
				GRCall(m, xt, yt, zt, vxt, vyt, vzt, ax, ay, az, aNx, aNy, aNz, N, Nperturbers, useHelio, GR, i);
			}
<<<<<<< HEAD
			//NonGrav(xt, yt, zt, vxt, vyt, vzt, ax, ay, az, A1, A2, A3, ALN, NK, NM, Nn, R0, i);
			//J2(m, xt, yt, zt, ax, ay, az, i);
=======
			NonGrav(xt, yt, zt, vxt, vyt, vzt, ax, ay, az, A1, A2, A3, ALN, NK, NM, Nn, R0, i);
			J2(m, xt, yt, zt, ax, ay, az, i);
>>>>>>> f28efb018735bbd62b2845e662c43faea3666aee

			kvx[i][3] = ax;
			kvy[i][3] = ay;
			kvz[i][3] = az;

		}
		//stage 5
		interpolate(Nperturbers, Ninterpolate, xp, yp, zp, vxp, vyp, vzp, timep, time + c5 * dt / dayUnit, xt, yt, zt, vxt, vyt, vzt);
		if(useHelio == 1){
			interpolate(1, Ninterpolate, xpbc, ypbc, zpbc, vxpbc, vypbc, vzpbc, timep, time + c5 * dt / dayUnit, xbc, ybc, zbc, vxbc, vybc, vzbc);
		}
		for(int i = Nperturbers; i < N; ++i){
			xt[i]  = x[i]  + dt * (a51 * kx[i][0]  + a52 * kx[i][1]  + a53 * kx[i][2]  + a54 * kx[i][3]);
			yt[i]  = y[i]  + dt * (a51 * ky[i][0]  + a52 * ky[i][1]  + a53 * ky[i][2]  + a54 * ky[i][3]);
			zt[i]  = z[i]  + dt * (a51 * kz[i][0]  + a52 * kz[i][1]  + a53 * kz[i][2]  + a54 * kz[i][3]);
			vxt[i] = vx[i] + dt * (a51 * kvx[i][0] + a52 * kvx[i][1] + a53 * kvx[i][2] + a54 * kvx[i][3]);
			vyt[i] = vy[i] + dt * (a51 * kvy[i][0] + a52 * kvy[i][1] + a53 * kvy[i][2] + a54 * kvy[i][3]);
			vzt[i] = vz[i] + dt * (a51 * kvz[i][0] + a52 * kvz[i][1] + a53 * kvz[i][2] + a54 * kvz[i][3]);
		}
		//Calculate Newtonian accelerations
		Newtonian(m, xt, yt, zt, aNx, aNy, aNz, aNxb, aNyb, aNzb, N, Nperturbers, useHelio);

		//convert to barycentric coordinates
		if(useHelio == 1 && GR > 0){	
/*
			comx = 0.0;
			comy = 0.0;
			comz = 0.0;
			vcomx = 0.0;
			vcomy = 0.0;
			vcomz = 0.0;
			mtot = 0.0;
			for(int i = 0; i < N; ++i){
				comx += m[i] * xt[i];
				comy += m[i] * yt[i];
				comz += m[i] * zt[i];
				vcomx += m[i] * vxt[i];
				vcomy += m[i] * vyt[i];
				vcomz += m[i] * vzt[i];
				mtot += m[i];
			}
			comx /= mtot;
			comy /= mtot;
			comz /= mtot;
			vcomx /= mtot;
			vcomy /= mtot;
			vcomz /= mtot;
		
			for(int i = 0; i < N; ++i){
				xb[i] = xt[i] - comx;
				yb[i] = yt[i] - comy;
				zb[i] = zt[i] - comz;
				vxb[i] = vxt[i] - vcomx;
				vyb[i] = vyt[i] - vcomy;
				vzb[i] = vzt[i] - vcomz;
			}
*/
			for(int i = 0; i < N; ++i){
				xb[i] = xt[i] - xbc[0];
				yb[i] = yt[i] - ybc[0];
				zb[i] = zt[i] - zbc[0];
				vxb[i] = vxt[i] - vxbc[0];
				vyb[i] = vyt[i] - vybc[0];
				vzb[i] = vzt[i] - vzbc[0];
			}
		}

		for(int i = Nperturbers; i < N; ++i){

			kx[i][4] = vxt[i];
			ky[i][4] = vyt[i];
			kz[i][4] = vzt[i];

			double ax = aNx[i];
			double ay = aNy[i];
			double az = aNz[i];

			if(useHelio == 1 && GR == 4){
				GRCall(m, xb, yb, zb, vxb, vyb, vzb, ax, ay, az, aNxb, aNyb, aNzb, N, Nperturbers, useHelio, GR, i);
			}
			else{
				GRCall(m, xt, yt, zt, vxt, vyt, vzt, ax, ay, az, aNx, aNy, aNz, N, Nperturbers, useHelio, GR, i);
			}
<<<<<<< HEAD
			//NonGrav(xt, yt, zt, vxt, vyt, vzt, ax, ay, az, A1, A2, A3, ALN, NK, NM, Nn, R0, i);
			//J2(m, xt, yt, zt, ax, ay, az, i);
=======
			NonGrav(xt, yt, zt, vxt, vyt, vzt, ax, ay, az, A1, A2, A3, ALN, NK, NM, Nn, R0, i);
			J2(m, xt, yt, zt, ax, ay, az, i);
>>>>>>> f28efb018735bbd62b2845e662c43faea3666aee

			kvx[i][4] = ax;
			kvy[i][4] = ay;
			kvz[i][4] = az;

		}

		//stage 6
		interpolate(Nperturbers, Ninterpolate, xp, yp, zp, vxp, vyp, vzp, timep, time + c6 * dt / dayUnit, xt, yt, zt, vxt, vyt, vzt);
		if(useHelio == 1){
			interpolate(1, Ninterpolate, xpbc, ypbc, zpbc, vxpbc, vypbc, vzpbc, timep, time + c6 * dt / dayUnit, xbc, ybc, zbc, vxbc, vybc, vzbc);
		}
		for(int i = Nperturbers; i < N; ++i){
			xt[i]  = x[i]  + dt * (a61 * kx[i][0]  + a62 * kx[i][1]  + a63 * kx[i][2]  + a64 * kx[i][3]  + a65 * kx[i][4]);
			yt[i]  = y[i]  + dt * (a61 * ky[i][0]  + a62 * ky[i][1]  + a63 * ky[i][2]  + a64 * ky[i][3]  + a65 * ky[i][4]);
			zt[i]  = z[i]  + dt * (a61 * kz[i][0]  + a62 * kz[i][1]  + a63 * kz[i][2]  + a64 * kz[i][3]  + a65 * kz[i][4]);
			vxt[i] = vx[i] + dt * (a61 * kvx[i][0] + a62 * kvx[i][1] + a63 * kvx[i][2] + a64 * kvx[i][3] + a65 * kvx[i][4]);
			vyt[i] = vy[i] + dt * (a61 * kvy[i][0] + a62 * kvy[i][1] + a63 * kvy[i][2] + a64 * kvy[i][3] + a65 * kvy[i][4]);
			vzt[i] = vz[i] + dt * (a61 * kvz[i][0] + a62 * kvz[i][1] + a63 * kvz[i][2] + a64 * kvz[i][3] + a65 * kvz[i][4]);
		}

		//Calculate Newtonian accelerations
		Newtonian(m, xt, yt, zt, aNx, aNy, aNz, aNxb, aNyb, aNzb, N, Nperturbers, useHelio);

		//convert to barycentric coordinates
		if(useHelio == 1 && GR > 0){	
/*
			comx = 0.0;
			comy = 0.0;
			comz = 0.0;
			vcomx = 0.0;
			vcomy = 0.0;
			vcomz = 0.0;
			mtot = 0.0;
			for(int i = 0; i < N; ++i){
				comx += m[i] * xt[i];
				comy += m[i] * yt[i];
				comz += m[i] * zt[i];
				vcomx += m[i] * vxt[i];
				vcomy += m[i] * vyt[i];
				vcomz += m[i] * vzt[i];
				mtot += m[i];
			}
			comx /= mtot;
			comy /= mtot;
			comz /= mtot;
			vcomx /= mtot;
			vcomy /= mtot;
			vcomz /= mtot;
		
			for(int i = 0; i < N; ++i){
				xb[i] = xt[i] - comx;
				yb[i] = yt[i] - comy;
				zb[i] = zt[i] - comz;
				vxb[i] = vxt[i] - vcomx;
				vyb[i] = vyt[i] - vcomy;
				vzb[i] = vzt[i] - vcomz;
			}
*/
			for(int i = 0; i < N; ++i){
				xb[i] = xt[i] - xbc[0];
				yb[i] = yt[i] - ybc[0];
				zb[i] = zt[i] - zbc[0];
				vxb[i] = vxt[i] - vxbc[0];
				vyb[i] = vyt[i] - vybc[0];
				vzb[i] = vzt[i] - vzbc[0];
			}
		}

		for(int i = Nperturbers; i < N; ++i){

			kx[i][5] = vxt[i];
			ky[i][5] = vyt[i];
			kz[i][5] = vzt[i];

			double ax = aNx[i];
			double ay = aNy[i];
			double az = aNz[i];

			if(useHelio == 1 && GR == 4){
				GRCall(m, xb, yb, zb, vxb, vyb, vzb, ax, ay, az, aNxb, aNyb, aNzb, N, Nperturbers, useHelio, GR, i);
			}
			else{
				GRCall(m, xt, yt, zt, vxt, vyt, vzt, ax, ay, az, aNx, aNy, aNz, N, Nperturbers, useHelio, GR, i);
			}
<<<<<<< HEAD
			//NonGrav(xt, yt, zt, vxt, vyt, vzt, ax, ay, az, A1, A2, A3, ALN, NK, NM, Nn, R0, i);
			//J2(m, xt, yt, zt, ax, ay, az, i);
=======
			NonGrav(xt, yt, zt, vxt, vyt, vzt, ax, ay, az, A1, A2, A3, ALN, NK, NM, Nn, R0, i);
			J2(m, xt, yt, zt, ax, ay, az, i);
>>>>>>> f28efb018735bbd62b2845e662c43faea3666aee

			kvx[i][5] = ax;
			kvy[i][5] = ay;
			kvz[i][5] = az;

		}
		double sc = 1.0e-15;

		/*	
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
		*/

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

		/*
		double ee = 1.0/5.0;	

		double s = pow( 1.0  / error, ee);
		*/

		//printf("%g %g\n", dt, s);			

		for(int i = Nperturbers; i < N; ++i){
		
			//RKF45
			x[i] += dt * (b1 * kx[i][0] + b3 * kx[i][2] + b4 * kx[i][3] + b5 * kx[i][4] + b6 * kx[i][5]);
			y[i] += dt * (b1 * ky[i][0] + b3 * ky[i][2] + b4 * ky[i][3] + b5 * ky[i][4] + b6 * ky[i][5]);
			z[i] += dt * (b1 * kz[i][0] + b3 * kz[i][2] + b4 * kz[i][3] + b5 * kz[i][4] + b6 * kz[i][5]);

			vx[i] += dt * (b1 * kvx[i][0] + b3 * kvx[i][2] + b4 * kvx[i][3] + b5 * kvx[i][4] + b6 * kvx[i][5]);
			vy[i] += dt * (b1 * kvy[i][0] + b3 * kvy[i][2] + b4 * kvy[i][3] + b5 * kvy[i][4] + b6 * kvy[i][5]);
			vz[i] += dt * (b1 * kvz[i][0] + b3 * kvz[i][2] + b4 * kvz[i][3] + b5 * kvz[i][4] + b6 * kvz[i][5]);
		}

		time = time0 + t * dt / dayUnit;

		//update table
		if(time - time1 >= dtime){
			int countNodes = Ninterpolate - 1;
			int er;
			for(int j = 0; j < Ninterpolate - 1; ++j){
				for(int i = 0; i < Nperturbers; ++i){
					xp[i * Ninterpolate + j] = xp[i * Ninterpolate + j + 1];
					yp[i * Ninterpolate + j] = yp[i * Ninterpolate + j + 1];
					zp[i * Ninterpolate + j] = zp[i * Ninterpolate + j + 1];
					timep[i * Ninterpolate + j] = timep[i * Ninterpolate + j + 1];
				}
			}

//printf("CountNodes %d\n", countNodes);
			for(int i = 0; i < Nperturbers; ++i){
				double skip;
				double timepp;
				int id;
				er = fscanf(XVfile, "%lf %d", &timepp, &id);
				fscanf(XVfile, "%lf %lf %lf", &xp[id * Ninterpolate + countNodes], &yp[id * Ninterpolate + countNodes], &zp[id * Ninterpolate + countNodes]);
				fscanf(XVfile, "%lf %lf %lf", &vxp[id * Ninterpolate + countNodes], &vyp[id * Ninterpolate + countNodes], &vzp[id * Ninterpolate + countNodes]);
				if(useHelio == 1 && i == 0){
					double bcTime;
					fscanf(BCfile, "%lf", &bcTime);
					fscanf(BCfile, "%lf %lf %lf", &xpbc[countNodes], &ypbc[countNodes], &zpbc[countNodes]);
					fscanf(BCfile, "%lf %lf %lf", &vxpbc[countNodes], &vypbc[countNodes], &vzpbc[countNodes]);
				}

				vxp[id * Ninterpolate + countNodes] /= dayUnit;
				vyp[id * Ninterpolate + countNodes] /= dayUnit;
				vzp[id * Ninterpolate + countNodes] /= dayUnit;

				if(useHelio == 1 && i == 0){
					vxpbc[countNodes] /= dayUnit;
					vypbc[countNodes] /= dayUnit;
					vzpbc[countNodes] /= dayUnit;
				}
				if(er < 0) break;
				timep[id * Ninterpolate + countNodes] = timepp;
	
//printf("%.20g %d %.20g %.20g %.20g %d\n", timep[id * Ninterpolate + countNodes], id, xp[id * Ninterpolate + countNodes], yp[id * Ninterpolate + countNodes], zp[id * Ninterpolate + countNodes], id * Ninterpolate + countNodes);

			}
			if(er < 0){
				printf("Error, time too large, not enough data after time\n");
				return 0;
			}
			time1 = time;
		}

		// ---------------------------------------
		//interpolate
		interpolate(Nperturbers, Ninterpolate, xp, yp, zp, vxp, vyp, vzp, timep, time, x, y, z, vx, vy, vz);
		if(useHelio == 1){
			interpolate(1, Ninterpolate, xpbc, ypbc, zpbc, vxpbc, vypbc, vzpbc, timep, time, xbc, ybc, zbc, vxbc, vybc, vzbc);
		}
		// ---------------------------------------


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
			mtot = 0.0;

			if(useHelio == 0 && outHelio == 1){
				//convert to heliocentric output
				comx = -x[0];
				comy = -y[0];
				comz = -z[0];
				vcomx = -vx[0];
				vcomy = -vy[0];
				vcomz = -vz[0];
			}
			if(useHelio == 1 && outHelio == 0){	
/*
				//convert to barycentric output
				for(int i = 0; i < N; ++i){
					comx += m[i] * x[i];
					comy += m[i] * y[i];
					comz += m[i] * z[i];
					vcomx += m[i] * vx[i];
					vcomy += m[i] * vy[i];
					vcomz += m[i] * vz[i];
					mtot += m[i];
				}
				comx /= mtot;
				comy /= mtot;
				comz /= mtot;
				vcomx /= mtot;
				vcomy /= mtot;
				vcomz /= mtot;

				comx = -comx;
				comy = -comy;
				comz = -comz;
				vcomx = -vcomx;
				vcomy = -vcomy;
				vcomz = -vcomz;
*/
				comx = -xbc[0];
				comy = -ybc[0];
				comz = -zbc[0];
				vcomx = -vxbc[0];
				vcomy = -vybc[0];
				vcomz = -vzbc[0];
			}

			
			for(int i = 0; i < N; ++i){
				fprintf(outfile, "%.10g %d %.40g %.40g %.40g %.40g %.40g %.40g %.40g\n", time, i, m[i], comx + x[i], comy + y[i], comz + z[i], vcomx + vx[i], vcomy + vy[i], vcomz + vz[i]);
				
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
	fclose(XVfile);
	if(useHelio == 1){
		fclose(BCfile);
	}
}
	

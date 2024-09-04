//Non Grav
inline void asteroid::NonGrav(double *x, double *y, double *z, double *vx, double *vy, double *vz){

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

		double f1 = A1_h[i] * gr / r;
		double f2 = A2_h[i] * gr / t;
		double f3 = A3_h[i] * gr / h;

//printf("NonGrav  %.20g %.20g %.20g %.20g |%.20g %.20g %.20g\n", gr, r, t, h, f1, f2, f3);
//printf("NonGrav a %.20g %.20g %.20g\n", ax_h[i], ay_h[i], az_h[i]);

		//ax_h[i] += f1 * xi + f2 * tx + f3 * hx;
		//ay_h[i] += f1 * yi + f2 * ty + f3 * hy;
		//az_h[i] += f1 * zi + f2 * tz + f3 * hz;
		ax_h[i] += A1_h[i] * gr * xi / r + A2_h[i] * gr * tx / t + A3_h[i] * gr * hx / h;
		ay_h[i] += A1_h[i] * gr * yi / r + A2_h[i] * gr * ty / t + A3_h[i] * gr * hy / h;
		az_h[i] += A1_h[i] * gr * zi / r + A2_h[i] * gr * tz / t + A3_h[i] * gr * hz / h;

//printf("NonGrav %.20g %.20g %.20g %.20g | %.20g %.20g %.20g\n", gr, r, t, h, f1, f2, f3);
//printf("NonGrav a %.20g %.20g %.20g\n", ax_h[i], ay_h[i], az_h[i]);

	}
}
inline void asteroid::GR(double *x, double *y, double *z, double *vx, double *vy, double *vz){
	//GR
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
		double vsq = vxi * vxi + vyi * vyi + vzi * vzi;

		double rv = xi * vxi + yi * vyi + zi * vzi;

		double f1 = GM[10] / (r * r * r * c2);
		double t1 = 4.0 * GM[10] / r - vsq;
		double t3 = 4.0 * rv;
		//printf("GRa %.20g %.20g %.20g %.20g %.20g\n", vsq, r, t1, t3, f1);
		//printf("a %d %.20g %.20g %.20g\n", i, ax_h, ay_h, az_h);

		//printf("A %d %.20g %.20g %.20g %.20g %.20g %.20g\n", i, xi, yi, zi, vxi, vyi, vzi);
		//printf("B %d %.20g %.20g %.20g\n", i, f1, t1, t3);

		double aax = f1 * (t1 * xi + t3 * vxi);
		double aay = f1 * (t1 * yi + t3 * vyi);
		double aaz = f1 * (t1 * zi + t3 * vzi);

		//printf("GR a %.20g %.20g %.20g\n", aax, aay, aaz);
		ax_h[i] += f1 * (t1 * xi + t3 * vxi);
		ay_h[i] += f1 * (t1 * yi + t3 * vyi);
		az_h[i] += f1 * (t1 * zi + t3 * vzi);
		//printf("GR a %d %.20g %.20g %.20g\n", i, ax_h[i], ay_h[i], az_h[i]);
	}
}

inline void asteroid::J2(double *x, double *y, double *z){
	//J2
	for(int i = Nperturbers; i < N; ++i){

//printf("J2 %.20g %.20g %.20g\n", J2E, RE, REAU);

		double xE = x[i] - x[2];
		double yE = y[i] - y[2];
		double zE = z[i] - z[2];
//printf("dz %.20g %.20g %.20g\n", xE, yE, zE);

		double rsq = xE * xE + yE * yE + zE * zE;
		double r = sqrt(rsq);
		//double r5 = rsq * rsq * r;

		//double t1 = GM[2] * 3.0 * J2E * REAU * REAU / (2.0 * r5);
		double t1 = 3.0 * J2E * REAU * REAU / rsq / rsq / r / 2.0;
		double t2 = 5.0 * (zE * zE / rsq) - 1.0;

//printf("%.20g %.20g %.20g\n", GM[2], t1, t2);

		double tx = GM[2] * t1 * t2 * xE;
		double ty = GM[2] * t1 * t2 * yE;
		double tz = GM[2] * t1 * (t2 - 2.0) * zE;

		ax_h[i] += tx;
		ay_h[i] += ty;
		az_h[i] += tz;

//printf("J2 a %.20g %.20g %.20g\n", tx, ty, tz);
	}
}

inline void asteroid::Gravity(double *x, double *y, double *z){
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

			ax_h[i] -= s*dx;
			ay_h[i] -= s*dy;
			az_h[i] -= s*dz;

//printf("position %d %d %.20g %.20g %.20g %.20g %.20g\n", pp, p, time, x[p], y[p], z[p], GM[p]);
		}
	}
}

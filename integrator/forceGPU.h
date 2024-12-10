//Non Grav
//Coordinates must be heliocentric
__device__ void NonGrav(double xi, double yi, double zi, double vxi, double vyi, double vzi, double A1i, double A2i, double A3i, double r, double &axi, double &ayi, double &azi){

//printf("r %.20g %.20g %.20g %.20g\n", xi, yi ,zi, r);
//printf("v %.20g %.20g %.20g %.20g %.20g\n", vxi, vyi, vzi, vx[i], vx[10]);

	//angular momenrum h = r x v
	double hx =  __dmul_rn(yi, vzi) - __dmul_rn(zi, vyi);
	double hy = -__dmul_rn(xi, vzi) + __dmul_rn(zi, vxi);
	double hz =  __dmul_rn(xi, vyi) - __dmul_rn(yi, vxi);

	double hsq = __dmul_rn(hx, hx) + __dmul_rn(hy, hy) + __dmul_rn(hz, hz);
	double h = sqrt(hsq);

	//Transverse velocity t = h x r
	double tx =  __dmul_rn(hy, zi) - __dmul_rn(hz, yi);
	double ty = -__dmul_rn(hx, zi) + __dmul_rn(hz, xi);
	double tz =  __dmul_rn(hx, yi) - __dmul_rn(hy, xi);

	double tsq = __dmul_rn(tx, tx) + __dmul_rn(ty,ty) + __dmul_rn(tz, tz);
	double t = sqrt(tsq);

	double gr = 0.0;

	if(cometFlag_c == 0){
		//const double alpha = 1.0;
		//const double nk = 0.0;
		//const double nm = 2.0;
		//const double nn = 5.093;
		//const double r0 = 1.0;
		//gr = alpha * pow(r / r0, -nm) * pow(1.0 + pow(r / r0, nn), -nk);
		gr = 1.0 / (r * r);
		//double gr = 1.0 / rsq;  //only valid for asteroids, not for comets 
	}
	else{
		double rr = r / nonGrav_r0_c;
		r = nonGrav_alpha_c * pow(rr, -nonGrav_nm_c) * pow(1.0 + pow(rr, nonGrav_nn_c), -nonGrav_nk_c);
	}

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

	double f1 = __dmul_rn(A1i, gr) / r;
	double f2 = __dmul_rn(A2i, gr) / t;
	double f3 = __dmul_rn(A3i, gr) / h;

//printf("NonGrav  %.20g %.20g %.20g %.20g |%.20g %.20g %.20g\n", gr, r, t, h, f1, f2, f3);
//printf("NonGrav a %.20g %.20g %.20g\n", ax_h[i], ay_h[i], az_h[i]);

	axi += __dmul_rn(f1, xi) + __dmul_rn(f2, tx) + __dmul_rn(f3, hx);
	ayi += __dmul_rn(f1, yi) + __dmul_rn(f2, ty) + __dmul_rn(f3, hy);
	azi += __dmul_rn(f1, zi) + __dmul_rn(f2, tz) + __dmul_rn(f3, hz);
	//axi += A1i * gr * xi / r + A2i * gr * tx / t + A3i * gr * hx / h;
	//ayi += A1i * gr * yi / r + A2i * gr * ty / t + A3i * gr * hy / h;
	//azi += A1i * gr * zi / r + A2i * gr * tz / t + A3i * gr * hz / h;

//printf("NonGrav %.20g %.20g %.20g %.20g | %.20g %.20g %.20g\n", gr, r, t, h, f1, f2, f3);
//printf("NonGrav a %.20g %.20g %.20g\n", ax_h[i], ay_h[i], az_h[i]);

}

//GR
//Coordinates must be heliocentric
__device__ void GR(double xi, double yi, double zi, double vxi, double vyi, double vzi, double r, double &axi, double &ayi, double &azi, double GMSun){
	//GR
	double vsq = __dmul_rn(vxi, vxi) + __dmul_rn(vyi, vyi) + __dmul_rn(vzi, vzi);

	double rv = __dmul_rn(xi, vxi) + __dmul_rn(yi, vyi) + __dmul_rn(zi, vzi);

	double f1 = GMSun / (r * r * r * c2_c);
	double t1 = 4.0 * GMSun / r - vsq;
	double t3 = 4.0 * rv;
	//printf("GRa %.20g %.20g %.20g %.20g %.20g\n", vsq, r, t1, t3, f1);
	//printf("a %d %.20g %.20g %.20g\n", i, ax_h, ay_h, az_h);

	//printf("A %.20g %.20g %.20g %.20g %.20g %.20g\n", xi, yi, zi, vxi, vyi, vzi);
	//printf("B %d %.20g %.20g %.20g\n", i, f1, t1, t3);

	double aax = f1 * (__dmul_rn(t1, xi) + __dmul_rn(t3, vxi));
	double aay = f1 * (__dmul_rn(t1, yi) + __dmul_rn(t3, vyi));
	double aaz = f1 * (__dmul_rn(t1, zi) + __dmul_rn(t3, vzi));

	//printf("GR a %.20g %.20g %.20g\n", aax, aay, aaz);

	axi = __dadd_rn(axi, aax);
	ayi = __dadd_rn(ayi, aay);
	azi = __dadd_rn(azi, aaz);

	//printf("GR a %d %.20g %.20g %.20g\n", i, ax_h[i], ay_h[i], az_h[i]);
}

//Coordinates must be Earth centric
__device__ void J2(double xE, double yE, double zE, double &axi, double &ayi, double &azi, double GMEarth){
	//J2

//printf("J2 %.20g %.20g\n", J2E_c, REAU_c);

	double rsq = xE * xE + yE * yE + zE * zE;
	double r = sqrt(rsq);
	//double r5 = rsq * rsq * r;

	//double t1 = GMEarth * 3.0 * J2E * REAU_c * REAU_c / (2.0 * r5);
	double t1 = 3.0 * J2E_c * REAU_c * REAU_c / rsq / rsq / r / 2.0;
	double t2 = 5.0 * (zE * zE / rsq) - 1.0;

//printf("%.20g %.20g %.20g\n", GMEarth, t1, t2);

	double tx = GMEarth * t1 * t2 * xE;
	double ty = GMEarth * t1 * t2 * yE;
	double tz = GMEarth * t1 * (t2 - 2.0) * zE;

	axi += tx;
	ayi += ty;
	azi += tz;

//printf("J2 a %.20g %.20g %.20g\n", tx, ty, tz);
}

__device__ void Gravity(double xi, double yi, double zi, double *xTable_s, double *yTable_s, double *zTable_s, double &axi, double &ayi, double &azi, double *GM_s, const int pp){

	int p = pp + 11;
	if(pp == 16) p = 8;	//Pluto
	if(pp == 17) p = 9;	//Moon
	if(pp == 18) p = 3;	//Mars
	if(pp == 19) p = 0;	//Mercur
	if(pp == 20) p = 7;	//Neptune
	if(pp == 21) p = 6;	//Uranus
	if(pp == 22) p = 2;	//Earth
	if(pp == 23) p = 1;	//Venus
	if(pp == 24) p = 5;	//Saturn
	if(pp == 25) p = 4;	//Jupiter
	if(pp == 26) p = 10;	//Sun

	double dx = xi - xTable_s[p];
	double dy = yi - yTable_s[p];
	double dz = zi - zTable_s[p];

	double rsq = __dmul_rn(dx, dx) + __dmul_rn(dy, dy) + __dmul_rn(dz, dz);
	double r = sqrt(rsq);
	double s = GM_s[p] / (r * rsq);

	double aax = -__dmul_rn(s, dx);
	double aay = -__dmul_rn(s, dy);
	double aaz = -__dmul_rn(s, dz);

	axi += aax;
	ayi += aay;
	azi += aaz;

//printf("position %d %d %.20g %.20g %.20g %.20g %.20g\n", pp, p, time, x[p], y[p], z[p], GM_h[p]);
}

__global__ void Gravity_kernel(double xi, double yi, double zi, double *xTable_s, double *yTable_s, double *zTable_s, double &axi, double &ayi, double &azi, double *GM_s, const int Nperturbers){

	int pp = threadIdx.x;
	if(pp < Nperturbers){

		int p = pp + 11;
		if(pp == 16) p = 8;	//Pluto
		if(pp == 17) p = 9;	//Moon
		if(pp == 18) p = 3;	//Mars
		if(pp == 19) p = 0;	//Mercur
		if(pp == 20) p = 7;	//Neptune
		if(pp == 21) p = 6;	//Uranus
		if(pp == 22) p = 2;	//Earth
		if(pp == 23) p = 1;	//Venus
		if(pp == 24) p = 5;	//Saturn
		if(pp == 25) p = 4;	//Jupiter
		if(pp == 26) p = 10;	//Sun

		double dx = xi - xTable_s[p];
		double dy = yi - yTable_s[p];
		double dz = zi - zTable_s[p];

		double rsq = __dmul_rn(dx, dx) + __dmul_rn(dy, dy) + __dmul_rn(dz, dz);
		double r = sqrt(rsq);
		double s = GM_s[p] / (r * rsq);

		double aax = -__dmul_rn(s, dx);
		double aay = -__dmul_rn(s, dy);
		double aaz = -__dmul_rn(s, dz);

		axi += aax;
		ayi += aay;
		azi += aaz;
	}

//printf("position %d %d %.20g %.20g %.20g %.20g %.20g\n", pp, p, time, x[p], y[p], z[p], GM_h[p]);
}


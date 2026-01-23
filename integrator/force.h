//Non Grav
//Coordinates must be heliocentric
//Marsden, Sekanina and Yeomans, 1973 Comets and nongravitational forces. V, https://adsabs.harvard.edu/full/1973AJ.....78..211M
inline void asteroid::NonGrav(double xi, double yi, double zi, double vxi, double vyi, double vzi, double A1i, double A2i, double A3i, double r, double &axi, double &ayi, double &azi){



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

	double gr = 0.0;

	if(cometFlag == 0){
		// *********Asteroids*************
		//const double alpha = 1.0;
		//const double nk = 0.0;
		//const double nm = 2.0;
		//const double nn = 5.093;
		//const double r0 = 1.0;
		//gr = alpha * pow(r / r0, -nm) * pow(1.0 + pow(r / r0, nn), -nk);
		gr = 1.0 / (r * r);
	}	
	else{
		// ************Comets**************
		//Larry uses constant values for these parameters, for all comets the same
		//add time delay for comets


		double tTau = time - nonGrav_tau;

		double rr = r / nonGrav_r0;
		double g1 = pow(rr, -nonGrav_nm);
		double g2 = pow(rr, nonGrav_nn);
		double g3 = pow(1.0 + g2, -nonGrav_nk);
		gr = nonGrav_alpha * g1 * g3;
	}


	double f1 = A1i * gr / r;
	double f2 = A2i * gr / t;
	double f3 = A3i * gr / h;


	axi += f1 * xi + f2 * tx + f3 * hx;
	ayi += f1 * yi + f2 * ty + f3 * hy;
	azi += f1 * zi + f2 * tz + f3 * hz;
	//axi += A1i * gr * xi / r + A2i * gr * tx / t + A3i * gr * hx / h;
	//ayi += A1i * gr * yi / r + A2i * gr * ty / t + A3i * gr * hy / h;
	//azi += A1i * gr * zi / r + A2i * gr * tz / t + A3i * gr * hz / h;

}

//GR
//Coordinates must be heliocentric
inline void asteroid::GR(double xi, double yi, double zi, double vxi, double vyi, double vzi, double r, double &axi, double &ayi, double &azi, double GMSun){
	//GR
	double vsq = vxi * vxi + vyi * vyi + vzi * vzi;

	double rv = xi * vxi + yi * vyi + zi * vzi;

	double f1 = GMSun / (r * r * r * c2);
	double t1 = 4.0 * GMSun / r - vsq;
	double t3 = 4.0 * rv;
	//printf("GRa %.20g %.20g %.20g %.20g %.20g\n", vsq, r, t1, t3, f1);
	//printf("a %d %.20g %.20g %.20g\n", i, ax_h, ay_h, az_h);

	//printf("A %.20g %.20g %.20g %.20g %.20g %.20g\n", xi, yi, zi, vxi, vyi, vzi);
	//printf("B %d %.20g %.20g %.20g\n", i, f1, t1, t3);

	double aax = f1 * (t1 * xi + t3 * vxi);
	double aay = f1 * (t1 * yi + t3 * vyi);
	double aaz = f1 * (t1 * zi + t3 * vzi);

	//printf("GR a %.20g %.20g %.20g\n", aax, aay, aaz);
	axi += aax;
	ayi += aay;
	azi += aaz;
	//printf("GR a %d %.20g %.20g %.20g\n", i, ax_h[i], ay_h[i], az_h[i]);
}

//Coordinates must be Earth centric
inline void asteroid::J2(double xE, double yE, double zE, double &axi, double &ayi, double &azi, double GMEarth){
	//J2

//printf("J2 %.20g %.20g %.20g\n", J2E, RE, REAU);

	double rsq = xE * xE + yE * yE + zE * zE;
	double r = sqrt(rsq);
	double r5 = rsq * rsq * r;

	double t1 = GMEarth * 3.0 * J2E * REAU * REAU / (2.0 * r5);
	double t2 = 5.0 * (zE * zE / rsq) - 1.0;

//printf("%.20g %.20g %.20g\n", GMEarth, t1, t2);

	double tx = t1 * t2 * xE;
	double ty = t1 * t2 * yE;
	double tz = t1 * (t2 - 2.0) * zE;

	axi += tx;
	ayi += ty;
	azi += tz;
//printf("J2      %.20g %.20g %.20g %.20g %.20g\n", GMEarth, xE, yE, zE, r); 

//printf("J2 a %.20g %.20g %.20g\n", tx, ty, tz);
}

inline void asteroid::Gravity(double xi, double yi, double zi, double *xTable_h, double *yTable_h, double *zTable_h, double &axi, double &ayi, double &azi, int i){

/*
//double k2 = (0.01720209895 * 0.01720209895);
double k2 = 0.017202098949957226 * 0.017202098949957226;	//Larrys value

//Masses from Larry:
GM_h[0] = 0.1660120825458948384934E-06 * k2;
GM_h[1] = 0.2447838287796943826673E-05 * k2;
GM_h[2] = 0.3003489615464968403002E-05 * k2;
GM_h[3] = 0.3227156082932277418834E-06 * k2;
GM_h[4] = 0.9547919099414246997537E-03 * k2;
GM_h[5] = 0.2858856700245945544128E-03 * k2;
GM_h[6] = 0.4366249613222119918066E-04 * k2;
GM_h[7] = 0.5151383772654273671045E-04 * k2;
GM_h[8] = 0.7350478973158631145947E-08 * k2;
GM_h[9] = 0.3694303349765110923667E-07 * k2;
GM_h[10] = 1.0 * k2;


//masses from de440 zip package AST_const.m
// masse are in m^3/s^2 need to be converted in AU^3/day^2
double AU = 149597870699.999988; //in m 
double day = 60 * 60 * 24.0;
double f = day * day / (AU * AU * AU);

GM_h[2]   = f * 398600.435507e9;
GM_h[10]  = f * 132712440041.279419e9;
GM_h[9]   = f * 4902.800118e9;
GM_h[0]   = f * 22031.868551e9;
GM_h[1]   = f * 324858.592000e9;
GM_h[3]   = f * 42828.375816e9;
GM_h[4]   = f * 126712764.100000e9;
GM_h[5]   = f * 37940584.841800e9;
GM_h[6]   = f * 5794556.400000e9;
GM_h[7]   = f * 6836527.100580e9;
GM_h[8]   = f * 975.500000e9;
*/

	for(int pp = 0; pp < Nperturbers; ++pp){


		int p = pp + 11;
		if(pp == 16) p = 8;	//Pluto
		if(pp == 17) p = 9;	//Moon
		if(pp == 18) p = 3;	//Mars
		if(pp == 19) p = 0;	//Mercury
		if(pp == 20) p = 7;	//Neptune
		if(pp == 21) p = 6;	//Uranus
		if(pp == 22) p = 2;	//Earth
		if(pp == 23) p = 1;	//Venus
		if(pp == 24) p = 5;	//Saturn
		if(pp == 25) p = 4;	//Jupiter
		if(pp == 26) p = 10;	//Sun


/*
int p = pp;
//order from Larry
if(pp ==  0) p = 0;	//Mercury
if(pp ==  1) p = 1;	//Venus
if(pp ==  2) p = 2;	//Earth
if(pp ==  3) p = 3;	//Mars
if(pp ==  4) p = 4;	//Jupiter
if(pp ==  5) p = 5;	//Saturn
if(pp ==  6) p = 6;	//Uranus
if(pp ==  7) p = 7;	//Neptune
if(pp ==  8) p = 8;	//Pluto
if(pp ==  9) p = 11 + 1;	 //Ceres
if(pp == 10) p = 11 + 11;	 //Pallas
if(pp == 11) p = 11 + 10;	 //Juno
if(pp == 12) p = 11 + 15;	 //Vesta
if(pp == 13) p = 11 + 7;	 //Hygiea
if(pp == 14) p = 11 + 4;	 //Eunomia
if(pp == 15) p = 11 + 5;	 //Euphrosyne
if(pp == 16) p = 11 + 6;	 //Europa
if(pp == 17) p = 11 + 3;	 //Davida
if(pp == 18) p = 11 + 8;	 //Interamnia
if(pp == 19) p = 11 + 12;	 //Psyche
if(pp == 20) p = 11 + 2;	 //Cybele
if(pp == 21) p = 11 + 14;	 //Thisbe
if(pp == 22) p = 11 + 0;	 //Camilla
if(pp == 23) p = 11 + 9;	 //Iris
if(pp == 24) p = 11 + 13;	 //Sylvia
if(pp == 25) p = 9;	//Moon
if(pp == 26) p = 10;	//Sun
*/

		double dx = xi - xTable_h[p];
		double dy = yi - yTable_h[p];
		double dz = zi - zTable_h[p];

		double rsq = dx * dx + dy * dy + dz * dz;
		double r = sqrt(rsq);
		double s = GM_h[p] / (r * rsq);
		//double s = GM_h[p] / (r * r * r);

		axi -= s * dx;
		ayi -= s * dy;
		azi -= s * dz;


//printf("position %d %d %25.20g %25.20g %25.20g | %25.20g\n", i, p, xTable_h[p], yTable_h[p], zTable_h[p], xi); 

	}
}

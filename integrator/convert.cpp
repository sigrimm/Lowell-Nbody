#include "asteroid.h"


void EccentricAnomaly(double M, double e, double &E){

	if(e < 1.0 - 1.0e-10){	
		//Eccentric Anomaly
		E = M + e * 0.5;
		double Eold = E;
		for(int j = 0; j < 32; ++j){
			E = E - (E - e * sin(E) - M) / (1.0 - e * cos(E));
			if(fabs(E - Eold) < 1.0e-15) break;
			Eold = E;
		}
	}
	else if(e > 1.0 + 1.0e-10){
		//hyperbolic
		//E is assumed to be the hyperbolic eccentricity 
		E = M;
		double Eold = E;
		for(int j = 0; j < 32; ++j){
			E = E + (E - e * sinh(E) + M) / (e * cosh(E) - 1.0);
			if(fabs(E - Eold) < 1.0e-15) break;
			Eold = E;
		}

	}
	else{
		//parabolic, solve Barkers equation 
		// M = D + D^3 / 3, 
		// use cot(s) = 1.5 * M  -> s = pi / 2 - atan(1.5 * M)

		//double s = M_PI * 0.5 - atan(1.5 * M);
		E = M;
		double Eold = E;
		for(int j = 0; j < 32; ++j){
			E = E - (E + E * E * E / 3.0 - M) / (1.0 + E * E);
			if(fabs(E - Eold) < 1.0e-15) break;
			Eold = E;
		}

	}
}

// *************************************
//This function converts Keplerian Elements into Cartesian Coordinates

//input a e inc Omega w M
//i is the index of the body
void asteroid::KepToCart_M(int i, double a, double e, double inc, double Omega, double w, double M){


	double mu = GM_h[10]; //Msun

	double E;
	if(e < 1.0 - 1.0e-10){	
		//Eccentric Anomaly
		E = M + e * 0.5;
		double Eold = E;
		for(int j = 0; j < 32; ++j){
			E = E - (E - e * sin(E) - M) / (1.0 - e * cos(E));
			if(fabs(E - Eold) < 1.0e-15) break;
			Eold = E;
		}
	}
	else if(e > 1.0 + 1.0e-10){
		//hyperbolic
		//E is assumed to be the hyperbolic eccentricity 
		E = M;
		double Eold = E;
		for(int j = 0; j < 32; ++j){
			E = E + (E - e * sinh(E) + M) / (e * cosh(E) - 1.0);
			if(fabs(E - Eold) < 1.0e-15) break;
			Eold = E;
		}

	}
	else{
		//parabolic, solve Barkers equation 
		// M = D + D^3 / 3, 
		// use cot(s) = 1.5 * M  -> s = pi / 2 - atan(1.5 * M)

		//double s = M_PI * 0.5 - atan(1.5 * M);
		E = M;
		double Eold = E;
		for(int j = 0; j < 32; ++j){
			E = E - (E + E * E * E / 3.0 - M) / (1.0 + E * E);
			if(fabs(E - Eold) < 1.0e-15) break;
			Eold = E;
		}

	}


	double cw = cos(w);
	double sw = sin(w);
	double cOmega = cos(Omega);
	double sOmega = sin(Omega);
	double ci = cos(inc);
	double si = sin(inc);

	double Px = cw * cOmega - sw * ci * sOmega;
	double Py = cw * sOmega + sw * ci * cOmega;
	double Pz = sw * si;

	double Qx = -sw * cOmega - cw * ci * sOmega;
	double Qy = -sw * sOmega + cw * ci * cOmega;
	double Qz = cw * si;

	double cE = cos(E);
	double sE = sin(E);

	double t0, t1, t2;

	if(e < 1.0 - 1.0e-10){
		//elliptic

		//double r = a * ( 1.0 - e * cE);
		//double r = a * (1.0 - e*e)/(1.0 + e *cos(Theta));
		//double t1 = r * cos(Theta); 
		//double t2 = r * sin(Theta); 
		t1 = a * (cE - e);
		t2 = a * sqrt(1.0 - e * e) * sE;
	}
	else if(e > 1.0 + 1.0e-10){
		//hyperbolic
		//double r = a * (1.0 - e*e)/(1.0 + e *cos(Theta));
		//or
		//double r = a * ( 1.0 - e * cosh(E));
		//t1 = r * cos(Theta); 
		//t2 = r * sin(Theta); 
		t1 = a * (cosh(E) - e);
		t2 = -a * sqrt(e * e - 1.0) * sinh(E);
	}
	else{
		//parabolic
		// a is assumed to be q, p = 2q, p = h^2/mu
		double Theta = 2.0 * atan(E);
		double r = 2 * a /(1.0 + cos(Theta));
		t1 = r * cos(Theta);
		t2 = r * sin(Theta);
	}

	x_h[i] =  t1 * Px + t2 * Qx;
	y_h[i] =  t1 * Py + t2 * Qy;
	z_h[i] =  t1 * Pz + t2 * Qz;

	if(e < 1.0 - 1.0e-10){
		//elliptic
		t0 = 1.0 / (1.0 - e * cE) * sqrt(mu / a);
		t1 = -sE;
		t2 = sqrt(1.0 - e * e) * cE;
	}
	else if(e > 1.0 + 1.0e-10){
		//hyperbolic
		//double r = a * (1.0 - e*e)/(1.0 + e *cos(Theta));
		double r = a * ( 1.0 - e * cosh(E));
		t0 = sqrt(-mu * a) / r;
		t1 = -sinh(E);
		t2 = sqrt(e * e - 1.0) * cosh(E);
	}
	else{
		//parabolic
		double Theta = 2.0 * atan(E);
		t0 = mu / sqrt(2.0 * a * mu);
		t1 = -sin(Theta);
		t2 = 1.0 +  cos(Theta);
	}

	vx_h[i] = t0 * (t1 * Px + t2 * Qx);
	vy_h[i] = t0 * (t1 * Py + t2 * Qy);
	vz_h[i] = t0 * (t1 * Pz + t2 * Qz);
}

// *************************************
//This function converts Keplerian Elements into Cartesian Coordinates

//input a e inc Omega w E
//i is the index of the body
void asteroid::KepToCart_E(int i, double a, double e, double inc, double Omega, double w, double E){

	double mu = GM_h[10]; //Msun

	double cw = cos(w);
	double sw = sin(w);
	double cOmega = cos(Omega);
	double sOmega = sin(Omega);
	double ci = cos(inc);
	double si = sin(inc);

	double Px = cw * cOmega - sw * ci * sOmega;
	double Py = cw * sOmega + sw * ci * cOmega;
	double Pz = sw * si;

	double Qx = -sw * cOmega - cw * ci * sOmega;
	double Qy = -sw * sOmega + cw * ci * cOmega;
	double Qz = cw * si;

	double cE = cos(E);
	double sE = sin(E);

	double t0, t1, t2;

	if(e < 1.0 - 1.0e-10){
		//elliptic

		//double r = a * ( 1.0 - e * cE);
		//double r = a * (1.0 - e*e)/(1.0 + e *cos(Theta));
		//double t1 = r * cos(Theta); 
		//double t2 = r * sin(Theta); 
		t1 = a * (cE - e);
		t2 = a * sqrt(1.0 - e * e) * sE;
	}
	else if(e > 1.0 + 1.0e-10){
		//hyperbolic
		//double r = a * (1.0 - e*e)/(1.0 + e *cos(Theta));
		//or
		//double r = a * ( 1.0 - e * cosh(E));
		//t1 = r * cos(Theta); 
		//t2 = r * sin(Theta); 
		t1 = a * (cosh(E) - e);
		t2 = -a * sqrt(e * e - 1.0) * sinh(E);
	}
	else{
		//parabolic
		// a is assumed to be q, p = 2q, p = h^2/mu
		double Theta = 2.0 * atan(E);
		double r = 2 * a /(1.0 + cos(Theta));
		t1 = r * cos(Theta);
		t2 = r * sin(Theta);
	}

	x_h[i] = t1 * Px + t2 * Qx;
	y_h[i] = t1 * Py + t2 * Qy;
	z_h[i] = t1 * Pz + t2 * Qz;

	if(e < 1.0 - 1.0e-10){
		//elliptic
		t0 = 1.0 / (1.0 - e * cE) * sqrt(mu / a);
		t1 = -sE;
		t2 = sqrt(1.0 - e * e) * cE;
	}
	else if(e > 1.0 + 1.0e-10){
		//hyperbolic
		//double r = a * (1.0 - e*e)/(1.0 + e *cos(Theta));
		double r = a * ( 1.0 - e * cosh(E));
		t0 = sqrt(-mu * a) / r;
		t1 = -sinh(E);
		t2 = sqrt(e * e - 1.0) * cosh(E);
	}
	else{
		//parabolic
		double Theta = 2.0 * atan(E);
		t0 = mu / sqrt(2.0 * a * mu);
		t1 = -sin(Theta);
		t2 = 1.0 +  cos(Theta);
	}

	vx_h[i] = t0 * (t1 * Px + t2 * Qx);
	vy_h[i] = t0 * (t1 * Py + t2 * Qy);
	vz_h[i] = t0 * (t1 * Pz + t2 * Qz);
//printf("B KtoC m:%g r:%g x:%g y:%g z:%g vx:%g vy:%g vz:%g\n", x.w, v.w, x.x, x.y, x.z, v.x ,v.y, v.z);
}

void asteroid::CartToKep(int i, double &a, double &e, double &inc, double &Omega, double &w, double &Theta, double &M, double &E){

	double mu = GM_h[10]; //Msun

	double x = x_h[i];
	double y = y_h[i];
	double z = z_h[i];
	double vx = vx_h[i];
	double vy = vy_h[i];
	double vz = vz_h[i];


	double rsq = x * x + y * y + z * z;
	double vsq = vx * vx + vy * vy + vz * vz;
	double u =  x * vx + y * vy + z * vz;
	double ir = 1.0 / sqrt(rsq);
	double ia = 2.0 * ir - vsq / mu;

	a = 1.0 / ia;

	//inclination
	double hx, hy, hz;

	hx = ( y * vz) - (z * vy);
	hy = (-x * vz) + (z * vx);
	hz = ( x * vy) - (y * vx);

	double h2 = hx * hx + hy * hy + hz * hz;
	double h = sqrt(h2);

	double t = hz / h;
	if(t < -1.0) t = -1.0;
	if(t > 1.0) t = 1.0;

	inc = acos(t);

	//longitude of ascending node
	double n = sqrt(hx * hx + hy * hy);
	Omega = acos(-hy / n);
	if(hx < 0.0){
		Omega = 2.0 * M_PI - Omega;
	}

	if(inc < 1.0e-10 || n == 0) Omega = 0.0;

	//argument of periapsis
	double ex, ey, ez;

	ex = ( vy * hz - vz * hy) / mu - x * ir;
	ey = (-vx * hz + vz * hx) / mu - y * ir;
	ez = ( vx * hy - vy * hx) / mu - z * ir;


	e = sqrt(ex * ex + ey * ey + ez * ez);

	t = (-hy * ex + hx * ey) / (n * e);
	if(t < -1.0) t = -1.0;
	if(t > 1.0) t = 1.0;
	w = acos(t);
	if(ez < 0.0) w = 2.0 * M_PI - w;
	if(n == 0) w = 0.0;

	//True Anomaly
	t = (ex * x + ey * y + ez * z) / e * ir;
	if(t < -1.0) t = -1.0;
	if(t > 1.0) t = 1.0;
	Theta = acos(t);

	if(u < 0.0){
		if(e < 1.0 - 1.0e-10){
			//elliptic
			Theta = 2.0 * M_PI - Theta;
		}
		else if(e > 1.0 + 1.0e-10){
			//hyperbolic
			Theta = -Theta;
		}
		else{
			//parabolic
			Theta = - Theta;
		}
	}

	//Non circular, equatorial orbit
	if(e > 1.0e-10 && inc < 1.0e-10){
		Omega = 0.0;
		w = acos(ex / e);
		if(ey < 0.0) w = 2.0 * M_PI - w;
	}

	//circular, inclinded orbit
	if(e <= 1.0e-10 && inc > 1.0e-11){
		w = 0.0;
	}

	//circular, equatorial orbit
	if(e <= 1.0e-10 && inc <= 1.0e-11){
		w = 0.0;
		Omega = 0.0;
	}

	if(w == 0 && Omega != 0.0){
		t = (-hy * x + hx * y) / n * ir;
		if(t < -1.0) t = -1.0;
		if(t > 1.0) t = 1.0;
		Theta = acos(t);
		if(z < 0.0){
			if(e < 1.0 - 1.0e-10){
				//elliptic
				Theta = 2.0 * M_PI - Theta;
			}
			else if(e > 1.0 + 1.0e-10){
				//hyperbolic
				Theta = -Theta;
			}
			else{
				//parabolic
				Theta = -Theta;
			}
		}
	}
	if(w == 0 && Omega == 0.0){
		Theta = acos(x * ir);
		if(y < 0.0){
			if(e < 1.0 - 1.0e-10){
				//elliptic
				Theta = 2.0 * M_PI - Theta;
			}
			else if(e > 1.0 + 1.0e-10){
				//hyperbolic
				Theta = -Theta;
			}
			else{
				//parabolic
				Theta = -Theta;
			}
		}
	}

	if(e < 1.0 - 1.0e-10){
		//Eccentric Anomaly
		E = acos((e + cos(Theta)) / (1.0 + e * cos(Theta)));
		if(M_PI < Theta && Theta < 2.0 * M_PI) E = 2.0 * M_PI - E;

		//Mean Anomaly
		M = E - e * sin(E);
//printf("%g %g %g %g\n", Theta, E, M, w);
	}
	else if(e > 1.0 + 1.0e-10){
		//Hyperbolic Anomaly
		//named still E instead of H or F
		E = acosh((e + t) / (1.0 + e * t));
		if(Theta < 0.0) E = - E;

		M = e * sinh(E) - E;
	}
	else{
		//Parabolic Anomaly
		E = tan(Theta * 0.5);
		if(E > M_PI) E = E - 2.0 * M_PI;

		M = E + E * E * E / 3.0;

		//use a to store q
		a = h * h / mu * 0.5;
	}
}


void asteroid::EcpliptictoEquatorial(){

	//Rotate around obliquity angle eps
	double eps = Obliquity / 3600.0 / 180.0 * M_PI;	// convert arcseconds to radians

	for(int i = 0; i < N; ++i){

		double x = x_h[i];
		double y = y_h[i];
		double z = z_h[i];
		double vx = vx_h[i];
		double vy = vy_h[i];
		double vz = vz_h[i];


		double ceps = cos(eps);
		double seps = sin(eps);

		x_h[i] = x;
		y_h[i] = ceps * y - seps * z;
		z_h[i] = seps * y + ceps * z;

		vx_h[i] = vx;
		vy_h[i] = ceps * vy - seps * vz;
		vz_h[i] = seps * vy + ceps * vz;

		printf("xyz %.40g %.40g %.40g %.40g %.40g %.40g %.40g %.40g %.40g\n", x_h[i], y_h[i], z_h[i], vx_h[i], vy_h[i], vz_h[i], A1_h[i], A2_h[i], A3_h[i]);

	}
}





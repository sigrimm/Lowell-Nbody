#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void aei(double xi, double yi, double zi, double vxi, double vyi, double vzi, double mu, double &a, double &e, double &inc, double &Omega, double &w, double &M){

	double Theta, E;
	double rsq = xi * xi + yi * yi + zi * zi;
	double vsq = vxi * vxi + vyi * vyi + vzi * vzi;
	double u =  xi * vxi + yi * vyi + zi * vzi;
	double ir = 1.0 / sqrt(rsq);
	double ia = 2.0 * ir - vsq / mu;

	a = 1.0 / ia;

	//inclination
	double hx, hy, hz;
	double h2, h, t;
	hx = ( yi * vzi) - (zi * vyi);
	hy = (-xi * vzi) + (zi * vxi);
	hz = ( xi * vyi) - (yi * vxi);

	h2 = hx * hx + hy * hy + hz * hz;
	h = sqrt(h2);

	t = hz / h;
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
	ex = ( vyi * hz - vzi * hy) / mu - xi * ir;
	ey = (-vxi * hz + vzi * hx) / mu - yi * ir;
	ez = ( vxi * hy - vyi * hx) / mu - zi * ir;


	e = sqrt(ex * ex + ey * ey + ez * ez);

	t = (-hy * ex + hx * ey) / (n * e);
	if(t < -1.0) t = -1.0;
	if(t > 1.0) t = 1.0;
	w = acos(t);
	if(ez < 0.0) w = 2.0 * M_PI - w;
	if(n == 0) w = 0.0;

	//True Anomaly
	t = (ex * xi + ey * yi + ez * zi) / e * ir;
	if(t < -1.0) t = -1.0;
	if(t > 1.0) t = 1.0;
	Theta = acos(t);
	if(u < 0.0) Theta = 2.0 * M_PI - Theta;


	//Non circular, equatorial orbit
	if(e > 1.0e-10 && inc < 1.0e-10){
		Omega = 0.0;
		w = acos(ex / e);
		if(ey < 0.0) w = 2.0 * M_PI - w;
	}

	//circular, inclinded orbit
		if(e < 1.0e-10 && inc > 1.0e-11){
		w = 0.0;
	}

	//circular, equatorial orbit
	if(e < 1.0e-10 && inc < 1.0e-11){
		w = 0.0;
		Omega = 0.0;
	}

	if(w == 0 && Omega != 0.0){
		t = (-hy * xi + hx * yi) / n * ir;
		if(t < -1.0) t = -1.0;
		if(t > 1.0) t = 1.0;
		Theta = acos(t);
		if(zi < 0.0) Theta = 2.0 * M_PI - Theta;
	}
	if(w == 0 && Omega == 0.0){
		Theta = acos(xi * ir);
		if(yi < 0.0) Theta = 2.0 * M_PI - Theta;

	}

	//Eccentric Anomaly
	E = acos((e + cos(Theta)) / (1.0 + e * cos(Theta)));
	if(M_PI < Theta && Theta < 2.0 * M_PI) E = 2.0 * M_PI - E;

	//Mean Anomaly
	M = E - e * sin(E);

	if(e >= 1){
		E = acosh((e + t) / (1.0 + e * t));
		if(M_PI < Theta && Theta < 2.0 * M_PI) E = 2.0 * M_PI - E;
		M = E - e * sinh(E);
	}

}


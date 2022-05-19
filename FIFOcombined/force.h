
// --------------------------------
//heliocentric coordinates
//sun part
__host__ __device__ void accS(double mu, double xi, double yi, double zi, double &ax, double &ay, double &az){

	double rx = -xi;
	double ry = -yi;
	double rz = -zi;
	double rsq = rx * rx + ry * ry + rz * rz;
	double r = sqrt(rsq);

	double s = mu / (r * rsq);

	ax += s * rx;
	ay += s * ry;
	az += s * rz;
}
//planet part
__host__ __device__ void accP(double mj, double xj, double yj, double zj, double xi, double yi, double zi, double &ax, double &ay, double &az){

	double rx = xj - xi;
	double ry = yj - yi;
	double rz = zj - zi;
	double rsq = rx * rx + ry * ry + rz * rz;
	double r = sqrt(rsq);

	double s = mj / (r * rsq);

	ax += s * rx;
	ay += s * ry;
	az += s * rz;
}
//planet part 2

__host__ __device__ void accP2(double mj, double xj, double yj, double zj, double &ax, double &ay, double &az){

	double rx = -xj;
	double ry = -yj;
	double rz = -zj;
	double rsq = rx * rx + ry * ry + rz * rz;
	double r = sqrt(rsq);

	double s = mj / (r * rsq);
	ax += s * rx;
	ay += s * ry;
	az += s * rz;
}
// --------------------------------------


//Sitarski 1982, Isotropic equation 5, heliocentric
//modified k2 to dayUnit
//should be equivalent to the Quinn et all function, assuming m[0] = 1.0
//heliocentric
__device__ __host__ void acchGR2(double xi, double yi, double zi, double vxi, double vyi, double vzi, double &ax, double &ay, double &az){
	
	double c2 = def_c * def_c;

	double rsq = xi * xi + yi * yi + zi * zi;
	double r = sqrt(rsq);
	double vsq = vxi * vxi + vyi * vyi + vzi * vzi;

	double rv = xi * vxi + yi * vyi + zi * vzi;

	double f1 = 1.0 / (r * rsq * c2);
	double t1 = 4.0 / r;
	double t2 = -vsq;
	double t3 = 4.0 * rv;
//printf("a %d %.20g %.20g %.20g\n", i, ax, ay, az);

//printf("A %d %.20g %.20g %.20g %.20g %.20g %.20g\n", i, xi, yi, zi, vxi, vyi, vzi);
//printf("B %d %.20g %.20g %.20g %.20g\n", i, f1, t1, t2, t3);
//printf("C %d %.20g %.20g %.20g %.20g\n", i, t1 + t2, (t1 + t2) * xi, ((t1 + t2) * xi + t3 * vxi), f1 * ((t1 + t2) * xi + t3 * vxi));

	ax += f1 * ((t1 + t2) * xi + t3 * vxi);
	ay += f1 * ((t1 + t2) * yi + t3 * vyi);
	az += f1 * ((t1 + t2) * zi + t3 * vzi);
//printf("D %d %.20g %.20g %.20g\n", i, ax, ay, az);

}



//A1, A2 and A3 terms for asteroids on heliocentric coordinates
__host__ __device__ void NonGrav(double xi, double yi, double zi, double vxi, double vyi, double vzi, double &ax, double &ay, double &az, double A1i, double A2i, double A3i){

	double rsq = xi * xi + yi * yi + zi * zi;
	double r = sqrt(rsq);

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

	double gr = 1.0 / rsq;	//only valid for asteroids, not for comets 
/*
	double rr = r / R0[i];
	double g1 = pow(rr, -NM[i]);
	double g2 = pow(rr, Nn[i]);
	double g3 = pow(1.0 + g2, -NK[i]);
	double gr = ALN[i] * g1 * g3;

printf("gr %.20g %.20g\n", gr1, gr);
*/

	double f1 = A1i * gr / r;
	double f2 = A2i * gr / t;
	double f3 = A3i * gr / h;
	
	
	ax += f1 * xi + f2 * tx + f3 * hx;
	ay += f1 * yi + f2 * ty + f3 * hy;
	az += f1 * zi + f2 * tz + f3 * hz;
//printf("NonGrav %d %.20g %.20g %.20g\n", i, (f1 * x[i] + f2 * tx + f3 * hx) * dayUnit * dayUnit, (f1 * y[i] + f2 * ty + f3 * hy) * dayUnit * dayUnit, (f1 * z[i] + f2 * tz + f3 * hz) * dayUnit * dayUnit);

}



//J2 perturbation from Earth
//Walter 2018, 12.2.10
__host__ __device__ void J2(double mj, double xj, double yj, double zj, double xi, double yi, double zi, double &ax, double &ay, double &az){

	double J2E = 0.00108262545; // J2 Earth from DE 430
	double RE = 6378136.3; // Earth radius in m from DE 430

	//double J2E = 1.08263e-3; //1.08262668e-3;
	//double RE = 6371.009; // Earth radius in km
	//double muE = 398600.44 // in km^3 /s^2	G * mEarth

	RE /= def_AU;	//Earth radius in AU

	//int iE = 3; 	//index of Earth

	//muE = 3.986004415e5 km^3s-2
	double muE = mj;

	double xE = xi - xj;
	double yE = yi - yj;
	double zE = zi - zj;

	double rsq = xE * xE + yE * yE + zE * zE;
	double r = sqrt(rsq);
	double r5 = rsq * rsq * r;

	double t1 = 3.0 * J2E * muE * RE * RE / (2.0 * r5);
	double t2 = 5.0 * zE * zE / rsq;

//printf("rm %.20g %.20g %.20g\n", RE, muE, t1);

	double tx = t1 * (t2 - 1.0) * xE;
	double ty = t1 * (t2 - 1.0) * yE;
	double tz = t1 * (t2 - 3.0) * zE;
	
	ax += tx;
	ay += ty;
	az += tz;

//printf("J2 %d %.20g %.20g %.20g %.20g | %.20g %.20g %.20g\n", i, r, tx * dayUnit * dayUnit, ty * dayUnit * dayUnit, tz * dayUnit * dayUnit, xE, yE, zE); 


}


#include <stdio.h>

#include "asteroid.h"
#include "Chebyshev.h"
#include "force.h"

int asteroid::HelioToBary(double *xx_h, double *yy_h, double *zz_h, double *vxx_h, double *vyy_h, double *vzz_h){
	//Update the Chebyshev coefficients if necessary
	int er;
	er = update_Chebyshev(timeStart);
	if(er <= 0){
		return 0;
	}
	update_perturbers(timeStart);


	//for(int i = 0; i < Nperturbers; ++i){
	//	printf("%d %.20g %.20g %.20g %.20g %.20g %.20g\n", i, xTable_h[i], yTable_h[i], zTable_h[i], vxTable_h[i], vyTable_h[i], vzTable_h[i]);
	//}

	//heliocentric coordinates
	for(int i = 0; i < N; ++i){
		xx_h[i] += xTable_h[10];
		yy_h[i] += yTable_h[10];
		zz_h[i] += zTable_h[10];

		vxx_h[i] += vxTable_h[10];
		vyy_h[i] += vyTable_h[10];
		vzz_h[i] += vzTable_h[10];


		//printf("H2B %.40g %.40g %.40g %.40g %.40g %.40g %.40g %.40g %.40g\n", xx_h[i], yy_h[i], zz_h[i], vxx_h[i], vyy_h[i], vzz_h[i], A1_h[i], A2_h[i], A3_h[i]);
	}
	return 1;
}

int asteroid::BaryToHelio(double *xx_h, double *yy_h, double *zz_h, double *vxx_h, double *vyy_h, double *vzz_h){
	//Update the Chebyshev coefficients if necessary
	int er;
	er = update_Chebyshev(time);
	if(er <= 0){
		return 0;
	}
	update_perturbers(time);


	//for(int i = 0; i < Nperturbers; ++i){
	//	printf("%d %.20g %.20g %.20g %.20g %.20g %.20g\n", i, xTable_h[i], yTable_h[i], zTable_h[i], vxTable_h[i], vyTable_h[i], vzTable_h[i]);
	//}

	//heliocentric coordinates
	for(int i = 0; i < N; ++i){
		xx_h[i] -= xTable_h[10];
		yy_h[i] -= yTable_h[10];
		zz_h[i] -= zTable_h[10];

		vxx_h[i] -= vxTable_h[10];
		vyy_h[i] -= vyTable_h[10];
		vzz_h[i] -= vzTable_h[10];
	}
	return 1;
}

int asteroid::BaryToGeo(double *xx_h, double *yy_h, double *zz_h, double *vxx_h, double *vyy_h, double *vzz_h){
	//Update the Chebyshev coefficients if necessary
	int er;
	er = update_Chebyshev(time);
	if(er <= 0){
		return 0;
	}
	update_perturbers(time);


	//for(int i = 0; i < Nperturbers; ++i){
	//	printf("%d %.20g %.20g %.20g %.20g %.20g %.20g\n", i, xTable_h[i], yTable_h[i], zTable_h[i], vxTable_h[i], vyTable_h[i], vzTable_h[i]);
	//}

	//heliocentric coordinates
	for(int i = 0; i < N; ++i){
		xx_h[i] -= xTable_h[2];
		yy_h[i] -= yTable_h[2];
		zz_h[i] -= zTable_h[2];

		vxx_h[i] -= vxTable_h[2];
		vyy_h[i] -= vyTable_h[2];
		vzz_h[i] -= vzTable_h[2];
//printf("Earth %.20g %.20g %.20g\n", xTable_h[2], yTable_h[2], zTable_h[2]);
	}
	return 1;
}

int asteroid::convertOutput(){
	
	int er;
	for(int i = 0; i < N; ++i){
		xout_h[i] = x_h[i];
		yout_h[i] = y_h[i];
		zout_h[i] = z_h[i];

		vxout_h[i] = vx_h[i];
		vyout_h[i] = vy_h[i];
		vzout_h[i] = vz_h[i];
	}

	if(Outheliocentric == 1){
		er = BaryToHelio(xout_h, yout_h, zout_h, vxout_h, vyout_h, vzout_h);
		if(er <= 0){
			return 0;
		}
	} 
	if(Outgeocentric == 1){
		er = BaryToGeo(xout_h, yout_h, zout_h, vxout_h, vyout_h, vzout_h);
		if(er <= 0){
			return 0;
		}
	} 


	//If needed, convert from equatorial coordinates to ecliptic coordinates
	if(Outecliptic == 1){
		EquatorialtoEcliptic(xout_h, yout_h, zout_h, vxout_h, vyout_h, vzout_h);
	}

	if(Outorbital == 1){
		for(int i = 0; i < N; ++i){
			double a, e, inc, Omega, w, Theta, M, E;
			CartToKep(xout_h, yout_h, zout_h, vxout_h, vyout_h, vzout_h, i, a, e, inc, Omega, w, Theta, M, E);


			inc = inc * 180.0 / M_PI;	//convert rad to deg
			Omega = Omega * 180.0 / M_PI;	//convert rad to deg
			w = w * 180.0 / M_PI;		//convert rad to deg
			M = M * 180.0 / M_PI;		//convert rad to deg

			xout_h[i] = a;
			yout_h[i] = e;
			zout_h[i] = inc;
			vxout_h[i] = Omega;
			vyout_h[i] = w;
			vzout_h[i] = M;

			//printf("aei %.40g %.40g %.40g %.40g %.40g %.40g %.40g %.40g %.40g\n", a, e, inc, Omega, w, M, A1_h[i], A2_h[i], A3_h[i]);

		}
	}
	return 1;
}



//Leapfrog step with fixed time step
inline int asteroid::leapfrog_step(){
	int er;
	//Drift
	for(int i = 0; i < N; ++i){
		x_h[i] += 0.5 * dt * vx_h[i];
		y_h[i] += 0.5 * dt * vy_h[i];
		z_h[i] += 0.5 * dt * vz_h[i];
	}

	time += dt * 0.5;
//printf("ta %.20g %.20g\n", time, dt);


	// ----------------------------------------------------------------------------
	//Update the Chebyshev coefficients if necessary

	er = update_Chebyshev(time);
	if(er <= 0){
		return 0;
	}

	update_perturbers(time);
	// ----------------------------------------------------------------------------


	// ----------------------------------------------------------------------------
	for(int i = 0; i < N; ++i){

		//compute force
		double ax = 0.0;
		double ay = 0.0;
		double az = 0.0;

		//heliocentric coordinates
		double xih = x_h[i] - xTable_h[10];
		double yih = y_h[i] - yTable_h[10];
		double zih = z_h[i] - zTable_h[10];

		double vxih = vx_h[i] - vxTable_h[10];
		double vyih = vy_h[i] - vyTable_h[10];
		double vzih = vz_h[i] - vzTable_h[10];

		//r is used in multiple forces, so reuse it
		double rsq = xih * xih + yih * yih + zih * zih;
		double r = sqrt(rsq);
		if(cometFlag > 0){
			if(i == 0){
				Tsave_h[timeStep] = time;
			}
			Rsave_h[i * Rbuffersize + timeStep] = r;
		}
		//Earth centric coordinates
		double xiE = x_h[i] - xTable_h[2];
		double yiE = y_h[i] - yTable_h[2];
		double ziE = z_h[i] - zTable_h[2];

		if(useNonGrav == 1){
			NonGrav(xih, yih, zih, vxih, vyih, vzih, A1_h[i], A2_h[i], A3_h[i], r, ax, ay, az);
		}
		if(useGR == 1){
			GR(xih, yih, zih, vxih, vyih, vzih, r, ax, ay, az, GM_h[10]);
		}
		if(useJ2 == 1){
			J2(xiE, yiE, ziE, ax, ay, az, GM_h[2]);
		}
		Gravity(x_h[i], y_h[i], z_h[i], xTable_h, yTable_h, zTable_h, ax, ay, az, i);


		//Kick
		vx_h[i] += dt * ax;
		vy_h[i] += dt * ay;
		vz_h[i] += dt * az;


		//Drift
		x_h[i] += 0.5 * dt * vx_h[i];
		y_h[i] += 0.5 * dt * vy_h[i];
		z_h[i] += 0.5 * dt * vz_h[i];
	}
	// ----------------------------------------------------------------------------


	time += dt * 0.5;
	++timeStep;
//printf("tb %.20g %.20g\n", time, dt); 
	return 1;
}


//Implicit Midpoint Method
// All bodies are integrated individually with an individual, adaptive time step
inline int asteroid::IMM_step(){

	int er;
	// ----------------------------------------------------------------------------
	//Update the Chebyshev coefficients if necessary
	er = update_Chebyshev(time + 0.5 * dt);
	if(er <= 0){
		return 0;
	}
	update_perturbers(time + 0.5 * dt);
	// ----------------------------------------------------------------------------


	for(int i = 0; i < N; ++i){

		double ax;
		double ay;
		double az;

		double xt = x_h[i];
		double yt = y_h[i];
		double zt = z_h[i];

		double vxt = vx_h[i];
		double vyt = vy_h[i];
		double vzt = vz_h[i];

		double xp;
		double yp;
		double zp;

		double vxp;
		double vyp;
		double vzp;



		for(int k = 0; k < 30; ++k){

			// ----------------------------------------------------------------------------
			//compute forces

			ax = 0.0;
			ay = 0.0;
			az = 0.0;

			//heliocentric coordinates
			double xih = xt - xTable_h[10];
			double yih = yt - yTable_h[10];
			double zih = zt - zTable_h[10];

			double vxih = vxt - vxTable_h[10];
			double vyih = vyt - vyTable_h[10];
			double vzih = vzt - vzTable_h[10];

			//r is used in multiple forces, so reuse it
			double rsq = xih * xih + yih * yih + zih * zih;
			double r = sqrt(rsq);

			if(cometFlag > 0 && k == 0){
				if(i == 0){
					Tsave_h[timeStep] = time;
				}
				Rsave_h[i * Rbuffersize + timeStep] = r;
			}

			//Earth centric coordinates
			double xiE = xt - xTable_h[2];
			double yiE = yt - yTable_h[2];
			double ziE = zt - zTable_h[2];

			if(useNonGrav == 1){
				NonGrav(xih, yih, zih, vxih, vyih, vzih, A1_h[i], A2_h[i], A3_h[i], r, ax, ay, az);
			}
			if(useGR == 1){
				GR(xih, yih, zih, vxih, vyih, vzih, r, ax, ay, az, GM_h[10]);
			}
			if(useJ2 == 1){
				J2(xiE, yiE, ziE, ax, ay, az, GM_h[2]);
			}
			Gravity(xt, yt, zt, xTable_h, yTable_h, zTable_h, ax, ay, az, i);

			// ----------------------------------------------------------------------------

			int stop = 1;
			//store old values to check for convergence
			xp = xt;
			yp = yt;
			zp = zt;

			vxp = vxt;
			vyp = vyt;
			vzp = vzt;

			xt = x_h[i] + 0.5 * dt * vxt;
			yt = y_h[i] + 0.5 * dt * vyt;
			zt = z_h[i] + 0.5 * dt * vzt;

			vxt = vx_h[i] + 0.5 * dt * ax;
			vyt = vy_h[i] + 0.5 * dt * ay;
			vzt = vz_h[i] + 0.5 * dt * az;

			if(fabs(xp - xt) >= 1.0e-20) stop = 0;
			if(fabs(yp - yt) >= 1.0e-20) stop = 0;
			if(fabs(zp - zt) >= 1.0e-20) stop = 0;

			if(fabs(vxp - vxt) >= 1.0e-20) stop = 0;
			if(fabs(vyp - vyt) >= 1.0e-20) stop = 0;
			if(fabs(vzp - vzt) >= 1.0e-20) stop = 0;

			if(stop == 1){
//if(k > 1) printf("k %d %d\n", i, k);
				break;
			}
			if(k >= 29){
				snew_h[0] = -1;
			}



		}// end of k loop

		
		//update
		x_h[i] += dt * vxt;
		y_h[i] += dt * vyt;
		z_h[i] += dt * vzt;

		vx_h[i] += dt * ax;
		vy_h[i] += dt * ay;
		vz_h[i] += dt * az;
//printf("%d %g %g %g\n", i, x_h[i], y_h[i], z_h[i]);




	} //end of i loop

	time += dt;
	++timeStep;
	return 1;
}


//Runge Kutta step with fixed time step
inline int asteroid::RK_step(){

	int er;
	for(int S = 0; S < RKFn; ++S){

		// ----------------------------------------------------------------------------
		//Update the Chebyshev coefficients if necessary
		er = update_Chebyshev(time + RKFc_h[S] * dt);
		if(er <= 0){
			return 0;
		}
		update_perturbers(time + RKFc_h[S] * dt);
		// ----------------------------------------------------------------------------

		for(int i = 0; i < N; ++i){
			double xt = x_h[i];
			double yt = y_h[i];
			double zt = z_h[i];

			double vxt = vx_h[i];
			double vyt = vy_h[i];
			double vzt = vz_h[i];

			for(int s = 0; s < S; ++s){
				double dtaa = dt * RKFa_h[S * RKFn + s];
				xt  += dtaa * kx_h[i + s * N];
				yt  += dtaa * ky_h[i + s * N];
				zt  += dtaa * kz_h[i + s * N];
				vxt += dtaa * kvx_h[i + s * N];
				vyt += dtaa * kvy_h[i + s * N];
				vzt += dtaa * kvz_h[i + s * N];
//printf("update 2 %d %d %g %g %g %g %g %g\n", S, i, xt, yt, zt, RKFa_h[S * RKFn + s], kx_h[s], dt);

			}

			kx_h[i + S * N] = vxt;
			ky_h[i + S * N] = vyt;
			kz_h[i + S * N] = vzt;


			// ----------------------------------------------------------------------------
			//compute forces

			double ax = 0.0;
			double ay = 0.0;
			double az = 0.0;


			//heliocentric coordinates
			double xih = xt - xTable_h[10];
			double yih = yt - yTable_h[10];
			double zih = zt - zTable_h[10];

			double vxih = vxt - vxTable_h[10];
			double vyih = vyt - vyTable_h[10];
			double vzih = vzt - vzTable_h[10];

			//r is used in multiple forces, so reuse it
			double rsq = xih * xih + yih * yih + zih * zih;
			double r = sqrt(rsq);

			if(cometFlag > 0 && S == 0){
				if(i == 0){
					Tsave_h[timeStep] = time;
				}
				Rsave_h[i * Rbuffersize + timeStep] = r;
			}

			//Earth centric coordinates
			double xiE = xt - xTable_h[2];
			double yiE = yt - yTable_h[2];
			double ziE = zt - zTable_h[2];

			if(useNonGrav == 1){
				NonGrav(xih, yih, zih, vxih, vyih, vzih, A1_h[i], A2_h[i], A3_h[i], r, ax, ay, az);
			}
			if(useGR == 1){
				GR(xih, yih, zih, vxih, vyih, vzih, r, ax, ay, az, GM_h[10]);
			}
			if(useJ2 == 1){
				J2(xiE, yiE, ziE, ax, ay, az, GM_h[2]);
			}
			Gravity(xt, yt, zt, xTable_h, yTable_h, zTable_h, ax, ay, az, i);

			kvx_h[i + S * N] = ax;
			kvy_h[i + S * N] = ay;
			kvz_h[i + S * N] = az;
			// ----------------------------------------------------------------------------
		}
	}

	//update
	for(int i = 0; i < N; ++i){

		double dx = 0.0;
		double dy = 0.0;
		double dz = 0.0;

		double dvx = 0.0;
		double dvy = 0.0;
		double dvz = 0.0;

		for(int S = 0; S < RKFn; ++S){
			double dtb = dt * RKFb_h[S];
			dx += dtb * kx_h[i + S * N];
			dy += dtb * ky_h[i + S * N];
			dz += dtb * kz_h[i + S * N];

			dvx += dtb * kvx_h[i + S * N];
			dvy += dtb * kvy_h[i + S * N];
			dvz += dtb * kvz_h[i + S * N];
		}

		x_h[i] += dx;
		y_h[i] += dy;
		z_h[i] += dz;

		vx_h[i] += dvx;
		vy_h[i] += dvy;
		vz_h[i] += dvz;
//printf("%d %g %g %g\n", i, x_h[i], y_h[i], z_h[i]);

	}
	time += dt;
	++timeStep;
	return 1;
}

// Runge Kutta Fehlberg step with adaptive time steps
// All bodies are integrated individually with an individual, adaptive time steps
inline int asteroid::RKF_individual_step(const int i){

	int er;

	double dt = dt_h[i];
	double ax;
	double ay;
	double az;

	double xt;
	double yt;
	double zt;

	double vxt;
	double vyt;
	double vzt;

	for(int S = 0; S < RKFn; ++S){

		// ----------------------------------------------------------------------------
		//Update the Chebyshev coefficients if necessary
		er = update_Chebyshev(time_h[i] + RKFc_h[S] * dt);
		if(er <= 0){
			return 0;
		}
 		update_perturbers(time_h[i] + RKFc_h[S] * dt);
		// ----------------------------------------------------------------------------

		ax = 0.0;
		ay = 0.0;
		az = 0.0;

		xt = x_h[i];
		yt = y_h[i];
		zt = z_h[i];

		vxt = vx_h[i];
		vyt = vy_h[i];
		vzt = vz_h[i];

		for(int s = 0; s < S; ++s){
			double dtaa = dt * RKFa_h[S * RKFn + s];
			int ii = i + s * N;
			xt  += dtaa * kx_h[ii];
			yt  += dtaa * ky_h[ii];
			zt  += dtaa * kz_h[ii];
			vxt += dtaa * kvx_h[ii];
			vyt += dtaa * kvy_h[ii];
			vzt += dtaa * kvz_h[ii];

		}
		int ii = i + S * N;
		kx_h[ii] = vxt;
		ky_h[ii] = vyt;
		kz_h[ii] = vzt;
//if(i < 2) printf("K %d %d %.20g %.20g %.20g\n", i, S, kx_h[i + S * N], ky_h[i + S * N], kz_h[i + S * N]);

		// ----------------------------------------------------------------------------
		//compute forces

		//heliocentric coordinates
		double xih = xt - xTable_h[10];
		double yih = yt - yTable_h[10];
		double zih = zt - zTable_h[10];

		double vxih = vxt - vxTable_h[10];
		double vyih = vyt - vyTable_h[10];
		double vzih = vzt - vzTable_h[10];

		//r is used in multiple forces, so reuse it
		double rsq = xih * xih + yih * yih + zih * zih;
		double r = sqrt(rsq);

//		if(cometFlag > 0 && S == 0){
//			if(i == 0){
//				Tsave_h[timeStep_h[i]] = time_h[i];
//			}
//			Rsave_h[i * Rbuffersize + timeStep_h[i]] = r;
//		}

		//Earth centric coordinates
		double xiE = xt - xTable_h[2];
		double yiE = yt - yTable_h[2];
		double ziE = zt - zTable_h[2];

		if(useNonGrav == 1){
			NonGrav(xih, yih, zih, vxih, vyih, vzih, A1_h[i], A2_h[i], A3_h[i], r, ax, ay, az);
		}
		if(useGR == 1){
			GR(xih, yih, zih, vxih, vyih, vzih, r, ax, ay, az, GM_h[10]);
		}
		if(useJ2 == 1){
			J2(xiE, yiE, ziE, ax, ay, az, GM_h[2]);
		}
		Gravity(xt, yt, zt, xTable_h, yTable_h, zTable_h, ax, ay, az, i);
		// ----------------------------------------------------------------------------

		kvx_h[ii] = ax;
		kvy_h[ii] = ay;
		kvz_h[ii] = az;
	}

	//update

	dx_h[i] = 0.0;
	dy_h[i] = 0.0;
	dz_h[i] = 0.0;

	dvx_h[i] = 0.0;
	dvy_h[i] = 0.0;
	dvz_h[i] = 0.0;

	for(int S = 0; S < RKFn; ++S){
		double dtb = dt * RKFb_h[S];
		int ii = i + S * N;
		dx_h[i] += dtb * kx_h[ii];
		dy_h[i] += dtb * ky_h[ii];
		dz_h[i] += dtb * kz_h[ii];

		dvx_h[i] += dtb * kvx_h[ii];
		dvy_h[i] += dtb * kvy_h[ii];
		dvz_h[i] += dtb * kvz_h[ii];
//if(i < 2) printf("dx %d %d %g %g %g %g %g %g\n", S, i, dx_h[i], dy_h[i], dz_h[i], RKFb_h[S], kx_h[i + S * N], dt);
	}



	//compute integration error

	double scalex;
	double scaley;
	double scalez;

	double scalevx;
	double scalevy;
	double scalevz;



	//Hairer
/*
	scalex	= atol + fmax(fabs(x_h[i]), fabs(x_h[i] + dx_h[i])) * rtol;
	scaley	= atol + fmax(fabs(y_h[i]), fabs(y_h[i] + dy_h[i])) * rtol;
	scalez	= atol + fmax(fabs(z_h[i]), fabs(z_h[i] + dz_h[i])) * rtol;

	scalevx = atol + fmax(fabs(vx_h[i]), fabs(vx_h[i] + dvx_h[i])) * rtol;
	scalevy = atol + fmax(fabs(vy_h[i]), fabs(vy_h[i] + dvy_h[i])) * rtol;
	scalevz = atol + fmax(fabs(vz_h[i]), fabs(vz_h[i] + dvz_h[i])) * rtol;
*/

	scalex	= atol + fabs(x_h[i]) * rtol;
	scaley	= atol + fabs(y_h[i]) * rtol;
	scalez	= atol + fabs(z_h[i]) * rtol;

	scalevx = atol + fabs(vx_h[i]) * rtol;
	scalevy = atol + fabs(vy_h[i]) * rtol;
	scalevz = atol + fabs(vz_h[i]) * rtol;

	//error estimation
	double errorkx = 0.0;
	double errorky = 0.0;
	double errorkz = 0.0;
	
	double errorkvx = 0.0;
	double errorkvy = 0.0;
	double errorkvz = 0.0;

	for(int S = 0; S < RKFn; ++S){
		double f = (RKFb_h[S] - RKFbb_h[S]) * dt;
		int ii = i + S * N;
		errorkx += f * kx_h[ii];
		errorky += f * ky_h[ii];
		errorkz += f * kz_h[ii];

		errorkvx += f * kvx_h[ii];
		errorkvy += f * kvy_h[ii];
		errorkvz += f * kvz_h[ii];
	}

	double errork = 0.0;
	errork += errorkx * errorkx / (scalex * scalex);
	errork += errorky * errorky / (scaley * scaley);
	errork += errorkz * errorkz / (scalez * scalez);
	errork += errorkvx * errorkvx / (scalevx * scalevx);
	errork += errorkvy * errorkvy / (scalevy * scalevy);
	errork += errorkvz * errorkvz / (scalevz * scalevz);

	errork = sqrt(errork / 6.0);	//6 is the number of dimensions

	double s = pow( 1.0  / errork, RKF_ee);

	if(s >= 1.0){
		//accept step
		x_h[i] += dx_h[i];
		y_h[i] += dy_h[i];
		z_h[i] += dz_h[i];

		vx_h[i] += dvx_h[i];
		vy_h[i] += dvy_h[i];
		vz_h[i] += dvz_h[i];

		time_h[i] += dt;
		++timeStep_h[i];

	}

	s = (RKF_fac * s > RKF_facmin) ? RKF_fac * s : RKF_facmin;
	snew_h[i] = (RKF_facmax < s) ? RKF_facmax : s;


//printf("dt %.20g %.20g %.20g\n", time_h[i], dt_h[i], snew_h[i]);

	return 1;
}

//Runge Kutta Fehlberg step with adaptive time step
inline int asteroid::RKF_step(const int level, double dtlimit){

	int er;
	double dt = dt_h[level];
	int N = Nlevel_h[level];
	int N0 = Nlevel_h[0];

	int stopFlag = 0;
	for(int l = 0; l <= level; ++l){
		if(stop_h[l] == 1){
			stopFlag = 1;
		}
	}

	if(stopFlag == 1){
		dtlimit = dt * dts;
	}

	for(int S = 0; S < RKFn; ++S){
		// ----------------------------------------------------------------------------
		//Update the Chebyshev coefficients if necessary
		er = update_Chebyshev(time_h[level] + RKFc_h[S] * dt);
		if(er <= 0){
			return 0;
		}
		update_perturbers(time_h[level] + RKFc_h[S] * dt);
		// ----------------------------------------------------------------------------

		for(int jj = 0; jj < N; ++jj){


			int i = jj;

			if(level > 0){
				i = index_h[(level - 1) * N0 + jj];
			}

			double xt = x_h[i];
			double yt = y_h[i];
			double zt = z_h[i];

			double vxt = vx_h[i];
			double vyt = vy_h[i];
			double vzt = vz_h[i];
			for(int s = 0; s < S; ++s){
				double dtaa = dt * RKFa_h[S * RKFn + s];
				int ii = jj + s * N;
				xt  += dtaa * kx_h[ii];
				yt  += dtaa * ky_h[ii];
				zt  += dtaa * kz_h[ii];
				vxt += dtaa * kvx_h[ii];
				vyt += dtaa * kvy_h[ii];
				vzt += dtaa * kvz_h[ii];

			}
			int ii = jj + S * N;
			kx_h[ii] = vxt;
			ky_h[ii] = vyt;
			kz_h[ii] = vzt;

			// ----------------------------------------------------------------------------
			//compute forces

			double ax = 0.0;
			double ay = 0.0;
			double az = 0.0;

			//heliocentric coordinates
			double xih = xt - xTable_h[10];
			double yih = yt - yTable_h[10];
			double zih = zt - zTable_h[10];

			double vxih = vxt - vxTable_h[10];
			double vyih = vyt - vyTable_h[10];
			double vzih = vzt - vzTable_h[10];

			//r is used in multiple forces, so reuse it
			double rsq = xih * xih + yih * yih + zih * zih;
			double r = sqrt(rsq);

//			if(cometFlag > 0 && S == 0){
//				if(i == 0){
//					Tsave_h[timeStep_h[level]] = time_h[level];
//				}
//				Rsave_h[i * Rbuffersize + timeStep_h[level]] = r;
//			}

			//Earth centric coordinates
			double xiE = xt - xTable_h[2];
			double yiE = yt - yTable_h[2];
			double ziE = zt - zTable_h[2];

			if(useNonGrav == 1){
				NonGrav(xih, yih, zih, vxih, vyih, vzih, A1_h[i], A2_h[i], A3_h[i], r, ax, ay, az);
			}
			if(useGR == 1){
				GR(xih, yih, zih, vxih, vyih, vzih, r, ax, ay, az, GM_h[10]);
			}
			if(useJ2 == 1){
				J2(xiE, yiE, ziE, ax, ay, az, GM_h[2]);
			}
			Gravity(xt, yt, zt, xTable_h, yTable_h, zTable_h, ax, ay, az, i);

			kvx_h[ii] = ax;
			kvy_h[ii] = ay;
			kvz_h[ii] = az;
			// ----------------------------------------------------------------------------
		}
	}

	//update
	for(int jj = 0; jj < N; ++jj){
		int i = jj;

		if(level > 0){
			i = index_h[(level - 1) * N0 + jj];
		}

		dx_h[i] = 0.0;
		dy_h[i] = 0.0;
		dz_h[i] = 0.0;

		dvx_h[i] = 0.0;
		dvy_h[i] = 0.0;
		dvz_h[i] = 0.0;

		for(int S = 0; S < RKFn; ++S){
			double dtb = dt * RKFb_h[S];
			int ii = jj + S * N;
			dx_h[i] += dtb * kx_h[ii];
			dy_h[i] += dtb * ky_h[ii];
			dz_h[i] += dtb * kz_h[ii];

			dvx_h[i] += dtb * kvx_h[ii];
			dvy_h[i] += dtb * kvy_h[ii];
			dvz_h[i] += dtb * kvz_h[ii];
//if(i < 2) printf("dx %d %d %g %g %g %g %g %g\n", S, i, dx_h[i], dy_h[i], dz_h[i], RKFb_h[S], kx_h[i + S * N], dt);
		}


		//compute integration error


		//Hairer
/*
		double scalex  = atol + fmax(fabs(x_h[i]), fabs(x_h[i] + dx_h[i])) * rtol;
		double scaley  = atol + fmax(fabs(y_h[i]), fabs(y_h[i] + dy_h[i])) * rtol;
		double scalez  = atol + fmax(fabs(z_h[i]), fabs(z_h[i] + dz_h[i])) * rtol;

		double scalevx = atol + fmax(fabs(vx_h[i]), fabs(vx_h[i] + dvx_h[i])) * rtol;
		double scalevy = atol + fmax(fabs(vy_h[i]), fabs(vy_h[i] + dvy_h[i])) * rtol;
		double scalevz = atol + fmax(fabs(vz_h[i]), fabs(vz_h[i] + dvz_h[i])) * rtol;
*/

		double scalex  = atol + fabs(x_h[i]) * rtol;
		double scaley  = atol + fabs(y_h[i]) * rtol;
		double scalez  = atol + fabs(z_h[i]) * rtol;

		double scalevx = atol + fabs(vx_h[i]) * rtol;
		double scalevy = atol + fabs(vy_h[i]) * rtol;
		double scalevz = atol + fabs(vz_h[i]) * rtol;
		//error estimation
		double errorkx = 0.0;
		double errorky = 0.0;
		double errorkz = 0.0;
		
		double errorkvx = 0.0;
		double errorkvy = 0.0;
		double errorkvz = 0.0;

		for(int S = 0; S < RKFn; ++S){
			double f = (RKFb_h[S] - RKFbb_h[S]) * dt;
			int ii = jj + S * N;
			errorkx += f * kx_h[ii];
			errorky += f * ky_h[ii];
			errorkz += f * kz_h[ii];

			errorkvx += f * kvx_h[ii];
			errorkvy += f * kvy_h[ii];
			errorkvz += f * kvz_h[ii];
//printf("error %d %d %.20g %.20g %.20g %.20g %.20g %.20g %.20g\n", i, S, f, kx_h[ii], ky_h[ii], kz_h[ii], kvx_h[ii], kvy_h[ii], kvz_h[ii]);
//printf("error %d %d %.20g %.20g %.20g %.20g %.20g %.20g %.20g\n", i, S, f, errorkx, errorky, errorkz, errorkvx, errorkvy, errorkvz);
		}

		double errork = 0.0;
		errork += errorkx * errorkx / (scalex * scalex);
		errork += errorky * errorky / (scaley * scaley);
		errork += errorkz * errorkz / (scalez * scalez);
		errork += errorkvx * errorkvx / (scalevx * scalevx);
		errork += errorkvy * errorkvy / (scalevy * scalevy);
		errork += errorkvz * errorkvz / (scalevz * scalevz);

		errork = sqrt(errork / 6.0);	//6 is the number of dimensions

		double s = pow( 1.0  / errork, RKF_ee);

		s = (RKF_fac * s > RKF_facmin) ? RKF_fac * s : RKF_facmin;
		s = (RKF_facmax < s) ? RKF_facmax : s;

//printf("snew %d %d %g %g %g \n", level, i, s, dt, s * dt);

		if(s * dt * dts >= dtlimit || level >= nL - 1){
			snew_h[i] = s;
		}
		else{
			snew_h[i] = 1.0e6;	//mark body for higher level integration
			time_h[level + 1] = time_h[level];
		}
	}


	double snew = 1.0e6;
	//Find minimum s factor
	for(int jj = 0; jj < N; ++jj){
		int i = jj;

		if(level > 0){
			i = index_h[(level - 1) * N0 + jj];
		}

		double s = snew_h[i];

		snew = (snew < s) ? snew : s;
	}
	if(snew == 1.0e6){
		snew = 1.0;
	}

	snewlevel_h[level] = snew;

	if(snew >= RKF_fac){
		//accept step
		for(int jj = 0; jj < N; ++jj){
			int i = jj;

			if(level > 0){
				i = index_h[(level - 1) * N0 + jj];
			}

			if(snew_h[i] < 1.0e6){
				x_h[i] += dx_h[i];
				y_h[i] += dy_h[i];
				z_h[i] += dz_h[i];

				vx_h[i] += dvx_h[i];
				vy_h[i] += dvy_h[i];
				vz_h[i] += dvz_h[i];

				if(stopFlag == 0){
					if(dts < 0){
						dtmin_h[i] = dt > dtmin_h[i] ? dt : dtmin_h[i];
					}
					else{
						dtmin_h[i] = dt < dtmin_h[i] ? dt : dtmin_h[i];
					}
				}
//printf("dtmin %d %d %d %g\n", i, level, stopFlag, dtmin_h[i]);

			}
			else{
				int j = Nlevel_h[level + 1]++;
				index_h[level * N0 + j] = i;
				//add to list
//printf("add to list %d %d %g\n", i, j, time_h[level + 1]);
			}

		}

		time_h[level] += dt;
		++timeStep_h[level];
	}
	else{
//printf("repeat level %d %g\n", level, snew);
	}


	return 1;
}

// Bulirsh-Stoer step with adaptive time step
inline int asteroid::BS_individual_step(const int i){

	int er;
	double dt = dt_h[i];

	double ax;
	double ay;
	double az;

	double xt;                      
	double yt;
	double zt;
		
	double vxt;
	double vyt;
	double vzt;

	double xp;                      
	double yp;
	double zp;
		
	double vxp;
	double vyp;
	double vzp;


	double scalex = atol + fabs(x_h[i]) * rtol;
	double scaley = atol + fabs(y_h[i]) * rtol;
	double scalez = atol + fabs(z_h[i]) * rtol;

	double scalevx = atol + fabs(vx_h[i]) * rtol;
	double scalevy = atol + fabs(vy_h[i]) * rtol;
	double scalevz = atol + fabs(vz_h[i]) * rtol;

	int f = 1;

	for(int n = 1; n <= 8; ++n){
		double dt2 = dt / (2.0 * n);
		double dt22 = dt2 * 2.0;

		// ----------------------------------------------------------------------------
		//Update the Chebyshev coefficients if necessary
		er = update_Chebyshev(time_h[i]);
		if(er <= 0){
			return 0;
		}
		update_perturbers(time_h[i]);
//printf("%d %d %g %g\n", n, 0, 0.0, dt);
		// ----------------------------------------------------------------------------

		ax = 0.0;
		ay = 0.0;
		az = 0.0;

		//compute forces
		// ----------------------------------------------------------------------------
		//heliocentric coordinates
		double xih = x_h[i] - xTable_h[10];
		double yih = y_h[i] - yTable_h[10];
		double zih = z_h[i] - zTable_h[10];

		double vxih = vx_h[i] - vxTable_h[10];
		double vyih = vy_h[i] - vyTable_h[10];
		double vzih = vz_h[i] - vzTable_h[10];

		//r is used in multiple forces, so reuse it
		double rsq = xih * xih + yih * yih + zih * zih;
		double r = sqrt(rsq);

		if(cometFlag > 0 && n == 1){
			if(i == 0){
				Tsave_h[timeStep_h[i]] = time_h[i];
			}
			Rsave_h[i * Rbuffersize + timeStep_h[i]] = r;
		}

		//Earth centric coordinates
		double xiE = x_h[i] - xTable_h[2];
		double yiE = y_h[i] - yTable_h[2];
		double ziE = z_h[i] - zTable_h[2];

		if(useNonGrav == 1){
			NonGrav(xih, yih, zih, vxih, vyih, vzih, A1_h[i], A2_h[i], A3_h[i], r, ax, ay, az);
		}
		if(useGR == 1){
			GR(xih, yih, zih, vxih, vyih, vzih, r, ax, ay, az, GM_h[10]);
		}
		if(useJ2 == 1){
			J2(xiE, yiE, ziE, ax, ay, az, GM_h[2]);
		}
		Gravity(x_h[i], y_h[i], z_h[i], xTable_h, yTable_h, zTable_h, ax, ay, az, i);
		// ----------------------------------------------------------------------------

		xp = x_h[i] + dt2 * vx_h[i];
		yp = y_h[i] + dt2 * vy_h[i];
		zp = z_h[i] + dt2 * vz_h[i];


		vxp = vx_h[i] + dt2 * ax;
		vyp = vy_h[i] + dt2 * ay;
		vzp = vz_h[i] + dt2 * az;

		// ----------------------------------------------------------------------------
		//Update the Chebyshev coefficients if necessary
		er = update_Chebyshev(time_h[i] + dt2);
		if(er <= 0){
			return 0;
		}
		update_perturbers(time_h[i] + dt2);
//printf("%d %d %g\n", n, 1, dt2 / dt);
		// ----------------------------------------------------------------------------

		ax = 0.0;
		ay = 0.0;
		az = 0.0;

		//compute forces
		// ----------------------------------------------------------------------------
		//heliocentric coordinates
		xih = xp - xTable_h[10];
		yih = yp - yTable_h[10];
		zih = zp - zTable_h[10];

		vxih = vxp - vxTable_h[10];
		vyih = vyp - vyTable_h[10];
		vzih = vzp - vzTable_h[10];

		//r is used in multiple forces, so reuse it
		rsq = xih * xih + yih * yih + zih * zih;
		r = sqrt(rsq);


		//Earth centric coordinates
		xiE = xp - xTable_h[2];
		yiE = yp - yTable_h[2];
		ziE = zp - zTable_h[2];

		if(useNonGrav == 1){
			NonGrav(xih, yih, zih, vxih, vyih, vzih, A1_h[i], A2_h[i], A3_h[i], r, ax, ay, az);
		}
		if(useGR == 1){
			GR(xih, yih, zih, vxih, vyih, vzih, r, ax, ay, az, GM_h[10]);
		}
		if(useJ2 == 1){
			J2(xiE, yiE, ziE, ax, ay, az, GM_h[2]);
		}
		Gravity(xp, yp, zp, xTable_h, yTable_h, zTable_h, ax, ay, az, i);
		// ----------------------------------------------------------------------------

		xt = x_h[i] + dt22 * vxp;
		yt = y_h[i] + dt22 * vyp;
		zt = z_h[i] + dt22 * vzp;

		vxt = vx_h[i] + dt22 * ax;
		vyt = vy_h[i] + dt22 * ay;
		vzt = vz_h[i] + dt22 * az;

		for(int m = 2; m <= n; ++m){

			// ----------------------------------------------------------------------------
			//Update the Chebyshev coefficients if necessary
			er = update_Chebyshev (time_h[i] + (m-1) * dt22);
			if(er <= 0){
				return 0;
			}
			update_perturbers(time_h[i] + (m-1) * dt22);
//printf("%d %d %g\n", n, m, (m-1) * dt22 / dt);
			// ----------------------------------------------------------------------------

			ax = 0.0;
			ay = 0.0;
			az = 0.0;

			//compute forces
			// ----------------------------------------------------------------------------
			//heliocentric coordinates
			xih = xt - xTable_h[10];
			yih = yt - yTable_h[10];
			zih = zt - zTable_h[10];

			vxih = vxt - vxTable_h[10];
			vyih = vyt - vyTable_h[10];
			vzih = vzt - vzTable_h[10];

			//r is used in multiple forces, so reuse it
			rsq = xih * xih + yih * yih + zih * zih;
			r = sqrt(rsq);

			//Earth centric coordinates
			xiE = xt - xTable_h[2];
			yiE = yt - yTable_h[2];
			ziE = zt - zTable_h[2];

			if(useNonGrav == 1){
				NonGrav(xih, yih, zih, vxih, vyih, vzih, A1_h[i], A2_h[i], A3_h[i], r, ax, ay, az);
			}
			if(useGR == 1){
				GR(xih, yih, zih, vxih, vyih, vzih, r, ax, ay, az, GM_h[10]);
			}
			if(useJ2 == 1){
				J2(xiE, yiE, ziE, ax, ay, az, GM_h[2]);
			}
			Gravity(xt, yt, zt, xTable_h, yTable_h, zTable_h, ax, ay, az, i);
			// ----------------------------------------------------------------------------

			xp += dt22 * vxt;
			yp += dt22 * vyt;
			zp += dt22 * vzt;

			vxp += dt22 * ax;
			vyp += dt22 * ay;
			vzp += dt22 * az;

			// ----------------------------------------------------------------------------
			//Update the Chebyshev coefficients if necessary
			er = update_Chebyshev(time_h[i] + (m-1) * dt22 + dt2);
			if(er <= 0){
				return 0;
			}
			update_perturbers(time_h[i] + (m-1) * dt22 + dt2);
//printf("%d %d %g\n", n, m, ((m-1) * dt22 + dt2) / dt);
			// ----------------------------------------------------------------------------

			ax = 0.0;
			ay = 0.0;
			az = 0.0;

			//compute forces
			// ----------------------------------------------------------------------------
			//heliocentric coordinates
			xih = xp - xTable_h[10];
			yih = yp - yTable_h[10];
			zih = zp - zTable_h[10];

			vxih = vxp - vxTable_h[10];
			vyih = vyp - vyTable_h[10];
			vzih = vzp - vzTable_h[10];

			//r is used in multiple forces, so reuse it
			rsq = xih * xih + yih * yih + zih * zih;
			r = sqrt(rsq);

			//Earth centric coordinates
			xiE = xp - xTable_h[2];
			yiE = yp - yTable_h[2];
			ziE = zp - zTable_h[2];

			if(useNonGrav == 1){
				NonGrav(xih, yih, zih, vxih, vyih, vzih, A1_h[i], A2_h[i], A3_h[i], r, ax, ay, az);
			}
			if(useGR == 1){
				GR(xih, yih, zih, vxih, vyih, vzih, r, ax, ay, az, GM_h[10]);
			}
			if(useJ2 == 1){
				J2(xiE, yiE, ziE, ax, ay, az, GM_h[2]);
			}
			Gravity(xp, yp, zp, xTable_h, yTable_h, zTable_h, ax, ay, az, i);
			// ----------------------------------------------------------------------------

			xt += dt22 * vxp;
			yt += dt22 * vyp;
			zt += dt22 * vzp;

			vxt += dt22 * ax;
			vyt += dt22 * ay;
			vzt += dt22 * az;

		} // end if m loop

		// ----------------------------------------------------------------------------
		//Update the Chebyshev coefficients if necessary
		er = update_Chebyshev(time_h[i] + dt);
		if(er <= 0){
			return 0;
		}
		update_perturbers(time_h[i] + dt);
//printf("%d %d %g\n", n, n+1, 1.0);
		// ----------------------------------------------------------------------------

		ax = 0.0;
		ay = 0.0;
		az = 0.0;

		//compute forces
		// ----------------------------------------------------------------------------
		//heliocentric coordinates
		xih = xt - xTable_h[10];
		yih = yt - yTable_h[10];
		zih = zt - zTable_h[10];

		vxih = vxt - vxTable_h[10];
		vyih = vyt - vyTable_h[10];
		vzih = vzt - vzTable_h[10];

		//r is used in multiple forces, so reuse it
		rsq = xih * xih + yih * yih + zih * zih;
		r = sqrt(rsq);

		//Earth centric coordinates
		xiE = xt - xTable_h[2];
		yiE = yt - yTable_h[2];
		ziE = zt - zTable_h[2];

		if(useNonGrav == 1){
			NonGrav(xih, yih, zih, vxih, vyih, vzih, A1_h[i], A2_h[i], A3_h[i], r, ax, ay, az);
		}
		if(useGR == 1){
			GR(xih, yih, zih, vxih, vyih, vzih, r, ax, ay, az, GM_h[10]);
		}
		if(useJ2 == 1){
			J2(xiE, yiE, ziE, ax, ay, az, GM_h[2]);
		}
		Gravity(xt, yt, zt, xTable_h, yTable_h, zTable_h, ax, ay, az, i);
		// ----------------------------------------------------------------------------

		xp += dt2 * vxt;
		yp += dt2 * vyt;
		zp += dt2 * vzt;

		vxp += dt2 * ax;
		vyp += dt2 * ay;
		vzp += dt2 * az;

		dx_h[i * 8 + (n-1)] = 0.5 * (xt + xp);
		dy_h[i * 8 + (n-1)] = 0.5 * (yt + yp);
		dz_h[i * 8 + (n-1)] = 0.5 * (zt + zp);

		dvx_h[i * 8 + (n-1)] = 0.5 * (vxt + vxp);
		dvy_h[i * 8 + (n-1)] = 0.5 * (vyt + vyp);
		dvz_h[i * 8 + (n-1)] = 0.5 * (vzt + vzp);


		//Extrapolation step

		for(int j = n - 1; j >= 1; --j){
			double t0 = BSt0_h[(n-1) * 8 + (j-1)];
			double t1 = t0 * BSddt_h[j];
			double t2 = t0 * BSddt_h[n-1];

			dx_h[i * 8 + (j-1)] = (t1 * dx_h[i * 8 + j]) - (t2 * dx_h[i * 8 + (j-1)]);
			dy_h[i * 8 + (j-1)] = (t1 * dy_h[i * 8 + j]) - (t2 * dy_h[i * 8 + (j-1)]);
			dz_h[i * 8 + (j-1)] = (t1 * dz_h[i * 8 + j]) - (t2 * dz_h[i * 8 + (j-1)]);

			dvx_h[i * 8 + (j-1)] = (t1 * dvx_h[i * 8 + j]) - (t2 * dvx_h[i * 8 + (j-1)]);
			dvy_h[i * 8 + (j-1)] = (t1 * dvy_h[i * 8 + j]) - (t2 * dvy_h[i * 8 + (j-1)]);
			dvz_h[i * 8 + (j-1)] = (t1 * dvz_h[i * 8 + j]) - (t2 * dvz_h[i * 8 + (j-1)]);

		}

		double error = 0.0;
		double error1 = 0.0;

		error1 = dx_h[i * 8] / scalex;
		error += error1 * error1;

		error1 = dy_h[i * 8] / scaley;
		error += error1 * error1;

		error1 = dz_h[i * 8] / scalez;
		error += error1 * error1;

		error1 = dvx_h[i * 8] / scalevx;
		error += error1 * error1;

		error1 = dvy_h[i * 8] / scalevy;
		error += error1 * error1;

		error1 = dvz_h[i * 8] / scalevz;
		error += error1 * error1;


		error = sqrt(error / 6.0);	//6 is the number of dimensions

		if(error < 1.0){
//printf("accept %lld %.20g %g\n", timeStep, time, dt);

			//accept step
			snew_h[i] = 1.0;

			xt = dx_h[i * 8];
			yt = dy_h[i * 8];
			zt = dz_h[i * 8];

			vxt = dvx_h[i * 8];
			vyt = dvy_h[i * 8];
			vzt = dvz_h[i * 8];

			for(int j = 1; j < n; ++j){
				xt += dx_h[i * 8 + j];
				yt += dy_h[i * 8 + j];
				zt += dz_h[i * 8 + j];

				vxt += dvx_h[i * 8 + j];
				vyt += dvy_h[i * 8 + j];
				vzt += dvz_h[i * 8 + j];


			}

			time_h[i] += dt;
			++timeStep_h[i];

			if(n >= 8){
				snew_h[i] = 0.55;
//printf("reduce time step %g\n", dt);
			}
			if(n < 7){
				snew_h[i] = 1.3;
//printf("increase time step %g\n", dt);
			}

			x_h[i] = xt;
			y_h[i] = yt;
			z_h[i] = zt;

			vx_h[i] = vxt;
			vy_h[i] = vyt;
			vz_h[i] = vzt;

			f = 0;
			break; //break n loop
		}
	} //end of n loop
	if(f == 1){
		snew_h[i] = 0.5;
//printf("repeat time step %g\n", dt);
	}
	return 1;
}
	

// Bulirsh-Stoer step with adaptive time step
inline int asteroid::BS_step(const int level, double dtlimit){

	int er;
	double dt = dt_h[level];
	int N = Nlevel_h[level];
	int N0 = Nlevel_h[0];

	int stopFlag = 0;
	for(int l = 0; l <= level; ++l){
		if(stop_h[l] == 1){
			stopFlag = 1;
		}               

	}       
	if(stopFlag == 1){
		dtlimit = dt * dts;
	}


	for(int jj = 0; jj < N; ++jj){

		int i = jj;

		if(level > 0){
			i = index_h[(level - 1) * N0 + jj];
		}


		scalex_h[i] = atol + fabs(x_h[i]) * rtol;
		scaley_h[i] = atol + fabs(y_h[i]) * rtol;
		scalez_h[i] = atol + fabs(z_h[i]) * rtol;

		scalevx_h[i] = atol + fabs(vx_h[i]) * rtol;
		scalevy_h[i] = atol + fabs(vy_h[i]) * rtol;
		scalevz_h[i] = atol + fabs(vz_h[i]) * rtol;

		snew_h[i] = -1000.0;
	}


	int f = 1;

	for(int n = 1; n <= 8; ++n){
		double dt2 = dt / (2.0 * n);
		double dt22 = dt2 * 2.0;

		// ----------------------------------------------------------------------------
		//Update the Chebyshev coefficients if necessary
		er = update_Chebyshev(time_h[level]);
		if(er <= 0){
			return 0;
		}
		update_perturbers(time_h[level]);
//printf("%d %d %g %g\n", n, 0, 0.0, dt);
		// ----------------------------------------------------------------------------

		for(int jj = 0; jj < N; ++jj){
			int i = jj;

			if(level > 0){
				i = index_h[(level - 1) * N0 + jj];
			}

			if(snew_h[i] > 0.0) continue;

			double ax = 0.0;
			double ay = 0.0;
			double az = 0.0;

			//compute forces
			// ----------------------------------------------------------------------------
			//heliocentric coordinates
			double xih = x_h[i] - xTable_h[10];
			double yih = y_h[i] - yTable_h[10];
			double zih = z_h[i] - zTable_h[10];

			double vxih = vx_h[i] - vxTable_h[10];
			double vyih = vy_h[i] - vyTable_h[10];
			double vzih = vz_h[i] - vzTable_h[10];

			//r is used in multiple forces, so reuse it
			double rsq = xih * xih + yih * yih + zih * zih;
			double r = sqrt(rsq);

//			if(cometFlag > 0 && n == 1){
//				if(i == 0){
//					Tsave_h[timeStep_h[level]] = time_h[level];
//				}
//				Rsave_h[i * Rbuffersize + timeStep_h[level]] = r;
//			}

			//Earth centric coordinates
			double xiE = x_h[i] - xTable_h[2];
			double yiE = y_h[i] - yTable_h[2];
			double ziE = z_h[i] - zTable_h[2];

			if(useNonGrav == 1){
				NonGrav(xih, yih, zih, vxih, vyih, vzih, A1_h[i], A2_h[i], A3_h[i], r, ax, ay, az);
			}
			if(useGR == 1){
				GR(xih, yih, zih, vxih, vyih, vzih, r, ax, ay, az, GM_h[10]);
			}
			if(useJ2 == 1){
				J2(xiE, yiE, ziE, ax, ay, az, GM_h[2]);
			}
			Gravity(x_h[i], y_h[i], z_h[i], xTable_h, yTable_h, zTable_h, ax, ay, az, i);
			// ----------------------------------------------------------------------------

			xp_h[i] = x_h[i] + dt2 * vx_h[i];
			yp_h[i] = y_h[i] + dt2 * vy_h[i];
			zp_h[i] = z_h[i] + dt2 * vz_h[i];

			vxp_h[i] = vx_h[i] + dt2 * ax;
			vyp_h[i] = vy_h[i] + dt2 * ay;
			vzp_h[i] = vz_h[i] + dt2 * az;
		}

		// ----------------------------------------------------------------------------
		//Update the Chebyshev coefficients if necessary
		er = update_Chebyshev(time_h[level] + dt2);
		if(er <= 0){
			return 0;
		}
		update_perturbers(time_h[level] + dt2);
//printf("%d %d %g\n", n, 1, dt2 / dt);
		// ----------------------------------------------------------------------------
		for(int jj = 0; jj < N; ++jj){
			int i = jj;

			if(level > 0){
				i = index_h[(level - 1) * N0 + jj];
			}

			if(snew_h[i] > 0.0) continue;

			double ax = 0.0;
			double ay = 0.0;
			double az = 0.0;

			//compute forces
			// ----------------------------------------------------------------------------
			//heliocentric coordinates
			double xih = xp_h[i] - xTable_h[10];
			double yih = yp_h[i] - yTable_h[10];
			double zih = zp_h[i] - zTable_h[10];

			double vxih = vxp_h[i] - vxTable_h[10];
			double vyih = vyp_h[i] - vyTable_h[10];
			double vzih = vzp_h[i] - vzTable_h[10];

			//r is used in multiple forces, so reuse it
			double rsq = xih * xih + yih * yih + zih * zih;
			double r = sqrt(rsq);


			//Earth centric coordinates
			double xiE = xp_h[i] - xTable_h[2];
			double yiE = yp_h[i] - yTable_h[2];
			double ziE = zp_h[i] - zTable_h[2];

			if(useNonGrav == 1){
				NonGrav(xih, yih, zih, vxih, vyih, vzih, A1_h[i], A2_h[i], A3_h[i], r, ax, ay, az);
			}
			if(useGR == 1){
				GR(xih, yih, zih, vxih, vyih, vzih, r, ax, ay, az, GM_h[10]);
			}
			if(useJ2 == 1){
				J2(xiE, yiE, ziE, ax, ay, az, GM_h[2]);
			}
			Gravity(xp_h[i], yp_h[i], zp_h[i], xTable_h, yTable_h, zTable_h, ax, ay, az, i);
			// ----------------------------------------------------------------------------

			xt_h[i] = x_h[i] + dt22 * vxp_h[i];
			yt_h[i] = y_h[i] + dt22 * vyp_h[i];
			zt_h[i] = z_h[i] + dt22 * vzp_h[i];

			vxt_h[i] = vx_h[i] + dt22 * ax;
			vyt_h[i] = vy_h[i] + dt22 * ay;
			vzt_h[i] = vz_h[i] + dt22 * az;
		}
		for(int m = 2; m <= n; ++m){

			// ----------------------------------------------------------------------------
			//Update the Chebyshev coefficients if necessary
			er = update_Chebyshev(time_h[level] + (m-1) * dt22);
			if(er <= 0){
				return 0;
			}
			update_perturbers(time_h[level] + (m-1) * dt22);
//printf("%d %d %g\n", n, m, (m-1) * dt22 / dt);
			// ----------------------------------------------------------------------------
			for(int jj = 0; jj < N; ++jj){
				int i = jj;

				if(level > 0){
					i = index_h[(level - 1) * N0 + jj];
				}

				if(snew_h[i] > 0.0) continue;

				double ax = 0.0;
				double ay = 0.0;
				double az = 0.0;

				//compute forces
				// ----------------------------------------------------------------------------
				//heliocentric coordinates
				double xih = xt_h[i] - xTable_h[10];
				double yih = yt_h[i] - yTable_h[10];
				double zih = zt_h[i] - zTable_h[10];

				double vxih = vxt_h[i] - vxTable_h[10];
				double vyih = vyt_h[i] - vyTable_h[10];
				double vzih = vzt_h[i] - vzTable_h[10];

				//r is used in multiple forces, so reuse it
				double rsq = xih * xih + yih * yih + zih * zih;
				double r = sqrt(rsq);

				//Earth centric coordinates
				double xiE = xt_h[i] - xTable_h[2];
				double yiE = yt_h[i] - yTable_h[2];
				double ziE = zt_h[i] - zTable_h[2];

				if(useNonGrav == 1){
					NonGrav(xih, yih, zih, vxih, vyih, vzih, A1_h[i], A2_h[i], A3_h[i], r, ax, ay, az);
				}
				if(useGR == 1){
					GR(xih, yih, zih, vxih, vyih, vzih, r, ax, ay, az, GM_h[10]);
				}
				if(useJ2 == 1){
					J2(xiE, yiE, ziE, ax, ay, az, GM_h[2]);
				}
				Gravity(xt_h[i], yt_h[i], zt_h[i], xTable_h, yTable_h, zTable_h, ax, ay, az, i);
				// ----------------------------------------------------------------------------

				xp_h[i] += dt22 * vxt_h[i];
				yp_h[i] += dt22 * vyt_h[i];
				zp_h[i] += dt22 * vzt_h[i];

				vxp_h[i] += dt22 * ax;
				vyp_h[i] += dt22 * ay;
				vzp_h[i] += dt22 * az;
			}
			// ----------------------------------------------------------------------------
			//Update the Chebyshev coefficients if necessary
			er = update_Chebyshev(time_h[level] + (m-1) * dt22 + dt2);
			if(er <= 0){
				return 0;
			}
			update_perturbers(time_h[level] + (m-1) * dt22 + dt2);
//printf("%d %d %g\n", n, m, ((m-1) * dt22 + dt2) / dt);
			// ----------------------------------------------------------------------------
			for(int jj = 0; jj < N; ++jj){
				int i = jj;

				if(level > 0){
					i = index_h[(level - 1) * N0 + jj];
				}

				if(snew_h[i] > 0.0) continue;

				double ax = 0.0;
				double ay = 0.0;
				double az = 0.0;

				//compute forces
				// ----------------------------------------------------------------------------
				//heliocentric coordinates
				double xih = xp_h[i] - xTable_h[10];
				double yih = yp_h[i] - yTable_h[10];
				double zih = zp_h[i] - zTable_h[10];

				double vxih = vxp_h[i] - vxTable_h[10];
				double vyih = vyp_h[i] - vyTable_h[10];
				double vzih = vzp_h[i] - vzTable_h[10];

				//r is used in multiple forces, so reuse it
				double rsq = xih * xih + yih * yih + zih * zih;
				double r = sqrt(rsq);

				//Earth centric coordinates
				double xiE = xp_h[i] - xTable_h[2];
				double yiE = yp_h[i] - yTable_h[2];
				double ziE = zp_h[i] - zTable_h[2];

				if(useNonGrav == 1){
					NonGrav(xih, yih, zih, vxih, vyih, vzih, A1_h[i], A2_h[i], A3_h[i], r, ax, ay, az);
				}
				if(useGR == 1){
					GR(xih, yih, zih, vxih, vyih, vzih, r, ax, ay, az, GM_h[10]);
				}
				if(useJ2 == 1){
					J2(xiE, yiE, ziE, ax, ay, az, GM_h[2]);
				}
				Gravity(xp_h[i], yp_h[i], zp_h[i], xTable_h, yTable_h, zTable_h, ax, ay, az, i);
				// ----------------------------------------------------------------------------

				xt_h[i] += dt22 * vxp_h[i];
				yt_h[i] += dt22 * vyp_h[i];
				zt_h[i] += dt22 * vzp_h[i];

				vxt_h[i] += dt22 * ax;
				vyt_h[i] += dt22 * ay;
				vzt_h[i] += dt22 * az;
			}
		} // end if m loop

		// ----------------------------------------------------------------------------
		//Update the Chebyshev coefficients if necessary
		er = update_Chebyshev(time_h[level] + dt);
		if(er <= 0){
			return 0;
		}
		update_perturbers(time_h[level] + dt);
//printf("%d %d %g\n", n, n+1, 1.0);
		// ----------------------------------------------------------------------------
		for(int jj = 0; jj < N; ++jj){
			int i = jj;

			if(level > 0){
				i = index_h[(level - 1) * N0 + jj];
			}

			if(snew_h[i] > 0.0) continue;

			double ax = 0.0;
			double ay = 0.0;
			double az = 0.0;

			//compute forces
			// ----------------------------------------------------------------------------
			//heliocentric coordinates
			double xih = xt_h[i] - xTable_h[10];
			double yih = yt_h[i] - yTable_h[10];
			double zih = zt_h[i] - zTable_h[10];

			double vxih = vxt_h[i] - vxTable_h[10];
			double vyih = vyt_h[i] - vyTable_h[10];
			double vzih = vzt_h[i] - vzTable_h[10];

			//r is used in multiple forces, so reuse it
			double rsq = xih * xih + yih * yih + zih * zih;
			double r = sqrt(rsq);

			//Earth centric coordinates
			double xiE = xt_h[i] - xTable_h[2];
			double yiE = yt_h[i] - yTable_h[2];
			double ziE = zt_h[i] - zTable_h[2];

			if(useNonGrav == 1){
				NonGrav(xih, yih, zih, vxih, vyih, vzih, A1_h[i], A2_h[i], A3_h[i], r, ax, ay, az);
			}
			if(useGR == 1){
				GR(xih, yih, zih, vxih, vyih, vzih, r, ax, ay, az, GM_h[10]);
			}
			if(useJ2 == 1){
				J2(xiE, yiE, ziE, ax, ay, az, GM_h[2]);
			}
			Gravity(xt_h[i], yt_h[i], zt_h[i], xTable_h, yTable_h, zTable_h, ax, ay, az, i);
			// ----------------------------------------------------------------------------

			xp_h[i] += dt2 * vxt_h[i];
			yp_h[i] += dt2 * vyt_h[i];
			zp_h[i] += dt2 * vzt_h[i];

			vxp_h[i] += dt2 * ax;
			vyp_h[i] += dt2 * ay;
			vzp_h[i] += dt2 * az;

			dx_h[i * 8 + (n-1)] = 0.5 * (xt_h[i] + xp_h[i]);
			dy_h[i * 8 + (n-1)] = 0.5 * (yt_h[i] + yp_h[i]);
			dz_h[i * 8 + (n-1)] = 0.5 * (zt_h[i] + zp_h[i]);

			dvx_h[i * 8 + (n-1)] = 0.5 * (vxt_h[i] + vxp_h[i]);
			dvy_h[i * 8 + (n-1)] = 0.5 * (vyt_h[i] + vyp_h[i]);
			dvz_h[i * 8 + (n-1)] = 0.5 * (vzt_h[i] + vzp_h[i]);
		}


		double errorMax = 0.0;
		//Extrapolation step
		for(int jj = 0; jj < N; ++jj){
			int i = jj;

			if(level > 0){
				i = index_h[(level - 1) * N0 + jj];
			}

			if(snew_h[i] > 0.0) continue;

			for(int j = n - 1; j >= 1; --j){
				double t0 = BSt0_h[(n-1) * 8 + (j-1)];
				double t1 = t0 * BSddt_h[j];
				double t2 = t0 * BSddt_h[n-1];

				dx_h[i * 8 + (j-1)] = (t1 * dx_h[i * 8 + j]) - (t2 * dx_h[i * 8 + (j-1)]);
				dy_h[i * 8 + (j-1)] = (t1 * dy_h[i * 8 + j]) - (t2 * dy_h[i * 8 + (j-1)]);
				dz_h[i * 8 + (j-1)] = (t1 * dz_h[i * 8 + j]) - (t2 * dz_h[i * 8 + (j-1)]);

				dvx_h[i * 8 + (j-1)] = (t1 * dvx_h[i * 8 + j]) - (t2 * dvx_h[i * 8 + (j-1)]);
				dvy_h[i * 8 + (j-1)] = (t1 * dvy_h[i * 8 + j]) - (t2 * dvy_h[i * 8 + (j-1)]);
				dvz_h[i * 8 + (j-1)] = (t1 * dvz_h[i * 8 + j]) - (t2 * dvz_h[i * 8 + (j-1)]);

			}

			double error = 0.0;
			double error1 = 0.0;

			error1 = dx_h[i * 8] / scalex_h[i];
			error += error1 * error1;

			error1 = dy_h[i * 8] / scaley_h[i];
			error += error1 * error1;

			error1 = dz_h[i * 8] / scalez_h[i];
			error += error1 * error1;

			error1 = dvx_h[i * 8] / scalevx_h[i];
			error += error1 * error1;

			error1 = dvy_h[i * 8] / scalevy_h[i];
			error += error1 * error1;

			error1 = dvz_h[i * 8] / scalevz_h[i];
			error += error1 * error1;


			error = sqrt(error / 6.0);	//6 is the number of dimensions

			errorMax = error > errorMax ? error : errorMax;


			if(error < 1.0){
				//accept step

				xt_h[i] = dx_h[i * 8];
				yt_h[i] = dy_h[i * 8];
				zt_h[i] = dz_h[i * 8];

				vxt_h[i] = dvx_h[i * 8];
				vyt_h[i] = dvy_h[i * 8];
				vzt_h[i] = dvz_h[i * 8];
	
				for(int j = 1; j < n; ++j){
					xt_h[i] += dx_h[i * 8 + j];
					yt_h[i] += dy_h[i * 8 + j];
					zt_h[i] += dz_h[i * 8 + j];

					vxt_h[i] += dvx_h[i * 8 + j];
					vyt_h[i] += dvy_h[i * 8 + j];
					vzt_h[i] += dvz_h[i * 8 + j];
				}

				snew_h[i] = 1.0;
				if(n >= 8){
					snew_h[i] = 0.55;
				}
				if(n < 7){
					snew_h[i] = 1.3;
				}
//printf("Accect %d %d %d %g %g\n", level, i, n, dt, snew_h[i]);
			}
			else if(n == 8){
				//repeat step
				snew_h[i] = 0.5;
//printf("Repeat %d %d %d %g %g\n", level, i, n, dt, snew_h[i]);
			}

		}

		if(errorMax < 1.0){
			f = 0;
			break;
		}
	}// end if n loop
	

	for(int jj = 0; jj < N; ++jj){
		int i = jj;

		if(level > 0){
			i = index_h[(level - 1) * N0 + jj];
		}
		double s = snew_h[i];

		if(s * dt * dts >= dtlimit || level >= nL - 1){
			//snew_h[i] = s;
		}
		else{
			snew_h[i] = 1.0e6;      //mark body for higher level integration
			time_h[level + 1] = time_h[level];
		}
	}


	double snew = 1.0e6;
	//Find maximum error
	for(int jj = 0; jj < N; ++jj){
		int i = jj;

		if(level > 0){
			i = index_h[(level - 1) * N0 + jj];
		}


		double s = snew_h[i];

		snew = (snew < s) ? snew : s;
	}
	if(snew == 1.0e6){
		snew = 1.0;
	}

//printf("snewMin %g %g\n", snew, dtlimit);
	snewlevel_h[level] = snew;

	if(snew > 0.5){
//printf("accept %lld %.20g %g\n", timeStep_h[level, time_h[level], dt);

		//accept step
		for(int jj = 0; jj < N; ++jj){
			int i = jj;

			if(level > 0){
				i = index_h[(level - 1) * N0 + jj];
			}

			if(snew_h[i] < 1.0e6){
				x_h[i] = xt_h[i];
				y_h[i] = yt_h[i];
				z_h[i] = zt_h[i];

				vx_h[i] = vxt_h[i];
				vy_h[i] = vyt_h[i];
				vz_h[i] = vzt_h[i];

//printf("update %d %g\n", i, time_h[level] + dt);

				if(stopFlag == 0){
					if(dts < 0){
						dtmin_h[i] = dt > dtmin_h[i] ? dt : dtmin_h[i];
					}
					else{
						dtmin_h[i] = dt < dtmin_h[i] ? dt : dtmin_h[i];
					}
				}
//printf("dtmin %d %d %d %g\n", i, level, stopFlag, dtmin_h[i]);

			}
			else{
				int j = Nlevel_h[level + 1]++;
				index_h[level * N0 + j] = i;
				//add to list
//printf("add to list %d %d %g\n", i, j, time_h[level + 1]);
			}


		}

		time_h[level] += dt;
		++timeStep_h[level];
	}
	else{
//printf("repeat level %d %g\n", level, snew);
	}


	return 1;
}
	
int asteroid::loop_individual(){
	int er;

	//If needed, convert from heliocentric equatorial coordinates to barycentric equatorial coordinates
	if(ICheliocentric == 1){
		printf("Convert heliocentric to barycentric coordinates\n");
		er = HelioToBary(x_h, y_h, z_h, vx_h, vy_h, vz_h);
		if(er <= 0){
			return 0;
		}
		printf("Convert heliocentric to barycentric coordinates OK\n");
	} 


	//At this point, the initial conditions coordinates are cartesian barycentric equatorial
	//The integration is also done in cartesian barycentric equatorial coordinates

	if(outBinary == 0){
		outputFile = fopen(outputFilename, "w");
	}
	else{
		outputFile = fopen(outputFilename, "wb");
	}
	
	if(printdt == 1){
		dtFile = fopen(dtFilename, "w");
	}
	printf("Start integration %.20g\n", timeStart + time_reference);

	for(int i = 0; i < N; ++i){
		dt_h[i] = dt;
		dtsave_h[i] = dt;
		dtmin_h[i] = dt;
		timeStep_h[i] = 0ll;
		time_h[i] = time;
	}
	dtminlevel_h[0] = dt;

	if(time_reference + time >= outStart){
		er = convertOutput();
		if(er <= 0){
			return 0;
		}
		printOutput();
	}


	for(int tt = 0; tt < MaxTimeSteps1; ++tt){


		//next output time
		double timett1 = timeStart + dts * (tt + 1) * outputInterval;

//printf("integrate %.20g %.20g %.20g\n", timeStart + dts * tt * 10.0, timett1, dt);



		for(int i = 0; i < N; ++i){

			//dtmin is the minimum time step of an output intervall, only used for diagnostics
			dtmin_h[i] = outputInterval * dts;

			//integrate until the next output interval
			for(int ttt = 0; ttt < MaxTimeSteps2; ++ttt){
				if(dts < 0){
					if(dt_h[i] < -outputInterval){
						dt_h[i] = -outputInterval;
					}

					//refine last time step of interval to match output time
					if(time_h[i] + dt_h[i] < timett1){
//printf("refine %d %.20g | %.20g %.20g %.20g\n", i, dt_h[i], time_h[i] + dt_h[i], timett1, timett1 - time_h[i]);

						dtsave_h[i] = dt_h[i];
						dt_h[i] = timett1 - time_h[i];
						stop = 1;
					}
					else{
						dtmin_h[i] = dt_h[i] > dtmin_h[i] ? dt_h[i] : dtmin_h[i];
					}
				}
				else{
					if(dt_h[i] > outputInterval){
						dt_h[i] = outputInterval;
					}

					//refine last time step of interval to match output time
					if(time_h[i] + dt_h[i] > timett1){
//printf("refine %d %.20g | %.20g %.20g %.20g\n", i, dt_h[i], time_h[i] + dt_h[i], timett1, timett1 - time_h[i]);

						dtsave_h[i] = dt_h[i];
						dt_h[i] = timett1 - time_h[i];
						stop = 1;
					}
					else{
						dtmin_h[i] = dt_h[i] < dtmin_h[i] ? dt_h[i] : dtmin_h[i];
					}
				}

				//do a time step of length dt_h[i]
				if(strcmp(integratorName, "RKF45") == 0){
					er = RKF_individual_step(i);
				}
				if(strcmp(integratorName, "DP54") == 0){
					er = RKF_individual_step(i);
				}
				if(strcmp(integratorName, "RKF78") == 0){
					er = RKF_individual_step(i);
				}
				if(strcmp(integratorName, "BS") == 0){
					er = BS_individual_step(i);
				}

				if(er <= 0){
					return 0;
				}

				if(printdt == 1){
					fprintf(dtFile, "%-25.20g %lld %-25.20g %d\n", time_h[i] + time_reference, timeStep_h[i], dt_h[i], i);
				}



//printf("dt     %d %lld %g %g %g\n", i, timeStep_h[i], snew_h[i], dt_h[i], dtsave_h[i]);


				dt_h[i] *= snew_h[i];



				if(dts < 0 && time_h[i] <= timett1){
					//set time step equal to the last accepted full time step
					
					if(snew_h[i] >= 1.0 && stop == 1){ 
	
						dt_h[i] = dtsave_h[i];
//printf("reset  %d %.20g\n", i, dt_h[i]);
					}

					stop = 0;
					break;
				}
				if(dts > 0 && time_h[i] >= timett1){
					//set time step equal to the last accepted full time step


					if(snew_h[i] >= 1.0 && stop == 1){

						dt_h[i] = dtsave_h[i];
//printf("reset  %d %.20g\n", i, dt_h[i]);
					}
					stop = 0;
					break;
				}



				if(ttt >= MaxTimeSteps2 - 1){

					printf("Error, time step loop2 did not finish\n");
					return 0;
				}

			}//end of ttt loop

		}//end if i loop

		time = time_h[0];

		if(time + time_reference > time1 || time + time_reference < time0){
			printf("Reached the end of the Chebyshev data file\n");
			return 0;
		}


		double dtmin = dtmin_h[0];
		for(int i = 1; i < N; ++i){
			if(dts < 0){
				dtmin = dtmin > dtmin_h[i] ? dtmin : dtmin_h[i];
			}
			else{
				dtmin = dtmin < dtmin_h[i] ? dtmin : dtmin_h[i];
			}
			
		}
		dtminlevel_h[0] = dtmin;

		if(time_reference + time >= outStart){
			er = convertOutput();
			if(er <= 0){
				return 0;
			}
			printOutput();
			fflush(outputFile);
		}

		if(dts < 0 && time <= timeEnd){
			printf("Reached the end of the integration %.20g\n", time + time_reference);
			return 0;
		}
		if(dts > 0 && time >= timeEnd){
			printf("Reached the end of the integration %.20g\n", time + time_reference);
			return 0;
		}

		if(tt >= MaxTimeSteps1 - 1){

			printf("Error, time step loop1 did not finish\n");
			return 0;
		}
	}//end of tt loop

	return 1;
}



int asteroid::loop(){
	int er;

	//If needed, convert from heliocentric equatorial coordinates to barycentric equatorial coordinates
	if(ICheliocentric == 1){
		printf("Convert heliocentric to barycentric coordinates\n");
		er = HelioToBary(x_h, y_h, z_h, vx_h, vy_h, vz_h);
		if(er <= 0){
			return 0;
		}
		printf("Convert heliocentric to barycentric coordinates OK\n");
	} 


	//At this point, the initial conditions coordinates are cartesian barycentric equatorial
	//The integration is also done in cartesian barycentric equatorial coordinates

	if(outBinary == 0){
		outputFile = fopen(outputFilename, "w");
	}
	else{
		outputFile = fopen(outputFilename, "wb");
	}
	
	if(printdt == 1){
		dtFile = fopen(dtFilename, "w");
	}
	printf("Start integration %.20g\n", timeStart + time_reference);


	dt_h[0] = dt;

	if(nL > 1 && dt_h[0] * dts < dtlimit[0]){
		dt_h[0] = dtlimit[0] * dts;
	}

	for(int i = 1; i < nL; ++i){
		dt_h[i] = dtlimit[i - 1] * dts;
	}
	for(int i = 0; i < nL; ++i){
		dtsave_h[i] = dt;
		time_h[i] = time;
		stop_h[i] = 0;
		dtminlevel_h[i] = dt_h[i];
	}

	//reduce loop to numer of levels
	for(int i = 0; i < N; ++i){
		dtmin_h[i] = dt;
		timeStep_h[i] = 0ll;
	}

	if(time_reference + time >= outStart){
		er = convertOutput();
		if(er <= 0){
			return 0;
		}
		printOutput();
	}

	for(int tt = 0; tt < MaxTimeSteps1; ++tt){

		//dtmin is the minimum time step of an output intervall, only used for diagnostics
		for(int i = 0; i < nL; ++i){
			dtminlevel_h[i] = 1000.0 * dts;
		}

		for(int i = 0; i < N; ++i){
			dtmin_h[i] = outputInterval * dts;
		}



		//next output time
		double timett1 = timeStart + dts * (tt + 1) * outputInterval;

		double snew = 10.0;
//printf("integrate %.20g %.20g %.20g\n", timeStart + dts * tt * 10.0, timett1, dt_h[0]);


		//integrate until the next output interval
		for(int ttt = 0; ttt < MaxTimeSteps2; ++ttt){

			if(dts < 0){
				if(dt_h[0] < -outputInterval){
					dt_h[0] = -outputInterval;
				}


				//refine last time step of interval to match output time
				if(time + dt_h[0] < timett1){
//printf("refine %d %.20g | %.20g %.20g %.20g\n", 0, dt_h[0], time + dt_h[0], timett1, timett1 - time);
					dtsave_h[0] = dt_h[0];
					dt_h[0] = timett1 - time;
					stop_h[0] = 1;
				}
				else{
					dtminlevel_h[0] = dt_h[0] > dtminlevel_h[0] ? dt_h[0] : dtminlevel_h[0];
				}
			}
			else{
				if(dt_h[0] > outputInterval){
					dt_h[0] = outputInterval;
				}

				//refine last time step of interval to match output time
				if(time + dt_h[0] > timett1){
//printf("refine %d %.20g | %.20g %.20g %.20g\n", 0, dt_h[0], time + dt_h[0], timett1, timett1 - time);
					dtsave_h[0] = dt_h[0];
					dt_h[0] = timett1 - time;
					stop_h[0] = 1;
				}
				else{
					dtminlevel_h[0] = dt_h[0] < dtminlevel_h[0] ? dt_h[0] : dtminlevel_h[0];
				}

			}

			Nlevel_h[0] = N;
			for(int i = 1; i < nL; ++i){
				Nlevel_h[i] = 0;
			}


			//do a time step of length dt
			if(strcmp(integratorName, "LF") == 0){
				er = leapfrog_step();
			}
			if(strcmp(integratorName, "RK4") == 0){
				er = RK_step();
			}
			if(strcmp(integratorName, "RK7") == 0){
				er = RK_step();
			}
			if(strcmp(integratorName, "RKF45") == 0){
				er = RKF_step(0, dtlimit[0]);

				if(Nlevel_h[1] > 0){

					loop_recursive(1);

				}

				snew = snewlevel_h[0];
				timeStep = timeStep_h[0]; 
				time = time_h[0]; 
			}
			if(strcmp(integratorName, "DP54") == 0){
				er = RKF_step(0, dtlimit[0]);

				if(Nlevel_h[1] > 0){

					loop_recursive(1);

				}

				snew = snewlevel_h[0];
				timeStep = timeStep_h[0]; 
				time = time_h[0]; 
			}
			if(strcmp(integratorName, "RKF78") == 0){
				er = RKF_step(0, dtlimit[0]);

				if(Nlevel_h[1] > 0){

					loop_recursive(1);

				}

				snew = snewlevel_h[0];
				timeStep = timeStep_h[0]; 
				time = time_h[0]; 
			}
			if(strcmp(integratorName, "BS") == 0){
				er = BS_step(0, dtlimit[0]);

				if(Nlevel_h[1] > 0){

					loop_recursive(1);

				}

				snew = snewlevel_h[0];
				timeStep = timeStep_h[0]; 
				time = time_h[0]; 
			}
			if(strcmp(integratorName, "IMM") == 0){
				er = IMM_step();
				if(snew_h[0] == -1){
					printf("Error, implicit midpoint method did not converge\n");
					return 0;
				}
			}
			if(er <= 0 ){
				return 0;
			}		
	
	
			if(printdt == 1){
				fprintf(dtFile, "%-25.20g %lld %-25.20g %d\n", time + time_reference, timeStep, dt_h[0] / snew, 0);
			}
//printf("dt %d %lld %g %g %g\n", 0, timeStep, snew, dt_h[0], dtsave_h[0]);



			dt_h[0] *= snew;

			if(dts < 0 && time <= timett1){
				//set time step equal to the last accepted full time step

				if(snew >= 1.0 && stop_h[0] == 1){

					for(int l = 0; l < nL; ++l){
						dt_h[l] = dtsave_h[l];
//printf("reset %d  %.20g\n", l, dt_h[l]);

					}
				}

				stop_h[0] = 0;
				break;
			}
			if(dts > 0 && time >= timett1){
				//set time step equal to the last accepted full time step


				if(snew >= 1.0 && stop_h[0] == 1){

					for(int l = 0; l < nL; ++l){
						dt_h[l] = dtsave_h[l];
//printf("reset %d  %.20g\n", l, dt_h[l]);

					}

				}
				stop_h[0] = 0;
				break;
			}

			if(time + time_reference > time1 || time + time_reference < time0){
				printf("Reached the end of the Chebyshev data file\n");
				return 0;
			}


			if(ttt >= MaxTimeSteps2 - 1){

				printf("Error, time step loop2 did not finish\n");
				return 0;
			}

		}//end of ttt loop

		if(time_reference + time >= outStart){
			er = convertOutput();
			if(er <= 0){
				return 0;
			}
			printOutput();
			fflush(outputFile);
		}
	
		if(dts < 0 && time <= timeEnd){
			printf("Reached the end of the integration %.20g\n", time + time_reference);
			return 0;
		}
		if(dts > 0 && time >= timeEnd){
			printf("Reached the end of the integration %.20g\n", time + time_reference);
			return 0;
		}


		if(tt >= MaxTimeSteps1 - 1){

			printf("Error, time step loop1 did not finish\n");
			return 0;
		}
	}//end of tt loop
	return 1;
}

int asteroid::loop_recursive(int level){
	int er;

//printf("\nStart level %d %g %g\n", level, dt_h[level - 1], dt_h[level]);
	for(int t = 0; t < MaxTimeSteps2; ++t){


		if(dts < 0){
			if(dt_h[level] < dt_h[level - 1]){
				dt_h[level] = dt_h[level - 1];
			}


			//refine last time step of interval to match output time
			if(time_h[level] + dt_h[level] < time_h[level - 1]){
//printf("refine %d %.20g | %.20g %.20g %.20g\n", level, dt_h[level - 1], time_h[level] + dt_h[level], time_h[level - 1], time_h[level - 1] - time_h[level]);
				dtsave_h[level] = dt_h[level];
				dt_h[level] = time_h[level - 1] - time_h[level];
				stop_h[level] = 1;
			}
			else{
				dtminlevel_h[level] = dt_h[level] > dtminlevel_h[level] ? dt_h[level] : dtminlevel_h[level];
			}
		}
		else{
			if(dt_h[level] > dt_h[level - 1]){
				dt_h[level] = dt_h[level - 1];
			}

			//refine last time step of interval to match output time
			if(time_h[level] + dt_h[level] > time_h[level - 1]){
//printf("refine %d %.20g | %.20g %.20g %.20g\n", 1, dt_h[level - 1], time_h[level] + dt_h[level], time_h[level - 1], time_h[level - 1] - time_h[level]);
				dtsave_h[level] = dt_h[level];
				dt_h[level] = time_h[level - 1] - time_h[level];
				stop_h[level] = 1;
			}
			else{
				dtminlevel_h[level] = dt_h[level] < dtminlevel_h[level] ? dt_h[level] : dtminlevel_h[level];
			}

		}

		if(strcmp(integratorName, "RKF45") == 0){
			er = RKF_step(level, dtlimit[level]);
		}
		if(strcmp(integratorName, "DP54") == 0){
			er = RKF_step(level, dtlimit[level]);
		}
		if(strcmp(integratorName, "RKF78") == 0){
			er = RKF_step(level, dtlimit[level]);
		}
		if(strcmp(integratorName, "BS") == 0){
			er = BS_step(level, dtlimit[level]);
		}

		if(Nlevel_h[level + 1] > 0){

			loop_recursive(level + 1);

		}

		if(er <= 0){
			return 0;
		}


		if(printdt == 1){
			fprintf(dtFile, "%-25.20g %lld %-25.20g %d\n", time_h[level] + time_reference, timeStep, dt_h[level] / snewlevel_h[level], level);
		}
//printf("dt %d %lld %g %g %g\n", level, timeStep, snewlevel_h[level], dt_h[level], dtsave_h[level]);

		dt_h[level] *= snewlevel_h[level];



		if(dts < 0 && time_h[level] <= time_h[level - 1]){
			//set time step equal to the last accepted full time step

			if(snewlevel_h[level] >= 1.0 && stop_h[level] == 1){

				for(int l = level; l < nL; ++l){
					dt_h[l] = dtsave_h[l];
//printf("reset %d  %.20g\n", l, dt_h[l]);

				}
			}

			stop_h[level] = 0;
			break;
		}
		if(dts > 0 && time_h[level] >= time_h[level - 1]){
			//set time step equal to the last accepted full time step


			if(snewlevel_h[level] >= 1.0 && stop_h[level] == 1){

				for(int l = level; l < nL; ++l){
					dt_h[l] = dtsave_h[l];
//printf("reset %d  %.20g\n", l, dt_h[l]);

				}


			}
			stop_h[level] = 0;
			break;
		}
		if(t >= MaxTimeSteps2 - 1){

			printf("Error, time step loop2 in level %d did not finish\n", level);
			return 0;
		}


	} // end of t loop
	Nlevel_h[level] = 0;

	return 1;
}

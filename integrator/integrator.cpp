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
//printf("Drift %d %.20g %.20g %.20g %.20g\n", i, x_h[i], vx_h[i], dt, 0.5* dt * vx_h[i]);
		x_h[i] += 0.5 * dt * vx_h[i];
		y_h[i] += 0.5 * dt * vy_h[i];
		z_h[i] += 0.5 * dt * vz_h[i];
	}
	time += dt * 0.5;
//printf("ta %.20g %.20g\n", time, dt);

	// ----------------------------------------------------------------------------
	for(int i = 0; i < N; ++i){
		ax_h[i] = 0.0;
		ay_h[i] = 0.0;
		az_h[i] = 0.0;
	}
	// ----------------------------------------------------------------------------
	//Update the Chebyshev coefficients if necessary

	er = update_Chebyshev(time);
	if(er <= 0){
		return 0;
	}

	update_perturbers(time);
	// ----------------------------------------------------------------------------

	// ----------------------------------------------------------------------------
	//compute force
	for(int i = 0; i < N; ++i){

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
			NonGrav(xih, yih, zih, vxih, vyih, vzih, A1_h[i], A2_h[i], A3_h[i], r, ax_h[i], ay_h[i], az_h[i]);
		}
		if(useGR == 1){
			GR(xih, yih, zih, vxih, vyih, vzih, r, ax_h[i], ay_h[i], az_h[i], GM_h[10]);
		}
		if(useJ2 == 1){
			J2(xiE, yiE, ziE, ax_h[i], ay_h[i], az_h[i], GM_h[2]);
		}
		Gravity(x_h[i], y_h[i], z_h[i], xTable_h, yTable_h, zTable_h, ax_h[i], ay_h[i], az_h[i], i);

	}
	// ----------------------------------------------------------------------------

	//Kick
	for(int i = 0; i < N; ++i){
//printf("Kick %d %.20g %.20g %.20g %.20g\n", i, vx_h[i], ax_h[i], dt, dt * ax_h[i]);
		vx_h[i] += dt * ax_h[i];
		vy_h[i] += dt * ay_h[i];
		vz_h[i] += dt * az_h[i];
	}
	//Drift
	for(int i = 0; i < N; ++i){
//printf("Drift %d %.20g %.20g %.20g %.20g\n", i, x_h[i], vx_h[i], dt, 0.5* dt * vx_h[i]);
		x_h[i] += 0.5 * dt * vx_h[i];
		y_h[i] += 0.5 * dt * vy_h[i];
		z_h[i] += 0.5 * dt * vz_h[i];
	}
	time += dt * 0.5;
	++timeStep;
//printf("tb %.20g %.20g\n", time, dt); 
	return 1;
}

//Implicit Midpoint Method
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
		xt_h[i] = x_h[i];
		yt_h[i] = y_h[i];
		zt_h[i] = z_h[i];

		vxt_h[i] = vx_h[i];
		vyt_h[i] = vy_h[i];
		vzt_h[i] = vz_h[i];
	}


	for(int k = 0; k < 30; ++k){

		// ----------------------------------------------------------------------------
		//compute forces
		for(int i = 0; i < N; ++i){

			ax_h[i] = 0.0;
			ay_h[i] = 0.0;
			az_h[i] = 0.0;

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

			if(cometFlag > 0 && k == 0){
				if(i == 0){
					Tsave_h[timeStep] = time;
				}
				Rsave_h[i * Rbuffersize + timeStep] = r;
			}

			//Earth centric coordinates
			double xiE = xt_h[i] - xTable_h[2];
			double yiE = yt_h[i] - yTable_h[2];
			double ziE = zt_h[i] - zTable_h[2];

			if(useNonGrav == 1){
				NonGrav(xih, yih, zih, vxih, vyih, vzih, A1_h[i], A2_h[i], A3_h[i], r, ax_h[i], ay_h[i], az_h[i]);
			}
			if(useGR == 1){
				GR(xih, yih, zih, vxih, vyih, vzih, r, ax_h[i], ay_h[i], az_h[i], GM_h[10]);
			}
			if(useJ2 == 1){
				J2(xiE, yiE, ziE, ax_h[i], ay_h[i], az_h[i], GM_h[2]);
			}
			Gravity(xt_h[i], yt_h[i], zt_h[i], xTable_h, yTable_h, zTable_h, ax_h[i], ay_h[i], az_h[i], i);

		}
		// ----------------------------------------------------------------------------

		int stop = 1;
		for(int i = 0; i < N; ++i){
			//store old values to check for convergence
			xp_h[i] = xt_h[i];
			yp_h[i] = yt_h[i];
			zp_h[i] = zt_h[i];

			vxp_h[i] = vxt_h[i];
			vyp_h[i] = vyt_h[i];
			vzp_h[i] = vzt_h[i];

			xt_h[i] = x_h[i] + 0.5 * dt * vxt_h[i];
			yt_h[i] = y_h[i] + 0.5 * dt * vyt_h[i];
			zt_h[i] = z_h[i] + 0.5 * dt * vzt_h[i];

			vxt_h[i] = vx_h[i] + 0.5 * dt * ax_h[i];
			vyt_h[i] = vy_h[i] + 0.5 * dt * ay_h[i];
			vzt_h[i] = vz_h[i] + 0.5 * dt * az_h[i];

			if(fabs(xp_h[i] - xt_h[i]) >= 1.0e-20) stop = 0;
			if(fabs(yp_h[i] - yt_h[i]) >= 1.0e-20) stop = 0;
			if(fabs(zp_h[i] - zt_h[i]) >= 1.0e-20) stop = 0;

			if(fabs(vxp_h[i] - vxt_h[i]) >= 1.0e-20) stop = 0;
			if(fabs(vyp_h[i] - vyt_h[i]) >= 1.0e-20) stop = 0;
			if(fabs(vzp_h[i] - vzt_h[i]) >= 1.0e-20) stop = 0;
		}
		if(stop == 1){
//if(k > 1) printf("k %d\n", k);
			break;
		}
		if(k >= 29){
			snew_h[0] = -1;
		}



	}// end of k loop

	//update
	for(int i = 0; i < N; ++i){
		x_h[i] += dt * vxt_h[i];
		y_h[i] += dt * vyt_h[i];
		z_h[i] += dt * vzt_h[i];

		vx_h[i] += dt * ax_h[i];
		vy_h[i] += dt * ay_h[i];
		vz_h[i] += dt * az_h[i];
//printf("%d %g %g %g\n", i, x_h[i], y_h[i], z_h[i]);

	}

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
			ax_h[i] = 0.0;
			ay_h[i] = 0.0;
			az_h[i] = 0.0;

			xt_h[i] = x_h[i];
			yt_h[i] = y_h[i];
			zt_h[i] = z_h[i];

			vxt_h[i] = vx_h[i];
			vyt_h[i] = vy_h[i];
			vzt_h[i] = vz_h[i];

			for(int s = 0; s < S; ++s){
				double dtaa = dt * RKFa_h[S * RKFn + s];
				xt_h[i]  += dtaa * kx_h[i + s * N];
				yt_h[i]  += dtaa * ky_h[i + s * N];
				zt_h[i]  += dtaa * kz_h[i + s * N];
				vxt_h[i] += dtaa * kvx_h[i + s * N];
				vyt_h[i] += dtaa * kvy_h[i + s * N];
				vzt_h[i] += dtaa * kvz_h[i + s * N];
//printf("update 2 %d %d %g %g %g %g %g %g\n", S, i, xt_h[i], yt_h[i], zt_h[i], RKFa_h[S * RKFn + s], kx_h[s], dt);

			}

			kx_h[i + S * N] = vxt_h[i];
			ky_h[i + S * N] = vyt_h[i];
			kz_h[i + S * N] = vzt_h[i];

		}

		// ----------------------------------------------------------------------------
		//compute forces
		for(int i = 0; i < N; ++i){

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

			if(cometFlag > 0 && S == 0){
				if(i == 0){
					Tsave_h[timeStep] = time;
				}
				Rsave_h[i * Rbuffersize + timeStep] = r;
			}

			//Earth centric coordinates
			double xiE = xt_h[i] - xTable_h[2];
			double yiE = yt_h[i] - yTable_h[2];
			double ziE = zt_h[i] - zTable_h[2];

			if(useNonGrav == 1){
				NonGrav(xih, yih, zih, vxih, vyih, vzih, A1_h[i], A2_h[i], A3_h[i], r, ax_h[i], ay_h[i], az_h[i]);
			}
			if(useGR == 1){
				GR(xih, yih, zih, vxih, vyih, vzih, r, ax_h[i], ay_h[i], az_h[i], GM_h[10]);
			}
			if(useJ2 == 1){
				J2(xiE, yiE, ziE, ax_h[i], ay_h[i], az_h[i], GM_h[2]);
			}
			Gravity(xt_h[i], yt_h[i], zt_h[i], xTable_h, yTable_h, zTable_h, ax_h[i], ay_h[i], az_h[i], i);

		}
		// ----------------------------------------------------------------------------
		for(int i = 0; i < N; ++i){
			kvx_h[i + S * N] = ax_h[i];
			kvy_h[i + S * N] = ay_h[i];
			kvz_h[i + S * N] = az_h[i];
		}
	}

	//update
	for(int i = 0; i < N; ++i){

		dx_h[i] = 0.0;
		dy_h[i] = 0.0;
		dz_h[i] = 0.0;

		dvx_h[i] = 0.0;
		dvy_h[i] = 0.0;
		dvz_h[i] = 0.0;

		for(int S = 0; S < RKFn; ++S){
			double dtb = dt * RKFb_h[S];
			dx_h[i] += dtb * kx_h[i + S * N];
			dy_h[i] += dtb * ky_h[i + S * N];
			dz_h[i] += dtb * kz_h[i + S * N];

			dvx_h[i] += dtb * kvx_h[i + S * N];
			dvy_h[i] += dtb * kvy_h[i + S * N];
			dvz_h[i] += dtb * kvz_h[i + S * N];
		}
	}

	for(int i = 0; i < N; ++i){
		x_h[i] += dx_h[i];
		y_h[i] += dy_h[i];
		z_h[i] += dz_h[i];

		vx_h[i] += dvx_h[i];
		vy_h[i] += dvy_h[i];
		vz_h[i] += dvz_h[i];
//printf("%d %g %g %g\n", i, x_h[i], y_h[i], z_h[i]);

	}
	time += dt;
	++timeStep;
	return 1;
}

//Runge Kutta Fehlberg step with adaptive time step
inline int asteroid::RKF_step(){

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
			ax_h[i] = 0.0;
			ay_h[i] = 0.0;
			az_h[i] = 0.0;

			xt_h[i] = x_h[i];
			yt_h[i] = y_h[i];
			zt_h[i] = z_h[i];

			vxt_h[i] = vx_h[i];
			vyt_h[i] = vy_h[i];
			vzt_h[i] = vz_h[i];

			for(int s = 0; s < S; ++s){
				double dtaa = dt * RKFa_h[S * RKFn + s];
				int ii = i + s * N;
				xt_h[i]  += dtaa * kx_h[ii];
				yt_h[i]  += dtaa * ky_h[ii];
				zt_h[i]  += dtaa * kz_h[ii];
				vxt_h[i] += dtaa * kvx_h[ii];
				vyt_h[i] += dtaa * kvy_h[ii];
				vzt_h[i] += dtaa * kvz_h[ii];

			}
			int ii = i + S * N;
			kx_h[ii] = vxt_h[i];
			ky_h[ii] = vyt_h[i];
			kz_h[ii] = vzt_h[i];
//if(i < 2) printf("K %d %d %.20g %.20g %.20g\n", i, S, kx_h[i + S * N], ky_h[i + S * N], kz_h[i + S * N]);
		}

		// ----------------------------------------------------------------------------
		//compute forces
		for(int i = 0; i < N; ++i){

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

			if(cometFlag > 0 && S == 0){
				if(i == 0){
					Tsave_h[timeStep] = time;
				}
				Rsave_h[i * Rbuffersize + timeStep] = r;
			}

			//Earth centric coordinates
			double xiE = xt_h[i] - xTable_h[2];
			double yiE = yt_h[i] - yTable_h[2];
			double ziE = zt_h[i] - zTable_h[2];

			if(useNonGrav == 1){
				NonGrav(xih, yih, zih, vxih, vyih, vzih, A1_h[i], A2_h[i], A3_h[i], r, ax_h[i], ay_h[i], az_h[i]);
			}
			if(useGR == 1){
				GR(xih, yih, zih, vxih, vyih, vzih, r, ax_h[i], ay_h[i], az_h[i], GM_h[10]);
			}
			if(useJ2 == 1){
				J2(xiE, yiE, ziE, ax_h[i], ay_h[i], az_h[i], GM_h[2]);
			}
			Gravity(xt_h[i], yt_h[i], zt_h[i], xTable_h, yTable_h, zTable_h, ax_h[i], ay_h[i], az_h[i], i);
		}
		// ----------------------------------------------------------------------------
		for(int i = 0; i < N; ++i){
			int ii = i + S * N;
			kvx_h[ii] = ax_h[i];
			kvy_h[ii] = ay_h[i];
			kvz_h[ii] = az_h[i];
		}
	}

	//update
	for(int i = 0; i < N; ++i){

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

	}


	//compute integration error
	double snew = 10.0;

	for(int i = 0; i < N; ++i){

		//Hairer
/*
		scalex_h[i]  = atol + fmax(fabs(x_h[i]), fabs(x_h[i] + dx_h[i])) * rtol;
		scaley_h[i]  = atol + fmax(fabs(y_h[i]), fabs(y_h[i] + dy_h[i])) * rtol;
		scalez_h[i]  = atol + fmax(fabs(z_h[i]), fabs(z_h[i] + dz_h[i])) * rtol;

		scalevx_h[i] = atol + fmax(fabs(vx_h[i]), fabs(vx_h[i] + dvx_h[i])) * rtol;
		scalevy_h[i] = atol + fmax(fabs(vy_h[i]), fabs(vy_h[i] + dvy_h[i])) * rtol;
		scalevz_h[i] = atol + fmax(fabs(vz_h[i]), fabs(vz_h[i] + dvz_h[i])) * rtol;
*/

		scalex_h[i]  = atol + fabs(x_h[i]) * rtol;
		scaley_h[i]  = atol + fabs(y_h[i]) * rtol;
		scalez_h[i]  = atol + fabs(z_h[i]) * rtol;

		scalevx_h[i] = atol + fabs(vx_h[i]) * rtol;
		scalevy_h[i] = atol + fabs(vy_h[i]) * rtol;
		scalevz_h[i] = atol + fabs(vz_h[i]) * rtol;

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
		errork += errorkx * errorkx / (scalex_h[i] * scalex_h[i]);
		errork += errorky * errorky / (scaley_h[i] * scaley_h[i]);
		errork += errorkz * errorkz / (scalez_h[i] * scalez_h[i]);
		errork += errorkvx * errorkvx / (scalevx_h[i] * scalevx_h[i]);
		errork += errorkvy * errorkvy / (scalevy_h[i] * scalevy_h[i]);
		errork += errorkvz * errorkvz / (scalevz_h[i] * scalevz_h[i]);

		errork = sqrt(errork / 6.0);	//6 is the number of dimensions

		double s = pow( 1.0  / errork, RKF_ee);


		//s = (RKF_fac * s > RKF_facmin) ? RKF_fac * s : RKF_facmin;
		//s = (RKF_facmax < s) ? RKF_facmax : s;

		//time steps of power of two
		if(s > 2.0) s = 2.0;
		else if (s < 1.0) s = 0.5;
		else s = 1;
		

		snew = (snew < s) ? snew : s;
		snew_h[i] = snew;
	}


	snew = 10.0;
	//Find minimum time step
	for(int i = 0; i < N; ++i){
		snew = (snew < snew_h[i]) ? snew : snew_h[i];

	}
	snew_h[0] = snew;
	if(snew >= 1.0){
		//accept step
		for(int i = 0; i < N; ++i){
			x_h[i] += dx_h[i];
			y_h[i] += dy_h[i];
			z_h[i] += dz_h[i];

			vx_h[i] += dvx_h[i];
			vy_h[i] += dvy_h[i];
			vz_h[i] += dvz_h[i];
		}
		time += dt;
		++timeStep;
		if(stop != 1){
			//only increase time step when stop == 0
			dt *= snew;
		}

		//set maximum time step to 1
		//if(abs(dt) > 1.0){
		//	dt = dts * 1.0;
		//}

	}
	else{
		//redo step
		dt *= snew;
	}
//printf("dt %.20g %.20g %.20g\n", time, dt, snew);
	return 1;
}

// Bulirsh-Stoer step with adaptive time step
inline int asteroid::BS_step(){
	int er;

	for(int i = 0; i < N; ++i){

		scalex_h[i] = atol + fabs(x_h[i]) * rtol;
		scaley_h[i] = atol + fabs(y_h[i]) * rtol;
		scalez_h[i] = atol + fabs(z_h[i]) * rtol;

		scalevx_h[i] = atol + fabs(vx_h[i]) * rtol;
		scalevy_h[i] = atol + fabs(vy_h[i]) * rtol;
		scalevz_h[i] = atol + fabs(vz_h[i]) * rtol;
	}


	int f = 1;

	for(int n = 1; n <= 8; ++n){
		double dt2 = dt / (2.0 * n);
		double dt22 = dt2 * 2.0;

		// ----------------------------------------------------------------------------
		//Update the Chebyshev coefficients if necessary
		er = update_Chebyshev(time);
		if(er <= 0){
			return 0;
		}
		update_perturbers(time);
//printf("%d %d %g %g\n", n, 0, 0.0, dt);
		// ----------------------------------------------------------------------------

		for(int i = 0; i < N; ++i){
			ax_h[i] = 0.0;
			ay_h[i] = 0.0;
			az_h[i] = 0.0;

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
					Tsave_h[timeStep] = time;
				}
				Rsave_h[i * Rbuffersize + timeStep] = r;
			}

			//Earth centric coordinates
			double xiE = x_h[i] - xTable_h[2];
			double yiE = y_h[i] - yTable_h[2];
			double ziE = z_h[i] - zTable_h[2];

			if(useNonGrav == 1){
				NonGrav(xih, yih, zih, vxih, vyih, vzih, A1_h[i], A2_h[i], A3_h[i], r, ax_h[i], ay_h[i], az_h[i]);
			}
			if(useGR == 1){
				GR(xih, yih, zih, vxih, vyih, vzih, r, ax_h[i], ay_h[i], az_h[i], GM_h[10]);
			}
			if(useJ2 == 1){
				J2(xiE, yiE, ziE, ax_h[i], ay_h[i], az_h[i], GM_h[2]);
			}
			Gravity(x_h[i], y_h[i], z_h[i], xTable_h, yTable_h, zTable_h, ax_h[i], ay_h[i], az_h[i], i);
			// ----------------------------------------------------------------------------

			xp_h[i] = x_h[i] + dt2 * vx_h[i];
			yp_h[i] = y_h[i] + dt2 * vy_h[i];
			zp_h[i] = z_h[i] + dt2 * vz_h[i];

			vxp_h[i] = vx_h[i] + dt2 * ax_h[i];
			vyp_h[i] = vy_h[i] + dt2 * ay_h[i];
			vzp_h[i] = vz_h[i] + dt2 * az_h[i];
		}

		// ----------------------------------------------------------------------------
		//Update the Chebyshev coefficients if necessary
		er = update_Chebyshev(time + dt2);
		if(er <= 0){
			return 0;
		}
		update_perturbers(time + dt2);
//printf("%d %d %g\n", n, 1, dt2 / dt);
		// ----------------------------------------------------------------------------
		for(int i = 0; i < N; ++i){
			ax_h[i] = 0.0;
			ay_h[i] = 0.0;
			az_h[i] = 0.0;

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
				NonGrav(xih, yih, zih, vxih, vyih, vzih, A1_h[i], A2_h[i], A3_h[i], r, ax_h[i], ay_h[i], az_h[i]);
			}
			if(useGR == 1){
				GR(xih, yih, zih, vxih, vyih, vzih, r, ax_h[i], ay_h[i], az_h[i], GM_h[10]);
			}
			if(useJ2 == 1){
				J2(xiE, yiE, ziE, ax_h[i], ay_h[i], az_h[i], GM_h[2]);
			}
			Gravity(xp_h[i], yp_h[i], zp_h[i], xTable_h, yTable_h, zTable_h, ax_h[i], ay_h[i], az_h[i], i);
			// ----------------------------------------------------------------------------

			xt_h[i] = x_h[i] + dt22 * vxp_h[i];
			yt_h[i] = y_h[i] + dt22 * vyp_h[i];
			zt_h[i] = z_h[i] + dt22 * vzp_h[i];

			vxt_h[i] = vx_h[i] + dt22 * ax_h[i];
			vyt_h[i] = vy_h[i] + dt22 * ay_h[i];
			vzt_h[i] = vz_h[i] + dt22 * az_h[i];
		}
		for(int m = 2; m <= n; ++m){

			// ----------------------------------------------------------------------------
			//Update the Chebyshev coefficients if necessary
			er = update_Chebyshev(time + (m-1) * dt22);
			if(er <= 0){
				return 0;
			}
			update_perturbers(time + (m-1) * dt22);
//printf("%d %d %g\n", n, m, (m-1) * dt22 / dt);
			// ----------------------------------------------------------------------------
			for(int i = 0; i < N; ++i){
				ax_h[i] = 0.0;
				ay_h[i] = 0.0;
				az_h[i] = 0.0;

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
					NonGrav(xih, yih, zih, vxih, vyih, vzih, A1_h[i], A2_h[i], A3_h[i], r, ax_h[i], ay_h[i], az_h[i]);
				}
				if(useGR == 1){
					GR(xih, yih, zih, vxih, vyih, vzih, r, ax_h[i], ay_h[i], az_h[i], GM_h[10]);
				}
				if(useJ2 == 1){
					J2(xiE, yiE, ziE, ax_h[i], ay_h[i], az_h[i], GM_h[2]);
				}
				Gravity(xt_h[i], yt_h[i], zt_h[i], xTable_h, yTable_h, zTable_h, ax_h[i], ay_h[i], az_h[i], i);
				// ----------------------------------------------------------------------------

				xp_h[i] += dt22 * vxt_h[i];
				yp_h[i] += dt22 * vyt_h[i];
				zp_h[i] += dt22 * vzt_h[i];

				vxp_h[i] += dt22 * ax_h[i];
				vyp_h[i] += dt22 * ay_h[i];
				vzp_h[i] += dt22 * az_h[i];
			}
			// ----------------------------------------------------------------------------
			//Update the Chebyshev coefficients if necessary
			er = update_Chebyshev(time + (m-1) * dt22 + dt2);
			if(er <= 0){
				return 0;
			}
			update_perturbers(time + (m-1) * dt22 + dt2);
//printf("%d %d %g\n", n, m, ((m-1) * dt22 + dt2) / dt);
			// ----------------------------------------------------------------------------
			for(int i = 0; i < N; ++i){
				ax_h[i] = 0.0;
				ay_h[i] = 0.0;
				az_h[i] = 0.0;

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
					NonGrav(xih, yih, zih, vxih, vyih, vzih, A1_h[i], A2_h[i], A3_h[i], r, ax_h[i], ay_h[i], az_h[i]);
				}
				if(useGR == 1){
					GR(xih, yih, zih, vxih, vyih, vzih, r, ax_h[i], ay_h[i], az_h[i], GM_h[10]);
				}
				if(useJ2 == 1){
					J2(xiE, yiE, ziE, ax_h[i], ay_h[i], az_h[i], GM_h[2]);
				}
				Gravity(xp_h[i], yp_h[i], zp_h[i], xTable_h, yTable_h, zTable_h, ax_h[i], ay_h[i], az_h[i], i);
				// ----------------------------------------------------------------------------

				xt_h[i] += dt22 * vxp_h[i];
				yt_h[i] += dt22 * vyp_h[i];
				zt_h[i] += dt22 * vzp_h[i];

				vxt_h[i] += dt22 * ax_h[i];
				vyt_h[i] += dt22 * ay_h[i];
				vzt_h[i] += dt22 * az_h[i];
			}
		} // end if m loop

		// ----------------------------------------------------------------------------
		//Update the Chebyshev coefficients if necessary
		er = update_Chebyshev(time + dt);
		if(er <= 0){
			return 0;
		}
		update_perturbers(time + dt);
//printf("%d %d %g\n", n, n+1, 1.0);
		// ----------------------------------------------------------------------------
		for(int i = 0; i < N; ++i){
			ax_h[i] = 0.0;
			ay_h[i] = 0.0;
			az_h[i] = 0.0;

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
				NonGrav(xih, yih, zih, vxih, vyih, vzih, A1_h[i], A2_h[i], A3_h[i], r, ax_h[i], ay_h[i], az_h[i]);
			}
			if(useGR == 1){
				GR(xih, yih, zih, vxih, vyih, vzih, r, ax_h[i], ay_h[i], az_h[i], GM_h[10]);
			}
			if(useJ2 == 1){
				J2(xiE, yiE, ziE, ax_h[i], ay_h[i], az_h[i], GM_h[2]);
			}
			Gravity(xt_h[i], yt_h[i], zt_h[i], xTable_h, yTable_h, zTable_h, ax_h[i], ay_h[i], az_h[i], i);
			// ----------------------------------------------------------------------------

			xp_h[i] += dt2 * vxt_h[i];
			yp_h[i] += dt2 * vyt_h[i];
			zp_h[i] += dt2 * vzt_h[i];

			vxp_h[i] += dt2 * ax_h[i];
			vyp_h[i] += dt2 * ay_h[i];
			vzp_h[i] += dt2 * az_h[i];
		}

		for(int i = 0; i < N; ++i){
			dx_h[i * 8 + (n-1)] = 0.5 * (xt_h[i] + xp_h[i]);
			dy_h[i * 8 + (n-1)] = 0.5 * (yt_h[i] + yp_h[i]);
			dz_h[i * 8 + (n-1)] = 0.5 * (zt_h[i] + zp_h[i]);

			dvx_h[i * 8 + (n-1)] = 0.5 * (vxt_h[i] + vxp_h[i]);
			dvy_h[i * 8 + (n-1)] = 0.5 * (vyt_h[i] + vyp_h[i]);
			dvz_h[i * 8 + (n-1)] = 0.5 * (vzt_h[i] + vzp_h[i]);
		}


		//Extrapolation step
		double errormax = 0.0;
		for(int i = 0; i < N; ++i){
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

			errormax = error > errormax ? error: errormax;

		}
//printf("error %d %g\n", n, errormax);
		if(errormax < 1.0){
//printf("accept %lld %.20g %g\n", timeStep, time, dt);

			//accept step
			snew_h[0] = 1.0;
			for(int i = 0; i < N; ++i){
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
			}

			time += dt;
			++timeStep;
			if(stop != 1){
				//only increase time step when stop == 0
				if(n < 6){
					dt *= 2.0;
					snew_h[0] = 2.0;
//printf("increase time step %g\n", dt);
				}
			}

			for(int i = 0; i < N; ++i){
				x_h[i] = xt_h[i];
				y_h[i] = yt_h[i];
				z_h[i] = zt_h[i];

				vx_h[i] = vxt_h[i];
				vy_h[i] = vyt_h[i];
				vz_h[i] = vzt_h[i];
			}

			f = 0;
			break; //break n loop
		}
	} //end of n loop
	if(f == 1){
		dt *= 0.5;
		snew_h[0] = 0.5;
//printf("repeat time step %g\n", dt);
	}
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

	if(time_reference + time >= outStart){
		er = convertOutput();
		if(er <= 0){
			return 0;
		}
		printOutput(dt);
	}

	//for(int tt = 0; tt < 2; ++tt){
	for(int tt = 0; tt < MaxTimeSteps1; ++tt){

		//dtmin is the minimum time step of an output intervall, only used for diagnostics
		double dtmin = dt;

		//next output time
		double timett1 = timeStart + dts * (tt + 1) * outputInterval;

		double snew = 10.0;
//printf("integrate %.20g %.20g\n", timeStart + dts * tt * 10.0, timett1);


		//integrate until the next output interval
		for(int ttt = 0; ttt < MaxTimeSteps2; ++ttt){

			//refine last time step of interval to match output time
			if(dts < 0){
				if(time + dt < timett1){
//printf("refine %.20g %.20g %.20g\n", time + dt, timett1, timett1 - time);
//fprintf(dtFile, "refine %.20g %.20g %.20g\n", time + dt, timett1, timett1 - time);
					dt1 = dt;
					dt = (timett1 - time);
					stop = 1;
				}
			}
			else{
				if(time + dt > timett1){
//printf("refine %.20g %.20g %.20g\n", time + dt, timett1, timett1 - time);
//fprintf(dtFile, "refine %.20g %.20g %.20g\n", time + dt, timett1, timett1 - time);
					dt1 = dt;
					dt = (timett1 - time);
					stop = 1;
				}

			}

			//do a time step of length dt
			if(strcmp(integratorName, "LF") == 0){
				er = leapfrog_step();
//				time = timeStart + timeStep * dt;
			}
			if(strcmp(integratorName, "RK4") == 0){
				er = RK_step();
			}
			if(strcmp(integratorName, "RK7") == 0){
				er = RK_step();
			}
			if(strcmp(integratorName, "RKF45") == 0){
				er = RKF_step();
				snew = snew_h[0];
			}
			if(strcmp(integratorName, "DP54") == 0){
				er = RKF_step();
				snew = snew_h[0];
			}
			if(strcmp(integratorName, "RKF78") == 0){
				er = RKF_step();
				snew = snew_h[0];
			}
			if(strcmp(integratorName, "BS") == 0){
				er = BS_step();
				snew = snew_h[0];
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
	
	
			if(stop == 0){
				dtmin = (abs(dt) < abs(dtmin)) ? dt : dtmin;
				if(printdt == 1){
					fprintf(dtFile, "%.20g %lld %.20g\n", time + time_reference, timeStep, dt);
				}
			}

			if(time + time_reference > time1 || time + time_reference < time0){
				printf("Reached the end of the Chebyshev data file\n");
				return 0;
			}

			if(dts < 0 && time < timeEnd){
				printf("Reached the end of the integration\n");
				return 0;
			}
			if(dts > 0 && time > timeEnd){
				printf("Reached the end of the integration\n");
				return 0;
			}
//convertOutput();
//printOutput(dtmin);

			if(stop == 1){
				stop = 0;
				if(snew >= 1.0){
					//set time step equal to the last accepted full time step
					dt = dt1;
//printf("reset %.20g\n", dt);
//fprintf(dtFile, "reset %.20g\n", dt);
					break;
				}
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
			printOutput(dtmin);
			fflush(outputFile);
		}

		if(tt >= MaxTimeSteps1 - 1){

			printf("Error, time step loop1 did not finish\n");
			return 0;
		}
	}//end of tt loop
	fclose(outputFile);
	if(printdt == 1){
		fclose(dtFile);
	}
	return 1;
}

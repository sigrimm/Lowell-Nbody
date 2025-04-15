#include <stdio.h>

#include "asteroid.h"
#include "Chebyshev.h"
#include "force.h"

void asteroid::HelioToBary(double *xx_h, double *yy_h, double *zz_h, double *vxx_h, double *vyy_h, double *vzz_h){
	//Update the Chebyshev coefficients if necessary
	update_Chebyshev(timeStart);
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


		//printf("xyz bary %.40g %.40g %.40g %.40g %.40g %.40g %.40g %.40g %.40g\n", x_h[i], y_h[i], z_h[i], vx_h[i], vy_h[i], vz_h[i], A1_h[i], A2_h[i], A3_h[i]);
	}
}

void asteroid::BaryToHelio(double *xx_h, double *yy_h, double *zz_h, double *vxx_h, double *vyy_h, double *vzz_h){
	//Update the Chebyshev coefficients if necessary
	update_Chebyshev(time);
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


		//printf("xyz bary %.40g %.40g %.40g %.40g %.40g %.40g %.40g %.40g %.40g\n", x_h[i], y_h[i], z_h[i], vx_h[i], vy_h[i], vz_h[i], A1_h[i], A2_h[i], A3_h[i]);
	}
}

void asteroid::convertOutput(){

	for(int i = 0; i < N; ++i){
		xout_h[i] = x_h[i];
		yout_h[i] = y_h[i];
		zout_h[i] = z_h[i];

		vxout_h[i] = vx_h[i];
		vyout_h[i] = vy_h[i];
		vzout_h[i] = vz_h[i];
	}

	if(Outheliocentric == 1){
		BaryToHelio(xout_h, yout_h, zout_h, vxout_h, vyout_h, vzout_h);
	} 


	//If needed, convert from equatorial coordinates to ecliptic coordinates
	if(Outecliptic == 1){
		EquatorialtoEcliptic(xout_h, yout_h, zout_h, vxout_h, vyout_h, vzout_h);        
	}

	if(Outorbital == 1){
		for(int i = 0; i < N; ++i){
			double a, e, inc, Omega, w, Theta, M, E;
			CartToKep(xout_h, yout_h, zout_h, vxout_h, vyout_h, vzout_h, i, a, e, inc, Omega, w, Theta, M, E);


			inc = inc * 180.0 / M_PI;       //convert rad to deg
			Omega = Omega * 180.0 / M_PI;   //convert rad to deg
			w = w * 180.0 / M_PI;           //convert rad to deg
			M = M * 180.0 / M_PI;           //convert rad to deg

			xout_h[i] = a;
			yout_h[i] = e;
			zout_h[i] = inc;
			vxout_h[i] = Omega;
			vyout_h[i] = w;
			vzout_h[i] = M;


		}
	}

}



//Leapfrog step with fixed time step
inline void asteroid::leapfrog_step(){
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
	update_Chebyshev(time);
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
}


//Runge Kutta step with fixed time step
inline void asteroid::RK_step(){

	for(int S = 0; S < RKFn; ++S){

		// ----------------------------------------------------------------------------
		//Update the Chebyshev coefficients if necessary
		update_Chebyshev(time + RKFc_h[S] * dt);
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
}

//Runge Kutta Fehlberg step with adaptive time step
inline void asteroid::RKF_step(){

	for(int S = 0; S < RKFn; ++S){

		// ----------------------------------------------------------------------------
		//Update the Chebyshev coefficients if necessary
		update_Chebyshev(time + RKFc_h[S] * dt);
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
//if(i < 2) printf("update 2 %d %d %g %g %g %g %g %g\n", S, i, xt_h[i], yt_h[i], zt_h[i], RKFa_h[S * RKFn + s], kx_h[s], dt);

			}

			kx_h[i + S * N] = vxt_h[i];
			ky_h[i + S * N] = vyt_h[i];
			kz_h[i + S * N] = vzt_h[i];
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
//if(i < 2) printf("dx %d %d %g %g %g %g %g %g\n", S, i, dx_h[i], dy_h[i], dz_h[i], RKFb_h[S], kx_h[i + S * N], dt);
		}

	}


	//compute integration error
	double snew = 10.0;

	for(int i = 0; i < N; ++i){
		double ym = 0.0;
		ym = (fabs(x_h[i]) > ym) ? fabs(x_h[i]) : ym;
		ym = (fabs(y_h[i]) > ym) ? fabs(y_h[i]) : ym;
		ym = (fabs(z_h[i]) > ym) ? fabs(z_h[i]) : ym;

		ym = (fabs(vx_h[i]) > ym) ? fabs(vx_h[i]) : ym;
		ym = (fabs(vy_h[i]) > ym) ? fabs(vy_h[i]) : ym;
		ym = (fabs(vz_h[i]) > ym) ? fabs(vz_h[i]) : ym;

		double isc = 1.0 / (RKF_atol + ym * RKF_rtol);
		isc *= isc;

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
//if(i < 2) printf("error %d %d %.20g %.20g %.20g %.20g %.20g %.20g %.20g\n", i, S, f, errorkx, errorky, errorkz, errorkvx, errorkvy, errorkvz);
		}

		double errork = 0.0;
		errork += errorkx * errorkx * isc;
		errork += errorky * errorky * isc;
		errork += errorkz * errorkz * isc;
		errork += errorkvx * errorkvx * isc;
		errork += errorkvy * errorkvy * isc;
		errork += errorkvz * errorkvz * isc;

		errork = sqrt(errork / 6.0);	//6 is the number of dimensions

		double s = pow( 1.0  / errork, RKF_ee);

		s = (RKF_fac * s > RKF_facmin) ? RKF_fac * s : RKF_facmin;
		s = (RKF_facmax < s) ? RKF_facmax : s;

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
}


	
int asteroid::loop(){

	//If needed, convert from heliocentric equatorial coordinates to barycentric equatorial coordinates
	if(ICheliocentric == 1){
		HelioToBary(x_h, y_h, z_h, vx_h, vy_h, vz_h);
	} 


	//At this point, the initial conditions coordinates are cartesian barycentric equatorial
	//The integration is also done in cartesian barycentric equatorial coordinates

	if(outBinary == 0){
		outputFile = fopen(outputFilename, "w");
	}
	else{
		outputFile = fopen(outputFilename, "wb");
	}
	printf("Start integration %.20g\n", timeStart + time_reference);

	if(time_reference + time >= outStart){
		convertOutput();
		printOutput(dt);
	}

	//for(int tt = 0; tt < 2; ++tt){
	for(int tt = 0; tt < 1000000; ++tt){
		double dtmin = dt;

		double timett1 = timeStart + dts * (tt + 1) * outputInterval;


		double snew = 10.0;
//printf("integrate %.20g %.20g\n", timeStart + dts * tt * 10.0, timett1);


		//integrate until the next output interval
		for(int ttt = 0; ttt < 1000000; ++ttt){

			//refine last time step of interval to match output time
			if(dts < 0){
				if(time + dt < timett1){
					dt1 = dt;
					dt = (timett1 - time);
					stop = 1;
//printf("refine %.20g\n", timett1 - time);
				}
			}
			else{
				if(time + dt > timett1){
					dt1 = dt;
					dt = (timett1 - time);
					stop = 1;
//printf("refine %.20g\n", timett1 - time);
				}

			}
			if(RKFn == 1){
				leapfrog_step();
			}
			if(RKFn == 4){
				RK_step();
			}
			if(RKFn == 6){
				RKF_step();
			}
			if(RKFn == 7){
				RKF_step();
			}
			if(RKFn == 13){
				RKF_step();
				snew = snew_h[0];
			}
			

			if(stop == 0){
				dtmin = (abs(dt) < abs(dtmin)) ? dt : dtmin;
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

			if(stop == 1){
				stop = 0;
				if(snew >= 1.0){
					//set time step equal to the last accepted full time step
					dt = dt1;
					break;
				}
			}

			if(ttt >= 1000000 - 1){

				printf("Error, time step loop did not finish\n");
				return 0;
			}

		}//end of ttt loop

		if(time_reference + time >= outStart){
			convertOutput();
			printOutput(dtmin);
		}

	}//end of tt loop
	fclose(outputFile);
	return 1;
}

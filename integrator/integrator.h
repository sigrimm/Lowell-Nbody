#include "asteroid.h"

//Leapfrog step with fixed time step
inline void asteroid::leapfrog_step(){
	//Drift
	for(int i = Nperturbers; i < N; ++i){
//printf("Drift %d %.20g %.20g %.20g %.20g\n", i, x_h[i], vx_h[i], dt, 0.5* dt * vx_h[i]);
		x_h[i] += 0.5* dt * vx_h[i];
		y_h[i] += 0.5* dt * vy_h[i];
		z_h[i] += 0.5* dt * vz_h[i];
	}
	time += dt / 2.0;
	//printf("ta %.20g\n", time);   

	// ----------------------------------------------------------------------------
	for(int i = Nperturbers; i < N; ++i){
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
	//compute forces
	NonGrav(x_h, y_h, z_h, vx_h, vy_h, vz_h);
	GR(x_h, y_h, z_h, vx_h, vy_h, vz_h);
	J2(x_h, y_h, z_h);
	Gravity(x_h, y_h, z_h);
	// ----------------------------------------------------------------------------

	//Kick
	for(int i = Nperturbers; i < N; ++i){
//printf("Kick %d %.20g %.20g %.20g %.20g\n", i, vx_h[i], ax_h[i], dt, dt * ax_h[i]);
		vx_h[i] += dt * ax_h[i];
		vy_h[i] += dt * ay_h[i];
		vz_h[i] += dt * az_h[i];
	}
	//Drift
	for(int i = Nperturbers; i < N; ++i){
//printf("Drift %d %.20g %.20g %.20g %.20g\n", i, x_h[i], vx_h[i], dt, 0.5* dt * vx_h[i]);
		x_h[i] += 0.5* dt * vx_h[i];
		y_h[i] += 0.5* dt * vy_h[i];
		z_h[i] += 0.5* dt * vz_h[i];
	}
	time += dt / 2.0;
	//printf("tb %.20g\n", time); 
}


//Runge Kutta step with fixed time step
inline void asteroid::RK_step(){

	for(int S = 0; S < RKFn; ++S){

		// ----------------------------------------------------------------------------
		//Update the Chebyshev coefficients if necessary
		update_Chebyshev(time + c_h[S] * dt);
		update_perturbers(time + c_h[S] * dt);

		for(int i = 0; i < Nperturbers; ++i){
			xt_h[i] = x_h[i];
			yt_h[i] = y_h[i];
			zt_h[i] = z_h[i];

			vxt_h[i] = vx_h[i];
			vyt_h[i] = vy_h[i];
			vzt_h[i] = vz_h[i];
		}
		// ----------------------------------------------------------------------------

		for(int i = Nperturbers; i < N; ++i){
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
				double dtaa = dt * a_h[S * RKFn + s];
				xt_h[i]  += dtaa * kx_h[i + s * N];
				yt_h[i]  += dtaa * ky_h[i + s * N];
				zt_h[i]  += dtaa * kz_h[i + s * N];
				vxt_h[i] += dtaa * kvx_h[i + s * N];
				vyt_h[i] += dtaa * kvy_h[i + s * N];
				vzt_h[i] += dtaa * kvz_h[i + s * N];
//printf("update 2 %d %d %g %g %g %g %g %g\n", S, i, xt_h[i], yt_h[i], zt_h[i], a_h[S * RKFn + s], kx_h[s], dt);

			}

			kx_h[i + S * N] = vxt_h[i];
			ky_h[i + S * N] = vyt_h[i];
			kz_h[i + S * N] = vzt_h[i];

		}

		// ----------------------------------------------------------------------------
		//compute forces
		NonGrav(xt_h, yt_h, zt_h, vxt_h, vyt_h, vzt_h);
		GR(xt_h, yt_h, zt_h, vxt_h, vyt_h, vzt_h);
		J2(xt_h, yt_h, zt_h);
		Gravity(xt_h, yt_h, zt_h);
		// ----------------------------------------------------------------------------
		for(int i = Nperturbers; i < N; ++i){
			kvx_h[i + S * N] = ax_h[i];
			kvy_h[i + S * N] = ay_h[i];
			kvz_h[i + S * N] = az_h[i];
		}
	}

	//update
	for(int i = Nperturbers; i < N; ++i){

		dx_h[i] = 0.0;
		dy_h[i] = 0.0;
		dz_h[i] = 0.0;

		dvx_h[i] = 0.0;
		dvy_h[i] = 0.0;
		dvz_h[i] = 0.0;

		for(int S = 0; S < RKFn; ++S){
			double dtb = dt * b_h[S];
			dx_h[i] += dtb * kx_h[i + S * N];
			dy_h[i] += dtb * ky_h[i + S * N];
			dz_h[i] += dtb * kz_h[i + S * N];

			dvx_h[i] += dtb * kvx_h[i + S * N];
			dvy_h[i] += dtb * kvy_h[i + S * N];
			dvz_h[i] += dtb * kvz_h[i + S * N];
		}
	}

	for(int i = Nperturbers; i < N; ++i){
		x_h[i] += dx_h[i];
		y_h[i] += dy_h[i];
		z_h[i] += dz_h[i];

		vx_h[i] += dvx_h[i];
		vy_h[i] += dvy_h[i];
		vz_h[i] += dvz_h[i];

	}
	time += dt;
}

//Runge Kutta Fehlberg step with adaptive time step
inline void asteroid::RKF_step(){

	for(int S = 0; S < RKFn; ++S){

		// ----------------------------------------------------------------------------
		//Update the Chebyshev coefficients if necessary
		update_Chebyshev(time + c_h[S] * dt);
		update_perturbers(time + c_h[S] * dt);

		for(int i = 0; i < Nperturbers; ++i){
			xt_h[i] = x_h[i];
			yt_h[i] = y_h[i];
			zt_h[i] = z_h[i];

			vxt_h[i] = vx_h[i];
			vyt_h[i] = vy_h[i];
			vzt_h[i] = vz_h[i];
		}
		// ----------------------------------------------------------------------------

		for(int i = Nperturbers; i < N; ++i){
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
				double dtaa = dt * a_h[S * RKFn + s];
				xt_h[i]  += dtaa * kx_h[i + s * N];
				yt_h[i]  += dtaa * ky_h[i + s * N];
				zt_h[i]  += dtaa * kz_h[i + s * N];
				vxt_h[i] += dtaa * kvx_h[i + s * N];
				vyt_h[i] += dtaa * kvy_h[i + s * N];
				vzt_h[i] += dtaa * kvz_h[i + s * N];
//printf("update 2 %d %d %g %g %g %g %g %g\n", S, i, xt_h[i], yt_h[i], zt_h[i], a_h[S * RKFn + s], kx_h[s], dt);

			}

			kx_h[i + S * N] = vxt_h[i];
			ky_h[i + S * N] = vyt_h[i];
			kz_h[i + S * N] = vzt_h[i];

		}

		// ----------------------------------------------------------------------------
		//compute forces
		NonGrav(xt_h, yt_h, zt_h, vxt_h, vyt_h, vzt_h);
		GR(xt_h, yt_h, zt_h, vxt_h, vyt_h, vzt_h);
		J2(xt_h, yt_h, zt_h);
		Gravity(xt_h, yt_h, zt_h);
		// ----------------------------------------------------------------------------
		for(int i = Nperturbers; i < N; ++i){
			kvx_h[i + S * N] = ax_h[i];
			kvy_h[i + S * N] = ay_h[i];
			kvz_h[i + S * N] = az_h[i];
		}
	}

	//update
	for(int i = Nperturbers; i < N; ++i){

		dx_h[i] = 0.0;
		dy_h[i] = 0.0;
		dz_h[i] = 0.0;

		dvx_h[i] = 0.0;
		dvy_h[i] = 0.0;
		dvz_h[i] = 0.0;

		for(int S = 0; S < RKFn; ++S){
			double dtb = dt * b_h[S];
			dx_h[i] += dtb * kx_h[i + S * N];
			dy_h[i] += dtb * ky_h[i + S * N];
			dz_h[i] += dtb * kz_h[i + S * N];

			dvx_h[i] += dtb * kvx_h[i + S * N];
			dvy_h[i] += dtb * kvy_h[i + S * N];
			dvz_h[i] += dtb * kvz_h[i + S * N];
		}
	}


	//compute integration error
	double snew = 10.0;

	for(int i = Nperturbers; i < N; ++i){
		double ym = 0.0;
		ym = (fabs(x_h[i]) > ym) ? fabs(x_h[i]) : ym;
		ym = (fabs(y_h[i]) > ym) ? fabs(y_h[i]) : ym;
		ym = (fabs(z_h[i]) > ym) ? fabs(z_h[i]) : ym;

		ym = (fabs(vx_h[i]) > ym) ? fabs(vx_h[i]) : ym;
		ym = (fabs(vy_h[i]) > ym) ? fabs(vy_h[i]) : ym;
		ym = (fabs(vz_h[i]) > ym) ? fabs(vz_h[i]) : ym;

		double isc = 1.0 / (RKF_atol + ym * RKF_rtol);

		//error estimation
		double errorkx = 0.0;
		double errorky = 0.0;
		double errorkz = 0.0;
		
		double errorkvx = 0.0;
		double errorkvy = 0.0;
		double errorkvz = 0.0;

		for(int S = 0; S < RKFn; ++S){
			double f = (b_h[S] - bb_h[S]) * dt;
			errorkx += f * kx_h[i + S * N];
			errorky += f * ky_h[i + S * N];
			errorkz += f * kz_h[i + S * N];

			errorkvx += f * kvx_h[i + S * N];
			errorkvy += f * kvy_h[i + S * N];
			errorkvz += f * kvz_h[i + S * N];
//printf("error %d %d %g %g\n", i, S, errorkx, kx[S]);
		}

		double errork = 0.0;
		errork += errorkx * errorkx * isc * isc;
		errork += errorky * errorky * isc * isc;
		errork += errorkz * errorkz * isc * isc;
		errork += errorkvx * errorkvx * isc * isc;
		errork += errorkvy * errorkvy * isc * isc;
		errork += errorkvz * errorkvz * isc * isc;

		errork = sqrt(errork / 6.0);    //6 is the number of dimensions

		double s = pow( 1.0  / errork, RKF_ee);
		s = fmax(RKF_facmin, RKF_fac * s);
		s = fmin(RKF_facmax, s);

		snew = (snew < s) ? snew : s;
	}

	if(stop == 1){
		//accept step
		for(int i = Nperturbers; i < N; ++i){
			x_h[i] += dx_h[i];
			y_h[i] += dy_h[i];
			z_h[i] += dz_h[i];

			vx_h[i] += dvx_h[i];
			vy_h[i] += dvy_h[i];
			vz_h[i] += dvz_h[i];

		}
		time += dt;
	}
	else if(snew >= 1.0){
		//accept step
		for(int i = Nperturbers; i < N; ++i){
			x_h[i] += dx_h[i];
			y_h[i] += dy_h[i];
			z_h[i] += dz_h[i];

			vx_h[i] += dvx_h[i];
			vy_h[i] += dvy_h[i];
			vz_h[i] += dvz_h[i];

		}
		time += dt;
		dt *= snew;

		//set maximum time step to 1
		if(abs(dt) > 1.0){
			dt = dts * 1.0;
		}

	}
	else{
		//redo step
		dt *= snew;
	}

//printf("dt %.20g %.20g %.20g\n", time, dt, snew);
}


	
inline int asteroid::loop(){


	outputFile = fopen("Out.dat", "w");
#if USEGPU == 1
		copyOutput();
#endif
	for(int p = Nperturbers; p < N; ++p){
		printf("start integration %.20g %.20g\n", time_reference + time, dt);
		fprintf(outputFile, "%.20g %d %.20g %.20g %.20g %.20g %.20g %.20g %.20g\n", time_reference + time, p, x_h[p], y_h[p], z_h[p], vx_h[p], vy_h[p], vz_h[p], dt);
	}
	//for(int tt = 0; tt < 2; ++tt){
	for(int tt = 0; tt < 1000000; ++tt){
		double dtmin = dt;

		double timett1 = timeStart + dts * (tt + 1) * outputInterval;

//printf("integrate %.20g %.20g\n", timeStart + dts * tt * 10.0, timett1);


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
//printf("refine %.20g\n", timett1 - A.time);
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
			}

			dtmin = (abs(dt) < abs(dtmin)) ? dt : dtmin;

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
				dt = dt1;
				break;
			}

			if(ttt >= 1000000 - 1){

				printf("Error time step loop did not finish\n");
				return 0;
			}

		}//end of ttt loop
#if USEGPU == 1
		copyOutput();
#endif
		for(int p = Nperturbers; p < N; ++p){
			printf("Reached time %.20g %.20g\n", time_reference + time, dt);
			fprintf(outputFile, "%.20g %d %.20g %.20g %.20g %.20g %.20g %.20g %.20g\n", time_reference + time, p, x_h[p], y_h[p], z_h[p], vx_h[p], vy_h[p], vz_h[p], dtmin);
		}

	}//end of tt loop
	fclose(outputFile);
	return 1;
}

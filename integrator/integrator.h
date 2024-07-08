#include "asteroid.h"

//Leapfrog step with fixed time step
inline void asteroid::leapfrog_step(){
	//Drift
	for(int p = Nperturbers; p < N; ++p){
//printf("Drift %d %.20g %.20g %.20g %.20g\n", p, x[p], vx[p], dt, 0.5* dt * vx[p]);
		x[p] += 0.5* dt * vx[p];
		y[p] += 0.5* dt * vy[p];
		z[p] += 0.5* dt * vz[p];
	}
	time += dt / 2.0;
	//printf("ta %.20g\n", time);   

	// ----------------------------------------------------------------------------
	for(int i = Nperturbers; i < N; ++i){
		ax[i] = 0.0;
		ay[i] = 0.0;
		az[i] = 0.0;
	}
	// ----------------------------------------------------------------------------
	//Update the Chebyshev coefficients if necessary
	update_Chebyshev(time);
	update_perturbers(time);
	// ----------------------------------------------------------------------------

	// ----------------------------------------------------------------------------
	//compute forces
	NonGrav(x, y, z, vx, vy, vz);
	GR(x, y, z, vx, vy, vz);
	J2(x, y, z);
	Gravity(x, y, z);
	// ----------------------------------------------------------------------------

	//Kick
	for(int i = Nperturbers; i < N; ++i){
//printf("Kick %d %.20g %.20g %.20g %.20g\n", i, vx[i], ax[i], dt, dt * ax[i]);
		vx[i] += dt * ax[i];
		vy[i] += dt * ay[i];
		vz[i] += dt * az[i];
	}
	//Drift
	for(int p = Nperturbers; p < N; ++p){
//printf("Drift %d %.20g %.20g %.20g %.20g\n", p, x[p], vx[p], dt, 0.5* dt * vx[p]);
		x[p] += 0.5* dt * vx[p];
		y[p] += 0.5* dt * vy[p];
		z[p] += 0.5* dt * vz[p];
	//printf("%.20g %d %.20g %.20g %.20g\n", time + dt / 2.0, p, x[p], y[p], z[p]);
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
			xt[i] = x[i];
			yt[i] = y[i];
			zt[i] = z[i];

			vxt[i] = vx[i];
			vyt[i] = vy[i];
			vzt[i] = vz[i];
		}
		// ----------------------------------------------------------------------------

		for(int i = Nperturbers; i < N; ++i){
			ax[i] = 0.0;
			ay[i] = 0.0;
			az[i] = 0.0;

			xt[i] = x[i];
			yt[i] = y[i];
			zt[i] = z[i];

			vxt[i] = vx[i];
			vyt[i] = vy[i];
			vzt[i] = vz[i];

			for(int s = 0; s < S; ++s){
				double dtaa = dt * a_h[S * RKFn + s];
				xt[i]  += dtaa * kx[i + s * N];
				yt[i]  += dtaa * ky[i + s * N];
				zt[i]  += dtaa * kz[i + s * N];
				vxt[i] += dtaa * kvx[i + s * N];
				vyt[i] += dtaa * kvy[i + s * N];
				vzt[i] += dtaa * kvz[i + s * N];
//printf("update 2 %d %d %g %g %g %g %g %g\n", S, i, xt_h[i], yt_h[i], zt_h[i], a_h[S * RKFn + s], kx[s], dt);

			}

			kx[i + S * N] = vxt[i];
			ky[i + S * N] = vyt[i];
			kz[i + S * N] = vzt[i];

		}

		// ----------------------------------------------------------------------------
		//compute forces
		NonGrav(xt, yt, zt, vxt, vyt, vzt);
		GR(xt, yt, zt, vxt, vyt, vzt);
		J2(xt, yt, zt);
		Gravity(xt, yt, zt);
		// ----------------------------------------------------------------------------
		for(int i = Nperturbers; i < N; ++i){
			kvx[i + S * N] = ax[i];
			kvy[i + S * N] = ay[i];
			kvz[i + S * N] = az[i];
		}
	}

	//update
	for(int i = Nperturbers; i < N; ++i){

		dx[i] = 0.0;
		dy[i] = 0.0;
		dz[i] = 0.0;

		dvx[i] = 0.0;
		dvy[i] = 0.0;
		dvz[i] = 0.0;

		for(int S = 0; S < RKFn; ++S){
			double dtb = dt * b_h[S];
			dx[i] += dtb * kx[i + S * N];
			dy[i] += dtb * ky[i + S * N];
			dz[i] += dtb * kz[i + S * N];

			dvx[i] += dtb * kvx[i + S * N];
			dvy[i] += dtb * kvy[i + S * N];
			dvz[i] += dtb * kvz[i + S * N];
		}
	}

	for(int i = Nperturbers; i < N; ++i){
		x[i] += dx[i];
		y[i] += dy[i];
		z[i] += dz[i];

		vx[i] += dvx[i];
		vy[i] += dvy[i];
		vz[i] += dvz[i];

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
			xt[i] = x[i];
			yt[i] = y[i];
			zt[i] = z[i];

			vxt[i] = vx[i];
			vyt[i] = vy[i];
			vzt[i] = vz[i];
		}
		// ----------------------------------------------------------------------------

		for(int i = Nperturbers; i < N; ++i){
			ax[i] = 0.0;
			ay[i] = 0.0;
			az[i] = 0.0;

			xt[i] = x[i];
			yt[i] = y[i];
			zt[i] = z[i];

			vxt[i] = vx[i];
			vyt[i] = vy[i];
			vzt[i] = vz[i];

			for(int s = 0; s < S; ++s){
				double dtaa = dt * a_h[S * RKFn + s];
				xt[i]  += dtaa * kx[i + s * N];
				yt[i]  += dtaa * ky[i + s * N];
				zt[i]  += dtaa * kz[i + s * N];
				vxt[i] += dtaa * kvx[i + s * N];
				vyt[i] += dtaa * kvy[i + s * N];
				vzt[i] += dtaa * kvz[i + s * N];
//printf("update 2 %d %d %g %g %g %g %g %g\n", S, i, xt_h[i], yt_h[i], zt_h[i], a_h[S * RKFn + s], kx[s], dt);

			}

			kx[i + S * N] = vxt[i];
			ky[i + S * N] = vyt[i];
			kz[i + S * N] = vzt[i];

		}

		// ----------------------------------------------------------------------------
		//compute forces
		NonGrav(xt, yt, zt, vxt, vyt, vzt);
		GR(xt, yt, zt, vxt, vyt, vzt);
		J2(xt, yt, zt);
		Gravity(xt, yt, zt);
		// ----------------------------------------------------------------------------
		for(int i = Nperturbers; i < N; ++i){
			kvx[i + S * N] = ax[i];
			kvy[i + S * N] = ay[i];
			kvz[i + S * N] = az[i];
		}
	}

	//update
	for(int i = Nperturbers; i < N; ++i){

		dx[i] = 0.0;
		dy[i] = 0.0;
		dz[i] = 0.0;

		dvx[i] = 0.0;
		dvy[i] = 0.0;
		dvz[i] = 0.0;

		for(int S = 0; S < RKFn; ++S){
			double dtb = dt * b_h[S];
			dx[i] += dtb * kx[i + S * N];
			dy[i] += dtb * ky[i + S * N];
			dz[i] += dtb * kz[i + S * N];

			dvx[i] += dtb * kvx[i + S * N];
			dvy[i] += dtb * kvy[i + S * N];
			dvz[i] += dtb * kvz[i + S * N];
		}
	}


	//compute integration error
	double snew = 10.0;

	for(int i = Nperturbers; i < N; ++i){
		double ym = 0.0;
		ym = (fabs(x[i]) > ym) ? fabs(x[i]) : ym;
		ym = (fabs(y[i]) > ym) ? fabs(y[i]) : ym;
		ym = (fabs(z[i]) > ym) ? fabs(z[i]) : ym;

		ym = (fabs(vx[i]) > ym) ? fabs(vx[i]) : ym;
		ym = (fabs(vy[i]) > ym) ? fabs(vy[i]) : ym;
		ym = (fabs(vz[i]) > ym) ? fabs(vz[i]) : ym;

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
			errorkx += f * kx[i + S * N];
			errorky += f * ky[i + S * N];
			errorkz += f * kz[i + S * N];

			errorkvx += f * kvx[i + S * N];
			errorkvy += f * kvy[i + S * N];
			errorkvz += f * kvz[i + S * N];
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
			x[i] += dx[i];
			y[i] += dy[i];
			z[i] += dz[i];

			vx[i] += dvx[i];
			vy[i] += dvy[i];
			vz[i] += dvz[i];

		}
		time += dt;
	}
	else if(snew >= 1.0){
		//accept step
		for(int i = Nperturbers; i < N; ++i){
			x[i] += dx[i];
			y[i] += dy[i];
			z[i] += dz[i];

			vx[i] += dvx[i];
			vy[i] += dvy[i];
			vz[i] += dvz[i];

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
	//for(int tt = 0; tt < 2; ++tt){
	for(int tt = 0; tt < 2000; ++tt){
		double dtmin = dt;

		double timett1 = timeStart + dts * (tt + 1) * 10.0;

//printf("integrate %.20g %.20g\n", timeStart + dts * tt * 10.0, timett1);


		for(int ttt = 0; ttt < 1000000; ++ttt){

			//refine last time step of interval to match output time
			if(dts < 0.0){
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

			if(time < timeEnd){
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
		for(int p = Nperturbers; p < N; ++p){
			printf("Reached time %.20g %.20g\n", time_reference + time, dt);
			fprintf(outputFile, "%.20g %d %.20g %.20g %.20g %.20g\n", time_reference + time, p, x[p], y[p], z[p], dtmin);
		}

	}//end of tt loop
	fclose(outputFile);
	return 1;
}

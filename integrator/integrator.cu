
#include "asteroid.h"
#include "ChebyshevGPU.h"
#include "forceGPU.h"


//Leapfrog step with fixed time step
__global__ void leapfrog_stepA_kernel(double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, const int Nperturbers, const int N, const double dt){

	int id = blockIdx.x * blockDim.x + threadIdx.x;

	//Drift
	if(id >= Nperturbers && id < N){
		x_d[id] += 0.5 * dt * vx_d[id];
		y_d[id] += 0.5 * dt * vy_d[id];
		z_d[id] += 0.5 * dt * vz_d[id];
	}

}

	// ----------------------------------------------------------------------------
	//Update the Chebyshev coefficients if necessary
	//update_Chebyshev(time);
	//update_perturbers(time);
	// ----------------------------------------------------------------------------

__global__ void leapfrog_stepB_kernel(double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *A1_d, double *A2_d, double *A3_d, double *GM_d, const int Nperturbers, const int N, const double dt, const double c2, const double J2E, const double REAU){

	int id = blockIdx.x * blockDim.x + threadIdx.x;

	double ax = 0.0;
	double ay = 0.0;
	double az = 0.0;

	//compute force

	double xi;
	double yi;
	double zi;

	double vxi;
	double vyi;
	double vzi;

	__shared__ double x_s[def_NP];
	__shared__ double y_s[def_NP];
	__shared__ double z_s[def_NP];

	__shared__ double vx_s[def_NP];
	__shared__ double vy_s[def_NP];
	__shared__ double vz_s[def_NP];

	__shared__ double GM_s[def_NP];

	if(id < Nperturbers){
		x_s[id] = x_d[id];
		y_s[id] = y_d[id];
		z_s[id] = z_d[id];

		vx_s[id] = vx_d[id];
		vy_s[id] = vy_d[id];
		vz_s[id] = vz_d[id];

		GM_s[id] = GM_d[id];
	}

	__syncthreads();


	if(id >= Nperturbers && id < N){

		xi = x_d[id];
		yi = y_d[id];
		zi = z_d[id];

		vxi = vx_d[id];
		vyi = vy_d[id];
		vzi = vz_d[id];

		//heliocentric coordinates
		double xih = xi - x_s[10];
		double yih = yi - y_s[10];
		double zih = zi - z_s[10];

		double vxih = vxi - vx_s[10];
		double vyih = vyi - vy_s[10];
		double vzih = vzi - vz_s[10];

		//r is used in multiple forces, so reuse it
		double rsq = xih * xih + yih * yih + zih * zih;
		double r = sqrt(rsq);

		//Earth centric coordinates
		double xiE = xi - x_s[2];
		double yiE = yi - y_s[2];
		double ziE = zi - z_s[2];

		NonGrav(xih, yih, zih, vxih, vyih, vzih, A1_d[id], A2_d[id], A3_d[id], r, ax, ay, az);
		GR(xih, yih, zih, vxih, vyih, vzih, r, ax, ay, az, GM_s[10], c2);
		J2(xiE, yiE, ziE, ax, ay, az, REAU, J2E, GM_s[2]);
		Gravity(xi, yi, zi, x_s, y_s, z_s, ax, ay, az, GM_s, Nperturbers);
	}
	// ----------------------------------------------------------------------------

	if(id >= Nperturbers && id < N){
		//Kick
		vx_d[id] += dt * ax;
		vy_d[id] += dt * ay;
		vz_d[id] += dt * az;

		//Drift
		x_d[id] += 0.5 * dt * vx_d[id];
		y_d[id] += 0.5 * dt * vy_d[id];
		z_d[id] += 0.5 * dt * vz_d[id];
	}
}

/*
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
		for(int i = Nperturbers; i < N; ++i){

			//heliocentric coordinates
			double xih = xt_h[i] - xt_h[10];
			double yih = yt_h[i] - yt_h[10];
			double zih = zt_h[i] - zt_h[10];

			double vxih = vxt_h[i] - vxt_h[10];
			double vyih = vyt_h[i] - vyt_h[10];
			double vzih = vzt_h[i] - vzt_h[10];

			//r is used in multiple forces, so reuse it
			double rsq = xih * xih + yih * yih + zih * zih;
			double r = sqrt(rsq);

			//Earth centric coordinates
			double xiE = xt_h[i] - xt_h[2];
			double yiE = yt_h[i] - yt_h[2];
			double ziE = zt_h[i] - zt_h[2];

			NonGrav(xih, yih, zih, vxih, vyih, vzih, A1_h[i], A2_h[i], A3_h[i], r, ax_h[i], ay_h[i], az_h[i]);
			GR(xih, yih, zih, vxih, vyih, vzih, r, ax_h[i], ay_h[i], az_h[i], GM_h[10]);
			J2(xiE, yiE, ziE, ax_h[i], ay_h[i], az_h[i], GM_h[2]);
			Gravity(xt_h, yt_h, zt_h, ax_h[i], ay_h[i], az_h[i], i);
		}
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
		for(int i = Nperturbers; i < N; ++i){

			//heliocentric coordinates
			double xih = xt_h[i] - xt_h[10];
			double yih = yt_h[i] - yt_h[10];
			double zih = zt_h[i] - zt_h[10];

			double vxih = vxt_h[i] - vxt_h[10];
			double vyih = vyt_h[i] - vyt_h[10];
			double vzih = vzt_h[i] - vzt_h[10];

			//r is used in multiple forces, so reuse it
			double rsq = xih * xih + yih * yih + zih * zih;
			double r = sqrt(rsq);

			//Earth centric coordinates
			double xiE = xt_h[i] - xt_h[2];
			double yiE = yt_h[i] - yt_h[2];
			double ziE = zt_h[i] - zt_h[2];

			NonGrav(xih, yih, zih, vxih, vyih, vzih, A1_h[i], A2_h[i], A3_h[i], r, ax_h[i], ay_h[i], az_h[i]);
			GR(xih, yih, zih, vxih, vyih, vzih, r, ax_h[i], ay_h[i], az_h[i], GM_h[10]);
			J2(xiE, yiE, ziE, ax_h[i], ay_h[i], az_h[i], GM_h[2]);
			Gravity(xt_h, yt_h, zt_h, ax_h[i], ay_h[i], az_h[i], i);
		}
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

*/
	
int asteroid::loop(){


	outputFile = fopen("Out.dat", "w");
	copyOutput();

	for(int p = Nperturbers; p < N; ++p){
		printf("Start integration\n");
		printf("Reached time %.20g dtmin %.8g\n", time_reference + time, dt);
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
/*
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
*/
			//if(RKFn == 1){
				leapfrog_stepA_kernel <<< (N + 255) / 256 , 256 >>> (x_d, y_d, z_d, vx_d, vy_d, vz_d, Nperturbers, N, dt);
				time += dt * 0.5;
				update_perturbers_kernel <<< 1, 32 >>>(x_d, y_d, z_d, vx_d, vy_d, vz_d, data_d, cdata_d, id_d, startTime_d, endTime_d, nChebyshev_d, offset0_d, time, time_reference, nCm, EM, AUtokm, Nperturbers);

				leapfrog_stepB_kernel <<< (N + 255) / 256 , 256 >>> (x_d, y_d, z_d, vx_d, vy_d, vz_d, A1_d, A2_d, A3_d, GM_d, Nperturbers, N, dt, c2, J2E, REAU);
				time += dt / 2.0;
			//}


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

				printf("Error, time step loop did not finish\n");
				return 0;
			}

		}//end of ttt loop
		copyOutput();
		for(int p = Nperturbers; p < N; ++p){
			printf("Reached time %.20g dtmin %.8g\n", time_reference + time, dt);
			fprintf(outputFile, "%.20g %d %.20g %.20g %.20g %.20g %.20g %.20g %.20g\n", time_reference + time, p, x_h[p], y_h[p], z_h[p], vx_h[p], vy_h[p], vz_h[p], dtmin);
		}

	}//end of tt loop
	fclose(outputFile);
	return 1;
}

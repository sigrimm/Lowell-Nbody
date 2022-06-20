#define dayUnit 0.01720209895
#include "interpolateGPU.h"

//check RKN12(10)17M

__host__ void Host::setRKF45(){

	a_h[1 * 6 + 0] = 1.0/4.0;       //21

	a_h[2 * 6 + 0] = 3.0/32.0;      //31
	a_h[2 * 6 + 1] = 9.0/32.0;      //32

	a_h[3 * 6 + 0] = 1932.0/2197.0; //41
	a_h[3 * 6 + 1] = -7200.0/2197.0;//42
	a_h[3 * 6 + 2] = 7296.0/2197.0; //43

	a_h[4 * 6 + 0] = 439.0/216.0;   //51
	a_h[4 * 6 + 1] = -8.0;          //52
	a_h[4 * 6 + 2] = 3680.0/513.0;  //53
	a_h[4 * 6 + 3] = -845.0/4104.0; //54

	a_h[5 * 6 + 0] = -8.0/27.0;     //61
	a_h[5 * 6 + 1] = 2.0;           //62
	a_h[5 * 6 + 2] = -3544/2565.0;  //63
	a_h[5 * 6 + 3] = 1859.0/4104.0; //64
	a_h[5 * 6 + 4] = -11.0/40.0;    //65


	b_h[0] = 25.0/216.0;
	b_h[1] = 0.0;
	b_h[2] = 1408.0/2565.0;
	b_h[3] = 2197.0/4104.0;
	b_h[4] = -1.0/5.0;
	b_h[5] = 0.0;

	bb_h[0] = 16.0/135.0;
	bb_h[1] = 0.0;
	bb_h[2] = 6656.0/12825.0;
	bb_h[3] = 28561.0/56430.0;
	bb_h[4] = -9.0/50.0;
	bb_h[5] = 2.0/55.0;

	c_h[0] = 0.0;
	c_h[1] = 0.25;
	c_h[2] = 3.0 / 8.0;
	c_h[3] = 12.0 / 13.0;
	c_h[4] = 1.0;
	c_h[5] = 0.5;

	ee = 1.0/5.0;

}

__host__ void Host::setDP54(){

	a_h[1 * 7 + 0] = 1.0/5.0;       //21

	a_h[2 * 7 + 0] = 3.0/40.0;      //31
	a_h[2 * 7 + 1] = 9.0/40.0;      //32

	a_h[3 * 7 + 0] = 44.0/45.0; //41
	a_h[3 * 7 + 1] = -56.0/15.0;//42
	a_h[3 * 7 + 2] = 32.0/9.0; //43

	a_h[4 * 7 + 0] = 19372.0/6561.0;   //51
	a_h[4 * 7 + 1] = -25360.0/2187.0;          //52
	a_h[4 * 7 + 2] = 64448.0/6561.0;  //53
	a_h[4 * 7 + 3] = -212.0/729.0; //54

	a_h[5 * 7 + 0] = 9017.0/3168.0;  //61
	a_h[5 * 7 + 1] = -355.0/33.0;           //62
	a_h[5 * 7 + 2] = 46732.0/5247.0;  //63
	a_h[5 * 7 + 3] = 49.0/176.0; //64
	a_h[5 * 7 + 4] = -5103.0/18656.0;    //65

	a_h[6 * 7 + 0] = 35.0/384.0;  //61
	a_h[6 * 7 + 1] = 0.0;           //62
	a_h[6 * 7 + 2] = 500.0/1113.0;  //63
	a_h[6 * 7 + 3] = 125.0/192.0; //64
	a_h[6 * 7 + 4] = -2187.0/6784.0;    //65
	a_h[6 * 7 + 5] = 11.0/84.0;    //65

	b_h[0] = 35.0/384.0;
	b_h[1] = 0.0;
	b_h[2] = 500.0/1113.0;
	b_h[3] = 125.0/192.0;
	b_h[4] = -2187.0/6784.0;
	b_h[5] = 11.0/84.0;
	b_h[6] = 0.0;

	bb_h[0] = 5179.0/57600.0;
	bb_h[1] = 0.0;
	bb_h[2] = 7571.0/16695.0;
	bb_h[3] = 393.0/640.0;
	bb_h[4] = -92097.0/339200.0;
	bb_h[5] = 187.0/2100.0;
	bb_h[6] = 1.0/40.0;

	c_h[0] = 0.0;
	c_h[1] = 1.0/5.0;
	c_h[2] = 3.0 / 10.0;
	c_h[3] = 4.0 / 5.0;
	c_h[4] = 8.0/9.0;
	c_h[5] = 1.0;
	c_h[6] = 1.0;

	ee = 1.0/5.0;

}


__host__ void Host::setRKF78(){


	a_h[1 * 13 + 0] = 2.0/27.0;

	a_h[2 * 13 + 0] = 1.0/36.0;
	a_h[2 * 13 + 1] = 1.0/12.0;

	a_h[3 * 13 + 0] = 1.0/24.0;
	a_h[3 * 13 + 1] = 0.0;
	a_h[3 * 13 + 2] = 1.0/8.0;

	a_h[4 * 13 + 0] = 5.0/12.0;
	a_h[4 * 13 + 1] = 0.0;
	a_h[4 * 13 + 2] = -25.0/16.0;
	a_h[4 * 13 + 3] = 25.0/16.0;

	a_h[5 * 13 + 0] = 1.0/20.0;
	a_h[5 * 13 + 1] = 0.0;
	a_h[5 * 13 + 2] = 0.0;
	a_h[5 * 13 + 3] = 1.0/4.0;
	a_h[5 * 13 + 4] = 1.0/5.0;

	a_h[6 * 13 + 0] = -25.0/ 108.0;
	a_h[6 * 13 + 1] = 0.0;
	a_h[6 * 13 + 2] = 0.0;
	a_h[6 * 13 + 3] = 125.0/ 108.0;
	a_h[6 * 13 + 4] = -65.0/ 27.0;
	a_h[6 * 13 + 5] = 125.0/ 54.0;

	a_h[7 * 13 + 0] = 31.0/300.0;
	a_h[7 * 13 + 1] = 0.0;
	a_h[7 * 13 + 2] = 0.0;
	a_h[7 * 13 + 3] = 0.0;
	a_h[7 * 13 + 4] = 61.0/225.0;
	a_h[7 * 13 + 5] = -2.0/9.0;
	a_h[7 * 13 + 6] = 13.0/900.0;

	a_h[8 * 13 + 0] = 2.0;
	a_h[8 * 13 + 1] = 0.0;
	a_h[8 * 13 + 2] = 0.0;
	a_h[8 * 13 + 3] = -53.0/6.0;
	a_h[8 * 13 + 4] = 704.0/45.0;
	a_h[8 * 13 + 5] = -107.0/9.0;
	a_h[8 * 13 + 6] = 67.0/90.0;
	a_h[8 * 13 + 7] = 3.0;

	a_h[9 * 13 + 0] = -91.0/108.0;
	a_h[9 * 13 + 1] = 0.0;
	a_h[9 * 13 + 2] = 0.0;
	a_h[9 * 13 + 3] = 23.0/108.0;
	a_h[9 * 13 + 4] = -976.0/135.0;
	a_h[9 * 13 + 5] = 311.0/54.0;
	a_h[9 * 13 + 6] = -19.0/60.0;
	a_h[9 * 13 + 7] = 17.0/6.0;
	a_h[9 * 13 + 8] = -1.0/12.0;

	a_h[10 * 13 + 0] = 2383.0/4100.0;
	a_h[10 * 13 + 1] = 0.0;
	a_h[10 * 13 + 2] = 0.0;
	a_h[10 * 13 + 3] = -341.0/164.0;
	a_h[10 * 13 + 4] = 4496.0/1025.0;
	a_h[10 * 13 + 5] = -301.0/82.0;
	a_h[10 * 13 + 6] = 2133.0/4100.0;
	a_h[10 * 13 + 7] = 45.0/82.0;
	a_h[10 * 13 + 8] = 45.0/164.0;
	a_h[10 * 13 + 9] = 18.0/41.0;

	a_h[11 * 13 + 0] = 3.0/205.0;
	a_h[11 * 13 + 1] = 0.0;
	a_h[11 * 13 + 2] = 0.0;
	a_h[11 * 13 + 3] = 0.0;
	a_h[11 * 13 + 4] = 0.0;
	a_h[11 * 13 + 5] = -6.0/41.0;
	a_h[11 * 13 + 6] = -3.0/205.0;
	a_h[11 * 13 + 7] = -3.0/41.0;
	a_h[11 * 13 + 8] = 3.0/41.0;
	a_h[11 * 13 + 9] = 6.0/41.0;
	a_h[11 * 13 + 10] = 0.0;

	a_h[12 * 13 + 0] = -1777.0/4100.0;
	a_h[12 * 13 + 1] = 0.0;
	a_h[12 * 13 + 2] = 0.0;
	a_h[12 * 13 + 3] = -341.0/164.0;
	a_h[12 * 13 + 4] = 4496.0/1025.0;
	a_h[12 * 13 + 5] = -289.0/82.0;
	a_h[12 * 13 + 6] = 2193.0/4100.0;
	a_h[12 * 13 + 7] = 51.0/82.0;
	a_h[12 * 13 + 8] = 33.0/164.0;
	a_h[12 * 13 + 9] = 12.0/41.0;		//19 / 41 or 12/41 ?
	a_h[12 * 13 + 10] = 0.0;
	a_h[12 * 13 + 11] = 1.0;		//0 or 1 ?



	b_h[0] = 41.0/840.0;
	b_h[1] = 0.0;
	b_h[2] = 0.0;
	b_h[3] = 0.0;
	b_h[4] = 0.0;
	b_h[5] = 34.0/105.0;
	b_h[6] = 9.0/35.0;
	b_h[7] = 9.0/35.0;
	b_h[8] = 9.0/280.0;
	b_h[9] = 9.0/280.0;
	b_h[10] = 41.0/840.0;
	b_h[11] = 0.0;
	b_h[12] = 0.0;

	bb_h[0] = 0.0;
	bb_h[1] = 0.0;
	bb_h[2] = 0.0;
	bb_h[3] = 0.0;
	bb_h[4] = 0.0;
	bb_h[5] = 34.0/105.0;
	bb_h[6] = 9.0/35.0;
	bb_h[7] = 9.0/35.0;
	bb_h[8] = 9.0/280.0;
	bb_h[9] = 9.0/280.0;
	bb_h[10] = 0.0;
	bb_h[11] = 41.0/840.0;
	bb_h[12] = 41.0/840.0;


	c_h[0] = 0.0;
	c_h[1] = 2.0/27.0;
	c_h[2] = 1.0/9.0;
	c_h[3] = 1.0/6.0;
	c_h[4] = 5.0/12.0;
	c_h[5] = 1.0/2.0;
	c_h[6] = 5.0/6.0; 
	c_h[7] = 1.0/6.0;
	c_h[8] = 2.0/3.0;
	c_h[9] = 1.0/3.0;
	c_h[10] = 1.0;
	c_h[11] = 0.0;
	c_h[12] = 1.0;

	ee = 1.0/8.0;
/*
	
	for(int i = 0; i < 13; ++i){
		for(int j = 0; j < i + 1; ++j){

			printf("%d %d %.20g\n", i, j, a_h[i * 13 + j]);
		}
	}
	for(int i = 0; i < 13; ++i){
		printf("%d %.20g\n", i, b_h[i]);
	}
	for(int i = 0; i < 13; ++i){
		printf("%d %.20g\n", i, bb_h[i]);
	}
	for(int i = 0; i < 13; ++i){
		printf("%d %.20g\n", i, c_h[i]);
	}
*/	
}

__host__ void Host::copyConst(){
	cudaMemcpyToSymbol(a_c, a_h, RKFn * RKFn * sizeof(double), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(b_c, b_h, RKFn * sizeof(double), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(bb_c, bb_h, RKFn * sizeof(double), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_c, c_h, RKFn * sizeof(double), 0, cudaMemcpyHostToDevice);
}

//only 1 particle
__host__ void Host::stageStep2(int i, double time, double dtiMin, double &snew){


	double kx[RKFn];
	double ky[RKFn];
	double kz[RKFn];
	double kvx[RKFn];
	double kvy[RKFn];
	double kvz[RKFn];

	for(int S = 0; S < RKFn; ++S){

		for(int p = 0; p < Nperturbers; ++p){
			interpolate(time + c_h[S] * dti, p);
			//interpolate2(time + c_h[S] * dti, p);
		}
		//update
		xt_h[i]  = x_h[i];
		yt_h[i]  = y_h[i];
		zt_h[i]  = z_h[i];
		vxt_h[i] = vx_h[i];
		vyt_h[i] = vy_h[i];
		vzt_h[i] = vz_h[i];

		for(int s = 0; s < S; ++s){
			double dtaa = dt * a_h[S * RKFn + s];
			xt_h[i]  += dtaa * kx[s];
			yt_h[i]  += dtaa * ky[s];
			zt_h[i]  += dtaa * kz[s];
			vxt_h[i] += dtaa * kvx[s];
			vyt_h[i] += dtaa * kvy[s];
			vzt_h[i] += dtaa * kvz[s];
//printf("update 2 %d %d %g %g %g %g %g %g\n", S, i, xt_h[i], yt_h[i], zt_h[i], a_h[S * RKFn + s], kx[s], dt);

		}

		kx[S] = vxt_h[i];
		ky[S] = vyt_h[i];
		kz[S] = vzt_h[i];

		double ax = 0.0;
		double ay = 0.0;
		double az = 0.0;

		if(useHelio == 0){
			for(int j = Nperturbers-1; j >= 0; --j){
				if(id_h[i] != id_h[j]){ 
					accP(m_h[j], xt_h[j], yt_h[j], zt_h[j], xt_h[i], yt_h[i], zt_h[i], ax, ay, az);
				}
			}
		}
		else{
			for(int j = Nperturbers-1; j >= 1; --j){
				if(id_h[i] != id_h[j]){
					accP(m_h[j], xt_h[j], yt_h[j], zt_h[j], xt_h[i], yt_h[i], zt_h[i], ax, ay, az);
	//if(i == 27) printf("Nij %d %d %llu %llu %.20g %.20g %.20g\n", i, j, id_h[i], id_h[j], ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit);
				}
			}
			accS(m_h[0] + m_h[i], xt_h[i], yt_h[i], zt_h[i], ax, ay, az);
	//if(i == 27) printf("Nij %d %d %.20g %.20g %.20g\n", i, 0, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit);
//printf("N0 %d %.20g %.20g %.20g\n", i, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit);
			for(int j = Nperturbers-1; j >= 1; --j){
				if(id_h[i] != id_h[j]){
					accP2(m_h[j], xt_h[j], yt_h[j], zt_h[j], ax, ay, az);
				}
			}
		//printf("Np %d %.20g %.20g %.20g %d\n", i, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit, S);
		}

		if(useGR == 2){
			acchGR2(xt_h[i], yt_h[i], zt_h[i], vxt_h[i], vyt_h[i], vzt_h[i], ax, ay, az);
		}

		if(useNonGrav == 1){
			NonGrav(xt_h[i], yt_h[i], zt_h[i], vxt_h[i], vyt_h[i], vzt_h[i], ax, ay, az, A1_h[i], A2_h[i], A3_h[i]);
		}

		if(useJ2 == 1){
			J2(m_h[3], xt_h[3], yt_h[3], zt_h[3], xt_h[i], yt_h[i], zt_h[i], ax, ay, az);
		}

		kvx[S] = ax;
		kvy[S] = ay;
		kvz[S] = az;
	}

	//update
	double dx = 0.0;
	double dy = 0.0;
	double dz = 0.0;
	double dvx = 0.0;
	double dvy = 0.0;
	double dvz = 0.0;

	for(int S = 0; S < RKFn; ++S){
		double dtb = dt * b_h[S];
		dx += dtb * kx[S];
		dy += dtb * ky[S];
		dz += dtb * kz[S];

		dvx += dtb * kvx[S];
		dvy += dtb * kvy[S];
		dvz += dtb * kvz[S];
	}

	dx_h[i] = dx;
	dy_h[i] = dy;
	dz_h[i] = dz;
	dvx_h[i] = dvx;
	dvy_h[i] = dvy;
	dvz_h[i] = dvz;

//if(i == 27) printf("dx %d %.20g %.20g %.20g %.20g %.20g %.20g\n", i, dx, dy, dz, dvx, dvy, dvz);

	if(useAdaptiveTimeSteps == 1){
		//compute error
		double ym = fabs(x_h[i]);
		ym = fmax(ym, fabs(y_h[i]));
		ym = fmax(ym, fabs(z_h[i]));
		ym = fmax(ym, fabs(vx_h[i]));
		ym = fmax(ym, fabs(vy_h[i]));
		ym = fmax(ym, fabs(vz_h[i]));

		ym = fmax(ym, fabs(x_h[i] + dx));
		ym = fmax(ym, fabs(y_h[i] + dy));
		ym = fmax(ym, fabs(z_h[i] + dz));
		ym = fmax(ym, fabs(vx_h[i] + dvx));
		ym = fmax(ym, fabs(vy_h[i] + dvy));
		ym = fmax(ym, fabs(vz_h[i] + dvz));

		double isc = 1.0 / (def_atol + ym * def_rtol);

		double s = 1.0e6;	//large number
		snew = 1.0e6;		//large number

		//error estimation
		double errorkx = 0.0;
		double errorky = 0.0;
		double errorkz = 0.0;
		double errorkvx = 0.0;
		double errorkvy = 0.0;
		double errorkvz = 0.0;

		for(int S = 0; S < RKFn; ++S){
			double f = (b_h[S] - bb_h[S]) * dt;
			errorkx += f * kx[S];
			errorky += f * ky[S];
			errorkz += f * kz[S];

			errorkvx += f * kvx[S];
			errorkvy += f * kvy[S];
			errorkvz += f * kvz[S];
//printf("error %d %d %g %g\n", i, S, errorkx, kx[S]);
		}

		double errork = 0.0;
		errork += errorkx * errorkx * isc * isc;
		errork += errorky * errorky * isc * isc;
		errork += errorkz * errorkz * isc * isc;
		errork += errorkvx * errorkvx * isc * isc;
		errork += errorkvy * errorkvy * isc * isc;
		errork += errorkvz * errorkvz * isc * isc;

		errork = sqrt(errork / 6.0);	//6 is the number of dimensions

		s = pow( 1.0  / errork, ee);
//printf("error %d %g %g %g\n", i, errork, s, errorkx);
		s = fmax(def_facmin, def_fac * s);
		s = fmin(def_facmax, s);

		// ****************************
		if(fabs(s * dti) < dtiMin){
			snew_h[i].y = fmin(snew_h[i].y, s);
		}
		if(snew_h[i].y >= 1.0){
			snew_h[i].x = s;
		}
		else{
			snew_h[i].x = 1.5;
		}

		snew = fmin(snew, snew_h[i].x);
//printf("snew %d %.20g %.20g %.20g %.20g %g\n", i, snew_h[i].y, snew_h[i].x, errork, snew, s);
	}

}

//N particles
__host__ void Host::stageStep(double time, double dtiMin, double &snew){


	for(int S = 0; S < RKFn; ++S){

		for(int p = 0; p < Nperturbers; ++p){
			interpolate(time + c_h[S] * dti, p);
			//interpolate2(time + c_h[S] * dti, p);
		}

		for(int i = Nperturbers; i < N; ++i){

			//update
			xt_h[i]  = x_h[i];
			yt_h[i]  = y_h[i];
			zt_h[i]  = z_h[i];
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
	//printf("update 2 %d %d %g %g %g %g %g %g\n", S, i, xt_h[i], yt_h[i], zt_h[i], a_h[S * RKFn + s], kx[s], dt);

			}

			kx_h[i + S * N] = vxt_h[i];
			ky_h[i + S * N] = vyt_h[i];
			kz_h[i + S * N] = vzt_h[i];

			double ax = 0.0;
			double ay = 0.0;
			double az = 0.0;

			if(useHelio == 0){
				for(int j = Nperturbers-1; j >= 0; --j){
					if(id_h[i] != id_h[j]){ 
						accP(m_h[j], xt_h[j], yt_h[j], zt_h[j], xt_h[i], yt_h[i], zt_h[i], ax, ay, az);
					}
				}
			}
			else{
				for(int j = Nperturbers-1; j >= 1; --j){
					if(id_h[i] != id_h[j]){
						accP(m_h[j], xt_h[j], yt_h[j], zt_h[j], xt_h[i], yt_h[i], zt_h[i], ax, ay, az);
//if(i == 27) printf("Nij %d %d %llu %llu %.20g %.20g %.20g\n", i, j, id_h[i], id_h[j], ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit);
					}
				}
				accS(m_h[0] + m_h[i], xt_h[i], yt_h[i], zt_h[i], ax, ay, az);
//if(i == 27) printf("Nij %d %d %.20g %.20g %.20g\n", i, 0, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit);
//printf("N0 %d %.20g %.20g %.20g\n", i, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit);
				for(int j = Nperturbers-1; j >= 1; --j){
					if(id_h[i] != id_h[j]){
						accP2(m_h[j], xt_h[j], yt_h[j], zt_h[j], ax, ay, az);
					}
				}
//printf("Np %d %.20g %.20g %.20g %d\n", i, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit, S);
			}

			if(useGR == 2){
				acchGR2(xt_h[i], yt_h[i], zt_h[i], vxt_h[i], vyt_h[i], vzt_h[i], ax, ay, az);
			}

			if(useNonGrav == 1){
				NonGrav(xt_h[i], yt_h[i], zt_h[i], vxt_h[i], vyt_h[i], vzt_h[i], ax, ay, az, A1_h[i], A2_h[i], A3_h[i]);
			}

			if(useJ2 == 1){
				J2(m_h[3], xt_h[3], yt_h[3], zt_h[3], xt_h[i], yt_h[i], zt_h[i], ax, ay, az);
			}

			kvx_h[i + S * N] = ax;
			kvy_h[i + S * N] = ay;
			kvz_h[i + S * N] = az;
		}
	}

	snew = 1.0e6;   //large number
	for(int i = Nperturbers; i < N; ++i){

		//update
		double dx = 0.0;
		double dy = 0.0;
		double dz = 0.0;
		double dvx = 0.0;
		double dvy = 0.0;
		double dvz = 0.0;

		for(int S = 0; S < RKFn; ++S){
			double dtb = dt * b_h[S];
			dx += dtb * kx_h[i + S * N];
			dy += dtb * ky_h[i + S * N];
			dz += dtb * kz_h[i + S * N];

			dvx += dtb * kvx_h[i + S * N];
			dvy += dtb * kvy_h[i + S * N];
			dvz += dtb * kvz_h[i + S * N];
		}

		dx_h[i] = dx;
		dy_h[i] = dy;
		dz_h[i] = dz;
		dvx_h[i] = dvx;
		dvy_h[i] = dvy;
		dvz_h[i] = dvz;

//if(i == 27) printf("dx %d %.20g %.20g %.20g %.20g %.20g %.20g\n", i, dx, dy, dz, dvx, dvy, dvz);

		if(useAdaptiveTimeSteps == 1){

			//compute error
			double ym = fabs(x_h[i]);
			ym = fmax(ym, fabs(y_h[i]));
			ym = fmax(ym, fabs(z_h[i]));
			ym = fmax(ym, fabs(vx_h[i]));
			ym = fmax(ym, fabs(vy_h[i]));
			ym = fmax(ym, fabs(vz_h[i]));

			ym = fmax(ym, fabs(x_h[i] + dx));
			ym = fmax(ym, fabs(y_h[i] + dy));
			ym = fmax(ym, fabs(z_h[i] + dz));
			ym = fmax(ym, fabs(vx_h[i] + dvx));
			ym = fmax(ym, fabs(vy_h[i] + dvy));
			ym = fmax(ym, fabs(vz_h[i] + dvz));

			double isc = 1.0 / (def_atol + ym * def_rtol);

			double s = 1.0e6;	//large number

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

			errork = sqrt(errork / 6.0);	//6 is the number of dimensions

			s = pow( 1.0  / errork, ee);
//printf("error %d %g %g %g\n", i, errork, s, errorkx);
			s = fmax(def_facmin, def_fac * s);
			s = fmin(def_facmax, s);

			// ****************************
			if(fabs(s * dti) < dtiMin){
				snew_h[i].y = fmin(snew_h[i].y, s);
			}
			if(snew_h[i].y >= 1.0){
				snew_h[i].x = s;
			}
			else{
				snew_h[i].x = 1.5;
			}

			snew = fmin(snew, snew_h[i].x);
//printf("snew %d %.20g %.20g %.20g %.20g %g\n", i, snew_h[i].y, snew_h[i].x, errork, snew, s);
		}
	}
}

//N particles, using interpolation table
__host__ void Host::stageStep1(double dtiMin, int iStart, int N, double &snew){


	double xt[Nperturbers];
	double yt[Nperturbers];
	double zt[Nperturbers];

	double kx[RKFn];
	double ky[RKFn];
	double kz[RKFn];
	double kvx[RKFn];
	double kvy[RKFn];
	double kvz[RKFn];

	snew = 1.0e6;	//large number

	for(int i = Nperturbers + iStart; i < N; ++i){
		for(int S = 0; S < RKFn; ++S){

			for(int p = 0; p < Nperturbers; ++p){
				int ii = p * RKFn + S;
				xt[p] = xTable_h[ii];
				yt[p] = yTable_h[ii];
				zt[p] = zTable_h[ii];
//printf("%d %d %.20g\n", i, S, xt[i]);
			}


			//update
			double xi = x_h[i];
			double yi = y_h[i];
			double zi = z_h[i];
			double vxi = vx_h[i];
			double vyi = vy_h[i];
			double vzi = vz_h[i];

			for(int s = 0; s < S; ++s){
				double dtaa = dt * a_h[S * RKFn + s];
				xi  += dtaa * kx[s];
				yi  += dtaa * ky[s];
				zi  += dtaa * kz[s];
				vxi += dtaa * kvx[s];
				vyi += dtaa * kvy[s];
				vzi += dtaa * kvz[s];
//printf("x %d %d %.20g %.20g\n", S, s, dtaa, xi);
			}


			kx[S] = vxi;
			ky[S] = vyi;
			kz[S] = vzi;

			double ax = 0.0;
			double ay = 0.0;
			double az = 0.0;

			if(useHelio == 0){
				for(int j = Nperturbers-1; j >= 0; --j){
					if(id_h[i] != id_h[j]){
						accP(m_h[j], xt[j], yt[j], zt[j], xi, yi, zi, ax, ay, az);
					}
				}
			}
			else{
				for(int j = Nperturbers-1; j >= 1; --j){
					if(id_h[i] != id_h[j]){
						accP(m_h[j], xt[j], yt[j], zt[j], xi, yi, zi, ax, ay, az);
//printf("Nij %d %d %llu %llu %.20g %.20g %.20g\n", i, j, id_h[i], id_h[j], ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit);
					}
				}
				accS(m_h[0] + m_h[i], xi, yi, zi, ax, ay, az);
//if(i == 27) printf("Nij %d %d %.20g %.20g %.20g\n", i, 0, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit);
//printf("N0 %d %.20g %.20g %.20g\n", i, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit);
				for(int j = Nperturbers-1; j >= 1; --j){
					if(id_h[i] != id_h[j]){
						accP2(m_h[j], xt[j], yt[j], zt[j], ax, ay, az);
					}
				}
//printf("Np %d %.20g %.20g %.20g %d\n", i, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit, S);
			}


			if(useGR == 2){
				acchGR2(xi, yi, zi, vxi, vyi, vzi, ax, ay, az);
			}

			if(useNonGrav == 1){
				NonGrav(xi, yi, zi, vxi, vyi, vzi, ax, ay, az, A1_h[i], A2_h[i], A3_h[i]);
			}

			if(useJ2 == 1){
				J2(m_h[3], xt[3], yt[3], zt[3], xi, yi, zi, ax, ay, az);
			}

			kvx[S] = ax;
			kvy[S] = ay;
			kvz[S] = az;
		}

		//udate
		double dx = 0.0;
		double dy = 0.0;
		double dz = 0.0;
		double dvx = 0.0;
		double dvy = 0.0;
		double dvz = 0.0;

		for(int S = 0; S < RKFn; ++S){
			double dtb = dt * b_h[S];
			dx += dtb * kx[S];
			dy += dtb * ky[S];
			dz += dtb * kz[S];

			dvx += dtb * kvx[S];
			dvy += dtb * kvy[S];
			dvz += dtb * kvz[S];
		}

		dx_h[i] = dx;
		dy_h[i] = dy;
		dz_h[i] = dz;
		dvx_h[i] = dvx;
		dvy_h[i] = dvy;
		dvz_h[i] = dvz;

		if(useAdaptiveTimeSteps == 1){
				
			//compute error
			double ym = fabs(x_h[i]);
			ym = fmax(ym, fabs(y_h[i]));
			ym = fmax(ym, fabs(z_h[i]));
			ym = fmax(ym, fabs(vx_h[i]));
			ym = fmax(ym, fabs(vy_h[i]));
			ym = fmax(ym, fabs(vz_h[i]));

			ym = fmax(ym, fabs(x_h[i] + dx));
			ym = fmax(ym, fabs(y_h[i] + dy));
			ym = fmax(ym, fabs(z_h[i] + dz));
			ym = fmax(ym, fabs(vx_h[i] + dvx));
			ym = fmax(ym, fabs(vy_h[i] + dvy));
			ym = fmax(ym, fabs(vz_h[i] + dvz));

			double isc = 1.0 / (def_atol + ym * def_rtol);

			double s = 1.0e6;	//large number

			double errorkx = 0.0;
			double errorky = 0.0;
			double errorkz = 0.0;
			double errorkvx = 0.0;
			double errorkvy = 0.0;
			double errorkvz = 0.0;

			for(int S = 0; S < RKFn; ++S){
				//this is y1i - y^1i
				double f = (b_h[S] - bb_h[S]) * dt;
				errorkx += f * kx[S];
				errorky += f * ky[S];
				errorkz += f * kz[S];

				errorkvx += f * kvx[S];
				errorkvy += f * kvy[S];
				errorkvz += f * kvz[S];
//printf("error %d %d %g %g\n", i, S, errorkx, kx[S]);

			}

			double errork = 0.0;
			errork += errorkx * errorkx * isc * isc;
			errork += errorky * errorky * isc * isc;
			errork += errorkz * errorkz * isc * isc;
			errork += errorkvx * errorkvx * isc * isc;
			errork += errorkvy * errorkvy * isc * isc;
			errork += errorkvz * errorkvz * isc * isc;

			errork = sqrt(errork / 6.0);	//6 is the number of dimensions

			s = pow( 1.0  / errork, ee);
//printf("error %d %g %g %g\n", i, errork, s, errorkx);
			s = fmax(def_facmin, def_fac * s);
			s = fmin(def_facmax, s);

			// ****************************
			if(fabs(s * dti) < dtiMin){
				snew_h[i].y = fmin(snew_h[i].y, s);
			}
			if(snew_h[i].y >= 1.0){
				snew_h[i].x = s;
			}
			else{
				snew_h[i].x = 1.5;
			}

			snew = fmin(snew, snew_h[i].x);
//printf("snew %d %.20g %.20g %.20g %.20g %g\n", i, snew_h[i].y, snew_h[i].x, errork, snew, s);
			// ****************************
		}
	}

}


__global__ void stageStep1_kernel(unsigned long long int *id_d, double *m_d, double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *dx_d, double *dy_d, double *dz_d, double *dvx_d, double *dvy_d, double *dvz_d, double *xTable_d, double *yTable_d, double *zTable_d, double *kx_d, double *ky_d, double *kz_d, double *kvx_d, double *kvy_d, double *kvz_d, double *A1_d, double *A2_d, double *A3_d, double2 *snew_d, double dt, double dti, double dtiMin, int N, int useHelio, int useGR, int useJ2, int useNonGrav, int useAdaptiveTimeSteps, double ee){

	int itx = threadIdx.x;
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	//shared memory contains only the perturbers
	//the particle idx is stored in registers
	__shared__ double x_s[Nperturbers];
	__shared__ double y_s[Nperturbers];
	__shared__ double z_s[Nperturbers];
	__shared__ double m_s[Nperturbers];
	__shared__ unsigned long long int id_s[Nperturbers];

	double xi0;
	double yi0;
	double zi0;
	double vxi0;
	double vyi0;
	double vzi0;

	double mi;
	double A1i;
	double A2i;
	double A3i;
	unsigned long long int idi;

	if(idx < N){
		xi0 = x_d[idx];
		yi0 = y_d[idx];
		zi0 = z_d[idx];
		vxi0 = vx_d[idx];
		vyi0 = vy_d[idx];
		vzi0 = vz_d[idx];

		mi = m_d[idx];
		A1i = A1_d[idx];
		A2i = A2_d[idx];
		A3i = A3_d[idx];
		idi = id_d[idx];
	}

	for(int S = 0; S < RKFn; ++S){


		if(itx < Nperturbers){
			int ii = itx * RKFn + S;
			x_s[itx] = xTable_d[ii];
			y_s[itx] = yTable_d[ii];
			z_s[itx] = zTable_d[ii];
			m_s[itx] = m_d[itx];
			id_s[itx] = id_d[itx];
		}
		__syncthreads();

		if(idx >= Nperturbers && idx < N){


			// ***********************
			// update
			double xi = xi0;
			double yi = yi0;
			double zi = zi0;
			double vxi = vxi0;
			double vyi = vyi0;
			double vzi = vzi0;

			for(int s = 0; s < S; ++s){
				double dtaa = dt * a_c[S * RKFn + s];
				xi  += dtaa * kx_d[idx + s * N];
				yi  += dtaa * ky_d[idx + s * N];
				zi  += dtaa * kz_d[idx + s * N];
				vxi += dtaa * kvx_d[idx + s * N];
				vyi += dtaa * kvy_d[idx + s * N];
				vzi += dtaa * kvz_d[idx + s * N];
//printf("x %d %d %.20g %.20g\n", S, s, dtaa, xi);
			}



			// end update
			// *****************************

			kx_d[idx + S * N] = vxi;
			ky_d[idx + S * N] = vyi;
			kz_d[idx + S * N] = vzi;

			double ax = 0.0;
			double ay = 0.0;
			double az = 0.0;

			if(useHelio == 0){
				for(int j = Nperturbers-1; j >= 0; --j){
					if(idi != id_s[j]){
						accP(m_s[j], x_s[j], y_s[j], z_s[j], xi, yi, zi, ax, ay, az);
					}
				}
			}
			else{
				for(int j = Nperturbers-1; j >= 1; --j){
					if(idi != id_s[j]){
						accP(m_s[j], x_s[j], y_s[j], z_s[j], xi, yi, zi, ax, ay, az);
//printf("Nij %d %d %.20g %.20g %.20g\n", idx, j, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit);
					}
				}
				accS(m_s[0] + mi, xi, yi, zi, ax, ay, az);

//if(idx == 27) printf("N0 %d %.20g %.20g %.20g\n", idx, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit);
				for(int j = Nperturbers-1; j >= 1; --j){
					if(idi != id_s[j]){
						accP2(m_s[j], x_s[j], y_s[j], z_s[j], ax, ay, az);
					}
				}
			}


			if(useGR == 2){
				acchGR2(xi, yi, zi, vxi, vyi, vzi, ax, ay, az);
			}

			if(useNonGrav == 1){
				NonGrav(xi, yi, zi, vxi, vyi, vzi, ax, ay, az, A1i, A2i, A3i);
			}

			if(useJ2 == 1){
				J2(m_s[3], x_s[3], y_s[3], z_s[3], xi, yi, zi, ax, ay, az);
			}

			kvx_d[idx + S * N] = ax;
			kvy_d[idx + S * N] = ay;
			kvz_d[idx + S * N] = az;
		}
		__syncthreads();
	}

	if(idx >= Nperturbers && idx < N){
		//update
		double dx = 0.0;
		double dy = 0.0;
		double dz = 0.0;
		double dvx = 0.0;
		double dvy = 0.0;
		double dvz = 0.0;
		
		for(int S = 0; S < RKFn; ++S){
			double dtb = dt * b_c[S];
			dx += dtb * kx_d[idx + S * N];
			dy += dtb * ky_d[idx + S * N];
			dz += dtb * kz_d[idx + S * N];

			dvx += dtb * kvx_d[idx + S * N];
			dvy += dtb * kvy_d[idx + S * N];
			dvz += dtb * kvz_d[idx + S * N];
		}
		
		dx_d[idx] = dx;
		dy_d[idx] = dy;
		dz_d[idx] = dz;
		dvx_d[idx] = dvx;
		dvy_d[idx] = dvy;
		dvz_d[idx] = dvz;

//if(idx == 27) printf("dx %d %.20g %.20g %.20g %.20g %.20g %.20g\n", idx, dx, dy, dz, dvx, dvy, dvz);

		if(useAdaptiveTimeSteps == 1){
			//compute error
			double ym = fabs(xi0);
			ym = fmax(ym, fabs(yi0));
			ym = fmax(ym, fabs(zi0));
			ym = fmax(ym, fabs(vxi0));
			ym = fmax(ym, fabs(vyi0));
			ym = fmax(ym, fabs(vzi0));

			ym = fmax(ym, fabs(xi0 + dx));
			ym = fmax(ym, fabs(yi0 + dy));
			ym = fmax(ym, fabs(zi0 + dz));
			ym = fmax(ym, fabs(vxi0 + dvx));
			ym = fmax(ym, fabs(vyi0 + dvy));
			ym = fmax(ym, fabs(vzi0 + dvz));


			double isc = 1.0 / (def_atol + ym * def_rtol);

			double s = 1.0e6;	//large number

			double errorkx = 0.0;
			double errorky = 0.0;
			double errorkz = 0.0;
			double errorkvx = 0.0;
			double errorkvy = 0.0;
			double errorkvz = 0.0;

			for(int S = 0; S < RKFn; ++S){
				//this is y1i - y^1i
				double f = (b_c[S] - bb_c[S]) * dt;
				errorkx += f * kx_d[S * N + idx];
				errorky += f * ky_d[S * N + idx];
				errorkz += f * kz_d[S * N + idx];

				errorkvx += f * kvx_d[S * N + idx];
				errorkvy += f * kvy_d[S * N + idx];
				errorkvz += f * kvz_d[S * N + idx];
			}

			double errork = 0.0;
			errork += errorkx * errorkx * isc * isc;
			errork += errorky * errorky * isc * isc;
			errork += errorkz * errorkz * isc * isc;
			errork += errorkvx * errorkvx * isc * isc;
			errork += errorkvy * errorkvy * isc * isc;
			errork += errorkvz * errorkvz * isc * isc;

			errork = sqrt(errork / 6.0);	//6 is the number of dimensions

			s = pow( 1.0  / errork, ee);
			s = fmax(def_facmin, def_fac * s);
			s = fmin(def_facmax, s);

			// ****************************
			double2 snew2 = snew_d[idx];
			if(fabs(s * dti) < dtiMin){
				snew2.y = fmin(snew2.y, s);
			}
			if(snew2.y >= 1.0){
				snew2.x = s;
			}
			else{
				snew2.x = 1.5;
			}
			snew_d[idx] = snew2;
	//printf("snew %d %.20g %.20g %.20g\n", idx, snew_d[idx].y, snew_d[idx].x, errork);
			// ****************************
		}
	}

}


//every particle runs on an own thread block
//the force calculation is parallelized along the threads and reduced in registers
__global__ void stageStep2_kernel(unsigned long long int *id_d, double *m_d, double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *dx_d, double *dy_d, double *dz_d, double *dvx_d, double *dvy_d, double *dvz_d, double *xTable_d, double *yTable_d, double *zTable_d, double *A1_d, double *A2_d, double *A3_d, double2 *snew_d, double dt, double dti, double dtiMin, int N, int useHelio, int useGR, int useJ2, int useNonGrav, int useAdaptiveTimeSteps, double ee){

	int itx = threadIdx.x;
	int idx = blockIdx.x + Nperturbers;

	//shared memory contains only the perturbers
	__shared__ double x_s[Nperturbers];
	__shared__ double y_s[Nperturbers];
	__shared__ double z_s[Nperturbers];
	__shared__ double m_s[Nperturbers];
	__shared__ unsigned long long int id_s[Nperturbers];


	__shared__ double kx_s[RKFn];
	__shared__ double ky_s[RKFn];
	__shared__ double kz_s[RKFn];
	__shared__ double kvx_s[RKFn];
	__shared__ double kvy_s[RKFn];
	__shared__ double kvz_s[RKFn];

	__shared__ double xi_s[1];
	__shared__ double yi_s[1];
	__shared__ double zi_s[1];
	__shared__ double vxi_s[1];
	__shared__ double vyi_s[1];
	__shared__ double vzi_s[1];

	if(idx < N){
		double xi0 = x_d[idx];
		double yi0 = y_d[idx];
		double zi0 = z_d[idx];
		double vxi0 = vx_d[idx];
		double vyi0 = vy_d[idx];
		double vzi0 = vz_d[idx];

		double mi = m_d[idx];
		double A1i = A1_d[idx];
		double A2i = A2_d[idx];
		double A3i = A3_d[idx];
		unsigned long long int idi = id_d[idx];

		if(itx < Nperturbers){
			m_s[itx] = m_d[itx];
			id_s[itx] = id_d[itx];
//if(idx == 27) printf("%d %.20g %.20g\n", itx, x_s[itx], m_s[itx]);
		}

		for(int S = 0; S < RKFn; ++S){

			if(threadIdx.x < Nperturbers){
				int ii = itx * RKFn + S;
				x_s[itx] = xTable_d[ii];
				y_s[itx] = yTable_d[ii];
				z_s[itx] = zTable_d[ii];
//printf("%d %d %.20g\n", itx, S, x_s[itx]);
			}
			__syncthreads();

			if(itx == 0){

				// ***********************
				// update
				xi_s[0] = xi0;
				yi_s[0] = yi0;
				zi_s[0] = zi0;
				vxi_s[0] = vxi0;
				vyi_s[0] = vyi0;
				vzi_s[0] = vzi0;

				for(int s = 0; s < S; ++s){
					double dtaa = dt * a_c[S * RKFn + s];
					xi_s[0]  += dtaa * kx_s[s];
					yi_s[0]  += dtaa * ky_s[s];
					zi_s[0]  += dtaa * kz_s[s];
					vxi_s[0] += dtaa * kvx_s[s];
					vyi_s[0] += dtaa * kvy_s[s];
					vzi_s[0] += dtaa * kvz_s[s];
	//printf("x %d %d %.20g %.20g\n", S, s, aa, xi);
				}


				// end update
				// *****************************

				kx_s[S] = vxi_s[0];
				ky_s[S] = vyi_s[0];
				kz_s[S] = vzi_s[0];

				if(useHelio > 0){
					m_s[0] += mi;
				}
			}
			__syncthreads();
		
			double ax = 0.0;
			double ay = 0.0;
			double az = 0.0;

			if(useHelio == 0){
				if(itx < Nperturbers){
					if(idi != id_s[itx]){
						accP(m_s[itx], x_s[itx], y_s[itx], z_s[itx], xi_s[0], yi_s[0], zi_s[0], ax, ay, az);
					}
				}
			}
			else{
				if(itx < Nperturbers){
					if(idi != id_s[itx]){
						accP(m_s[itx], x_s[itx], y_s[itx], z_s[itx], xi_s[0], yi_s[0], zi_s[0], ax, ay, az);
//if(idx == 27) printf("Nij %d %d %.20g %.20g %.20g %g %g %g\n", idx, itx, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit, x_s[itx], xi_s[0], m_s[itx]);
					}
				}

//if(idx == 27) printf("N0 %d %d %.20g %.20g %.20g\n", idx, itx, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit);
				if(itx > 0 && itx < Nperturbers){
					if(idi != id_s[itx]){
						accP2(m_s[itx], x_s[itx], y_s[itx], z_s[itx], ax, ay, az);
	//printf("Npi %d %d %.20g %.20g %.20g %d\n", idx, j, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit, S);
					}
				}
			}
			__syncthreads();
			//reduce

			for(int i = 1; i < warpSize; i*=2){
			#if def_OldShuffle == 0
				ax += __shfl_xor_sync(0xffffffff, ax, i, warpSize);
				ay += __shfl_xor_sync(0xffffffff, ay, i, warpSize);
				az += __shfl_xor_sync(0xffffffff, az, i, warpSize);
			#else
				ax += __shfld_xor(ax, i);
				ay += __shfld_xor(ay, i);
				az += __shfld_xor(az, i);
			#endif
			//printf("s reduce  %d %d %d %.20g\n", blockIdx.x, i, threadIdx.x, s);
			}

			__syncthreads();
	//if(itx == 0){
	//printf("Np %d %.20g %.20g %.20g %d\n", idx, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit, S);
	//}

			if(itx == 0){

				if(useGR == 2){
					acchGR2(xi_s[0], yi_s[0], zi_s[0], vxi_s[0], vyi_s[0], vzi_s[0], ax, ay, az);
				}

				if(useNonGrav == 1){
					NonGrav(xi_s[0], yi_s[0], zi_s[0], vxi_s[0], vyi_s[0], vzi_s[0], ax, ay, az, A1i, A2i, A3i);
				}

				if(useJ2 == 1){
					J2(m_s[3], x_s[3], y_s[3], z_s[3], xi_s[0], yi_s[0], zi_s[0], ax, ay, az);
				}

				kvx_s[S] = ax;
				kvy_s[S] = ay;
				kvz_s[S] = az;

			}
			__syncthreads();
		}

		if(itx == 0){
			//update

			double dx = 0.0;
			double dy = 0.0;
			double dz = 0.0;
			double dvx = 0.0;
			double dvy = 0.0;
			double dvz = 0.0;
			
			for(int S = 0; S < RKFn; ++S){
				double dtb = dt * b_c[S];
				dx += dtb * kx_s[S];
				dy += dtb * ky_s[S];
				dz += dtb * kz_s[S];

				dvx += dtb * kvx_s[S];
				dvy += dtb * kvy_s[S];
				dvz += dtb * kvz_s[S];
			}
			
			dx_d[idx] = dx;
			dy_d[idx] = dy;
			dz_d[idx] = dz;
			dvx_d[idx] = dvx;
			dvy_d[idx] = dvy;
			dvz_d[idx] = dvz;

//if(idx == 27) printf("dx %d %.20g %.20g %.20g %.20g %.20g %.20g\n", idx, dx, dy, dz, dvx, dvy, dvz);

			if(useAdaptiveTimeSteps == 1){

				//compute error
				double ym = fabs(xi0);
				ym = fmax(ym, fabs(yi0));
				ym = fmax(ym, fabs(zi0));
				ym = fmax(ym, fabs(vxi0));
				ym = fmax(ym, fabs(vyi0));
				ym = fmax(ym, fabs(vzi0));

				ym = fmax(ym, fabs(xi0 + dx));
				ym = fmax(ym, fabs(yi0 + dy));
				ym = fmax(ym, fabs(zi0 + dz));
				ym = fmax(ym, fabs(vxi0 + dvx));
				ym = fmax(ym, fabs(vyi0 + dvy));
				ym = fmax(ym, fabs(vzi0 + dvz));


				double isc = 1.0 / (def_atol + ym * def_rtol);

				double s = 1.0e6;	//large number

				double errorkx = 0.0;
				double errorky = 0.0;
				double errorkz = 0.0;
				double errorkvx = 0.0;
				double errorkvy = 0.0;
				double errorkvz = 0.0;

				for(int S = 0; S < RKFn; ++S){
					//this is y1i - y^1i
					double f = (b_c[S] - bb_c[S]) * dt;
					errorkx += f * kx_s[S];
					errorky += f * ky_s[S];
					errorkz += f * kz_s[S];

					errorkvx += f * kvx_s[S];
					errorkvy += f * kvy_s[S];
					errorkvz += f * kvz_s[S];
				}
				
				double errork = 0.0;
				errork += errorkx * errorkx * isc * isc;
				errork += errorky * errorky * isc * isc;
				errork += errorkz * errorkz * isc * isc;
				errork += errorkvx * errorkvx * isc * isc;
				errork += errorkvy * errorkvy * isc * isc;
				errork += errorkvz * errorkvz * isc * isc;

				errork = sqrt(errork / 6.0);	//6 is the number of dimensions

				s = pow( 1.0  / errork, ee);
				s = fmax(def_facmin, def_fac * s);
				s = fmin(def_facmax, s);

				// ****************************
				double2 snew2 = snew_d[idx];
				if(fabs(s * dti) < dtiMin){
					snew2.y = fmin(snew2.y, s);
				}
				if(snew2.y >= 1.0){
					snew2.x = s;
				}
				else{
					snew2.x = 1.5;
				}
				snew_d[idx] = snew2;
//printf("snew %d %.20g %.20g %.20g\n", idx, snew_d[idx].y, snew_d[idx].x, errork);
				// ****************************
			}
		}
	}
}


__global__ void computeError_d1_kernel(double2 *snew_d, int N){

	int id = blockIdx.x * blockDim.x + threadIdx.x;

	double s = 1.0e6;	//large number

	extern __shared__ double se_s[];
	double *s_s = se_s;

	int lane = threadIdx.x % warpSize;
	int warp = threadIdx.x / warpSize;

	if(warp == 0){
		s_s[threadIdx.x] = 1.0e6;
	}
	__syncthreads();

	if(id >= Nperturbers && id < N){	
		s = snew_d[id].x;
	}
//printf("snew %d %g\n", id, s);

	__syncthreads();

	double sr = s;

	for(int i = 1; i < warpSize; i*=2){
#if def_OldShuffle == 0
		sr = __shfl_xor_sync(0xffffffff, s, i, warpSize);
#else
		sr = __shfld_xor(s, i);
#endif
		s = fmin(s, sr);
//printf("s reduce  %d %d %d %.20g\n", blockIdx.x, i, threadIdx.x, s);
	}
//printf("s reduce  %d %d %d %d %.20g\n", blockIdx.x, 0, id, threadIdx.x, s);

	__syncthreads();
	if(blockDim.x > warpSize){
		//reduce across warps
		if(lane == 0){
			s_s[warp] = s;
		}
		__syncthreads();
		//reduce previous warp results in the first warp
		if(warp == 0){
			s = s_s[threadIdx.x];
//printf("r reduce 1  %d %d %d %.20g %d %d\n", blockIdx.x, 0, threadIdx.x, s, int(blockDim.x), warpSize);
			for(int i = 1; i < warpSize; i*=2){
#if def_OldShuffle == 0
				sr = __shfl_xor_sync(0xffffffff, s, i, warpSize);
#else
				sr = __shfld_xor(s, i);
#endif
				s = fmin(s, sr);
//printf("r reduce 2  %d %d %d %.20g\n", blockIdx.x, i, threadIdx.x, s);
			}
		}
	}
	__syncthreads();
	if(threadIdx.x == 0){

		snew_d[blockIdx.x].x = s;
//printf("r reduce 3  %d %d %.20g\n", blockIdx.x, threadIdx.x, s);
	}

}

__global__ void computeError_d2_kernel(double2 *snew_d, const int N){

	int idy = threadIdx.x;
	int idx = blockIdx.x;

	double s = 1.0e6;

	extern __shared__ double se2_s[];
	double *s_s = se2_s;

	int lane = threadIdx.x % warpSize;
	int warp = threadIdx.x / warpSize;

	if(warp == 0){
		s_s[threadIdx.x] = 1.0e6;
	}


	if(idy < N){
		s = snew_d[idy].x;
	}
//printf("s reduce d2  %d %.20g\n", idy, s);

	__syncthreads();

	double sr = s;

	for(int i = 1; i < warpSize; i*=2){
#if def_OldShuffle == 0
		sr = __shfl_xor_sync(0xffffffff, s, i, warpSize);
#else
		sr = __shfld_xor(s, i);
#endif
		s = fmin(s, sr);
//if(idx == 0) printf("HC2bx %d %d %.20g\n", i, idy, a);
	}
	__syncthreads();

	if(blockDim.x > warpSize){
		//reduce across warps
		if(lane == 0){
			s_s[warp] = s;
		}
		__syncthreads();
		//reduce previous warp results in the first warp
		if(warp == 0){
			s = s_s[threadIdx.x];
			//if(idx == 0) printf("HC2cx %d %d %.20g %d %d\n", 0, idy, a, int(blockDim.x), warpSize);
			for(int i = 1; i < warpSize; i*=2){
#if def_OldShuffle == 0
				sr = __shfl_xor_sync(0xffffffff, s, i, warpSize);
#else
				sr = __shfld_xor(s, i);
#endif
				s = fmin(s, sr);
//printf("s reduce d2b %d %d %.20g\n", i, idy, s);
			}
		}
	}
	__syncthreads();
	if(threadIdx.x == 0){
		if(idx == 0){
//printf("s reduce d2c %d %.20g\n", idy, s);
			snew_d[0].x = s;
		}
	}
}


__host__ void Host::update(int i){

	x_h[i] += dx_h[i];
	y_h[i] += dy_h[i];
	z_h[i] += dz_h[i];

	vx_h[i] += dvx_h[i];
	vy_h[i] += dvy_h[i];
	vz_h[i] += dvz_h[i];

//printf("dx %d %.20g %.20g %.20g %.20g %.20g %.20g\n", i, dx_h[i], dy_h[i], dz_h[i], dvx_h[i], dvy_h[i], dvz_h[i]);
//printf("dx %d %.20g %.20g %.20g %.20g %.20g %.20g\n", i, x_h[i], y_h[i], z_h[i], vx_h[i], vy_h[i], vz_h[i]);
}

__global__ void update_kernel(double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *dx_d, double *dy_d, double *dz_d, double *dvx_d, double *dvy_d, double *dvz_d, int N){

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if(idx >= Nperturbers && idx < N){
		x_d[idx] += dx_d[idx];
		y_d[idx] += dy_d[idx];
		z_d[idx] += dz_d[idx];

		vx_d[idx] += dvx_d[idx];
		vy_d[idx] += dvy_d[idx];
		vz_d[idx] += dvz_d[idx];
	}
}


__host__ int Host::preIntegration(){
	if(useGPU > 0){
		//copy perturbers data to CPU, because its not yet here
		copyPerturbersCPU();
	}


	//Do pre-Integration to synchronize all initial conditions
	printf("Start pre-Integration up to time %.20g\n", outStart);

	double dtiMin = 1.0e-8; //This is a local copy of dtiMin
	for(int i = Nperturbers; i < N; ++i){
		time = jd_init_h[i];
		double snew = 1.0;
		dti = 10.0;
		dt = dti * dayUnit;
		int stop = 0;

		printf("preIntegration %d %.20g %.20g\n", i, time, outStart);

		if(time == outStart){
			//These particles are already at the correct epoch
			continue;
		}
		if(time > outStart){
			dti = -dti;
			dt = -dt;
		}

		for(long long int t = 1; t <= Nsteps; ++t){
		//for(long long int t = 1; t <= 10; ++t){
//printf("A %d %lld %.20g %.20g %.20g %.20g\n", i, id_h[i], time, x_h[i], y_h[i], z_h[i]);

			//interpolateTable(time);
			//stageStep1(dtiMin, i - Nperturbers, i + 1, snew);

			stageStep2(i, time, dtiMin, snew);

			if(useAdaptiveTimeSteps == 0){
				snew = 1.0;
			}

			if(snew >= 1.0){
				update(i);

				time += dti;

				dti *= snew;
				dt = dti * dayUnit;

				if(stop == 1){
					break;
				}
				else{
					dtmin_h[i] = fmin(fabs(dt / dayUnit), dtmin_h[i]);
				}
			}
			else{
				dti *= snew;
				dt = dti * dayUnit;
				stop = 0;
			}

//printf("B %d %lld %.20g %.20g %.20g %.20g %.20g %.20g\n", i, id_h[i], time, dti, snew, x_h[i], y_h[i], z_h[i]);


			if(dti > 0.0 && time + dti > outStart){
				dti = (outStart - time);
				dt = dti * dayUnit;
				stop = 1;
			//printf("Final time step %.20g\n", dti);
			}
			if(dti < 0.0 && time + dti < outStart){
				dti = (outStart - time);
				dt = dti * dayUnit;
				stop = 1;
			//printf("Final time step %.20g\n", dti);
			}

		} //en of t loop
//printf("C %d %lld %.20g %.20g %.20g %.20g %.20g %.20g\n", i, id_h[i], time, dti, snew, x_h[i], y_h[i], z_h[i]);

		if(fabs(time - outStart) > 1.0e-10){
			printf("Error in pre-integration of particle %d, start time not reached\n", i);
			return 0;
		}
		if(snew_h[i].y < 1.0){
			printf("Error, integration time step smaller than %g. Integration stopped\n", dtiMin);
			return 0;
		}

	}//end of i loop

	if(useGPU > 0){
		cudaMemcpy(x_d, x_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(y_d, y_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(z_d, z_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(vx_d, vx_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(vy_d, vy_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(vz_d, vz_h, N * sizeof(double), cudaMemcpyHostToDevice);
	}

	return 1;
}

__host__ void Host::IntegrationLoop(int S){
	//save initial values for j loop (individual time steps)
	double dti0 = dti;
	double dts0 = dts;
	long long int outI0 = outI;
	double time00 = time;

	int Nj = Nperturbers + 1;
	if(useIndividualTimeSteps == 1){
		Nj = N;
		//dtiMin = 1.0e-5;
	}

	//the loop is used for individual time step integrations
	for(int j = Nperturbers; j < Nj; ++j){
		dti = dti0;
		dts = dts0;
		dt = dti * dayUnit;
		outI = outI0;
		time = time00;

		unsigned long long int cOut = 0llu;	//counter for output
		int ci = 0;
		double dtiOld = dti;

printf("J %d time: %.20g, %.20g %.20g %.20g dti: %g, dts: %g, outI: %llu\n", j, time, outStart, time0, time1, dti, dts, outI);

		for(long long int t = 1; t <= Nsteps; ++t){
		//for(long long int t = 1; t < 5000; ++t){
		//for(long long int t = 1; t < 10; ++t){
			
			//cudaDeviceSynchronize();
			//printf("%lld %d | %d %d\n", t, 0, N, Nperturbers);	
			
			double snew = 1.0;
			
			if(useGPU == 0){

				if(useIndividualTimeSteps == 1){
					//single particle
					stageStep2(j, time, dtiMin[S], snew);
				}
				else{
					//interpolateTable(time);	
					//stageStep1(dtiMin[S], 0, N, snew);
				
					stageStep(time, dtiMin[S], snew);
				}
			}
			else{
				if(useGPU == 1){
					interpolateTable_kernel <<< dim3(Nperturbers, RKFn, 1), Ninterpolate >>> (xp_d, yp_d, zp_d, timep_d, timep0, dtimep, time, dti, xTable_d, yTable_d, zTable_d);
					if(N > 300){
						stageStep1_kernel <<< (N + 127) / 128, 128 >>> (id_d, m_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, dx_d, dy_d, dz_d, dvx_d, dvy_d, dvz_d, xTable_d, yTable_d, zTable_d, kx_d, ky_d, kz_d, kvx_d, kvy_d, kvz_d, A1_d, A2_d, A3_d, snew_d, dt, dti, dtiMin[S], N, useHelio, useGR, useJ2, useNonGrav, useAdaptiveTimeSteps, ee);
					}
					else{
						stageStep2_kernel <<< (N - Nperturbers), 32 >>> (id_d, m_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, dx_d, dy_d, dz_d, dvx_d, dvy_d, dvz_d, xTable_d, yTable_d, zTable_d, A1_d, A2_d, A3_d, snew_d, dt, dti, dtiMin[S], N, useHelio, useGR, useJ2, useNonGrav, useAdaptiveTimeSteps, ee);
						//stageStep3_kernel < 2 > <<< ((N - Nperturbers) + 1) / 2, dim3(32, 2, 1) >>> (id_d, m_d, x_d, y_d, z_d, vx_d, vy_d, vz_d, dx_d, dy_d, dz_d, dvx_d, dvy_d, dvz_d, xTable_d, yTable_d, zTable_d, A1_d, A2_d, A3_d, snew_d, dt, dti, dtiMin[S], N, useHelio, useGR, useJ2, useNonGrav, useAdaptiveTimeSteps, ee);
					}
					
				}

				if(useAdaptiveTimeSteps == 1){			
					int nct = 512;
					int ncb = min((N + nct - 1) / nct, 1024);
					computeError_d1_kernel <<< ncb, nct, WarpSize * sizeof(double)  >>> (snew_d, N);
					if(ncb > 1){
						computeError_d2_kernel <<< 1, ((ncb + WarpSize - 1) / WarpSize) * WarpSize, WarpSize * sizeof(double)  >>> (snew_d, ncb);
					}
					cudaMemcpy(snew_h, snew_d, sizeof(double2), cudaMemcpyDeviceToHost);
					cudaDeviceSynchronize();
					snew = snew_h[0].x;
				}
			}
			//Only used for testing
			//if(useGPU > 0){
				//cudaMemcpy(x_h, x_d, N * sizeof(double), cudaMemcpyDeviceToHost);
				//cudaMemcpy(y_h, y_d, N * sizeof(double), cudaMemcpyDeviceToHost);
				//cudaMemcpy(z_h, z_d, N * sizeof(double), cudaMemcpyDeviceToHost);
			//}
			printf("%.20g %llu dt: %.20g %.g %g %d\n", time, cOut, dti, dts, snew, S);
			fprintf(dtfile, "%.20g %llu dt: %.20g %.g %g %d\n", time, cOut, dti, dts, snew, S);
			
			if(useAdaptiveTimeSteps == 0){
				snew = 1.0;
			}
			
			if(snew >= 1.0){
				ci = (fabs(dti) + 0.5 * dts) / dts;
				printf("---- accept step %llu %llu %llu %.20g %.20g %g----\n", cOut, cOut + ci, outI, time, time + dti, dti);		
				if(useGPU == 0){
					if(useIndividualTimeSteps == 1){
						update(j);
					}
					else{
						for(int i = Nperturbers; i < N; ++i){
							update(i);
						}
					}
				}
				else{
					update_kernel <<< (N + 127) / 128, 128 >>> (x_d, y_d, z_d, vx_d, vy_d, vz_d, dx_d, dy_d, dz_d, dvx_d, dvy_d, dvz_d, N);	
				}
				cOut += ci;		
				
				time += dti;
				
				if(cOut >= outI && ((dt > 0 && time >= outStart) || (dt < 0 && time <= outStart))){
				//if(t % 10 == 0){
					if(useGPU > 0){
						cudaMemcpy(snew_h, snew_d, N * sizeof(double2), cudaMemcpyDeviceToHost);
						cudaMemcpy(x_h, x_d, N * sizeof(double), cudaMemcpyDeviceToHost);
						cudaMemcpy(y_h, y_d, N * sizeof(double), cudaMemcpyDeviceToHost);
						cudaMemcpy(z_h, z_d, N * sizeof(double), cudaMemcpyDeviceToHost);
						cudaMemcpy(vx_h, vx_d, N * sizeof(double), cudaMemcpyDeviceToHost);
						cudaMemcpy(vy_h, vy_d, N * sizeof(double), cudaMemcpyDeviceToHost);
						cudaMemcpy(vz_h, vz_d, N * sizeof(double), cudaMemcpyDeviceToHost);

						cudaMemcpy(xTable_h, xTable_d, Nperturbers * RKFn * sizeof(double), cudaMemcpyDeviceToHost);
						cudaMemcpy(yTable_h, yTable_d, Nperturbers * RKFn * sizeof(double), cudaMemcpyDeviceToHost);
						cudaMemcpy(zTable_h, zTable_d, Nperturbers * RKFn * sizeof(double), cudaMemcpyDeviceToHost);


					}
					output(t, time, S);

					dti = dtiOld;
					dt = dti * dayUnit;
					snew = 1.0;
					cOut = 0;
					outI = (outInterval + 0.5 * dts) / dts; //needed only at the first time

//reduceCall();
				}
			}
			else{
				printf(" ---- repeat step %llu %llu %.20g ----\n", cOut, outI, time);		
				
			}
			
			
			dti *= snew;
			
			if(fabs(dtiOld) >= 65.0 * fabs(dts) && cOut % 10 == 0 && snew >= 1.0){
				//synchronize with coarser time steps
				dts *= 10;
				outI /= 10;
				cOut /= 10;
				dti = 5.0 * dts;
				if(dt < 0) dti = - dti;
printf("increase time step C %g %g\n", dti, dts);
			}
			
			//round dti to dts intervals
			int dtt = dti / dts;
printf("dta %.20g %d %g %llu\n", dti, dtt, dts, cOut);
			dti = dtt * dts;
			//if(dti < dts) dti = dts;
			
			
printf("dtb %.20g %d %g %llu\n", dti, dtt, dts, cOut);
			ci = (fabs(dti) + 0.5 * dts) / dts;
			
			if(dti > 0 && dti < dtiMin[S]) dti = dtiMin[S];
			if(dti < 0 && -dti < dtiMin[S]) dti = -dtiMin[S];
			
			
			dtiOld = dti;
			
			dt = dti * dayUnit;
			//printf("%llu %llu %.20g, %.20g %.20g\n", cOut + ci, outI, time, time + dti, outStart);
			
/*			if(cOut + ci > outI && ((dt > 0 && time + dti >= outStart) || (dt < 0 && time + dti <= outStart))){
				dti = (outI - cOut) * dts;
				if(dt < 0) dti = -dti;
				
				dt = dti * dayUnit;
printf("   correct %.20g %.20g %.20g %.20g %llu %llu\n", time, time + dti, dti, dtiOld, cOut, outI);
printf("dtc %.20g %g\n", dti, dts);
			}
			
			
			else{
*/{				
				if(fabs(dti) >= 65.0 * dts){
					//synchronize with coarser time steps
printf("increase time step A %g %g\n", dti, dts);
					dti = (((cOut + ci) / 10) * 10 - cOut) * dts;
					if(dt < 0) dti = -dti;
					dt = dti * dayUnit;
					//dts *= 10;
					//outI /= 10;
					//cOut /= 10;
printf("increase time step B %g %g\n", dti, dts);
				}
				if(fabs(dti) <= 2.0 * dts){
printf("refine time steps %g %g\n", dti, dts);
					//refine interpolation points
					dts *= 0.1;
					outI *= 10;
					cOut *= 10;
				}
			}
			
			
			//add condition for backward integration
			if((dt > 0 && time >= time1) || (dt < 0 && time <= time1)){
				printf("Reached the end\n");
				break;
			}
			
		}	// end of time step loop
	}// end of j loop

}

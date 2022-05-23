#define dayUnit 0.01720209895

//check RKN12(10)17M

void setRKF45(double *a_h, double *b_h, double *bb_h, double *c_h){

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

}

void setDP54(double *a_h, double *b_h, double *bb_h, double *c_h){

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

}


void setRKF78(double *a_h, double *b_h, double *bb_h, double *c_h){


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


__host__ void update2(double *xt, double *yt, double *zt, double *vxt, double *vyt, double *vzt, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *kx, double *ky, double *kz, double *kvx, double *kvy, double *kvz, int i, int N, double dt, int S, int RKFn, double *a){

        xt[i]  = x[i];
        yt[i]  = y[i];
        zt[i]  = z[i];
        vxt[i] = vx[i];
        vyt[i] = vy[i];
        vzt[i] = vz[i];

        for(int s = 0; s < S; ++s){
                double dtaa = dt * a[S * RKFn + s];
                xt[i]  += dtaa * kx[i + s * N];
                yt[i]  += dtaa * ky[i + s * N];
                zt[i]  += dtaa * kz[i + s * N];
                vxt[i] += dtaa * kvx[i + s * N];
                vyt[i] += dtaa * kvy[i + s * N];
                vzt[i] += dtaa * kvz[i + s * N];
//printf("update 2 %d %d %g %g %g %g %g %g\n", S, i, xt[i], yt[i], zt[i], a[S * RKFn + s], kx[i + s * N], dt);

        }
//printf("update 2 %d %g %g %g %g %g %g\n", i, xt[i], yt[i], zt[i], vxt[i], vyt[i], vzt[i]);
}


__host__ void stageStep(unsigned long long int *id_h, double *m_h, double *xt, double *yt, double *zt, double *vxt, double *vyt, double *vzt, double *kx_h, double *ky_h, double *kz_h, double *kvx_h, double *kvy_h, double *kvz_h, double *A1_h, double *A2_h, double *A3_h, int S, int i, int Nperturbers, int N, int useHelio, int useGR, int useJ2, int useNonGrav){

	kx_h[i + S * N] = vxt[i];
	ky_h[i + S * N] = vyt[i];
	kz_h[i + S * N] = vzt[i];

	double ax = 0.0;
	double ay = 0.0;
	double az = 0.0;

	if(useHelio == 0){
		for(int j = Nperturbers-1; j >= 0; --j){
			if(id_h[i] != id_h[j]){ 
				accP(m_h[j], xt[j], yt[j], zt[j], xt[i], yt[i], zt[i], ax, ay, az);
			}
		}
	}
	else{
		for(int j = Nperturbers-1; j >= 1; --j){
			if(id_h[i] != id_h[j]){
				accP(m_h[j], xt[j], yt[j], zt[j], xt[i], yt[i], zt[i], ax, ay, az);
//if(i == 27) printf("Nij %d %d %llu %llu %.20g %.20g %.20g\n", i, j, id_h[i], id_h[j], ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit);
			}
		}
		accS(m_h[0] + m_h[i], xt[i], yt[i], zt[i], ax, ay, az);
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
		acchGR2(xt[i], yt[i], zt[i], vxt[i], vyt[i], vzt[i], ax, ay, az);
	}

	if(useNonGrav == 1){
		NonGrav(xt[i], yt[i], zt[i], vxt[i], vyt[i], vzt[i], ax, ay, az, A1_h[i], A2_h[i], A3_h[i]);
	}

	if(useJ2 == 1){
		J2(m_h[3], xt[3], yt[3], zt[3], xt[i], yt[i], zt[i], ax, ay, az);
	}

	kvx_h[i + S * N] = ax;
	kvy_h[i + S * N] = ay;
	kvz_h[i + S * N] = az;
}

__host__ void stageStep1(unsigned long long int *id_h, double *m_h, double *x_h, double *y_h, double *z_h, double *vx_h, double *vy_h, double *vz_h, double *dx_h, double *dy_h, double *dz_h, double *dvx_h, double *dvy_h, double *dvz_h, double *xTable_h, double *yTable_h, double *zTable_h, double *A1_h, double *A2_h, double *A3_h, double2 *snew_h, double *a_h, double *b_h, double *bb_h, double dt, double dti, double dtiMin, int RKFn, int Nperturbers, int N, int iStart, int useHelio, int useGR, int useJ2, int useNonGrav, double ee, double &snew){


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

		double s = 1.0e6;       //large number

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

		errork = sqrt(errork / 6.0);    //6 is the number of dimensions

		s = pow( 1.0  / errork, ee);
//printf("error %d %g %g %g\n", i, errork, s, errorkx);
		s = fmax(def_facmin, def_fac * s);
		s = fmin(def_facmax, s);

		// ****************************
		if(snew_h[i].y >= 1.0){
			snew_h[i].x = s;
		}
		else{
			snew_h[i].x = 1.5;
		}
		if(fabs(s * dti) < dtiMin){
			snew_h[i].y = fmin(snew_h[i].y, s);
		}

		snew = fmin(snew, snew_h[i].x);
//printf("snew %d %.20g %.20g %.20g %.20g %g\n", i, snew_h[i].y, snew_h[i].x, errork, snew, s);
		// ****************************

	}

}


template < int Nperturbers >
__global__ void stageStep1_kernel(unsigned long long int *id_d, double *m_d, double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *dx_d, double *dy_d, double *dz_d, double *dvx_d, double *dvy_d, double *dvz_d, double *xTable_d, double *yTable_d, double *zTable_d, double *kx_d, double *ky_d, double *kz_d, double *kvx_d, double *kvy_d, double *kvz_d, double *A1_d, double *A2_d, double *A3_d, double2 *snew_d, double dt, double dti, double dtiMin, int RKFn, int N, int useHelio, int useGR, int useJ2, int useNonGrav, double ee){

	int itx = threadIdx.x;
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	//shared memory contains only the perturbers
	//the particle idx is stored in registers
	__shared__ double x_s[Nperturbers];
	__shared__ double y_s[Nperturbers];
	__shared__ double z_s[Nperturbers];
	__shared__ double m_s[Nperturbers];
	__shared__ unsigned long long int id_s[Nperturbers];

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

//printf("N0 %d %.20g %.20g %.20g\n", idx, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit);
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

		double s = 1.0e6;       //large number

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

		errork = sqrt(errork / 6.0);    //6 is the number of dimensions

		s = pow( 1.0  / errork, ee);
		s = fmax(def_facmin, def_fac * s);
		s = fmin(def_facmax, s);

		// ****************************
		if(snew_d[idx].y >= 1.0){
			snew_d[idx].x = s;
		}
		else{
			snew_d[idx].x = 1.5;
		}
		if(fabs(s * dti) < dtiMin){
			snew_d[idx].y = fmin(snew_d[idx].y, s);
		}
//printf("snew %d %.20g %.20g %.20g\n", idx, snew_d[idx].y, snew_d[idx].x, errork);
		// ****************************
	}


}


//every particle runs on an own thread block
//the force calculation is parallelized along the threads and reduced in registers
template < int Nperturbers, int RKFn >
__global__ void stageStep2_kernel(unsigned long long int *id_d, double *m_d, double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *dx_d, double *dy_d, double *dz_d, double *dvx_d, double *dvy_d, double *dvz_d, double *xTable_d, double *yTable_d, double *zTable_d, double *A1_d, double *A2_d, double *A3_d, double2 *snew_d, double dt, double dti, double dtiMin, int N, int useHelio, int useGR, int useJ2, int useNonGrav, double ee){

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
//printf("%d %d %.20g\n", itx, S, x_s[itx]);
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
	//printf("Nij %d %d %.20g %.20g %.20g\n", idx, j, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit);
					}
				}

	//printf("N0 %d %.20g %.20g %.20g\n", idx, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit);
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

			double s = 1.0e6;       //large number

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
			if(snew2.y >= 1.0){
				snew2.x = s;
			}
			else{
				snew2.x = 1.5;
			}
			if(fabs(s * dti) < dtiMin){
				snew2.y = fmin(snew2.y, s);
			}
			snew_d[idx] = snew2;
//printf("snew %d %.20g %.20g %.20g\n", idx, snew_d[idx].y, snew_d[idx].x, errork);
			// ****************************

		}
	}
}


//compute adaptive time step for only 1 particle
void computeError1(double2 *snew_h, double *x_h, double *y_h, double *z_h, double *vx_h, double *vy_h, double *vz_h,double *kx_h, double *ky_h, double *kz_h, double *kvx_h, double *kvy_h, double *kvz_h, double *dx_h, double *dy_h, double *dz_h, double *dvx_h, double *dvy_h, double *dvz_h, double *b_h, double *bb_h, int RKFn, int i, int N, double &snew, double dt, double dti, double dtiMin, double ee){


	//update
	double dx = 0.0;
	double dy = 0.0;
	double dz = 0.0;
	double dvx = 0.0;
	double dvy = 0.0;
	double dvz = 0.0;

	for(int S = 0; S < RKFn; ++S){
		double dtb = dt * b_h[S];
		dx += dtb * kx_h[S * N + i];
		dy += dtb * ky_h[S * N + i];
		dz += dtb * kz_h[S * N + i];

		dvx += dtb * kvx_h[S * N + i];
		dvy += dtb * kvy_h[S * N + i];
		dvz += dtb * kvz_h[S * N + i];
	}

	dx_h[i] = dx;
	dy_h[i] = dy;
	dz_h[i] = dz;
	dvx_h[i] = dvx;
	dvy_h[i] = dvy;
	dvz_h[i] = dvz;

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

	double s = 1.0e6;       //large number
	snew = 1.0e6;   //large number

	//error estimation
	double errorkx = 0.0;
	double errorky = 0.0;
	double errorkz = 0.0;
	double errorkvx = 0.0;
	double errorkvy = 0.0;
	double errorkvz = 0.0;

	for(int S = 0; S < RKFn; ++S){
		double f = (b_h[S] - bb_h[S]) * dt;
		errorkx += f * kx_h[S * N + i];
		errorky += f * ky_h[S * N + i];
		errorkz += f * kz_h[S * N + i];

		errorkvx += f * kvx_h[S * N + i];
		errorkvy += f * kvy_h[S * N + i];
		errorkvz += f * kvz_h[S * N + i];
//printf("error %d %d %g %g\n", i, S, errorkx, kx_h[S * N + i]);
	}

	double errork = 0.0;
	errork += errorkx * errorkx * isc * isc;
	errork += errorky * errorky * isc * isc;
	errork += errorkz * errorkz * isc * isc;
	errork += errorkvx * errorkvx * isc * isc;
	errork += errorkvy * errorkvy * isc * isc;
	errork += errorkvz * errorkvz * isc * isc;

	errork = sqrt(errork / 6.0);    //6 is the number of dimensions

	s = pow( 1.0  / errork, ee);
//printf("error %d %g %g %g\n", i, errork, s, errorkx);
	s = fmax(def_facmin, def_fac * s);
	s = fmin(def_facmax, s);

	// ****************************
	if(snew_h[i].y >= 1.0){
		snew_h[i].x = s;
	}
	else{
		snew_h[i].x = 1.5;
	}
	if(fabs(s * dti) < dtiMin){
		snew_h[i].y = fmin(snew_h[i].y, s);
	}

	snew = fmin(snew, snew_h[i].x);
//printf("snew %d %.20g %.20g %.20g %.20g %g\n", i, snew_h[i].y, snew_h[i].x, errork, snew, s);

}


void computeError(double2 *snew_h, double *x_h, double *y_h, double *z_h, double *vx_h, double *vy_h, double *vz_h,double *kx_h, double *ky_h, double *kz_h, double *kvx_h, double *kvy_h, double *kvz_h, double *dx_h, double *dy_h, double *dz_h, double *dvx_h, double *dvy_h, double *dvz_h, double *b_h, double *bb_h, int RKFn, int Nperturbers, int N, double &snew, double dt, double dti, double dtiMin, double ee){


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
			dx += dtb * kx_h[S * N + i];
			dy += dtb * ky_h[S * N + i];
			dz += dtb * kz_h[S * N + i];

			dvx += dtb * kvx_h[S * N + i];
			dvy += dtb * kvy_h[S * N + i];
			dvz += dtb * kvz_h[S * N + i];
		}

		dx_h[i] = dx;
		dy_h[i] = dy;
		dz_h[i] = dz;
		dvx_h[i] = dvx;
		dvy_h[i] = dvy;
		dvz_h[i] = dvz;

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

		double s = 1.0e6;       //large number

		//error estimation
		double errorkx = 0.0;
		double errorky = 0.0;
		double errorkz = 0.0;
		double errorkvx = 0.0;
		double errorkvy = 0.0;
		double errorkvz = 0.0;

		for(int S = 0; S < RKFn; ++S){
			double f = (b_h[S] - bb_h[S]) * dt;
			errorkx += f * kx_h[S * N + i];
			errorky += f * ky_h[S * N + i];
			errorkz += f * kz_h[S * N + i];

			errorkvx += f * kvx_h[S * N + i];
			errorkvy += f * kvy_h[S * N + i];
			errorkvz += f * kvz_h[S * N + i];
	//printf("error %d %d %g %g\n", i, S, errorkx, kx_h[S * N + i]);
		}

		double errork = 0.0;
		errork += errorkx * errorkx * isc * isc;
		errork += errorky * errorky * isc * isc;
		errork += errorkz * errorkz * isc * isc;
		errork += errorkvx * errorkvx * isc * isc;
		errork += errorkvy * errorkvy * isc * isc;
		errork += errorkvz * errorkvz * isc * isc;

		errork = sqrt(errork / 6.0);    //6 is the number of dimensions

		s = pow( 1.0  / errork, ee);
	//printf("error %d %g %g %g\n", i, errork, s, errorkx);
		s = fmax(def_facmin, def_fac * s);
		s = fmin(def_facmax, s);

		// ****************************
		if(snew_h[i].y >= 1.0){
			snew_h[i].x = s;
		}
		else{
			snew_h[i].x = 1.5;
		}
		if(fabs(s * dti) < dtiMin){
			snew_h[i].y = fmin(snew_h[i].y, s);
		}

		snew = fmin(snew, snew_h[i].x);
//printf("snew %d %.20g %.20g %.20g %.20g %g\n", i, snew_h[i].y, snew_h[i].x, errork, snew, s);
	}
}

__global__ void computeError_d1_kernel(double2 *snew_d, int Nperturbers, int N){

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


void update(double *x_h, double *y_h, double *z_h, double *vx_h, double *vy_h, double *vz_h, double *dx_h, double *dy_h, double *dz_h, double *dvx_h, double *dvy_h, double *dvz_h, int i){

	x_h[i] += dx_h[i];
	y_h[i] += dy_h[i];
	z_h[i] += dz_h[i];

	vx_h[i] += dvx_h[i];
	vy_h[i] += dvy_h[i];
	vz_h[i] += dvz_h[i];

//printf("dx %d %.20g %.20g %.20g %.20g %.20g %.20g\n", i, dx_h[i], dy_h[i], dz_h[i], dvx_h[i], dvy_h[i], dvz_h[i]);
//printf("dx %d %.20g %.20g %.20g %.20g %.20g %.20g\n", i, x_h[i], y_h[i], z_h[i], vx_h[i], vy_h[i], vz_h[i]);
}

__global__ void update_kernel(double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *dx_d, double *dy_d, double *dz_d, double *dvx_d, double *dvy_d, double *dvz_d, int N, int Nperturbers){

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



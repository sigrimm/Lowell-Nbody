//every particle runs on an own thread block
//the force calculation is parallelized along the threads and reduced in registers
template < int ny >
__global__ void stageStep3_kernel(unsigned long long int *id_d, double *m_d, double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *dx_d, double *dy_d, double *dz_d, double *dvx_d, double *dvy_d, double *dvz_d, double *xTable_d, double *yTable_d, double *zTable_d, double *A1_d, double *A2_d, double *A3_d, double2 *snew_d, double dt, double dti, double dtiMin, int N, int useHelio, int useGR, int useJ2, int useNonGrav, int useAdaptiveTimeSteps, double ee){

	int itx = threadIdx.x;
	int ity = threadIdx.y;
	int idx = blockIdx.x * blockDim.y + ity + Nperturbers;

	//shared memory contains only the perturbers
	__shared__ double x_s[Nperturbers];
	__shared__ double y_s[Nperturbers];
	__shared__ double z_s[Nperturbers];
	__shared__ double m_s[Nperturbers];
	__shared__ unsigned long long int id_s[Nperturbers];


	__shared__ double kx_s[RKFn][ny];
	__shared__ double ky_s[RKFn][ny];
	__shared__ double kz_s[RKFn][ny];
	__shared__ double kvx_s[RKFn][ny];
	__shared__ double kvy_s[RKFn][ny];
	__shared__ double kvz_s[RKFn][ny];

	if(idx < N){
		double xi0 = x_d[idx];
		double yi0 = y_d[idx];
		double zi0 = z_d[idx];
		double vxi0 = vx_d[idx];
		double vyi0 = vy_d[idx];
		double vzi0 = vz_d[idx];
	
		double xi;
		double yi;
		double zi;
		double vxi;
		double vyi;
		double vzi;

		double mi = m_d[idx];
		double A1i = A1_d[idx];
		double A2i = A2_d[idx];
		double A3i = A3_d[idx];
		unsigned long long int idi = id_d[idx];

		if(itx < Nperturbers && ity == 0){
			m_s[itx] = m_d[itx];
			id_s[itx] = id_d[itx];
//printf("%d %d %.20g\n", itx, S, x_s[itx]);
		}
		__syncthreads();

		double mH = 0.0;
		for(int S = 0; S < RKFn; ++S){

			if(threadIdx.x < Nperturbers && ity == 0){
				int ii = itx * RKFn + S;
				x_s[itx] = xTable_d[ii];
				y_s[itx] = yTable_d[ii];
				z_s[itx] = zTable_d[ii];
//printf("%d %d %.20g\n", itx, S, x_s[itx]);
			}
			__syncthreads();


			// ***********************
			// update
			xi = xi0;
			yi = yi0;
			zi = zi0;
			vxi = vxi0;
			vyi = vyi0;
			vzi = vzi0;

			for(int s = 0; s < S; ++s){
				double dtaa = dt * a_c[S * RKFn + s];
				xi  += dtaa * kx_s[s][ity];
				yi  += dtaa * ky_s[s][ity];
				zi  += dtaa * kz_s[s][ity];
				vxi += dtaa * kvx_s[s][ity];
				vyi += dtaa * kvy_s[s][ity];
				vzi += dtaa * kvz_s[s][ity];
//printf("x %d %d %.20g %.20g\n", S, s, aa, xi);
			}


			// end update
			// *****************************

			if(itx == 0){
				kx_s[S][ity] = vxi;
				ky_s[S][ity] = vyi;
				kz_s[S][ity] = vzi;
	
				if(useHelio > 0){
					mH = mi;
					//m_s[0] += mi;
				}
			}
			__syncthreads();
		
			double ax = 0.0;
			double ay = 0.0;
			double az = 0.0;

			if(useHelio == 0){
				if(itx < Nperturbers){
					if(idi != id_s[itx]){
						accP(m_s[itx], x_s[itx], y_s[itx], z_s[itx], xi, yi, zi, ax, ay, az);
					}
				}
			}
			else{
				if(itx < Nperturbers){
					if(idi != id_s[itx]){
						accP(m_s[itx] + mH, x_s[itx], y_s[itx], z_s[itx], xi, yi, zi, ax, ay, az);
//printf("Nij %d %d %.20g %.20g %.20g\n", idx, j, ax * dayUnit * dayUnit, ay * dayUnit * dayUnit, az * dayUnit * dayUnit);
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
					acchGR2(xi, yi, zi, vxi, vyi, vzi, ax, ay, az);
				}

				if(useNonGrav == 1){
					NonGrav(xi, yi, zi, vxi, vyi, vzi, ax, ay, az, A1i, A2i, A3i);
				}

				if(useJ2 == 1){
					J2(m_s[3], x_s[3], y_s[3], z_s[3], xi, yi, zi, ax, ay, az);
				}

				kvx_s[S][ity] = ax;
				kvy_s[S][ity] = ay;
				kvz_s[S][ity] = az;

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
				dx += dtb * kx_s[S][ity];
				dy += dtb * ky_s[S][ity];
				dz += dtb * kz_s[S][ity];

				dvx += dtb * kvx_s[S][ity];
				dvy += dtb * kvy_s[S][ity];
				dvz += dtb * kvz_s[S][ity];
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
					errorkx += f * kx_s[S][ity];
					errorky += f * ky_s[S][ity];
					errorkz += f * kz_s[S][ity];

					errorkvx += f * kvx_s[S][ity];
					errorkvy += f * kvy_s[S][ity];
					errorkvz += f * kvz_s[S][ity];
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
}


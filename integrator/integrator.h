#include "asteroid.h"

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

	//Kick
	// ----------------------------------------------------------------------------
	for(int i = Nperturbers; i < N; ++i){
		ax[i] = 0.0;
		ay[i] = 0.0;
		az[i] = 0.0;
	}
	// ----------------------------------------------------------------------------
	//Update the Chebyshev coefficients if necessary
	update_Chebyshev(time);
	// ----------------------------------------------------------------------------
	update_perturbers(time);

	NonGrav();
	GR();
	J2();
	Gravity();

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


inline int asteroid::loop(){


        for(int tt = 0; tt < 2000; ++tt){

                double timett1 = timeStart + dts * (tt + 1) * 10.0;

//printf("integrate %.20g %.20g\n", timeStart + dts * tt * 10.0, timett1);


                for(int ttt = 0; ttt < 10 * 100 + 1; ++ttt){

                        //refine last time step of interval to match output time
                        if(time + dt < timett1){
                                dt = (timett1 - time);
                                stop = 1;
//printf("refine %.20g %.20g\n", timett1 - A.time, (timett1 - A.time) / 2.0);

                        }

                        leapfrog_step();

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
                }//end of ttt loop
                for(int p = Nperturbers; p < N; ++p){
printf("%.20g %d %.20g %.20g %.20g\n", time_reference + time, p, x[p], y[p], z[p]);

                }

        }//end of tt loop
	return 1;
}

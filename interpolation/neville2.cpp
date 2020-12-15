#include <stdlib.h>
#include <stdio.h>
#include <math.h>



int main(){
		
	
	const int NN = 10;

	double timet[21 * NN];
	double xt[21 * NN];
	double yt[21 * NN];
	double zt[21 * NN];

	for(int tt = 0; tt < 39200; ++tt){
		double time = 2433290.5 + tt* 0.01;
		//double time = 2433290.5;
		double dtime = 1.0;	//interval between stored time steps

		FILE *XVfile;
		XVfile = fopen("xv.dat", "r");


		int countNodes = 0;
		for(int t = 0; t < 10000; ++t){
			int er;
//printf("CountNodes %d\n", countNodes);
			for(int i = 0; i < 21; ++i){
				double skip;
				double timett;
				int id;
				er = fscanf(XVfile, "%lf %d %lf %lf", &timett, &id, &skip, &skip);
				fscanf(XVfile, "%lf %lf %lf", &xt[id * NN + countNodes], &yt[id * NN + countNodes], &zt[id * NN + countNodes]);
				fscanf(XVfile, "%lf %lf %lf", &skip, &skip, &skip); //v
				fscanf(XVfile, "%lf %lf %lf", &skip, &skip, &skip);
				fscanf(XVfile, "%lf %lf %lf", &skip, &skip, &skip);
				fscanf(XVfile, "%lf %lf %lf", &skip, &skip, &skip);
				er = fscanf(XVfile, "%lf %lf", &skip, &skip);
				if(er < 0) break;
				timet[id * NN + countNodes] = timett;
//printf("%.20g %d %.20g %.20g %.20g %d\n", timet[id * NN + countNodes], id, xt[id * NN + countNodes], yt[id * NN + countNodes], zt[id * NN + countNodes], id * NN + countNodes);


				
				if(i == 0 && t == 0 && timet[id * NN + countNodes] > time - (NN/2 - 1) * dtime){
					printf("Error, time too small, not enough data before time\n");
					return 0;
				}
				if(i == 20 && timet[id * NN + countNodes] > time - NN/2 * dtime){
					++countNodes;
				}

			}
			if(er < 0) break;
			if(countNodes >= NN){
				break;
			}
		}
		fclose(XVfile);
		if(countNodes < NN){
			printf("Error, time too large, not enough data after time\n");
			return 0;
		}



		for(int p = 0; p < 21; ++p){
			double r3[3];

			for(int k = 0; k < 3; ++k){
				
				double P[NN][NN];
				double tn[NN];

				for(int i = 0; i < NN; ++i){
					if(k == 0) P[0][i] = xt[p * NN + i];
					if(k == 1) P[0][i] = yt[p * NN + i];
					if(k == 2) P[0][i] = zt[p * NN + i];
					tn[i] = timet[p * NN + i];
					

//printf("%d %.20g %.20g\n", i, tn[i], P[0][i]);
				}

				for(int j = 1; j < NN; ++j){
//printf("****\n");
					for(int i = 0; i < NN - j; ++i){
						P[j][i] = ((time - tn[i+j]) * P[j-1][i] + (tn[i] - time) * P[j-1][i+1]) / (tn[i] - tn[i+j]);
//printf("%d %d %g %g %g %g %.20g\n", i, i+j, tn[i], tn[i+j], P[j-1][i], P[j-1][i+1], P[j][i]);
				
					}
				}
				r3[k] = P[NN-1][0];

			}

			printf("%.20g %d %.20g %.20g %.20g\n", time, p, r3[0], r3[1], r3[2]);	

		}

	}
	return 0;
}

#include <stdlib.h>
#include <stdio.h>
#include <math.h>


double fx(double x){

	double f = 1.0 / (1.0 + x * x);

}


int main(){


		
	//for(double x = - 5.0; x < 5.0; x+=0.01){
	//	printf("%g %g\n", x, fx(x));
	//
	//}


		
	double x = 0.0;
	const int NN = 6;

	double xn[NN];
	double yn[NN];
	int i = 0;
	for(double x = -5; x <= 5.0; x += 2.0){
		xn[i] = x;
		yn[i] = fx(x);
//printf("%d %g %g\n", i, x, yn[i]);
		++i;

	}
	/*	
	double x = 8.4;
	const int NN = 4;

	double xn[NN] = {8.1, 8.3, 8.6, 8.7};
	double yn[NN] = {16.9446, 17.56492, 18.50515, 18.82091};
	*/

	double P[NN][NN];

	for(int i = 0; i < NN; ++i){
		P[0][i] = yn[i];

//printf("%d %g\n", i, P[0][i]);
	}

	for(int i = 0; i < NN - 1; ++i){
		P[1][i] = ((x - xn[i+1]) * P[0][i] + (xn[i] - x) * P[0][i+1]) / (xn[i] - xn[i+1]);
printf("%d %d %g %g %g %g %g\n", i, i+1, xn[i], xn[i+1], P[0][i], P[0][i+1], P[1][i]);
	
	}

	for(int j = 2; j < 6; ++j){
printf("****\n");
		for(int i = 0; i < NN - j; ++i){
			P[j][i] = ((x - xn[i+j]) * P[j-1][i] + (xn[i] - x) * P[j-1][i+1]) / (xn[i] - xn[i+j]);
printf("%d %d %g %g %g %g %.20g\n", i, i+j, xn[i], xn[i+j], P[j-1][i], P[j-1][i+1], P[j][i]);
	
		}
	}

	return 0;
}

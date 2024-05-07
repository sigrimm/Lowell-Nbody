#include "define.h"

class planets{

public:
	int Nplanets;
	
	int *c1050;     //coefficients in group 1050, The size will be determined later
	int nc1050 = 0; //number of columns in group 1050
	char **cname;	//array with constant names

	

	int *id;		//id of the perturbers
	int *nChebyshev;	//Number of Chebyshev polynomials
	int *p_offset;		//Start point of individual perturbers in data array
	int *p_N;		//Number of records in desired time interval

	double *pertdata;


	int alloc();
	
	int readHeader(FILE *);
	int readPlanets(FILE *, double, double);
	int printPlanets(FILE *);

private:

};

#include "define.h"

class perturbers{

public:
	int Npert;


        int ND = 0;             //Number of double precision components in summaries
        int NI = 0;             //Number of integer components in summaries
        int NEXT = 10000000;    //record number of the next block
        int PREV = 0;           //record number of the last block

	int *id;		//id of the perturbers
	int *nChebyshev;	//Number of Chebyshev polynomials
	int *p_offset;		//Start point of individual perturbers in data array
	int *p_N;		//Number of records in desired time interval

	double *pertdata;


	int alloc();

	int readPerturbers1(FILE *);
	int readPerturbers2(FILE *, double, double);
	int printPerturbers(FILE *);

private:

};

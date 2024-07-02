#include "perturbers.h"
//https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html

int perturbers::alloc(){

	nChebyshev = (int*)malloc(Npert * sizeof(int));
	p_offset0 = (int*)malloc(Npert * sizeof(int));
	p_offset1 = (int*)malloc(Npert * sizeof(int));
	p_N = (int*)malloc(Npert * sizeof(int));
	id = (int*)malloc(Npert * sizeof(int));
	GM = (double*)malloc(Npert * sizeof(double));

	for(int i = 0; i < Npert; ++i){
		nChebyshev[i] = 0;
		GM[i] = 0.0;
	}

	return 1;
}


int perturbers::readPerturbers1(FILE *infile){

	//-------------------------------------------------------
	//Read the first file record, containing global information
	//The record size is 1024
	//The first records contains the information listed below:
	//See https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html
	//-------------------------------------------------------
	//+1 char
	char LOCIDW[9];			//Identification word
	//int ND;			//Number of double precision components in summaries
	//int NI;			//Number of integer components in summaries
	char LOCIFN[61];		//File name
	int FWARD;			//The record number of the first summary record
	int BWARD;			//The record number of the last summary record
	int FREE;			//First free address
	char LOCFMT[9];			//Binary format string
	char PRENUL[604];		//Block of Nulls
	char FTPSTR[29];		//FTP validation string
	char PSTNUL[298];		//Block of Nulls

	fread(LOCIDW, 8 * sizeof(char), 1, infile); 
	fread(&ND, sizeof(int), 1, infile); 
	fread(&NI, sizeof(int), 1, infile); 
	fread(LOCIFN, 60 * sizeof(char), 1, infile); 
	fread(&FWARD, sizeof(int), 1, infile); 
	fread(&BWARD, sizeof(int), 1, infile); 
	fread(&FREE, sizeof(int), 1, infile); 
	fread(LOCFMT, 8 * sizeof(char), 1, infile); 
	fread(PRENUL, 603 * sizeof(char), 1, infile); 
	fread(FTPSTR, 28 * sizeof(char), 1, infile); 
	fread(PSTNUL, 297 * sizeof(char), 1, infile); 

	printf("LOCIDW %s\n", LOCIDW);
	printf("ND %d\n", ND);
	printf("NI %d\n", NI);
	printf("LOCIFN %s\n", LOCIFN);
	printf("FWARD %d\n", FWARD);
	printf("BWARD %d\n", BWARD);
	printf("FREE %d\n", FREE);
	printf("LOCFMT %s\n", LOCFMT);
	printf("PRENUL |%s|\n", PRENUL);
	//printf("FTPSTR ");
	//for(int i = 0; i < 10; ++i){
	//	printf("%c", FTPSTR[i]);
	//}
	//printf("\n");
	printf("PSTNUL |%s|\n", PSTNUL);
	//-------------------------------------------------------


	//-------------------------------------------------------
	//Read comment area
	//-------------------------------------------------------
	char CDATA[DE_recordsize + 1];
	fread(CDATA, DE_recordsize * sizeof(char), 1, infile); 
	printf("START COMMENT\n");
	for(int i = 0; i < DE_recordsize; ++i){
		printf("%c", CDATA[i]);
	}
	printf("\nEND COMMENT\n");

	NEXT = FWARD;
	PREV = BWARD;

	return 1;
}


int perturbers::readPerturbers2(FILE *infile, FILE *outfile, double time0, double time1, int planetoffset){

	int aStart[Npert];	//start of the arrays
	int aEnd[Npert];	//end
	double Dtrecord[Npert];


	double data[DE_recordsize];


	//-------------------------------------------------------
	for(int k = 0; k < 100; ++k){

		//-------------------------------------------------------
		//Read Summary Records
		//A summary consists of ND = 6 double precision numbers 
		//and NI = 6 integer numbers.
		//These are:
		//1. Initial epoch of the interval in segment
		//2. Final epoch of the interval in segment
		//1. NAIF integer code for the target
		//2. NAIF integer code for the center (10 = Sun)
		//3. NAIF integer code for the reference frame (1 = J2000.0)
		//4. Integer code for the representation (data type, 2 = Chebyshev polynomials, position only, x y z components)
		//5. Initial address of the array
		//6. Final address of the array

		//https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html
		//https://www.gb.nrao.edu/ovlbi/spk.req

		//-------------------------------------------------------
		fseek(infile, (NEXT - 1) * DE_recordsize, SEEK_SET);
		double ddata[128];
		fread(ddata, 128 * sizeof(double), 1, infile); 

		int nSummaries = int(ddata[2]);
		int summarySize = ND + (NI + 1) / 2;

		NEXT = int(ddata[0]);
		PREV = int(ddata[1]);


		printf("NEXT %d\n", NEXT);
		printf("PREV %d\n", PREV);
		printf("Number of Summaries %d\n", nSummaries);
		printf("Summary Size %d\n", summarySize);
		printf("%d\n", nSummaries * summarySize);


		for(int i = 0; i < nSummaries; ++i){
			printf("%d ", i);
			double *d;
			int *idata;


			d = &ddata[i * summarySize + 3];
			double Time0 = d[0] / 86400.0 + 2451545.0; 
			double Time1 = d[1] / 86400.0 + 2451545.0; 
			printf("%.20g %.20g ", Time0, Time1);
			//for(int j = 0; j < ND; ++j){
			//	printf("%.20g ", d[j]);
			//}
			idata = (int *)&ddata[i * summarySize + 3 + ND];
			for(int j = 0; j < NI; ++j){
				printf("%d ", idata[j]);
			}
			//for(int j = 0; j < summarySize; ++j){
			//	printf("%.20g ", ddata[i * summarySize + j + 3]);
			//}
			printf("\n");

			if(Time0 <= time0 && Time1 > time0){
				for(int j = 0; j < Npert; ++j){
					if(id[j] == idata[0]){
						aStart[j] = idata[4];	//start of the array
						aEnd[j] = idata[5];	//end of the array
						printf("**** %d %d %d %d\n", i, j, id[j], aStart[j]);
					}
				}
			}

		}
		//-------------------------------------------------------

		//-------------------------------------------------------
		//Read Name Records
		//-------------------------------------------------------
		char CDATA[DE_recordsize + 1];
		fread(CDATA, DE_recordsize * sizeof(char), 1, infile); 
		int NC = 8 * (ND + (NI + 1) / 2);
		for(int i = 0; i < nSummaries; ++i){
			printf("%d ", i);
			for(int j = 0; j < NC; ++j){
				printf("%c", CDATA[i * NC + j]);
			}
			printf("\n");
		}
		//-------------------------------------------------------

		if(NEXT == 0){
			printf("\nReached the end of %s\n\n", "sb441-n16.bsp");
			break;
		}
	}	

	//-------------------------------------------------------
	//Scan through the data to extract the number of Chebyshev coefficients and 
	//count the number of epochs
	//-------------------------------------------------------
	dataSize = 0;
	for(int i = 0; i < Npert; ++i){

		//-------------------------------------------------------
		//Read Element Records
		//Elements contain:

		//Record 1
		//Record 2
		//...
		//Record N
		//INIT		Initial epoch of the first epoch
		//INTLEN	length of interval covered by each record
		//RSIZE		total size in each record, all records have the same size
		//N		Number of records contained in the segment


		//Each record contains:
		// Midtime of segment
		// Time Radius -> (MidTime - Radius) to (MidTime + Radius)
		// n Chebyshev coefficients for x component
		// n Chebyshev coefficients for y component
		// n Chebyshev coefficients for z component

		// The degree of the polynomial is n = (RSIZE - 2) / 3 - 1
		//-------------------------------------------------------

		//Read first the information at the end of the segment
		fseek(infile, (aEnd[i] - 4) * 8, SEEK_SET);
		double ddata[128];
		fread(ddata, 4 * sizeof(double), 1, infile);
		double INIT = ddata[0] / 86400.0 + 2451545.0; 
		double INTLEN = ddata[1] / 86400.0;
		int RSIZE = int(ddata[2]); 
		int N = int(ddata[3]); 

		//calculate the number of Chebyshev coefficients
		int NC = (RSIZE - 2) / 3;

		nChebyshev[i] = (nChebyshev[i] > NC) ? nChebyshev[i] : NC;
		Dtrecord[i] = INTLEN;

		int offset = floor((time0 - INIT) / INTLEN);

		printf("%d %d %.20g %.20g %d %d | %d\n", id[i], aStart[i], INIT, INTLEN, RSIZE, N, offset);

		fseek(infile, (aStart[i] + offset * RSIZE - 1) * 8, SEEK_SET);

		for(int j = 0; j < N; ++j){
			fread(ddata, RSIZE * sizeof(double), 1, infile); 

			double start = (ddata[0] - ddata[1])  / 86400.0 + 2451545.0;
			double end = (ddata[0] + ddata[1])  / 86400.0 + 2451545.0;
			double dt = ddata[1]  / 86400.0 * 2.0;

			if(start > time1){
				printf("%d %d\n", id[i],  j);
				p_offset0[i] = dataSize;
				p_N[i] = j;
				dataSize += j * RSIZE;
				p_offset1[i] = dataSize;
				//fprintf(outfile, "%d %d %d %d\n", id[i], nChebyshev[i][i], p_offset0[i] + planetoffset, p_offset1[i] + planetoffset);
				fwrite(&id[i], sizeof(int), 1, outfile);
				fwrite(&nChebyshev[i], sizeof(int), 1, outfile);
				int o0 = p_offset0[i] + planetoffset;
				int o1 = p_offset1[i] + planetoffset;
				fwrite(&o0, sizeof(int), 1, outfile);
				fwrite(&o1, sizeof(int), 1, outfile);
				fwrite(&GM[i], sizeof(double), 1, outfile);
				break;
			}

			//printf("%d %d %.20g %.20g | %.20g %.20g %.20g\n", id[i], j, time0, time1, start, end, dt);

			if(j >= N -1){

				printf("Error, reached the end of a block, Code needs to be improved to read in next block\n");
				return 0;
			}

		}
	}

	pertdata = (double*)malloc(dataSize * sizeof(double));

	for(int i = 0; i < Npert; ++i){
		//Read first the information at the end of the segment
		fseek(infile, (aEnd[i] - 4) * 8, SEEK_SET);
		double ddata[128];
		fread(ddata, 4 * sizeof(double), 1, infile);
		double INIT = ddata[0] / 86400.0 + 2451545.0; 
		double INTLEN = ddata[1] / 86400.0;
		int RSIZE = int(ddata[2]); 
		int N = int(ddata[3]); 

		int offset = floor((time0 - INIT) / INTLEN);

		//printf("%d %d %.20g %.20g %d %d | %d\n", id[i], aStart[i], INIT, INTLEN, RSIZE, N, offset);

		fseek(infile, (aStart[i] + offset * RSIZE - 1) * 8, SEEK_SET);

		for(int j = 0; j < p_N[i]; ++j){
			fread(pertdata + p_offset0[i] + j * RSIZE, RSIZE * sizeof(double), 1, infile); 
			
			double d0 = pertdata[p_offset0[i] + j * RSIZE];
			double d1 = pertdata[p_offset0[i] + j * RSIZE + 1];

			double start = (d0 - d1)  / 86400.0 + 2451545.0;
			double end = (d0 + d1)  / 86400.0 + 2451545.0;

			//overwrite midpoint time and time radius with start time and end time
			pertdata[p_offset0[i] + j * RSIZE] = start;
			pertdata[p_offset0[i] + j * RSIZE + 1] = end;

			//printf("%d %d %.20g %.20g %d \n", id[i], j, start, end, p_offset0[i] + j * RSIZE);
		}
	}
	return 1;
}

int perturbers::printPerturbers(FILE *outfile){

	//-------------------------------------------------------
	//Print all Chebyshev polynomial in a new file
	//-------------------------------------------------------
	for(int i = 0; i < Npert; ++i){
		int RSIZE = nChebyshev[i] * 3 + 2;
		for(int j = 0; j < p_N[i]; ++j){
			for(int k = 0; k < RSIZE; ++k){
				//fprintf(outfile, "%.20g ", pertdata[p_offset0[i] + j * RSIZE + k]);
				fwrite(&pertdata[p_offset0[i] + j * RSIZE + k], sizeof(double), 1, outfile);
			}
			//fprintf(outfile, "\n");
		}
	}
	return 1;
}


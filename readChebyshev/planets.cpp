#include "planets.h"


int planets::alloc(){
                
        nChebyshev = (int*)malloc(Nplanets * sizeof(int));
        p_offset0 = (int*)malloc(Nplanets * sizeof(int));
        p_offset1 = (int*)malloc(Nplanets * sizeof(int));
        p_N = (int*)malloc(Nplanets * sizeof(int));
        id = (int*)malloc((Nplanets  + Npert)* sizeof(int));
        GM = (double*)malloc((Nplanets + Npert) * sizeof(double));
                        
        for(int i = 0; i < Nplanets; ++i){
                nChebyshev[i] = 0;
        }               
        for(int i = 0; i < Nplanets + Npert; ++i){
		GM[i] = 0.0;
        }               
                
        return 1;
}


//search for block index
int getBlockIndex(FILE *infile, double time, int &blockIndex){
	for(int i = 2; i < 1000000; ++i){
		fseek(infile, i * DE_blocksize8, SEEK_SET);
		//read start and end time of the block
		double cdata[2];
		int er = fread(cdata, sizeof(double), 2, infile);
		if(er < 2){
			printf("Error in reading data\n");
			return 0;
		}
		//printf("Block %d %.20g %.20g %.20g\n", i, time, cdata[0], cdata[1]);
		if(time >= cdata[0] && time < cdata[1]){
			printf("Starting Block %d %.20g %.20g %.20g\n", i, time, cdata[0], cdata[1]);
			blockIndex = i;
			break;
		}
	}
	if(blockIndex == -1){
		printf("Error, start block index not found\n");
		return 0;
	}
	return 1;
}


int planets::readHeader(FILE *hfile){


	if(hfile){
		char line[256];
		std::string str;		
		//------------------------------------------
		//read until GROUP 1040
		//------------------------------------------
		while(fgets(line, 256, hfile)){
			str = line;
			//std::cout << "|" << str << "|\n";
			if (str.find("GROUP   1040") == 0){
				break;
			}

		}
		printf("reading group 1040\n");
		//read empty line
		fgets(line, 256, hfile);
		
		//read number of constants
		int nconst = 0;
		fscanf(hfile, "%d", &nconst);
		printf("Number of constants: %d\n", nconst);


		cname = (char**)malloc(nconst * sizeof(char *));
		for(int i = 0; i < nconst; ++i){
			cname[i] = (char*)malloc(16);
			fscanf(hfile, "%s", cname[i]);
			//printf("%s ", cname[i]);
		}
		//printf("/n");
		//------------------------------------------

		
		/* /------------------------------------------
		//read until GROUP 1041
		//reading the numbers from the header file causes a loss of precision
		//So we skip this part and read the numbers from the binary file
		//------------------------------------------
		while(fgets(line, 256, hfile)){
			str = line;
			//std::cout << "|" << str << "|\n";
			if (str.find("GROUP   1041") == 0){
				break;
			}

		}
		printf("reading group 1041\n");
		//read empty line
		fgets(line, 256, hfile);
		
		//read number of constants
		nconst = 0;
		fscanf(hfile, "%d", &nconst);
		printf("%d\n", nconst);


		double cvalue[nconst];

		char svalue[32];
		double gvalue;
		for(int i = 0; i < nconst; ++i){
			
			fscanf(hfile, "%s", svalue);
			std::string s = svalue;
			std::replace(s.begin(), s.end(), 'D', 'E');
			gvalue = std::stod(s);
			if(i<10){
				printf("%.20e %s\n", gvalue, svalue);
				std::cout << s << std::endl;
			}
		}
		*/ //------------------------------------------

		//------------------------------------------
		//read until GROUP 1050
		//the number of columns of group 1050 is determined and stored in nc1050
		//The group 1050 consists of 3 rows:
		//row 1: offset of planet in data array
		//row 2: number of chebyshev coeffitients per entry
		//row 3: numer of subintervals per entry
		//------------------------------------------
		while(fgets(line, 256, hfile)){
			str = line;
			//std::cout << "|" << str << "|\n";
			if (str.find("GROUP   1050") == 0){
				break;
			}

		}
		printf("reading group 1050\n");
		//read empty line
		fgets(line, 256, hfile);
		//read first line and count the number of columns
		int offset = 0;
		int t, n;
		fgets(line, 256, hfile);
		//Determine the number of columns of group 1050
		//n indicates the number of character read, 
		//this is then used to offset the starting point of the string to read
		while(sscanf(line + offset, "%d%n", &t, &n)>0){
			//printf("%d %d\n", t, n);
			offset += n;
			++nc1050;
		}
		printf("number of columns in group 1050 %d\n", nc1050);

		c1050 = (int*)malloc(nc1050 * 3 * sizeof(int));
		
		for(int i = 0; i < 3; ++i){
			//read until the end of the line
			//read 3 lines from file
			offset = 0;
			int nc = 0;
			while(sscanf(line + offset, "%d%n", &t, &n)>0){
				c1050[nc + i * nc1050] = t;
				//printf("%d %d\n", t, c1050[nc + i * nc1050]);
				offset += n;
				++nc;
			}
			fgets(line, 256, hfile);
		}
		//------------------------------------------
		
		printf("reading header file finished\n");
	}

	return 0;
}	

int planets::readPlanets(FILE *infile, FILE *outfileT, FILE *outfile, double time0, double time1){


	double startTime, endTime, daysPerBlock;
	int Nc;
	double AUtokm, EM;

	//The offsets are explained in
	//https://www.celestialprogramming.com/jpl-ephemeris-format/jpl-ephemeris-format.html

	printf("seek %d\n", 0x0A5C);	//2652 offset to StartJD
	fseek(infile, 0x0A5C, SEEK_SET);

	fread(&startTime, sizeof(double), 1, infile); 
	fread(&endTime, sizeof(double), 1, infile); 
	fread(&daysPerBlock, sizeof(double), 1, infile); 
	fread(&Nc, sizeof(int), 1, infile); 
	fread(&AUtokm, sizeof(double), 1, infile);		//AU 
	fread(&EM, sizeof(double), 1, infile);			//EMRAT 


	printf("Start time %.20g\n", startTime);
	printf("End time %.20g\n", endTime);
	printf("days per block %.20g\n", daysPerBlock);
	printf("Number of constants %d\n", Nc);
	printf("AUtokm %.20e\n", AUtokm);
	printf("EM %.20e\n", EM);


	printf("read coefficients of group 1050\n");
	//read 12 columns of group 1050
	int nSub[12];		//number of sub intervals per planet
	int nChebyshevMax = 0;
	int body_offset[12];

	for(int i = 0; i < 12; ++i){
		int c1, c2, c3;
		fread(&c1, sizeof(int), 1, infile);
		fread(&c2, sizeof(int), 1, infile);
		fread(&c3, sizeof(int), 1, infile);
		body_offset[i] = c1 - 1;	//minus 1 to start at 0
		if(i < 12){
			nChebyshev[i] = c2;
			nChebyshevMax = (c2 > nChebyshevMax) ? c2 : nChebyshevMax;
		}
		nSub[i] = c3;
		printf("g1050 column %d %d %d %d\n", i, c1, c2, c3);
		//printf("%d %d %d %d\n", i, c1050[i + 0 * nc1050], c1050[i + 1 * nc1050], c1050[i + 2 * nc1050]);
		if(c1 != c1050[i + 0 * nc1050]){
			printf("Error, coefficients in group 1050 do not agree\n");
			return 0;
		}
		if(c2 != c1050[i + 1 * nc1050]){
			printf("Error, coefficients in group 1050 do not agree\n");
			return 0;
		}
		if(c3 != c1050[i + 2 * nc1050]){
			printf("Error, coefficients in group 1050 do not agree\n");
			return 0;
		}
	}
	printf("Maximum number of Chebyshev coefficients: %d\n", nChebyshevMax);

	//read version number
	int version;
	fread(&version, sizeof(int), 1, infile);
	printf("Version of data file: %d\n", version);

	if(DE_version != version){
		printf("Error, DE version do not agree! %d %d\n", DE_version, version);
		return 0;
	}
	//read  column 13 of group 1050
	for(int i = 12; i < 13; ++i){
		int c1, c2, c3;
		fread(&c1, sizeof(int), 1, infile);
		fread(&c2, sizeof(int), 1, infile);
		fread(&c3, sizeof(int), 1, infile);
		printf("g1050 column %d %d %d %d\n", i, c1, c2, c3);
		//printf("%d %d %d %d\n", i, c1050[i + 0 * nc1050], c1050[i + 1 * nc1050], c1050[i + 2 * nc1050]);
		if(c1 != c1050[i + 0 * nc1050]){
			printf("Error, coefficients in group 1050 do not agree\n");
			return 0;
		}
		if(c2 != c1050[i + 1 * nc1050]){
			printf("Error, coefficients in group 1050 do not agree\n");
			return 0;
		}
		if(c3 != c1050[i + 2 * nc1050]){
			printf("Error, coefficients in group 1050 do not agree\n");
			return 0;
		}
	}

	//read constants
	fseek(infile, DE_blocksize8, SEEK_SET);
	double dd;

	double AU, EMRAT, CLIGHT, RE, J2E;
	for(int i = 0; i < Nc; ++i){
		fread(&dd, sizeof(double), 1, infile);
		//Number of kilometers per astronomical unit
		if(strcmp(cname[i], "AU") == 0){
			AU = dd;
			printf("%s %.20e\n", cname[i], dd);
		}
		//speed of light
		if(strcmp(cname[i], "CLIGHT") == 0){
			CLIGHT = dd;
			printf("%s %.20e\n", cname[i], dd);
		}
		//Radius of Earth
		if(strcmp(cname[i], "RE") == 0){
			RE = dd;
			printf("%s %.20e\n", cname[i], dd);
		}
		//J2 of Earth
		if(strcmp(cname[i], "J2E") == 0){
			J2E = dd;
			printf("%s %.20e\n", cname[i], dd);
		}
		//Earth-Moon mass ratio
		if(strcmp(cname[i], "EMRAT") == 0){
			EMRAT = dd;
			printf("%s %.20e\n", cname[i], dd);
		}
		//Mass of Sun
		if(strcmp(cname[i], "GMS") == 0){
			GM[10] = dd;
			printf("%s %.20e\n", cname[i], dd / def_k2);
		}
		//Mass of Earth Moon barycenter
		if(strcmp(cname[i], "GMB") == 0){
			double GMB = dd;
			printf("%s %.20e\n", cname[i], dd / def_k2);

			// E + M = B
			//EM = E/M

			//Mass of Earth
			// M = E/EM
			//E + E/EM = B
			//E(EM+1) = B * EM
			GM[2] = (EMRAT / (1.0 + EMRAT)) * GMB;
			printf("%s %.20e %.20g\n", "GME", GM[2] / def_k2, GM[2]);
		
			//Mass of Moon
			//E = EM * M
			//EM * M + M = B
			//M(EM + 1) = B
			GM[9] = 1.0 / (1.0 + EMRAT) * GMB;
			printf("%s %.20e\n", "GM9", GM[9] / def_k2);
		}
		//Mass of Mercury
		if(strcmp(cname[i], "GM1") == 0){
			GM[0] = dd;
			printf("%s %.20e\n", cname[i], dd / def_k2);
		}
		//Mass of Venus
		if(strcmp(cname[i], "GM2") == 0){
			GM[1] = dd;
			printf("%s %.20e\n", cname[i], dd / def_k2);
		}
		//Mass of Mars
		if(strcmp(cname[i], "GM4") == 0){
			GM[3] = dd;
			printf("%s %.20e\n", cname[i], dd / def_k2);
		}
		//Mass of Jupiter
		if(strcmp(cname[i], "GM5") == 0){
			GM[4] = dd;
			printf("%s %.20e\n", cname[i], dd / def_k2);
		}
		//Mass of Saturn
		if(strcmp(cname[i], "GM6") == 0){
			GM[5] = dd;
			printf("%s %.20e\n", cname[i], dd / def_k2);
		}
		//Mass of Uranus
		if(strcmp(cname[i], "GM7") == 0){
			GM[6] = dd;
			printf("%s %.20e\n", cname[i], dd / def_k2);
		}
		//Mass of Neptune
		if(strcmp(cname[i], "GM8") == 0){
			GM[7] = dd;
			printf("%s %.20e\n", cname[i], dd / def_k2);
		}
		//Mass of Pluto
		if(strcmp(cname[i], "GM9") == 0){
			GM[8] = dd;
			printf("%s %.20e\n", cname[i], dd / def_k2);
		}
		for(int j = 0; j < Npert; ++j){
			char iid[7];
			sprintf(iid, "MA%04d", id[Nplanets + j]);
			if(strcmp(cname[i], iid) == 0){
				GM[Nplanets + j] = dd;
				printf("Mass %d %s %g\n", j, iid, dd);
			}
		}

	
	}

#if def_printT == 1
	fprintf(outfileT,"time0 %.20g\n", time0);
	fprintf(outfileT,"time1 %.20g\n", time1);
	fprintf(outfileT,"AUtokm %.20g\n", AUtokm);
	fprintf(outfileT,"EM %.20g\n", EM);
	fprintf(outfileT,"CLIGHT %.20g\n", CLIGHT);
	fprintf(outfileT,"RE %.20g\n", RE);
	fprintf(outfileT,"J2E %.20g\n", J2E);
#endif

	fwrite(&time0, sizeof(double), 1, outfile);
	fwrite(&time1, sizeof(double), 1, outfile);
	fwrite(&AUtokm, sizeof(double), 1, outfile);
	fwrite(&EM, sizeof(double), 1, outfile);
	fwrite(&CLIGHT, sizeof(double), 1, outfile);
	fwrite(&RE, sizeof(double), 1, outfile);
	fwrite(&J2E, sizeof(double), 1, outfile);

	double *cdata;
	cdata = (double*)malloc(DE_blocksize * sizeof(double));	
	int er;

	//search for start block
	int blockIndex = -1;
	er = getBlockIndex(infile, time0, blockIndex);
	if(er <= 0){
		printf("Error in finding start block\n");
		return 0;
	}

	int dblock = (time1 - time0) / daysPerBlock + 1;
	printf("dblock %d\n", dblock);


	dataSize = 0;

	for(int i = 0; i < Nplanets; ++i){
		p_offset0[i] = dataSize;
		p_N[i] = 0;
		dataSize += nSub[i] *  (nChebyshev[i] * 3 + 2) * dblock;
		p_offset1[i] = dataSize;
		printf("planet offset %d %d %d %d %d %d\n", i, nSub[i], nChebyshev[i], p_offset0[i], p_offset1[i], nSub[i] *  (nChebyshev[i] * 3 + 2) * dblock);
#if def_printT == 1
		fprintf(outfileT, "%d %d %d %d %g\n", id[i], nChebyshev[i], p_offset0[i], p_offset1[i], GM[i]);
#endif
//printf("%d %d %d %d %.20g\n", id[i], nChebyshev[i], p_offset0[i], p_offset1[i], GM[i]);
		fwrite(&id[i], sizeof(int), 1, outfile);
		fwrite(&nChebyshev[i], sizeof(int), 1, outfile);
		fwrite(&p_offset0[i], sizeof(int), 1, outfile);
		fwrite(&p_offset1[i], sizeof(int), 1, outfile);
		fwrite(&GM[i], sizeof(double), 1, outfile);

	}

	pertdata = (double*)malloc(dataSize * sizeof(double));

	for(int i = 0; i < dblock; ++i){
		//read data block
		fseek(infile, blockIndex * DE_blocksize8, SEEK_SET);
		er = fread(cdata, sizeof(double), DE_blocksize, infile);
		++blockIndex;

		if(er < DE_blocksize){
			printf("Error in reading data\n");
			return 0;
		}

		if(cdata[0] >= time1){
			printf("Reached the end of the time range\n");
			break;
		}

		if(cdata[0] >= endTime){
			printf("Reached the end of the data file\n");
			return 0;
		}

		printf("Reading Block %d %.20g %.20g | %.20g %.20g\n", i, time0, time1, cdata[0], cdata[1]);
	
		
		for(int j = 0; j < Nplanets; ++j){

			int sizeSubInterval = daysPerBlock / nSub[j];
			for(int iSub = 0; iSub < nSub[j]; ++iSub){

				int offset = body_offset[j] + iSub * nChebyshev[j] * 3;

				double subStart = cdata[0] + iSub * sizeSubInterval;	//starting time of sub interval
				double subEnd = subStart + sizeSubInterval;	//ending time of sub interval
				int RSIZE = nChebyshev[j] * 3 + 2;

				//printf("%d %.20g %.20g\n", j, subStart, subEnd);

				pertdata[p_offset0[j] + p_N[j] * RSIZE] = subStart;
				pertdata[p_offset0[j] + p_N[j] * RSIZE + 1] = subEnd;

				for(int k = 0; k < nChebyshev[j] * 3; ++k){
					pertdata[p_offset0[j] + p_N[j] * RSIZE + 2 + k] = cdata[offset + k];
				}
				++p_N[j]; 
			}

		}
	
	}


	return 0;
}


int planets::printPlanets(FILE *outfileT, FILE *outfile){

	//-------------------------------------------------------
	//Print all Chebyshev polynomial in a new file
	//-------------------------------------------------------
	for(int i = 0; i < Nplanets; ++i){
		int RSIZE = nChebyshev[i] * 3 + 2;
		for(int j = 0; j < p_N[i]; ++j){
			for(int k = 0; k < RSIZE; ++k){
#if def_printT == 1
				fprintf(outfileT, "%.20g ", pertdata[p_offset0[i] + j * RSIZE + k]);
#endif
				fwrite(&pertdata[p_offset0[i] + j * RSIZE + k], sizeof(double), 1, outfile);
			}
#if def_printT == 1
			fprintf(outfileT, "\n");
#endif
		}
	}
	return 1;
}

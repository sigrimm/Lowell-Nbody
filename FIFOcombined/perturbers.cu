#include "Host.h"
__host__ void Host::perturbersMass(){

	/*
	double mass1[] = {
		//masses form DE440
		//in units of GM
		1.3271244004127942E+11, //Sun           10
		2.2031868551400003E+04, //Mercury       1
		3.2485859200000000E+05, //Venus         2
		3.9860043550702266E+05, //Earth         399
		4.2828375815756102E+04, //Mars          4
		1.2671276409999998E+08, //Jupiter       5
		3.7940584841799997E+07, //Saturn        6
		5.7945563999999985E+06, //Uranus        7
		6.8365271005803989E+06, //Neptune       8
		9.7550000000000000E+02, //Pluto         9
		4.9028001184575496E+03, //Moon          301

		6.2628888644409933E+01, //Ceres         1
		1.3665878145967422E+01, //Pallas        2
		1.9205707002025889E+00, //Juno          3
		1.7288232879171513E+01, //Vesta         4
		1.1398723232184107E+00, //Iris          7
		5.6251476453852289E+00, //Hygiea        10
		2.0230209871098284E+00, //Eunomia       15
		1.5896582441709424E+00, //Psyche        16
		1.0793714577033560E+00, //Euphrosyne    31
		2.6830359242821795E+00, //Europa        52
		9.3810575639151328E-01, //Cybele        65
		2.1682320736996910E+00, //Sylvia        87
		1.1898077088121908E+00, //Thisbe        88
		1.4437384031866001E+00, //Camilla       107
		3.8944831481705644E+00, //Davida        511
		2.8304096393299849E+00, //Interamnia    704
	};
	for(int i = 0; i < Nperturbers; ++i){
		m_h[i] = mass1[i] / mass1[0];
//printf("m %d %.20g\n", i, m_h[i]);
	}
	*/

	//Units are 1/(mass of object in solar masses)
	double pmass[] = {
		//This values come from the Larry Wasserman aiset.f
		//From DE-440
		//The order must agree with the perturbers coordinate file
		1.000000000000000e0,      	// Sun        10
		6.02365794492933620e+06,	// Mercury    1
		4.08523718656268204e+05,	// Venus      2
		3.32946048773048213e+05,	// Earth      399
		3.09870354673725273e+06,	// Mars       4
		1.04734863124400454e+03,	// Jupiter    5
		3.49790180079320044e+03,	// Saturn     6
		2.29029507834766191e+04,	// Uranus     7
		1.94122597758754673e+04,	// Neptune    8
		1.36045556167380244e+08,	// Pluto      9
		2.70687029523511454e+07,	// Moon       301

		2.11902913996742201e+09,	// Ceres      1
		9.71122664959812355e+09,	// Pallas     2
		6.91005231035131531e+10,	// Juno       3
		7.67646068680440331e+09,	// Vesta      4
		1.16427460635738602e+11,	// Iris       7
		2.35927034111103439e+10,	// Hygiea     10
		6.56011187658897781e+10,	// Eunomia    15
		8.34848877284898376e+10,	// Psyche     16
		1.22953445817123718e+11,	// Euphrosyne 31
		4.94635345133462448e+10,	// Europa     52
		1.41468527548287018e+11,	// Cybele     65
		6.12076731319770355e+10,	// Sylvia     87
		1.11541082696269424e+11,	// Thisbe     88
		9.19227747550098419e+10,	// Camilla    107
		3.40770353836609001e+10,	// Davida     511
		4.68880681429191055e+10,	// Interamnia 704
	};
	for(int i = 0; i < Nperturbers; ++i){
		m_h[i] = 1.0/pmass[i];
printf("m %d %.20g\n", i, m_h[i]);
//printf("m %d %.20g %.20g\n", i, m_h[i], mass1[i] / mass1[0]);
	}

	// *********************************************
	//Old masses and perturbers
	// *********************************************
	//Units are 1/(mass of object in solar masses)
	/*
	double pmass[] = {
		//The sun has to be at the first position
		1.000000000000000e0,      // Sun        (0)
		6023682.155592479e0,      // Mercury    (1)
		408523.7186582996e0,      // Venus      (2)
		332946.0488339480e0,      // Earth      (3)
		3098703.590290707e0,      // Mars       (4)
		1047.348625463337e0,      // Jupiter    (5)
		3497.901767786633e0,      // Saturn     (6)
		22902.98161308703e0,      // Uranus     (7)
		19412.25977597307e0,      // Neptune    (8)
		135836683.7686175e0,      // Pluto      (9)
		2112939391.8192508,                  // Ceres      (10)
		9531877787.0654011,                  // Pallas     (11)
		81799329362.428986,                  // Juno       (12)
		7676559929.1351004,                  // Vesta      (13)
		23944976514.662392,                  // Hygiea     (14)
		63251980219.354561,                  // Eunomia    (15)
		46649712166.264168,                  // Euphrosyne (16)
		119474172269.94408,                  // Europa     (17)
		56926698684.931702,                  // Davida     (18)
		56298080671.641434,                  // Interamnia (19)
		27068703.24120323e0,      // Moon       (20)

		86737410876.841156,               // Psyche
		93034865412.812271,               // Cybele
		114823090351.20033,               // Thisbe
		116910898662.48077,               // Doris
		128906361339.41116,               // Patientia
		134548655333.38321,               // Sylvia

		0.0                       // test particle

	};
	for(int i = 0; i < Nperturbers; ++i){
		m_h[i] = 1.0/pmass[i];
//printf("m %d %.20g\n", i, m_h[i]);
	}
	*/


}


__host__ void Host::perturbersIDs(){
	//perturbers indices
	unsigned long long int id[] = {
		10,
		1,
		2,
		399,
		4,
		5,
		6,
		7,
		8,
		72057594038083239,      // 9,
		301,
		72057594037948899,	//Ceres      1
		72057594037948900,	//Pallas     2
		72057594037948901,	//Juno       3
		72057594037948902,	//Vesta      4
		72057594037948905,	//Iris       7
		72057594037948908,	//Hygiea     10
		72057594037948913,	//Eunomia    15
		72057594037948914,	//Psyche     16
		72057594037948929,	//Euphrosyne 31
		72057594037948950,	//Europa     52
		72057594037948963,	//Cybele     65
		72057594037948985,	//Sylvia     87
		72057594037948986,	//Thisbe     88
		72057594037949005,	//Camilla    107
		72057594037949409,	//Davida     511
		72057594037949602,	//Interamnia 704
	};

	for(int i = 0; i < Nperturbers; ++i){
		id_h[i] = id[i];
	}
}


__global__ void bufferToX_kernel(double *XYdata_d, double *timep_d, double *xp_d, double *yp_d, double *zp_d, int N){

	int id = threadIdx.x + blockDim.x * blockIdx.x;

	if(id < N){
		timep_d[id] = XYdata_d[id * 4];
		xp_d[id] = XYdata_d[id * 4 + 1];
		yp_d[id] = XYdata_d[id * 4 + 2];
		zp_d[id] = XYdata_d[id * 4 + 3];
if(id < 30 || id > N - 10 ) printf("buffer %d %.20g %.20g %.20g %.20g\n", id, timep_d[id], xp_d[id], yp_d[id], zp_d[id]);
	}
}
__host__ void bufferToX(double *XYdata_h, double *timep_h, double *xp_h, double *yp_h, double *zp_h, int N){


	for(int id = 0; id < N; ++id){
		timep_h[id] = XYdata_h[id * 4];
		xp_h[id] = XYdata_h[id * 4 + 1];
		yp_h[id] = XYdata_h[id * 4 + 2];
		zp_h[id] = XYdata_h[id * 4 + 3];
if(id < 30 || id > N - 10) printf("buffer %d %.20g %.20g %.20g %.20g\n", id, timep_h[id], xp_h[id], yp_h[id], zp_h[id]);
	}
}



//read the perturbers table
__host__ int Host::readTable(){

	if(useHelio == 1){
		XVfile = fopen("All_h.bin", "rb");
	}
	else{
		XVfile = fopen("All_b.bin", "rb");
	}

	if(XVfile == NULL){
		printf("Error, perturber file not found\n");
		return 0;
	}

	int NTableC = 0;
	//use two data buffers for interleaved data transfer
	for(int t = 0; t < 1000000; ++t){

	//for(int t = 0; t < 10; ++t){
//printf("t %d\n", t);
		int er;

		if(t % 2 == 0){
//printf("start read A\n");
			er = fread(readBufferA_h, Nperturbers * 4 * sizeof(double), 1, XVfile);
//printf("end read A\n");
		}
		else{
//printf("start read B\n");
			er = fread(readBufferB_h, Nperturbers * 4 * sizeof(double), 1, XVfile);
//printf("end read B\n");
		}

		/*
		//only here for checking
		if(t < 4){
			for(int i = 0; i < Nperturbers; ++i){
				if(t % 2 == 0) printf("XYa %d %.20g %g %g %g\n", i, readBufferA_h[i * 4 + 0], readBufferA_h[i * 4 + 1], readBufferA_h[i * 4 + 2], readBufferA_h[i * 4 + 3]);
				if(t % 2 == 1) printf("XYb %d %.20g %g %g %g\n", i, readBufferB_h[i * 4 + 0], readBufferB_h[i * 4 + 1], readBufferB_h[i * 4 + 2], readBufferB_h[i * 4 + 3]);
			}
		}
		*/

		if(t == 0){
			//set start time of perturbers file
			timep0 = readBufferA_h[0];
		}

		if(useGPU == 0){
			if(t % 2 == 0){
				memcpy(XYdata_h + t * Nperturbers * 4, readBufferA_h, Nperturbers * 4 * sizeof(double));
			}
			else{
				memcpy(XYdata_h + t * Nperturbers * 4, readBufferB_h, Nperturbers * 4 * sizeof(double));
			}
		}
		else{
			cudaDeviceSynchronize(); //this must be here

			//both buffers A and B use the same stream, so copy can overlap with the next read, but not with the
			//next copy. 
			if(t % 2 == 0){
//printf("start copy A\n");
				cudaMemcpyAsync(XYdata_d + t * Nperturbers * 4, readBufferA_h, Nperturbers * 4 * sizeof(double), cudaMemcpyHostToDevice);
			}
			else{
//printf("start copy B\n");
				cudaMemcpyAsync(XYdata_d + t * Nperturbers * 4, readBufferB_h, Nperturbers * 4 * sizeof(double), cudaMemcpyHostToDevice);
			}
		}

		if(er <= 0){
//printf("readbuffer %d %d %d %.20g %g %g %g\n", er, t, bSwap, readBuffer_h[0], readBuffer_h[1], readBuffer_h[2], readBuffer_h[3]);
			NTableC = t;
			break;
		}

//if(t < 4) printf("readbuffer %d %d %d %.20g %g %g %g\n", er, t, bSwap, readBuffer_h[0], readBuffer_h[1], readBuffer_h[2], readBuffer_h[3]);
	}

	printf("NTableC: %d\n", NTableC);

	fclose(XVfile);
	cudaDeviceSynchronize();

	if(useGPU == 0){
		bufferToX (XYdata_h, timep_h, xp_h, yp_h, zp_h, NTableC * Nperturbers);
		free(XYdata_h);
		free(readBufferA_h);
		free(readBufferB_h);
	}
	else{
		bufferToX_kernel <<< (NTableC * Nperturbers + 127) / 128, 128 >>> (XYdata_d, timep_d, xp_d, yp_d, zp_d, NTableC * Nperturbers);
		cudaFree(XYdata_d);
		cudaFreeHost(readBufferA_h);
		cudaFreeHost(readBufferB_h);
	}

	return 1;
}

//needed with useGPU > 0 but still using the CPU for the preintegration
__host__ void Host::copyPerturbersCPU(){
	cudaMemcpy(timep_h, timep_d, Nperturbers * NTable * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(xp_h, xp_d, Nperturbers * NTable * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(yp_h, yp_d, Nperturbers * NTable * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(zp_h, zp_d, Nperturbers * NTable * sizeof(double), cudaMemcpyDeviceToHost);
}

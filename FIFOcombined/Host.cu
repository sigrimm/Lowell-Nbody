#include "Host.h"

//initialize variables
__host__ Host::Host(){

	timep0 = 0.0;

	comet = 0;

	useGR = 2;
	useJ2 = 1;
	useNonGrav = 1;

	useGPU = 1;
	useAdaptiveTimeSteps = 1;
	useIndividualTimeSteps = 0;

	useFIFO = 2;
	InVersion= 0.0;
	DoPreIntegration = 0;	

	useHelio = 1;
	outHelio = 1;

	outBinary = 0;

	Nsteps = 1e9;
	outInterval = 1;


	time0 = 0.0;
	time1 = 0.0;
	outStart = 0.0;
	time = 0.0;

	dti = 0.05;
	dt = dti * dayUnit;

	dts = 0.01;
	outI = 1llu;

	N = Nperturbers;
	NMax = 10000000;

	nRuns = 1;
	dtiMin = new double[10];
	runsN = new int[10];		//assume maximally 10 levels
	runsdt = new double[10];
	for(int i = 0; i < 10; ++i){
		runsN[i] = Nperturbers;
		runsdt[i] = dti;
		dtiMin[i] = 0.01;
	}

	outfilename = new char[160];
	dtfilename = new char[160];
	sprintf(dtfilename, "timesteps.dat");

	//Erase Outbinary file
	if(outBinary == 1){
		if(outHelio == 1){
			sprintf(outfilename, "Outhelio.bin");
		}
		else{
			sprintf(outfilename, "Outbary.bin");
		}
		outfile = fopen(outfilename, "wb");
		fclose(outfile);
	}

}

__host__ int Host::readparam(int argc, char*argv[]){

	FILE *paramfile;
	paramfile = fopen("param.dat", "r");


        char sp[160];
        int er;

	for(int j = 0; j < 1000; ++j){ //loop around all lines in the param.dat file
		int c;
		for(int i = 0; i < 50; ++i){
			c = fgetc(paramfile);
			if(c == EOF) break;
			sp[i] = char(c);
			if(c == '=' || c == ':'){
				sp[i + 1] = '\0';
				break;
			}
		}
		if(c == EOF) break;

		if(strcmp(sp, "Number of time step levels =") == 0){
			er = fscanf (paramfile, "%d", &nRuns);
			if(er <= 0){
				printf("Error: Number of time step levels is not valid!\n");
				return 0;
			}
			fgets(sp, 3, paramfile);
		}
		else if(strcmp(sp, "dt =") == 0){
			for(int t = 0; t < nRuns; ++t){
				er = fscanf (paramfile, "%lf", &runsdt[t]);
				if(er <= 0){
					printf("Error: dtMin value is not valid!\n");
					return 0;
				}
			}
			fgets(sp, 3, paramfile);
		}
		else if(strcmp(sp, "dtMin =") == 0){
			for(int t = 0; t < nRuns; ++t){
				er = fscanf (paramfile, "%lf", &dtiMin[t]);
				if(er <= 0){
					printf("Error: dtMin value is not valid!\n");
					return 0;
				}
			}
			fgets(sp, 3, paramfile);
		}
	}




	fclose(paramfile);
	return 1;
}


__host__ void Host::Alloc1(){
	//allocate data on host
	id_h = (unsigned long long int*)malloc(N * sizeof(unsigned long long int));
	index_h = (unsigned int*)malloc(N * sizeof(unsigned int));
	m_h = (double*)malloc(N * sizeof(double));
	x_h = (double*)malloc(N * sizeof(double));
	y_h = (double*)malloc(N * sizeof(double));
	z_h = (double*)malloc(N * sizeof(double));
	vx_h = (double*)malloc(N * sizeof(double));
	vy_h = (double*)malloc(N * sizeof(double));
	vz_h = (double*)malloc(N * sizeof(double));
	A1_h = (double*)malloc(N * sizeof(double));
	A2_h = (double*)malloc(N * sizeof(double));
	A3_h = (double*)malloc(N * sizeof(double));
	jd_init_h = (double*)malloc(N * sizeof(double));

	timep_h = (double*)malloc(Nperturbers * NTable * sizeof(double));
	xp_h = (double*)malloc(Nperturbers * NTable * sizeof(double));
	yp_h = (double*)malloc(Nperturbers * NTable * sizeof(double));
	zp_h = (double*)malloc(Nperturbers * NTable * sizeof(double));

	x0_h = (double*)malloc(N * sizeof(double));
	y0_h = (double*)malloc(N * sizeof(double));
	z0_h = (double*)malloc(N * sizeof(double));
	vx0_h = (double*)malloc(N * sizeof(double));
	vy0_h = (double*)malloc(N * sizeof(double));
	vz0_h = (double*)malloc(N * sizeof(double));
	A10_h = (double*)malloc(N * sizeof(double));
	A20_h = (double*)malloc(N * sizeof(double));
	A30_h = (double*)malloc(N * sizeof(double));
	m0_h = (double*)malloc(N * sizeof(double));
	id0_h = (unsigned long long int*)malloc(N * sizeof(unsigned long long int));
	index0_h = (unsigned int*)malloc(N * sizeof(unsigned int));

	xb_h = (double*)malloc(N * sizeof(double));
	yb_h = (double*)malloc(N * sizeof(double));
	zb_h = (double*)malloc(N * sizeof(double));
	vxb_h = (double*)malloc(N * sizeof(double));
	vyb_h = (double*)malloc(N * sizeof(double));
	vzb_h = (double*)malloc(N * sizeof(double));
	A1b_h = (double*)malloc(N * sizeof(double));
	A2b_h = (double*)malloc(N * sizeof(double));
	A3b_h = (double*)malloc(N * sizeof(double));
	mb_h = (double*)malloc(N * sizeof(double));
	idb_h = (unsigned long long int*)malloc(N * sizeof(unsigned long long int));

	dx_h = (double*)malloc(N * sizeof(double));
	dy_h = (double*)malloc(N * sizeof(double));
	dz_h = (double*)malloc(N * sizeof(double));
	dvx_h = (double*)malloc(N * sizeof(double));
	dvy_h = (double*)malloc(N * sizeof(double));
	dvz_h = (double*)malloc(N * sizeof(double));

	if(useGPU > 0){
		//allocate data on the device
		cudaMalloc((void **) &id_d, N * sizeof(unsigned long long int));
		cudaMalloc((void **) &index_d, N * sizeof(unsigned int));
		cudaMalloc((void **) &m_d, N * sizeof(double));
		cudaMalloc((void **) &x_d, N * sizeof(double));
		cudaMalloc((void **) &y_d, N * sizeof(double));
		cudaMalloc((void **) &z_d, N * sizeof(double));
		cudaMalloc((void **) &vx_d, N * sizeof(double));
		cudaMalloc((void **) &vy_d, N * sizeof(double));
		cudaMalloc((void **) &vz_d, N * sizeof(double));
		cudaMalloc((void **) &A1_d, N * sizeof(double));
		cudaMalloc((void **) &A2_d, N * sizeof(double));
		cudaMalloc((void **) &A3_d, N * sizeof(double));
		cudaMalloc((void **) &jd_init_d, N * sizeof(double));

		cudaMalloc((void **) &id0_d, N * sizeof(unsigned long long int));
		cudaMalloc((void **) &index0_d, N * sizeof(unsigned int));
		cudaMalloc((void **) &m0_d, N * sizeof(double));
		cudaMalloc((void **) &x0_d, N * sizeof(double));
		cudaMalloc((void **) &y0_d, N * sizeof(double));
		cudaMalloc((void **) &z0_d, N * sizeof(double));
		cudaMalloc((void **) &vx0_d, N * sizeof(double));
		cudaMalloc((void **) &vy0_d, N * sizeof(double));
		cudaMalloc((void **) &vz0_d, N * sizeof(double));
		cudaMalloc((void **) &A10_d, N * sizeof(double));
		cudaMalloc((void **) &A20_d, N * sizeof(double));
		cudaMalloc((void **) &A30_d, N * sizeof(double));

		cudaMalloc((void **) &timep_d, Nperturbers * NTable * sizeof(double));
		cudaMalloc((void **) &xp_d, Nperturbers * NTable * sizeof(double));
		cudaMalloc((void **) &yp_d, Nperturbers * NTable * sizeof(double));
		cudaMalloc((void **) &zp_d, Nperturbers * NTable * sizeof(double));

		cudaMalloc((void **) &dx_d, N * sizeof(double));
		cudaMalloc((void **) &dy_d, N * sizeof(double));
		cudaMalloc((void **) &dz_d, N * sizeof(double));
		cudaMalloc((void **) &dvx_d, N * sizeof(double));
		cudaMalloc((void **) &dvy_d, N * sizeof(double));
		cudaMalloc((void **) &dvz_d, N * sizeof(double));
	}

	//perturbers data
	if(useGPU == 0){
		XYdata_h = (double*)malloc(Nperturbers * NTable * 4 * sizeof(double));
		readBufferA_h = (double*)malloc(Nperturbers * 4 * sizeof(double));
		readBufferB_h = (double*)malloc(Nperturbers * 4 * sizeof(double));
	}
	else{
		//allocate data for interleaved data transfer
		cudaHostAlloc((void **) &readBufferA_h, Nperturbers * 4 * sizeof(double), cudaHostAllocDefault);
		cudaHostAlloc((void **) &readBufferB_h, Nperturbers * 4 * sizeof(double), cudaHostAllocDefault);
		cudaMalloc((void **) &XYdata_d, Nperturbers * NTable * 4 * sizeof(double));
	}

}

__host__ void Host::initialize1(){
	for(int i = 0; i < N; ++i){
		A1_h[i] = 0.0;
		A2_h[i] = 0.0;
		A3_h[i] = 0.0;
		jd_init_h[i] = 0.0;
	}

	//Sun
	id_h[0] = 10;
	index_h[0] = 0;
	m_h[0] = 1.0;
	x_h[0] = 0.0;
	y_h[0] = 0.0;
	z_h[0] = 0.0;
	vx_h[0] = 0.0;
	vy_h[0] = 0.0;
	vz_h[0] = 0.0;

	for(int i = 1; i < N; ++i){
		id_h[i] = Nperturbers + i;
		index_h[i] = i;
		m_h[i] = 0.0;
		x_h[i] = 0.0;
		y_h[i] = 0.0;
		z_h[i] = 0.0;
		vx_h[i] = 0.0;
		vy_h[i] = 0.0;
		vz_h[i] = 0.0;
	}
}

__host__ void Host::copy1(){

	cudaMemcpy(m_d, m_h, N * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(id_d, id_h, N * sizeof(unsigned long long int), cudaMemcpyHostToDevice);
	cudaMemcpy(index_d, index_h, N * sizeof(unsigned int), cudaMemcpyHostToDevice);
	cudaMemcpy(x_d, x_h, N * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(y_d, y_h, N * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(z_d, z_h, N * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(vx_d, vx_h, N * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(vy_d, vy_h, N * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(vz_d, vz_h, N * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(A1_d, A1_h, N * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(A2_d, A2_h, N * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(A3_d, A3_h, N * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(jd_init_d, jd_init_h, N * sizeof(double), cudaMemcpyHostToDevice);

}

__host__ void Host::Alloc2(){
	xt_h = (double*)malloc(N * sizeof(double));
	yt_h = (double*)malloc(N * sizeof(double));
	zt_h = (double*)malloc(N * sizeof(double));
	vxt_h = (double*)malloc(N * sizeof(double));
	vyt_h = (double*)malloc(N * sizeof(double));
	vzt_h = (double*)malloc(N * sizeof(double));

	kx_h = (double*)malloc(N * RKFn * sizeof(double));
	ky_h = (double*)malloc(N * RKFn * sizeof(double));
	kz_h = (double*)malloc(N * RKFn * sizeof(double));
	kvx_h = (double*)malloc(N * RKFn * sizeof(double));
	kvy_h = (double*)malloc(N * RKFn * sizeof(double));
	kvz_h = (double*)malloc(N * RKFn * sizeof(double));

	snew_h = (double2*)malloc(N * sizeof(double2));
	dtmin_h = (double*)malloc(N * sizeof(double));

	xTable_h = (double*)malloc(Nperturbers * RKFn * sizeof(double));
	yTable_h = (double*)malloc(Nperturbers * RKFn * sizeof(double));
	zTable_h = (double*)malloc(Nperturbers * RKFn * sizeof(double));

	if(useGPU > 0){
		cudaMalloc((void **) &kx_d, N * RKFn * sizeof(double));
		cudaMalloc((void **) &ky_d, N * RKFn * sizeof(double));
		cudaMalloc((void **) &kz_d, N * RKFn * sizeof(double));
		cudaMalloc((void **) &kvx_d, N * RKFn * sizeof(double));
		cudaMalloc((void **) &kvy_d, N * RKFn * sizeof(double));
		cudaMalloc((void **) &kvz_d, N * RKFn * sizeof(double));

		cudaMalloc((void **) &snew_d, N * sizeof(double2));
		cudaMalloc((void **) &scan_d, N * sizeof(int2));
		cudaMalloc((void **) &N_d, sizeof(int));

		cudaMalloc((void **) &xTable_d, Nperturbers * RKFn * sizeof(double));
		cudaMalloc((void **) &yTable_d, Nperturbers * RKFn * sizeof(double));
		cudaMalloc((void **) &zTable_d, Nperturbers * RKFn * sizeof(double));
	}

	a_h = (double*)malloc(RKFn * RKFn * sizeof(double));
	b_h = (double*)malloc(RKFn * sizeof(double));
	bb_h = (double*)malloc(RKFn * sizeof(double));
	c_h = (double*)malloc(RKFn * sizeof(double));

}
__host__ void Host::initialize2(){
	for(int i = 0; i < N; ++i){
		dtmin_h[i] = 1.0e6;
	}

	for(int i = 0; i < RKFn; ++i){
		for(int j = 0; j < RKFn; ++j){
			a_h[i * RKFn + j] = 0.0;
		}
		b_h[i] = 0.0;
		bb_h[i] = 0.0;
		c_h[i] = 0.0;
	}
}


__host__ void Host::initialize3(){
	for(int i = 0; i < N; ++i){
		dtmin_h[i] = 1.0e6;
	}

	//save coordinates for backward integrations
	for(int i = 0; i < N; ++i){
		x0_h[i] = x_h[i];
		y0_h[i] = y_h[i];
		z0_h[i] = z_h[i];
		vx0_h[i] = vx_h[i];
		vy0_h[i] = vy_h[i];
		vz0_h[i] = vz_h[i];
		A10_h[i] = A1_h[i];
		A20_h[i] = A2_h[i];
		A30_h[i] = A3_h[i];
		m0_h[i] = m_h[i];
		id0_h[i] = id_h[i];

		xb_h[i] = x_h[i];
		yb_h[i] = y_h[i];
		zb_h[i] = z_h[i];
		vxb_h[i] = vx_h[i];
		vyb_h[i] = vy_h[i];
		vzb_h[i] = vz_h[i];
		A1b_h[i] = A1_h[i];
		A2b_h[i] = A2_h[i];
		A3b_h[i] = A3_h[i];
		mb_h[i] = m_h[i];
		idb_h[i] = id_h[i];
//if(id_h[i] == 72057594038045489) printf("S %d %llu %.20g %.20g %.20g %.20g\n", i, id_h[i], m_h[i], x0_h[i], A1_h[i], snew_h[i].y);
	}
	if(useGPU > 0){
		cudaMemcpy(m0_d, m0_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(id0_d, id0_h, N * sizeof(unsigned long long int), cudaMemcpyHostToDevice);
		cudaMemcpy(x0_d, x0_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(y0_d, y0_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(z0_d, z0_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(vx0_d, vx0_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(vy0_d, vy0_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(vz0_d, vz0_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(A10_d, A10_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(A20_d, A20_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(A30_d, A30_h, N * sizeof(double), cudaMemcpyHostToDevice);
	}
}
__host__ void Host::restore3(){

	//restore coordinates for backward integrations
	for(int i = 0; i < N; ++i){
		dtmin_h[i] = 1.0e6;

		x_h[i] = xb_h[i];
		y_h[i] = yb_h[i];
		z_h[i] = zb_h[i];
		vx_h[i] = vxb_h[i];
		vy_h[i] = vyb_h[i];
		vz_h[i] = vzb_h[i];
		A1_h[i] = A1b_h[i];
		A2_h[i] = A2b_h[i];
		A3_h[i] = A3b_h[i];
		m_h[i] = mb_h[i];
		id_h[i] = idb_h[i];
		index_h[i] = i;

		x0_h[i] = x_h[i];
		y0_h[i] = y_h[i];
		z0_h[i] = z_h[i];
		vx0_h[i] = vx_h[i];
		vy0_h[i] = vy_h[i];
		vz0_h[i] = vz_h[i];
		A10_h[i] = A1_h[i];
		A20_h[i] = A2_h[i];
		A30_h[i] = A3_h[i];
		m0_h[i] = m_h[i];
		id0_h[i] = id_h[i];
		index0_h[i] = i;

//if(id_h[i] == 72057594038045489) printf("R %d %llu %.20g %.20g %.20g %.20g\n", i, id_h[i], m_h[i], x0_h[i], A1_h[i], snew_h[i].y);
	}
	if(useGPU > 0){
		cudaMemcpy(m_d, m_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(id_d, id_h, N * sizeof(unsigned long long int), cudaMemcpyHostToDevice);
		cudaMemcpy(index_d, index_h, N * sizeof(unsigned int), cudaMemcpyHostToDevice);
		cudaMemcpy(x_d, x_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(y_d, y_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(z_d, z_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(vx_d, vx_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(vy_d, vy_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(vz_d, vz_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(A1_d, A1_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(A2_d, A2_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(A3_d, A3_h, N * sizeof(double), cudaMemcpyHostToDevice);

		cudaMemcpy(m0_d, m0_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(id0_d, id0_h, N * sizeof(unsigned long long int), cudaMemcpyHostToDevice);
		cudaMemcpy(index0_d, index0_h, N * sizeof(unsigned int), cudaMemcpyHostToDevice);
		cudaMemcpy(x0_d, x0_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(y0_d, y0_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(z0_d, z0_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(vx0_d, vx0_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(vy0_d, vy0_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(vz0_d, vz0_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(A10_d, A10_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(A20_d, A20_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(A30_d, A30_h, N * sizeof(double), cudaMemcpyHostToDevice);
	}
}

__global__ void setSnew_kernel(double2 *snew_d, int N){

	int id = blockIdx.x * blockDim.x + threadIdx.x;
	if(id < N){
		snew_d[id].x = 1.5;
		snew_d[id].y = 1.5;
	}
}

__host__ void Host::setSnew(){

	if(useGPU == 0){
		for(int i = 0; i < N; ++i){
			snew_h[i].x = 1.5;
			snew_h[i].y = 1.5;
		}
	}
	else{
		setSnew_kernel <<< (N + 255) / 256, 256 >>> (snew_d, N);
	}
}

__host__ void Host::reduce(int S){
	//reduce arrays for repeated integration with a smaller time step
	int k = Nperturbers;

	for(int i = Nperturbers; i < N; ++i){
		index0_h[i] = index_h[i]; 
	}
	
	for(int i = Nperturbers; i < N; ++i){
		if(snew_h[i].y < 1.0){
			int ii = index0_h[i];

			x_h[k] = x0_h[ii];
			y_h[k] = y0_h[ii];
			z_h[k] = z0_h[ii];
			vx_h[k] = vx0_h[ii];
			vy_h[k] = vy0_h[ii];
			vz_h[k] = vz0_h[ii];
			A1_h[k] = A10_h[ii];
			A2_h[k] = A20_h[ii];
			A3_h[k] = A30_h[ii];
			m_h[k] = m0_h[ii];
			id_h[k] = id0_h[ii];
			index_h[k] = index0_h[ii];

//printf("%d %d %u %llu\n", i, k, ii, id_h[k]);
printf("%d %d %u %llu %.20g %.20g\n", i, k, ii, id_h[k], x_h[k], A1_h[k]);
			++k; 
		}
	}
	N = k;

	runsN[S + 1] = N;
}

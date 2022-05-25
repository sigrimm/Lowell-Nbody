#include "Host.h"

//initialize variables
__host__ Host::Host(){

	timep0 = 0.0;

	comet = 0;

	useGR = 2;
	useJ2 = 1;
	useNonGrav = 1;

	useGPU = 1;

	useFIFO = 2;
	InVersion= 0.0;
	DoPreIntegration = 0;	

	useHelio = 1;
	outHelio = 1;

	outBinary = 1;

	Nsteps = 1e9;
	outInterval = 1;


	time0 = 0.0;
	time1 = 0.0;
	outStart = 0.0;
	time = 0.0;

	dti = 0.05;
	dt = dti * dayUnit;
	dtiMin = 1.0e-5;

	dts = 0.01;
	outI = 1llu;

	N = Nperturbers;

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

__host__ void Host::Alloc1(){
	//allocate data on host
	id_h = (unsigned long long int*)malloc(N * sizeof(unsigned long long int));
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

	dx_h = (double*)malloc(N * sizeof(double));
	dy_h = (double*)malloc(N * sizeof(double));
	dz_h = (double*)malloc(N * sizeof(double));
	dvx_h = (double*)malloc(N * sizeof(double));
	dvy_h = (double*)malloc(N * sizeof(double));
	dvz_h = (double*)malloc(N * sizeof(double));

	if(useGPU > 0){
		//allocate data on the device
		cudaMalloc((void **) &id_d, N * sizeof(unsigned long long int));
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

__host__ void Host::Initialize1(){
	for(int i = 0; i < N; ++i){
		A1_h[i] = 0.0;
		A2_h[i] = 0.0;
		A3_h[i] = 0.0;
		jd_init_h[i] = 0.0;
	}

	//Sun
	id_h[0] = 10;
	m_h[0] = 1.0;
	x_h[0] = 0.0;
	y_h[0] = 0.0;
	z_h[0] = 0.0;
	vx_h[0] = 0.0;
	vy_h[0] = 0.0;
	vz_h[0] = 0.0;

	for(int i = 1; i < N; ++i){
		id_h[i] = Nperturbers + i;
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

		cudaMalloc((void **) &xTable_d, Nperturbers * RKFn * sizeof(double));
		cudaMalloc((void **) &yTable_d, Nperturbers * RKFn * sizeof(double));
		cudaMalloc((void **) &zTable_d, Nperturbers * RKFn * sizeof(double));
	}

	a_h = (double*)malloc(RKFn * RKFn * sizeof(double));
	b_h = (double*)malloc(RKFn * sizeof(double));
	bb_h = (double*)malloc(RKFn * sizeof(double));
	c_h = (double*)malloc(RKFn * sizeof(double));

}
__host__ void Host::Initialize2(){
	for(int i = 0; i < N; ++i){
		dtmin_h[i] = 1.0e6;
		snew_h[i].x = 1.5;
		snew_h[i].y = 1.5;
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


__host__ void Host::Initialize3(){
	for(int i = 0; i < N; ++i){
		dtmin_h[i] = 1.0e6;
		snew_h[i].x = 1.5;
		snew_h[i].y = 1.5;
	}
	if(useGPU > 0){
		cudaMemcpy(snew_d, snew_h, N * sizeof(double2), cudaMemcpyHostToDevice);
	}

	//save coordinates for repeated integrations
	for(int i = 0; i < N; ++i){
		x0_h[i] = x_h[i];
		y0_h[i] = y_h[i];
		z0_h[i] = z_h[i];
		vx0_h[i] = vx_h[i];
		vy0_h[i] = vy_h[i];
		vz0_h[i] = vz_h[i];
	}

}


__host__ void Host::reduce(int S){
	//reduce arrays for repeated integration with a smaller time step
	int k = Nperturbers;
	for(int i = Nperturbers; i < N; ++i){
		if(snew_h[i].y < 1.0){
			m_h[k] = m_h[i];
			id_h[k] = id_h[i];
			x_h[k] = x0_h[i];
			y_h[k] = y0_h[i];
			z_h[k] = z0_h[i];
			vx_h[k] = vx0_h[i];
			vy_h[k] = vy0_h[i];
			vz_h[k] = vz0_h[i];
			A1_h[k] = A1_h[i];
			A2_h[k] = A2_h[i];
			A3_h[k] = A3_h[i];

			x0_h[k] = x0_h[i];
			y0_h[k] = y0_h[i];
			z0_h[k] = z0_h[i];
			vx0_h[k] = vx0_h[i];
			vy0_h[k] = vy0_h[i];
			vz0_h[k] = vz0_h[i];
			snew_h[k].x = 1.5;
			snew_h[k].y = 1.5;
printf("%d %d %llu\n", i, k, id_h[k]);
			++k; 
		}
	}
	N = k;
	Sn[S + 1] = N;
	if(useGPU > 0){
		cudaMemcpy(m_d, m_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(id_d, id_h, N * sizeof(unsigned long long int), cudaMemcpyHostToDevice);
		cudaMemcpy(x_d, x_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(y_d, y_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(z_d, z_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(vx_d, vx_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(vy_d, vy_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(vz_d, vz_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(A1_d, A1_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(A2_d, A2_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(A3_d, A3_h, N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(snew_d, snew_h, N * sizeof(double2), cudaMemcpyHostToDevice);
	}
	if(S == 0){
		dti = 1.0;
		dts = 0.01;
		dtiMin = 0.1;
		dt = dti * dayUnit;
		outI = (outInterval + 0.5 * dts) / dts;
	}
	if(S == 1){
		dti = 0.01;
		dts = 1.0e-5;
		dtiMin = 1.0e-4;
		dt = dti * dayUnit;
		outI = (outInterval + 0.5 * dts) / dts;
	}

}

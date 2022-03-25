void output(double *x_d, double *y_d, double *z_d, double *vx_d, double *vy_d, double *vz_d, double *x_h, double *y_h, double *z_h, double *vx_h, double *vy_h, double *vz_h, double *m_h, unsigned long long int *id_h, long long int t, double time, int N, int Nperturbers, int NN, int useGPU, int useHelio, int outHelio, int outBinary){

	printf("Output %.20g %lld\n", time, t);

	if(useGPU == 1){
		cudaMemcpy(x_h, x_d, NN * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(y_h, y_d, NN * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(z_h, z_d, NN * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(vx_h, vx_d, NN * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(vy_h, vy_d, NN * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(vz_h, vz_d, NN * sizeof(double), cudaMemcpyDeviceToHost);
	}

	FILE *outfile;
	char outfilename[160];
	
	if(outHelio == 1){
		if(outBinary == 0){	
			sprintf(outfilename, "Outhelio10_%.12lld.dat", t);
		}
		else{
			sprintf(outfilename, "Outhelio10.bin");
		}
	}
	else{
		if(outBinary == 0){
			sprintf(outfilename, "Outbary10_%.12lld.dat", t);
		}
		else{
			sprintf(outfilename, "Outbary10.bin");
		}
	}
	if(outBinary == 0){
		outfile = fopen(outfilename, "w");
	}
	else{
		if(t == 0){
			outfile = fopen(outfilename, "wb");
		}
		else{
			outfile = fopen(outfilename, "ab");
		}
	}

//	printf("%s\n", outfilename);

	double comx = 0.0;
	double comy = 0.0;
	double comz = 0.0;
	double vcomx = 0.0;
	double vcomy = 0.0;
	double vcomz = 0.0;

	if(useHelio == 0 && outHelio == 1){
		//convert to heliocentric output
		comx = -x_h[0];
		comy = -y_h[0];
		comz = -z_h[0];
		vcomx = -vx_h[0];
		vcomy = -vy_h[0];
		vcomz = -vz_h[0];
	}
	
	if(outBinary == 0){
		for(int i = Nperturbers; i < N; ++i){
			fprintf(outfile, "%.10g %llu %.40g %.40g %.40g %.40g %.40g %.40g %.40g\n", time, id_h[i], m_h[i], comx + x_h[i], comy + y_h[i], comz + z_h[i], (vcomx + vx_h[i]) * dayUnit, (vcomy + vy_h[i]) * dayUnit, (vcomz + vz_h[i]) * dayUnit/*, dtmin_h[i]*/);
		}
	}
	else{
		for(int i = Nperturbers; i < N; ++i){

			//unsigned long long int id = id_h[i];
			unsigned long long int id = __builtin_bswap64 (id_h[i]);
			double xx = comx + x_h[i];
			double yy = comy + y_h[i];
			double zz = comz + z_h[i];
			double vxx = (vcomx + vx_h[i]) * dayUnit;
			double vyy = (vcomy + vy_h[i]) * dayUnit;
			double vzz = (vcomz + vz_h[i]) * dayUnit;

			fwrite(&id, 1, sizeof(unsigned long long int), outfile);
			fwrite(&time, 1, sizeof(double), outfile);
			fwrite(&xx, 1, sizeof(double), outfile);
			fwrite(&yy, 1, sizeof(double), outfile);
			fwrite(&zz, 1, sizeof(double), outfile);
			fwrite(&vxx, 1, sizeof(double), outfile);
			fwrite(&vyy, 1, sizeof(double), outfile);
			fwrite(&vzz, 1, sizeof(double), outfile);

			//printf("%llu %g %g %g %g %g %g %g\n", id, time, xx, yy, zz, vxx, vyy, vzz);
		}
		
	}
	fclose(outfile);

}


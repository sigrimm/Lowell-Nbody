void output(double2 *snew_h, double *dtmin_h, double *x_h, double *y_h, double *z_h, double *vx_h, double *vy_h, double *vz_h, double *xt_h, double *yt_h, double *zt_h, double *m_h, unsigned long long int *id_h, long long int t, double time, int N, int Nperturbers, int useGPU, int useHelio, int outHelio, int outBinary, int S){

	printf("Output %d %.20g %lld\n", S, time, t);


	FILE *outfile;
	char outfilename[160];
	
	if(outHelio == 1){
		if(outBinary == 0){	
			sprintf(outfilename, "Outhelio_%d_%.12lld.dat", S, t);
		}
		else{
			sprintf(outfilename, "Outhelio_%d.bin", S);
		}
	}
	else{
		if(outBinary == 0){
			sprintf(outfilename, "Outbary_%d_%.12lld.dat", S, t);
		}
		else{
			sprintf(outfilename, "Outbary_%d.bin", S);
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
		//for(int i = 0; i < Nperturbers; ++i){
		//	fprintf(outfile, "%.10g %llu %.40g %.40g %.40g %.40g %.40g %.40g %.40g %g\n", time, id_h[i], m_h[i], xt_h[i], yt_h[i], zt_h[i], 0.0, 0.0, 0.0, 0.0);
		//}
		for(int i = Nperturbers; i < N; ++i){
			//fprintf(outfile, "%.10g %llu %.40g %.40g %.40g %.40g %.40g %.40g %.40g %g\n", time, id_h[i], m_h[i], comx + x_h[i], comy + y_h[i], comz + z_h[i], (vcomx + vx_h[i]) * dayUnit, (vcomy + vy_h[i]) * dayUnit, (vcomz + vz_h[i]) * dayUnit, dtmin_h[i]);
//			if(snew_h[i].y >= 1.0){
				fprintf(outfile, "%.10g %llu %.40g %.40g %.40g %.40g %.40g %.40g %.40g %g\n", time, id_h[i], m_h[i], comx + x_h[i], comy + y_h[i], comz + z_h[i], (vcomx + vx_h[i]) * dayUnit, (vcomy + vy_h[i]) * dayUnit, (vcomz + vz_h[i]) * dayUnit, snew_h[i].y);
//			}
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


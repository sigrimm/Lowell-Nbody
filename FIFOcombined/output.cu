#include "Host.h"
__host__ void Host::output(long long int t, double time){

	printf("Output %.20g %lld\n", time, t);


	FILE *outfile;
	char outfilename[160];
	
	if(outHelio == 1){
		if(outBinary == 0){	
			sprintf(outfilename, "Outhelio_%.12lld.dat", t);
		}
		else{
			sprintf(outfilename, "Outhelio.bin");
		}
	}
	else{
		if(outBinary == 0){
			sprintf(outfilename, "Outbary_%.12lld.dat", t);
		}
		else{
			sprintf(outfilename, "Outbary.bin");
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


	if(useGPU > 0 && t > 0){
		cudaMemcpy(snew_h, snew_d, N * sizeof(double2), cudaMemcpyDeviceToHost);
		cudaMemcpy(x0_h, x0_d, N * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(y0_h, y0_d, N * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(z0_h, z0_d, N * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(vx0_h, vx0_d, N * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(vy0_h, vy0_d, N * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(vz0_h, vz0_d, N * sizeof(double), cudaMemcpyDeviceToHost);

		cudaMemcpy(xTable_h, xTable_d, Nperturbers * RKFn * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(yTable_h, yTable_d, Nperturbers * RKFn * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(zTable_h, zTable_d, Nperturbers * RKFn * sizeof(double), cudaMemcpyDeviceToHost);
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
		//for(int p = 0; p < Nperturbers; ++p){
		//	int ii = p * RKFn + 12;	
		//	fprintf(outfile, "%.10g %llu %.40g %.40g %.40g %.40g %.40g %.40g %.40g %g\n", time, id_h[p], m_h[p], xTable_h[ii], yTable_h[ii], zTable_h[ii], 0.0, 0.0, 0.0, 0.0);
		//}
		for(int i = Nperturbers; i < N; ++i){

			fprintf(outfile, "%.10g %llu %.40g %.40g %.40g %.40g %.40g %.40g %.40g %g\n", time, id0_h[i], m0_h[i], comx + x0_h[i], comy + y0_h[i], comz + z0_h[i], (vcomx + vx0_h[i]) * dayUnit, (vcomy + vy0_h[i]) * dayUnit, (vcomz + vz0_h[i]) * dayUnit, dtmin_h[i]);

			//if(snew_h[i].y >= 1.0){
			//	fprintf(outfile, "%.10g %llu %.40g %.40g %.40g %.40g %.40g %.40g %.40g %g\n", time, id_h[i], m_h[i], comx + x_h[i], comy + y_h[i], comz + z_h[i], (vcomx + vx_h[i]) * dayUnit, (vcomy + vy_h[i]) * dayUnit, (vcomz + vz_h[i]) * dayUnit, snew_h[i].y);
			//}
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


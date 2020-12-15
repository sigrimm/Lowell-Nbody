#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define dayUnit 0.01720209895
//#define dayUnit 0.01720209894846
//#define dayUnit 1.0
//Erics value:
//#define dayUnit 0.01720412578565452474399499749324604636058


void acc(double *m, double *x, double *y, double *z, double &ax, double &ay, double &az, int i, int j){
	double rx = x[j] - x[i];
	double ry = y[j] - y[i];
	double rz = z[j] - z[i];
	double rsq = rx * rx + ry * ry + rz * rz;
	double r = sqrt(rsq);
	
	double k = 1.0;
	//use tis if dayUnit = 1.0
	//double k = 0.01720209895 * 0.01720209895;
	//double k = 0.000295981944048623;

	ax += k * m[j] / (r * rsq) * rx;
	ay += k * m[j] / (r * rsq) * ry;
	az += k * m[j] / (r * rsq) * rz;
}


int main(){

	//Number of planets
	
	int N = 2;

	int outHelio = 1;
	//1 print output in heliocentric coordinates
	//0 print output in barycentric  coordinates

	long long int Nsteps = 4000;	
	long long int outInterval = 1;
	double dt = 1.0 * dayUnit;
	
	double m[N];
	double rad[N];
	double x[N];
	double y[N];
	double z[N];
	double vx[N];
	double vy[N];
	double vz[N];

	double Spinx[N];
	double Spiny[N];
	double Spinz[N];

	double love[N];
	double tau[N];

	FILE *outfile;
	char outfilename[160];	

	FILE *infile;
	char infilename[160];

	sprintf(infilename, "initial.dat");
	infile = fopen(infilename, "r");

	double time = 0.0;

	for(int i = 1; i < N; ++i){
		int skip;
		fscanf(infile, "%lf", &time);
		fscanf(infile, "%lf", &m[i]);
		fscanf(infile, "%lf", &x[i]);
		fscanf(infile, "%lf", &y[i]);
		fscanf(infile, "%lf", &z[i]);
		fscanf(infile, "%lf", &vx[i]);
		fscanf(infile, "%lf", &vy[i]);
		fscanf(infile, "%lf", &vz[i]);

		vx[i] /= dayUnit;
		vy[i] /= dayUnit;
		vz[i] /= dayUnit;
	}

	//Sun
	m[0] = 1.0;
	rad[0] = 1.0e-3; //?
	x[0] = 0.0;
	y[0] = 0.0;
	z[0] = 0.0;
	vx[0] = 0.0;
	vy[0] = 0.0;
	vz[0] = 0.0;


	fclose(infile);
	
	
	
	//barycentric, all particles
	
	//convert to barycentric coordinates
	double comx = 0.0;
	double comy = 0.0;
	double comz = 0.0;
	double vcomx = 0.0;
	double vcomy = 0.0;
	double vcomz = 0.0;
	double mtot = 0.0;
	
	for(int i = 0; i < N; ++i){
		comx += m[i] * x[i];
		comy += m[i] * y[i];
		comz += m[i] * z[i];
		vcomx += m[i] * vx[i];
		vcomy += m[i] * vy[i];
		vcomz += m[i] * vz[i];
		mtot += m[i];
	}

	for(int i = 0; i < N; ++i){
		x[i] -= comx / mtot;
		y[i] -= comy / mtot;
		z[i] -= comz / mtot;
		vx[i] -= vcomx / mtot;
		vy[i] -= vcomy / mtot;
		vz[i] -= vcomz / mtot;
	}

	//first output
	
	//convert to heliocentric output
	comx = 0.0;
	comy = 0.0;
	comz = 0.0;
	vcomx = 0.0;
	vcomy = 0.0;
	vcomz = 0.0;
	
	if(outHelio == 1){
		for(int i = 1; i < N; ++i){
			comx += m[i] * x[i];
			comy += m[i] * y[i];
			comz += m[i] * z[i];
			vcomx += m[i] * vx[i];
			vcomy += m[i] * vy[i];
			vcomz += m[i] * vz[i];
		}
		comx /= m[0];
		comy /= m[0];
		comz /= m[0];
		vcomx /= m[0];
		vcomy /= m[0];
		vcomz /= m[0];
	}
	
	if(outHelio == 1){
		sprintf(outfilename, "Outhelio10_%.12d.dat", 0);
	}
	else{
		sprintf(outfilename, "Outbary10_%.12d.dat", 0);
	}
	outfile = fopen(outfilename, "w");
	printf("%s\n", outfilename);
	for(int i = 0; i < N; ++i){
		fprintf(outfile, "%.10g %d %.40g %.40g %.40g %.40g %.40g %.40g %.40g %.40g 0 0 0 0 0 0 0 0 0 0 0\n", time, i, m[i], rad[i], comx + x[i], comy + y[i], comz + z[i], vcomx + vx[i], vcomy + vy[i], vcomz + vz[i]);

	}
	fclose(outfile);

	double kx[N][13];
	double ky[N][13];
	double kz[N][13];
	double kvx[N][13];
	double kvy[N][13];
	double kvz[N][13];

	double xt[N];
	double yt[N];
	double zt[N];
	double vxt[N];
	double vyt[N];
	double vzt[N];

	double errorkx[N];
	double errorky[N];
	double errorkz[N];
	double errorkvx[N];
	double errorkvy[N];
	double errorkvz[N];

	double a21 = 2.0/27.0;

	double a31 = 1.0/36.0;
	double a32 = 1.0/12.0;

	double a41 = 1.0/24.0;
	double a43 = 1.0/8.0;

	double a51 = 5.0/12.0;
	double a53 = -25.0/16.0;
	double a54 = 25.0/16.0;

	double a61 = 1.0/20.0;
	double a64 = 1.0/4.0;
	double a65 = 1.0/5.0;

	double a71 = -25.0/ 108.0;
	double a74 =  125.0/ 108.0;
	double a75 = -65.0/ 27.0;
	double a76 =  125.0/ 54.0;

	double a81 = 31.0/300.0;
	double a85 = 61.0/225.0;
	double a86 = -2.0/9.0;
	double a87 = 13.0/900.0;

	double a91 = 2.0;
	double a94 = -53.0/6.0;
	double a95 = 704.0/45.0;
	double a96 = -107.0/9.0;
	double a97 = 67.0/90.0;
	double a98 = 3.0;

	double a101 = -91.0/108.0;
	double a104 = 23.0/108.0;
	double a105 = -976.0/135.0;
	double a106 = 311.0/54.0;
	double a107 = -19.0/60.0;
	double a108 = 17.0/6.0;
	double a109 = -1.0/12.0;

	double a111 = 2383.0/4100.0;
	double a114 = -341.0/164.0;
	double a115 = 4496.0/1025.0;
	double a116 = -301.0/82.0;
	double a117 = 2133.0/4100.0;
	double a118 = 45.0/82.0;
	double a119 = 45.0/164.0;
	double a1110 = 18.0/41.0;

	double a121 = 3.0/205.0;
	double a126 = - 6.0/41.0;
	double a127 = - 3.0/205.0;
	double a128 = - 3.0/41.0;
	double a129 = 3.0/41.0;
	double a1210 = 6.0/41.0;


	//double a131 = -1777.0/4100.0, a134 = -341.0/164.0, a135 = 4496.0/1025.0, a136 = -289.0/82.0, a137 = 2193.0/4100.0, a138 = 51.0/82.0, a139 = 33.0/164.0, a1310 = 12.0/41.0;
	double a131 = -1777.0/4100.0;
	double a134 = -341.0/164.0;
	double a135 = 4496.0/1025.0;
	double a136 = -289.0/82.0;
	double a137 = 2193.0/4100.0;
	double a138 = 51.0/82.0;
	double a139 = 33.0/164.0;
	double a1310 = 19.0/41.0;
	//a1310 = 19.0 / 41.0

	double b1 = 41.0/840.0;
	double b6 = 34.0/105.0;
	double b7 = 9.0/35.0;
	double b8 = 9.0/35.0;
	double b9 = 9.0/280.0;
	double b10 = 9.0/280.0;
	double b11 =41.0/840.0;

	double er = -41.0/840.0;


	for(long long int t = 1; t < Nsteps; ++t){
		double ax = 0.0;
		double ay = 0.0;
		double az = 0.0;
		
		//stage 1	
		for(int i = 0; i < N; ++i){
			xt[i] = x[i];
			yt[i] = y[i];
			zt[i] = z[i];
			vxt[i] = vx[i];
			vyt[i] = vy[i];
			vzt[i] = vz[i];
		}
		for(int i = 0; i < N; ++i){
			kx[i][0] = vxt[i];
			ky[i][0] = vyt[i];
			kz[i][0] = vzt[i];

			ax = 0.0;
			ay = 0.0;
			az = 0.0;

			for(int j = 0; j < N; ++j){
				if(i != j){
					acc(m, xt, yt, zt, ax, ay, az, i, j);
				}
			}

			kvx[i][0] = ax;
			kvy[i][0] = ay;
			kvz[i][0] = az;
		}

		
		//stage 2
		for(int i = 0; i < N; ++i){
			xt[i]  = x[i]  + dt * a21 * kx[i][0];
			yt[i]  = y[i]  + dt * a21 * ky[i][0];
			zt[i]  = z[i]  + dt * a21 * kz[i][0];
			vxt[i] = vx[i] + dt * a21 * kvx[i][0];
			vyt[i] = vy[i] + dt * a21 * kvy[i][0];
			vzt[i] = vz[i] + dt * a21 * kvz[i][0];
		}

		for(int i = 0; i < N; ++i){

			kx[i][1] = vxt[i];
			ky[i][1] = vyt[i];
			kz[i][1] = vzt[i];

			ax = 0.0;
			ay = 0.0;
			az = 0.0;

			for(int j = 0; j < N; ++j){
				if(i != j){
					acc(m, xt, yt, zt, ax, ay, az, i, j);
				}
			}

			kvx[i][1] = ax;
			kvy[i][1] = ay;
			kvz[i][1] = az;

		}
		//stage 3
		for(int i = 0; i < N; ++i){
			xt[i]  = x[i]  + dt * (a31 * kx[i][0]  + a32 * kx[i][1]);
			yt[i]  = y[i]  + dt * (a31 * ky[i][0]  + a32 * ky[i][1]);
			zt[i]  = z[i]  + dt * (a31 * kz[i][0]  + a32 * kz[i][1]);
			vxt[i] = vx[i] + dt * (a31 * kvx[i][0] + a32 * kvx[i][1]);
			vyt[i] = vy[i] + dt * (a31 * kvy[i][0] + a32 * kvy[i][1]);
			vzt[i] = vz[i] + dt * (a31 * kvz[i][0] + a32 * kvz[i][1]);
		}

		for(int i = 0; i < N; ++i){

			kx[i][2] = vxt[i];
			ky[i][2] = vyt[i];
			kz[i][2] = vzt[i];

			ax = 0.0;
			ay = 0.0;
			az = 0.0;

			for(int j = 0; j < N; ++j){
				if(i != j){
					acc(m, xt, yt, zt, ax, ay, az, i, j);
				}
			}

			kvx[i][2] = ax;
			kvy[i][2] = ay;
			kvz[i][2] = az;

		}
		//stage 4
		for(int i = 0; i < N; ++i){
			xt[i]  = x[i]  + dt * (a41 * kx[i][0]  + a43 * kx[i][2]);
			yt[i]  = y[i]  + dt * (a41 * ky[i][0]  + a43 * ky[i][2]);
			zt[i]  = z[i]  + dt * (a41 * kz[i][0]  + a43 * kz[i][2]);
			vxt[i] = vx[i] + dt * (a41 * kvx[i][0] + a43 * kvx[i][2]);
			vyt[i] = vy[i] + dt * (a41 * kvy[i][0] + a43 * kvy[i][2]);
			vzt[i] = vz[i] + dt * (a41 * kvz[i][0] + a43 * kvz[i][2]);
		}

		for(int i = 0; i < N; ++i){

			kx[i][3] = vxt[i];
			ky[i][3] = vyt[i];
			kz[i][3] = vzt[i];

			ax = 0.0;
			ay = 0.0;
			az = 0.0;

			for(int j = 0; j < N; ++j){
				if(i != j){
					acc(m, xt, yt, zt, ax, ay, az, i, j);
				}
			}

			kvx[i][3] = ax;
			kvy[i][3] = ay;
			kvz[i][3] = az;

		}
		//stage 5
		for(int i = 0; i < N; ++i){
			xt[i]  = x[i]  + dt * (a51 * kx[i][0]  + a53 * kx[i][2]  + a54 * kx[i][3]);
			yt[i]  = y[i]  + dt * (a51 * ky[i][0]  + a53 * ky[i][2]  + a54 * ky[i][3]);
			zt[i]  = z[i]  + dt * (a51 * kz[i][0]  + a53 * kz[i][2]  + a54 * kz[i][3]);
			vxt[i] = vx[i] + dt * (a51 * kvx[i][0] + a53 * kvx[i][2] + a54 * kvx[i][3]);
			vyt[i] = vy[i] + dt * (a51 * kvy[i][0] + a53 * kvy[i][2] + a54 * kvy[i][3]);
			vzt[i] = vz[i] + dt * (a51 * kvz[i][0] + a53 * kvz[i][2] + a54 * kvz[i][3]);
		}

		for(int i = 0; i < N; ++i){

			kx[i][4] = vxt[i];
			ky[i][4] = vyt[i];
			kz[i][4] = vzt[i];

			ax = 0.0;
			ay = 0.0;
			az = 0.0;

			for(int j = 0; j < N; ++j){
				if(i != j){
					acc(m, xt, yt, zt, ax, ay, az, i, j);
				}
			}

			kvx[i][4] = ax;
			kvy[i][4] = ay;
			kvz[i][4] = az;

		}
		//stage 6
		for(int i = 0; i < N; ++i){
			xt[i]  = x[i]  + dt * (a61 * kx[i][0]  + a64 * kx[i][3]  + a65 * kx[i][4]);
			yt[i]  = y[i]  + dt * (a61 * ky[i][0]  + a64 * ky[i][3]  + a65 * ky[i][4]);
			zt[i]  = z[i]  + dt * (a61 * kz[i][0]  + a64 * kz[i][3]  + a65 * kz[i][4]);
			vxt[i] = vx[i] + dt * (a61 * kvx[i][0] + a64 * kvx[i][3] + a65 * kvx[i][4]);
			vyt[i] = vy[i] + dt * (a61 * kvy[i][0] + a64 * kvy[i][3] + a65 * kvy[i][4]);
			vzt[i] = vz[i] + dt * (a61 * kvz[i][0] + a64 * kvz[i][3] + a65 * kvz[i][4]);
		}

		for(int i = 0; i < N; ++i){

			kx[i][5] = vxt[i];
			ky[i][5] = vyt[i];
			kz[i][5] = vzt[i];

			ax = 0.0;
			ay = 0.0;
			az = 0.0;

			for(int j = 0; j < N; ++j){
				if(i != j){
					acc(m, xt, yt, zt, ax, ay, az, i, j);
				}
			}

			kvx[i][5] = ax;
			kvy[i][5] = ay;
			kvz[i][5] = az;

		}
		//stage 7
		for(int i = 0; i < N; ++i){
			xt[i]  = x[i]  + dt * (a71 * kx[i][0]  + a74 * kx[i][3]  + a75 * kx[i][4]  + a76 * kx[i][5]);
			yt[i]  = y[i]  + dt * (a71 * ky[i][0]  + a74 * ky[i][3]  + a75 * ky[i][4]  + a76 * ky[i][5]);
			zt[i]  = z[i]  + dt * (a71 * kz[i][0]  + a74 * kz[i][3]  + a75 * kz[i][4]  + a76 * kz[i][5]);
			vxt[i] = vx[i] + dt * (a71 * kvx[i][0] + a74 * kvx[i][3] + a75 * kvx[i][4] + a76 * kvx[i][5]);
			vyt[i] = vy[i] + dt * (a71 * kvy[i][0] + a74 * kvy[i][3] + a75 * kvy[i][4] + a76 * kvy[i][5]);
			vzt[i] = vz[i] + dt * (a71 * kvz[i][0] + a74 * kvz[i][3] + a75 * kvz[i][4] + a76 * kvz[i][5]);
		}

		for(int i = 0; i < N; ++i){

			kx[i][6] = vxt[i];
			ky[i][6] = vyt[i];
			kz[i][6] = vzt[i];

			ax = 0.0;
			ay = 0.0;
			az = 0.0;

			for(int j = 0; j < N; ++j){
				if(i != j){
					acc(m, xt, yt, zt, ax, ay, az, i, j);
				}
			}

			kvx[i][6] = ax;
			kvy[i][6] = ay;
			kvz[i][6] = az;

		}
		//stage 8
		for(int i = 0; i < N; ++i){
			xt[i]  = x[i]  + dt * (a81 * kx[i][0]  + a85 * kx[i][4]  + a86 * kx[i][5]  + a87 * kx[i][6]);
			yt[i]  = y[i]  + dt * (a81 * ky[i][0]  + a85 * ky[i][4]  + a86 * ky[i][5]  + a87 * ky[i][6]);
			zt[i]  = z[i]  + dt * (a81 * kz[i][0]  + a85 * kz[i][4]  + a86 * kz[i][5]  + a87 * kz[i][6]);
			vxt[i] = vx[i] + dt * (a81 * kvx[i][0] + a85 * kvx[i][4] + a86 * kvx[i][5] + a87 * kvx[i][6]);
			vyt[i] = vy[i] + dt * (a81 * kvy[i][0] + a85 * kvy[i][4] + a86 * kvy[i][5] + a87 * kvy[i][6]);
			vzt[i] = vz[i] + dt * (a81 * kvz[i][0] + a85 * kvz[i][4] + a86 * kvz[i][5] + a87 * kvz[i][6]);
		}

		for(int i = 0; i < N; ++i){

			kx[i][7] = vxt[i];
			ky[i][7] = vyt[i];
			kz[i][7] = vzt[i];

			ax = 0.0;
			ay = 0.0;
			az = 0.0;

			for(int j = 0; j < N; ++j){
				if(i != j){
					acc(m, xt, yt, zt, ax, ay, az, i, j);
				}
			}

			kvx[i][7] = ax;
			kvy[i][7] = ay;
			kvz[i][7] = az;

		}
		//stage 9
		for(int i = 0; i < N; ++i){
			xt[i]  = x[i]  + dt * (a91 * kx[i][0]  + a94 * kx[i][3]  + a95 * kx[i][4]  + a96 * kx[i][5]  + a97 * kx[i][6]  + a98 * kx[i][7]);
			yt[i]  = y[i]  + dt * (a91 * ky[i][0]  + a94 * ky[i][3]  + a95 * ky[i][4]  + a96 * ky[i][5]  + a97 * ky[i][6]  + a98 * ky[i][7]);
			zt[i]  = z[i]  + dt * (a91 * kz[i][0]  + a94 * kz[i][3]  + a95 * kz[i][4]  + a96 * kz[i][5]  + a97 * kz[i][6]  + a98 * kz[i][7]);
			vxt[i] = vx[i] + dt * (a91 * kvx[i][0] + a94 * kvx[i][3] + a95 * kvx[i][4] + a96 * kvx[i][5] + a97 * kvx[i][6] + a98 * kvx[i][7]);
			vyt[i] = vy[i] + dt * (a91 * kvy[i][0] + a94 * kvy[i][3] + a95 * kvy[i][4] + a96 * kvy[i][5] + a97 * kvy[i][6] + a98 * kvy[i][7]);
			vzt[i] = vz[i] + dt * (a91 * kvz[i][0] + a94 * kvz[i][3] + a95 * kvz[i][4] + a96 * kvz[i][5] + a97 * kvz[i][6] + a98 * kvz[i][7]);
		}

		for(int i = 0; i < N; ++i){

			kx[i][8] = vxt[i];
			ky[i][8] = vyt[i];
			kz[i][8] = vzt[i];

			ax = 0.0;
			ay = 0.0;
			az = 0.0;

			for(int j = 0; j < N; ++j){
				if(i != j){
					acc(m, xt, yt, zt, ax, ay, az, i, j);
				}
			}

			kvx[i][8] = ax;
			kvy[i][8] = ay;
			kvz[i][8] = az;

		}
		//stage 10
		for(int i = 0; i < N; ++i){
			xt[i]  = x[i]  + dt * (a101 * kx[i][0]  + a104 * kx[i][3]  + a105 * kx[i][4]  + a106 * kx[i][5]  + a107 * kx[i][6]  + a108 * kx[i][7]);
			yt[i]  = y[i]  + dt * (a101 * ky[i][0]  + a104 * ky[i][3]  + a105 * ky[i][4]  + a106 * ky[i][5]  + a107 * ky[i][6]  + a108 * ky[i][7]);
			zt[i]  = z[i]  + dt * (a101 * kz[i][0]  + a104 * kz[i][3]  + a105 * kz[i][4]  + a106 * kz[i][5]  + a107 * kz[i][6]  + a108 * kz[i][7]);
			vxt[i] = vx[i] + dt * (a101 * kvx[i][0] + a104 * kvx[i][3] + a105 * kvx[i][4] + a106 * kvx[i][5] + a107 * kvx[i][6] + a108 * kvx[i][7]);
			vyt[i] = vy[i] + dt * (a101 * kvy[i][0] + a104 * kvy[i][3] + a105 * kvy[i][4] + a106 * kvy[i][5] + a107 * kvy[i][6] + a108 * kvy[i][7]);
			vzt[i] = vz[i] + dt * (a101 * kvz[i][0] + a104 * kvz[i][3] + a105 * kvz[i][4] + a106 * kvz[i][5] + a107 * kvz[i][6] + a108 * kvz[i][7]);
		}

		for(int i = 0; i < N; ++i){

			kx[i][9] = vxt[i];
			ky[i][9] = vyt[i];
			kz[i][9] = vzt[i];

			ax = 0.0;
			ay = 0.0;
			az = 0.0;

			for(int j = 0; j < N; ++j){
				if(i != j){
					acc(m, xt, yt, zt, ax, ay, az, i, j);
				}
			}

			kvx[i][9] = ax;
			kvy[i][9] = ay;
			kvz[i][9] = az;

		}
		//stage 11
		for(int i = 0; i < N; ++i){
			xt[i]  = x[i]  + dt * (a111 * kx[i][0]  + a114 * kx[i][3]  + a115 * kx[i][4]  + a116 * kx[i][5]  + a117 * kx[i][6]  + a118 * kx[i][7]  + a119 * kx[i][8]);
			yt[i]  = y[i]  + dt * (a111 * ky[i][0]  + a114 * ky[i][3]  + a115 * ky[i][4]  + a116 * ky[i][5]  + a117 * ky[i][6]  + a118 * ky[i][7]  + a119 * ky[i][8]);
			zt[i]  = z[i]  + dt * (a111 * kz[i][0]  + a114 * kz[i][3]  + a115 * kz[i][4]  + a116 * kz[i][5]  + a117 * kz[i][6]  + a118 * kz[i][7]  + a119 * kz[i][8]);
			vxt[i] = vx[i] + dt * (a111 * kvx[i][0] + a114 * kvx[i][3] + a115 * kvx[i][4] + a116 * kvx[i][5] + a117 * kvx[i][6] + a118 * kvx[i][7] + a119 * kvx[i][8]);
			vyt[i] = vy[i] + dt * (a111 * kvy[i][0] + a114 * kvy[i][3] + a115 * kvy[i][4] + a116 * kvy[i][5] + a117 * kvy[i][6] + a118 * kvy[i][7] + a119 * kvy[i][8]);
			vzt[i] = vz[i] + dt * (a111 * kvz[i][0] + a114 * kvz[i][3] + a115 * kvz[i][4] + a116 * kvz[i][5] + a117 * kvz[i][6] + a118 * kvz[i][7] + a119 * kvz[i][8]);
		}

		for(int i = 0; i < N; ++i){

			kx[i][10] = vxt[i];
			ky[i][10] = vyt[i];
			kz[i][10] = vzt[i];

			ax = 0.0;
			ay = 0.0;
			az = 0.0;

			for(int j = 0; j < N; ++j){
				if(i != j){
					acc(m, xt, yt, zt, ax, ay, az, i, j);
				}
			}

			kvx[i][10] = ax;
			kvy[i][10] = ay;
			kvz[i][10] = az;

		}
		//stage 12
		for(int i = 0; i < N; ++i){
			xt[i]  = x[i]  + dt * (a121 * kx[i][0]  + a126 * kx[i][5]  + a127 * kx[i][6]  + a128 * kx[i][7]  + a129 * kx[i][8]  + a1210 * kx[i][9]);
			yt[i]  = y[i]  + dt * (a121 * ky[i][0]  + a126 * ky[i][5]  + a127 * ky[i][6]  + a128 * ky[i][7]  + a129 * ky[i][8]  + a1210 * ky[i][9]);
			zt[i]  = z[i]  + dt * (a121 * kz[i][0]  + a126 * kz[i][5]  + a127 * kz[i][6]  + a128 * kz[i][7]  + a129 * kz[i][8]  + a1210 * kz[i][9]);
			vxt[i] = vx[i] + dt * (a121 * kvx[i][0] + a126 * kvx[i][5] + a127 * kvx[i][6] + a128 * kvx[i][7] + a129 * kvx[i][8] + a1210 * kvx[i][9]);
			vyt[i] = vy[i] + dt * (a121 * kvy[i][0] + a126 * kvy[i][5] + a127 * kvy[i][6] + a128 * kvy[i][7] + a129 * kvy[i][8] + a1210 * kvy[i][9]);
			vzt[i] = vz[i] + dt * (a121 * kvz[i][0] + a126 * kvz[i][5] + a127 * kvz[i][6] + a128 * kvz[i][7] + a129 * kvz[i][8] + a1210 * kvz[i][9]);
		}

		for(int i = 0; i < N; ++i){

			kx[i][11] = vxt[i];
			ky[i][11] = vyt[i];
			kz[i][11] = vzt[i];

			ax = 0.0;
			ay = 0.0;
			az = 0.0;

			for(int j = 0; j < N; ++j){
				if(i != j){
					acc(m, xt, yt, zt, ax, ay, az, i, j);
				}
			}

			kvx[i][11] = ax;
			kvy[i][11] = ay;
			kvz[i][11] = az;

		}
		//stage 13
		for(int i = 0; i < N; ++i){
			xt[i]  = x[i]  + dt * (a131 * kx[i][0]  + a134 * kx[i][3]  + a135 * kx[i][4]  + a136 * kx[i][5]  + a137 * kx[i][6]  + a138 * kx[i][7]  + a139 * kx[i][8]  + a1310 * kx[i][9] + kx[i][11]);
			yt[i]  = y[i]  + dt * (a131 * ky[i][0]  + a134 * ky[i][3]  + a135 * ky[i][4]  + a136 * ky[i][5]  + a137 * ky[i][6]  + a138 * ky[i][7]  + a139 * ky[i][8]  + a1310 * ky[i][9] + ky[i][11]);
			zt[i]  = z[i]  + dt * (a131 * kz[i][0]  + a134 * kz[i][3]  + a135 * kz[i][4]  + a136 * kz[i][5]  + a137 * kz[i][6]  + a138 * kz[i][7]  + a139 * kz[i][8]  + a1310 * kz[i][9] + kz[i][11]);
			vxt[i] = vx[i] + dt * (a131 * kvx[i][0] + a134 * kvx[i][3] + a135 * kvx[i][4] + a136 * kvx[i][5] + a137 * kvx[i][6] + a138 * kvx[i][7] + a139 * kvx[i][8] + a1310 * kvx[i][9] + kvx[i][11]);
			vyt[i] = vy[i] + dt * (a131 * kvy[i][0] + a134 * kvy[i][3] + a135 * kvy[i][4] + a136 * kvy[i][5] + a137 * kvy[i][6] + a138 * kvy[i][7] + a139 * kvy[i][8] + a1310 * kvy[i][9] + kvy[i][11]);
			vzt[i] = vz[i] + dt * (a131 * kvz[i][0] + a134 * kvz[i][3] + a135 * kvz[i][4] + a136 * kvz[i][5] + a137 * kvz[i][6] + a138 * kvz[i][7] + a139 * kvz[i][8] + a1310 * kvz[i][9] + kvz[i][11]);
		}

		for(int i = 0; i < N; ++i){

			kx[i][12] = vxt[i];
			ky[i][12] = vyt[i];
			kz[i][12] = vzt[i];

			ax = 0.0;
			ay = 0.0;
			az = 0.0;

			for(int j = 0; j < N; ++j){
				if(i != j){
					acc(m, xt, yt, zt, ax, ay, az, i, j);
				}
			}

			kvx[i][12] = ax;
			kvy[i][12] = ay;
			kvz[i][12] = az;

		}
	
		//error estimation
		for(int i = 0; i < N; ++i){
			errorkx[i] = abs( er * (kx[i][0] + kx[i][10] - kx[i][11] - kx[i][12]));
			errorky[i] = abs( er * (ky[i][0] + ky[i][10] - ky[i][11] - ky[i][12]));
			errorkz[i] = abs( er * (kz[i][0] + kz[i][10] - kz[i][11] - kz[i][12]));
			errorkvx[i] = abs( er * (kvx[i][0] + kvx[i][10] - kvx[i][11] - kvx[i][12]));
			errorkvy[i] = abs( er * (kvy[i][0] + kvy[i][10] - kvy[i][11] - kvy[i][12]));
			errorkvz[i] = abs( er * (kvz[i][0] + kvz[i][10] - kvz[i][11] - kvz[i][12]));
		}

		double error = 0.0;
		for(int i = 0; i < N; ++i){
			error = fmax(error, errorkx[i]);
			error = fmax(error, errorky[i]);
			error = fmax(error, errorkz[i]);
			error = fmax(error, errorkvx[i]);
			error = fmax(error, errorkvy[i]);
			error = fmax(error, errorkvz[i]);
		}


		double e = 1.0e-10;
		double ee = 1.0/7.0;
		double s = 1.0;

		s = 0.64 * pow( e  / error, ee);
	
		for(int i = 0; i < N; ++i){
		
			//RKF78
			x[i] += dt * (b1 * kx[i][0] + b6 * kx[i][5] + b7 * kx[i][6] + b8 * kx[i][7] + b9 * kx[i][8] + b10 * kx[i][9] + b11 * kx[i][10]);
			y[i] += dt * (b1 * ky[i][0] + b6 * ky[i][5] + b7 * ky[i][6] + b8 * ky[i][7] + b9 * ky[i][8] + b10 * ky[i][9] + b11 * ky[i][10]);
			z[i] += dt * (b1 * kz[i][0] + b6 * kz[i][5] + b7 * kz[i][6] + b8 * kz[i][7] + b9 * kz[i][8] + b10 * kz[i][9] + b11 * kz[i][10]);

			vx[i] += dt * (b1 * kvx[i][0] + b6 * kvx[i][5] + b7 * kvx[i][6] + b8 * kvx[i][7] + b9 * kvx[i][8] + b10 * kvx[i][9] + b11 * kvx[i][10]);
			vy[i] += dt * (b1 * kvy[i][0] + b6 * kvy[i][5] + b7 * kvy[i][6] + b8 * kvy[i][7] + b9 * kvy[i][8] + b10 * kvy[i][9] + b11 * kvy[i][10]);
			vz[i] += dt * (b1 * kvz[i][0] + b6 * kvz[i][5] + b7 * kvz[i][6] + b8 * kvz[i][7] + b9 * kvz[i][8] + b10 * kvz[i][9] + b11 * kvz[i][10]);

		}




		if(t % outInterval == 0){
			
			if(outHelio == 1){
				sprintf(outfilename, "Outhelio10_%.12lld.dat", t);
			}
			else{
				sprintf(outfilename, "Outbary10_%.12lld.dat", t);
			}
			outfile = fopen(outfilename, "w");
			printf("%s\n", outfilename);
			
			//convert to heliocentric output
			comx = 0.0;
			comy = 0.0;
			comz = 0.0;
			vcomx = 0.0;
			vcomy = 0.0;
			vcomz = 0.0;
			
			if(outHelio == 1){
				for(int i = 1; i < N; ++i){
					comx += m[i] * x[i];
					comy += m[i] * y[i];
					comz += m[i] * z[i];
					vcomx += m[i] * vx[i];
					vcomy += m[i] * vy[i];
					vcomz += m[i] * vz[i];
				}
				comx /= m[0];
				comy /= m[0];
				comz /= m[0];
				vcomx /= m[0];
				vcomy /= m[0];
				vcomz /= m[0];
			}
			for(int i = 0; i < N; ++i){
				fprintf(outfile, "%.10g %d %.40g %.40g %.40g %.40g %.40g %.40g %.40g %.40g 0 0 0 0 0 0 0 0 0 0 0\n", time + t * dt / dayUnit, i, m[i], rad[i], comx + x[i], comy + y[i], comz + z[i], vcomx + vx[i], vcomy + vy[i], vcomz + vz[i]);
				
			}
			fclose(outfile);
		}
	//dt *= s;

	}
	
}
	

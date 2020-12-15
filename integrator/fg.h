void fgfull(double &xi, double &yi, double &zi, double &vxi, double &vyi, double &vzi, double dt, double mu, int GR){

	double f,g,fd,gd;                               /* Gauss's f, g, fdot and gdot */
	double rsq,vsq,ir;
	double u;                                       /* r v cos(phi) */
	double ia;                                      //inverse of a
	double ria;
	double air_;
	double e;                                        /* eccentricity */
	double ec,es;                             	  /* e cos(E), e sin(E) */
	double ien;                                   /* inverse mean motion */
	double en;                                    /* mean motion */
	double dec;                                      /* delta E */
	double dm;
	double mw;                                       /* minus function to zero */
	double wp;                                       /* first derivative */
	double iwp;
	double wpp;                                      /* second derivative */
	double wppp;                                     /* third derivative */
	double dx,s,c;
	double t1;
	double tx, ty, tz;
	const double DOUBLE_EPS = 1.2e-16;
	double converge;
	double UP = 2*M_PI;
	double LOW = -2*M_PI;
	double next;
	int i;
	/*
	* Evaluate some orbital quantites.
	*/

	rsq = xi*xi + yi*yi + zi*zi;
	vsq = vxi*vxi + vyi*vyi + vzi*vzi;
	u =  xi*vxi + yi*vyi + zi*vzi;
	ir = 1.0 / sqrt(rsq);
	ia = 2.0*ir-vsq/mu;

		
	if(GR == 1){// GR time rescale (Saha & Tremaine 1994)
		double c = 10065.3201686;//c in AU / day * 0.0172020989
		double c2 = c * c;
		dt *= 1.0 - 1.5 * mu * ia / c2;
	}
	

//printf("%g %g %g %g %g %g\n", rsq, vsq, u, ir, ia, mu);
	if(ia > 0){

		t1 = ia*ia;
		ria = rsq*ir*ia;
		ien = 1.0 / sqrt(mu*t1*ia);
		en = 1.0/ien;
		ec = 1.0-ria;
		es = u*t1*ien;
		e = sqrt(ec*ec + es*es);
		dm = en * dt - es;
		if ((es*cos(dm)+ec*sin(dm)) > 0){
			dec = dm + 0.85*e;
		}
		else dec = dm - 0.85*e;
		converge = fabs(en * dt *DOUBLE_EPS);

		for(i = 0; i < 128; ++i) {

			s = sin(dec);
			c = cos(dec);
			//sincos(dec, &s, &c);
			wpp = ec*s + es*c;
			wppp = ec*c - es*s;
			mw = dm - dec + wpp;
//printf("%.10g %.10g %.10g\n", mw, dec, wpp);
			if(mw < 0.0){
				UP = dec;
			}
			else LOW = dec;
			wp = 1.0 - wppp;
			wpp *= 0.5;
			dx = mw/wp;
			dx = mw/(wp + dx*wpp);
			dx = mw/(wp + dx*(wpp + (1.0/6.0)*dx*wppp));
			next = dec + dx;
//printf("dx %d %.20g %g %g %g %.20g\n", i, dx, wp, mw, dm, dec);
			if (fabs(dx) <= converge) break;
			if(next > LOW && next < UP){
				dec = next;
			}
			else dec = 0.5*(LOW + UP);
			if (dec==LOW || dec==UP) break;
		}
		if(i < 127){
			iwp = 1.0/wp;
			air_ = -1.0/ria;
			t1 = (1.0-c);
			f = 1.0 + air_*t1;
			g = dt + (s-dec)*ien;
			fd = air_*iwp*s*en;
			gd = 1.0 - iwp*t1;

			tx = f*xi+g*vxi;
			ty = f*yi+g*vyi;
			tz = f*zi+g*vzi;
			vxi = fd*xi+gd*vxi;
			vyi = fd*yi+gd*vyi;
			vzi = fd*zi+gd*vzi;
			xi = tx;
			yi = ty;
			zi = tz;
//printf("%g %g\n", xi, vyi);
		}
		else{
printf("not converged\n");
		}

	}
	else{
printf("Negative a %g %g\n", ia, 1.0/ia);
	}
}

//store dec
void fgfull2(double &xi, double &yi, double &zi, double &vxi, double &vyi, double &vzi, double &dec, double dt, double mu, int GR, int t){

	double f,g,fd,gd;                               /* Gauss's f, g, fdot and gdot */
	double rsq,vsq,ir;
	double u;                                       /* r v cos(phi) */
	double ia;                                      //inverse of a
	double ria;
	double air_;
	double e;                                        /* eccentricity */
	double ec,es;                             	  /* e cos(E), e sin(E) */
	double ien;                                   /* inverse mean motion */
	double en;                                    /* mean motion */
	//double dec;                                      /* delta E */
	double dm;
	double mw;                                       /* minus function to zero */
	double wp;                                       /* first derivative */
	double iwp;
	double wpp;                                      /* second derivative */
	double wppp;                                     /* third derivative */
	double dx,s,c;
	double t1;
	double tx, ty, tz;
	const double DOUBLE_EPS = 1.2e-16;
	double converge;
	double UP = 2*M_PI;
	double LOW = -2*M_PI;
	double next;
	int i;
	/*
	* Evaluate some orbital quantites.
	*/

	rsq = xi*xi + yi*yi + zi*zi;
	vsq = vxi*vxi + vyi*vyi + vzi*vzi;
	u =  xi*vxi + yi*vyi + zi*vzi;
	ir = 1.0 / sqrt(rsq);
	ia = 2.0*ir-vsq/mu;
	
	if(GR == 1){// GR time rescale (Saha & Tremaine 1994)
		double c = 10065.3201686;//c in AU / day * 0.0172020989
		double c2 = c * c;
		dt *= 1.0 - 1.5 * mu * ia / c2;
	}
	
	if(ia > 0){

		t1 = ia*ia;
		ria = rsq*ir*ia;
		ien = 1.0 / sqrt(mu*t1*ia);
		en = 1.0/ien;
		ec = 1.0-ria;
		es = u*t1*ien;
		e = sqrt(ec*ec + es*es);
		dm = en * dt - es;
		if(t == 1){
			if ((es*cos(dm)+ec*sin(dm)) > 0){
				dec = dm + 0.85*e;
			}
			else dec = dm - 0.85*e;
		}
		converge = fabs(en * dt *DOUBLE_EPS);
		for(i = 0; i < 128; ++i) {

			s = sin(dec);
			c = cos(dec);
			//sincos(dec, &s, &c);
			wpp = ec*s + es*c;
			wppp = ec*c - es*s;
			mw = dm - dec + wpp;
//printf("%.10g %.10g %.10g\n", mw, dec, wpp);
			if(mw < 0.0){
				UP = dec;
			}
			else LOW = dec;
			wp = 1.0 - wppp;
			wpp *= 0.5;
			dx = mw/wp;
			dx = mw/(wp + dx*wpp);
			dx = mw/(wp + dx*(wpp + (1.0/6.0)*dx*wppp));
			next = dec + dx;
//printf("dx %d %.20g %g %g %g %.20g\n", i, dx, wp, mw, dm, dec);
			if (fabs(dx) <= converge) break;
			if(next > LOW && next < UP){
				dec = next;
			}
			else dec = 0.5*(LOW + UP);
			if (dec==LOW || dec==UP) break;
		}
		if(i < 127){
			iwp = 1.0/wp;
			air_ = -1.0/ria;
			t1 = (1.0-c);
			f = 1.0 + air_*t1;
			g = dt + (s-dec)*ien;
			fd = air_*iwp*s*en;
			gd = 1.0 - iwp*t1;

			tx = f*xi+g*vxi;
			ty = f*yi+g*vyi;
			tz = f*zi+g*vzi;
			vxi = fd*xi+gd*vxi;
			vyi = fd*yi+gd*vyi;
			vzi = fd*zi+gd*vzi;
			xi = tx;
			yi = ty;
			zi = tz;
//printf("%g %g\n", xi, vyi);
		}
		else{
printf("not converged\n");
		}

	}
	else{
printf("Negative a %g %g\n", ia, 1.0/ia);
	}
}


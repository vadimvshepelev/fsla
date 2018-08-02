
#ifndef EOSANALYTIC_H
#define EOSANALYTIC_H

#include "eos.h"

// Аналитическое уравнение состояния для газообразного алюминия

class EOSAnalytic : public EOS {
public:
	EOSAnalytic();

	EOSType getType() { return analytic; }
	
	double    getpi(double ro, double ti);
	double    getpe(double ro, double ti, double te);
	
	double    getei(double ro, double ti);
	double    getee(double ro, double ti, double te);
		
	double    getti(double ro, double ei);
	double    gette(double ro, double ti, double ee);

	double     getC(double ro, double ti, double te);

	double    getci(double ro, double ti);
	double    getce(double ro, double te);
	double getkappa(double ro, double ti, double te);
	double getAlpha(double ro, double ti, double te);

	double getphase(double ro, double ti);
	double   getmix(double ro, double ti);
	double  getnuWR(double ro, double ti, double te, double b);
	double getEntropy(double ro, double ti, double te);

	double getGamma(void);
/*
	double getdpdro  (double ro, double ti, double te);
	double getdpdroe (double ro, double ti, double te);
	double getdpdroei(double ro, double ti, double te);

	double getdpdt   (double ro, double ti, double te);
	double getdedt   (double ro, double ti, double te);

	double getdpdro_rov_roE(double ro, double ti, double te, double v);
	double getdpdrov_ro_roE(double ro, double ti, double te, double v);
	double getdpdroE_ro_rov(double ro, double ti, double te, double v);
	*/
protected:
	double __ee(double ro, double teta);
	double __pe(double ro, double teta);
	double solve_ti(double ro, double ei, double low_border, double high_border);
	double solve_te(double ro, double ti, double ee, double low_border, double high_border);
	double dpde(double ro, double ti, double te);
	double dpdei(double ro, double ti, double te);
	double dpdti(double ro, double ti, double te);
	double dedti(double ro, double ti, double te);
	double deidti(double ro, double ti, double te);
	double dpdte(double ro, double ti, double te);
	double dedte(double ro, double ti, double te);
};


class EOSPyrexGlass: public EOS {
	double  getnuWR(double ro, double ti, double te, double b) {return 0.;}
public:
	EOSPyrexGlass() {}
	EOSType getType() { return EOSType::analyticPyrexGlass; }
	double    getpi(double ro, double ti);
	double    getpe(double ro, double ti, double te) {return 0.;}
	double    getei(double ro, double ti);
	double    getee(double ro, double ti, double te) {return 0.;}
	double    getti(double ro, double ei);
	double    gette(double ro, double ti, double ee) {return 0.;}
	double     getC(double ro, double ti, double te) {return 0.;}
	double    getci(double ro, double ti) {return 1.; };
	double    getce(double ro, double te) {return 1.;}
	double getkappa(double ro, double ti, double te) {return 0.;}
	double getAlpha(double ro, double ti, double te) {return 0.;}
	double getphase(double ro, double ti) {return 0.;}
	double   getmix(double ro, double ti) {return 0.;}
	double getEntropy(double ro, double ti, double te) {return 0.;}
	double getGamma(void) {return 0.;};
/*	double getdpdro  (double ro, double ti, double te) {return 0.;}
	double getdpdroe (double ro, double ti, double te) {return 0.;}
	double getdpdroei(double ro, double ti, double te) {return 0.;}
	double getdpdt   (double ro, double ti, double te) {return 0.;}
	double getdedt   (double ro, double ti, double te) {return 0.;}
*/
protected:
	double __ee(double ro, double teta) {return 0.;}
	double __pe(double ro, double teta) {return 0.;}
	double solve_ti(double ro, double ei, double low_border, double high_border);
/*	double dpde(double ro, double ti, double te) {return 0.;}
	double dpdei(double ro, double ti, double te) {return 0.;}
	double dpdti(double ro, double ti, double te) {return 0.;}
	double dedti(double ro, double ti, double te) {return 0.;}
	double deidti(double ro, double ti, double te) {return 0.;}
	double dpdte(double ro, double ti, double te) {return 0.;}
	double dedte(double ro, double ti, double te) {return 0.;}*/
};


class EOSSimpleWater: public EOS {
	double  getnuWR(double ro, double ti, double te, double b) {return 0.;}
public:
	EOSSimpleWater() {}
	EOSType getType() { return EOSType::simpleWater; }
	double    getpi(double ro, double ti);
	double    getpe(double ro, double ti, double te) {return 0.;}
	double    getei(double ro, double ti);
	double    getee(double ro, double ti, double te) {return 0.;}
	double    getti(double ro, double ei);
	double    gette(double ro, double ti, double ee) {return 0.;}
	double     getC(double ro, double ti, double te) {return 0.;}
	double    getci(double ro, double ti);
	double    getce(double ro, double te) {return 1.;}
	double getkappa(double ro, double ti, double te) {return 0.;}
	double getAlpha(double ro, double ti, double te) {return 0.;}
	double getphase(double ro, double ti) {return 0.;}
	double   getmix(double ro, double ti) {return 0.;}
	double getEntropy(double ro, double ti, double te) {return 0.;}
	double getGamma(void) {return 0.;}
/*	double getdpdro  (double ro, double ti, double te) {return 0.;}
	double getdpdroe (double ro, double ti, double te) {return 0.;}
	double getdpdroei(double ro, double ti, double te) {return 0.;}
	double getdpdt   (double ro, double ti, double te) {return 0.;}
	double getdedt   (double ro, double ti, double te) {return 0.;}
*/
protected:
	double __ee(double ro, double teta) {return 0.;}
	double __pe(double ro, double teta) {return 0.;}
	double solve_ti(double ro, double ei, double low_border, double high_border);
/*	double dpde(double ro, double ti, double te) {return 0.;}
	double dpdei(double ro, double ti, double te) {return 0.;}
	double dpdti(double ro, double ti, double te) {return 0.;}
	double dedti(double ro, double ti, double te) {return 0.;}
	double deidti(double ro, double ti, double te) {return 0.;}
	double dpdte(double ro, double ti, double te) {return 0.;}
	double dedte(double ro, double ti, double te) {return 0.;}*/
};


class EOSSimpleCr : public EOSAnalytic {
public:
	EOSType getType() { return simpleCr; }
	double    getei(double ro, double ti) { return (1.95*ti/1000. + 3.)*1.e9/7190.; }	// [J/kg]
	double    getee(double ro, double ti, double te) { return (__ee(ro, te) - __ee(ro, ti)); } // [J/kg]
	double    gette(double ro, double ti, double ee) { return solve_te(ro, ti, ee, 0.00001, 100000.0); } // [K]
	double     getC(double ro, double ti, double te) { return 4.e3;}	 // [m/s]
	double    getci(double ro, double ti) {
		double M = 52.e-3, NA = 6.0e23, kB = 1.38e-23, m0 = M/NA, n = 7190./m0;
		double ci = 3.0*n*kB;
		return ci; 
	} // [J/K/m^3]   
	double    getce(double ro, double te) { return 1./(0.034/(te/1000.) + .14/sqrt(te/1000.)); }	// [J/K/m^3]
	double getkappa(double ro, double ti, double te) {
		double EFkK = 100.;
		double t = 6.*(te/1000.)/EFkK, tm = 6*1.3/EFkK;
		double Sse = 3.8e-4 * t / (1.-2.*sqrt(t) + 1.4*t + 0.04*t*t);
		double TckK = 7., beta = 2.4/(1.+.35*(ti/1000.)/TckK);
		double CC = t*(1+11.2*t*t)/(1.+3.35*pow(t,2.05));
		double CC13 =  tm*(1+11.2*tm*tm)/(1.+3.35*pow(tm,2.05));
		double w = .0043 * (1. + 126.*(ti/1000.)/TckK) / (1.+0.86*(ti/1000.)/TckK);
		double TmkK = 1.3;
		double wTmkK = .0043 * (1. + 126.*TmkK/TckK) / (1.+0.86*TmkK/TckK);
		double kapSIliq = 140.*CC/CC13 * wTmkK/w;
		double kappa = 1./(Sse + 1./kapSIliq);
		return 0.6*kappa;
	} // [W/K/m]
	double getAlpha(double ro, double ti, double te) { return -3.e17;} // [W/K/m^3]
protected:
	double __ee(double ro, double te) {return (1.e8*(0.84*sqrt(te/1000.) - 1.73*(te/1000.) + 4.76*pow(te/1000., 3./2.))/7190.); } // [J/kg]
};

class EOSSimpleAu : public EOSAnalytic {
	EOSType getType() { return simpleAu; }
	double    getei(double ro, double ti) { return (1.315*ti/1000. + 2.7)*1.e9/19320.;} // [J/kg]
	double    getee(double ro, double ti, double te) { return (__ee(ro, te) - __ee(ro, ti)); } // [J/kg]
	double    gette(double ro, double ti, double ee) { return solve_te(ro, ti, ee, 0.00001, 100000.0); } // [K]
	double     getC(double ro, double ti, double te) { return 4.e3;}	 // [m/s]
	double    getci(double ro, double ti) {
		double M = 197.e-3, NA = 6.0e23, kB = 1.38e-23, m0 = M/NA, n = 19320./m0;
		double ci = 3.0*n*kB;
		return ci; 
	} // [J/K/m^3]                                           
	double    getce(double ro, double te) { return 1./(0.034/(te/1000.) + .14/sqrt(te/1000.)); }  // [J/K/m^3]
	double getkappa(double ro, double ti, double te) {
		double EFkK = 100.;
		double t = 6.*(te/1000.)/EFkK, tm = 6*1.3/EFkK;
		double Sse = 3.8e-4 * t / (1.-2.*sqrt(t) + 1.4*t + 0.04*t*t);
		double TckK = 7., beta = 2.4/(1.+.35*(ti/1000.)/TckK);
		double CC = t*(1+11.2*t*t)/(1.+3.35*pow(t,2.05));
		double CC13 =  tm*(1+11.2*tm*tm)/(1.+3.35*pow(tm,2.05));
		double w = .0043 * (1. + 126.*(ti/1000.)/TckK) / (1.+0.86*(ti/1000.)/TckK);
		double TmkK = 1.3;
		double wTmkK = .0043 * (1. + 126.*TmkK/TckK) / (1.+0.86*TmkK/TckK);
		double kapSIliq = 140.*CC/CC13 * wTmkK/w;
		double kappa = 1./(Sse + 1./kapSIliq);
		return kappa;
	} // [W/K/m]
	double getAlpha(double ro, double ti, double te) {
		double Ka = 4.;
		double eV = 11600.; // [K]
		double tEv = te/eV;
		double alpha = -1.e17 * (.2+4.3/Ka*pow(tEv,3.6) / (1. + pow(tEv,3.5) + .9*pow(tEv,4.1)));
		return alpha;
	}// [W/K/m^3]
protected:
	double __ee(double ro, double te) {return (1.e8*(0.84*sqrt(te/1000.) - 1.73*(te/1000.) + 4.76*pow(te/1000., 3./2.))/19320.); } // [J/kg]

};

class EOSSimpleSi : public EOSAnalytic {
	EOSType getType() { return simpleSi; }
	double    getei(double ro, double ti) {return (1.17*ti/1000. + 2.)*1.e9/2330.;} // [J/kg]
	double    getee(double ro, double ti, double te) { return (__ee(ro, te) - __ee(ro, ti)); } // [J/kg]
	double    gette(double ro, double ti, double ee) { return solve_te(ro, ti, ee, 0.00001, 100000.0); } // [K]
	double     getC(double ro, double ti, double te) { return 4.e3;}	 // [m/s]
	double    getci(double ro, double ti) {
		double M = 28.e-3, NA = 6.0e23, kB = 1.38e-23, m0 = M/NA, n = 2330./m0;
		double ci = 3.0*n*kB;
		return ci; 
	} // [J/K/m^3]   
	double    getce(double ro, double te) { return 1./(0.034/(te/1000.) + .14/sqrt(te/1000.)); }	// [J/K/m^3]
	double getkappa(double ro, double ti, double te) { return 0.; /*150.;*/ }	 // [W/m/K]
	double getAlpha(double ro, double ti, double te) { return -3.e17;; } // [W/K/m^3]
protected:
	double __ee(double ro, double te) {return (1.e8*(0.84*sqrt(te/1000.) - 1.73*(te/1000.) + 4.76*pow(te/1000., 3./2.))/2330.); } // [J/kg]
};


class EOSMieGruneisenRu : public EOS {
	const double A, a, b, ro0, v0;
	double al1, al2, bet1, bet2, c1, c2;
	double G(double x) {
		//return .5*(1. + c2)*(pow(x, al2)/(1. + c2*pow(x, bet1)))/((1. + c1)*pow(x, al1))*(a + b + 5./3); 
		//1/3 + ( alpha + ( alpha - betaa ) cc x^betaa )/2/( 1 + cc x^betaa )
		return 1./3. + ( al1 + (al1 - bet1) * c1 * pow(x, bet1)/2./( 1 + c1*pow(x, bet1)));
	}
public:
//	EOSMieGruneisenRu() : A(1.9642003675331705), a(1.5426884720681746), b(0.9792335085217448), ro0(12410.), v0(13.524/0.148) {
//	EOSMieGruneisenRu() : A(2.6063032461372684), a(1.4467314356992462), b(1.007457482946519), ro0(12470.), v0(13.4587/0.148) {
//	EOSMieGruneisenRu() : A(2.024682220294127), a(1.6569036734344376), b(1.1736399604338026), ro0(12470.), v0(90.6514925212472) {
	EOSMieGruneisenRu() : A(2.024682220294127), 
		                  a(1.6569036734344376), 
						  b(1.1736399604338026), 
						  ro0(12506.32423124619377), 
						  v0(13.416420893148256 / 0.148) {
		al1  = 2.*a + 1.;
		bet1 = a + 1.;
		c1   = (a - b)/(b + 1.);
		bet2 = (a + 1.)*(a + 2./3.)/(a + b + 5./3.);
		al2  = bet2 + a; 
		c2   = (a + 1.)*(a + 2./3.)/(b + 1.)/(b + 2./3.) - 1.;	
	}
	EOSMieGruneisenRu(const double _A, 
			          const double _a, const double _b, 
					  const double _ro0, const double _v0,
					  double _al1, double _al2, double _bet1, double _bet2, double _c1, double _c2) : 
	                  A(_A), a(_a), b(_b), ro0(_ro0), v0(_v0), al1(_al1), al2(_al2), bet1(_bet1), bet2(_bet2), c1(_c1), c2(_c2) {}
	EOSType getType() { return MieGruneisen; }	
	double    getpi(double ro, double ti);					// [Pa]
	double    getpe(double ro, double ti, double te);		// [Pa]	
	double    getei(double ro, double ti);					// [J/kg]
	double    getee(double ro, double ti, double te);		// [J/kg]	
	double    getti(double ro, double ei);					// [K]
	double    gette(double ro, double ti, double ee);		// [K]
	double     getC(double ro, double ti, double te);
	double    getci(double ro, double ti);
	double    getce(double ro, double te);
	double getkappa(double ro, double ti, double te);
	double getAlpha(double ro, double ti, double te);
	double getphase(double ro, double ti);
	double   getmix(double ro, double ti);
	double  getnuWR(double ro, double ti, double te, double b);
	double getEntropy(double ro, double ti, double te);
	double getGamma(void);
protected:
	double solve_ti(double ro, double ei, double low_border, double high_border);
	double solve_te(double ro, double ti, double ee, double low_border, double high_border);
};


#endif
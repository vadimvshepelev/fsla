#include "eosAnalyticAu.h"
#include "defines.h"
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

// Заглушки

double EOSAnalyticAu::getpi(double ro, double ti)    { return 1.0; }
double EOSAnalyticAu::getei(double ro, double ti)    { return 1.0;	}
double EOSAnalyticAu::getci(double ro, double ti)    {	return 1.0; }
double EOSAnalyticAu::getC(double ro, double ti, double te) {	return 1.0; }
double EOSAnalyticAu::getphase(double ro, double ti) { return 1.0; }
double EOSAnalyticAu::getmix(double ro, double ti)   { return 1.0; }
//double  EOSAnalyticAu::getnuWR(double ro, double ti, double te, double b) { return 0.0; }
double EOSAnalyticAu::getkappa(double ro, double ti, double te) { return 1.0; }
double EOSAnalyticAu::getEntropy(double ro, double ti, double te) { return 0.0; }
double EOSAnalyticAu::getGamma(void) {return 0.0;}

// Основные функции

double EOSAnalyticAu::getee(double ro, double ti, double te) {
	double PI = 3.14159;
	double tn = atan((te-3500.)/700.)*2./PI;
    double ee = 1.e-10/19.3*338*te*te* ((1. - tn)/2. + 1.4*(1. + tn)/2.);  // [kJ/g]
	// [kJ/g] = SI [10^3 J / 10^-3 kg] = [10^6 J/kg]
	return ee * 1.e6; //[J/kg]
}

double EOSAnalyticAu::getpe(double ro, double ti, double te) {
	double pe = 2./3.*ro*getee(ro, ti, te);
	return pe; // [Pa]
}

double EOSAnalyticAu::getce(double ro, double te) {
	double PI = 3.14159265359;
	double tn = atan((te-3500.)/700.)*2./PI;
//    double ce = 1.e-10/19.3*338.*
//		      ( 0.4/700./PI*te*te / ( 1. + (te-3500.)/700.*(te-3500.)/700. )  + 
//			    2.*te* ( (1. - tn)/2. + 1.4*(1. + tn)/2. ) );  // [kJ/g/K]
	double A = 0.4/700./PI*te*te / ( 1. + (te-3500.)*(te-3500.)/700./700. );
	double B = 2.*te* ( (1.-tn)/2 + 1.4*(1.+tn)/2);

	double ce = 1.e-10/19.3*338. * (A + B); // [kJ/g/K]
	// [kJ/g/K] = [10^3 J/10^-3 kg/K] = [10^6 J/kg/K]
	// Для [J/m3/K] умножаем на плотность ro
	return ce*1.e6*ro; // [J/m3/K]
}

double EOSAnalyticAu::getAlpha(double ro, double ti, double te) {
	double kBeV = 8.817343e-5; // [eV/K]
	double teeV = kBeV*te;     // [eV]
	double k_alpha = 4.;       
	double alphaei = (0.2 + 4.3*pow(teeV, 3.6)/(1.+pow(teeV, 3.5)+0.9*pow(teeV, 4.1))/k_alpha)*1.e8/(ro0/1000.);
	// [kW/g/K] = SI [10^3 J/s/10^-3 kg/K] = SI [10^6 J/s/kg/K] 
	// Для [W/K/m^3] нужно умножить на плотность ro
	return -ro*alphaei*1.e6;   // [W/K/m^3]
}

// Обратные функции

double EOSAnalyticAu::getti(double ro, double ei) {
	double ti = solve_ti(ro, ei, 0.00001, 60000.0);
	return ti;
}

double EOSAnalyticAu::gette(double ro, double ti, double ee) {
	double te=0.0;
	double eps = 1.0e-10;
	double te_dee, low_dee, high_dee;
	double low_border = 0.001, high_border = 500000.0;
	if(ee==0.0)
		return ti;
	for(;;) {
		te = (low_border + high_border)/2.0;
		if(fabs(low_border-high_border) < eps)
			return te;
		te_dee   = getee(ro, ti, te)-ee;
		low_dee  = getee(ro, ti, low_border)-ee;
		high_dee = getee(ro, ti, high_border)-ee;
		if(fabs(te_dee) < eps)
			return te;
		if(low_dee * high_dee > 0) 	{
			printf("\nsolve_te: Error in initial data -- no equation root on the interval\n");
			printf("ro=%e, ee=%e, F(a)=%e, F(b)=%e\n", ro, ee, low_dee, high_dee);
			exit(1);
		}
		if(te_dee * low_dee > 0)
			low_border  += (high_border-low_border)/2.0;
		else
			high_border -= (high_border-low_border)/2.0;
	}
}

double EOSAnalyticAu::solve_ti(double ro, double ei, double low_border, double high_border) {
	double ti;
	double eps = 0.01;
	double ti_dei, low_dei, high_dei;
	for(;;) {
		ti = (low_border + high_border)/2.0;
		if(fabs(low_border-high_border) < eps)
			return ti;
		ti_dei   = getei(ro, ti)-ei;
		low_dei  = getei(ro, low_border)-ei;
		high_dei = getei(ro, high_border)-ei;
		if(fabs(ti_dei) < eps)
			return ti;
		if(low_dei * high_dei > 0) {
			printf("\nsolve_ti: Error in initial data -- no equation root on the interval\n");
			printf("ro=%e, ei=%e, F(a)=%e, F(b)=%e\n", ro, ei, low_dei, high_dei);
			return -1.;
		}
		if(ti_dei * low_dei > 0)
			low_border  += (high_border-low_border)/2.0;
		else  
			high_border -= (high_border-low_border)/2.0;
	}
}

// Производные
/*
double EOSAnalyticAu::getdpdro(double ro, double ti, double te) {
	return ( getpi(ro + EOS_EPS, ti) + getpe(ro + EOS_EPS, ti, te) 
		   - getpi(ro, ti)			 - getpe(ro, ti, te) ) / EOS_EPS;
}

double EOSAnalyticAu::getdpdroe(double ro, double ti, double te) {
	return ( getpi(ro + EOS_EPS, ti) + getpe(ro + EOS_EPS, ti, te) 
		   - getpi(ro, ti)			 - getpe(ro, ti, te) ) / EOS_EPS;
}

double EOSAnalyticAu::getdpdroei(double ro, double ti, double te) {
	return ( getpi(ro + EOS_EPS, ti) + getpe(ro + EOS_EPS, ti, te) 
		   - getpi(ro, ti)			 - getpe(ro, ti, te) ) / EOS_EPS;
}

double EOSAnalyticAu::getdpdro_rov_roE(double ro, double ti, double te, double v)
{ return 0;}

double EOSAnalyticAu::getdpdrov_ro_roE(double ro, double ti, double te, double v)
{ return 0;}

double EOSAnalyticAu::getdpdroE_ro_rov(double ro, double ti, double te, double v)
{ return 0;}

double EOSAnalyticAu::getdpdt(double ro, double ti, double te) {
	return (getpi(ro, ti+EOS_EPS)-getpi(ro, ti))/EOS_EPS;
}

double EOSAnalyticAu::getdedt(double ro, double ti, double te) {
	return (getei(ro, ti+EOS_EPS)-getei(ro, ti))/EOS_EPS;
}

////////////// PRIVATE ////////////////

double EOSAnalyticAu::dpdei(double ro, double ti, double te) {
	return dpdti(ro, ti, te) / deidti(ro, ti, te) - 
		   dpdte(ro, ti, te) / dedte(ro, ti, te);
}

double EOSAnalyticAu::dpde(double ro, double ti, double te) {
	return dpdte(ro, ti, te) / dedte(ro, ti, te);
}

double EOSAnalyticAu::dpdti(double ro, double ti, double te) {
	return ( getpi(ro, ti + EOS_EPS) + getpe(ro, EOS_EPS, te) 
		   - getpi(ro, ti)           - getpe(ro, ti, te) ) / EOS_EPS;
}

double EOSAnalyticAu::dedti(double ro, double ti, double te) {
	return ( getei(ro, ti + EOS_EPS) + getee(ro, ti + EOS_EPS, te) 
		   - getei(ro, ti)           - getee(ro, ti, te) ) / EOS_EPS;
}

double EOSAnalyticAu::deidti(double ro, double ti, double te) {
	return ( getei(ro, ti + EOS_EPS) - getei(ro, ti) ) / EOS_EPS;
}

double EOSAnalyticAu::dpdte(double ro, double ti, double te) {
	return ( getpe(ro, ti, te + EOS_EPS) - getpe(ro, ti, te) ) / EOS_EPS;

}

double EOSAnalyticAu::dedte(double ro, double ti, double te) {
	return ( getee(ro, ti, te + EOS_EPS) - getee(ro, ti, te) ) / EOS_EPS;
}*/

EOSAnalyticAu::EOSAnalyticAu() {	
/*	MAX_T=100000.0; 
	MIN_T=300.0; 
	MAX_RO=5000.0; 
	MIN_RO=0.001;*/
	ro0 = 19300.;
}

/*complex<double> EOSAnalyticAu::geteps(double ro, double ti, double te, double Z=3.0) {
	return complex<double>(0.0);
}*/

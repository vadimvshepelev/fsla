
#include"eosAnalyticNi.h"

#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>


double EOSAnalyticNi::getpi(double ro, double ti)
{
	return 1.0;
}

double EOSAnalyticNi::getpe(double ro, double ti, double te) {
	double pe = __pe(ro, te) - __pe(ro, ti);
	return pe;
}

double EOSAnalyticNi::__pe(double ro, double teta) {
	double ro0   = 8900.0;
	double K = 1077.0*teta*teta/2.0;
	double S = 1.275e5*pow(teta, 1.3);
	double pe = 1.1*pow(ro/ro0, 1.1)*K*S/(K+S);
	return pe;
}

double EOSAnalyticNi::getei(double ro, double ti) 
{
	return 1.0;	
}

double EOSAnalyticNi::getee(double ro, double ti, double te) {
	double ee = (__ee(ro, te) - __ee(ro, ti)) ; // [J/m^3]
	return ee/ro; // [J/kg]
}

double EOSAnalyticNi::__ee(double ro, double teta)
{
	double ro0 = 8900.0;
	
	// Previous formula
	// double ee  = 1.0e5*pow(ro/ro0, 0.33)*teta*teta/(185.2+0.5882*pow(teta, 0.7));

	double K = 1077.0*teta*teta/2.0;
	double S = 1.275e5*pow(teta, 1.3);
	double ee = pow(ro/ro0, 0.6)*K*S/(K+S);
	return ee; // [J/m^3]
}

double EOSAnalyticNi::getti(double ro, double ei)
{
	double ti = solve_ti(ro, ei, 0.00001, 60000.0);
	return ti;
}

double EOSAnalyticNi::gette(double ro, double ti, double ee) {
	double te=0.0;
	double eps = 1.0e-10;
	double te_dee, low_dee, high_dee;
	double low_border = 0.001, high_border = 500000.0;
	if(ee==0.0)
		return ti;
	for(;;)
	{
		te = (low_border + high_border)/2.0;

		if(fabs(low_border-high_border) < eps)
			return te;

		te_dee   = getee(ro, ti, te)-ee;
		low_dee  = getee(ro, ti, low_border)-ee;
		high_dee = getee(ro, ti, high_border)-ee;

		if(fabs(te_dee) < eps)
			return te;

		if(low_dee * high_dee > 0)
		{
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

double EOSAnalyticNi::solve_ti(double ro, double ei, double low_border, double high_border)
{
	double ti;
	double eps = 0.01;
	double ti_dei, low_dei, high_dei;

	for(;;)
	{
		ti = (low_border + high_border)/2.0;

		if(fabs(low_border-high_border) < eps)
			return ti;

		ti_dei   = getei(ro, ti)-ei;
		low_dei  = getei(ro, low_border)-ei;
		high_dei = getei(ro, high_border)-ei;

		if(fabs(ti_dei) < eps)
			return ti;

		if(low_dei * high_dei > 0)
		{
			printf("\nsolve_ti: Error in initial data -- no equation root on the interval\n");
			printf("ro=%e, ei=%e, F(a)=%e, F(b)=%e\n", ro, ei, low_dei, high_dei);
			exit(1);
		}

		if(ti_dei * low_dei > 0)
			low_border  += (high_border-low_border)/2.0;
		else
			high_border -= (high_border-low_border)/2.0;
	}
}

double EOSAnalyticNi::getC(double ro, double ti, double te)
{	
	return /*sqrt( getdpdro(ro, ti, te) + 
			    (getpi(ro, ti)+getpe(ro, ti, te)) * dpde(ro, ti, te) / ro / ro +
				 getpi(ro, ti) * dpdei(ro, ti, te) / ro / ro );*/ 3000.;
}


double EOSAnalyticNi::getci(double ro, double ti)
{
	double NA    = 6.0e23;
	double kB    = 1.38e-16;
	double Matom = 58.7/NA;
	double nat   = ro/1000.0/Matom;
	double ci = 3.0*nat*kB;

	//return ci*1.0e-1;

	ci = ( getei(ro, 1.05*ti) - getei(ro, ti) ) / (0.05*ti) ;

	return ci;
}

double EOSAnalyticNi::getce(double ro, double te) {
	double ro0   = 8900.0;
	double K = 1077.0*te*te/2.0;
	double S = 1.275e5*pow(te, 1.3);
	double K_prime = 1077.0 * te;
	double S_prime = 1.6575e5*pow(te, 0.3);
	double ce = pow(ro/ro0, 0.6)*(K_prime*S*S + S_prime*K*K)/(K+S)/(K+S); // [J/K/m^3]
	return ce; // [J/K/m^3]
	// (если поделить на /ro, то будет [J/K/kg])
}

double EOSAnalyticNi::getkappa(double ro, double ti, double te) {
	return 1.0;
}


double EOSAnalyticNi::getAlpha(double ro, double ti, double te)
{
	double ro0 = 8900.0;
	
	// (test one) return -36.0e17 * 1.0e-1 * ro/ro0 ;
	
	// (actual one)     
	
	return -ro/ro0 * 3.72e19 / pow(te, 0.61); // [W/K/m^3]

	// (very old one)	return -ro/ro0 * 3.9e18 / pow(te, 0.25);

}

double EOSAnalyticNi::getEntropy(double ro, double ti, double te)
{
	return 0.0;
}

double EOSAnalyticNi::getGamma(void) {
	return 0.0;
}
/*
double EOSAnalyticNi::getdpdro(double ro, double ti, double te)
{
	return ( getpi(ro + EOS_EPS, ti) + getpe(ro + EOS_EPS, ti, te) 
		   - getpi(ro, ti)			 - getpe(ro, ti, te) ) / EOS_EPS;
}

double EOSAnalyticNi::getdpdroe(double ro, double ti, double te)
{
	return ( getpi(ro + EOS_EPS, ti) + getpe(ro + EOS_EPS, ti, te) 
		   - getpi(ro, ti)			 - getpe(ro, ti, te) ) / EOS_EPS;
}

double EOSAnalyticNi::getdpdroei(double ro, double ti, double te)
{
	return ( getpi(ro + EOS_EPS, ti) + getpe(ro + EOS_EPS, ti, te) 
		   - getpi(ro, ti)			 - getpe(ro, ti, te) ) / EOS_EPS;
}


double EOSAnalyticNi::getdpdro_rov_roE(double ro, double ti, double te, double v)
{ return 0;}

double EOSAnalyticNi::getdpdrov_ro_roE(double ro, double ti, double te, double v)
{ return 0;}

double EOSAnalyticNi::getdpdroE_ro_rov(double ro, double ti, double te, double v)
{ return 0;}

double EOSAnalyticNi::getdpdt(double ro, double ti, double te)
{
	return (getpi(ro, ti+EOS_EPS)-getpi(ro, ti))/EOS_EPS;
}

double EOSAnalyticNi::getdedt(double ro, double ti, double te)
{
	return (getei(ro, ti+EOS_EPS)-getei(ro, ti))/EOS_EPS;
}



////////////// PRIVATE ////////////////

double EOSAnalyticNi::dpdei(double ro, double ti, double te)
{
	return dpdti(ro, ti, te) / deidti(ro, ti, te) - 
		   dpdte(ro, ti, te) / dedte(ro, ti, te);
}


double EOSAnalyticNi::dpde(double ro, double ti, double te)
{
	return dpdte(ro, ti, te) / dedte(ro, ti, te);
}


double EOSAnalyticNi::dpdti(double ro, double ti, double te)
{
	return ( getpi(ro, ti + EOS_EPS) + getpe(ro, EOS_EPS, te) 
		   - getpi(ro, ti)           - getpe(ro, ti, te) ) / EOS_EPS;
}


double EOSAnalyticNi::dedti(double ro, double ti, double te)
{
	return ( getei(ro, ti + EOS_EPS) + getee(ro, ti + EOS_EPS, te) 
		   - getei(ro, ti)           - getee(ro, ti, te) ) / EOS_EPS;
}


double EOSAnalyticNi::deidti(double ro, double ti, double te)
{
	return ( getei(ro, ti + EOS_EPS) - getei(ro, ti) ) / EOS_EPS;
}


double EOSAnalyticNi::dpdte(double ro, double ti, double te)
{
	return ( getpe(ro, ti, te + EOS_EPS) - getpe(ro, ti, te) ) / EOS_EPS;

}


double EOSAnalyticNi::dedte(double ro, double ti, double te)
{
	return ( getee(ro, ti, te + EOS_EPS) - getee(ro, ti, te) ) / EOS_EPS;
}
*/

EOSAnalyticNi::EOSAnalyticNi()
{	
/*	MAX_T=100000.0; 
	MIN_T=300.0; 
	MAX_RO=5000.0; 
	MIN_RO=0.001;*/
	ro0 = 8900.;
}

double EOSAnalyticNi::getphase(double ro, double ti)
{
	return 1.0;
}


double EOSAnalyticNi::getmix(double ro, double ti)
{
	return 1.0;
}

/*double  EOSAnalyticNi::getnuWR(double ro, double ti, double te, double b)
{
	return 0.0;
}
*/


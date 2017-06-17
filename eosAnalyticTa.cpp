
#include"eosAnalyticTa.h"

#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define PI 3.14159

double EOSAnalyticTa::getpi(double ro, double ti)
{
	return 1.0;
}

double EOSAnalyticTa::getpe(double ro, double ti, double te, double Z=3.0)
{
	double pe = __pe(ro, te, Z) - __pe(ro, ti, Z);
	return pe;
}

double EOSAnalyticTa::__pe(double ro, double teta, double Z)
{
	double ro0   = 16654.0;
	double a	 = .05;
	double teKK  = teta/1000.;
	double at	 = atan((teKK-13.)/2.) / PI;
	double pe    = 1.09*pow(ro/ro0, a)* ( .09*teKK*teKK*(-at + .5) + (1.4*pow(teKK, 1.15) - 11.)*(at + .5) ) + .58259;// [GPa]
	return pe*1.e9; //[Pa]
}

double EOSAnalyticTa::getei(double ro, double ti) 
{
	return 1.0;	
}

double EOSAnalyticTa::getee(double ro, double ti, double te, double Z=3.0) 
{
	double ee = (__ee(ro, te, Z) - __ee(ro, ti, Z)) ; // [J/m^3]
	return ee/ro; // [J/kg]
}

double EOSAnalyticTa::__ee(double ro, double teta, double Z)
{
	double ro0  = 16654.0;
	double b    = .1;
	double teKK = teta/1000.;
	double ee = pow(ro/ro0, -b)*0.09*teKK*teKK/(1.+5.e-4*teKK*teKK)*(1.+1.5e-5*pow(teKK, 2.7))*1.e9; //[J/m3]
	return ee; // [J/m^3]
}

double EOSAnalyticTa::getti(double ro, double ei)
{
	double ti = solve_ti(ro, ei, 0.00001, 60000.0);
	return ti;
}

double EOSAnalyticTa::gette(double ro, double ti, double ee, double Z)
{
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

		te_dee   = getee(ro, ti, te, Z)-ee;
		low_dee  = getee(ro, ti, low_border, Z)-ee;
		high_dee = getee(ro, ti, high_border, Z)-ee;

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

double EOSAnalyticTa::solve_ti(double ro, double ei, double low_border, double high_border)
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
			return -1.;
		}

		if(ti_dei * low_dei > 0)
			low_border  += (high_border-low_border)/2.0;
		else
			high_border -= (high_border-low_border)/2.0;
	}
}

double EOSAnalyticTa::getC(double ro, double ti, double te, double v)
{	
	return sqrt( getdpdro(ro, ti, te) + 
			    (getpi(ro, ti)+getpe(ro, ti, te)) * dpde(ro, ti, te) / ro / ro +
				 getpi(ro, ti) * dpdei(ro, ti, te) / ro / ro );
}


double EOSAnalyticTa::getci(double ro, double ti)
{
	double NA    = 6.0e23;
	double kB    = 1.38e-23;
	double Matom = 181.0e-3/NA;
	double nat   = ro/Matom;
	double ci = 3.0*nat*kB;
	return ci; // [J/m3/K]
}

double EOSAnalyticTa::getce(double ro, double te, double Z)
{
	double ro0   = 16654.;
	double b	 = .1;
	double teKK  = te/1000.;
	double ce = 5.3* pow(ro/ro0, -b)*(0.5*teKK + 9.e-8*teKK*teKK*teKK*teKK*teKK)/(1.+3.5e-3*pow(teKK, 1.8)) * 1.e5; // [J/K/m^3]
	return ce; 
}

double EOSAnalyticTa::getkappa(double ro, double ti, double te, double Z)
{
	return 1.0;
}


double EOSAnalyticTa::getAlpha(double ro, double ti, double te)
{
	double ro0 = 16654.;
	return -1.27e17*ro/ro0; //[J/s/m3/K]
}

double EOSAnalyticTa::getEntropy(double ro, double ti, double te)
{
	return 0.0;
}

double EOSAnalyticTa::getGamma(void) {
	return 0.0;
}


double EOSAnalyticTa::getdpdro(double ro, double ti, double te)
{
	return ( getpi(ro + EOS_EPS, ti) + getpe(ro + EOS_EPS, ti, te) 
		   - getpi(ro, ti)			 - getpe(ro, ti, te) ) / EOS_EPS;
}

double EOSAnalyticTa::getdpdroe(double ro, double ti, double te)
{
	return ( getpi(ro + EOS_EPS, ti) + getpe(ro + EOS_EPS, ti, te) 
		   - getpi(ro, ti)			 - getpe(ro, ti, te) ) / EOS_EPS;
}

double EOSAnalyticTa::getdpdroei(double ro, double ti, double te)
{
	return ( getpi(ro + EOS_EPS, ti) + getpe(ro + EOS_EPS, ti, te) 
		   - getpi(ro, ti)			 - getpe(ro, ti, te) ) / EOS_EPS;
}

double EOSAnalyticTa::getdpdro_rov_roE(double ro, double ti, double te, double v)
{ return 0;}

double EOSAnalyticTa::getdpdrov_ro_roE(double ro, double ti, double te, double v)
{ return 0;}

double EOSAnalyticTa::getdpdroE_ro_rov(double ro, double ti, double te, double v)
{ return 0;}

double EOSAnalyticTa::getdpdt(double ro, double ti, double te)
{
	return (getpi(ro, ti+EOS_EPS)-getpi(ro, ti))/EOS_EPS;
}

double EOSAnalyticTa::getdedt(double ro, double ti, double te)
{
	return (getei(ro, ti+EOS_EPS)-getei(ro, ti))/EOS_EPS;
}

////////////// PRIVATE ////////////////

double EOSAnalyticTa::dpdei(double ro, double ti, double te)
{
	return dpdti(ro, ti, te) / deidti(ro, ti, te) - 
		   dpdte(ro, ti, te) / dedte(ro, ti, te);
}


double EOSAnalyticTa::dpde(double ro, double ti, double te)
{
	return dpdte(ro, ti, te) / dedte(ro, ti, te);
}


double EOSAnalyticTa::dpdti(double ro, double ti, double te)
{
	return ( getpi(ro, ti + EOS_EPS) + getpe(ro, EOS_EPS, te) 
		   - getpi(ro, ti)           - getpe(ro, ti, te) ) / EOS_EPS;
}


double EOSAnalyticTa::dedti(double ro, double ti, double te)
{
	return ( getei(ro, ti + EOS_EPS) + getee(ro, ti + EOS_EPS, te) 
		   - getei(ro, ti)           - getee(ro, ti, te) ) / EOS_EPS;
}


double EOSAnalyticTa::deidti(double ro, double ti, double te)
{
	return ( getei(ro, ti + EOS_EPS) - getei(ro, ti) ) / EOS_EPS;
}


double EOSAnalyticTa::dpdte(double ro, double ti, double te)
{
	return ( getpe(ro, ti, te + EOS_EPS) - getpe(ro, ti, te) ) / EOS_EPS;

}


double EOSAnalyticTa::dedte(double ro, double ti, double te)
{
	return ( getee(ro, ti, te + EOS_EPS) - getee(ro, ti, te) ) / EOS_EPS;
}


EOSAnalyticTa::EOSAnalyticTa()
{	
/*	MAX_T=100000.0; 
	MIN_T=300.0; 
	MAX_RO=5000.0; 
	MIN_RO=0.001;*/
	ro0 = 16654.;
}

double EOSAnalyticTa::getphase(double ro, double ti)
{
	return 1.0;
}


double EOSAnalyticTa::getmix(double ro, double ti)
{
	return 1.0;
}

double  EOSAnalyticTa::getnuWR(double ro, double ti, double te, double b, double Z)
{
	return 0.0;
}

complex<double> EOSAnalyticTa::geteps(double ro, double ti, double te, double Z=3.0)
{
	return complex<double>(0.0);
}


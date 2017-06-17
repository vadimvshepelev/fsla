
#include "eosIdeal.h"

#include <math.h>


double EOSIdeal::gette(double ro, double ti, double ee) { return 0.0; }
double EOSIdeal::getti(double ro, double e)				{ return (gamma-1.0)*e*M/R; }

double EOSIdeal::getpe(double ro, double ti, double te) { return 0.0; }
double EOSIdeal::getpi(double ro, double ti)			{ return ro*R*ti/M; }

double EOSIdeal::getee(double ro, double ti, double te) { return 0.0; }
double EOSIdeal::getei(double ro, double ti)			{ return R*ti/M/(gamma-1.0); }

double EOSIdeal::getce(double ro, double te)			{ return 0.0; }
double EOSIdeal::getci(double ro, double ti)			{ return R/M/(gamma-1.); } // [J/kg/K]

double EOSIdeal::getC(double ro, double ti, double te)	
{ 
	if (ro==0.0)
		return 0.0;
	else
		return sqrt(gamma*R*ti/M); 
}

double EOSIdeal::getAlpha(double ro, double ti, double te)	   { return 0.0; }
double EOSIdeal::getkappa(double ro, double ti, double te) { return 0.0; } 

double EOSIdeal::getphase(double ro, double ti) { return 5.0; }
double EOSIdeal::getmix(double ro, double ti) { return 1.0; }

double EOSIdeal::getGamma(void) {return gamma;}

//double EOSIdeal::getnuWR(double ro, double ti, double te, double b, double Z) {return 0.0; }

/*complex<double> EOSIdeal::geteps(double ro, double ti, double te, double Z) 
{ return complex<double>(0.0);}*/

/*
double EOSIdeal::getdpdro(double ro, double ti, double te)   {return R*ti/M;}
double EOSIdeal::getdpdroe(double ro, double ti, double te)  { return gamma-1.0; }
double EOSIdeal::getdpdroei(double ro, double ti, double te) { return 0.0; }

double EOSIdeal::getdpdt(double ro, double ti, double te){return ro*R/M;}
double EOSIdeal::getdedt(double ro, double ti, double te){return R/(gamma-1)/M;}

double EOSIdeal::getdpdro_rov_roE(double ro, double ti, double te, double v) {return 0.5*(gamma-1.0)*v*v;}
double EOSIdeal::getdpdrov_ro_roE(double ro, double ti, double te, double v) {return -(gamma-1.0)*v;}
double EOSIdeal::getdpdroE_ro_rov(double ro, double ti, double te, double v) {return (gamma-1.0);} */

double EOSIdeal::getdpdro(double ro, double ti, double te)
{
	return (getpi(ro+EOS_EPS, ti)-getpi(ro, ti))/EOS_EPS;	
}

double EOSIdeal::getdedro(double ro, double ti, double te)
{
	return (getei(ro+EOS_EPS, ti)-getei(ro, ti))/EOS_EPS;	
}

double EOSIdeal::getdpdt(double ro, double ti, double te)
{
	return (getpi(ro, ti+EOS_EPS)-getpi(ro, ti))/EOS_EPS;
}

double EOSIdeal::getdedt(double ro, double ti, double te)
{
	return (getei(ro, ti+EOS_EPS)-getei(ro, ti))/EOS_EPS;
}





/*double EOSIdeal::getdpdro_rov_roE(double ro, double ti, double te, double v)
{	
	double dpde = getdpdt(ro, ti, te)/getdedt(ro, ti, te); 
	return (getdpdro(ro, ti, te) + dpde/ro*(v*v/2.0 - getei(ro, ti) - ro*getdedro(ro, ti, te)));		
}

double EOSIdeal::getdpdrov_ro_roE(double ro, double ti, double te, double v) 
{
	return (- getdpdt(ro, ti, te)/getdedt(ro, ti, te)*v/ro);
}

double EOSIdeal::getdpdroE_ro_rov(double ro, double ti, double te, double v) 
{ 
	return (getdpdt(ro, ti, te)/getdedt(ro, ti, te)/ro);
}
*/

EOSIdeal::EOSIdeal(double _gamma) : gamma(_gamma) { R=8.31; M=27.0e-3; MAX_T=100000.0; MIN_T=0.001; MAX_RO=5000.0; MIN_RO=0.001;}



/*
	double dpdro  = 0.5*(5.0/3.0 - 1.0)*n.v*n.v;
	double dpdrov =    -(5.0/3.0 - 1.0)*n.v;
	double dpdroE =     (5.0/3.0 - 1.0);
*/
#ifndef _FEOSMGLIF_H_
#define _FEOSMGLIF_H_

#include "../FEOS.h"


class FEOSMGLiF : public FEOS {
	const double rho0, G, c, s;
	double pCold(double rho);
	double eCold(double rho);
	double pColdPrimeX(double rho);
	double eColdPrimeX(double rho);
public:
	FEOSMGLiF() : rho0(2650.), G(.71), c(5150.), s(1.35) {}
	double getp(double rho, double e);
	double gete(double rho, double p);
	double getc(double rho, double p);
	double gets(double rho, double p) { return 1.; }
	double getdpdrho(double rho, double e) { return 0.; }
	double getdpde(double rho, double e) { return 0.; }
	string gettype(void) { return string("mg");}
};

#endif // ! _FEOSMGLIF_H_


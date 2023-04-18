#include<math.h>

#include "FEOSMGLiF.h"


double FEOSMGLiF::pCold(double rho) { 
	double x = rho0 / rho; 
	return rho0 * c * c * (1. - x) / ((1. - s * (1. - x)) * (1. - s * (1. - x))); 
}

double FEOSMGLiF::eCold(double rho) { 
	double x = rho0 / rho; 
	return pCold(rho) / rho0 * (1. - x) / 2.; 
}

double FEOSMGLiF::pColdPrimeX(double rho) {
	double x = rho0 / rho;
	double denom = (1. - s * (1. - x)) * (1. - s * (1. - x)) * (1. - s * (1. - x));
	return -rho0 * c * c * (1. + s * (1. - x)) / denom;
}
double FEOSMGLiF::eColdPrimeX(double rho) {
	double x = rho0 / rho;
	return 1. / (2. * rho0) * ((1. - x) * pColdPrimeX(rho) - pCold(rho));
}

double FEOSMGLiF::getp(double rho, double e) {
	return pCold(rho) + G * rho * (e - eCold(rho)); 
}

double FEOSMGLiF::gete(double rho, double p) {
	return eCold(rho) + (p - pCold(rho)) / (rho * G); 
}

double FEOSMGLiF::getc(double rho, double p) {
	double pRho = -rho0 / (rho * rho) * pColdPrimeX(rho) - G * eCold(rho) + G * rho0 / rho * eColdPrimeX(rho);
	double pE = G * rho;
	return sqrt(pRho + p * pE / (rho * rho)); // [m/s]
}
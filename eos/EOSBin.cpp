#include "EOSBin.h"

double EOSBin::getp(double ro, double e) {
	return (gamma-1.)*ro*e + c0*c0*(ro-ro0);
}

double EOSBin::gete(double ro, double p) {
	return (p-c0*c0*(ro-ro0))/(gamma-1.)/ro;
}

double EOSBin::getC(double ro, double e) {
	double p = (gamma-1.)*ro*e + c0*c0*(ro-ro0);
	double p0 = 1./gamma*ro0*c0*c0; //????????????????????????????????? точно ли делить на гамму
	return sqrt(gamma*(p+p0)/ro);
}
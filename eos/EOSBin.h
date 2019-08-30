#ifndef _EOSBIN_H_
#define _EOSBIN_H_

#include "..\\eosold.h"

class EOSBin {
	// p = (gamma-1)*ro*e + ñ0^2*(ro-ro0)
	// e = ...(ro,T)
	//
public:
	const double gamma, ro0, c0;
	EOSBin();
	EOSBin(double _gamma, double _ro0, double _c0): gamma(_gamma), ro0(_ro0), c0(_c0) {}
	double getp(double ro, double e);
	double gete(double ro, double p);
    double getC(double ro, double e);	
};




#endif
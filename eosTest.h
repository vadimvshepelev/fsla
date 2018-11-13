#ifndef __EOSTEST_H_
#define __EOSTEST_H_

#include "eosAnalytic.h"

class EOSTest : public EOSAnalytic
{
public:

	EOSType getType() { return test; }

	double    getce(double ro, double te, double Z) { return 1.0; }
	double    getci(double ro, double ti) { return 1.0; }
	double    getkappa(double ro, double ti, double te) { return 1.0e-8; }
	double    getAlpha(double ro, double ti, double te) { return -1.0; };
	double    getEntropy(double ro, double T) { return 0.; }
};


#endif
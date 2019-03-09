
#ifndef EOSIDEAL_H
#define EOSIDEAL_H


#include "eosold.h"


// Уравнение состояния идеального газа

class EOSIdeal : public EOSOld
{

public:
	const double gamma;

	EOSIdeal(double _gamma);

	EOSType getType() { return ideal; }
	
	double    getpi(double ro, double ti);
	double    getpe(double ro, double ti, double te);
	
	double    getei(double ro, double ti);
	double    getee(double ro, double ti, double te);
		
	double    getti(double ro, double ei);
	double    gette(double ro, double ti, double ee);

	double    getci(double ro, double ti);
	double    getce(double ro, double te);

	double     getC(double ro, double ti, double te);
	double getkappa(double ro, double ti, double te);
	double     getAlpha(double ro, double ti, double te);

	double getphase(double ro, double ti);
	double   getmix(double ro, double ti);

	double getEntropy(double ro, double ti);
	double getGamma(void);

	double getdpdro  (double ro, double ti, double te);
	double getdedro  (double ro, double ti, double te);
	double getdpdt   (double ro, double ti, double te);
	double getdedt   (double ro, double ti, double te);
private:
	double R;				
	double M;
	
};



#endif

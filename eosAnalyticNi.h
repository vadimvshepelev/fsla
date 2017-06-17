
#ifndef EOSANALYTICNI_H
#define EOSANALYTICNI_H


#include "eos.h"


// јналитическое уравнение состо€ни€ дл€ никел€

class EOSAnalyticNi : public EOS {
public:
	EOSAnalyticNi(); 
	EOSType getType()  { return analyticNi; } 	
	double    getpi(double ro, double ti);
	double    getpe(double ro, double ti, double te);	
	double    getei(double ro, double ti);
	double    getee(double ro, double ti, double te);		
	double    getti(double ro, double ei);
	double    gette(double ro, double ti, double ee);
	double     getC(double ro, double ti, double te);
	double    getci(double ro, double ti);
	double    getce(double ro, double te);
	double getkappa(double ro, double ti, double te);
	double getAlpha(double ro, double ti, double te);
	double getphase(double ro, double ti);
	double   getmix(double ro, double ti);
	double getEntropy(double ro, double ti, double te);
	double getGamma(void);
	// Partial derivatives
/*	double getdpdro  (double ro, double ti, double te);
	double getdpdroe (double ro, double ti, double te);
	double getdpdroei(double ro, double ti, double te);
	double getdpdt   (double ro, double ti, double te);
	double getdedt   (double ro, double ti, double te);
	*/
protected:
	double __ee(double ro, double teta);
	double __pe(double ro, double teta);
	double solve_ti(double ro, double ei, double low_border, double high_border);
	// Auxilary partial derivatives
/*	double dpde(double ro, double ti, double te);
	double dpdei(double ro, double ti, double te);
	double dpdti(double ro, double ti, double te);
	double dedti(double ro, double ti, double te);
	double deidti(double ro, double ti, double te);
	double dpdte(double ro, double ti, double te);
	double dedte(double ro, double ti, double te);
	*/
};


#endif


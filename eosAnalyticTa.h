#ifndef EOSANALYTICTA_H
#define EOSANALYTICTA_H

#include "eosold.h"

// јналитическое уравнение состо€ни€ дл€ тантала

class EOSAnalyticTa : public EOSOld
{
public:
	EOSAnalyticTa();
	EOSType getType() { return analyticTa; }

	double    getpi(double ro, double ti);
	double    getpe(double ro, double ti, double te, double Z);
	
	double    getei(double ro, double ti);
	double    getee(double ro, double ti, double te, double Z);
		
	double    getti(double ro, double ei);
	double    gette(double ro, double ti, double ee, double Z);
	

	double     getC(double ro, double ti, double te, double v=0);

	double    getci(double ro, double ti);
	double    getce(double ro, double te, double Z=3.0);
	double getkappa(double ro, double ti, double te, double Z=3.0);
	double getAlpha(double ro, double ti, double te);

	double getphase(double ro, double ti);
	double   getmix(double ro, double ti);

	double  getnuWR(double ro, double ti, double te, double b, double Z=3.0);

	double getEntropy(double ro, double ti, double te);
	double getGamma(void);

	complex<double> geteps(double ro, double ti, double te, double Z);

	// Partial derivatives

	double getdpdro  (double ro, double ti, double te);
	double getdpdroe (double ro, double ti, double te);
	double getdpdroei(double ro, double ti, double te);

	double getdpdt   (double ro, double ti, double te);
	double getdedt   (double ro, double ti, double te);

	double getdpdro_rov_roE(double ro, double ti, double te, double v);
	double getdpdrov_ro_roE(double ro, double ti, double te, double v);
	double getdpdroE_ro_rov(double ro, double ti, double te, double v);

protected:

	//double z;


	double __ee(double ro, double teta, double Z=3.0);
	double __pe(double ro, double teta, double Z=3.0);
	double solve_ti(double ro, double ei, double low_border, double high_border);

	// Auxilary partial derivatives

	double dpde(double ro, double ti, double te);
	double dpdei(double ro, double ti, double te);
	double dpdti(double ro, double ti, double te);
	double dedti(double ro, double ti, double te);
	double deidti(double ro, double ti, double te);
	double dpdte(double ro, double ti, double te);
	double dedte(double ro, double ti, double te);

};


#endif


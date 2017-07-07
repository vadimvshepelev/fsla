#ifndef _EOSBIN_H_
#define _EOSBIN_H_


class EOSBin {
	// p = (gamma-1)*ro*e + C0^2*(ro-ro0)
	// e = ...(ro,T)
	//
	
public:
	const double gamma, ro0, C0;
	EOSBin();
	double getp(double ro, double e);
	double gete(double ro, double p);
    double getC(double ro, double e);	
	double getT(double ro, double e);	
    double getcv(double ro, double ti);
};






#endif
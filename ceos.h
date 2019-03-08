#ifndef _CEOS_H_
#define _EOS_H_

class CEOS {
public:
	virtual double getp(double ro, double e);
	virtual double gete(double ro, double p);
};


#endif


class CEOSIdeal : public CEOS {
public:
	const double gamma, R;
	CEOSIdeal(double _gamma) : gamma(_gamma), R(8.31) {}
	double getp(double ro, double e) { return (gamma-1.)*ro*e; }
	double gete(double ro, double p) { return p/(gamma-1.)/ro; }
};

class CEOSMieGruneisen : public CEOS {
public:
	const double ro0;
	CEOSMieGruneisen(double _ro0) : ro0(_ro0) {}
	double getp(double ro, double e);
	double gete(double ro, double p);
};



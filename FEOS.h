#ifndef _CEOS_H_
#define _CEOS_H_

#include<string>

using namespace std;

class FEOS {
public:
	virtual double getp(double ro, double e) = 0;
	virtual double gete(double ro, double p) = 0;
	virtual double getc(double ro, double p) = 0;
	virtual string gettype(void) = 0;
	virtual double getdpdrho(double rho, double e) = 0;
	virtual double getdpde(double rho, double e) = 0;
};


class FEOSIdeal : public FEOS {
public:
	const double gamma, R;
	FEOSIdeal(double _gamma) : gamma(_gamma), R(8.31) {}
	double getp(double ro, double e) { return (gamma-1.)*ro*e; }
	double gete(double ro, double p) { if(ro!=0.) return p/(gamma-1.)/ro; else return 0.; }
	double getc(double ro, double p) { if(ro!=0.) return sqrt(gamma*p/ro); else return 0.;}
	string gettype(void) {return string("ideal"); }
	double getdpdrho(double rho, double e) {return (gamma-1.)*e;}
	double getdpde(double rho, double e) {return (gamma-1.)*rho;}
};


class FEOSMieGruneisen : public FEOS {
public:
	const double ro0;
	FEOSMieGruneisen(): ro0(998.2) {}     // ro0 = [kg/m3]
	FEOSMieGruneisen(double _ro0, double _e0) : ro0(_ro0) {}
	double getp(double ro, double e);
	double gete(double ro, double p);
	double getc(double ro, double p);
	double getG(double ro);
	// Derivatives of G, p0, KS
	double getGPrime(double ro);
	double getp0Prime(double ro);
	double getKSPrime(double ro, double e);
	// Cold components, Born-Meyer potential 
	double getp0(double ro);
	double gete0(double ro);
	string gettype(void) {return string("mg"); }
	double getdpdrho(double rho, double e) {return 0.;}
	double getdpde(double rho, double e) {return 0.;}
};


class FEOSMieGruneisenAl : public FEOS {
public:
	const double ro0;
	FEOSMieGruneisenAl(): ro0(2700.) {}     // ro0 = [kg/m3]
	double getp(double ro, double e);
	double gete(double ro, double p);
	double getc(double ro, double p);
	string gettype(void) {return string("mg"); }
};


class FEOSMGAlPrecise : public FEOS {
	const double rho0, p0, a, b;
	double G(double x);
	double GPrime(double x);
	double Gx1(double x); 
	double Gx2(double x);
	double pCold(double rho);
	double eCold(double rho);
	double pColdPrime(double rho);
public:
	FEOSMGAlPrecise() : rho0(2750.), p0(560.964e9), a(1.12657), b(0.975511) {}
	double getp(double rho, double e);
	double gete(double rho, double p);
	double getc(double rho, double p);
	string gettype(void) {return string("mg"); }
};


class FEOSMGAlPrecise6 : public FEOS {
	const double rho0, p0, a, b, ph;
	double G(double x);
	double GPrime(double x);
	double Gx1(double x); 
	double Gx2(double x);
	double pCold(double rho);
	double pColdPrime(double rho);
	double eCold(double rho);
	double eColdPrime(double rho);
public:	
	FEOSMGAlPrecise6() : rho0(2750.), p0(560.964e9), a(1.12657), b(0.975511), ph(15.e9) {}
	double getp(double rho, double e);
	double gete(double rho, double p);
	double getc(double rho, double p);
	double gets(double rho, double p); // Entropy
	double getdpdrho(double rho, double e);
	double getdpde(double rho, double e);
	string gettype(void) {return string("mg"); }
};


// Constants set for 22 materials of Lomonosov (presumeably wide-range, not pure Mie-Gruneisen) EOS
extern double __V0[22], __E0[22], __DX[22], __GM[22], __CMN[22], __GN[22], __ES[22], __GC[22], __QS[22], __SM[22], __RS[22], 
	          __A1[22], __A2[22], __A3[22], __A4[22], __A5[22], __TA[22], __EA[22], __GI[22], __GF[22], __CR[22], __C1R[22];

class CEOSLomonosov {
public:
	double V0, E0, DX, GM, CMN, GN, ES, GC, QS, SM, RS, A1, A2, A3, A4, A5, TA, EA, GI, GF, CR, C1R;
	CEOSLomonosov(int _id);
	double getp(double ro, double e);
	double gete(double ro, double p);
	double getc(double ro, double e);
	string gettype(void) {return string("lomonosov"); }
private:
	void EOSE5(double ro, double e, double& P, double& C, double& g, bool& nonPhysical);

};


#endif

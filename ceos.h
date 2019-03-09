#ifndef _CEOS_H_
#define _EOS_H_

class CEOS {
public:
	virtual double getp(double ro, double e) = 0;
	virtual double gete(double ro, double p) = 0;
	virtual double getc(double ro, double e) = 0;
};


#endif


class CEOSIdeal : public CEOS {
public:
	const double gamma, R;
	CEOSIdeal(double _gamma) : gamma(_gamma), R(8.31) {}
	double getp(double ro, double e) { return (gamma-1.)*ro*e; }
	double gete(double ro, double p) { return p/(gamma-1.)/ro; }
	double getc(double ro, double e) { return sqrt(gamma*(gamma-1)*e);}
};


class CEOSMieGruneisen : public CEOS {
public:
	const double ro0;
	//CEOSMieGruneisen(): ro0(0.) {}
	CEOSMieGruneisen(double _ro0) : ro0(_ro0) {}
	double getp(double ro, double e);
	double gete(double ro, double p);
	double getc(double ro, double e);
};


// Constants set for 22 materials of Lomonosov (presumeably wide-range, not pure Mie-Gruneisen) EOS
extern double __V0[22], __E0[22], __DX[22], __GM[22], __CMN[22], __GN[22], __ES[22], __GC[22], __QS[22], __SM[22], __RS[22], 
	          __A1[22], __A2[22], __A3[22], __A4[22], __A5[22], __TA[22], __EA[22], __GI[22], __GF[22], __CR[22], __C1R[22];

class CEOSLomonosov : public CEOSMieGruneisen {
public:
	double V0, E0, DX, GM, CMN, GN, ES, GC, QS, SM, RS, A1, A2, A3, A4, A5, TA, EA, GI, GF, CR, C1R;
	CEOSLomonosov(int _id);
	double getp(double ro, double e);
	double gete(double ro, double p);
	double getc(double ro, double e);
private:
	void EOSE5(double ro, double e, double& P, double& C, double& g, bool& nonPhysical);

};
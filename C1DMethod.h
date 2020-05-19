#ifndef _C1DMETHOD_H_
#define _C1DMETHOD_H_

#include "_vector4.h"
#include "C1Dfield.h"
#include "F1DReconstruction.h"


enum RPWaveConfig {nothing, swrw, rwsw, swsw, rwrw, vacrw, rwvac, rwvacrw};


// Structure for the Riemann problem solution vector of primitive variables -- (roL, roR, u ,p)
struct RPValues {
	RPValues() : roL(0.), roR(0.), u(0.), p(0.), type(RPWaveConfig::nothing) {}
	double roL, roR, u, p;
	RPWaveConfig type;
};


struct C1DVectorPrimitive {
	double ro;
	double u;
	double p;
};

Vector4 calcPhysicalFlux(CEOS& eos, double ro, double u, double p);

class CRiemannSolver {
public: 
	virtual Vector4 calcFlux(CEOS& eos, double roL, double rouL, double roEL, double roR, double rouR, double roER) = 0;
	virtual int isSupported(CEOSIdeal& eos) = 0;
	virtual int isSupported(CEOSMieGruneisen& eos) = 0;
};


class CExactRiemannSolver : public CRiemannSolver {
public:
	CExactRiemannSolver() {}
	Vector4 calcFlux(CEOS& eos, double roL, double rouL, double roEL, double roR, double rouR, double roER);
	int isSupported(CEOSIdeal& eos) { return 1; }
	int isSupported(CEOSMieGruneisen& eos) {return 0; }
private:
	C1DVectorPrimitive calcSolution(CEOS& eos, double roL, double uL, double pL, double roR, double uR, double pR, double x, double t);
	RPValues calcValues(CEOS& eos, double roL, double uL, double pL, double roR, double uR, double pR);
	double fL(CEOS& eos, double p, double roL, double uL, double pL);
	double dfLdp(CEOS& eos, double p, double roL, double uL, double pL);
	double fR(CEOS& eos, double p, double roR, double uR, double pR);
	double dfRdp(CEOS& eos, double p, double roR, double uR, double pR);
};


class CHLLRiemannSolver : public CRiemannSolver {
public:
	CHLLRiemannSolver() {}
	Vector4 calcFlux(CEOS& eos, double roL, double rouL, double roEL, double roR, double rouR, double roER);
	int isSupported(CEOSIdeal& eos) { return 1; }
	int isSupported(CEOSMieGruneisen& eos) {return 1; }
};


class CHLLCRiemannSolver : public CRiemannSolver {
public:
	CHLLCRiemannSolver() {}
	Vector4 calcFlux(CEOS& eos, double roL, double rouL, double roEL, double roR, double rouR, double roER);
	int isSupported(CEOSIdeal& eos) { return 1; }
	int isSupported(CEOSMieGruneisen& eos) {return 1; }
};


// Godunov-Prokhorov-Safronov simple linearized entropy non-decreasing Riemann solver
class CGPSRiemannSolver : public CRiemannSolver {
public: 
	CGPSRiemannSolver() {}
	Vector4 calcFlux(CEOS& eos, double roL, double rouL, double roEL, double roR, double rouR, double roER);
	int isSupported(CEOSIdeal& eos) { return 1; }
	int isSupported(CEOSMieGruneisen& eos) {return 0; }
};


class C1DMethod {
public:
	virtual void calc(C1DProblem& pr, CEOS& eos, C1DField& fld)=0;
	virtual double calcdt(C1DProblem& pr, CEOS& eos, C1DField& fld)=0;	
};


class C1DGodunovMethodMillerPuckett : public C1DMethod {
public:
	void calc(C1DProblem& pr, CEOSMieGruneisen& eos, C1DField& fld);
	double calcdt(C1DProblem& pr, CEOSMieGruneisen& eos, C1DField& fld);
	C1DVectorPrimitive calcRPExactMillerPuckett(CEOSMieGruneisen& eos, double roL, double vL, double pL, double roR, double vR, double pR);
};


class C1DGodunovTypeMethod : public C1DMethod {
public:
	CRiemannSolver& rslv;
	C1DGodunovTypeMethod();  
	C1DGodunovTypeMethod(CRiemannSolver& _rslv) : rslv(_rslv) {}
	void calc(C1DProblem& pr, CEOS& eos, C1DField& fld);
	double calcdt(C1DProblem& pr, CEOS& eos, C1DField& fld);	
};


class C1DGodunovTypeMethodVacuum : public C1DMethod {
public:
	CRiemannSolver& rslv;
	C1DGodunovTypeMethodVacuum();  
	C1DGodunovTypeMethodVacuum(CRiemannSolver& _rslv, double _xbnd) : rslv(_rslv), xbnd(_xbnd) {}
	double xbnd;
	void calc(C1DProblem& pr, CEOS& eos, C1DField& fld);
	double calcdt(C1DProblem& pr, CEOS& eos, C1DField& fld);	
};

class C1D2ndOrderMethod : public C1DGodunovTypeMethod {
public:
	C1D2ndOrderMethod(CRiemannSolver& _rslv, F1DReconstruction& _rec) : C1DGodunovTypeMethod(_rslv), rec(_rec) {} 
	F1DReconstruction& rec;
	void calc(C1DProblem& pr, CEOS& eos, C1DField& fld);
};


#endif
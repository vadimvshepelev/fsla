#ifndef _C1DMETHOD_H_
#define _C1DMETHOD_H_

#include"C1Dfield.h"
#include"solver.h"


class CRiemannSolver {
public: 
	virtual Vector4 calcFlux(CEOS& eos, double roL, double rouL, double roEL, double roR, double rouR, double roER) = 0;
};


class CExactRiemannSolver : public CRiemannSolver {
public:
	CExactRiemannSolver();
	Vector4 calcFlux(CEOS& eos, double roL, double rouL, double roEL, double roR, double rouR, double roER);
};


class CHLLRiemannSolver : public CRiemannSolver {
public:
	CHLLRiemannSolver() {}
	Vector4 calcFlux(CEOS& eos, double roL, double rouL, double roEL, double roR, double rouR, double roER);
};


class CHLLCRiemannSolver : public CRiemannSolver {
public:
	CHLLCRiemannSolver() {}
	Vector4 calcFlux(CEOS& eos, double roL, double rouL, double roEL, double roR, double rouR, double roER);
};


// Godunov-Prokhorov-Safronov simple linearized entropy non-decreasing Riemann solver
class CGPSRiemannSolver : public CRiemannSolver {
public: 
	CGPSRiemannSolver() {}
	Vector4 calcFlux(CEOS& eos, double roL, double rouL, double roEL, double roR, double rouR, double roER);
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
	CVectorPrimitive calcRPExactMillerPuckett(CEOSMieGruneisen& eos, double roL, double vL, double pL, double roR, double vR, double pR);
};


class C1DGodunovTypeMethod : public C1DMethod {
public:
	CRiemannSolver& rslv;
	C1DGodunovTypeMethod();  
	C1DGodunovTypeMethod(CRiemannSolver& _rslv) : rslv(_rslv) {}
	void calc(C1DProblem& pr, CEOS& eos, C1DField& fld);
	double calcdt(C1DProblem& pr, CEOS& eos, C1DField& fld);	
};


#endif
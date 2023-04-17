#ifndef _C1DMETHOD_H_
#define _C1DMETHOD_H_

#include "_vector4.h"
#include "_matrix4.h"
#include "C1DField.h"
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

Vector4 calcPhysicalFlux(FEOS& eos, double ro, double u, double p);

class CRiemannSolver {
public: 
	virtual Vector4 calcFlux(FEOS& eos, double roL, double rouL, double roEL, double roR, double rouR, double roER) = 0;
	virtual int isSupported(FEOSIdeal& eos) = 0;
	virtual int isSupported(FEOSMieGruneisen& eos) = 0;
};


class CExactRiemannSolver : public CRiemannSolver {
public:
	CExactRiemannSolver() {}
	Vector4 calcFlux(FEOS& eos, double roL, double rouL, double roEL, double roR, double rouR, double roER);
	int isSupported(FEOSIdeal& eos) { return 1; }
	int isSupported(FEOSMieGruneisen& eos) {return 0; }
};

extern CExactRiemannSolver exrslv;


C1DVectorPrimitive calcSolution(FEOS& eos, double roL, double uL, double pL, double roR, double uR, double pR, double x, double t);
RPValues calcValues(FEOS& eos, double roL, double uL, double pL, double roR, double uR, double pR);
double fL(FEOS& eos, double p, double roL, double uL, double pL);
double dfLdp(FEOS& eos, double p, double roL, double uL, double pL);
double fR(FEOS& eos, double p, double roR, double uR, double pR);
double dfRdp(FEOS& eos, double p, double roR, double uR, double pR);



class CHLLRiemannSolver : public CRiemannSolver {
public:
	CHLLRiemannSolver() {}
	Vector4 calcFlux(FEOS& eos, double roL, double rouL, double roEL, double roR, double rouR, double roER);
	int isSupported(FEOSIdeal& eos) { return 1; }
	int isSupported(FEOSMieGruneisen& eos) {return 1; }
};


class CHLLCRiemannSolver : public CRiemannSolver {
public:
	CHLLCRiemannSolver() {}
	Vector4 calcFlux(FEOS& eos, double roL, double rouL, double roEL, double roR, double rouR, double roER);
	int isSupported(FEOSIdeal& eos) { return 1; }
	int isSupported(FEOSMieGruneisen& eos) {return 1; }
};

// Godunov-Prokhorov-Safronov simple linearized entropy non-decreasing Riemann solver
class CGPSRiemannSolver : public CRiemannSolver {
public: 
	CGPSRiemannSolver() {}
	Vector4 calcFlux(FEOS& eos, double roL, double rouL, double roEL, double roR, double rouR, double roER);
	int isSupported(FEOSIdeal& eos) { return 1; }
	int isSupported(FEOSMieGruneisen& eos) {return 0; }
};

// Lax-Friedrichs Riemann solver
class CLFRiemannSolver : public CRiemannSolver {
public: 
	CLFRiemannSolver() {}
	Vector4 calcFlux(
			FEOS& eos,
			double roL, double rouL, double roEL,
			double roR, double rouR, double roER,
			double dx, double dt);
	Vector4 calcFlux(
			FEOS& eos,
			double roL, double rouL, double roEL,
			double roR, double rouR, double roER) { return Vector4::ZERO; }
	int isSupported(FEOSIdeal& eos) { return 1; }
	int isSupported(FEOSMieGruneisen& eos) {return 1; }
};


// Global Lax-Friedrichs Riemann solver
class CLFGlobalRiemannSolver : public CRiemannSolver {
public:
	CLFGlobalRiemannSolver() {}
	Vector4 calcFlux(
			FEOS& eos,
			double roL, double rouL, double roEL,
			double roR, double rouR, double roER,
			double lambda);
	Vector4 calcFlux(
			FEOS& eos,
			double roL, double rouL, double roEL,
			double roR, double rouR, double roER) { return Vector4::ZERO; }
	int isSupported(FEOSIdeal& eos) { return 1; }
	int isSupported(FEOSMieGruneisen& eos) {return 1; }
};


// BGK Riemann solver including 2nd order monotone reconstruction 
class CBGKRiemannSolver : public CRiemannSolver {
	Vector4 La, Lb, Lg, Ld;
	Matrix4 getOmegaInv(FEOS& eos, Vector4 U);
	Matrix4 getOmega(FEOS& eos, Vector4 U);
	void fillLambda(Vector4 Fm, Vector4 F, Vector4 Fp, Vector4 Fpp, Vector4 L, double step);
	void fillLambdaComponent(int i, double lambda, double criteria, double curant);
public: 
	CBGKRiemannSolver() {}
	Vector4 calcFlux(FEOS& eos, Vector4 Um, Vector4 U, Vector4 Up, Vector4 Upp, double dx, double dt);
	Vector4 calcFlux(FEOS& eos, double roL, double rouL, double roEL, double roR, double rouR, double roER) { return Vector4::ZERO; }
	int isSupported(FEOSIdeal& eos) { return 1; }
	int isSupported(FEOSMieGruneisen& eos) {return 1; }
};


class CRoeRiemannSolver : public CRiemannSolver {
public: 
	CRoeRiemannSolver() {}
	Vector4 calcFlux(FEOS& eos, double roL, double rouL, double roEL, double roR, double rouR, double roER);
	int isSupported(FEOSIdeal& eos) { return 1; }
	int isSupported(FEOSMieGruneisen& eos) {return 0; }
};


class CRoeGeneralRiemannSolver : public CRiemannSolver {
public: 
	CRoeGeneralRiemannSolver() {}
	Vector4 calcFlux(FEOS& eos, double roL, double rouL, double roEL, double roR, double rouR, double roER);
	int isSupported(FEOSIdeal& eos) { return 1; }
	int isSupported(FEOSMieGruneisen& eos) {return 1; }
};


class C1DMethodSamarskii {
public:
	void calc(C1DProblem& pr, FEOS& eos, C1DFieldPrimitive& fld);
	double calcdt(C1DProblem& pr, FEOS& eos, C1DFieldPrimitive& fld);
};


class C1DMethod {
public:
	virtual void calc(C1DProblem& pr, FEOS& eos, C1DField& fld)=0;
	virtual double calcdt(C1DProblem& pr, FEOS& eos, C1DField& fld)=0;
	virtual void calcFluxField(C1DProblem& pr, FEOS& eos, C1DField& fld)=0;
};


class C1DGodunovMethodMillerPuckett : public C1DMethod {
public:
	void calc(C1DProblem& pr, FEOSMieGruneisen& eos, C1DField& fld);
	double calcdt(C1DProblem& pr, FEOSMieGruneisen& eos, C1DField& fld);
	C1DVectorPrimitive calcRPExactMillerPuckett(FEOSMieGruneisen& eos, double roL, double vL, double pL, double roR, double vR, double pR);
};


class C1DGodunovTypeMethod : public C1DMethod {
public:
	CRiemannSolver& rslv;
	C1DGodunovTypeMethod() : rslv(exrslv) {}
	C1DGodunovTypeMethod(CRiemannSolver& _rslv) : rslv(_rslv) {}
	virtual void calc(C1DProblem& pr, FEOS& eos, C1DField& fld) override;
	double calcdt(C1DProblem& pr, FEOS& eos, C1DField& fld) override;
	virtual void calcFluxField(C1DProblem& pr, FEOS& eos, C1DField& fld) override;
};


class C1DLFMethod : public C1DGodunovTypeMethod {
	CLFRiemannSolver& lfrslv;
public:
	C1DLFMethod();  
	C1DLFMethod(CLFRiemannSolver& _lfrslv) : lfrslv(_lfrslv) {}
	void calc(C1DProblem& pr, FEOS& eos, C1DField& fld);
};


class C1DLFGlobalMethod : public C1DGodunovTypeMethod {
	CLFGlobalRiemannSolver& lfrslv;
public:
	C1DLFGlobalMethod();
	C1DLFGlobalMethod(CLFGlobalRiemannSolver& _lfrslv) : lfrslv(_lfrslv) {}
	void calc(C1DProblem& pr, FEOS& eos, C1DField& fld) override;
};



class C1DGodunovTypeMethodVacuum : public C1DMethod {
public:
	CRiemannSolver& rslv;
	C1DGodunovTypeMethodVacuum();  
	C1DGodunovTypeMethodVacuum(CRiemannSolver& _rslv, double _xbnd) : rslv(_rslv), xbnd(_xbnd) {}
	double xbnd;
	void calc(C1DProblem& pr, FEOS& eos, C1DField& fld);
	double calcdt(C1DProblem& pr, FEOS& eos, C1DField& fld);	
};

class C1D2ndOrderMethod : public C1DGodunovTypeMethod {
public:
	C1D2ndOrderMethod(CRiemannSolver& _rslv, F1DReconstruction& _rec)
		: C1DGodunovTypeMethod(_rslv), rec(_rec) {}

	F1DReconstruction& rec;
	virtual void calc(C1DProblem& pr, FEOS& eos, C1DField& fld) override;
	virtual void calcFluxField(C1DProblem& pr, FEOS& eos, C1DField& fld) override;
};


class C1D2ndOrderLFGlobalMethod : public C1D2ndOrderMethod {
public:
	CLFGlobalRiemannSolver& lfrslv;
	// F1DReconstruction& rec;

	C1D2ndOrderLFGlobalMethod(CLFGlobalRiemannSolver& _rslv,
							  F1DReconstruction& _rec)
		: C1D2ndOrderMethod(_rslv, _rec), lfrslv(_rslv) {}

//	void calc(C1DProblem& pr, FEOS& eos, C1DField& fld) override;
	void calcFluxField(C1DProblem& pr, FEOS& eos, C1DField& fld) override;
};


class C1DBGKMethod : public C1DGodunovTypeMethod {
	CBGKRiemannSolver &bgk;	
public:
	C1DBGKMethod();
	C1DBGKMethod(CBGKRiemannSolver& _bgk): bgk(_bgk) {}
	void calc(C1DProblem& pr, FEOS& eos, C1DField& fld);
};


#endif

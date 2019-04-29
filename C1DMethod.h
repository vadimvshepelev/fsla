#ifndef _C1DMETHOD_H_
#define _C1DMETHOD_H_

#include"C1Dfield.h"
#include"solver.h"


class C1DMethod {
public:
	virtual void calc(C1DProblem& pr, CEOSMieGruneisen& eos, C1DField& fld)=0;
	virtual double calcdt(C1DProblem& pr, CEOSMieGruneisen& eos, C1DField& fld)=0;
};


class C1DGodunovMethodMillerPuckett : public C1DMethod {
	void calc(C1DProblem& pr, CEOSMieGruneisen& eos, C1DField& fld);
	double calcdt(C1DProblem& pr, CEOSMieGruneisen& eos, C1DField& fld);
	CVectorPrimitive calcRPExactMillerPuckett(CEOSMieGruneisen& eos, double roL, double vL, double pL, double roR, double vR, double pR);
};



#endif
#ifndef _C1DSIMULATION_H_
#define _C1DSIMULATION_H_

#include"c1dproblem.h"
#include"ceos.h"
#include"c1dfield.h"
#include"c1dmethod.h"
#include"coutput.h"


class C1DSimulation {
	C1DProblem &pr;
	CEOS &eos;
	C1DField &fld;	
	C1DMethod &mtd;	
	COutput &outp;
public:
	C1DSimulation(C1DProblem& _pr, CEOS& _eos, C1DField& _fld, C1DMethod& _mtd, COutput& _outp) : pr(_pr), eos(_eos), fld(_fld), mtd(_mtd), outp(_outp) {}
	void run();
};

#endif
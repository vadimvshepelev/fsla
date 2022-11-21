#ifndef _F1DSIMULATION_H_
#define _F1DSIMULATION_H_

#include "C1DProblem.h"
#include "FEOS.h"
#include "C1DField.h"
#include "C1DMethod.h"
#include "COutput.h"
#include "odesolver.h"


class F1DSimulation {
	C1DProblem& pr;
	FEOS& eos;
	C1DField& fld;
	C1DMethod& mtd;
	COutput& outp;
	ODESolver& solve;
public:
	F1DSimulation(
			C1DProblem& _pr,
			FEOS& _eos,
			C1DField& _fld,
			C1DMethod& _mtd,
			COutput& _outp,
			ODESolver& _odemtd)
		: pr(_pr), eos(_eos), fld(_fld),
		  mtd(_mtd), outp(_outp), solve(_odemtd) {
	}

	void run();
};

#endif

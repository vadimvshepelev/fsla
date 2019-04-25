#ifndef _C1DSIMULATION_H_
#define _C1DSIMULATION_H_

#include "C1DProblem.h"
#include "ceos.h"
#include "C1DField.h"
#include "C1DMethod.h"

class C1DSimulation {
	C1DProblem &pr;
	CEOS &eos;
	C1DField &fld;
	C1DMethod &mtd;
public:
	void run();

};


#endif
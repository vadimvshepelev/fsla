#ifndef _F1DRECONSTRUCTION_H_
#define _F1DRECONSTRUCTION_H_

#include"C1DProblem.h"
#include"C1DField.h"
#include<vector>

using namespace std;

class F1DReconstruction {
public:
	vector<vector<double>> ULx, URx;
	F1DReconstruction(C1DField& fld);
	virtual void calc(C1DField& fld)=0;
};


class F1DENO2Reconstruction : public F1DReconstruction {
	vector<double> rm1r, r0r, rm1l, r0l;
public:
	F1DENO2Reconstruction(C1DField& fld);
	void calc(C1DField& fld);
};



#endif
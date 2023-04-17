#ifndef _C1DFIELD_H_
#define _C1DFIELD_H_

#include <vector>
#include "C1DProblem.h"
#include "FEOS.h"

using namespace std;


class C1DField {
public:
	const int imin, imax;
	vector<double> x;
	const double dx;
	double t, dt;
	vector<Vector4> U, newU;
	vector<Vector4> F;
	C1DField(C1DProblem& pr);
	~C1DField();
};


class C1DFieldPrimitive {
public:
	const int imin, imax;
	vector<double> x, newx;
	const double dm;
	double t, dt;
	vector<Vector4> W, newW, prevW;
	C1DFieldPrimitive(C1DProblem& pr);
	~C1DFieldPrimitive();
};


#endif

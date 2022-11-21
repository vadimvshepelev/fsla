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
//	vector<vector<double>> U, newU, F; 	// TODO: попробовать вариант с double ** или vector<CVector4> и сравнить по производительности
	vector<Vector4> U;
	vector<Vector4> newU;
	vector<Vector4> F;

	C1DField(C1DProblem& pr);
	~C1DField();
};


#endif

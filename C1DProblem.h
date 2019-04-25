#ifndef _C1DPROBLEM_H_
#define _C1DPROBLEM_H_

#include<string>

using namespace std;


class C1DProblem {	
public:
	C1DProblem(string _name, 
		       double _roL, double _uL, double _pL, double _roR, double _uR, double _pR, 
			   double _xMin, double _xMax, double _tMin, double _tMax, double _x0,
			   int _NX, double _CFL) : roL(_roL), uL(_uL), pL(_pL), roR(_roR), uR(_uR), pR(_pR),
			   xMin(_xMin), xMax(_xMax), tMin(_tMin), tMax(_tMax), x0(_x0), NX(_NX), CFL(_CFL) {}
	string name;
	const double roL, uL, pL, roR, uR, pR;
	const double xMin, xMax, tMin, tMax, x0;
	const int NX;
	double CFL;	
};


extern C1DProblem prNBtest;

#endif
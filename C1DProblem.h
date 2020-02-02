#ifndef _C1DPROBLEM_H_
#define _C1DPROBLEM_H_

#include<string>

using namespace std;


class C1DProblem {	
public:
//	C1DProblem() : name(0), rol(0.), ul(0.), pl(0.), ror(0.), ur(0.), pr(0.),
//			   xmin(0.), xmax(0.), tmin(0.), tmax(0.), x0(0.), nx(0), cfl(0.), bcs(0) {}
	C1DProblem(string _name, 
		       double _rol, double _ul, double _pl, double _ror, double _ur, double _pr, 
			   double _xmin, double _xmax, double _tmin, double _tmax, double _x0,
			   int _nx, double _cfl, string _bcs) : name(_name), rol(_rol), ul(_ul), pl(_pl), ror(_ror), ur(_ur), pr(_pr),
			   xmin(_xmin), xmax(_xmax), tmin(_tmin), tmax(_tmax), x0(_x0), nx(_nx), cfl(_cfl), bcs(_bcs) {}
	string name;
	const double rol, ul, pl, ror, ur, pr;
	const double xmin, xmax, tmin, tmax, x0;
	const int nx;
	double cfl;	
	string bcs;
	void setics(CEOS& eos, vector<vector<Vector4>> U);
	void setbcs(vector<vector<Vector4>> U);
};

// !!! А вот бы к каждому типу задачи (наследование) приделать функцию setICs, берущую как параметры EOS и CField, 
// и проставляющую начальные условия, как ей надо !!!

extern C1DProblem prNBtest;
extern C1DProblem prToro1Idealtest, prToro2Idealtest, prToro3Idealtest, prToro4Idealtest, prToro5Idealtest;
extern C1DProblem prDenisenko1, prDenisenko2, prDenisenko3;

#endif
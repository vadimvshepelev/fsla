#ifndef _C1DPROBLEM_H_
#define _C1DPROBLEM_H_

#include<string>
#include<vector>

#include"_vector4.h"
#include"FEOS.h"

using namespace std;


class C1DProblem {
public:
//	C1DProblem()
//		: name(0),
//		  rol(0.), ul(0.), pl(0.),
//		  ror(0.), ur(0.), pr(0.),
//		  xmin(0.), xmax(0.),
//		  tmin(0.), tmax(0.), x0(0.),
//		  nx(0), cfl(0.), bcs(0), order(1) {}

	C1DProblem(string _name,
			   double _rhol, double _ul, double _pl,
			   double _rhor, double _ur, double _pr,
			   double _xmin, double _xmax,
			   double _tmin, double _tmax, double _x0,
			   int _nx, double _cfl, string _bcs, int _order)
		: name(_name),
		  rol(_rhol), ul(_ul), pl(_pl),
		  ror(_rhor), ur(_ur), pr(_pr),
		  xmin(_xmin), xmax(_xmax),
		  tmin(_tmin), tmax(_tmax), x0(_x0),
		  nx(_nx), cfl(_cfl), bcs(_bcs), order(_order) {}

	// Default order = 2 to support old ctors
	C1DProblem(string _name,
			   double _rhol, double _ul, double _pl,
			   double _rhor, double _ur, double _pr,
			   double _xmin, double _xmax,
			   double _tmin, double _tmax, double _x0,
			   int _nx, double _cfl, string _bcs)
		: 	C1DProblem(_name,
					   _rhol, _ul, _pl,
					   _rhor, _ur, _pr,
					   _xmin, _xmax,
					   _tmin, _tmax, _x0,
					   _nx, _cfl, _bcs, 2) {}

	string name;
	const double rol, ul, pl, ror, ur, pr;
	const double xmin, xmax, tmin, tmax, x0;
	const int nx;
	double cfl;
	string bcs;
	const int order;
//	virtual void setics(FEOS& eos, vector<double>& x,
//						vector<vector<double>>& U);
	virtual void setics(FEOS& eos, vector<double>& x, vector<Vector4>& U);
//	 void setbcs(vector<vector<double>>& U);
	void setbcs(vector<Vector4>& U);
	constexpr int get_order() { return order; }
};

class C1DLaserProblem : public C1DProblem {
public:
	C1DLaserProblem(string _name,
					double _rhol, double _ul, double _pl,
					double _rhom, double _um, double _pm,
					double _rhor, double _ur, double _pr,
					double _xmin, double _xmax, double _tmin,
					double _tmax, double _x0, double _x1,
					int _nx, double _cfl, string _bcs, int _order)
		: C1DProblem(_name,
					 _rhol, _ul, _pl,
					 _rhor, _ur, _pr,
					 _xmin, _xmax,
					 _tmin, _tmax, _x0,
					 _nx, _cfl, _bcs, _order),
		  rom(_rhom), um(_um), pm(_pm), x1(_x1) {}

	C1DLaserProblem(string _name,
					double _rhol, double _ul, double _pl,
					double _rhom, double _um, double _pm,
					double _rhor, double _ur, double _pr,
					double _xmin, double _xmax, double _tmin,
					double _tmax, double _x0, double _x1,
					int _nx, double _cfl, string _bcs)
		: C1DProblem(_name,
					 _rhol, _ul, _pl,
					 _rhor, _ur, _pr,
					 _xmin, _xmax,
					 _tmin, _tmax, _x0,
					 _nx, _cfl, _bcs, 2),
		  rom(_rhom), um(_um), pm(_pm), x1(_x1) {}

	const double rom, um, pm;
	const double x1;
//	void setics(FEOS& eos, vector<double>& x,
//				vector<vector<double>>& U);
	void setics(FEOS& eos, vector<double>& x, vector<Vector4>& U);
};


extern C1DProblem prNBtest;
extern C1DProblem prToro1Idealtest;
extern C1DProblem prToro2Idealtest;
extern C1DProblem prToro3Idealtest;
extern C1DProblem prToro4Idealtest;
extern C1DProblem prToro5Idealtest;

extern C1DProblem prDenisenko1;
extern C1DProblem prDenisenko2;
extern C1DProblem prDenisenko3;

extern C1DProblem prLaserVTAlIdealTest1;
extern C1DProblem prLaserVTAlIdealTest2;

extern C1DProblem prLaserVTAlMGTest1;
extern C1DProblem prLaserVTAlMGTestNum1;
extern C1DProblem prLaserVTAlMGTest2;

extern C1DProblem prIdealVacTest;

extern C1DProblem prVTAlMGTest1;
extern C1DProblem prVTAlMGTest2;
extern C1DProblem prVTAlMGTest2_2;

extern C1DProblem prDenisenko;

extern C1DLaserProblem prHoles;

extern C1DProblem prFedorAl;

#endif

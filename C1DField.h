#ifndef _C1DFIELD_H_
#define _C1DFIELD_H_

#include<vector>
#include"C1DProblem.h"
#include"ceos.h"

using namespace std;


class C1DField {
public:
	const int iMin, iMax;
	vector<double> x;
	const double dx;
	double t, dt;
	vector<vector<double>> U, newU, F; 	// TODO: попробовать вариант с double ** или vector<CVector4> и сравнить по производительности
	C1DField(C1DProblem const& pr);
	~C1DField();
	void setICs(C1DProblem const& pr, CEOS const& eos);
};


#endif




/*const int nGhostCells;
	const int nComp;
	const int iMin, iMax, jMin, jMax, kMin, kMax;
	vector<double> x, y, z;
	const double dx, dy, dz;	
	double t, dt;	
	vector<vector<vector<CVector5>>> U, UNew;	
	vector<vector<vector<CVector5>>> F, G, H, S;
	CField() : nGhostCells(0), nComp(0), iMin(0), iMax(0), jMin(0), jMax(0), kMin(0), kMax(0), 
		       x(0), y(0), z(0), dx(0.), dy(0.), dz(0.),  t(0.), dt(0.), U(0), UNew(0), F(0), G(0), H(0), S(0) {}
	CField(CProblem _problem, CEOSIdeal _eos);
	int setICs(const CProblem& _problem, CEOSIdeal& _eos);
	int setBCs(const CProblem& _problem);
};*/
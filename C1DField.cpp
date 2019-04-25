#include"C1DField.h"

C1DField::C1DField(C1DProblem const& pr): iMin(2), iMax(2+pr.NX), x(vector<double>(2+pr.NX+2)), dx((pr.xMax-pr.xMin)/pr.NX),
	                                      U(vector<vector<double>>()), newU(vector<vector<double>>()), F(vector<vector<double>>()) {
	U.resize(2+pr.NX+2);
	newU.resize(2+pr.NX+2);
	F.resize(2+pr.NX+2);	
}

C1DField::~C1DField() {
	x.clear();
	U.clear();
	newU.clear();
	F.clear();
}

void C1DField::setICs(C1DProblem const& pr, CEOS const& eos) {
	
}
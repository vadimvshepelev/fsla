#include<assert.h>
#include"C1DField.h"

C1DField::C1DField(C1DProblem& pr): imin(2), imax(2+pr.nx), x(vector<double>(2+pr.nx+2+1)), dx((pr.xmax-pr.xmin)/pr.nx), t(0.), dt(0.),
	                                      U(vector<vector<double>>()), newU(vector<vector<double>>()), F(vector<vector<double>>()) {
	
	U.resize(2+pr.nx+2);
	newU.resize(2+pr.nx+2);
	F.resize(2+pr.nx+2);
	for(int i=0; i<imax+2; i++) {
		U[i].resize(3);
		newU[i].resize(3);
		F[i].resize(3);
	}
}

C1DField::~C1DField() {
	x.clear();
	U.clear();
	newU.clear();
	F.clear();
}

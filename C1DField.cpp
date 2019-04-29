#include<assert.h>
#include"C1DField.h"

C1DField::C1DField(C1DProblem& pr): imin(2), imax(2+pr.nx), x(vector<double>(2+pr.nx+2)), dx((pr.xmax-pr.xmin)/pr.nx),
	                                      U(vector<vector<double>>()), newU(vector<vector<double>>()), F(vector<vector<double>>()) {
	U.resize(2+pr.nx+2);
	newU.resize(2+pr.nx+2);
	F.resize(2+pr.nx+2);	
}

C1DField::~C1DField() {
	x.clear();
	U.clear();
	newU.clear();
	F.clear();
}

void C1DField::setics(C1DProblem& pr, CEOSMieGruneisen& eos) {
	
}

void C1DField::setbcs(C1DProblem& pr) {
	assert(pr.bcs[0]=='t');
	switch(pr.bcs[0]) {
	case 't': 
		for(int counter=0; counter<3; counter++) {
			U[imin-1][counter] = U[imin][counter];
			U[imin-2][counter] = U[imin][counter];
		}
		break;
	}
	assert(pr.bcs[1]=='t');
	switch(pr.bcs[0]) {
	case 't': 
		for(int counter=0; counter<3; counter++) {
			U[imax][counter] = U[imax][counter];
			U[imax+1][counter] = U[imax][counter];
		}
		break;
	}
}
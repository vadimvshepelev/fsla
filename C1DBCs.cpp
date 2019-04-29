#include"C1DBCs.h"

void CLeftTransmissiveBoundary::set(C1DField& fld) {
	vector<vector<double>> U = fld.U;
	int imin = fld.imin;
	for(int counter=0; counter<3; counter++) {
		fld.U[imin-1][counter] = fld.U[imin][counter];
		fld.U[imin-2][counter] = fld.U[imin][counter];
	}
}

void CRightTransmissiveBoundary::set(C1DField& fld) {
	vector<vector<double>> U = fld.U;
	int imax = fld.imax;
	for(int counter=0; counter<3; counter++) {
		fld.U[imax][counter] = fld.U[imax][counter];
		fld.U[imax+1][counter] = fld.U[imax][counter];
	}
}



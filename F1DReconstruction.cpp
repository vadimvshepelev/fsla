#include"F1DReconstruction.h"

// ENO2 stencils
double _rm1r[] = {-1./2., 3./2.,  0.   },  
	    _r0r[] = {    0., 1./2.,  1./2.},
	   _rm1l[] = { 1./2., 1./2.,  0.   },
	    _r0l[] = {    0., 3./2., -1./2.};

F1DReconstruction::F1DReconstruction(C1DField& _fld) {
	int i=0;
	double _v[] = {0., 0., 0., 0.};
	vector<double> tempVect = vector<double>(_v, _v+3);
	for(i=0; i<_fld.imax+2; i++) {		
		ULx.push_back(tempVect);
		URx.push_back(tempVect);
	}
}



F1DENO2Reconstruction::F1DENO2Reconstruction(C1DField& _fld) : F1DReconstruction(_fld), 
	                                                           rm1r(vector<double>(_rm1r, _rm1r+3)),
															   r0r(vector<double>(_r0r, _r0r+3)),
															   rm1l(vector<double>(_rm1l, _rm1l+3)),
															   r0l(vector<double>(_r0l, _r0l+3)) {
/*	rm1r = vector<double>(_rm1r, _rm1r+3);		// Для типа vector<T> при инициализация через два параметра первый это указатель на первый элемент, 
	r0r = vector<double>(_r0r, _r0r+3);			// а второй -- указатель на последний элемент. Таким образом очень удобно инициализировать вектор
	rm1l = vector<double>(_rm1l, _rm1l+3);		// через статический массив, который, в свою очередь можно инициализировать списком значений {..., ..., ... etc.}
	r0l = vector<double>(_r0l, _r0l+3);			// Пока не дошли до C11, едем на этом. :) */
}

void F1DENO2Reconstruction::calc(C1DField& fld) {
	int i=0, n=0;
	const int nComp=4;
	int imin=fld.imin, imax=fld.imax;
	vector<vector<double>>& U = fld.U;
	Vector4 diffPlus=Vector4::ZERO, diffMinus=Vector4::ZERO;			 
	int stencilIndex[]={0, 0, 0, 0};
	for(i=imin-1; i<imax+1; i++) {


//		cout << U[i][0];






	    diffPlus = Vector4(U[i+1][0], U[i+1][1], U[i+1][2], U[i+1][3])-Vector4(U[i][0], U[i][1], U[i][2], U[i][3]);
        diffMinus = Vector4(U[i][0], U[i][1], U[i][2], U[i][3])-Vector4(U[i-1][0], U[i-1][1], U[i-1][2], U[i-1][3]);
		for(n=0; n<nComp; n++) {
		    // Choose stencil
			if(fabs(diffMinus[n])<=fabs(diffPlus[n])) stencilIndex[n] = i-1; else stencilIndex[n] = i;
			// Calculate ENO-2 approximations themselves
			if(stencilIndex[n]==i-1) {
				URx[i][n] = rm1r[0]*U[i-1][n] + rm1r[1]*U[i][n] + rm1r[2]*U[i+1][n];
				ULx[i][n] = rm1l[0]*U[i-1][n] + rm1l[1]*U[i][n] + rm1l[2]*U[i+1][n];
			} else {
				URx[i][n] = r0r[0]*U[i-1][n] + r0r[1]*U[i][n] + r0r[2]*U[i+1][n];
				ULx[i][n] = r0l[0]*U[i-1][n] + r0l[1]*U[i][n] + r0l[2]*U[i+1][n];
			}
		}
	}
}
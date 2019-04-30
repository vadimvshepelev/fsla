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

void C1DField::setics(C1DProblem& pr, CEOSMieGruneisen& eos) {
	// TODO
	// Идея сделать "задача" = "начальные условия" + "граничные условия". 
	// Граничные условия реалзиовать тем крутым способом, что уже здесь сделан, но убран из-за пересечений с C1DField, там где они задаются функицей bcs.set(C1DField& fld).
	// Начальные условия -- инкапсулируют постановку начальных условий, тоже через C1DField. У каждого типа задачи своя функиця set, чтобы избежать ветвлений. 
	// (Сейчас не так.) (Пока легко, пока делаем просто под задачу Римана.) (А еще кстати каждая задача может и сама себя печатать по своему своей функцией dump().)
	// Вопрос, зачем тогда нужен класс "задача"? Зачем он нужен в HyperSolver?
	int i = 0;
	double E = 0.;
	for(i=0; i<imax+2+1; i++) {
		x[i] = (pr.xmin - 2.*dx) + (double)i*dx;
	}	
	for(i=imin; i<imax; i++) {
		if(x[i]<pr.x0) {
			U[i][0] = pr.rol;
			U[i][1] = pr.rol*pr.ul;
			E = eos.gete(pr.rol, pr.pl) + .5*pr.ul*pr.ul;
			U[i][2] = pr.rol*E;
		} else {
			U[i][0] = pr.ror;
			U[i][1] = pr.ror*pr.ur;
			E = eos.gete(pr.ror, pr.pr) + .5*pr.ur*pr.ur;
			U[i][2] = pr.ror*E;
		}
	}
	setbcs(pr);
	return;
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
	switch(pr.bcs[1]) {
	case 't': 
		for(int counter=0; counter<3; counter++) {
			U[imax][counter] = U[imax-1][counter];
			U[imax+1][counter] = U[imax-1][counter];
		}
		break;
	}
}
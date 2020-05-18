#include "C1DProblem.h"

void C1DProblem::setics(CEOS& eos, vector<double>& x, vector<vector<double>>& U) {
	// TODO
	// Идея сделать "задача" = "начальные условия" + "граничные условия". 
	// Граничные условия реалзиовать тем крутым способом, что уже здесь сделан, но убран из-за пересечений с C1DField, там где они задаются функицей bcs.set(C1DField& fld).
	// Начальные условия -- инкапсулируют постановку начальных условий, тоже через C1DField. У каждого типа задачи своя функиця set, чтобы избежать ветвлений. 
	// (Сейчас не так.) (Пока легко, пока делаем просто под задачу Римана.) (А еще кстати каждая задача может и сама себя печатать по своему своей функцией dump().)
	// Вопрос, зачем тогда нужен класс "задача"? Зачем он нужен в HyperSolver?
	int i = 0;
	int imin = 2, imax = imin + nx;
	double dx = (xmax - xmin)/nx;
	double E = 0.;
	for(i=0; i<imax+2+1; i++) {
		x[i] = (xmin - 2.*dx) + (double)i*dx;
	}
	for(i=imin; i<imax; i++) {
		if(x[i]<x0) {
			U[i][0] = rol;
			U[i][1] = rol*ul;
			E = eos.gete(rol, pl) + .5*ul*ul;
			U[i][2] = rol*E;
		} else {
			U[i][0] = ror;
			U[i][1] = ror*ur;
			E = eos.gete(ror, pr) + .5*ur*ur;
			U[i][2] = ror*E;
		}
	}
	setbcs(U);
	return;
}

void C1DProblem::setbcs(vector<vector<double>>& U) {
	int imin = 2, imax = imin + nx;
	assert(bcs[0]=='t');
	switch(bcs[0]) {
	case 't': 
		for(int counter=0; counter<3; counter++) {
			U[imin-1][counter] = U[imin][counter];
			U[imin-2][counter] = U[imin][counter];
		}
		break;
	}
	assert(bcs[1]=='t');
	switch(bcs[1]) {
	case 't': 
		for(int counter=0; counter<3; counter++) {
			U[imax][counter] = U[imax-1][counter];
			U[imax+1][counter] = U[imax-1][counter];
		}
		break;
	}
}





// Test for non-ideal (Mie-Gruneisen) EOS of Bolotova-Nigmatullin
C1DProblem prNBtest = C1DProblem("NBtestHLL", 100., 0., 1.e9, 1000., 0., 1.e5, 0., 1., 0., 100.e-6, .7, 1000, .9, "tt"); 
// 5 Toro tests
C1DProblem prToro1Idealtest = C1DProblem("Toro-1-ex", 1.,           .75,      1.,    .125,       0.,      .1,  0., 1., 0.,   .2, .3, 100, .3, "tt");
C1DProblem prToro2Idealtest = C1DProblem("Toro-2", 1.,           -2.,      .4,      1.,       2.,      .4,  0., 1., 0.,  .15, .5, 100, .9, "tt");
C1DProblem prToro3Idealtest = C1DProblem("Toro-3", 1.,            0.,   1000.,      1.,       0.,     .01,  0., 1., 0., .012, .5, 100, .9, "tt");
C1DProblem prToro4Idealtest = C1DProblem("Toro-4", 5.99924,  19.5975, 460.894, 5.99242, -6.19633, 46.0950,  0., 1., 0., .035, .4, 100, .9, "tt");
C1DProblem prToro5Idealtest = C1DProblem("Toro-5-ex", 1.,     -19.59745,   1000.,      1., -19.59745,     .01,  0., 1., 0., .012, .8, 100, .3, "tt");
// 3 tests for V. Denisenko article about Godunov method
C1DProblem prDenisenko1 = C1DProblem("Denisenko-1-8100", 2., 0., 2., 1., 0., 1., 0., 1., 0., .225, .5, 8100, .9, "tt");
C1DProblem prDenisenko2 = C1DProblem("Denisenko-2-8100", 1., -1., 1., 1., 1., 1., 0., 1., 0., .15, .5, 8100, .9, "tt");
C1DProblem prDenisenko3 = C1DProblem("Denisenko-3-8100", 3., 4., 2., 2., 2., 1., 0., 1., 0., .09, .5, 8100, .9, "tt");
// 08.02.2020 First Riemann problem test for laser volume target problem for Al with ideal EOS
C1DProblem prLaserVTAlIdealTest1 = C1DProblem("LaserVTAl1-ideal-2000", 2700., 0., 300.e9, 2., 0., 19.4872e9, -50.e-9, 50e-9, 0., 2.e-12, 0., 2000, .9, "tt");
// 08.02.2020 Second Riemann problem test for laser volume target problem for Al with ideal EOS
C1DProblem prLaserVTAlIdealTest2 = C1DProblem("LaserVTAl2-ideal-2000", 2700., 0., 19.4872e9, 2700., 0., 300.e9, -100.e-9, 0., 0., 2.e-12, -50.e-9, 2000, .9, "tt");
// 15.02.2020 First Riemann problem test for laser VT problem with Mie-Gruneisen EOS
C1DProblem prLaserVTAlMGTest1 = C1DProblem("LaserVTAl1-MG-1nm", 2700., 0., 300.e9, 2., 0., 194872.e9, -50.e-9, 50.e-9, 0., 1.e-13, 0., 100, .9, "tt");
// 16.02.2020 Second Riemann problem test for laserVT problem with Mie-Gruneisen EOS
C1DProblem prLaserVTAlMGTest2 = C1DProblem("LaserVTAl2-MG-1nm", 2700., 0., 19.4872e9, 2700., 0., 300.e9, -100.e-9, 0., 0., 1.e-12, -50.e-9, 100, .9, "tt");
// 16.05.2020 We restarted the test of Riemann problem for Mie-Gruneisen EOS -- just shifted test prLaserVTAlMGTest2
C1DProblem prLaserVTAlMGTestNum1 = C1DProblem("LaserVTAlNum1-MG-1nm-hll", 2700., 0., 19.4872e9, 2700., 0., 300.e9, 0., 100.e-9, 0., 2.e-12, 50.e-9, 100, .9, "tt");
// 17.05.2020 I modified Godunov method in C1DGodunovTypeMethodVacuum class for metal-vacuum border tracking
C1DProblem prIdealVacTest = C1DProblem("IdealVac", 0., 0., 0., 1., 0., 1., 0., 1., 0., .1, .5, 100, .9, "tt");
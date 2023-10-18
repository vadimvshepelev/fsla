#include "C1DProblem.h"

#include "_vector4.h"

/*void C1DProblem::setics(FEOS& eos, vector<double>& x,
						vector<Vector4>& U) {
	// TODO
	// Идея сделать "задача" = "начальные условия" + "граничные условия".
	// Граничные условия реалзиовать тем крутым способом, что уже здесь сделан,
	// но убран из-за пересечений с C1DField, там где они задаются функицей
	// bcs.set(C1DField& fld).
	// Начальные условия -- инкапсулируют постановку начальных условий,
	// тоже через C1DField. У каждого типа задачи своя функиця set,
	// чтобы избежать ветвлений.
	// (Сейчас не так.) (Пока легко, пока делаем просто под задачу Римана.)
	// (А еще кстати каждая задача может и сама себя печатать по-своему
	// своей функцией dump().)
	// Вопрос, зачем тогда нужен класс "задача"?
	// Зачем он нужен в HyperSolver?
	std::size_t i = 0;
	const std::size_t imin = get_order();
	const std::size_t imax = get_order() + nx;
	double dx = (xmax - xmin) / nx;
	double E = 0.;
	for (i = 0; i < x.size(); ++ i) {
		x[i] = (xmin - get_order() * dx) + static_cast<double>(i) * dx;
	}
	for (i = imin; i < imax; ++ i) {
		if (x[i] < x0 && fabs(x[i] - x0) > 1.e-5*dx ) {
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
}*/

void C1DProblem::setics(FEOS& eos, vector<double>& x, vector<Vector4>& U) {
	std::size_t i = 0;
	const std::size_t imin = get_order();
	const std::size_t imax = get_order() + nx;
	double dx = (xmax - xmin) / nx;
	double E = 0.;
	for (i = 0; i < x.size(); ++ i) {
		x[i] = (xmin - get_order() * dx) + static_cast<double>(i) * dx;
	}
	// Проектируем генеральную таблицу соответствия разных методов и классов постановкам задач
	// Здесь используем систему, основанную на флагах
	// Даем флаг primitiveVariablesUsedFlag специально для заполнения примитивными переменными
	bool primitiveVariablesUsedFlag = false;
	if (name == "LiF" ||
		name.substr(name.size()-8, 8) == "lagrange")
		primitiveVariablesUsedFlag = true;

	for (i = imin; i < imax; ++ i) {
		if (x[i] < x0 && fabs(x[i] - x0) > 1.e-5*dx ) {
			if (!primitiveVariablesUsedFlag) {
				U[i][0] = rol;
				U[i][1] = rol * ul;
				E = eos.gete(rol, pl) + .5 * ul * ul;
				U[i][2] = rol * E;
			}
			else {
				U[i][0] = rol;
				U[i][1] = ul;
				U[i][2] = pl;
			}

		} else {
			if (!primitiveVariablesUsedFlag) {
				U[i][0] = ror;
				U[i][1] = ror * ur;
				E = eos.gete(ror, pr) + .5 * ur * ur;
				U[i][2] = ror * E;
			}
			else {
				U[i][0] = ror;
				U[i][1] = ur;
				U[i][2] = pr;
			}
		}
	}
	setbcs(U);
	return;
}


void C1DProblem::setbcs(vector<Vector4>& U) {
	const std::size_t imin = get_order();
	const std::size_t imax = imin + nx;
	assert(bcs[0] == 't' || bcs[0] == 'v');
	switch (bcs[0]) {
	case 't':
		for (std::size_t counter = 1; counter <= imin; ++ counter) {
			for (std::size_t component = 0; component < 3; ++ component) {
				U[imin - counter][component] = U[imin][component];
			}
		}
		break;
	case 'v':  // Warning! This means density rho = 0!
		for (std::size_t counter = 1; counter <= imin; ++counter) {
			for (std::size_t component = 0; component < 3; ++component) {
				U[imin - counter][component] = 0.;
			}
		}
		break;
	}
	assert(bcs[1] == 't' || bcs[0] == 'v');
	switch (bcs[1]) {
	case 't':
		for (std::size_t counter = 0; counter < imin; ++ counter) {
			for (std::size_t component = 0; component < 3; ++ component) {
				U[imax + counter][component] = U[imax - 1][component];
			}
		}
		break;
	case 'v':
		for (std::size_t counter = 0; counter < imin; ++counter) {
			for (std::size_t component = 0; component < 3; ++component) {
				U[imax + counter][component] = 0.;
			}
		}
		break;
	}
}


void C1DLaserProblem::setics(
		FEOS& eos,
		vector<double>& x,
		vector<Vector4>&/*vector<vector<double>>&*/ U) {
	std::size_t i = 0;
	const std::size_t imin = get_order();
	const std::size_t imax = imin + nx;
	double dx = (xmax - xmin) / nx;
	double E = 0.;
	for (i = 0; i < x.size(); ++ i) {
		x[i] = (xmin - get_order() * dx) + static_cast<double>(i) * dx;
	}
	for (i = imin; i < imax; ++ i) {
		if (x[i] < x0 && fabs(x[i] - x0) > 1.e-5*dx ) {
			U[i][0] = rol;
			U[i][1] = rol*ul;
			E = eos.gete(rol, pl) + .5*ul*ul;
			U[i][2] = rol*E;
		} else if (x[i] < x1 && fabs(x[i] - x1) > 1.e-5*dx ) {
			U[i][0] = rom;
			U[i][1] = rom*um;
			E = eos.gete(rom, pm) + .5*um*um;
			U[i][2] = rom*E;
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



// Test for non-ideal (Mie-Gruneisen) EOS of Bolotova-Nigmatullin
C1DProblem prNBtest = C1DProblem("NBtestHLL", 100., 0., 1.e9, 1000., 0., 1.e5, 0., 1., 0., 100.e-6, .7, 1000, .9, "tt");
// 5 Toro tests
C1DProblem prToro1Idealtest = C1DProblem("Toro-1-lagrange", 1.,           .75,      1.,    .125,       0.,      .1,  0., 1., 0.,   .2, .3, 100, .1, "tt", 1);
C1DProblem prToro2Idealtest = C1DProblem("Toro-2", 1.,           -2.,      .4,      1.,       2.,      .4,  0., 1., 0.,  .15, .5, 100, .9, "tt");
C1DProblem prToro3Idealtest = C1DProblem("Toro-3", 1.,            0.,   1000.,      1.,       0.,     .01,  0., 1., 0., .012, .5, 100, .9, "tt");
C1DProblem prToro4Idealtest = C1DProblem("Toro-4", 5.99924,  19.5975, 460.894, 5.99242, -6.19633, 46.0950,  0., 1., 0., .035, .4, 100, .9, "tt");
C1DProblem prToro5Idealtest = C1DProblem("Toro-5-ex", 1.,     -19.59745,   1000.,      1., -19.59745,     .01,  0., 1., 0., .012, .8, 100, .3, "tt");
// 3 tests for V. Denisenko article about Godunov method
C1DProblem prDenisenko1 = C1DProblem("Denisenko-1-8100", 2., 0., 2., 1., 0., 1., 0., 1., 0., .225, .5, 8100, .9, "tt");
C1DProblem prDenisenko2 = C1DProblem("Denisenko-2-8100", 1., -1., 1., 1., 1., 1., 0., 1., 0., .15, .5, 8100, .9, "tt");
C1DProblem prDenisenko3 = C1DProblem("Denisenko-3-8100", 3., 4., 2., 2., 2., 1., 0., 1., 0., .09, .5, 8100, .9, "tt");
// 08.02.2020 First Riemann problem test for laser volume target problem for Al with ideal EOS
C1DProblem prLaserVTAlIdealTest1 = C1DProblem("LaserVTAl1-ideal", 2700., 0., 50.e9, 2., 0., 19.4872e9, -50.e-9, 50e-9, 0., 2.e-12, 0., 500, .9, "tt");
// 08.02.2020 Second Riemann problem test for laser volume target problem for Al with ideal EOS
C1DProblem prLaserVTAlIdealTest2 = C1DProblem("LaserVTAl2-ideal-100-corr20201226", 2700., 0., 19.4872e9, 2700., 0., 300.e9, -100.e-9, 0., 0., 2.e-12, -50.e-9, 100, .9, "tt");
// 15.02.2020 First Riemann problem test for laser VT problem with Mie-Gruneisen EOS
C1DProblem prLaserVTAlMGTest1 = C1DProblem("LaserVTAl1-MG-1nm",
										   2700., 0., 300.e9,
										   2., 0., 194872.e9,
										   -50.e-9, 50.e-9, 0., 1.e-13, 0., 100, .1, "tt");
// 16.02.2020 Second Riemann problem test for laserVT problem with Mie-Gruneisen EOS
C1DProblem prLaserVTAlMGTest2 = C1DProblem("LaserVTAl2-MG-1nm", 2700., 0., 19.4872e9, 2700., 0., 300.e9, -100.e-9, 0., 0., 1.e-12, -50.e-9, 100, .9, "tt");
// 16.05.2020 We restarted the test of Riemann problem for Mie-Gruneisen EOS -- just shifted test prLaserVTAlMGTest2
C1DProblem prLaserVTAlMGTestNum1 = C1DProblem("vtAlNum1-MG-1nm-hll", 2700., 0., 19.4872e9, 2700., 0., 300.e9, 0., 100.e-9, 0., 2.e-12, 50.e-9, 100, .9, "tt");
// 17.05.2020 I modified Godunov method in C1DGodunovTypeMethodVacuum class for metal-vacuum border tracking
C1DProblem prIdealVacTest = C1DProblem("IdealVac", 0., 0., 0., 1., 0., 1., 0., 1., 0., .1, .5, 100, .9, "tt");
// 30.12.2020 modified tests with new precise MG EOS
// C1DProblem prVTAlMGTest1 = C1DProblem("vtAl1MG-1nm-hll", 2700., 0., 0., 2700., 0., 300.e9, -100.e-9, 0., 0., 20.e-12, -50.e-9, 1000, .5, "tt", 3);
C1DProblem prVTAlMGTest1 = C1DProblem("vtAl1MG-1nm-lf-weno5", 2700., 0., 0., 2700., 0., 300.e9, -100.e-9, 0., 0., 2.e-12, -50.e-9, 100, .9, "tt", 3);
C1DProblem prVTAlMGTest2 = C1DProblem("vtAl2MG-1nm-hllc", 2700., 0., 300.e9, 2., 0., 0., -50.e-9, 50.e-9, 0., 1.5e-13, 0., 100, .5, "tt");
// 29.03.2021 modified test#2 with cold metal and vacuum -- result is non-moving (stable) contact boundary because of equal pressures
C1DProblem prVTAlMGTest2_2 = C1DProblem("vtAl2MG-1nm-roegen", 2413., 0., 20.e9, 2., 0., 20.e9, -50.e-9, 50.e-9, 0., 1.e-12, 0., 100, .9, "tt");

// 30.05.2022 problem for synchronizing with Baer-Nunziato model and Petr Chuprov
C1DLaserProblem prHoles = C1DLaserProblem("holes",
										  2230., 0., 100.e9, 19300., 0., 1000.e9, 2., 0., 100.e9,
										  -1250.e-9, 1000.e-9, 0., 10.e-12, -50.e-9, 0.,
										  2250, .9, "tt");

// 30.06.2021 Vova shockwave
C1DProblem prDenisenko = C1DProblem("denisenko-roe-eno2", 1.2714, .2928, 1.4017, 1., 0., 1., 0., 1., 0., .4, .2, 100, .3, "tt");

C1DProblem prFedorAl = C1DProblem("FedorAl", 1., 1., 1., .1, 0., .125, 0., 1., 0., .2, .3, 100, .9, "tt");
C1DProblem prVTAlBel = C1DLaserProblem("vtAlBel",
										2700., 0., 0.,
										2700., 0., 1000.,
										2., 0., 0.,
										0., 1250.,
										0., .01,
										1000., 1050.,
										2000, 0.4,
										"tt", 3);

// 16.04.2023 LiF problem
C1DProblem prLiF = C1DProblem("LiF", 2650., 0., 25.e9, 2650., 0., 0., 0., 50.e-6, 0., 10.e-9, 5.e-6, 500, .05, "tt", 2);
// C1DProblem prLiF = C1DProblem("LiF", 2650., 0., 25.e9, 2650., 0., 25.e9, 0., 50.e-6, 0., .5e-9, 5.e-6, 500, .1, "vt", 1);
C1DProblem prTestLagrange1D = C1DProblem("LiF", 1., 0., 1., 1, 0., 1., 0., 1., 0., .2, .5, 100, .5, "vt", 1);

// 18.10.2023 LiF slice (LiF 2D problem simulation doesn't work, so we make 1d simulation on the symmetry axis
// Eulerian version
C1DProblem prLiFSliceEuler = C1DProblem("LiF", 2650., 0., 10.e9, 2., 0., 0., -10.e-6, 2.e-6, 0., 1.e-9, 0., 100, .1, "tt", 2);
//Lagrangian version
// C1DProblem prTestLagrange1D = C1DProblem("LiF", 1., 0., 1., 1, 0., 1., 0., 1., 0., .2, .5, 100, .5, "vt", 1);
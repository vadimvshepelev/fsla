#include <cassert>
#include <ctime>
#include <iomanip>

#include "F1DSimulation.h"



void F1DSimulation::run() {
	pr.setics(eos, fld.x, fld.U);
	int counter = 0;

	clock_t tStartGlobal = 0;
	clock_t tEndGlobal = 0;
	clock_t tStart = 0;
	clock_t tEnd = 0;

	double cfl = 0.;
	double tCalc = 0.;

	outp.manageFileOutput(pr, fld, eos);
	cout << "Starting simulation..." << endl;

	tStartGlobal = clock();
	while (fld.t < pr.tmax) {
		cfl = pr.cfl;
		fld.dt = mtd.calcdt(pr, eos, fld);
		if (counter < 5) {
			cfl = cfl * .2;
			fld.dt *= .2;
		}
		if (fld.t+fld.dt > pr.tmax) fld.dt = pr.tmax-fld.t;
		tStart = clock();
		mtd.calc(pr, eos, fld);
		// solve.solve();
		tEnd = clock();
		tCalc = static_cast<double>(tEnd - tStart) / CLOCKS_PER_SEC;
		fld.t += fld.dt;
		outp.manageScreenOutput(pr, counter, fld.t, fld.dt, cfl, tCalc);
		outp.manageFileOutput(pr, fld, eos);
		++ counter;
	}
	tEndGlobal = clock();

	tCalc = static_cast<double>(
				tEndGlobal - tStartGlobal) / CLOCKS_PER_SEC;
	std::cout << std::setprecision(3) << "total time=" << tCalc << "s" << "\n";
	cout << "...done!" << "\n";
	return;
}


void F1DSimulationLagrange::run() {
	pr.setics(eos, fld.x, fld.W);
	int counter = 0;

	clock_t tStartGlobal = 0;
	clock_t tEndGlobal = 0;
	clock_t tStart = 0;
	clock_t tEnd = 0;

	double cfl = 0.;
	double tCalc = 0.;

	outp.manageFileOutput(pr, fld, eos);
	cout << "Starting simulation..." << endl;

	tStartGlobal = clock();
	while (fld.t < pr.tmax) {
		cfl = pr.cfl;
		fld.dt = mtd.calcdt(pr, eos, fld);
		if (counter < 5) {
			cfl = cfl * .2;
			fld.dt *= .2;
		}
		if (fld.t + fld.dt > pr.tmax) fld.dt = pr.tmax - fld.t;
		tStart = clock();




		if (counter == 331) {

			double pp = 0.;
		}





		int itNum = mtd.calc(pr, eos, fld);
		tEnd = clock();
		tCalc = static_cast<double>(tEnd - tStart) / CLOCKS_PER_SEC;
		fld.t += fld.dt;
		outp.manageScreenOutput(pr, counter, fld.t, fld.dt, cfl, tCalc, itNum);
		outp.manageFileOutput(pr, fld, eos);
		++counter;


		/*
		// Определяем, чему равны давление, скорость и плотность в точке x = -.5 нм
		double x0 = -.5e-6, rho0=0., u0=0., p0=0., alpha = 0.;
		int i0 = 0;
		for (auto i = 0; i < fld.imax-1; i++) {
			double xL = fld.x[i] + .5 * (fld.x[i + 1] - fld.x[i]), xR = fld.x[i+1] + .5*(fld.x[i+2]-fld.x[i+1]);
			if ( xL < x0 && x0 < xR) {
				i0 = i;
				alpha = x0 - xL / (xR - xL);
				rho0 = alpha * fld.W[i0][0] + (1. - alpha) * fld.W[i0 + 1][0];
				u0 = alpha * fld.W[i0][1] + (1. - alpha) * fld.W[i0 + 1][1];
				p0 = alpha * fld.W[i0][1] + (1. - alpha) * fld.W[i0 + 1][2];
				break;
			}
			rho0 = .01;
			u0 = 0.;
			p0 = 0.;
		}

		if (p0 != 0.)
			cout << endl << fld.t << " " << p0 << endl;
		*/


	}
	tEndGlobal = clock();

	tCalc = static_cast<double>(
		tEndGlobal - tStartGlobal) / CLOCKS_PER_SEC;
	std::cout << std::setprecision(3) << "total time=" << tCalc << "s" << "\n";
	cout << "...done!" << "\n";
	return;
}

#include <cassert>
#include <ctime>
#include <iomanip>

#include "F1DSimulation.h"



void F1DSimulation::run() {
	pr.setics(eos, fld.x, fld.U);
	int counter = 0;

	clock_t tStart = 0;
	clock_t tEnd = 0;
	clock_t tStartGlobal = 0;
	clock_t tEndGlobal = 0;

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
		// mtd.calc(pr, eos, fld);
		solve.solve();
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

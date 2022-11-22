#include<time.h>
#include<assert.h>
#include "F1DSimulation.h"



void F1DSimulation::run() {
	clock_t tStart = 0, tEnd = 0;
	tStart = clock(); 

	pr.setics(eos, fld.x, fld.U);
	int counter=0;	
	double cfl=0., tCalc=0.;		
	outp.manageFileOutput(pr, fld, eos);
	cout << "Starting simulation..." << endl;
	
	while(fld.t < pr.tmax) {
		cfl = pr.cfl;
		fld.dt = mtd.calcdt(pr, eos, fld);
		if(counter < 5) {
			cfl = cfl*.2;
			fld.dt *= .2;
		}
		if(fld.t+fld.dt > pr.tmax) fld.dt = pr.tmax-fld.t; 						
		mtd.calc(pr, eos, fld); 
		fld.t += fld.dt;				
		// outp.manageScreenOutput(pr, counter, fld.t, fld.dt, cfl, 0./*tCalc*/);
		outp.manageFileOutput(pr, fld, eos);
		counter++;
	}	
	cout << "...done!" << endl;

	tEnd = clock(); 
	tCalc = (double)(tEnd - tStart) / CLOCKS_PER_SEC;		
	cout << "Totally " << tCalc << " secs taken" << endl;
	return;
}

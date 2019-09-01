#ifndef _COUTPUT_H_
#define _COUTPUT_H_

#include "c1dproblem.h"
#include "c1dfield.h"

class COutput {
	string subDir;
	vector<double> dtt; // Vector of time moments to write .dat files to disk
	string tUnit;
	double tMul;
	int nDump, tPrecision, dtPrecision;	
	string getProgressBar(C1DProblem& _pr, double t);
	int print(void);
	//int dump(C1DProblem& pr, C1DField& fld, CEOSMieGruneisen& eos, string fName);	
	int dump(C1DProblem& pr, C1DField& fld, CEOS& eos, string fName);	
public:
	COutput(C1DProblem& _pr, string _subdir, vector<double> _dtt);
	//int manageFileOutput(C1DProblem& pr, C1DField& fld, CEOSMieGruneisen& eos);
	int manageFileOutput(C1DProblem& pr, C1DField& fld, CEOS& eos);
	int manageScreenOutput(C1DProblem& pr, int iteration, double t, double dt, double cfl, double tCalc);
};



#endif
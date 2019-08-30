#include <assert.h>

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>

#include"coutput.h"
#include"solver.h"
#include<boost\\filesystem.hpp>

#define BOOST_LIB_DIAGNOSTIC
#define BOOST_SYSTEM_NO_DEPRECATED

using namespace std;

 COutput::COutput(C1DProblem& pr, string _subdir, vector<double> _dtt) : 
	              subDir(_subdir), dtt(_dtt), tUnit(""), tMul(1.), nDump(0), tPrecision(4), dtPrecision(2) {
    // Разбираемся с каталогом и создаем новый, если его нет	
	cout << "Creating subdirectory ";
	boost::filesystem::path outputPath(subDir);
		if(!exists(outputPath)) {
			cout << "'" << subDir << "'" << "...";
			if(!boost::filesystem::create_directory(outputPath)) {
				cerr << endl << "Error: cannot create 'output' subdirectory.";
				exit(1);
			}
			cout << "done!" << endl;
		}		
		subDir += "\\" + pr.name;
		boost::filesystem::path p(subDir);
	if(exists(p)) {
		cout << "'" << p.string() << "'" << "...already exists!" << endl;
	} else if(!boost::filesystem::create_directory(p)) {
		cerr << "Error: cannot make a subdirectory.";
		exit(1);
	} else 
		cout << p << "...subdirectory successfully created!" << endl;	
 }

 string COutput::getProgressBar(C1DProblem& _pr, double t) {
	 int i = 0;
	 ostringstream oss;	 
	 oss << "[";
	 double progressRate = (t - _pr.tmin)/(_pr.tmax - _pr.tmin);        
	 for(i = 0; i <= 10; i++) 
		 if ((int)progressRate >= i) 
			 oss << "."; 
		 else 
			 oss << " ";
	 oss << "] " << (int)progressRate*100 << "% ";
	 return oss.str();
 }

int COutput::manageScreenOutput(C1DProblem& _prm, int iteration, double t, double dt, double CFL, double tCalc) {
	ostringstream oss;
	string buf;
	oss << /*"\r" << getProgressBar(_problem, t) <<*/ "iter=" << iteration << 
		   setprecision(tPrecision)  << " t="  << t        << 
		   setprecision(dtPrecision) << " dt=" << dt       << " CFL=" << CFL << " time=" << tCalc << "s" << endl;
	cout << oss.str();
	return 1;
}

int COutput::manageFileOutput(C1DProblem& pr, C1DField& fld, CEOSMieGruneisen& eos) {	
	assert(!dtt.empty());		
	if (fld.t>=dtt[0]) {
		ostringstream oss1;
		oss1 << subDir << "\\" << pr.name << "-" << nDump++ << ".dat"; string fName1 = oss1.str();
		cout << "Writing to file '" << fName1 << "'...";		
		dump(pr, fld, eos, fName1);
		cout << "done!" << endl;		
		dtt.erase(dtt.begin());		
	}
	return 1;
}

int COutput::manageFileOutput(C1DProblem& pr, C1DField& fld, CEOSIdeal& eos) {	
	assert(!dtt.empty());		
	if (fld.t>=dtt[0]) {
		ostringstream oss1;
		oss1 << subDir << "\\" << pr.name << "-" << nDump++ << ".dat"; string fName1 = oss1.str();
		cout << "Writing to file '" << fName1 << "'...";		
		dump(pr, fld, eos, fName1);
		cout << "done!" << endl;		
		dtt.erase(dtt.begin());		
	}
	return 1;
}

int COutput::dump(C1DProblem& prb, C1DField& fld, CEOSMieGruneisen& eos, string fName) {
	int i = 0, imin = fld.imin, imax = fld.imax;
	vector<vector<double>> U = fld.U;
	vector<double> x = fld.x;
	double t = fld.t, dx = fld.dx; 
	double rol = prb.rol, ul = prb.ul, pl = prb.pl, ror = prb.ror, ur = prb.ur, pr = prb.pr, x0 = prb.x0;
	ofstream ofs(fName);
	CVectorPrimitive res = CVectorPrimitive();	
	double e_ex = 0., q=0.;
	if(!ofs) {
		cout << "COutput::dump1D() reports error: cannot open output file." << endl;		
		exit(1);
	}
	ofs << "TITLE=\"Riemann Problem 1D slice t=" << t << "\"" << endl;
	ofs << "VARIABLES=\"x\",\"ro\",\"u\",\"p\",\"e\",\"ro_ex\",\"u_ex\",\"p_ex\",\"e_ex\"" << endl;	
	ofs << "ZONE T=\"Numerical\", I=" << imax - imin << ", F=POINT" << endl;
	double mul_x=1., mul_u=1., mul_p=1., mul_e=1.;
	double _ro = 0., _u = 0., _v = 0., _w = 0., _e = 0., _p = 0.;	
	for(i = imin; i < imax; i++) {						
		_ro = U[i][0];
		_u  = U[i][1]/_ro; 
		_e  = U[i][2]/_ro - .5*_u*_u;
		_p  = eos.getp(_ro, _e); 		
		//res = calcRPAnalyticalSolution (eos, rol, ul, pl, ror, ur, pr, x-x0, t);
		e_ex = eos.gete(res.ro, res.p) ? res.ro!=0 : 0.;
		ofs << (fld.x[imin+i]+.5*fld.dx)*mul_x << " " << _ro    << " " << _u*mul_u << " " << _p*mul_p << " " << _e*mul_e << " " << 
			   res.ro << " " << res.v*mul_u << " " << res.p*mul_p << " " << e_ex*mul_e << endl;				
	}	
	ofs.close();	
	return 1;
}

int COutput::dump(C1DProblem& prb, C1DField& fld, CEOSIdeal& eos, string fName) {
	int i = 0, imin = fld.imin, imax = fld.imax;
	vector<vector<double>> U = fld.U;
	vector<double> x = fld.x;
	double t = fld.t, dx = fld.dx; 
	double rol = prb.rol, ul = prb.ul, pl = prb.pl, ror = prb.ror, ur = prb.ur, pr = prb.pr, x0 = prb.x0;
	ofstream ofs(fName);
	CVectorPrimitive res = CVectorPrimitive();	
	double e_ex = 0., q=0.;
	if(!ofs) {
		cout << "COutput::dump1D() reports error: cannot open output file." << endl;		
		exit(1);
	}
	ofs << "TITLE=\"Riemann Problem 1D slice t=" << t << "\"" << endl;
	ofs << "VARIABLES=\"x\",\"ro\",\"u\",\"p\",\"e\",\"ro_ex\",\"u_ex\",\"p_ex\",\"e_ex\"" << endl;	
	ofs << "ZONE T=\"Numerical\", I=" << imax - imin << ", F=POINT" << endl;
	double mul_x=1., mul_u=1., mul_p=1., mul_e=1.;
	double _ro = 0., _u = 0., _v = 0., _w = 0., _e = 0., _p = 0.;	
	for(i = imin; i < imax; i++) {						
		_ro = U[i][0];
		_u  = U[i][1]/_ro; 
		_e  = U[i][2]/_ro - .5*_u*_u;
		_p  = eos.getp(_ro, _e); 		
		//res = calcRPAnalyticalSolution (eos, rol, ul, pl, ror, ur, pr, x-x0, t);
		e_ex = eos.gete(res.ro, res.p) ? res.ro!=0 : 0.;
		ofs << (fld.x[imin+i]+.5*fld.dx)*mul_x << " " << _ro    << " " << _u*mul_u << " " << _p*mul_p << " " << _e*mul_e << " " << 
			   res.ro << " " << res.v*mul_u << " " << res.p*mul_p << " " << e_ex*mul_e << endl;				
	}	
	ofs.close();	
	return 1;
}
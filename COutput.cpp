#include <assert.h>

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>

#include"COutput.h"
#include <boost/filesystem.hpp>

#define BOOST_LIB_DIAGNOSTIC
#define BOOST_SYSTEM_NO_DEPRECATED

using namespace std;

 COutput::COutput(C1DProblem& pr, string _subdir, vector<double> _dtt) : 
	              subDir(_subdir), dtt(_dtt), tUnit("ps"), tMul(1.e9), nDump(0), tPrecision(4), dtPrecision(2) {
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
		subDir += "/" + pr.name;
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

/*int COutput::manageFileOutput(C1DProblem& pr, C1DField& fld, CEOSMieGruneisen& eos) {	
	assert(!dtt.empty());		
	if (fld.t>=dtt[0]) {
		ostringstream oss1;
		oss1 << subDir << "/" << pr.name << "-" << nDump++ << ".dat"; string fName1 = oss1.str();
		cout << "Writing to file '" << fName1 << "'...";		
		dump(pr, fld, eos, fName1);
		cout << "done!" << endl;		
		dtt.erase(dtt.begin());		
	}
	return 1;
}
*/
int COutput::manageFileOutput(C1DProblem& pr, C1DField& fld, FEOS& eos) {	
	assert(!dtt.empty());	
	if (fld.t>=dtt[0]) {
		ostringstream oss1;
		oss1 << subDir << "/" << pr.name << "-" << nDump++ << ".dat"; string fName1 = oss1.str();
		cout << "Writing to file '" << fName1 << "'...";		
		dump(pr, fld, eos, fName1);
		cout << "done!" << endl;		
		dtt.erase(dtt.begin());		
	}



	
	
	
/*	ostringstream oss2;
	oss2 << subDir << "/" << pr.name << "-" << 100 << ".dat"; string fName1 = oss2.str();
	cout << "Writing to file '" << fName1 << "'...";		
	dump(pr, fld, eos, fName1);
	cout << "done!" << endl;		*/
	


	return 1;
}

int COutput::dump(C1DProblem& prb, C1DField& fld, FEOS& eos, string fName) {
	int i = 0, imin = fld.imin, imax = fld.imax;
	auto U = fld.U;
	vector<double> x = fld.x;
	double t = fld.t, dx = fld.dx; 
	double rol = prb.rol, ul = prb.ul, pl = prb.pl, ror = prb.ror, ur = prb.ur, pr = prb.pr, x0 = prb.x0;
	ofstream ofs(fName);
	CVectorPrimitive res = CVectorPrimitive();	
	double ro_ex = 0., p_ex = 0., e_ex = 0., q=0.;
	if(!ofs) {
		cout << "COutput::dump1D() reports error: cannot open output file." << endl;		
		exit(1);
	}
	if(eos.gettype() == "ideal" && prb.name != "holes") {
		ofs << "TITLE=\"Riemann Problem 1D slice t=" << t << "\"" << endl;
		ofs << "VARIABLES=\"x\",\"rho\",\"u\",\"p\",\"e\",\"rho_ex\",\"u_ex\",\"p_ex\",\"e_ex\"" << endl;	
		ofs << "ZONE T=\"Numerical\", I=" << imax - imin << ", F=POINT" << endl;
		double mul_x=1., mul_u=1., mul_p=1., mul_e=1.;
		double _ro = 0., _u = 0., _v = 0., _w = 0., _e = 0., _p = 0.;	
		for(i = imin; i < imax; i++) {						
			_ro = U[i][0];
			if(_ro!=0) {
				_u  = U[i][1]/_ro; 
				_e  = U[i][2]/_ro - .5*_u*_u;
			} else {
				_u = 0.;
				_e = 0.;
			}
			_p  = eos.getp(_ro, _e); 		
			res = calcRPAnalyticalSolution(eos, rol, ul, pl, ror, ur, pr, x[i]+.5*dx-x0, t);
			ro_ex = res.ro;
			p_ex = res.p;
			e_ex = eos.gete(ro_ex, p_ex);
			ofs << (fld.x[i]+.5*dx)*mul_x << " " << _ro    << " " << _u*mul_u << " " << _p*mul_p << " " << _e*mul_e << " " << 
				   res.ro << " " << res.v*mul_u << " " << res.p*mul_p << " " << e_ex*mul_e << endl;				
		}	
	} else {
		ofs << "TITLE=\"Laser problem, t=" << t << "\"" << endl;
		ofs << "VARIABLES=\"x[nm]\",\"ro[kg/m3]\",\"u[m/s]\",\"p[GPa]\",\"e[MJ/kg]\"" << endl;	
		ofs << "ZONE T=\"Numerical\", I=" << imax - imin << ", F=POINT" << endl;
		double mul_x=1.e9, mul_u=1., mul_p=1.e-9, mul_e=1.e-6;
		double _ro = 0., _u = 0., _v = 0., _w = 0., _e = 0., _p = 0.;	
		for(i = imin; i < imax; i++) {						
			_ro = U[i][0];
			if(_ro!=0) {
				_u  = U[i][1]/_ro; 
				_e  = U[i][2]/_ro - .5*_u*_u;
			} else {
				_u = 0.;
				_e = 0.;
			}
			_p  = eos.getp(_ro, _e); 		
			ofs << (fld.x[i]+.5*dx)*mul_x << " " << _ro    << " " << _u*mul_u << " " << _p*mul_p << " " << _e*mul_e << " " << endl;				
	    }
	}
	ofs.close();	
	return 1;
}
/*
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
}*/

CVectorPrimitive COutput::calcRPAnalyticalSolution(FEOS& eos, double roL, double vL, double pL, double roR, double vR, double pR, double x, double t){
	RPSolutionPrimitive res = solveRP(eos, roL, vL, pL, roR, vR, pR);
	// V = (ro, v, p)T
	CVectorPrimitive V;
	double xi = 0.;
    if (t != 0.) xi = x/t; else xi = 0.;
	double cL = eos.getc(roL, pL), cR = eos.getc(roR, pR);
	const double gamma = eos.getc(1., 1.)*eos.getc(1., 1.);
	double xiFront=0., xiHead=0., xiTail=0., xiHeadL=0., xiTailL=0., xiHeadR=0., xiTailR=0.;
	// Если вакуум
	if(res.type == VacRW) { 
		xiHead = vR + cR;
		xiTail = vR - 2.*cR/(gamma-1.);
		if(xi<=xiTail) {
			V.ro = 0.;
			V.v  = vR - 2.*cR/(gamma-1.);
			V.p  = 0.;
		} else if (xi<xiHead) {
			V.ro = roR*pow(2./(gamma+1.) - (gamma-1.)/(gamma+1.)/cR*(vR-xi), 2./(gamma-1.));
			V.v  = 2./(gamma+1)*(-cR + (gamma-1.)/2.*vR + xi);
			V.p  = pR*pow(2./(gamma+1.) - (gamma-1.)/(gamma+1.)/cR*(vR-xi), 2.*gamma/(gamma-1.)); 
		} else {
			V.ro = roR;
			V.v  = vR;
			V.p  = pR;
		}
		return V;
	}
	if(res.type == RWVac) {
		xiHead = vL - cL;
		xiTail = vL + 2.*cL/(gamma-1.);
		if(xi>=xiTail) {
			V.ro = 0.;
			V.v  = 0.;
			V.p  = 0.;
		} else if (xi>xiHead) {
			V.ro = roL*pow(2./(gamma+1.)+(gamma-1.)/(gamma+1.)/cL*(vL-xi), 2./(gamma-1.));
			V.v  = 2./(gamma+1)*(cL + (gamma-1.)/2.*vL + xi);
			V.p  = pL*pow(2./(gamma+1.)+(gamma-1.)/(gamma+1.)/cL*(vL-xi), 2.*gamma/(gamma-1.));
		} else {
			V.ro = roL;
			V.v = vL;
			V.p = pL;
		}
		return V;
	}
	if(res.type == RWVacRW) {
		xiHeadL = vL - cL;
		xiTailL = vL + 2.*cL/(gamma-1.);
		xiHeadR = vR + cR;
		xiTailR = vR - 2.*cR/(gamma-1.);
		if(xi<=xiHeadL) {
			V.ro = roL;
			V.v  = vL;
			V.p  = pL;
		} else if (xi<xiTailL) {
			V.ro = roL*pow(2./(gamma+1.)+(gamma-1.)/(gamma+1.)/cL*(vL-xi), 2./(gamma-1.));
			V.v  = 2./(gamma+1)*(cL + (gamma-1.)/2.*vL + xi);
			V.p  = pL*pow(2./(gamma+1.)+(gamma-1.)/(gamma+1.)/cL*(vL-xi), 2.*gamma/(gamma-1.));
		} else if (xi<=xiTailR) {
			V.ro = 0.;
			V.v  = 0.;
			V.p  = 0.;
		} else if (xi<xiHeadR) {
			V.ro = roR*pow(2./(gamma+1.) - (gamma-1.)/(gamma+1.)/cR*(vR-xi), 2./(gamma-1.));
			V.v  = 2./(gamma+1)*(-cR + (gamma-1.)/2.*vR + xi);
			V.p  = pR*pow(2./(gamma+1.) - (gamma-1.)/(gamma+1.)/cR*(vR-xi), 2.*gamma/(gamma-1.)); 
		} else {
			V.ro = roR;
			V.v  = vR;
			V.p  = pR;
		}
		return V;
	}
	double cLLocal = sqrt(gamma*res.p/res.roL), cRLocal = sqrt(gamma*res.p/res.roR);
	// Если не вакуум. Пусть точка слева от контактного разрыва (xiContact = res.v)
	if(xi<res.v) {
		if(res.type == SWSW || res.type == SWRW) { 
			xiFront = vL - cL*sqrt((gamma+1.)/2./gamma*res.p/pL + (gamma-1.)/2./gamma);
			if(xi<xiFront) {
				V.ro = roL;
				V.v  = vL;
				V.p  = pL;
			} else {
				V.ro = res.roL;
				V.v = res.v;
				V.p = res.p;
			}
		} else if (res.type == RWSW || res.type == RWRW) {
			xiHead = vL-cL;
			xiTail = res.v-cLLocal;
			if(xi<=xiHead) {
				V.ro = roL;
				V.v  = vL;
				V.p  = pL;
			} else if(xi>=xiTail) {
				V.ro = res.roL;
				V.v  = res.v;
				V.p  = res.p;
			} else {
				V.ro = roL*pow(2./(gamma+1.)+(gamma-1.)/(gamma+1.)/cL*(vL-xi), 2./(gamma-1.));
				V.v  = 2./(gamma+1)*(cL + (gamma-1.)/2.*vL + xi);
				V.p  = pL*pow(2./(gamma+1.)+(gamma-1.)/(gamma+1.)/cL*(vL-xi), 2.*gamma/(gamma-1.));
			}
		} 
	//Пусть точка справа от контактного разрыва (xiContact = res.v)
	} else {
		if(res.type == RWSW || res.type == SWSW) {
			xiFront = vR + cR*sqrt((gamma+1.)/2./gamma*res.p/pR + (gamma-1.)/2./gamma);
			if(xi>xiFront) {
				V.ro = roR;
				V.v  = vR;
				V.p  = pR;
			} else {
				V.ro = res.roR;
				V.v  = res.v;
				V.p  = res.p;
			}
		} else if(res.type == RWRW || res.type == SWRW) {
			xiHead = vR + cR;
			xiTail = res.v + cRLocal;
			if(xi >= xiHead) {
				V.ro = roR;
				V.v  = vR;
				V.p  = pR;
			} else if (xi <= xiTail) {
				V.ro = res.roR;
				V.v  = res.v;
				V.p  = res.p;
			} else {
				V.ro = roR*pow(2./(gamma+1.) - (gamma-1.)/(gamma+1.)/cR*(vR-xi), 2./(gamma-1.));
				V.v  = 2./(gamma+1)*(-cR + (gamma-1.)/2.*vR + xi);
				V.p  = pR*pow(2./(gamma+1.) - (gamma-1.)/(gamma+1.)/cR*(vR-xi), 2.*gamma/(gamma-1.));
			}
		}
	}
	return V;
}

RPSolutionPrimitive COutput::solveRP(FEOS& eos, double roL, double vL, double pL, double roR, double vR, double pR) {
	// Решаем нелинейное уравнение относительно давления методом касательных Ньютона
	RPSolutionPrimitive res; res.roL = 0.; res.roR=0.; res.v = 0.; res.p = 0.;
	double p = 0., pPrev = 0.;
	double TOL = 1.e-6;
    double cL = eos.getc(roL, pL), cR = eos.getc(roR, pR);
	const double gamma = eos.getc(1., 1.)*eos.getc(1., 1.);
	int itCounter = 0;	
	// Пытаюсь определить возможную конфигурацию решения, чтобы вернее выставить начальное приближение
	// Похоже, итерации нужны только в случаях "УВ+УВ" и "УВ + ВР", т.к. в случае ВР+ВР и ВР+вакуум есть 
	// аналитические решения для идеального газа
	//
	// Также вызывает вопрос последний тест Торо, где полученное решение отличается от его решения 
	// во втором знаке после запятой
	if(roL==roR && vL==vR && pL==pR) {
		res.type = RWRW;
		res.roL  = roL;
		res.roR  = roL;
		res.p    = pL;
		res.v	 = vL;
		return res;
	}
	if(roL==0.) {
		res.type = VacRW;
		res.roL  = 0.;
		res.roR  = roR;
		res.p	 = 0.;
		res.v    = vR - 2.*cR/(gamma-1.);
		return res;
	}
	if(roR==0.) {
		res.type = RWVac;
		res.roL  = roL;
		res.roR  = 0.;
		res.p	 = 0.;
		res.v    = vL + 2.*cL/(gamma-1.);
		return res;
	}
	if(2.*cL/(gamma-1) + 2*cR/(gamma-1.) < fabs(vL-vR)){
		res.type  = RWVacRW;
		res.roL	  = 0.;
		res.roR   = 0.;
		res.v     = 0.;
		res.p	  = 0.;
		return res;
	}

	double fLmin = fL(eos, pL, roL, vL, pL) + fR(eos, pL, roR, vR, pR) + vR-vL;
	double fRMax = fL(eos, pR, roL, vL, pL) + fR(eos, pR, roR, vR, pR) + vR-vL;
	// Начальное приближение
	//p = 0.5*(pL+pR);
	p=pL/2.;
	do {
		pPrev = p;
		p = pPrev - (fL(eos, pPrev, roL, vL, pL) + fR(eos, pPrev, roR, vR, pR) + vR - vL )/
			        (dfLdp(eos, pPrev, roL, vL, pL) + dfRdp(eos, pPrev, roR, vR, pR)); 
		if (p<=0.)
			p = TOL;
		itCounter++;
	} while (fabs(2*(p-pPrev)/(p+pPrev))>TOL);
	res.p   = p;
	res.v   = 0.5*(vL + vR) + 0.5*(fR(eos, p, roR, vR, pR) - fL(eos, p, roL, vL, pL));
	if( p<pL && p>pR) {
		res.type = RWSW;
		res.roL  = roL*pow(res.p/pL, 1./gamma); 
		res.roR  = roR*(res.p/pR + (gamma-1.)/(gamma+1.))/((gamma-1.)/(gamma+1.)*res.p/pR + 1.);
	} else if(p<=pL && p<=pR) {
		res.type = RWRW;
		res.roL  = roL*pow(res.p/pL, 1./gamma); 
		res.roR  = roR*pow(res.p/pR, 1./gamma); 
	} else if(p>pL && p<pR) {
		res.type = SWRW;
		res.roL  = roL*(res.p/pL + (gamma-1.)/(gamma+1.))/((gamma-1.)/(gamma+1.)*res.p/pL + 1.);
		res.roR  = roR*pow(res.p/pR, 1./gamma); 
	} else {
		res.type = SWSW;
		res.roL  = roL*(res.p/pL + (gamma-1.)/(gamma+1.))/((gamma-1.)/(gamma+1.)*res.p/pL + 1.);
		res.roR  = roR*(res.p/pR + (gamma-1.)/(gamma+1.))/((gamma-1.)/(gamma+1.)*res.p/pR + 1.);
	}
	return res;
}


double COutput::fL(FEOS& eos, double p, double roL, double vL, double pL) {
	double cL = eos.getc(roL, pL);	
	const double gamma = roL*cL*cL/pL;
	double f = 0.;
	if(p>pL) {
		double AL = 2./(gamma+1)/roL;
		double BL = (gamma-1.)/(gamma+1.)*pL;
		f = (p-pL) * sqrt(AL/(p+BL));
		return f;
	} else {		
		f = 2.*cL/(gamma-1.) * ( (pow(p/pL, (gamma-1.)/2./gamma)) - 1. );
		return f;	
	}
}

double COutput::dfLdp(FEOS& eos, double p, double roL, double vL, double pL) {
	double cL = eos.getc(roL, pL);	
	const double gamma = roL*cL*cL/pL;
	double dfdp = 0.;
	if (p>pL) {
		double AL = 2./(gamma+1)/roL;
		double BL = (gamma-1.)/(gamma+1.)*pL;
		dfdp = sqrt(AL/(p+BL)) * (1. - (p-pL)/2./(p+BL));
		return dfdp;
	}
	else {		
		dfdp = cL/pL/gamma*pow(p/pL, -(gamma+1)/2./gamma); 
		return dfdp;
	}
}

double COutput::fR(FEOS& eos, double p, double roR, double vR, double pR) {
	double cR = eos.getc(roR, pR);	
	const double gamma = roR*cR*cR/pR;
	double f = 0.;
	if(p>pR) {
		double AR = 2./(gamma+1)/roR;
		double BR = (gamma-1.)/(gamma+1.)*pR;
		f = (p-pR) * sqrt(AR/(p+BR));
		return f;
	} else {
		double cR = sqrt(gamma*pR/roR);
		f = 2.*cR/(gamma-1.) * ( (pow(p/pR, (gamma-1.)/2./gamma)) - 1. );
		return f;
	}
}

double COutput::dfRdp(FEOS& eos, double p, double roR, double vR, double pR) {
	double cR = eos.getc(roR, pR);	
	const double gamma = roR*cR*cR/pR;
	double dfdp = 0.;
	if (p>pR) {
		double AR = 2./(gamma+1)/roR;
		double BR = (gamma-1.)/(gamma+1.)*pR;
		dfdp = sqrt(AR/(p+BR)) * (1. - (p-pR)/2./(p+BR));
		return dfdp;
	} else {
		double cR = sqrt(gamma*pR/roR);
		dfdp = cR/pR/gamma*pow(p/pR, -(gamma+1)/2./gamma); 
		return dfdp;
	}
;}

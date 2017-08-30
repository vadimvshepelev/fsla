#include "..\\solver.h"
#include "..\\EOS\\EOSBin.h"


// Function finds the solution of Riemann problem for ideal gas with binary EOS
/*
RPSolutionPrimitive CSolver::solveRPEOSBin(double roL, double vL, double pL, double roR, double vR, double pR) {
	// Решаем нелинейное уравнение относительно давления методом касательных Ньютона
	RPSolutionPrimitive res; res.roL = 0.; res.roR=0.; res.v = 0.; res.p = 0.;
	double p = 0., pPrev = 0.;
	double TOL = 1.e-6;
	double gamma = 1.4;
	int itCounter = 0;
	double cL = 0., cR = 0.;
	if(roL!=0.) cL = sqrt(gamma*pL/roL);
	if(roR!=0.) cR = sqrt(gamma*pR/roR);
	// Пытаюсь определить возможную конфигурацию решения, чтобы вернее выставить начальное приближение
	// Похоже, итерации нужны только в случаях "УВ+УВ" и "УВ + ВР", т.к. в случае ВР+ВР и ВР+вакуум есть 
	// аналитические решения для идеального газа
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

	double fLmin = fL(pL, roL, vL, pL) + fR(pL, roR, vR, pR) + vR-vL;
	double fRMax = fL(pR, roL, vL, pL) + fR(pR, roR, vR, pR) + vR-vL;
	// Начальное приближение
	//p = 0.5*(pL+pR);
	p=pL/2.;
	do {
		pPrev = p;
		p = pPrev - (fL(pPrev, roL, vL, pL) + fR(pPrev, roR, vR, pR) + vR - vL )/
			        (dfLdp(pPrev, roL, vL, pL) + dfRdp(pPrev, roR, vR, pR)); 
		if (p<=0.)
			p = TOL;
		itCounter++;
	} while (fabs(2*(p-pPrev)/(p+pPrev))>TOL);
	res.p   = p;
	res.v   = 0.5*(vL + vR) + 0.5*(fR(p, roR, vR, pR) - fL(p, roL, vL, pL));
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
	//Тестирование (только без вакуума)
	//double L =  - fL(p, roL, vL, pL);
	//double R = fR(p, roR, vR, pR) + vR - vL;
	//double delta = fabs((L-R)/0.5/(L+R));
	//double LToro = - fL(res.p, roL, vL, pL);
	//double RToro = fR(res.p, roR, vR, pR) + vR - vL;
	//double deltaToro = fabs((LToro-RToro)/0.5/(LToro+RToro));
	return res;
}*/

/*
double CSolver::fLEOSBin(double p, double roL, double vL, double pL) {
	double gamma = 1.4;
	double f = 0.;
	if(p>pL) {
		double AL = 2./(gamma+1)/roL;
		double BL = (gamma-1.)/(gamma+1.)*pL;
		f = (p-pL) * sqrt(AL/(p+BL));
		return f;
	} else {
		double cL = sqrt(gamma*pL/roL);
		f = 2.*cL/(gamma-1.) * ( (pow(p/pL, (gamma-1.)/2./gamma)) - 1. );
		return f;	
	}
}*/

/*double CSolver::dfLdpEOSBin(double p, double roL, double vL, double pL) {
	 double gamma = 1.4;
	double dfdp = 0.;
	if (p>pL) {
		double AL = 2./(gamma+1)/roL;
		double BL = (gamma-1.)/(gamma+1.)*pL;
		dfdp = sqrt(AL/(p+BL)) * (1. - (p-pL)/2./(p+BL));
		return dfdp;
	}
	else {
		double cL = sqrt(gamma*pL/roL);
		dfdp = cL/pL/gamma*pow(p/pL, -(gamma+1)/2./gamma); 
		return dfdp;
	}
}*/

/*double CSolver::fREOSBin(double p, double roR, double vR, double pR) {
	double gamma = 1.4;
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
}*/

/*double CSolver::dfRdpEOSBin(double p, double roR, double vR, double pR) {
	double gamma = 1.4;
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
;}*/

/*
CVectorPrimitive CSolver::calcRPAnalyticalSolutionEOSBin(double roL, double vL, double pL, double roR, double vR, double pR, double x, double t){
	RPSolutionPrimitive res = solveRP(roL, vL, pL, roR, vR, pR);
	// V = (ro, v, p)T
	CVectorPrimitive V;
	double xi = x/t;
	double gamma = 1.4;
	double cL = 0., cR = 0.;
	if(roL!=0.) cL = sqrt(gamma*pL/roL);
	if(roR!=0.) cR = sqrt(gamma*pR/roR);
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
}*/


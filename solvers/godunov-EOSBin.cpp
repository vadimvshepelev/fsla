#include "..\\solver.h"
#include "..\\EOS\\EOSBin.h"


// Function finds the solution of Riemann problem for ideal gas with binary EOS
RPSolutionPrimitive CSolver::solveRPEOSBin(double roL, double vL, double pL, double roR, double vR, double pR) {
	// Решаем нелинейное уравнение относительно давления методом касательных Ньютона
	RPSolutionPrimitive res; 
	EOSBin* eos = task.eosBin;
	const double TOL = 1.e-6, gamma = eos->gamma, ro0 = eos->ro0, c0 = eos->c0, p0 = 1./gamma*ro0*ro0*c0;
	int itCounter = 0;
	double p = 0., pPrev = 0., cL = 0., cR = 0.;
	if(roL!=0.) cL = sqrt(gamma*(pL+p0)/roL);
	if(roR!=0.) cR = sqrt(gamma*(pR+p0)/roR);
	// Определяем возможную конфигурацию решения, чтобы вернее выставить начальное приближение
	// Похоже, итерации нужны только в случаях "УВ+УВ" и "УВ+ВР", т.к. в случае ВР+ВР и ВР+вакуум есть 
	// аналитические решения
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

	double fLmin = fLEOSBin(pL, roL, vL, pL) + fREOSBin(pL, roR, vR, pR) + vR-vL;
	double fRMax = fLEOSBin(pR, roL, vL, pL) + fREOSBin(pR, roR, vR, pR) + vR-vL;
	// Начальное приближение
	//p = 0.5*(pL+pR);
	p=pL/2.;
	do {
		pPrev = p;
		p = pPrev - (fLEOSBin(pPrev, roL, vL, pL) + fREOSBin(pPrev, roR, vR, pR) + vR - vL )/
			        (dfLdpEOSBin(pPrev, roL, vL, pL) + dfRdpEOSBin(pPrev, roR, vR, pR)); 
		// Отрицательные давления все же допускаются в случае с двучленным УРС!
		/*if (p<=0.)
			p = TOL;*/
		// Но если p+p0 < 0., это уже что-то неверное
		if (p+p0 <= 0.) {
			cout << "CSolver::SolveRPEOSBin error: sufficiently negative pressure occured!" << endl;
			exit(1);
		}			

		itCounter++;
	} while (fabs(2*(p-pPrev)/(p+pPrev))>TOL);
	res.p   = p;
	res.v   = 0.5*(vL + vR) + 0.5*(fREOSBin(p, roR, vR, pR) - fLEOSBin(p, roL, vL, pL));
	if( p<pL && p>pR) {
		res.type = RWSW;
		res.roL  = roL*pow((res.p+p0)/(pL+p0), 1./gamma); 
		res.roR  = roR*((res.p+p0)/(pR+p0) + (gamma-1.)/(gamma+1.))/((gamma-1.)/(gamma+1.)*(res.p+p0)/(pR+p0) + 1.);
	} else if(p<=pL && p<=pR) {
		res.type = RWRW;
		res.roL  = roL*pow((res.p+p0)/(pL+p0), 1./gamma); 
		res.roR  = roR*pow((res.p+p0)/(pR+p0), 1./gamma); 
	} else if(p>pL && p<pR) {
		res.type = SWRW;
		res.roL  = roL*((res.p+p0)/(pL+p0) + (gamma-1.)/(gamma+1.))/((gamma-1.)/(gamma+1.)*(res.p+p0)/(pL+p0) + 1.);
		res.roR  = roR*pow((res.p+p0)/(pR+p0), 1./gamma); 
	} else {
		res.type = SWSW;
		res.roL  = roL*((res.p+p0)/(pL+p0) + (gamma-1.)/(gamma+1.))/((gamma-1.)/(gamma+1.)*(res.p+p0)/(pL+p0) + 1.);
		res.roR  = roR*((res.p+p0)/(pR+p0) + (gamma-1.)/(gamma+1.))/((gamma-1.)/(gamma+1.)*(res.p+p0)/(pR+p0) + 1.);
	}
	//Тестирование (только без вакуума)
	//double L =  - fL(p, roL, vL, pL);
	//double R = fR(p, roR, vR, pR) + vR - vL;
	//double delta = fabs((L-R)/0.5/(L+R));
	//double LToro = - fL(res.p, roL, vL, pL);
	//double RToro = fR(res.p, roR, vR, pR) + vR - vL;
	//double deltaToro = fabs((LToro-RToro)/0.5/(LToro+RToro));
	return res;
}


double CSolver::fLEOSBin(double p, double roL, double vL, double pL) {
	EOSBin* eos = task.eosBin;
	const double gamma = eos->gamma, ro0 = eos->ro0, c0 = eos->c0, p0 = 1./gamma*ro0*ro0*c0;
	double f = 0.;
	if(p>pL) {
		double AL = 2./(gamma+1)/roL;
		double BL = (gamma-1.)/(gamma+1.)*(pL+p0);
		f = (p-pL) * sqrt(AL/(p+p0+BL));
		return f;
	} else {
		double cL = sqrt(gamma*(pL+p0)/roL);
		f = 2.*cL/(gamma-1.) * ( (pow((p+p0)/(pL+p0), (gamma-1.)/2./gamma)) - 1. );
		return f;	
	}
}

double CSolver::dfLdpEOSBin(double p, double roL, double vL, double pL) {
    EOSBin* eos = task.eosBin;
	const double gamma = eos->gamma, ro0 = eos->ro0, c0 = eos->c0, p0 = 1./gamma*ro0*ro0*c0;
	double dfdp = 0.;
	if (p>pL) {
		double AL = 2./(gamma+1)/roL;
		double BL = (gamma-1.)/(gamma+1.)*(pL+p0);
		dfdp = sqrt(AL/(p+p0+BL)) * (1. - (p-pL)/2./(p+p0+BL));
		return dfdp;
	}
	else {
		double cL = sqrt(gamma*(pL+p0)/roL);
		dfdp = cL/(pL+p0)/gamma*pow((p+p0)/(pL+p0), -(gamma+1)/2./gamma); 
		return dfdp;
	} 
}

double CSolver::fREOSBin(double p, double roR, double vR, double pR) {
	const double gamma = task.eosBin->gamma, ro0 = task.eosBin->ro0, c0 = task.eosBin->c0, p0 = 1./gamma*ro0*ro0*c0;
	double f = 0.;
	if(p>pR) {
		double AR = 2./(gamma+1)/roR;
		double BR = (gamma-1.)/(gamma+1.)*(pR+p0);
		f = (p-pR) * sqrt(AR/(p+p0+BR));
		return f;
	} else {
		double cR = sqrt(gamma*(pR+p0)/roR);
		f = 2.*cR/(gamma-1.) * ( (pow((p+p0)/(pR+p0), (gamma-1.)/2./gamma)) - 1. );
		return f;
	} 
}

double CSolver::dfRdpEOSBin(double p, double roR, double vR, double pR) {
	const double gamma = task.eosBin->gamma, ro0 = task.eosBin->ro0, c0 = task.eosBin->c0, p0 = 1./gamma*ro0*ro0*c0;
	double dfdp = 0.;
	if (p>pR) {
		double AR = 2./(gamma+1)/roR;
		double BR = (gamma-1.)/(gamma+1.)*(pR+p0);
		dfdp = sqrt(AR/(p+p0+BR)) * (1. - (p-pR)/2./(p+p0+BR));
		return dfdp;
	} else {
		double cR = sqrt(gamma*(pR+p0)/roR);
		dfdp = cR/(pR+p0)/gamma*pow((p+p0)/(pR+p0), -(gamma+1)/2./gamma); 
		return dfdp;
	} 
;}

CVectorPrimitive CSolver::calcRPAnalyticalSolutionEOSBin(double roL, double vL, double pL, double roR, double vR, double pR, double x, double t){
	
	
	
	



	/// DEBUG
	//RPSolutionPrimitive res0 = solveRP(roL, vL, pL, roR, vR, pR);
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	RPSolutionPrimitive res = solveRPEOSBin(roL, vL, pL, roR, vR, pR);
	// V = (ro, v, p)T
	CVectorPrimitive V;
	const double TOL = 1.e-6, gamma = task.eosBin->gamma, ro0 = task.eosBin->ro0, c0 = task.eosBin->c0, p0 = 1./gamma*ro0*ro0*c0;
	double xi = x/t, cL = 0., cR = 0.;
	if(roL!=0.) cL = sqrt(gamma*(pL+p0)/roL);
	if(roR!=0.) cR = sqrt(gamma*(pR+p0)/roR);
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
			V.p  = (pR+p0)*pow(2./(gamma+1.) - (gamma-1.)/(gamma+1.)/cR*(vR-xi), 2.*gamma/(gamma-1.)); 
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
			V.p  = (pL+p0)*pow(2./(gamma+1.)+(gamma-1.)/(gamma+1.)/cL*(vL-xi), 2.*gamma/(gamma-1.));
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
			V.p  = (pL+p0)*pow(2./(gamma+1.)+(gamma-1.)/(gamma+1.)/cL*(vL-xi), 2.*gamma/(gamma-1.));
		} else if (xi<=xiTailR) {
			V.ro = 0.;
			V.v  = 0.;
			V.p  = 0.;
		} else if (xi<xiHeadR) {
			V.ro = roR*pow(2./(gamma+1.) - (gamma-1.)/(gamma+1.)/cR*(vR-xi), 2./(gamma-1.));
			V.v  = 2./(gamma+1)*(-cR + (gamma-1.)/2.*vR + xi);
			V.p  = (pR+p0)*pow(2./(gamma+1.) - (gamma-1.)/(gamma+1.)/cR*(vR-xi), 2.*gamma/(gamma-1.)); 
		} else {
			V.ro = roR;
			V.v  = vR;
			V.p  = pR;
		}
		return V;
	}
	double cLLocal = sqrt(gamma*(res.p+p0)/res.roL), cRLocal = sqrt(gamma*(res.p+p0)/res.roR);
	// Если не вакуум. Пусть точка слева от контактного разрыва (xiContact = res.v)
	if(xi<res.v) {
		if(res.type == SWSW || res.type == SWRW) { 
			xiFront = vL - cL*sqrt((gamma+1.)/2./gamma*(res.p+p0)/(pL+p0) + (gamma-1.)/2./gamma);
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



				double c = 2./(gamma+1.) * (cL + (gamma-1.)/2.*(vL-xi));





				V.p  = (pL+p0)*pow(2./(gamma+1.)+(gamma-1.)/(gamma+1.)/cL*(vL-xi), 2.*gamma/(gamma-1.)) - p0;
			}
		} 
	//Пусть точка справа от контактного разрыва (xiContact = res.v)
	} else {
		if(res.type == RWSW || res.type == SWSW) {
			xiFront = vR + cR*sqrt((gamma+1.)/2./gamma*(res.p+p0)/(pR+p0) + (gamma-1.)/2./gamma);
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
				V.p  = (pR+p0)*pow(2./(gamma+1.) - (gamma-1.)/(gamma+1.)/cR*(vR-xi), 2.*gamma/(gamma-1.)) - p0;
			}
		}
	}
	return V; 
}

Vector4 CSolver::calcGodunovFluxEOSBin(double roL, double rouL, double roEL, double roR, double rouR, double roER) {
	EOSBin eos = *(task.eosBin);	
	double uL = rouL/roL, uR = rouR/roR;
	double eL = roEL/roL - .5*uL*uL, eR = roER/roR - .5*uR*uR; 
	double pL = eos.getp(roL, eL), pR = eos.getp(roR, eR);
	CVectorPrimitive res = calcRPAnalyticalSolutionEOSBin(roL, uL, pL, roR, uR, pR, 0., .01);
	
	
	
	
	
	
	
	//CVectorPrimitive res = calcRPAnalyticalSolutionEOSBin(1., -2., .4, 1., 2., .4, 0., .2);
	
	
	
	
	
	
	
	
	double E = eos.gete(res.ro, res.p) + .5*res.v*res.v;
	Vector4 FGodunov = Vector4(res.ro*res.v, res.ro*res.v*res.v + res.p, res.v*(res.ro*E+res.p), 0.);
	return FGodunov;
}

void CSolver::calcHydroStageGodunovEOSBin(double t, double tau) {
	EOSBin eos = *(task.eosBin);
	const double gamma = eos.gamma;
	double E=0., h=ms[1].x-ms[0].x;
	const int nSize = ms.getSize(); 
	int i=0; 	
	Vector4 L = Vector4::ZERO, R = Vector4::ZERO, D = Vector4::ZERO, V = Vector4::ZERO;	
	// Значения на границе (в фиктивных ячейках)
	double roLB=0., vLB = 0., ELB=0., roRB=0., vRB=0., ERB=0.;
	// Задаем граничные условия
	// Transmissive left boundary
	roLB = ms[0].ro; vLB = ms[0].v;	ELB = ms[0].e + .5*ms[0].v*ms[0].v;
	// Reflective left boundary
	//roLB = ms[0].ro; vLB = -ms[0].v;	ELB = ms[0].e + .5*ms[0].v*ms[0].v;
	// Transmissive right boundary
	roRB = ms[nSize-1].ro; vRB = ms[nSize-1].v; ERB = ms[nSize-1].e + .5*ms[nSize-1].v*ms[nSize-1].v;
	for(i=0; i<ms.getSize(); i++) {
		Node &n = ms[i];



		if(i == 55) {
			double qq = 0.;
		}






		if(i==0)
			//n.F = calcGodunovFluxEOSBin(roLB, roLB*vLB, roLB*ELB, n.W[0], n.W[1], n.W[2]);
			n.F = calcHLLFluxEOSBin(roLB, roLB*vLB, roLB*ELB, n.W[0], n.W[1], n.W[2]);
		else
			//n.F = calcGodunovFluxEOSBin(ms[i-1].W[0], ms[i-1].W[1], ms[i-1].W[2], n.W[0], n.W[1], n.W[2]);	
			n.F = calcHLLFluxEOSBin(ms[i-1].W[0], ms[i-1].W[1], ms[i-1].W[2], n.W[0], n.W[1], n.W[2]);	
	}
	//ms[nSize].F = calcGodunovFluxEOSBin(ms[nSize-1].W[0], ms[nSize-1].W[1], ms[nSize-1].W[2], roRB, roRB*vRB, roRB*ERB);
	ms[nSize].F = calcHLLFluxEOSBin(ms[nSize-1].W[0], ms[nSize-1].W[1], ms[nSize-1].W[2], roRB, roRB*vRB, roRB*ERB);
	// Main cycle
	for(i=0; i<nSize; i++) {				
		ms[i].W_temp = ms[i].W - tau/h*(ms[i+1].F-ms[i].F);
	}
	for(i=0; i<nSize; i++) {
		// Обновляем консервативные переменные
		Node& n=ms[i]; 
		n.W[0] = n.W_temp[0];
		n.W[1] = n.W_temp[1];
		n.W[2] = n.W_temp[2];
		// Обновляем все переменные
		n.ro = n.W[0];
		n.v  = n.W[1]/n.ro;
		n.e  = n.W[2]/n.ro - 0.5*n.v*n.v;
		n.p  = eos.getp(n.ro, n.e);
		n.C  = eos.getC(n.ro, n.e);
	}
}
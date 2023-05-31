#include "C1DMethod.h"

#include <algorithm>
#include <ranges>
#include <utility>

#include "_vector4.h"
#include "solver.h"


CExactRiemannSolver exrslv;

void C1DGodunovMethodMillerPuckett::calc(C1DProblem& pr, FEOSMieGruneisen& eos, C1DField& fld) {
	double roL = 0., uL = 0., eL = 0., pL = 0.,
		   roR = 0., uR = 0., eR = 0., pR = 0., 
		   E = 0.;
	double dx = fld.dx, t = fld.t, dt = fld.dt;
	int imin = fld.imin, imax = fld.imax;
//	auto U = fld.U;
//	auto newU = fld.newU;
//	auto F = fld.F;
	auto&& U = fld.U;
	auto&& newU = fld.newU;
	auto&& F = fld.F;
	// TODO: проверить на скорость выполнения операций, сравнить с реализацией через тип Vector4 -- если не медленнее, то в дальнейшем избавиться от Vector4 везде
	int i=0;
	C1DVectorPrimitive res;
	// Потоки считаем по алгоритму решения задаче о распаде разрыва для УРС Ми-Грюнайзена из работы [Miller, Puckett]
	for(i=imin; i<=imax; i++) {
		roL = U[i-1][0]; uL = U[i-1][1]/roL; eL = U[i-1][2]/roL - .5*uL*uL; pL = eos.getp(roL, eL);
		roR = U[i][0];   uR = U[i][1]/roR;   eR = U[i][2]/roR - .5*uR*uR;   pR = eos.getp(roR, eR);								
		res = calcRPExactMillerPuckett(eos, roL, uL, pL, roR, uR, pR);
		assert(res.ro!=0.);
		E = eos.gete(res.ro, res.p) + .5*res.u*res.u;
		F[i] = {
			res.ro * res.u,
			res.ro * res.u * res.u + res.p,
			res.u * (res.ro * E + res.p),
			0.
		};
		//n.W_temp = n.W - dt/h*(Fp-Fm);
	}
	for(i=imin; i<imax; i++) {
		for(int counter=0; counter<3; counter++) newU[i][counter] = U[i][counter] - dt/dx*(F[i+1][counter] - F[i][counter]);
	}
	for(i=imin; i<imax; i++) {
		for(int counter=0; counter<3; counter++) U[i][counter] = newU[i][counter];
	}	
	pr.setbcs(fld.U);
}

double C1DGodunovMethodMillerPuckett::calcdt(C1DProblem& pr, FEOSMieGruneisen& eos, C1DField& fld) {
	auto&& U = fld.U;
	int imin = fld.imin, imax = fld.imax;
	double ro = U[imin][0], u = U[imin][1]/ro, e = U[imin][2]/ro-.5*u*u, p=eos.getp(ro,e), c = eos.getc(ro, p);
	vector<double> x = fld.x;	
	double umax = max(fabs(u), max(fabs(u-c), fabs(u+c))); 	
	double dt1 = (x[imin+1]-x[imin])/umax;
	double dt2 = 0.;
	for(int i=imin; i<imax; i++) {		
		ro = U[i][0]; u = U[i][1]/ro; e = U[i][2]/ro-.5*u*u, p=eos.getp(ro,e), c = eos.getc(ro, p); 
		umax = fabs(u) + c; 
		dt2 = (x[i+1]-x[i])/umax;
		if(dt1 > dt2) dt1 = dt2;
	}
	return pr.cfl*dt1;
}

C1DVectorPrimitive C1DGodunovMethodMillerPuckett::calcRPExactMillerPuckett(FEOSMieGruneisen& eos, double roL, double uL, double pL, double roR, double uR, double pR) {
	C1DVectorPrimitive V;
	const double gammaL = eos.getG(roL), gammaR = eos.getG(roR);
	// Первое замечание по физике: почему нет "нормального" нулевого давления? То давление, которое "должно" быть нулевым по логике вещей, существенно отрицательно.
	// const double K0S = eos.ro0*eos.getc(eos.ro0, eos.getp0(eos.ro0));
	// Исправляем:
	const double K0S = eos.ro0*eos.getc(eos.ro0, eos.getp(eos.ro0, 2.e7));
	const double eL = eos.gete(roL, pL), eR = eos.gete(roR, pR), cL = eos.getc(roL, pL), cR = eos.getc(roR, pR); 
	double _e = 0.;
	const double KSL = roL*cL*cL, KSR = roL*cR*cR;
	const double KSPrimeL = eos.getKSPrime(roL, eL), KSPrimeR= eos.getKSPrime(roR, eR); // Аккуратно посчитать через первые и вторые производные, просто это надо чуть времени
	double a0L = cL, a1L = (KSPrimeL + 1.)/4., a0R = cR, a1R = (KSPrimeR + 1.)/4.; 
	double b0 = pL-pR + roL*uL*(a0L+a1L*uL) + roR*uR*(a0L-a1R*uR),
		   b1 = -roL*(a0L+2.*a1L*uL) - roR*(a0R-2.*a1R*uR),
		   b2 = roL*a1L - roR*a1R;
	double uCD = (-b1 - sqrt(b1*b1-4.*b0*b2))/2./b2,
		   pCD = pL + roR*a0R*(uCD - uR) + roR*a1R*(uCD - uR)*(uCD - uR),
		   uS = 0.;		   
	double pLHat = 0., pRHat = 0.;
	if(uL>uR) {
		pLHat = pL + roL*(uL-uR)*(a0L + a1L*(uL-uR)); 
		pRHat = pR + roR*(uL-uR)*(a0R + a1R*(uL-uR));
	} else {
		pLHat = pL + roL*a0L*(uL-uR);
		pRHat = pR + roR*a0R*(uL-uR);
	}
	if(pL>pRHat) a1L = 0.; 
	if(pR>pLHat) a1R = 0.;	
	if(uCD>0.) {
        // U(0) belongs 'L' part of material, initial discontinuity moves in the right direction from x=0 point 
		if(pCD>pL) {
			// Left shock
			uS = uL - a0L + a1L*(uCD - uL);
			if(uS>0.) {
				V.ro = roL* (-a0L + a1L*(uCD - uL)) / (uL - a0L + a1L*(uCD - uL) - uCD);
				V.u = uCD;
				//_e = eL + .5*(pL + pCD)*(1./roL-1./V.ro);
				//V.p = eos.getp(V.ro, _e);
				V.p = pCD;			
			} else {
				V.ro = roL;
				V.u = uL;
				V.p = pL;
			}
		} else {
		    // Left rarefaction
			double roCDL = roL/(1.+(pL-pCD)/KSL);
			double cCDL = a0L*roL/roCDL;
			double sigmaLWave = (a0L-uL)/(a0L-uL-uCD-cCDL);
			double sigmaL = min(1., max(0., sigmaLWave));
			V.ro = (1.-sigmaL)*roL + sigmaL*roCDL;
			V.u = (1.-sigmaL)*uL + sigmaL*uCD;
			V.p = (1.-sigmaL)*pL + sigmaL*pCD;


			// И тут вопрос, а какой будет энергия? Если мы получим V.e (энергию в нуле х) интерполяцией и сравним с V.p, что мы получим? Это повод для теста.


		}
	} else {
		// U(0) belongs 'R' part of material, initial discontinuity moves in the left direction from x=0 point 
		// Right shock
		if(pCD>pR) {
			// Left shock
			uS = uR + a0R + a1R*(uCD - uR);
			if(uS>0.) {
				V.ro = roR* (a0R + a1R*(uCD - uR)) / (a0R + a1R*(uCD - uR) + uR - uCD);
				V.u = uCD;
				//_e = eL + .5*(pL + pCD)*(1./roL-1./V.ro);
				//V.p = eos.getp(V.ro, _e);
				V.p = pCD;			
			} else {
				V.ro = roR;
				V.u = uR;
				V.p = pR;
			}
		} else {
		    // Right rarefaction
			double roCDR = roL/(1.+(pR-pCD)/KSR);
			double cCDR = a0R*roR/roCDR;
			double sigmaRWave = (a0R-uR)/(a0R-uR-uCD-cCDR);
			double sigmaR = min(1., max(0., sigmaRWave));
			V.ro = (1.-sigmaR)*roR + sigmaR*roCDR;
			V.u = (1.-sigmaR)*uL + sigmaR*uCD;
			V.p = (1.-sigmaR)*pL + sigmaR*pCD;
		}
	}
	return V;
}


// Function finds resulting ro, u, p and wave structure of Riemann problem
RPValues calcValues(FEOS& eos, double roL, double uL, double pL, double roR, double uR, double pR) {
	// Решаем нелинейное уравнение относительно давления методом касательных Ньютона
	RPValues res; 
	double p = 0., pPrev = 0.;
	const double TOL = 1.e-6;
	const double gamma = eos.getc(1., 1.)*eos.getc(1., 1.);
	int itCounter = 0;
	double cL = 0., cR = 0.;
	if(roL!=0.) cL = eos.getc(roL, pL);
	if(roR!=0.) cR = eos.getc(roR, pR);;
	// Пытаюсь определить возможную конфигурацию решения, чтобы вернее выставить начальное приближение
	// Похоже, итерации нужны только в случаях "УВ+УВ" и "УВ + ВР", т.к. в случае ВР+ВР и ВР+вакуум есть 
	// аналитические решения для идеального газа
	//
	// Также вызывает вопрос последний тест Торо, где полученное решение отличается от его решения 
	// во втором знаке после запятой
	if(roL==roR && uL==uR && pL==pR) {
		res.type = RPWaveConfig::rwrw;
		res.roL  = roL;
		res.roR  = roL;
		res.p    = pL;
		res.u	 = uL;
		return res;
	}
	if(roL==0.) {
		res.type = RPWaveConfig::vacrw;
		res.roL  = 0.;
		res.roR  = roR;
		res.p	 = 0.;
		res.u    = uR - 2.*cR/(gamma-1.);
		return res;
	}
	if(roR==0.) {
		res.type = RPWaveConfig::rwvac;
		res.roL  = roL;
		res.roR  = 0.;
		res.p	 = 0.;
		res.u    = uL + 2.*cL/(gamma-1.);
		return res;
	}
	if(2.*cL/(gamma-1) + 2*cR/(gamma-1.) < fabs(uL-uR)){
		res.type  = RPWaveConfig::rwvacrw;
		res.roL	  = 0.;
		res.roR   = 0.;
		res.u     = 0.;
		res.p	  = 0.;
		return res;
	}

	double fLmin = fL(eos, pL, roL, uL, pL) + fR(eos, pL, roR, uR, pR) + uR-uL;
	double fRMax = fL(eos, pR, roL, uL, pL) + fR(eos, pR, roR, uR, pR) + uR-uL;
	// Начальное приближение
	//p = 0.5*(pL+pR);
	p=pL/2.;
	do {
		pPrev = p;
		p = pPrev - (fL(eos, pPrev, roL, uL, pL) + fR(eos, pPrev, roR, uR, pR) + uR - uL )/
			        (dfLdp(eos, pPrev, roL, uL, pL) + dfRdp(eos, pPrev, roR, uR, pR)); 
		if (p<=0.)
			p = TOL;
		itCounter++;
	} while (fabs(2*(p-pPrev)/(p+pPrev))>TOL);
	res.p   = p;
	res.u   = 0.5*(uL + uR) + 0.5*(fR(eos, p, roR, uR, pR) - fL(eos, p, roL, uL, pL));
	if( p<pL && p>pR) {
		res.type = RPWaveConfig::rwsw;
		res.roL  = roL*pow(res.p/pL, 1./gamma); 
		res.roR  = roR*(res.p/pR + (gamma-1.)/(gamma+1.))/((gamma-1.)/(gamma+1.)*res.p/pR + 1.);
	} else if(p<=pL && p<=pR) {
		res.type = RPWaveConfig::swsw;
		res.roL  = roL*pow(res.p/pL, 1./gamma); 
		res.roR  = roR*pow(res.p/pR, 1./gamma); 
	} else if(p>pL && p<pR) {
		res.type = RPWaveConfig::swsw;
		res.roL  = roL*(res.p/pL + (gamma-1.)/(gamma+1.))/((gamma-1.)/(gamma+1.)*res.p/pL + 1.);
		res.roR  = roR*pow(res.p/pR, 1./gamma); 
	} else {
		res.type = RPWaveConfig::swsw;
		res.roL  = roL*(res.p/pL + (gamma-1.)/(gamma+1.))/((gamma-1.)/(gamma+1.)*res.p/pL + 1.);
		res.roR  = roR*(res.p/pR + (gamma-1.)/(gamma+1.))/((gamma-1.)/(gamma+1.)*res.p/pR + 1.);
	}
	return res;
}

double fL(FEOS& eos, double p, double roL, double uL, double pL) {
	const double gamma = eos.getc(1., 1.)*eos.getc(1.,1.);
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
}

double dfLdp(FEOS& eos, double p, double roL, double vL, double pL) {
	const double gamma = eos.getc(1., 1.)*eos.getc(1.,1.);
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
}

double fR(FEOS& eos, double p, double roR, double vR, double pR) {
	const double gamma = eos.getc(1., 1.)*eos.getc(1.,1.);
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

double dfRdp(FEOS& eos, double p, double roR, double vR, double pR) {
	const double gamma = eos.getc(1., 1.)*eos.getc(1.,1.);
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
}

C1DVectorPrimitive calcSolution(FEOS& eos, double roL, double uL, double pL, double roR, double uR, double pR, double x, double t){
	RPValues res = calcValues(eos, roL, uL, pL, roR, uR, pR);
	// V = (ro, v, p)T
	C1DVectorPrimitive V;
	double xi = x/t;
	const double gamma = eos.getc(1., 1.)*eos.getc(1., 1.);
	double cL = 0., cR = 0.;
	if(roL!=0.) cL = eos.getc(roL, pL);
	if(roR!=0.) cR = eos.getc(roR, pR);;
	double xiFront=0., xiHead=0., xiTail=0., xiHeadL=0., xiTailL=0., xiHeadR=0., xiTailR=0.;
	// Если вакуум
	if(res.type == RPWaveConfig::vacrw) { 
		xiHead = uR + cR;
		xiTail = uR - 2.*cR/(gamma-1.);
		if(xi<=xiTail) {
			V.ro = 0.;
			V.u  = uR - 2.*cR/(gamma-1.);
			V.p  = 0.;
		} else if (xi<xiHead) {
			V.ro = roR*pow(2./(gamma+1.) - (gamma-1.)/(gamma+1.)/cR*(uR-xi), 2./(gamma-1.));
			V.u  = 2./(gamma+1)*(-cR + (gamma-1.)/2.*uR + xi);
			V.p  = pR*pow(2./(gamma+1.) - (gamma-1.)/(gamma+1.)/cR*(uR-xi), 2.*gamma/(gamma-1.)); 
		} else {
			V.ro = roR;
			V.u  = uR;
			V.p  = pR;
		}
		return V;
	}
	if(res.type == RPWaveConfig::rwvac) {
		xiHead = uL - cL;
		xiTail = uL + 2.*cL/(gamma-1.);
		if(xi>=xiTail) {
			V.ro = 0.;
			V.u  = 0.;
			V.p  = 0.;
		} else if (xi>xiHead) {
			V.ro = roL*pow(2./(gamma+1.)+(gamma-1.)/(gamma+1.)/cL*(uL-xi), 2./(gamma-1.));
			V.u  = 2./(gamma+1)*(cL + (gamma-1.)/2.*uL + xi);
			V.p  = pL*pow(2./(gamma+1.)+(gamma-1.)/(gamma+1.)/cL*(uL-xi), 2.*gamma/(gamma-1.));
		} else {
			V.ro = roL;
			V.u = uL;
			V.p = pL;
		}
		return V;
	}
	if(res.type ==  RPWaveConfig::rwvacrw) {
		xiHeadL = uL - cL;
		xiTailL = uL + 2.*cL/(gamma-1.);
		xiHeadR = uR + cR;
		xiTailR = uR - 2.*cR/(gamma-1.);
		if(xi<=xiHeadL) {
			V.ro = roL;
			V.u  = uL;
			V.p  = pL;
		} else if (xi<xiTailL) {
			V.ro = roL*pow(2./(gamma+1.)+(gamma-1.)/(gamma+1.)/cL*(uL-xi), 2./(gamma-1.));
			V.u  = 2./(gamma+1)*(cL + (gamma-1.)/2.*uL + xi);
			V.p  = pL*pow(2./(gamma+1.)+(gamma-1.)/(gamma+1.)/cL*(uL-xi), 2.*gamma/(gamma-1.));
		} else if (xi<=xiTailR) {
			V.ro = 0.;
			V.u  = 0.;
			V.p  = 0.;
		} else if (xi<xiHeadR) {
			V.ro = roR*pow(2./(gamma+1.) - (gamma-1.)/(gamma+1.)/cR*(uR-xi), 2./(gamma-1.));
			V.u  = 2./(gamma+1)*(-cR + (gamma-1.)/2.*uR + xi);
			V.p  = pR*pow(2./(gamma+1.) - (gamma-1.)/(gamma+1.)/cR*(uR-xi), 2.*gamma/(gamma-1.)); 
		} else {
			V.ro = roR;
			V.u  = uR;
			V.p  = pR;
		}
		return V;
	}
	double cLLocal = eos.getc(res.roL, res.p), cRLocal = eos.getc(res.roR, res.p);
	// Если не вакуум. Пусть точка слева от контактного разрыва (xiContact = res.u)
	if(xi<res.u) {
		if(res.type ==  RPWaveConfig::swsw || res.type ==  RPWaveConfig::swrw) { 
			xiFront = uL - cL*sqrt((gamma+1.)/2./gamma*res.p/pL + (gamma-1.)/2./gamma);
			if(xi<xiFront) {
				V.ro = roL;
				V.u  = uL;
				V.p  = pL;
			} else {
				V.ro = res.roL;
				V.u = res.u;
				V.p = res.p;
			}
		} else if (res.type ==  RPWaveConfig::rwsw || res.type == RPWaveConfig::rwrw) {
			xiHead = uL-cL;
			xiTail = res.u-cLLocal;
			if(xi<=xiHead) {
				V.ro = roL;
				V.u  = uL;
				V.p  = pL;
			} else if(xi>=xiTail) {
				V.ro = res.roL;
				V.u  = res.u;
				V.p  = res.p;
			} else {
				V.ro = roL*pow(2./(gamma+1.)+(gamma-1.)/(gamma+1.)/cL*(uL-xi), 2./(gamma-1.));
				V.u  = 2./(gamma+1)*(cL + (gamma-1.)/2.*uL + xi);
				V.p  = pL*pow(2./(gamma+1.)+(gamma-1.)/(gamma+1.)/cL*(uL-xi), 2.*gamma/(gamma-1.));
			}
		} 
	//Пусть точка справа от контактного разрыва (xiContact = res.v)
	} else {
		if(res.type ==  RPWaveConfig::rwsw || res.type ==  RPWaveConfig::swsw) {
			xiFront = uR + cR*sqrt((gamma+1.)/2./gamma*res.p/pR + (gamma-1.)/2./gamma);
			if(xi>xiFront) {
				V.ro = roR;
				V.u  = uR;
				V.p  = pR;
			} else {
				V.ro = res.roR;
				V.u  = res.u;
				V.p  = res.p;
			}
		} else if(res.type ==  RPWaveConfig::rwrw || res.type ==  RPWaveConfig::swrw) {
			xiHead = uR + cR;
			xiTail = res.u + cRLocal;
			if(xi >= xiHead) {
				V.ro = roR;
				V.u  = uR;
				V.p  = pR;
			} else if (xi <= xiTail) {
				V.ro = res.roR;
				V.u  = res.u;
				V.p  = res.p;
			} else {
				V.ro = roR*pow(2./(gamma+1.) - (gamma-1.)/(gamma+1.)/cR*(uR-xi), 2./(gamma-1.));
				V.u  = 2./(gamma+1)*(-cR + (gamma-1.)/2.*uR + xi);
				V.p  = pR*pow(2./(gamma+1.) - (gamma-1.)/(gamma+1.)/cR*(uR-xi), 2.*gamma/(gamma-1.));
			}
		}
	}
	return V;
}


double C1DGodunovTypeMethod::calcdt(C1DProblem& pr, FEOS& eos, C1DField& fld) {
	auto&& U = fld.U;
	int imin = fld.imin, imax = fld.imax;
	double ro = U[imin][0], u = U[imin][1]/ro, e = U[imin][2]/ro-.5*u*u, p=eos.getp(ro,e), c = eos.getc(ro, p);
	vector<double> x = fld.x;	
	double umax = max(fabs(u), max(fabs(u-c), fabs(u+c))); 	
	double dt1 = (x[imin+1]-x[imin])/umax;
	double dt2 = 0.;
	for (int i = imin; i < imax; ++ i) {
		ro = U[i][0]; u = U[i][1]/ro; e = U[i][2]/ro-.5*u*u, p=eos.getp(ro,e), c = eos.getc(ro, p); 
		umax = fabs(u) + c; 
		dt2 = (x[i+1]-x[i])/umax;
		if (dt1 > dt2) dt1 = dt2;
	}
	return pr.cfl * dt1;
}


void C1DGodunovTypeMethod::calcFluxField(
		C1DProblem& pr, FEOS& eos, C1DField& fld) {
//	double roL = 0., uL = 0., eL = 0., pL = 0.,
//		   roR = 0., uR = 0., eR = 0., pR = 0.,
//		   E = 0.;
	int imin = fld.imin, imax = fld.imax;
	auto&& U = fld.U;
	auto&& F = fld.F;
	// TODO: проверить на скорость выполнения операций, сравнить
	// с реализацией через тип Vector4 -- если не медленнее,
	// то в дальнейшем избавиться от Vector4 везде
	int i = 0;
	// Потоки считаем по алгоритму решения задаче о распаде разрыва
	// для УРС Ми-Грюнайзена из работы [Miller, Puckett]
	for (i = imin; i <= imax; ++ i) {
		if (i == 53) {
			double q = 1.;
		}

		F[i] = rslv.calcFlux(
					eos, U[i-1][0], U[i-1][1], U[i-1][2],
					U[i][0], U[i][1], U[i][2]);
		// n.W_temp = n.W - dt/h*(Fp-Fm);
	}
}


void C1DGodunovTypeMethod::calc(C1DProblem& pr, FEOS& eos, C1DField& fld) {
	double dx = fld.dx, t = fld.t, dt = fld.dt;
	int imin = fld.imin, imax = fld.imax;
	auto&& U = fld.U;
	auto&& newU = fld.newU;
	auto&& F = fld.F;
	int i = 0;

	calcFluxField(pr, eos, fld);

	for (i = imin; i < imax; ++ i) {
		for(int counter=0; counter<3; counter++)
			newU[i][counter] = U[i][counter]
					- dt/dx*(F[i+1][counter] - F[i][counter]);
	}
	for (i = imin; i < imax; ++ i) {
		for(int counter=0; counter<3; counter++)
			U[i][counter] = newU[i][counter];
	}	
	pr.setbcs(fld.U);
}


double C1DGodunovTypeMethodVacuum::calcdt(C1DProblem& pr, FEOS& eos, C1DField& fld) {
	auto&& U = fld.U;
	int imin = fld.imin, imax = fld.imax;
	double ro = U[imax-1][0], u = U[imax-1][1]/ro, e = U[imax-1][2]/ro-.5*u*u, p=eos.getp(ro,e), c = eos.getc(ro, p);
	vector<double> x = fld.x;	
	const double gamma = eos.getc(1., 1.)*eos.getc(1., 1.);
	double umax = max(fabs(u), max(fabs(u-c), fabs(u+c))); 	
	double dt1 = (x[imax]-x[imax-1])/umax;
	double dt2 = 0.;
	for(int i=imin; i<imax; i++) {
		ro = U[i][0]; u = U[i][1]/ro; e = U[i][2]/ro-.5*u*u, p=eos.getp(ro,e), c = eos.getc(ro, p); 
		if(ro == 0.) 
			continue;
		umax = fabs(u) + c; 
		dt2 = (x[i+1]-x[i])/umax;
		if(dt1 > dt2) dt1 = dt2;
	}
	return pr.cfl*dt1;
}

void C1DGodunovTypeMethodVacuum::calc(C1DProblem& pr, FEOS& eos, C1DField& fld) {
	double roL = 0., uL = 0., eL = 0., pL = 0.,
		   roR = 0., uR = 0., eR = 0., pR = 0., 
		   E = 0.;
	double dx = fld.dx, t = fld.t, dt = fld.dt;
	int imin = fld.imin, imax = fld.imax;
	auto&& U = fld.U;
	auto&& newU = fld.newU;
	auto&& F = fld.F;
	vector<double>& x = fld.x;
	int i=0;
	// Considering only "gas-vacuum" configuration, i.e. gas fills the left side of the simulation segment
	const double eps = .01;
	double uvac = 0.;
	// Calculating the n+1-time-layer x coordinate of gas-vacuum boundary
	for(i=imin; i<=imax; i++) {
		// Is it interface cell?
		Vector4 flux = Vector4::ZERO;
		if(U[i-1][0] == 0.) {
			if(U[i][0] != 0.) {
				double gamma = eos.getc(1., 1)*eos.getc(1., 1),
					   _ro = U[i][0],
					   _u = U[i][1]/_ro, 
					   _e = U[i][2]/_ro - .5*_u*_u,
					   _p = eos.getp(_ro, _e),
					   _c = eos.getc(_ro, _p);
				uvac = _u - 2./(gamma - 1.)*_c;
				xbnd += uvac*dt;
				break;
			}
		}
	}
	for(i=imin; i<=imax; i++) {
		
		
		
		
		
		
		if (i==52) {

			double qq = 0.;
		}
		
		
		
		
		
		
		
		
		
		
		// Is it interface cell?
		Vector4 flux = Vector4::ZERO;
		if(U[i-1][0] == 0.) {
			if(U[i][0] != 0.) {
				double _ro = U[i][0],
						_u = U[i][1]/_ro, 
						_e = U[i][2]/_ro - .5*_u*_u,
						_p = eos.getp(_ro, _e),
						_c = eos.getc(_ro, _p);
				if(x[i]-xbnd > eps*dx) {
					if (fabs(_u) > _c) {
						/*// поток от начального условия F(U[0](t))
						double gamma = eos.getc(1., 1)*eos.getc(1., 1);
						double _ro = U[i][0];
						double _u = U[i][1]/_ro;
						double _e = U[i][2]/_ro-.5*_u*_u;
						double _p = eos.getp(_ro, _e);
						double _c = eos.getc(_ro, _p);
						double _x = (x[i] + x[i+1])/2.;
						double _u_res = ((gamma-1.)*pr.ur + 2.*(_x/t + _c)) / (gamma+1.);
						double _ro_res = pow((_x/t - _u_res)*(_x/t-_u_res)*pow(pr.ror, gamma)/gamma/pr.pr, 1./(gamma-1.));
						double _p_res = pow(_ro_res, gamma)/pr.ror*pr.pr;
						// flux = calcPhysicalFlux(eos, pr.ror, pr.ur, pr.pr);
						flux = calcPhysicalFlux(eos, _ro_res, _u_res, _p_res);*/
					
						
						
						
					//	double ro0 = pr.ror, u0 = pr.ur, E0 = eos.gete(pr.ror, pr.pr) + .5*u0*u0;
					//	flux = rslv.calcFlux(eos, 0., 0., 0., ro0, ro0*u0, ro0*E0);

						// double q = xbnd-x[i];
						// double qq = dt;
						// double qqq =q/qq;
						
						//flux = rslv.calcFlux(eos, 0., 0., 0., U[i][0], U[i][1], U[i][2]); 

						double _ro = U[i][0];
						double _u = U[i][1]/_ro;
						double _e = U[i][2]/_ro-.5*_u*_u;
						double _p = eos.getp(_ro, _e);
						flux = calcPhysicalFlux(eos, _ro, _u, _p);						
						//C1DVectorPrimitive res = calcSolution(eos, 0., 0., 0., _ro, _u, _p, xbnd-pr.x0, dt);
						//flux = calcPhysicalFlux(eos, res.ro, res.u, res.p);						
						


					} else {
					    // поток F(U[i])
						//double _ro = U[i][0], 
						//	    _u = U[i][1]/_ro,
						//		_e = U[i][2]/_ro - .5*_u*_u,
						//		_p = eos.getp(_ro, _e);
						// flux = calcPhysicalFlux(eos, _ro, _u, _p);
						// flux = rslv.calcFlux(eos, 0., 0., 0., U[i][0], U[i][1], U[i][2]); 
						double ro0 = pr.ror, u0 = pr.ur, p0 = pr.pr, E0 = eos.gete(pr.ror, pr.pr) + .5*u0*u0;
						double c0 = eos.getc(ro0, p0);
						//Vector4 _flux = rslv.calcFlux(eos, 0., 0., 0., ro0, ro0*u0, ro0*E0);
						double gamma = eos.getc(1., 1)*eos.getc(1., 1);
						double _u = ((gamma-1.)*u0 - 2*c0)/(gamma+1.);
						double _ro = pow(_u*_u*pow(ro0, gamma)/(gamma/p0), 1./(gamma-1.));
						double _p = pow(_ro/ro0, gamma)*p0;
						flux = calcPhysicalFlux(eos, _ro, _u, _p);
						// flux = calcPhysicalFlux(eos, ro0, u0, p0);
						double qq = 0.;
				    }
			    } else {
					flux = Vector4::ZERO;
				}

		    } else {
		        flux = Vector4::ZERO;
		    }
		} else {
			flux = rslv.calcFlux(eos, U[i-1][0], U[i-1][1], U[i-1][2], U[i][0], U[i][1], U[i][2]); 
		}
		F[i] = flux;
	}
	for(i=imin; i<imax; i++) {
		for(int counter=0; counter<3; counter++) newU[i][counter] = U[i][counter] - dt/dx*(F[i+1][counter] - F[i][counter]);
	}
	for(i=imin; i<imax; i++) {
		for(int counter=0; counter<3; counter++) U[i][counter] = newU[i][counter];
	}	
	pr.setbcs(fld.U);
}


Vector4 calcPhysicalFlux(FEOS& eos, double rho, double u, double p) {
	if (rho == 0.) return Vector4::ZERO;
	double e = eos.gete(rho, p);
	return Vector4(rho*u, p + rho*u*u, u*(p + rho*(e + .5*u*u)), 0.);
}

// Godunov-exact solver flux
Vector4 CExactRiemannSolver::calcFlux(FEOS& eos, double roL, double rouL, double roEL, double roR, double rouR, double roER) {	
	double uL = 0., uR = 0., pL = 0., pR = 0.;
	if(roL!=0.) {
	    uL = rouL/roL;
		pL = eos.getp(roL, roEL/roL - .5*uL*uL);
	}
	if(roR!=0.) {
		uR = rouR/roR;
		pR = eos.getp(roR, roER/roR - .5*uR*uR);
	}
	C1DVectorPrimitive res = calcSolution(eos, roL, uL, pL, roR, uR, pR, 0., .01);
	double E = eos.gete(res.ro, res.p) + .5*res.u*res.u;
	Vector4 FGodunov = Vector4(res.ro*res.u, res.ro*res.u*res.u + res.p, res.u*(res.ro*E+res.p), 0.);
	return FGodunov;
}

// HLL flux for general EOS
Vector4 CHLLRiemannSolver::calcFlux(FEOS& eos, double roL, double rouL, double roEL, double roR, double rouR, double roER) {
	double _ro = 0., _u = 0., _e = 0., _p = 0.;
	double uL = rouL/roL, uR = rouR/roR;	
	double EL = roEL/roL, ER = roER/roR;
	double eL = roEL/roL - .5*uL*uL, eR = roER/roR - .5*uR*uR; 
	if(roL == 0. || roR == 0.) {
		cout << "Error: CSolver::calcHLLFluxEOSIdeal(): vacuum is present." << endl;
		exit(1);	}
	double pL = eos.getp(roL, eL), pR = eos.getp(roR, eR);
	double HL = EL + pL/roL, HR = ER + pR/roR;
	double cL = eos.getc(roL, pL), cR = eos.getc(roR, pR);	
	Vector4 UL = Vector4(roL, rouL, roEL, 0.), UR = Vector4(roR, rouR, roER, 0.);
	_ro = roL, _u = rouL/roL, _e = roEL/roL-.5*_u*_u, _p=eos.getp(_ro, _e);
	Vector4 FL = Vector4(rouL, _p + _ro*_u*_u, _u*(_p + roEL), 0.);
	_ro = roR, _u = rouR/roR, _e = roER/roR-.5*_u*_u, _p=eos.getp(_ro, _e);
	Vector4 FR = Vector4(rouR, _p + _ro*_u*_u, _u*(_p + roER), 0.); 

	// Roe averaging
	double roAv = sqrt(roL*roR);
	double uAv = (uL*sqrt(roL) + uR*sqrt(roR))/(sqrt(roL)+sqrt(roR));
	double HAv = (HL*sqrt(roL) + HR*sqrt(roR))/(sqrt(roL)+sqrt(roR));
	double EAv = (EL*sqrt(roL) + ER*sqrt(roR))/(sqrt(roL)+sqrt(roR));
	double pAv = (HAv - EAv) * roAv;
	double cAv = eos.getc(roAv, pAv);
	double SL = uAv-cAv, SR = uAv+cAv;

	// Vector4 Uhll = (SR*UR - SL*UL + FL - FR)/(SR - SL);
	Vector4 Fhll = Vector4::ZERO;
    if (0. <= SL)  
		Fhll = FL;
	else if (SL <= 0. && 0 <= SR)
		Fhll = (SR*FL - SL*FR + SL*SR*(UR - UL))/(SR - SL);
	else 
		Fhll = FR;
	return Fhll;
}

// HLLC flux for general EOS
Vector4 CHLLCRiemannSolver::calcFlux(FEOS& eos, double roL, double rouL, double roEL, double roR, double rouR, double roER) {
	double _ro = 0., _u = 0., _e = 0., _p = 0.;
	double uL = rouL/roL, uR = rouR/roR;	
	double EL = roEL/roL, ER = roER/roR;
	double eL = roEL/roL - .5*uL*uL, eR = roER/roR - .5*uR*uR; 
	if(roL == 0. || roR == 0.) {
		cout << "Error: CSolver::calcHLLFluxEOSIdeal(): vacuum is present." << endl;
		exit(1);	}
	double pL = eos.getp(roL, eL), pR = eos.getp(roR, eR);
	double HL = EL + pL/roL, HR = ER + pR/roR;
	double cL = eos.getc(roL, pL), cR = eos.getc(roR, pR);	
	Vector4 UL = Vector4(roL, rouL, roEL, 0.), UR = Vector4(roR, rouR, roER, 0.);
	_ro = roL, _u = rouL/roL, _e = roEL/roL-.5*_u*_u, _p=eos.getp(_ro, _e);
	Vector4 FL = Vector4(rouL, _p + _ro*_u*_u, _u*(_p + roEL), 0.);
	_ro = roR, _u = rouR/roR, _e = roER/roR-.5*_u*_u, _p=eos.getp(_ro, _e);
	Vector4 FR = Vector4(rouR, _p + _ro*_u*_u, _u*(_p + roER), 0.); 
	// Step 1 -- pressure estimate from primitive-variable Riemann solver (PVRS)
	double roAv = .5*(roL + roR), cAv = .5*(cL + cR), pPVRS = .5*(pL + pR) - .5*(uR - uL)*roAv*cAv;
	//double pStar = max(0., pPVRS);
	double pStar = max(-15.e9, pPVRS);
	
	// Step 2 -- wave speed estimates	
	
	// Uncomment for Roe averaging. Works quite well for Mie-Gruneisen EOS but fails on Toro test 1 at rarefaction wave (imposes the discontinuity)
	double _roAv = sqrt(roL*roR);
	double _uAv = (uL*sqrt(roL) + uR*sqrt(roR))/(sqrt(roL)+sqrt(roR));
	double _HAv = (HL*sqrt(roL) + HR*sqrt(roR))/(sqrt(roL)+sqrt(roR));
	double _EAv = (EL*sqrt(roL) + ER*sqrt(roR))/(sqrt(roL)+sqrt(roR));
	double _pAv = (_HAv - _EAv) * _roAv;
	double _cAv = eos.getc(_roAv, _pAv);

	double SL = min(uL - cL, _uAv - _cAv), SR = max(uR + cR, _uAv + _cAv);

	// Uncomment for Einfeldt estimates (is it possibly just HLLE solver?)
	/*
	double uAv = (uL*sqrt(roL) + uR*sqrt(roR))/(sqrt(roL)+sqrt(roR));
	double eta2 = .5*sqrt(roL*roR)/(sqrt(roL) + sqrt(roR))/(sqrt(roL) + sqrt(roR));
	double dSqAv = (cL*cL*sqrt(roL) + cR*cR*sqrt(roR))/(sqrt(roL)+sqrt(roR)) + eta2*(uR-uL)*(uR-uL);
	double SL = min(uL - cL, uAv - sqrt(dSqAv)), SR = max(uR + cR, uAv + sqrt(dSqAv));
	*/

	// Uncomment for naive estimate from Toro. Works pretty on Toro tests but thigs go not so well with Mie-Gruneien EOS
	// double SL = min(min(uL - cL, uR - cR), 0.), SR = max(max(uL + cL, uR + cR), 0.);
	
	
/*  double _roAv = sqrt(roL*roR);
	double sqroL = sqrt(roL), sqroR = sqrt(roR);
	double _uAv = (sqroL*uL + sqroR*uR)/(sqroL+sqroR), 
		   _HAv = (sqroL*HL + sqroR*HR)/(sqroL+sqroR),
		   _eAv = (sqroL*eL + sqroR*eR)/(sqroL+sqroR);	
	double _pAv = (_HAv - _eAv - .5*_uAv*_uAv) * _roAv;
	double dpdrhoAv = 0., dpdeAv = 0.;
	if(eL != eR)
		dpdeAv = (.5*(pR + eos.getp(roL, eR)) - .5*(eos.getp(roR, eL)+pL))/(eR-eL);
	else 
		dpdeAv = .5*(eos.getdpde(roL, eL)+eos.getdpde(roR, eR));
	if(roL!=roR)
		dpdrhoAv = (.5*(pR + eos.getp(roR, eL)) - .5*(eos.getp(roL, eR)+pL))/(roR-roL);
	else
		dpdrhoAv = .5*(eos.getdpdrho(roL, eL)+eos.getdpdrho(roR, eR));
    double _cAv = sqrt(_pAv*dpdeAv/_roAv/_roAv + dpdrhoAv); */



	// double SL = _uAv-_cAv, SR = _uAv+_cAv;








	//double SL = min(min(uL - cL, uR - cR), 0.), SR = max(max(uL + cL, uR + cR), 0.);

	double SStar = (pR - pL + roL*uL*(SL - uL) - roR*uR*(SR - uR))/(roL*(SL - uL) - roR*(SR - uR));
	// Vector4 Uhll = (SR*UR - SL*UL + FL - FR)/(SR - SL);
	Vector4 Fhllc = Vector4::ZERO;
	Vector4 D = Vector4(0., 1., SStar, 0.);
    Vector4 UStarL = (SL*UL - FL + pStar*D)/(SL - SStar), UStarR = (SR*UR - FR + pStar*D)/(SR - SStar);
	Vector4 FStarL = FL + SL*(UStarL - UL), FStarR = FR + SR*(UStarR - UR);	
	if (0. <= SL)  
		Fhllc = FL;
	else if (SL <= 0. && 0 <= SStar)
		Fhllc = FStarL;
	else if (SStar <= 0. && 0 <= SR)
		Fhllc = FStarR;
	else 
		Fhllc = FR;
	return Fhllc;
}


Vector4 CGPSRiemannSolver::calcFlux(FEOS& eos, double roL, double rouL, double roEL, double roR, double rouR, double roER) {
	double _ro=0., _u=0., _p=0., _e=0., _E=0., _sigma = 0., sigmaL = 0., sigmaR = 0.;
	double uL = rouL/roL, uR = rouR/roR, EL = roEL/roL, ER = roER/roR, eL = EL - 0.5*uL*uL, eR = ER - 0.5*uR*uR, 
		   pL = eos.getp(roL, eL), pR = eos.getp(roR, eR), cL = eos.getc(roL, pL), cR = eos.getc(roR, pR);
	double gamma = cL*cL/(pL/roL);
	sigmaL = pL/pow(roL, gamma);
	sigmaR = pR/pow(roR, gamma);
	if(uL > cL) 
		return Vector4(roL*uL, pL+roL*uL*uL, uL*(pL+roL*EL), 0.);
	else if (uR < -cR)
		return Vector4(roR*uR, pR+roR*uR*uR, uR*(pR+roR*ER), 0.);
	_p = (pL / roL / cL + pR / roR / cR + uL - uR) / (1. / roL / cL + 1 / roR / cR);
	_u = (roL * cL * uL + roR * cR * uR + pL - pR) / (roL * cL + roR * cR);
	if(_u > 0.)
		_sigma = sigmaL; 
	else 
		_sigma = sigmaR;
	_ro = pow(_p/_sigma, 1./gamma);
	_e  = eos.gete(_ro, _p);
	_E  = _e + 0.5*_u*_u;
	return Vector4(_ro*_u, _p+_ro*_u*_u, _u*(_p+_ro*_E), 0.);
}


Vector4 CLFRiemannSolver::calcFlux(
		FEOS& eos,
		double roL, double rouL, double roEL,
		double roR, double rouR, double roER,
		double dx, double dt) {
	Vector4 UL = Vector4(roL, rouL, roEL, 0.);
	Vector4 UR = Vector4(roR, rouR, roER, 0.);

	double uL = rouL / roL;
	double eL = roEL/roL - .5*uL*uL;
	double uR = rouR / roR;
	double eR = roER/roR - .5*uR*uR;
	Vector4 FL = calcPhysicalFlux(eos, roL, rouL / roL, eos.getp(roL, eL));
	Vector4 FR = calcPhysicalFlux(eos, roR, rouR / roR, eos.getp(roR, eR));

	Vector4 F = 0.5 * (FL + FR) - 0.5 * dx / dt * (UR - UL);

	return F;
}


//double calcSquareSoundSpeed(
//		double rho, double rho_v, double rho_E, double gamma = 1.4) {
//	/* Compute the square of sound speed. */

//	return std::abs(gamma * (gamma - 1.)
//					* (rho_E - rho_v * rho_v * 0.5 / rho) / rho);
//}


double calcMaxWaveSpeedD(
		const std::ranges::common_range auto& u_arr,
		FEOS& eos) {
	/* Calculate max (sound speed + abs(v)) ~ |df/du|
	 * for 1D Euler eq'ns. */

//	return std::ranges::max(std::ranges::transform_view(
//		std::as_const(u_arr),
//		[gamma](const auto& u_arr_vec_pt) -> double {
//			return (std::sqrt(std::abs(calcSquareSoundSpeed(
//								u_arr_vec_pt[0],
//								u_arr_vec_pt[1],
//								u_arr_vec_pt[2], gamma)))
//					+ std::abs(u_arr_vec_pt[1] / u_arr_vec_pt[0]));
//		}
//	));
	auto&& U = u_arr;
	double rho = U[0][0];
	double u = U[0][1]/rho;
	double e = U[0][2]/rho - .5*u*u;
	double p = eos.getp(rho, e);
	double c = eos.getc(rho, p);
	double max_lambda = std::max(
				std::abs(u), std::max(std::abs(u-c), std::abs(u+c)));
	double umax = max_lambda;
	for (std::size_t i = 0; i < std::ranges::size(U); ++ i) {
		rho = U[i][0];
		u = U[i][1]/rho;
		e = U[i][2]/rho - .5*u*u;
		p = eos.getp(rho, e);
		c = eos.getc(rho, p);
		umax = std::abs(u) + std::abs(c);
		if (max_lambda < umax) max_lambda = umax;
	}

	return max_lambda;
}


Vector4 CLFGlobalRiemannSolver::calcFlux(
		FEOS& eos,
		double roL, double rouL, double roEL,
		double roR, double rouR, double roER,
		double lambda) {
	Vector4 UL = Vector4(roL, rouL, roEL, 0.);
	Vector4 UR = Vector4(roR, rouR, roER, 0.);

	double uL = rouL / roL;
	double eL = roEL/roL - .5*uL*uL;
	double uR = rouR / roR;
	double eR = roER/roR - .5*uR*uR;
	Vector4 FL = calcPhysicalFlux(eos, roL, rouL / roL, eos.getp(roL, eL));
	Vector4 FR = calcPhysicalFlux(eos, roR, rouR / roR, eos.getp(roR, eR));

	double alpha = /* 1.1 * */ lambda;

//	Vector4 F(0., 0., 0., 0.);

//	for (std::size_t cmpnt = 0; cmpnt < 3; ++ cmpnt)
//		 F[cmpnt] = 0.5 * (FL[cmpnt] + FR[cmpnt])
//				 - 0.5 * alpha * (UR[cmpnt] - UL[cmpnt]);

	Vector4 F = 0.5 * (FL + FR) - 0.5 * alpha * (UR - UL);

	return F;
}


void C1D2ndOrderMethod::calcFluxField(
		C1DProblem& pr, FEOS& eos, C1DField& fld) {
//	double roL = 0., uL = 0., eL = 0., pL = 0.,
//			roR = 0., uR = 0., eR = 0., pR = 0., E = 0.;
	int i=0, imin = fld.imin, imax = fld.imax;
	auto&& F = fld.F;
	auto&& ULx = rec.ULx;
	auto&& URx = rec.URx;
	// TODO: проверить на скорость выполнения операций,
	// сравнить с реализацией через тип Vector4 --
	// если не медленнее, то в дальнейшем избавиться от Vector4 везде
	rec.calc(fld);
	for (i = imin; i <= imax; ++ i) {
		F[i] = rslv.calcFlux(
					eos,
					URx[i-1][0], URx[i-1][1], URx[i-1][2],
					ULx[i][0], ULx[i][1], ULx[i][2]);
	}
}


void C1D2ndOrderMethod::calc(
		C1DProblem& pr, FEOS& eos, C1DField& fld) {
	double dx = fld.dx;
	double t = fld.t;
	double dt = fld.dt;

	int i = 0;
	int imin = fld.imin;
	int imax = fld.imax;

	auto&& ULx = rec.ULx;
	auto&& URx = rec.URx;
	pr.setbcs(ULx);
	pr.setbcs(URx);

	auto&& F = fld.F;
	auto&& U = fld.U;
	auto&& newU = fld.newU;

	calcFluxField(pr, eos, fld);

	for (i = imin; i < imax; ++ i) {
//		for (int counter = 0; counter < 3; ++ counter)
//			newU[i][counter] = U[i][counter] - dt / dx * (
//						F[i+1][counter] - F[i][counter]);
		newU[i] = U[i] - dt / dx * (F[i+1] - F[i]);
	}

	for (i = imin; i < imax; ++ i) {
//		for (int counter = 0; counter < 3; ++ counter)
//			U[i][counter] = newU[i][counter];
		U[i] = newU[i];
	}	
	pr.setbcs(fld.U);
}


void C1D2ndOrderLFGlobalMethod::calcFluxField(
		C1DProblem& pr, FEOS& eos, C1DField& fld) {
//	double roL = 0., uL = 0., eL = 0., pL = 0.,
//			roR = 0., uR = 0., eR = 0., pR = 0., E = 0.;
	int i=0, imin = fld.imin, imax = fld.imax;
	auto&& F = fld.F;
	auto&& ULx = rec.ULx;
	auto&& URx = rec.URx;
	// TODO: проверить на скорость выполнения операций,
	// сравнить с реализацией через тип Vector4 --
	// если не медленнее, то в дальнейшем избавиться от Vector4 везде
	rec.calc(fld);
	double lambda = calcMaxWaveSpeedD(fld.U, eos);
	for (i = imin - 1; i <= imax; ++ i) {
		F[i] = lfrslv.calcFlux(
					eos,
					URx[i-1][0], URx[i-1][1], URx[i-1][2],
					ULx[i][0], ULx[i][1], ULx[i][2], lambda);
	}
}


void C1DBGKMethod::calc(C1DProblem& pr, FEOS& eos, C1DField& fld) {
	double roL = 0., uL = 0., eL = 0., pL = 0., roR = 0., uR = 0., eR = 0., pR = 0., E = 0.;
	double dx = fld.dx, t = fld.t, dt = fld.dt;
	int i=0, imin = fld.imin, imax = fld.imax;
	auto&& U = fld.U;
	auto&& newU = fld.newU;
	auto&& F = fld.F;
	for (i = imin - 1; i < imax; ++ i) {







		if (i == 31)
		{ 
			
			
			double q = 1.;
		
		
		}




		Vector4 _Um = Vector4(U[i-1][0],U[i-1][1],U[i-1][2],0.);
		Vector4 _U = Vector4(U[i][0],U[i][1],U[i][2],0.);
		Vector4 _Up = Vector4(U[i+1][0],U[i+1][1],U[i+1][2],0.);
		Vector4 _Upp = Vector4(U[i+2][0],U[i+2][1],U[i+2][2],0.);
		F[i] = bgk.calcFlux(eos, _Um, _U, _Up, _Upp, dx, dt);
	}
	for(i=imin; i<imax; i++) {
		for(int counter=0; counter<3; counter++) newU[i][counter] = U[i][counter] - dt/dx*(F[i][counter] - F[i-1][counter]);
	}
	for(i=imin; i<imax; i++) {
		for(int counter=0; counter<3; counter++) U[i][counter] = newU[i][counter];
	}	
	pr.setbcs(fld.U);
}

void C1DLFMethod::calc(C1DProblem& pr, FEOS& eos, C1DField& fld) {
	double roL = 0., uL = 0., eL = 0., pL = 0.,
		   roR = 0., uR = 0., eR = 0., pR = 0., 
		   E = 0.;
	double dx = fld.dx, t = fld.t, dt = fld.dt;
	int imin = fld.imin, imax = fld.imax;
	auto&& U = fld.U;
	auto&& newU = fld.newU;
	auto&& F = fld.F;
	// TODO: проверить на скорость выполнения операций,
	// сравнить с реализацией через тип Vector4 -- если не медленнее,
	// то в дальнейшем избавиться от Vector4 везде
	int i=0;	
	// Потоки считаем по алгоритму решения задаче о распаде разрыва для УРС
	// Ми-Грюнайзена из работы [Miller, Puckett]
	for (i = imin; i <= imax; i++) {
		F[i] = lfrslv.calcFlux(
					eos,
					U[i-1][0], U[i-1][1], U[i-1][2],
					U[i][0], U[i][1], U[i][2],
					dx, dt);
		//n.W_temp = n.W - dt/h*(Fp-Fm);
	}
	for (i = imin; i < imax; ++ i) {
		for (int counter = 0; counter < 3; ++ counter)
			newU[i][counter]
					= U[i][counter] - dt / dx * (
						F[i+1][counter] - F[i][counter]);
	}
	for (i = imin; i < imax; ++ i) {
		for (int counter = 0; counter < 3; ++ counter)
			U[i][counter] = newU[i][counter];
	}
	pr.setbcs(fld.U);
}


void C1DLFGlobalMethod::calc(C1DProblem& pr, FEOS& eos, C1DField& fld) {
	double roL = 0., uL = 0., eL = 0., pL = 0.,
		   roR = 0., uR = 0., eR = 0., pR = 0.,
		   E = 0.;
	double dx = fld.dx, t = fld.t, dt = fld.dt;
	int imin = fld.imin, imax = fld.imax;
	auto&& U = fld.U;
	auto&& newU = fld.newU;
	auto&& F = fld.F;
	// TODO: проверить на скорость выполнения операций,
	// сравнить с реализацией через тип Vector4 -- если не медленнее,
	// то в дальнейшем избавиться от Vector4 везде
	int i=0;
	// Потоки считаем по алгоритму решения задаче о распаде разрыва для УРС
	// Ми-Грюнайзена из работы [Miller, Puckett]
	double lambda = calcMaxWaveSpeedD(fld.U, eos);
	for (i = imin; i <= imax; i++) {
		F[i] = lfrslv.calcFlux(
					eos,
					U[i-1][0], U[i-1][1], U[i-1][2],
					U[i][0], U[i][1], U[i][2],
					lambda);
		//n.W_temp = n.W - dt/h*(Fp-Fm);
	}
	for (i = imin; i < imax; ++ i) {
		for (int counter = 0; counter < 3; ++ counter)
			newU[i][counter]
					= U[i][counter] - dt / dx * (
						F[i+1][counter] - F[i][counter]);
	}
	for (i = imin; i < imax; ++ i) {
		for (int counter = 0; counter < 3; ++ counter)
			U[i][counter] = newU[i][counter];
	}
	pr.setbcs(fld.U);
}


Vector4 CRoeRiemannSolver::calcFlux(FEOS& eos, double roL, double rouL, double roEL, double roR, double rouR, double roER) {
	double uL = rouL/roL, eL = roEL/roL - 0.5*uL*uL, uR = rouR/roR, eR = roER/roR - 0.5*uR*uR;
	double pL = eos.getp(roL,eL), pR = eos.getp(roR,eR), cL = eos.getc(roL,pL), cR = eos.getc(roR,pR);
	double EL = roEL/roL, ER = roER/roR, HL = (roEL + pL)/roL, HR = (roER + pR)/roR;
	double roAv = sqrt(roL*roR);
	double sqroL = sqrt(roL), sqroR = sqrt(roR);
	double uAv = (sqroL*uL + sqroR*uR)/(sqroL+sqroR), 
		   HAv = (sqroL*HL + sqroR*HR)/(sqroL+sqroR),
		   EAv = (sqroL*EL + sqroR*ER)/(sqroL+sqroR);	
	double pAv = (HAv - EAv) * roAv;
	double cAv = eos.getc(roAv,pAv);			
	// Ideal version to compare:	
	/*double _gam = 5./3.;
	cAv = sqrt((HAv-.5*uAv*uAv)*(_gam-1.));*/
	Vector4 Lambda = Vector4(uAv - cAv, uAv, uAv + cAv, 0.),
		        K0 = Vector4(1., uAv-cAv, HAv-uAv*cAv, 0.),
			    K1 = Vector4(1., uAv,     uAv*uAv/2.,  0.),
				K2 = Vector4(1., uAv+cAv, HAv+uAv*cAv, 0.);
	Vector4 Alpha = Vector4::ZERO;
	// For ideal gas: 
	// Alpha[1] = (gamma-1.)/cAv/cAv * ((HAv-uAv*uAv)*(roR-roL) + uAv*(rouR-rouL) - (roER-roEL));
	// Using cAv2 = (gamma-1) * (HAv-.5*uAv2), (gamma-1)/cAv2 = 1/(HAv-.5uAv2) 
	Alpha[1] = 1./(HAv-.5*uAv*uAv) * ((HAv-uAv*uAv)*(roR-roL) + uAv*(rouR-rouL) - (roER-roEL));
	Alpha[0] = 1./2./cAv * ((roR-roL)*(uAv+cAv) - (rouR-rouL) - cAv*Alpha[1]);
	Alpha[2] = (roR-roL) - Alpha[0] - Alpha[1];
	Vector4 FL = Vector4(roL*uL, pL+roL*uL*uL, uL*(pL+roEL), 0.);
	Vector4 FR = Vector4(roR*uR, pR+roR*uR*uR, uR*(pR+roER), 0.);
	double delta=0.;
	delta = max(0., 4.*((uR-cR)-(uL-cL)) );
    if (fabs(Lambda[0]) < 0.5*delta) Lambda[0] = Lambda[0]*Lambda[0]/delta + .25*delta;
    delta = max(0., 4.*((uR+cR)-(uL+cL)) );
    if (fabs(Lambda[2]) < 0.5*delta) Lambda[2] = Lambda[2]*Lambda[2]/delta + .25*delta;
	Vector4 FRoe = 0.5*(FL + FR) - 0.5*(Alpha[0]*fabs(Lambda[0])*K0 + 
		                                Alpha[1]*fabs(Lambda[1])*K1 + 
										Alpha[2]*fabs(Lambda[2])*K2);
	return FRoe;
}


// According to [T.Glaister, JoCP, 72, 382-408 (1988)]
Vector4 CRoeGeneralRiemannSolver::calcFlux(FEOS& eos, double roL, double rouL, double roEL, double roR, double rouR, double roER) {
	double uL = rouL/roL, eL = roEL/roL - 0.5*uL*uL, uR = rouR/roR, eR = roER/roR - 0.5*uR*uR;
	double pL = eos.getp(roL,eL), pR = eos.getp(roR,eR), cL = eos.getc(roL,pL), cR = eos.getc(roR,pR);
	double EL = roEL/roL, ER = roER/roR, HL = (roEL + pL)/roL, HR = (roER + pR)/roR;
	double roAv = sqrt(roL*roR);
	double sqroL = sqrt(roL), sqroR = sqrt(roR);
	double uAv = (sqroL*uL + sqroR*uR)/(sqroL+sqroR), 
		   HAv = (sqroL*HL + sqroR*HR)/(sqroL+sqroR),
		   eAv = (sqroL*eL + sqroR*eR)/(sqroL+sqroR);	
	double pAv = (HAv - eAv - .5*uAv*uAv) * roAv;
	double dpdrhoAv = 0., dpdeAv = 0.;
	if(eL != eR)
		dpdeAv = (.5*(pR + eos.getp(roL, eR)) - .5*(eos.getp(roR, eL)+pL))/(eR-eL);
	else 
		dpdeAv = .5*(eos.getdpde(roL, eL)+eos.getdpde(roR, eR));
	if(roL!=roR)
		dpdrhoAv = (.5*(pR + eos.getp(roR, eL)) - .5*(eos.getp(roL, eR)+pL))/(roR-roL);
	else
		dpdrhoAv = .5*(eos.getdpdrho(roL, eL)+eos.getdpdrho(roR, eR));
    double cAv = sqrt(pAv*dpdeAv/roAv/roAv + dpdrhoAv);
	
	
	
	
//	double _c = eos.getc(roAv, pAv);
	
	


	// Ideal version to compare:	
	/*double _gam = 5./3.;
	cAv = sqrt((HAv-.5*uAv*uAv)*(_gam-1.));*/
	Vector4 Lambda = Vector4(uAv - cAv, uAv, uAv + cAv, 0.),
		        K0 = Vector4(1., uAv-cAv, HAv-uAv*cAv, 0.),
			    K1 = Vector4(1., uAv,     uAv*uAv/2.,  0.),
				K2 = Vector4(1., uAv+cAv, HAv+uAv*cAv, 0.);
	Vector4 Alpha = Vector4::ZERO;
	// For ideal gas: 
	// Alpha[1] = (gamma-1.)/cAv/cAv * ((HAv-uAv*uAv)*(roR-roL) + uAv*(rouR-rouL) - (roER-roEL));
	// Using cAv2 = (gamma-1) * (HAv-.5*uAv2), (gamma-1)/cAv2 = 1/(HAv-.5uAv2) 
	Alpha[1] = 1./(HAv-.5*uAv*uAv) * ((HAv-uAv*uAv)*(roR-roL) + uAv*(rouR-rouL) - (roER-roEL));
	Alpha[0] = 1./2./cAv * ((roR-roL)*(uAv+cAv) - (rouR-rouL) - cAv*Alpha[1]);
	Alpha[2] = (roR-roL) - Alpha[0] - Alpha[1];
	Vector4 FL = Vector4(roL*uL, pL+roL*uL*uL, uL*(pL+roEL), 0.);
	Vector4 FR = Vector4(roR*uR, pR+roR*uR*uR, uR*(pR+roER), 0.);
	double delta=0.;
	delta = max(0., 4.*((uR-cR)-(uL-cL)) );
    if (fabs(Lambda[0]) < 0.5*delta) Lambda[0] = Lambda[0]*Lambda[0]/delta + .25*delta;
    delta = max(0., 4.*((uR+cR)-(uL+cL)) );
    if (fabs(Lambda[2]) < 0.5*delta) Lambda[2] = Lambda[2]*Lambda[2]/delta + .25*delta;
	Vector4 FRoe = 0.5*(FL + FR) - 0.5*(Alpha[0]*fabs(Lambda[0])*K0 + 
		                                Alpha[1]*fabs(Lambda[1])*K1 + 
										Alpha[2]*fabs(Lambda[2])*K2);
	return FRoe;
}


Matrix4 CBGKRiemannSolver::getOmega(FEOS& eos, Vector4 parameter) {
	const double gamma = eos.getc(1.,1.)*eos.getc(1.,1.);
	//double _rho = U[0], _u = U[1]/_rho, _E = U[2]/_rho, _e = _E-_u*_u/2., _p=eos.getp(_rho,_e), _c = eos.getc(_rho,_p);
	double _rho = parameter[0], _u = parameter[1], _H = parameter[2], _c = parameter[3];	
	
	/*double _H
	double _c = sqrt((gamma-1.)*(_H-.5*_u*_u));*/
	
	
	
	
	
	double dpdroE = gamma-1.;
/*	double _Eprho = _E + _p/_rho; 
	return Matrix4(
		1.,			  1.,		           1.,			   0.,
		_u-_c,		  _u,				   _u+_c,		   0.,
		_Eprho-_u*_c, _Eprho-_c*_c/dpdroE, _Eprho+_u*_c,   0.,
		0.,			  0.,				   0.0,			   1. ); */

	
	return Matrix4(
		1.,		  1.,	    1.,		  0.,
		_u-_c,	  _u,		_u+_c,	  0.,
		_H-_u*_c, .5*_u*_u, _H+_u*_c, 0.,
		0.,		  0.,		0.,		  1. );
}

Matrix4 CBGKRiemannSolver::getOmegaInv(FEOS& eos, Vector4 U) {
	double _rho = U[0], _u = U[1]/_rho, _E = U[2]/_rho, _e = _E-_u*_u/2., _p=eos.getp(_rho,_e), _c = eos.getc(_rho,_p);
	
/* // Uncomment for ideal gas values
	const double gamma = eos.getc(1.,1.)*eos.getc(1.,1.);
	double dpdro  = .5*(gamma-1.)*_u*_u;
	double dpdrov = -(gamma-1.)*_u;
	double dpdroE = gamma-1.;
	*/

	double _prho = eos.getdpdrho(_rho,_e);
	double _pe = eos.getdpde(_rho,_e);
	double dpdro = _prho-(_e-.5*_u*_u)/_rho*_pe;
	double dpdrov = -_u/_rho*_pe;
	double dpdroE = 1./_rho*_pe;

	double cc = 1./2./_c/_c;
	return Matrix4(
			cc*(dpdro + _u*_c),	   cc*(dpdrov - _c), cc*(dpdroE),	   0.0,
			cc*(2.*(_c*_c-dpdro)), cc*(-2.*dpdrov),  cc*(-2.*dpdroE),  0.0,
			cc*(dpdro - _u*_c),	   cc*(dpdrov +_c),  cc*(dpdroE),	   0.0,
			0.0,				   0.0,				 0.0,		       1.0
		);
}

void CBGKRiemannSolver::fillLambda(Vector4 Fm, Vector4 F, Vector4 Fp, Vector4 Fpp, Vector4 L, double step) {
	double curant = 0.;
	double criteria = 0.;  
	for(int i=0; i<4; i++) 	{
		curant   = fabs(L[i]) * step;
		criteria = L[i] * (fabs(Fpp[i]-Fp[i]) - fabs(F[i]-Fm[i]));
		fillLambdaComponent(i, L[i], criteria, curant);
	}
}

void CBGKRiemannSolver::fillLambdaComponent(int i, double lambda, double criteria, double curant) {
	double alpha, beta, gamma, delta;
	if( criteria >= 0 && lambda >= 0 ) {
		alpha = -0.5 * (1.0-curant);
		 beta =  1.0 - alpha;
		gamma =  0.0; 
		delta =  0.0;
	}
	else if( criteria >= 0 && lambda < 0 ) {
		alpha =  0.0;
		 beta =  0.0;
		delta = -0.5 * (1.0-curant);
		gamma =  1.0 - delta; 
	}
	else if( criteria < 0 && lambda >= 0 ) {
		alpha =  0.0;
		gamma =  0.5 * (1.0-curant);
		 beta =  1.0 - gamma;
		delta =  0.0;
	}
	else {
		alpha =  0.0;
		 beta =  0.5 * (1.0-curant);
		gamma =  1.0 - beta;
		delta =  0.0;
	} 
	La[i] = alpha*lambda;
	Lb[i] =  beta*lambda;
	Lg[i] = gamma*lambda;
	Ld[i] = delta*lambda;
}

Vector4 CBGKRiemannSolver::calcFlux(FEOS& eos, Vector4 Um, Vector4 U, Vector4 Up, Vector4 Upp, double dx, double dt) {
	const double gamma = eos.getc(1.,1.)*eos.getc(1.,1.);
	Vector4 Fm  = getOmegaInv(eos, Um)*Um;
	Vector4 F   = getOmegaInv(eos, U)*U;
	Vector4 Fp  = getOmegaInv(eos, Up)*Up;
	Vector4 Fpp = getOmegaInv(eos, Upp)*Upp;
	Matrix4 O, OI;
	Vector4 L;
    Vector4 V1 = Vector4::ZERO, V2 = Vector4::ZERO, V3 = Vector4::ZERO, V4 = Vector4::ZERO;
	//Averaging
	// First use aryphmetic mean, then use Roe averaging, it'll be cool
	double rhoL = U[0], rhoR = Up[0]; 
	double uL = U[1]/rhoL, uR = Up[1]/rhoR;
	double EL = U[2]/rhoL, eL = EL - .5*uL*uL, ER = Up[2]/rhoR, eR = ER - .5*uR*uR;
	double pL = eos.getp(rhoL,eL), pR = eos.getp(rhoR,eR);
	double HL = EL +pL/rhoL, HR = ER + pR/rhoR;
	
	
	// Uncomment for aryphmetical mean averaging -- simply doesn't work on large discontinuities
/*	double _rho=.5*(rhoL+rhoR), _u=.5*(rhoL*uL+rhoR*uR)/_rho, _e=.5*(rhoL*eL+rhoR*eR), _p=eos.getp(_rho,_e), _c=eos.getc(_rho,_p),
		   _E = _rho*(_e+.5*_u*_u);
	Vector4 UAv = Vector4(_rho,_rho*_u,_rho*(_e+.5*_u*_u),0.); */
// Uncomment for Roe averaging	
	double _rho = sqrt(rhoL*rhoR);
	double _u = (uL*sqrt(rhoL) + uR*sqrt(rhoR))/(sqrt(rhoL)+sqrt(rhoR));
	double _H = (HL*sqrt(rhoL) + HR*sqrt(rhoR))/(sqrt(rhoL)+sqrt(rhoR));
	double _E = (EL*sqrt(rhoL) + ER*sqrt(rhoR))/(sqrt(rhoL)+sqrt(rhoR));
	double _p = (_H - _E) * _rho;	
	double _c = sqrt((gamma-1.)*(_H-.5*_u*_u));
//	double _c = eos.getc(_rho, _p);
	Vector4 UAv = Vector4(_rho,_rho*_u,_rho*_E,0.); 
	Vector4 parameter = Vector4(_rho, _u, _H, _c);

	O = getOmega(eos, parameter);
	OI = getOmegaInv(eos, UAv);
	L = Vector4(_u-_c,_u,_u+_c,0.);
	fillLambda(Fm, F, Fp, Fpp, L, dt/dx);
	Vector4 FBGK = Vector4::ZERO;
	FBGK = O*(La*(OI*Um)) + O*(Lb*(OI*U)) + O*(Lg*(OI*Up)) + O*(Ld*(OI*Upp));
    return FBGK;
}


int C1DMethodSamarskii::calc(C1DProblem& pr, FEOS& eos, C1DFieldPrimitive& fld) {
	int i = 0, counter = 0;
	double p_next_plus = 0., p_next_minus = 0., p_plus = 0., p_minus = 0.;
	auto&& W = fld.W, && newW = fld.newW, &&prevW = fld.prevW;
	auto&& x = fld.x, && newx = fld.newx;
	auto imin = fld.imin, imax = fld.imax;
	auto dm = fld.dm, dt = fld.dt;

	vector<double> g(fld.imax + 1);
	for (i = imin; i < imax; i++)  {
		double du_dm = (W[i+1][1] - W[i][1])/dm;
		g[i] = -.5*.000005*W[i][0]*fabs(du_dm)*(du_dm-fabs(du_dm));



		/*if (du_dm < 0)
			// g[i] = 6000.0 * W[i][0] * du * du; // Al;
			g[i] = .00001 * W[i][0] * du_dm * du_dm;  // -0.01 * W[i][0] * du + .001 * W[i][0] * du * du; // Al;
		else
			g[i] = 0;*/


		if (i == 500) {
			double pp = 0.;

		}



	}
	std::copy(W.begin(), W.end(), newW.begin());
	int itCounter = 0, maxIt = 30;
	const double eps = .01; // .0001;
	vector<double> diff(fld.imax + 1);
	do {
		itCounter++;
		if (itCounter > maxIt) {
			cout << endl <<
				"C1DMethodLagrange::calc() error: no convergence in 30 iterations" << std::endl;
			// dumpToFile(t);
			exit(1);
		}
		std::copy(newW.begin(), newW.end(), prevW.begin());
		for (i = imin; i < imax; i++) {
			p_next_plus = prevW[i][2] + g[i];
			p_next_minus = prevW[i-1][2] + g[i-1];
			p_plus = W[i][2] + g[i];
			p_minus = W[i-1][2] + g[i-1];
			/*
			if (i == 0) {
				p_next_plus = ms_prev[i].p;
				p_next_minus = 0.0;
				p_plus = ms[i].p;
				p_minus = 0.0;
			}
			else if (i == nSize) {
				p_next_plus = 0.0;
				p_next_minus = ms_prev[i - 1].p;
				p_plus = 0.0;
				p_minus = ms[i - 1].p;
			}
			else {
				p_next_plus = ms_prev[i].p + g[i];
				p_next_minus = ms_prev[i - 1].p + g[i - 1];
				p_plus = ms[i].p + g[i];
				p_minus = ms[i - 1].p + g[i - 1];
			}
			

			ms_temp[i].v = ms[i].v - tau / 2.0 / h * (p_next_plus - p_next_minus + p_plus - p_minus);
			ms_temp[i].x = ms[i].x + 0.5 * tau * (ms_temp[i].v + ms[i].v);

			*/
			newW[i][1] = W[i][1] - .5*dt/dm*(p_next_plus - p_next_minus + p_plus - p_minus);
			newx[i] = x[i] + .5*dt*(newW[i][1] + W[i][1]);
		}
		newx[imax] = x[imax] + .5*dt*(newW[imax][1] + W[imax][1]);
		newx[0] = newx[1] - (x[1] - x[0]);

		for (int i = imin; i < imax; i++) {
		/*	ms_temp[i].ro = 1.0 / (1.0 / ms[i].ro + tau / 2.0 / h *
				(ms_temp[i + 1].v + ms[i + 1].v - ms_temp[i].v - ms[i].v));
			ms_temp[i].ei = ms[i].ei - tau / 4.0 / h * (ms_prev[i].pi + ms[i].pi + g[i]) *
				(ms_temp[i + 1].v + ms[i + 1].v - ms_temp[i].v - ms[i].v);
			ms_temp[i].ee = ms[i].ee - tau / 4.0 / h * (ms_prev[i].pe + ms[i].pe + g[i]) *
				(ms_temp[i + 1].v + ms[i + 1].v - ms_temp[i].v - ms[i].v);
			ms_temp[i].e = ms_temp[i].ei + ms_temp[i].ee;*/
			newW[i][0] = 1. / (1. / W[i][0] + .5 * dt / dm * (newW[i+1][1] + W[i+1][1] - newW[i][1] - W[i][1]));
			double e = eos.gete(W[i][0], W[i][2]);
			double newe = e - .25 * dt / dm * (prevW[i][2] + W[i][2] + g[i]) * (newW[i+1][1] + W[i+1][1] - newW[i][1] - W[i][1]);
			newW[i][2] = eos.getp(newW[i][0], newe);
		}
		for (int i = 1; i < imax; i++) {
			const double M = 27.e-3, R = 8.31;
			double newti = newW[i][2] * M / (newW[i][0] * R);
			double prevti = prevW[i][2] * M / (prevW[i][0] * R);
			diff[i] = fabs(newti - prevti);
		}	
		diff[0] = 0.; diff[imax] = 0;
	} while (*std::max_element(diff.begin(), diff.end()) > eps);
	std::copy(newW.begin(), newW.end(), W.begin());
	std::copy(newx.begin(), newx.end(), x.begin());





	for(i = imin; i < imax; i++) {
		double c = eos.getc(W[i][0], W[i][2]);
		/*if (std::isnan(c)) {
			std::cout << "Acoustic disaster in the cell " << i << std::endl;	
			exit(1);
		}*/

		if (std::isnan(c)) {
			W[i][2] = 0.;
			std::cout << "\n\n\n\n\nSpallation in the cell " << i << "!!!\n\n\n\n\n";
			exit(1);
		}
	}





	pr.setbcs(fld.W);
	g.clear();
	diff.clear();
	return itCounter;
}

double C1DMethodSamarskii::calcdt(C1DProblem& pr, FEOS& eos, C1DFieldPrimitive& fld) {
	auto&& W = fld.W;
	auto&& x = fld.x;
	int imin = fld.imin, imax = fld.imax;
	double rho = W[imin][0], u = W[imin][1], p = W[imin][2], c = eos.getc(rho, p);
	double umax = max(fabs(u), max(fabs(u - c), fabs(u + c)));
	double dt1 = (x[imin + 1] - x[imin]) / umax;
	double dt2 = 0.;
	for (int i = imin; i < imax; ++i) {
		rho = W[i][0]; u = W[i][1]; p = W[i][1]; c = eos.getc(rho, p);
		umax = fabs(u) + c;
		dt2 = (x[i+1] - x[i]) / umax;
		if (dt1 > dt2) dt1 = dt2;
	}
	return pr.cfl * dt1;
}
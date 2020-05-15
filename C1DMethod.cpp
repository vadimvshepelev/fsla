#include"C1DMethod.h"
#include"solver.h"


void C1DGodunovMethodMillerPuckett::calc(C1DProblem& pr, CEOSMieGruneisen& eos, C1DField& fld) {
	double roL = 0., uL = 0., eL = 0., pL = 0.,
		   roR = 0., uR = 0., eR = 0., pR = 0., 
		   E = 0.;
	double dx = fld.dx, t = fld.t, dt = fld.dt;
	int imin = fld.imin, imax = fld.imax;
	vector<vector<double>> U = fld.U, newU = fld.newU, F = fld.F;
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
		double _F[] = {res.ro*res.u, res.ro*res.u*res.u + res.p, res.u*(res.ro*E+res.p)};
		F[i] = vector<double>(_F, _F+sizeof(_F)/sizeof(_F[0]));		
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

double C1DGodunovMethodMillerPuckett::calcdt(C1DProblem& pr, CEOSMieGruneisen& eos, C1DField& fld) {
	vector<vector<double>> U = fld.U; 
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

C1DVectorPrimitive C1DGodunovMethodMillerPuckett::calcRPExactMillerPuckett(CEOSMieGruneisen& eos, double roL, double uL, double pL, double roR, double uR, double pR) {
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
RPValues CExactRiemannSolver::calcValues(CEOS& eos, double roL, double uL, double pL, double roR, double uR, double pR) {
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

double CExactRiemannSolver::fL(CEOS& eos, double p, double roL, double uL, double pL) {
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

double CExactRiemannSolver::dfLdp(CEOS& eos, double p, double roL, double vL, double pL) {
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

double CExactRiemannSolver::fR(CEOS& eos, double p, double roR, double vR, double pR) {
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

double CExactRiemannSolver::dfRdp(CEOS& eos, double p, double roR, double vR, double pR) {
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

C1DVectorPrimitive CExactRiemannSolver::calcSolution(CEOS& eos, double roL, double uL, double pL, double roR, double uR, double pR, double x, double t){
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









double C1DGodunovTypeMethod::calcdt(C1DProblem& pr, CEOS& eos, C1DField& fld) {
	vector<vector<double>> U = fld.U; 
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

void C1DGodunovTypeMethod::calc(C1DProblem& pr, CEOS& eos, C1DField& fld) {
	double roL = 0., uL = 0., eL = 0., pL = 0.,
		   roR = 0., uR = 0., eR = 0., pR = 0., 
		   E = 0.;
	double dx = fld.dx, t = fld.t, dt = fld.dt;
	int imin = fld.imin, imax = fld.imax;
	vector<vector<double>> &U = fld.U, &newU = fld.newU, &F = fld.F;
	// TODO: проверить на скорость выполнения операций, сравнить с реализацией через тип Vector4 -- если не медленнее, то в дальнейшем избавиться от Vector4 везде
	int i=0;	
	// Потоки считаем по алгоритму решения задаче о распаде разрыва для УРС Ми-Грюнайзена из работы [Miller, Puckett]
	for(i=imin; i<=imax; i++) {
		Vector4 flux = rslv.calcFlux(eos, U[i-1][0], U[i-1][1], U[i-1][2], U[i][0], U[i][1], U[i][2]); 
		double _F[] = {flux[0], flux[1], flux[2]};
		F[i] = vector<double>(_F, _F+sizeof(_F)/sizeof(_F[0]));		
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

// Godunov-exact solver flux
Vector4 CExactRiemannSolver::calcFlux(CEOS& eos, double roL, double rouL, double roEL, double roR, double rouR, double roER) {	
	double uL = rouL/roL, uR = rouR/roR;
	double pL = eos.getp(roL, roEL/roL - .5*uL*uL), pR = eos.getp(roR, roER/roR - .5*uR*uR);
	C1DVectorPrimitive res = calcSolution(eos, roL, uL, pL, roR, uR, pR, 0., .01);
	double E = eos.gete(res.ro, res.p) + .5*res.u*res.u;
	Vector4 FGodunov = Vector4(res.ro*res.u, res.ro*res.u*res.u + res.p, res.u*(res.ro*E+res.p), 0.);
	return FGodunov;
}

// HLL flux for general EOS
Vector4 CHLLRiemannSolver::calcFlux(CEOS& eos, double roL, double rouL, double roEL, double roR, double rouR, double roER) {
	double _ro = 0., _u = 0., _e = 0., _p = 0.;
	double uL = rouL/roL, uR = rouR/roR;	
	double eL = roEL/roL - .5*uL*uL, eR = roER/roR - .5*uR*uR; 
	if(roL == 0. || roR == 0.) {
		cout << "Error: CSolver::calcHLLFluxEOSIdeal(): vacuum is present." << endl;
		exit(1);	}
	double pL = eos.getp(roL, eL), pR = eos.getp(roR, eR);
	double cL = eos.getc(roL, pL), cR = eos.getc(roR, pR);	
	Vector4 UL = Vector4(roL, rouL, roEL, 0.), UR = Vector4(roR, rouR, roER, 0.);
	_ro = roL, _u = rouL/roL, _e = roEL/roL-.5*_u*_u, _p=eos.getp(_ro, _e);
	Vector4 FL = Vector4(rouL, _p + _ro*_u*_u, _u*(_p + roEL), 0.);
	_ro = roR, _u = rouR/roR, _e = roER/roR-.5*_u*_u, _p=eos.getp(_ro, _e);
	Vector4 FR = Vector4(rouR, _p + _ro*_u*_u, _u*(_p + roER), 0.); 
	double SL = min(min(uL - cL, uR - cR), 0.), 
		   SR = max(max(uL + cL, uR + cR), 0.);
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
Vector4 CHLLCRiemannSolver::calcFlux(CEOS& eos, double roL, double rouL, double roEL, double roR, double rouR, double roER) {
	double _ro = 0., _u = 0., _e = 0., _p = 0.;
	double uL = rouL/roL, uR = rouR/roR;	
	double eL = roEL/roL - .5*uL*uL, eR = roER/roR - .5*uR*uR; 
	if(roL == 0. || roR == 0.) {
		cout << "Error: CSolver::calcHLLFluxEOSIdeal(): vacuum is present." << endl;
		exit(1);	}
	double pL = eos.getp(roL, eL), pR = eos.getp(roR, eR);
	double cL = eos.getc(roL, pL), cR = eos.getc(roR, pR);	
	Vector4 UL = Vector4(roL, rouL, roEL, 0.), UR = Vector4(roR, rouR, roER, 0.);
	_ro = roL, _u = rouL/roL, _e = roEL/roL-.5*_u*_u, _p=eos.getp(_ro, _e);
	Vector4 FL = Vector4(rouL, _p + _ro*_u*_u, _u*(_p + roEL), 0.);
	_ro = roR, _u = rouR/roR, _e = roER/roR-.5*_u*_u, _p=eos.getp(_ro, _e);
	Vector4 FR = Vector4(rouR, _p + _ro*_u*_u, _u*(_p + roER), 0.); 
	// Step 1 -- pressure estimate from primitive-variable Riemann solver (PVRS)
	double roAv = .5*(roL + roR), cAv = .5*(cL + cR), pPVRS = .5*(pL + pR) - .5*(uR - uL)*roAv*cAv;
	double pStar = max(0., pPVRS);
	// Step 2 -- wave speed estimates	
	double SL = min(min(uL - cL, uR - cR), 0.), SR = max(max(uL + cL, uR + cR), 0.);
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


Vector4 CGPSRiemannSolver::calcFlux(CEOS& eos, double roL, double rouL, double roEL, double roR, double rouR, double roER) {
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


void C1D2ndOrderMethod::calc(C1DProblem& pr, CEOS& eos, C1DField& fld) {
	double roL = 0., uL = 0., eL = 0., pL = 0., roR = 0., uR = 0., eR = 0., pR = 0., E = 0.;
	double dx = fld.dx, t = fld.t, dt = fld.dt;
	int i=0, imin = fld.imin, imax = fld.imax;
	vector<vector<double>> &U = fld.U, &newU = fld.newU, &F = fld.F;
	vector<vector<double>> &ULx = rec.ULx, &URx = rec.URx;
	// TODO: проверить на скорость выполнения операций, сравнить с реализацией через тип Vector4 -- если не медленнее, то в дальнейшем избавиться от Vector4 везде
	rec.calc(fld);
	for(i=imin; i<=imax; i++) {								
		Vector4 flux = rslv.calcFlux(eos, URx[i-1][0], URx[i-1][1], URx[i-1][2], ULx[i][0], ULx[i][1], ULx[i][2]); 
		double _F[] = {flux[0], flux[1], flux[2]};
		F[i] = vector<double>(_F, _F+sizeof(_F)/sizeof(_F[0]));		
	}
	for(i=imin; i<imax; i++) {
		for(int counter=0; counter<3; counter++) newU[i][counter] = U[i][counter] - dt/dx*(F[i+1][counter] - F[i][counter]);
	}
	for(i=imin; i<imax; i++) {
		for(int counter=0; counter<3; counter++) U[i][counter] = newU[i][counter];
	}	
	pr.setbcs(fld.U);
}
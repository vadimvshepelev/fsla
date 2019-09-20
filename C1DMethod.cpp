#include"C1DMethod.h"


void C1DGodunovMethodMillerPuckett::calc(C1DProblem& pr, CEOSMieGruneisen& eos, C1DField& fld) {
	double roL = 0., uL = 0., eL = 0., pL = 0.,
		   roR = 0., uR = 0., eR = 0., pR = 0., 
		   E = 0.;
	double dx = fld.dx, t = fld.t, dt = fld.dt;
	int imin = fld.imin, imax = fld.imax;
	vector<vector<double>> U = fld.U, newU = fld.newU, F = fld.F;
	// TODO: проверить на скорость выполнения операций, сравнить с реализацией через тип Vector4 -- если не медленнее, то в дальнейшем избавиться от Vector4 везде
	int i=0;
	CVectorPrimitive res;
	// Потоки считаем по алгоритму решения задаче о распаде разрыва для УРС Ми-Грюнайзена из работы [Miller, Puckett]
	for(i=imin; i<=imax; i++) {
		roL = U[i-1][0]; uL = U[i-1][1]/roL; eL = U[i-1][2]/roL - .5*uL*uL; pL = eos.getp(roL, eL);
		roR = U[i][0];   uR = U[i][1]/roR;   eR = U[i][2]/roR - .5*uR*uR;   pR = eos.getp(roR, eR);								
		res = calcRPExactMillerPuckett(eos, roL, uL, pL, roR, uR, pR);
		assert(res.ro!=0.);
		E = eos.gete(res.ro, res.p) + .5*res.v*res.v;
		double _F[] = {res.ro*res.v, res.ro*res.v*res.v + res.p, res.v*(res.ro*E+res.p)};
		F[i] = vector<double>(_F, _F+sizeof(_F)/sizeof(_F[0]));		
		//n.W_temp = n.W - dt/h*(Fp-Fm);
	}
	for(i=imin; i<imax; i++) {
		for(int counter=0; counter<3; counter++) newU[i][counter] = U[i][counter] - dt/dx*(F[i+1][counter] - F[i][counter]);
	}
	for(i=imin; i<imax; i++) {
		for(int counter=0; counter<3; counter++) U[i][counter] = newU[i][counter];
	}	
	fld.setbcs(pr);
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

CVectorPrimitive C1DGodunovMethodMillerPuckett::calcRPExactMillerPuckett(CEOSMieGruneisen& eos, double roL, double uL, double pL, double roR, double uR, double pR) {
	CVectorPrimitive V;
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
				V.v = uCD;
				//_e = eL + .5*(pL + pCD)*(1./roL-1./V.ro);
				//V.p = eos.getp(V.ro, _e);
				V.p = pCD;			
			} else {
				V.ro = roL;
				V.v = uL;
				V.p = pL;
			}
		} else {
		    // Left rarefaction
			double roCDL = roL/(1.+(pL-pCD)/KSL);
			double cCDL = a0L*roL/roCDL;
			double sigmaLWave = (a0L-uL)/(a0L-uL-uCD-cCDL);
			double sigmaL = min(1., max(0., sigmaLWave));
			V.ro = (1.-sigmaL)*roL + sigmaL*roCDL;
			V.v = (1.-sigmaL)*uL + sigmaL*uCD;
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
				V.v = uCD;
				//_e = eL + .5*(pL + pCD)*(1./roL-1./V.ro);
				//V.p = eos.getp(V.ro, _e);
				V.p = pCD;			
			} else {
				V.ro = roR;
				V.v = uR;
				V.p = pR;
			}
		} else {
		    // Right rarefaction
			double roCDR = roL/(1.+(pR-pCD)/KSR);
			double cCDR = a0R*roR/roCDR;
			double sigmaRWave = (a0R-uR)/(a0R-uR-uCD-cCDR);
			double sigmaR = min(1., max(0., sigmaRWave));
			V.ro = (1.-sigmaR)*roR + sigmaR*roCDR;
			V.v = (1.-sigmaR)*uL + sigmaR*uCD;
			V.p = (1.-sigmaR)*pL + sigmaR*pCD;
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
		roL = U[i-1][0]; uL = U[i-1][1]/roL; eL = U[i-1][2]/roL - .5*uL*uL; pL = eos.getp(roL, eL);
		roR = U[i][0];   uR = U[i][1]/roR;   eR = U[i][2]/roR - .5*uR*uR;   pR = eos.getp(roR, eR);								
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
	fld.setbcs(pr);
}

// HLL flux for ideal EOS
Vector4 CHLLRiemannSolver::calcFlux(CEOS& eos, double roL, double rouL, double roEL, double roR, double rouR, double roER) {
	double _ro = 0., _u = 0., _e = 0., _p = 0.;
	double uL = rouL/roL, uR = rouR/roR;	
	double eL = roEL/roL - .5*uL*uL, eR = roER/roR - .5*uR*uR; 
	if(roL == 0. || roR == 0.) {
		cout << "Error: CSolver::calcHLLFluxEOSIdeal(): vacuum is present." << endl;
		exit(1);	}
	double pL = eos.getp(roL, eL), pR = eos.getp(roR, eR);
	double cL = eos.getc(roL, pL), cR = eos.getc(roR, pR);	
/*	double pMin = min(pL, pR), pMax = max(pL, pR), pStar = 1./(cL+cR)*(cR*pL + cL*pR + cL*cR*(uL-uR));
	double Q = pMax/pMin, QUser = 2.;
	double roLRes = 0., roRRes = 0., uRes = 0., pRes = 0.;
	double qL = 0., qR = 0.;
	double z = 0., pLR = 0.;
	double AL = 0., AR = 0., BL = 0., BR = 0., gL = 0., gR = 0., p_0 = 0.;
	if(Q < QUser && pStar > pMin && pStar > pMax) {
		// Primitive-variable noniteration solver (PVRS)
		pRes = pStar;
		uRes = 1./(cL+cR)*(cL*uL + cL*cR + pL-pR);
		roLRes  = roL + (pStar-pL)/cL/cL;
		roRRes  = roR + (pStar-pR)/cR/cR;
	} else if (pStar < pMin) {
		// Two-rarefaction Riemann solver (TRRS)	
		z = (gamma-1.)/2./gamma;
		pLR = pow(pL/pR, z);
		pRes = pow( (cL+cR - (gamma-1.)/2.*(uR-uL)) / (cL/pow(pL, z) + cR/pow(pR, z)), 1./z);
		uRes = (pLR*uL/cL + uR/cR + 2*(pLR-1)/(gamma-1.))/(pLR/cL + 1./cR);
		roLRes = roL*pow(pRes/pL, 1./gamma);
		roRRes = roR*pow(pRes/pR, 1./gamma);
	} else {
		// Two-shock Riemann solver (TSRS)
		AL = 2./(gamma+1)/roL; AR = 2./(gamma+1)/roR;
		BL = (gamma-1.)/(gamma+1)*pL; BR = (gamma-1.)/(gamma+1.)*pR;
		p_0 = max(0., pStar);
		gL = sqrt(AL/(p_0+BL)); gR = sqrt(AR/(p_0+BR));
		pRes = (gL*pL + gR*pR - (uR-uL))/(gL+gR);
		uRes = 0.5*(uR+uL) + 0.5*((pRes-pR)*gR - (pRes-pL)*gL);
		roLRes = (pRes/pL + (gamma-1.)/(gamma+1)) / ((gamma-1.)/(gamma+1.)*pRes/pL + 1.);
		roRRes = (pRes/pR + (gamma-1.)/(gamma+1)) / ((gamma-1.)/(gamma+1.)*pRes/pR + 1.);
	}
	if(pRes <= pL) qL = 1.; else qL = sqrt(1. + (gamma+1.)/2./gamma*(pRes/pL-1.));
	if(pRes <= pR) qR = 1.; else qR = sqrt(1. + (gamma+1.)/2./gamma*(pRes/pR-1.));
	//double SL = min(uL-cL, uR-cR), SR = max(uL+cL, uR+cR);
	double SL = uL - cL*qL, SR = uR - cR*qR;
	double _ro=0, _u=0., _e=0., _p=0.;
	Vector4 UL = Vector4(roL, rouL, roEL, 0.), UR = Vector4(roR, rouR, roER, 0.);
	_ro = roL, _u = rouL/roL, _e = roEL/roL-.5*_u*_u, _p=eos.getp(_ro, _e);
	Vector4 FL = Vector4(rouL, _p + _ro*_u*_u, _u*(_p + roEL), 0.);
	_ro = roR, _u = rouR/roR, _e = roER/roR-.5*_u*_u, _p=eos.getp(_ro, _e);
	Vector4 FR = Vector4(rouR, _p + _ro*_u*_u, _u*(_p + roER), 0.); 
	if(0 <= SL) 
		return FL;
	else if (0 >= SR )
		return FR;
	else if (SL<=0 && SR >=0) {		
		return 1./(SR-SL)*(SR*FL - SL*FR - SL*SL*(UR-UL));
	} else {
	   cout << "Error: C1DGodunovMethod::calcFlux(): unexpected wave configuration." << endl;
	   exit(1);	
	}*/
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

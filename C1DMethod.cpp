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
	const double K0S = eos.ro0*eos.getc(eos.ro0, eos.getp0(eos.ro0));
	const double eL = eos.gete(roL, pL), eR = eos.gete(roR, pR), cL = eos.getc(roL, eL), cR = eos.getc(roR, eR); 
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
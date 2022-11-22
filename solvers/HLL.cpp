#include "../solver.h"
#include "../eos/EOSBin.h"

Vector4 CSolver::calcHLLFluxEOSBin(double roL, double rouL, double roEL, double roR, double rouR, double roER) {
	EOSBin &eos = (task.eosBin);	
	const double gamma = eos.gamma;
	double uL = rouL/roL, uR = rouR/roR;	
	double eL = roEL/roL - .5*uL*uL, eR = roER/roR - .5*uR*uR; 
	if(roL == 0. || roR == 0.) {
		cout << "Error: CSolver::calcHLLFluxEOSBin(): vacuum is present." << endl;
		exit(1);
	}
	double cL = eos.getC(roL, eL), cR = eos.getC(roR, eR);
	double pL = eos.getp(roL, eL), pR = eos.getp(roR, eR);
	double pMin = min(pL, pR), pMax = max(pL, pR), pStar = 1./(cL+cR)*(cR*pL + cL*pR + cL*cR*(uL-uR));
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
	Vector4 UL = Vector4(roL, rouL, roEL, 0.), 
		    UR = Vector4(roR, rouR, roER, 0.),
	        FL = calcPhysicalFluxEOSBin(roL, rouL, roEL), 
			FR = calcPhysicalFluxEOSBin(roR, rouR, roER); 
	if(0 <= SL) 
		return FL;
	else if (0 >= SR )
		return FR;
	else if (SL<=0 && SR >=0) {		
		return 1./(SR-SL)*(SR*FL - SL*FR - SL*SL*(UR-UL));
	} else {
	   cout << "Error: CSolver::calcHLLFluxEOSBin(): unexpected wave configuration." << endl;
	   exit(1);	
	}
}


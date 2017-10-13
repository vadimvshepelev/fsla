#include "..\\solver.h"
#include "..\\EOS\\EOSBin.h"

Vector4 CSolver::calcHLLFluxEOSBin(double roL, double rouL, double roEL, double roR, double rouR, double roER) {
	EOSBin eos = *(task.eosBin);	
	double uL = rouL/roL, uR = rouR/roR;	
	double eL = roEL/roL - .5*uL*uL, eR = roER/roR - .5*uR*uR; 
	double cL = eos.getC(roL, eL), cR = eos.getC(roR, eR);
	/*double pL = eos.getp(roL, eL), pR = eos.getp(roR, eR);*/
	double SL = min(uL-cL, uR-cR), SR = max(uL+cL, uR+cR);
	Vector4 UL = Vector4(roL, rouL, roEL, 0.), 
		    UR = Vector4(roR, rouR, roER, 0.),
	        FL = calcPhysicalFluxEOSBin(roL, rouL, roEL), 
			FR = calcPhysicalFluxEOSBin(roR, rouR, roER); 
	if(0 <= SL) 
		return FL;
	else if (0 >= SR )
		return FR;
	else if (SL<=0 && SR >=0) 
		return 1./(SR-SL)*(SR*FL - SL*FR - SL*SL*(UR-UL));
	else 
	   exit(1);	
}


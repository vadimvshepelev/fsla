#include "../solver.h"
#include "../eos/EOSBin.h"

Vector4 CSolver::calcHLLCFluxEOSBin(double roL, double rouL, double roEL, double roR, double rouR, double roER) {
	EOSBin &eos = task.eosBin;	
	double uL = rouL/roL, uR = rouR/roR;
	double eL = roEL/roL - .5*uL*uL, eR = roER/roR - .5*uR*uR; 
	double pL = eos.getp(roL, eL), pR = eos.getp(roR, eR);
	double SL = 0., SR = 0.;


	Vector4 FHLLC = Vector4::ZERO;//Vector4(res.ro*res.v, res.ro*res.v*res.v + res.p, res.v*(res.ro*E+res.p), 0.);
	return FHLLC;
}

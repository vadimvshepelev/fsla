#include <math.h>

#include "ceos.h"


double CEOSMieGruneisen::getp(double ro, double e) {
	const double A = .7626e9; // [Pa]
	const double b = 11.55;   	
	const double K = 1.15e9;  // [Pa]
	const double beta = .3333; 
	const double xi = .85; 
	// Cold component, Born-Meyer potential 
	double pp = A*pow(ro/ro0, -beta+1.) * exp(b*(1.-pow(ro/ro0, -beta))) - K*(pow(ro/ro0, xi+1.)); 



}

double CEOSMieGruneisen::gete(double ro, double p) {



}
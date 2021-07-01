#include <math.h>
#include <iostream>
#include <string>
#include <sstream>

using namespace std;

#include "feos.h"


double FEOSMieGruneisen::getG(double ro) {
	double x = ro/ro0;
	const double a0 = 2.95, a1 = 2.408, a2 = 12.151;
	const double M = 18.;     // [g/mole]
	const double CVLiq = 4150.;   // [J/kg/K] 
	const double CVGas = 1430.;	  // [J/kg/K]
	const double R = 8.31;        // [J/mole/K]
	// Здесь используем только жидкую фазу! Можно будет дифференцировать при усложнении
	const double CV = CVLiq;
	double G = R/CV/M*(a0 + (1.-a0)*exp(-pow(x/.5273, 1.7)) + a1*exp(-pow(x/1.0904, -3.5)) + a2*exp(-pow(x/1.3927, -5.)));
 	return G;
}

double FEOSMieGruneisen::getp0(double ro) {
	const double A = .6726e9; // [Pa]
	const double b = 11.55;   	
	const double K = 1.15e9;  // [Pa]
	const double beta = .3333; 
	const double xi = .85; 
	// Cold component, Born-Meyer potential 
	double x = ro/ro0;
	double p0 =  x*(A*pow(x, -beta) * exp(b*(1.-pow(x, -beta))) - K*(pow(x, xi))); 
	return p0;
}

double FEOSMieGruneisen::gete0(double ro) {
	const double A = .6726e9; // [Pa]
	const double b = 11.55;   	
	const double K = 1.15e9;  // [Pa]
	const double beta = .3333; 
	const double xi = .85; 
	// Cold component, Born-Meyer potential 
	double x = ro/ro0;
	double e0 = 1./ro0*(A/beta/b * exp(b*(1.-pow(x, -beta))) - K/xi*pow(x, xi));
	return e0;
}

double FEOSMieGruneisen::getp(double ro, double e) {
	double x = ro/ro0;
	double p0 = getp0(ro);
	double e0 = gete0(ro);
	double G = getG(ro);
	double p = p0 + ro*G*(e-e0);
	return p;
}

double FEOSMieGruneisen::gete(double ro, double p) {
	double x = ro/ro0;
	double p0 = getp0(ro);
	double e0 = gete0(ro);
	double G = getG(ro);
	double e = e0 + 1./ro/G*(p-p0);
	return e;
}

double FEOSMieGruneisen::getc(double ro, double p) {
	double x = ro/ro0;
	double G = getG(ro);
	double p0 = getp0(ro);
	double p0Prime = getp0Prime(ro);
	double GPrime = getGPrime(ro);
	double c2 = 1/ro0*(p0Prime + (GPrime/G + 1./x + G/x)*(p-p0));
	return sqrt(c2);
}


double FEOSMieGruneisen::getGPrime(double ro) {
	double x = ro/ro0;
	const double a0 = 2.95, a1 = 2.408, a2 = 12.151;
	const double M = 18.;     // [g/mole]
	const double CVLiq = 4150.;   // [J/kg/K] 
	const double CVGas = 1430.;	  // [J/kg/K]
	const double R = 8.31;        // [J/mole/K]
	// Здесь используем только жидкую фазу! Можно будет дифференцировать при усложнении
	const double CV = CVLiq;
	/* double GPrime = R/CV/M*((1.-a0)*exp(-pow(x/.5273, 1.7))*pow(x/.5273, -2.7)*(-1.7/.5273) + 
		                        a1 *exp(-pow(x/1.0904, -3.5))*pow(x/1.0904, -4.5)*(3.5/1.0904) +
								a2 *exp(-pow(x/1.3927, -5.0))*pow(x/1.3927, -6.0)*(5.0/1.3927)); */
	double GPrime = .035412*exp(-5.23948*pow(x, -5.))*pow(x, -6) 
		          + .00126927*exp(-1.35379*pow(x, -3.5))*pow(x, -4.5) 
				  + .00109463*exp(-2.96826*pow(x, 1.7))*pow(x, .7);
	return GPrime;
}
	
double FEOSMieGruneisen::getp0Prime(double ro) {
	const double A = .7626e9; // [Pa]
	const double b = 11.55;   	
	const double K = 1.15e9;  // [Pa]
	const double beta = .3333; 
	const double xi = .85; 
	double x = ro/ro0;
	double p0Prime = /*A*pow(x, -beta+1.) * exp(b*(1.-pow(x, -beta))) - K*(pow(x, xi+1.)) + 
		             A*(-beta+1.)*pow(x, -beta)*exp(b*(1.-pow(x, -beta))) + A*pow(x, -beta+1.) * exp(b*(1.-pow(x, -beta))) * 		
		x*(         A*pow(x, -beta+1.) * exp(b*(1.-pow(x, -beta))) - K*(pow(x, xi+1.))); 
		*/
		/*A*pow(x, -beta)*exp(b*(1.-pow(x, -beta))) - K*pow(x, xi) + 
		A*beta*exp(b*(1.-pow(x, -beta)))*(-pow(x, -beta) + pow(x, -2.*beta)) - K*xi*pow(x, xi);*/

		exp(-11.55/pow(x,.3333))*(2.68705e14/pow(x, .6666) + 4.65359e13/pow(x,.3333)) - 2.1275e9*pow(x,.85);
	return p0Prime;
}

double FEOSMieGruneisen::getKSPrime(double ro, double e) {
	double p = getp(ro, e), c = getc(ro, p), KS = ro*c*c, G = getG(ro), x = ro/ro0; 
	// Partial derivatives
	// G'(ro)
	double a0 = 2.95, a1 = 2.408, a2 = 12.151;
	double b0 = .5273, b1 = 1.0904, b2 = 1.3927;
	double c0 = 1.7, c1 = -3.5, c2 = -5.;
	double dGdx = - c0*(1.-a0)/b0*exp(-pow(x/b0, c0))*pow(x/b0, c0-1.) - a1*c1/b1*exp(-pow(x/b1, c1))*pow(x/b1, c1-1.) - a2*c2/b2*exp(-pow(x/b2, c2))*pow(x/b2, c2-1.);
	double dGdro = ro0*dGdx;
	// G"(ro)
	double d2Gdx2 = - c0*(1.-a0)/b0/b0*exp(-pow(x/b0, c0))*(-c0*pow(x/b0, 2.*c0-2.)+(c0-1.)*pow(x/b0, c0-2.))
		                 - a1*c1/b1/b1*exp(-pow(x/b1, c1))*(-c1*pow(x/b1, 2.*c1-2.)+(c1-1.)*pow(x/b1, c1-2.))
						 - a2*c2/b2/b2*exp(-pow(x/b2, c2))*(-c2*pow(x/b2, 2.*c2-2.)+(c2-1.)*pow(x/b1, c2-2.));
	double d2Gdro2 = ro0*ro0*d2Gdx2;
	// p'(ro)
	double A = .6726e9, beta = .3333, b = 11.55, K = 1.15e9, xi = .85;
	double dp0dx = A*pow(x, -beta)*exp(b*(1.-pow(x, -beta))) - K*pow(x, xi) + A*beta*exp(b*(1.-pow(x, -beta)))*(-pow(x, -beta) + pow(x, -2.*beta)) - K*xi*pow(x, xi);
	double dp0dro = ro*dp0dx;
	double e0 = 1000./ro0*(A/beta/b * exp(b*(1.-pow(x, -beta))) - K/xi*pow(x, xi));
	double de0dx = 1./ro0*(A*exp(b*1.-pow(x, -beta))*pow(x, -beta-1.) - K*pow(x, xi-1));
	double de0dro = ro0*de0dx;	
	double dpdro = dp0dro + (G+ro*dGdro)*(e-e0) - ro*G*de0dro;
    // p'(e) 
	double dpde = ro*G;
	double d2p0dx2 = A*exp(b*(1.-pow(x, -beta)))*(-beta*(1.-beta)*pow(x, -beta-1.) + (b*beta-(2.+b)*beta*beta)*pow(x, -2.*beta-1.) + b*beta*beta*pow(x, -3.*beta-1.)) - K*xi*(xi+1.)*pow(x, xi-1.);
	double d2p0dro2 = ro0*ro0*d2p0dx2;
	double d2e0dx2 = 1./ro0*(A*exp(b*(1-pow(x, -beta)))*(b*beta*pow(x, 2.*beta-2.) - (beta+1.)*pow(x, -beta-2.)) - K*pow(x, xi-2.));
	double d2e0dro2 = ro0*ro0*d2e0dx2;
	double d2pdro2 = d2p0dro2 + (2.*G+ro*d2Gdro2)*(e-e0) - (G+ro*dGdro)*de0dro - ro*G*d2e0dro2;
	double d2pdrode = G + ro*dGdro;
	double d2pde2 = 0.;
	double KSPrime = (ro*dpdro + ro*ro*d2pdro2 + dpdro*dpde + 2.*d2pdrode - p/ro*dpde - p/ro/ro*dpde*dpde + p*p/ro/ro*d2pde2)/KS;
	return KSPrime;
}


double FEOSMieGruneisenAl::getp(double ro, double e) {
	const double gamma = 3.9, G = 2., 
		         B = 76.e9; // Bulk modulus of Al, [Pa];
	double x = ro/ro0;
	// p = (pc - G*ec) + G*E
	// e c = (B/ (gamma (gamma -1 ) )) x^ gamma – (B /( gamma -1 )) x + B/ gamma
	double p_c = B/gamma * (pow(x, gamma)-1.);
	double e_c = (B/gamma/(gamma-1.)*pow(x, gamma) - B*x/(gamma-1.) + B/gamma)/ro;
	double p  = (p_c-G*ro*e_c) + G*ro*e;
	return p;
}

double FEOSMieGruneisenAl::gete(double ro, double p) {
	const double gamma = 3.9, G = 2., 
		         B = 76.e9; // Bulk modulus of Al, [Pa];
	double x = ro/ro0;
	double p_c = B/gamma * (pow(x, gamma)-1.);
	double e_c = (B/gamma/(gamma-1.)*pow(x, gamma) - B*x/(gamma-1.) + B/gamma)/ro;
	double e = (p-p_c+G*ro*e_c)/G/ro;
	return e;
}

double FEOSMieGruneisenAl::getc(double ro, double p) {
	const double gamma = 3.9, G = 2., 
		         B = 76.e9; // Bulk modulus of Al, [Pa];
	double x = ro/ro0;
	double p_c = B/gamma * (pow(x, gamma)-1.);
	double c  = sqrt((G+1)*(p-p_c)/ro + B/ro0*pow(x, gamma-1));
	return c;
}

double FEOSMGAlPrecise::G(double x) {
	if(x > 1.) 
		return 1./3. + .5*(a*(a+1.)*pow(x, a) - b*(b+1.)*pow(x, b))/((a+1.)*pow(x, a) - (b+1.)*pow(x, b));
	double B0 = a+b+1./3., 
		   B1 = -(a+1.)*(b+1.)/(a+b+1.), 
		   B2 = a*b+a+b+2.-a*a*b*b/(a+b+1.)/(a+b+1.);
	double g1 = 3.*B0 + (B0 + 2./3.)*(-2.*B1+(B1*B1 + B2)/2.),
		   g2 = -3.*B0 + (B0+2./3.)*(3.*B1-B1*B1-B2),
		   g3 = B0 - (B0+2./3.)*(B1-(B1*B1+B2)/2.);
	return 2./3. + .5*(g1*x + g2*x*x + g3*x*x*x);
}

double FEOSMGAlPrecise::Gx1(double x) {
	double B0 = a+b+1./3., 
		   B1 = -(a+1.)*(b+1.)/(a+b+1.), 
		   B2 = a*b+a+b+2.-a*a*b*b/(a+b+1.)/(a+b+1.);
	double g1 = 3.*B0 + (B0 + 2./3.)*(-2.*B1+(B1*B1 + B2)/2.),
		   g2 = -3.*B0 + (B0+2./3.)*(3.*B1-B1*B1-B2),
		   g3 = B0 - (B0+2./3.)*(B1-(B1*B1+B2)/2.);
	double Gx=1./2.*(3.*g3*x*x+2.*g2*x+g1);
	return Gx;
}

double FEOSMGAlPrecise::Gx2(double x) {
	double Gx = 1./(2.*x)* (
		        (a*a*(a+1.)*pow(x,a)-b*b*(b+1.)*pow(x,b)) / ((a+1.)*pow(x,a)-(b+1.)*pow(x,b)) -
				( (a*(a+1.)*pow(x,a)-b*(b+1.)*pow(x,b)) / ((a+1.)*pow(x,a)-(b+1.)*pow(x,b)) ) *
				( (a*(a+1.)*pow(x,a)-b*(b+1.)*pow(x,b)) / ((a+1.)*pow(x,a)-(b+1.)*pow(x,b)) ) 
				);
	return Gx;
}

double FEOSMGAlPrecise::GPrime(double x) {
	if(x<1.)
		return Gx1(x);
	else
		return Gx2(x);
}

double FEOSMGAlPrecise::pCold(double rho) {
	double x = rho/rho0;
	double pc = 0.;
	if(x>=1.)
	    pc=p0*x*(pow(x,a)-pow(x,b));
	else
		pc=(-17.1387+42.1454*x*x*x*x*x*x*x-30.2134*x*x*x*x*x*x*x*x*x+5.20668*pow(x,11.8389))*1.e9;
	return pc;
}

double FEOSMGAlPrecise::pColdPrime(double rho) {
	double x = rho/rho0;
	double pcx = 0.;
	if(x>=1.)
	    pcx = p0*((a+1.)*pow(x,a)-(b+1.)*pow(x,b));
	else
		pcx = (42.1454*7.*x*x*x*x*x*x-30.2134*9.*x*x*x*x*x*x*x*x+5.20668*11.8389*pow(x,10.8389))*1.e9;
	return pcx;
}

double FEOSMGAlPrecise::eCold(double rho) {  // [J/kg]
	double x = rho/rho0;
	double ec = 0.;
	if(x>=1.) 
		ec = 2.03987e8*((pow(x,a)/a-pow(x,b)/b-(1./a-1./b))); 		
	else
		ec = 2.03987e17/p0*(-17.1387*(-1./x + 1.) + 42.1454/6.*(x*x*x*x*x*x - 1.) - 30.2134/8.*(x*x*x*x*x*x*x*x - 1.) + 5.20668/10.8389*(pow(x,10.8389)-1.));
	return ec;
}


double FEOSMGAlPrecise::getp(double rho, double e) {
	double x = rho/rho0;
	return pCold(rho)+G(x)*rho*(e-eCold(rho));
}

double FEOSMGAlPrecise::gete(double rho, double p) {
	double x = rho/rho0;
	return eCold(rho)+1./rho/G(x)*(p-pCold(rho));
}

double FEOSMGAlPrecise::getc(double rho, double p) {
	double x = rho/rho0;
	double _G = G(x);
	double G_x = GPrime(x);
	double pc = pCold(rho)/p0;
	double pcx = pColdPrime(rho)/p0;
	double c2 = p0/rho0*((p/p0 - pc)/(x*G(x))*(_G+_G*_G+x*G_x) + pcx);
	return sqrt(c2);
}


double FEOSMGAlPrecise6::G(double x) {
	return 2.;
}

double FEOSMGAlPrecise6::GPrime(double x) {
	return 0.;
}

double FEOSMGAlPrecise6::pCold(double rho) {
	double x = rho/rho0;
	double pc = 0.;
	if(x>=1.)
	    pc=p0*x*(pow(x,a)-pow(x,b));
	else {
		double n = p0*(a-b)/ph;
		pc=ph*(pow(x,n)-1.);
	}
	return pc;
}

double FEOSMGAlPrecise6::pColdPrime(double rho) {
	double x = rho/rho0;
	double pcx = 0.;
	if(x>=1.)
	    pcx = p0*((a+1.)*pow(x,a)-(b+1.)*pow(x,b));
	else {
		double n = p0*(a-b)/ph;
		pcx = p0*(a-b)*pow(x,n-1.);
	}
	return pcx;
}

double FEOSMGAlPrecise6::eCold(double rho) {  // [J/kg]
	double x = rho/rho0;
	double ec = 0.;
	if(x>=1.) 
		ec = 2.03987e8*(pow(x,a)/a-pow(x,b)/b); 		
	else {
		double n = p0*(a-b)/ph;
		ec = 2.03987e8*(a-b)*(1./n*(pow(x,n-1.)/(n-1.)+1./x) - 1./(n-1.) - 1./a/b);    ;
	}
	return ec;
}


double FEOSMGAlPrecise6::getp(double rho, double e) {
	double x = rho/rho0;
	return pCold(rho)+G(x)*rho*(e-eCold(rho));
}

double FEOSMGAlPrecise6::gete(double rho, double p) {
	double x = rho/rho0;
	double _ec = eCold(rho);
	double _pc = pCold(rho);
	return _ec + 1./rho/G(x)*(p-_pc);
}

double FEOSMGAlPrecise6::getc(double rho, double p) {
/*	double x = rho/rho0;
	double _G = G(x);
	double G_x = GPrime(x);
	double pc = pCold(rho)/1.e9;
	double pcx = pColdPrime(rho)/1.e9;
	double c2 = 10.e9/rho0*((p/1.e9 - pc)/(x*G(x))*(_G+_G*_G+x*G_x) + pcx); */
	double _e = gete(rho,p);
	double _prho = getdpdrho(rho,_e); 
	double _pe = getdpde(rho,_e);
	double c2 = p*_pe/rho/rho + _prho;
	return sqrt(c2);
}

double FEOSMGAlPrecise6::gets(double rho, double p) {
	// cv = 1. (хотя не важно), x1 = 1., p1 = 0. -- референсное значение для энтропии s1
	const double x1 = 1., p1 = 0., cv = 1.;
	double x = rho/rho0;
	return cv*log( (p-pCold(x))/(p1-pCold(1.))*pow(x1/x,G(x)+1.) );
}

double FEOSMGAlPrecise6::eColdPrime(double rho) {
	double x = rho/rho0;
	double ecx = 0.;
	if(x>=1.) 
		ecx = 2.03987e8*(pow(x,a-1.)-pow(x,b-1.)); 		
	else {
		double n = p0*(a-b)/ph;
		ecx = 2.03987e17/p0*(a-b)*(1./n*(pow(x,n-2.)-1./x/x)); 
	}
	return ecx;
}

double FEOSMGAlPrecise6::getdpdrho(double rho, double e) {
	double x = rho/rho0;
	double _G = G(x);
	double _pc = pCold(rho);
	double _pcx = pColdPrime(rho);
	double _ec = eCold(rho);
	// double _ecx = eColdPrime(rho);
	// return _pcx/rho0 + _G*(e-_ec) - rho*_G*_ecx/rho0; 
	return 1./rho0*(_pcx + _G*rho0*(e-_ec) - _G/x*_pc); 
}


double FEOSMGAlPrecise6::getdpde(double rho, double e) {
	double _G = G(rho/rho0);
	return rho*_G;
}


















// Значения констант для широкодиапазонных (по-видимому) УРС Ломоносова
double __V0[] = { 0.12700000E+00, 0.88200003E-01, 
	            0.54200000E+00, 0.57700002E+00, 0.36899999E+00, 0.12000000E+00,
			    0.13900000E+00, 0.11300000E+00, 0.11200000E+00, 0.11200000E+00,
	            0.88200003E-01, 0.97800002E-01, 0.11400000E+00, 0.59700001E-01,
			    0.51899999E-01, 0.91200002E-01, 0.22200000E+00, 0.51800001E-01,
			    0.37900001E+00, 0.25099999E+00, 0.23400000E+00, 0.11200000E+01  };

double __E0[] = { 0.13392337E+00, 0.36263671E-01, 
	            0.82284731E+00, 0.30740389E+00, 0.27874166E+00, 0.13305190E+00,
				0.14209883E+00, 0.12565036E+00, 0.12638120E+00, 0.11601043E+00,
	            0.36263671E-01, 0.76777481E-01, 0.66890948E-01, 0.40647842E-01,
				0.40003933E-01, 0.90129100E-01, 0.15371029E+00, 0.37791826E-01,
				0.29866642E+00, 0.71235187E-01,	0.12139247E+00, 0.84046292E+00 }; 

double __DX[] = { 0.98664743E+00, 0.97725165E+00, 
	            0.98145735E+00, 0.97675198E+00, 0.97923738E+00, 0.98654467E+00,
				0.98867989E+00, 0.98776722E+00, 0.98788631E+00, 0.98512793E+00,
	            0.97725165E+00, 0.99516839E+00, 0.97772300E+00, 0.99392861E+00,
				0.99559391E+00, 0.97596389E+00, 0.99185765E+00, 0.98891979E+00,
				0.98592889E+00, 0.99903792E+00, 0.99866945E+00, 0.95384234E+00 };

double __GM[] = { 0.30000000E+01, 0.20000000E+01,
                0.30000000E+01, 0.30000000E+01, 0.30000000E+01, 0.30000000E+01,
				0.30000000E+01, 0.30000000E+01, 0.20000000E+01, 0.20000000E+01,
				0.20000000E+01, 0.30000000E+01, 0.30000000E+01, 0.30000000E+01,
				0.30000000E+01, 0.90000000E+01, 0.30000000E+01, 0.30000000E+01,
				0.30000000E+01, 0.30000000E+01, 0.10000000E+01, 0.30000000E+01 };

double __CMN[] = { 0.95555325E+01,-0.14261333E+02,
	             0.42994564E+02, 0.11484624E+02, 0.13236353E+02, 0.93908024E+01,
				 0.17328915E+02, 0.13348364E+02, 0.40190891E+02, 0.32090862E+02,
				 -.14261333E+02, 0.15142752E+02, 0.16921934E+02, 0.55232458E+01,
				 0.88132391E+01, 0.94326097E+00, 0.10463772E+02, 0.80316305E+01,
				 0.95852518E+01, 0.27481421E+02, -.65462183E+03, 0.10416588E+02 };

double __GN[]  = { 0.90187567E+00, 0.23021955E+01,
				 0.14658144E+01, 0.11745121E+01, 0.79678905E+00, 0.88841671E+00,
				 0.14138776E+01, 0.12774221E+01, 0.14686730E+01, 0.15083531E+01,
				 0.23021955E+01, 0.12716897E+01, 0.25520797E+01, 0.89647335E+00,
				 0.11707672E+01, 0.39387769E+00, 0.79043901E+00, 0.18141516E+01,
				 0.11521416E+00, 0.30822426E+00, 0.11392220E+01, 0.38477069E+00	};

double __ES[] = { 0.74099998E+01, 0.93599999E+00,
				0.15000000E+02, 0.59499998E+01, 0.12200000E+02, 0.74400001E+01,
				0.64800000E+01, 0.60000000E+01, 0.72700000E+01, 0.52300000E+01,
				0.93599999E+00, 0.68600001E+01, 0.99000001E+00, 0.43200002E+01,
				0.45900002E+01, 0.22900000E+01, 0.97500000E+01, 0.17500000E+01,
				0.80000000E+02, 0.80000000E+02, 0.80000000E+02, 0.23600000E+02 };

double __GC[] = { 0.21484993E+02, 0.39365627E+02,
                0.11733573E+02, 0.13373021E+02, 0.27289089E+02, 0.21815765E+02,
				0.28064201E+02, 0.11755576E+02, 0.21993205E+03, 0.20048616E+02,
				0.39365627E+02, 0.29287760E+02, 0.37389256E+02, 0.31843267E+02,
			    0.18399178E+02, 0.98581100E+00, 0.76616817E+01, 0.45784122E+02,
				0.98853817E+01, 0.47292509E+01, 0.61693435E+01, 0.11861769E+02 };

double __QS[] = { 0.22000000E+01, 0.20000000E+01,
                0.20000000E+01, 0.20000000E+01, 0.20000000E+01, 0.22000000E+01,
				0.20000000E+01, 0.20000000E+01, 0.30000000E+01, 0.20000000E+01,
				0.20000000E+01, 0.30000000E+01, 0.20000000E+01, 0.30000000E+01,
				0.20000000E+01, 0.50000000E+00, 0.20000000E+01, 0.20000000E+01,
				0.20000000E+01, 0.20000000E+01, 0.20000000E+01, 0.20000000E+01 };

double __SM[] = { 0.88413364E+00, 0.78492051E+00,
				0.84969723E+00, 0.84408921E+00, 0.79834670E+00, 0.88117069E+00,
				0.79869562E+00, 0.84020197E+00, 0.86620516E+00, 0.84974420E+00,
				0.78492051E+00, 0.90510607E+00, 0.77372539E+00, 0.90547490E+00,
				0.82473016E+00, 0.76297575E+00, 0.87047338E+00, 0.77146745E+00,
				0.83532584E+00, 0.86597520E+00, 0.80329478E+00, 0.84098595E+00 }; 

double __RS[] = { 0.50000000E+00, 0.50000000E+00, 
			     .50000000E+00,  .50000000E+00, .50000000E+00, .50000000E+00,
			     .50000000E+00,  .15000000E+01, .10000000E+00, .50000000E+00,
				 .50000000E+00,  .50000000E+00, .50000000E+00, .50000000E+00,
				 .50000000E+00,  .50000000E+01, .50000000E+00, .50000000E+00,
				 .50000000E+00,  .50000000E+00, .50000000E+00, .50000000E+00 };

double __A1[] = { -.10226000E+04, 0.16000900E+03,
                -.13212900E+04, 0.35834999E+02, 0.39950000E+03, -.10416100E+04,
				0.69557001E+03, 0.70087701E+03, 0.18928400E+04, -.65371002E+03,
				0.16000900E+03, 0.45810999E+03, 0.84567999E+03, -.25452299E+02,
				-.40144601E+03, 0.74053003E+03, -.29326501E+03, 0.18746400E+04,
				0.45814001E+03, -.81706201E+03, -.23182000E+04, 0.13791600E+03 };

double __A2[] = { 0.37132200E+04, -.28083600E+03,
                0.31578201E+04, -.19723599E+03, -.11932000E+04, 0.42018701E+04,
			    -.20732300E+04, -.20388700E+04, -.45849902E+04, 0.28425300E+04,
		        -.28083600E+03, -.19048800E+04, -.19349200E+04, -.31007999E+03,
				0.34389801E+03, -.16957900E+04, 0.34880301E+03, -.42837700E+04,
				-.12964900E+04, 0.67171399E+03, 0.37505901E+04, -.41737799E+03 };

double __A3[] = { -.50760498E+04, -.83797997E+02,
                -.27733201E+04, 0.17582001E+03, 0.95741998E+03, -.60544502E+04,
                0.14612500E+04, 0.13551400E+04, 0.29190200E+04, -.44766201E+04,
                -.83797997E+02, 0.16349200E+04, 0.11284700E+04, 0.97382004E+02,
                -.44822000E+03, 0.89040997E+03, -.13732100E+03, 0.23809600E+04,
                0.10074200E+04, 0.23222800E+03, -.17300400E+04, 0.35004999E+03 };

double __A4[] = { 0.26228401E+04, 0.22006700E+03,
                0.98634003E+03, -.14870000E+02, -.17536000E+03, 0.31672300E+04,
                -.85099998E+02, -.15187000E+02, -.23344000E+03, 0.26115000E+04,
                0.22006700E+03, -.19434000E+03, -.38910000E+02, 0.24988499E+03,
                0.53452899E+03, 0.70820000E+02, 0.86291000E+02, 0.33049999E+02,
                -.17992999E+03, -.96350998E+02, 0.31112000E+03, -.76129997E+02 };

double __A5[] = { -.23741000E+03, -.15442000E+02,
				-.49549999E+02, 0.45100001E+00, 0.11640000E+02, -.27304001E+03,
				0.15100000E+01, -.19600000E+01, 0.65700002E+01, -.32370001E+03,
                -.15442000E+02, 0.61900001E+01, -.31999999E+00, -.11734700E+02,
                -.28761000E+02, -.59699998E+01, -.45079999E+01, -.48800001E+01,
                0.10860000E+02, 0.94709997E+01, -.13470000E+02, 0.55419998E+01 };

double __TA[] = { 0.62708032E+00, -.33601224E+00,
                -.41694292E+00, 0.40492105E+00, -.14242107E+00, 0.11158217E+01,
				-.25900945E+00, -.29996026E+00, 0.12563297E+01, 0.88399196E+00,
			    -.33601224E+00, 0.29947156E+00, 0.12490201E+01, 0.10827350E+00,
		        -.15773018E+00, 0.13764737E+01, 0.44122058E+00, 0.17943431E+00,
                0.17006929E+01, 0.66973019E+00, 0.17068748E+00, 0.11692852E+01 };

double __EA[] = { 0.42000000E+02, 0.93000002E+01,
				0.18000000E+02, 0.14000000E+02, 0.14000000E+02 ,0.42000000E+02,
                0.20000000E+02, 0.19000000E+02, 0.30000000E+02, 0.18000000E+02,
                0.93000002E+01, 0.25000000E+02, 0.10000000E+02, 0.25000000E+02,
                0.24000000E+02, 0.12000000E+02, 0.20000000E+02, 0.25000000E+02,
                0.60000000E+02, 0.10000000E+03, 0.60000000E+02, 0.30000000E+02 };

double __GI[] = { 0.50000000E+00, 0.46000001E+00,
                0.56000000E+00, 0.50000000E+00, 0.56000000E+00, 0.47999999E+00,
				0.50000000E+00, 0.44999999E+00, 0.60000002E+00, 0.50000000E+00,
                0.46000001E+00, 0.55000001E+00, 0.50000000E+00, 0.55000001E+00,
		        0.38000000E+00, 0.50000000E+00, 0.34999999E+00, 0.55000001E+00,
			    0.50000000E+00, 0.44999999E+00, 0.50000000E+00, 0.50000000E+00 };

double __GF[] = { 0.10031745E+01, 0.10038391E+01,
                0.10456556E+01, 0.10218971E+01, 0.10197377E+01, 0.10031537E+01,
                0.10070708E+01, 0.10065769E+01, 0.10041945E+01, 0.10064120E+01,
                0.10038391E+01, 0.10030688E+01, 0.10066030E+01, 0.10016240E+01,
                0.10016652E+01, 0.10074447E+01, 0.10076891E+01, 0.10014995E+01,
                0.10049790E+01, 0.10007124E+01, 0.10020235E+01, 0.10279195E+01 };

double __CR[] = { 0.34523895E+01, 0.16488581E+01,
                0.45034027E+01, 0.29177928E+01, 0.41500096E+01, 0.34477754E+01,
                0.37076824E+01, 0.35676460E+01, 0.38610494E+01, 0.31175354E+01,
                0.16488581E+01, 0.35281672E+01, 0.17087947E+01, 0.26236989E+01,
                0.28201730E+01, 0.13170183E+01, 0.31558700E+01, 0.23904498E+01,
                0.41986074E+01, 0.61198459E+01, 0.77948456E+01, 0.39223223E+01 };

double __C1R[]= { 0.32979987E+01, 0.81879050E+00, 
                0.16976673E+02, 0.44948950E+01, 0.52029748E+01, 0.31223445E+01,
				0.52075500E+01, 0.39344752E+01, 0.39021666E+01, 0.25285969E+01,
			    0.81879050E+00, 0.54939213E+01, 0.18226281E+01, 0.20840724E+01,
				0.22750535E+01, 0.13973175E+01, 0.30427899E+01, 0.19867337E+01,
				0.78179336E+01, 0.19883865E+02, 0.40253441E+02, 0.52068253E+01 };

CEOSLomonosov::CEOSLomonosov(int _id) {
	 V0 = __V0[_id]; E0 = __E0[_id]; DX = __DX[_id]; GM = __GM[_id]; CMN = __CMN[_id]; 
	 GN = __GN[_id]; ES = __ES[_id]; GC = __GC[_id]; QS = __QS[_id]; SM = __SM[_id]; 
	 RS = __RS[_id]; A1 = __A1[_id]; A2 = __A2[_id]; A3 = __A3[_id]; A4 = __A4[_id]; 
	 A5 = __A5[_id]; TA = __TA[_id]; EA = __EA[_id]; GI = __GI[_id]; GF = __GF[_id]; 
	 CR = __CR[_id]; C1R = __C1R[_id];
}

double CEOSLomonosov::getc(double ro, double e) {
	double P = 0., C = 0, g = 0.;
	bool nonPhysical = false;
	EOSE5(ro/1000., e/1.e6, P, C, g, nonPhysical);   // C = [km/s]
	return C*1.e3;                                   // P = [Pa]
}

double CEOSLomonosov::getp(double ro, double e) {
	double _e = 1.e-6*e;                             // [J/kg] = [1e-3 kJ / 1e3 g = 1e-6 kJ/g]
	double P = 0., C = 0, g = 0.;
	bool nonPhysical = false;
	EOSE5(ro/1000., _e, P, C, g, nonPhysical);       // P = [GPa]
	return P*1.e9;                                   // P = [Pa]
}

double CEOSLomonosov::gete(double ro, double p) {
	double _e = 0., eps = .01, p_dp = 0., low_dp = 0., high_dp = 0., eMin = 1.e-6, eMax = 1.e20;	
	for(;;) {
		_e = (eMin + eMax)/2.0;
		if(fabs(eMax-eMin) < eps)
			return _e;
		p_dp = getp(ro, _e)-p;
		low_dp  = getp(ro, eMin) - p;
		high_dp = getp(ro, eMax) - p;
		if(fabs(p_dp) < eps)
			return _e;
		if(low_dp * high_dp > 0) {
			cerr << endl << "EOSLomonosov::gete(): Error in initial data -- no equation root on the interval" << endl;
			cerr << "ro=" << ro << ", p=" << p << ", F(a)=" << low_dp << ", F(b)=" << high_dp << endl;
			return -1.;
		}
		if(p_dp * low_dp > 0)
			eMin += (eMax - eMin)/2.0;
		else
			eMax -= (eMax - eMin)/2.0;
	}
}

void CEOSLomonosov::EOSE5(double ro, double e, double& P, double& C, double& g, bool& nonPhysical) {
	//int nmet = _nmat;
	const double GGI = 2.0/3.0;
	const double R13 = 1.0/3.0;
	double x = log(ro * V0);
	e = e + E0;
	double vr = exp(x);
	double v = DX * vr;
	double v0 = V0; //specific material volume at T = 0
	double px = 0;
	double gx = 0;
	double c2x = 0;
	double et = 0;
	double g1x = 0;
	if ((1.0 - v) > 0) {
		//Cold lattice contribution for liquid
		double cm = 0;
		if (x > -100.0/GM) cm = CMN*pow(v, GM);
		double cn = 0;
		if (x > -100.0/GN) cn = CMN*pow(v, GN);
		px = cm - cn;
		double ex = cm / GM - cn / GN + ES;
		c2x = px + cm * GM - cn * GN;
		px = px * v / (v0 * DX);
		double sq = GC * pow(v, QS);
		double sr = sq*pow(v/SM, RS);
		double g = sq/QS;
		double gr = sr/(QS + RS);
		gx = GGI + g - gr;
		et = e - ex;
		double g1x = sq - sr;
	} else {
		//Cold lattice contribution for solid
		double v13 = pow(v,R13);
		double v23 = v13*v13;
		double v33 = v;
		double v43 = v23*v23;
		double v53 = v23*v33;
		double v1 = A1*v13;
		double v2 = A2*v23;
		double v3 = A3*v33;
		double v4 = A4*v43;
		double v5 = A5*v53;
		double ex = 3.0 * (v1-A1 
		+ 0.5*(v2-A2) 
		+ (v3 - A3)/3.0
		+ 0.25*(v4-A4) 
		+ 0.2*(v5-A5)
		);
		double r0 = v1 + v2 + v3 + v4 + v5;
		double r1 = 4./3.*v1+5./3.*v2+2.*v3+7./3.*v4+8./3.*v5;
		double r2 = 4./9.*v1+10./9.*v2+2.*v3+28./9.*v4+40./9.*v5;
		double r3 = -8./27.*v1-10./27.*v2+28./27.*v4+80./27.*v5;
		ex = ex*V0*DX;
		px = r0*v;
		double c2x = r1*V0*DX;
		r3 = r3 + TA * ((r1*(TA + 1.0) - r2) * 3.
		   - r0 * (TA + 1.) * (TA + 2.0));
		r2 = r2 + TA * (r0*(TA + 1.0) - 2.*r1);
		r1 = r1 - r0*TA;
		double gr = r2/r1;
		gx = 0.5 * (TA + gr + GGI);
		et = e - ex;
		g1x = 0.5 * (gr*(1.-gr)+r3/r1);
	}
	// Full pressure, sound speed and g calculations ---
	double eav = EA * pow(v, GGI);
	v = V0/vr;
	g = GI + (gx - GI) * GF/(1.0+et/eav);
	double pt = g*et/v;
	P = pt + px;
	e = e - E0;
	double c2s = 0;
	// Check pressure position  ---
	if (et < 0) {
		//Non physical state
		nonPhysical = true;
		P=1.E-8; //dummy value for pressure
		C=1.E-4; //dummy value for sound speed
		std::stringstream msg;
		msg.str("");
		msg<<"Unphysical state detected while computing pressure.\n";
		msg<<"e = "<< e << "kJ/g\n";
		msg<<"ro = "<<  ro  << "g/cm^3\n";
		// -- _logger->WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::FATAL_ERROR, msg.str());
		cout << msg.str() << endl;
		// -- throw new Exception("Non physical state obtained");
		return;
	} else {
		//Correct physical state
		nonPhysical = false;
		if ((P < 1.E-4) && (vr < 0.7)) { //TO DO why hardcode?
			c2s = CR*CR*GI*(et+pt*v)/C1R;
			//c2s = sign(c2s) * sqrt(abs(c2s));
			if (c2s>=0) c2s = sqrt(fabs(c2s)); else c2s = -sqrt(fabs(c2s));
		} else {
			double e1x = GGI*EA;
			double ga = g - GF * (gx - GI)/(1.+et/eav)/(1.+et/eav)*et/eav;
			c2s = c2x 
				+ ga*P*v 
				+ g*(et-px*v) 
				+ et*GF/(1.+et/eav)*(g1x+(gx-GI)/(et+eav)*(px*v+et*e1x/eav));   









			/*double q1 = sign(c2x);
			double q2 = abs(c2x);
			double q3 = sqrt(abs(c2x));
			double q4 = sign(c2x) * sqrt(fabs(c2x));*/

			if (c2x>=0) c2x = sqrt(fabs(c2x)); else c2x = -sqrt(fabs(c2x));
			if (c2s>=0) c2s = sqrt(fabs(c2s)); else c2s = -sqrt(fabs(c2s));


			//c2x = sign(c2x) * sqrt(fabs(c2x)); 
			//c2s = sign(c2s) * sqrt(fabs(c2s));
		}
	}
	C = c2s;
	return;
}
#include"eosAnalytic.h"

#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

double EOSAnalytic::getpi(double ro, double ti) {
	double t   = ti/27.2/11605.0; double x   = 0.64/(ro/1000.0);
	double m1  = 0.00001; double m2  = 0.000004; double m3  = 0.00002; double m4  = 0.000002;
	double tau = (m1/pow(x, 1.0/3.0) + m2/pow(x, 2.0/3.0) + m3/x + m4/pow(x, 4.0/3.0)) / t;
	double tau1 =  -1.0/3.0 / t * (m1/pow(x, 4.0/3.0) + 2*m2/pow(x, 5.0/3.0) + 3*m3/(x*x) + 4*m4/pow(x, 7.0/3.0));
	// p1
	double alpha = 2.98; double beta  = 2.2; double A     = 4.483e-5;double B     = 1.378e-4;	
	double p1 = A/pow(x, alpha) - B/pow(x, beta);
	// pi
	double dzeta = 3.69; double d1    = 95.0163; double d2    = 440.63; double nc    = 0.237/111.8;
	double pi = (-dzeta*nc*t*(1.0/tau + d1 - 2*d2*tau) * tau1 + p1) * 29400.0e10;
	return pi / 10.0;
}

double EOSAnalytic::getpe(double ro, double ti, double te) {
	double pe = __pe(ro, te) - __pe(ro, ti); return pe / 10.0;
}

double EOSAnalytic::__pe(double ro, double teta) {
	double NA    = 6.0e23; double kB    = 1.38e-16; double Matom = 27.0/NA; double nat   = ro/1000.0/Matom;
	double eF    = 115000.0; double pF    = 2.0/5.0*3.*nat*kB*eF; double pCl   = 3.*nat*kB*teta;
	double pe = sqrt(pF*pF + pCl*pCl);
	return pe;
}

double EOSAnalytic::getei(double ro, double ti) {
	// tau
	double t   = ti/27.2/11605.0; double x   = 0.64/(ro/1000.0);
	double m1  = 0.00001; double m2  = 0.000004; double m3  = 0.00002; double m4  = 0.000002;
	double tau = (m1/pow(x, 1.0/3.0) + m2/pow(x, 2.0/3.0) + m3/x + m4/pow(x, 4.0/3.0)) / t;
	// en0
	double en0 = -1.05835;
	// en1
	double alpha = 2.98; double beta  = 2.2; double A     = 4.483e-5; double B     = 1.378e-4;
	double en1 = A/(alpha-1)/pow(x, alpha-1) - B/(beta-1)/pow(x, beta-1);
	// ei
	double dzeta = 3.69; double d1    = 95.0163; double d2    = 440.63; double nc    = 0.237/111.8;
	double ei = (en0 + t*dzeta*(1.0+d1*tau-2.0*d2*tau*tau) * 27.2 +	en1/nc*27.2) * 1.6e-12 * 6.0e23/27.0 * ro/1000.0;
	return ei / (ro/1000.0) * 1.0e-4;
}

double EOSAnalytic::getee(double ro, double ti, double te) {
	double ee = __ee(ro, te) - __ee(ro, ti);
	return ee / (ro/1000.0) * 1.0e-4;
}

double EOSAnalytic::__ee(double ro, double teta) {
	double NA    = 6.0e23; double kB    = 1.38e-16; double Matom = 27.0/NA; double nat   = ro/1000.0/Matom;
//	Old gamma!	
//	double gamma = 1450.0*pow(ro/2700.0, 1.0/3.0);
	double   TF  = 115000.0; double   EF  = kB * TF; double gamma = 3.14159*3.14159*3.*nat*kB*kB/2.0/EF;
	double ceCl  = 3.0/2.0*3.*nat*kB; double ee = ceCl/gamma * (sqrt(ceCl*ceCl + gamma*gamma*teta*teta) - ceCl);
	return ee;
}

double EOSAnalytic::getti(double ro, double ei) {
	double ti = solve_ti(ro, ei, 0.00001, 60000.0); 
	return ti;
}

double EOSAnalytic::gette(double ro, double ti, double ee)
{
	double NA    = 6.0e23; double kB    = 1.38e-16; double Matom = 27.0/NA; double nat   = ro/1000.0/Matom;
// Old gamma
// double gamma = 1450.0*pow(ro/2700.0, 1.0/3.0);
	double   TF  = 115000.0; double   EF  = kB * TF;
	double gamma = 3.14159*3.14159*3.*nat*kB*kB/2.0/EF; double ceCl  = 3.0/2.0*3.*nat*kB;
	double A = sqrt(ceCl*ceCl+gamma*gamma*ti*ti) + gamma/ceCl*ee*(ro/1000.0)/1.0e-4;;
	double te = 1.0/gamma * sqrt( A * A - ceCl * ceCl );
	return te;
}

double EOSAnalytic::solve_ti(double ro, double ei, double low_border, double high_border) {
	double ti; double eps = 0.01; double ti_dei, low_dei, high_dei;
	for(;;) {
		ti = (low_border + high_border)/2.0;
		if(fabs(low_border-high_border) < eps)
			return ti;
		ti_dei   = getei(ro, ti)-ei;
		low_dei  = getei(ro, low_border)-ei;
		high_dei = getei(ro, high_border)-ei;
		if(fabs(ti_dei) < eps)
			return ti;
		if(low_dei * high_dei > 0) {
			printf("\nsolve_ti: Error in initial data -- no equation root on the interval\n");
			printf("ro=%e, ei=%e, F(a)=%e, F(b)=%e\n", ro, ei, low_dei, high_dei);
			return -1.;
		}
		if(ti_dei * low_dei > 0)
			low_border  += (high_border-low_border)/2.0;
		else
			high_border -= (high_border-low_border)/2.0;
	}
}


double EOSAnalytic::solve_te(double ro, double ti, double ee, double low_border, double high_border) {
	double te; double eps = 0.01, te_dee = 0., low_dee = 0., high_dee = 0.;
	for(;;) {
		te = (low_border + high_border)/2.0;
		if(fabs(low_border-high_border) < eps)
			return te;
		te_dee   = getee(ro, ti, te)-ee;
		low_dee  = getee(ro, ti, low_border)-ee;
		high_dee = getee(ro, ti, high_border)-ee;
		if(fabs(te_dee) < eps)
			return te;
		if(low_dee * high_dee > 0) {
			printf("\nsolve_te: Error in initial data -- no equation root on the interval\n");
			printf("ro=%e, ee=%e, F(a)=%e, F(b)=%e\n", ro, ee, low_dee, high_dee);
			return -1.;
		}
		if(te_dee * low_dee > 0)
			low_border  += (high_border-low_border)/2.0;
		else
			high_border -= (high_border-low_border)/2.0;
	}
}



////////"Условная" скорость звука для случая одной температуры
double EOSAnalytic::getC(double ro, double ti, double te) {

	return 0.;
	/*
	return sqrt( getdpdro(ro, ti, te) + 
			    (getpi(ro, ti)+getpe(ro, ti, te)) * dpde(ro, ti, te) / ro / ro +
				 getpi(ro, ti) * dpdei(ro, ti, te) / ro / ro );*/
}

double EOSAnalytic::getci(double ro, double ti) {
	// return 3.2e7 * 1.0e-1;
	double NA    = 6.0e23;
	double kB    = 1.38e-16;
	double Matom = 27.0/NA;
	double nat   = ro/1000.0/Matom;
	double ci = 3.0*nat*kB;
	return ci*1.0e-1;
	// 4 testExchangeStage() : return 1.0;
}

double EOSAnalytic::getce(double ro, double te) {
	//	ce(rho, T)=( ( gam T)^( - 2 ) + cCle^( - 2) )^( - 1/2)
	double NA    = 6.0e23; double kB    = 1.38e-16; double Matom = 27.0/NA; double nat   = ro/1000.0/Matom;
	//	double gamma = 1450.0*pow(ro/2700.0, 1.0/3.0);
	double   TF  = 115000.0; double   EF  = kB * TF;
	double gamma = 3.14159*3.14159*3.*nat*kB*kB/2.0/EF;
	double ceCl  = 3.0/2.0*3.*nat*kB;
	double ce = pow(1.0/gamma/gamma/te/te + 1.0/ceCl/ceCl, -0.5);
	return ce * 1.0e-1; 
	// 4 testHeatStage()    : return 1.0;
	// 4 testExchangeStage(): return 1.0;
}

double EOSAnalytic::getkappa(double ro, double ti, double te) {
	double C      = 7.25; double beta   = 0.504; double TF     = 115000.0;
	double theta  = te/TF; double thetai = ti/TF;
	double sk54   = pow(theta*theta + 0.16, 1.25); double sk     = theta*theta + 0.44;
	double sk12   = pow(theta*theta + 0.092, 0.5);
	double kappa  = C*sk54*sk/sk12*theta / (theta*theta+beta*thetai) * 1.0e7;
	return kappa * 1.0e-5;
	// 4 heatTest(): return 1.0e-8;
}

double EOSAnalytic::getAlpha(double ro, double ti, double te) {

// Алюминий используется как стекло -- т.е. теплопроводность и обмен в нем отсутствуют
	
	return 0.;


	//return -36.0e17 * 1.0e-1 * ro/ro0 ;
	// 4 testExchangeStage() : return -1.0;
}

/*
double EOSAnalytic::getdpdro(double ro, double ti, double te)
{
	return ( getpi(ro + EOS_EPS, ti) + getpe(ro + EOS_EPS, ti, te) 
		   - getpi(ro, ti)			 - getpe(ro, ti, te) ) / EOS_EPS;
}

double EOSAnalytic::getdpdroe(double ro, double ti, double te)
{
	return ( getpi(ro + EOS_EPS, ti) + getpe(ro + EOS_EPS, ti, te) 
		   - getpi(ro, ti)			 - getpe(ro, ti, te) ) / EOS_EPS;
}

double EOSAnalytic::getdpdroei(double ro, double ti, double te)
{
	return ( getpi(ro + EOS_EPS, ti) + getpe(ro + EOS_EPS, ti, te) 
		   - getpi(ro, ti)			 - getpe(ro, ti, te) ) / EOS_EPS;
}
*/
/*

double EOSAnalytic::getdpdro_rov_roE(double ro, double ti, double te, double v)
{ return 0;}

double EOSAnalytic::getdpdrov_ro_roE(double ro, double ti, double te, double v)
{ return 0;}

double EOSAnalytic::getdpdroE_ro_rov(double ro, double ti, double te, double v)
{ return 0;}

*/
////////////// PRIVATE ////////////////

///////////////////////////
////////// УБИТЬ ЛИШНЕЕ !!!
///////////////////////////

double EOSAnalytic::dpdei(double ro, double ti, double te)
{
	return dpdti(ro, ti, te) / deidti(ro, ti, te) - 
		   dpdte(ro, ti, te) / dedte(ro, ti, te);
}


double EOSAnalytic::dpde(double ro, double ti, double te)
{
	return dpdte(ro, ti, te) / dedte(ro, ti, te);
}


double EOSAnalytic::dpdti(double ro, double ti, double te)
{
	return ( getpi(ro, ti + EOS_EPS) + getpe(ro, EOS_EPS, te) 
		   - getpi(ro, ti)           - getpe(ro, ti, te) ) / EOS_EPS;
}


double EOSAnalytic::dedti(double ro, double ti, double te)
{
	return ( getei(ro, ti + EOS_EPS) + getee(ro, ti + EOS_EPS, te) 
		   - getei(ro, ti)           - getee(ro, ti, te) ) / EOS_EPS;
}


double EOSAnalytic::deidti(double ro, double ti, double te)
{
	return ( getei(ro, ti + EOS_EPS) - getei(ro, ti) ) / EOS_EPS;
}


double EOSAnalytic::dpdte(double ro, double ti, double te)
{
	return ( getpe(ro, ti, te + EOS_EPS) - getpe(ro, ti, te) ) / EOS_EPS;

}


double EOSAnalytic::dedte(double ro, double ti, double te)
{
	return ( getee(ro, ti, te + EOS_EPS) - getee(ro, ti, te) ) / EOS_EPS;
}


EOSAnalytic::EOSAnalytic()
{	
	MAX_T=100000.0; 
	MIN_T=300.0; 
	MAX_RO=5000.0; 
	MIN_RO=0.001;
	ro0 = 2700.;
}

double EOSAnalytic::getphase(double ro, double ti)
{
	return 1.0;
}


double EOSAnalytic::getmix(double ro, double ti)
{
	return 1.0;
}

double  EOSAnalytic::getnuWR(double ro, double ti, double te, double b)
{
	return 0.0;
}

double EOSAnalytic::getEntropy(double ro, double ti, double te)
{
	return 0.0;
}

double EOSAnalytic::getGamma(void) {
	return 0.0;
}
/*
complex<double> EOSAnalytic::geteps(double ro, double ti, double te, double Z=3.0)
{
	return complex<double>(0.0);
}

double EOSAnalytic::getdpdt(double ro, double ti, double te)
{
	return (getpi(ro, ti+EOS_EPS)-getpi(ro, ti))/EOS_EPS;
}

double EOSAnalytic::getdedt(double ro, double ti, double te)
{
	return (getei(ro, ti+EOS_EPS)-getei(ro, ti))/EOS_EPS;
}*/

/////////////////////////
//// EOSPyrexGlass()
/////////////////////////


double EOSPyrexGlass::getpi(double ro, double ti) {
	double _pi = (ro/1000. - 2.23)*(10.3 + 5.047*(ro/1000. - 3.52)*(ro/1000. - 3.52))*exp(0.13*ro/1000.) * 1.e9; // [Pa]
	return _pi;
}

double EOSPyrexGlass::getei(double ro, double ti) {
	double _ei = 830. * ti;  //[J/kg]
	return _ei;
}

double EOSPyrexGlass::getti(double ro, double ei) {
	double _ti = ei / 830.;
	return _ti;
}


// Mie-Gruneisen EOS

double EOSMieGruneisenRu::getpi(double ro, double ti) {
	//p[x_,TiK_]=( A/v0*x*(x^a - x^b) + 3/v0*Ti*x*G1[x] )*29400 
	double rof = 12410.;
	double x = ro/rof;
	double tiEnergyUnits = ti/11605./27.2;   
	return (A/v0*x*(pow(x,a) - pow(x,b)) + 3./v0*tiEnergyUnits*x*G(x))*29400.*1.e9; // Что такое 29400 и 29.2??? [Pa]
}

double EOSMieGruneisenRu::getei(double ro, double ti) {
	//en[x_,TiK_]=( A/a*(  x^a - a/b*(x^b)  ) + 3*Ti )*27.2*0.9533       ((1Е)
	//(*внутренняя энергия в kJ/g при температуре Ti в Кельвинах и сжатии  . Т.е. теплоемкость принята равной 3*кБ.    *)
	double rof = 12410.;
	double x = ro/rof;
	double tiEnergyUnits = ti/11605./27.2;   
	return ( A/a*(pow(x,a) - a/b*pow(x,b) ) + 3*tiEnergyUnits )*27.2*0.9533*1.e6; // [J/kg]
}

double EOSMieGruneisenRu::getpe(double ro, double ti, double te) { 
	//Pe[Te_] = 0.00127947*Te + 1.696*10^-7*Te^2 - 4.12488*10^-12 *Te^3 + 3.41851*10^-17*Te^4
	//return (0.00127947*te + 1.696e-7*te*te - 4.12488e-12*te*te*te + 3.41851e-17*te*te*te*te)*1.e9; 
	//Pe[ x_, Te_ ] = xi ( 0.00127947*Te + 1.696*10^-7*Te^2 - 4.12488*10^-12 *Te^3 + 3.41851*10^-17*Te^4), 
	//где xi = ( 1 + a )*( 1/x + a/x^0.2) )^(-1), здесь a = 5.154
	/*
	P1 = 0.00127947; P2 = 1.696E-7;
    P3 = -4.12488E-12; P4 = 3.41851E-17;
    E2 = 2.11384E-7; E3 = -4.45326E-12;
    E4 = 4.20893E-17;

    (Pe = Te*(P1 + Te*(P2 + Te*(P3 + Te*P4)))*f1(x), GPa; Ee = Te*Te*(E2 + Te*(E3 + Te*E4))/rho*f1(x), kJ/g;)
	*/


	double rof = 12410.;
	double x = ro/rof;
	double af = 5.154;
	double xi = (1. + af)/(1./x + af/pow(x, .2));
	return xi*(0.00127947*te + 1.696e-7*te*te - 4.12488e-12*te*te*te + 3.41851e-17*te*te*te*te)*1.e9; // [Pa]
}


double EOSMieGruneisenRu::getee(double ro, double ti, double te) { 
	//Ee[Te_] = 2.11384*10^-7*Te^2 - 4.45326*10^-12*Te^3 + 4.20893*10^-17*Te^4
	double rof = 12410.;
	double x = ro/rof;
	double af = 5.154;
	double xi = (1. + af)/(1./x + af/pow(x, .2));
	return (2.11384e-7*te*te - 4.45326e-12*te*te*te + 4.20893e-17*te*te*te*te)*1.e9/ro; // [J/kg]  -- !!! Возможно, ошибка тут, надо домножать на 10^6, 
																						// а не на 10^9? или эти 10^3 наигрываются при делении на 
																					    // плотность в граммах у В.А?
}

double EOSMieGruneisenRu::getAlpha(double ro, double ti, double te) {
	//18 – 12.5 Te/(50000 + Te)
	return -(18. - 12.5*te/(50000. + te))*1.e17; 
}

double EOSMieGruneisenRu::getkappa(double ro, double ti, double te) {
	//kappas-e = x^2 ( 0.8e6 / Te + pow( Te, 1.76 ) / 3.13e5 ). 
	//kappas-i  =  x^2 60*Te/Ti W/m/K
	//kappa = ( 1/kappas-e + 1/kappas-i )^-1
	//kappaSE = 0.8 10^6/Te + Te^1.76/(3.13 10^5);   kappaSI = 60 Te/Ti;    kappa[Te_,Ti_] = ( 1/kappaSE + 1/kappaSI )^-1;


	/*     
Теплопроводность
    krt = 117;
    Tr = 298;

    a = 1.6569036734344376; b = 1.1736399604338026;
    am2 = 0.27046824537560393; am = 0.1991637915011424;
    bm = 1.937141548654876;
    a0 = 25.123052588773383; a1 = 0.2524851528775459;
    b1 = 0.40172793052485617; b2 = 1.787755373843094;
    b3 = 0.3725238860898682;

	***

	В другом варианте имеем
t = Te*6/(   9.9*11605*x^(2/3)   );   a0=25.123;  a1=0.252485;  b1=0.4017279;  b2=1.787755;  b3=0.37252; 
kappaSE = 10^3*x*(1+ b1*t+b2*t^2+b3*t^3) / t / (a1*t+a0)       (*Ю*)

t = Te*6/(   8*11605*x^(2/3)   ); am2=0.270468;     am=0.19916379;       bm=1.9371415;
cv[t1_] =  t1*(1.+am2*t1*t1)/(1.+am*(t1^bm)) ;
a=1.44673;    b=1.007457;     cab=(a-b)/(b+1);      y[x_] = (1+cab)*(x^(2*a+1.))/(1+cab*(x^(a+1.)))

Trt=298;     xrt=12.407/12.47;            trt=Trt*6/(   8*11605.*xrt^(2/3)   );         krt=117.
kappaSI = krt*(x/xrt)*(y[x]/y[xrt])*Trt/Ti*cv[t]/cv[trt]
kappa[ Te_, Ti_, x_ ] = ( 1/kappaSE + 1/kappaSI )^-1;  


*/


	double rof = 12410.;
	double x = ro/rof;
/*
	Вычисление каппа_с-е по методу Ю.В. Петрова.
t = Te*6/(   8*11605*x^(2/3)   )    (*вопрос 1 Юре  исправлено:  почему здесь ЕФ при норм.условиях равно 9.9 эВ ? Должно быть 8 эВ*)
 (* keem - это не теплопроводность! Не каппа, а тепл.сопротивление=величина обратная к каппа.*)
(*коэфф для s-e от 19.5.2017*)
keem[Te_,x_]=(10^(-3))/x*t*(a1*t+a0)/(1+ b1*t+b2*t^2+b3*t^3)    (*Ю*) !!!!              = 1/kappase*/
	double a0 = 25.123052588773383;
	double a1 = 0.2524851528775459; 
	double b1 = 0.40172793052485617;
	double b2=1.787755373843094;
    double b3=0.3725238860898682;
	double t = te*6./(8.*11605.*pow(x, 2./3.));
	double kappase = 1./ ( 1.e-3/x*t*(a1*t+a0)/(1.+b1*t+b2*t*t+b3*t*t*t) );
	double am2=0.27046824537560393, am=0.1991637915011424, bm=1.937141548654876;
	// cv[t1_] =  t1*(1.+am2*t1*t1)/(1.+am*(t1^bm)) ; t1=t!!!
	double krt = 117.;
	double Trt = 298.;
	double xrt=12.407/12.47;
	double trt = Trt*6./(8.*11605.*pow(xrt, 2./3.));
	double cvt = t*(1.+am2*t*t)/(1.+am*(pow(t, bm)));
	double cvtrt = trt*(1.+am2*trt*trt)/(1.+am*(pow(trt, bm)));
	double cab = (a-b)/(b+1);
	//y[x_]=(1+cab)*(x^(2*a+1.))/(1+cab*(x^(a+1.)))
	double y = (1+cab)*(pow(x, 2*a+1.))/(1+cab*(pow(x, a+1.)));
	double yxrt = (1+cab)*(pow(xrt, 2*a+1.))/(1+cab*(pow(xrt, a+1.)));
	/*	Trt=298;     xrt=12.407/12.47;            trt=Trt*6/(   8*11605.*xrt^(2/3)   );         krt=117.
	(* 117 Вт/м/К – это табличное значение коэффициента теплопроводности при комнатных условиях*)
	ksei[Te_,Ti_,x_] = krt*(x/xrt)*(y[x]/y[xrt])*Trt/Ti*cv[t]/cv[trt]*/
	double kappasi = krt * (x/xrt) * (y/yxrt) * Trt/ti * cvt/cvtrt;
	return 1./(1./kappase + 1./kappasi); // [W/m/K]
}

double EOSMieGruneisenRu::getC(double ro, double ti, double te) { 
/*	double x = ro/ro0;
	double c = c1;
	double s = 5970 * sqrt( (c+1.)*pow(x, 2.*a+1.)/(1. + c*pow(x,a+1.)));
	return s;*/
	return 20.e4;
}

double EOSMieGruneisenRu::getci(double ro, double ti) { 
	//return 25.*1.e5; 
	double _ci = 3.*.9533*1.e6/11605.*ro; return _ci; // [J/m3/K]
	//return 3.*8.31*ro/101.07e-3;  //[J/m3/K]
}

double EOSMieGruneisenRu::getce(double ro, double te) {
	//Ce[Te_ ] = 0.00422767 Te - 1.33598*10^-7*Te^2 + 1.68357*10^-12*Te^3
	//return (0.00422767*te - 1.33598e-7*te*te + 1.68357e-12*te*te*te)*1.e5; //[J/kg/K]
	// посмотреть систему и единицы
	// в pe исправлено, а в ce осталось то что было
	// xi(0.00422767 Te - 1.33598*10^-7*Te^2 + 1.68357*10^-12*Te^3)   (форма маth)
	// в единицах (10^5 J/m^3/K).
	double x = ro/ro0;
	double a = 5.154;
	double xi = (1. + a)/(1./x + a/pow(x, .2));
	return xi*(0.00422767*te - 1.33598e-7*te*te + 1.68357e-12*te*te*te) * 1.e5; // [J/m3/K]
}

double EOSMieGruneisenRu::getphase(double ro, double ti) {return 0.; }
double EOSMieGruneisenRu::getmix(double ro, double ti) {return 0.; }

double EOSMieGruneisenRu::getti(double ro, double ei) {
	double ti = solve_ti(ro, ei, 0.001, 100000.0); 
	return ti;
}

double EOSMieGruneisenRu::gette(double ro, double ti, double ee) {
	double te = solve_te(ro, ti, ee, 0.001, 100000.0); 
	return te; 
}
	
double EOSMieGruneisenRu::solve_ti(double ro, double ei, double low_border, double high_border) {
	double ti; double eps = 0.01; double ti_dei, low_dei, high_dei;
	for(;;) {
		ti = (low_border + high_border)/2.0;
		if(fabs(low_border-high_border) < eps)
			return ti;
		ti_dei   = getei(ro, ti)-ei;
		low_dei  = getei(ro, low_border)-ei;
		high_dei = getei(ro, high_border)-ei;
		if(fabs(ti_dei) < eps)
			return ti;
		if(low_dei * high_dei > 0) {
			printf("\nsolve_ti: Error in initial data -- no equation root on the interval\n");
			printf("ro=%e, ei=%e, F(a)=%e, F(b)=%e\n", ro, ei, low_dei, high_dei);
			return -1.;
		}
		if(ti_dei * low_dei > 0)
			low_border  += (high_border-low_border)/2.0;
		else
			high_border -= (high_border-low_border)/2.0;
	}
}

double EOSMieGruneisenRu::solve_te(double ro, double ti, double ee, double low_border, double high_border) {
	double te; double eps = 0.01, te_dee = 0., low_dee = 0., high_dee = 0.;
	for(;;) {
		te = (low_border + high_border)/2.0;
		if(fabs(low_border-high_border) < eps)
			return te;
		te_dee   = getee(ro, ti, te)-ee;
		low_dee  = getee(ro, ti, low_border)-ee;
		high_dee = getee(ro, ti, high_border)-ee;
		if(fabs(te_dee) < eps)
			return te;
		if(low_dee * high_dee > 0) {
			printf("\nsolve_te: Error in initial data -- no equation root on the interval\n");
			printf("ro=%e, ee=%e, F(a)=%e, F(b)=%e\n", ro, ee, low_dee, high_dee);
			return -1.;
		}
		if(te_dee * low_dee > 0)
			low_border  += (high_border-low_border)/2.0;
		else
			high_border -= (high_border-low_border)/2.0;
	}
}


double EOSMieGruneisenRu::getGamma(void) { return 0.; }
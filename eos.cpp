#include <math.h>
#include "eos.h"

/*double EOS::getC_an(double ro, double ti, double te)
{
	double	
			   pi =        getpi(ro, ti),
			pi_ro = getdpidro_ti(ro ,ti),
			pi_ti = getdpidti_ro(ro ,ti),
			ei_ro = getdeidro_ti(ro ,ti),
			ei_ti = getdeidti_ro(ro ,ti);
	
	//////////////DEBUG/////////////////
	double      c = sqrt(fabs((pi_ro - pi_ti * (ei_ro - pi/ro/ro) / ei_ti)));
	
	
	double		c2=sqrt(5.3 * getpi(ro, ti)/ro);

	////////////////////////////////////

	 

	return c;
}

double EOS::getdpidro_ti(double ro, double ti)
{
	double h = 0.1;
	
	return ( getpi(ro*(1.0+h), ti)-getpi(ro, ti) ) / (ro*h);
}

double EOS::getdpidti_ro(double ro, double ti)
{
	double h = 0.1;
	return ( getpi(ro, ti*(1.0+h))-getpi(ro, ti) ) / (ti*h);
}

double EOS::getdeidro_ti(double ro, double ti)
{
	double h = 0.1;
	
	return ( getei(ro*(1.0+h), ti)-getei(ro, ti) ) / (ro*h);
}

double EOS::getdeidti_ro(double ro, double ti)
{
	double h = 0.1;
	return ( getei(ro, ti*(1.0+h))-getei(ro, ti) ) / (ti*h);
}
*/
double EOS::getdeedte_ro(double ro, double ti, double te)
{
	double h = 0.1;
	return ( getee(ro, ti, te*(1.0+h))-getee(ro, ti, te) ) / (te*h);
}


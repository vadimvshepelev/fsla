
#include "methodold.h"

#include <math.h>
#include <stdio.h>


CMethodOld::CMethodOld(EOSOld *_eos) : eos(*_eos),
							  La(Vector4::ZERO),
							  Lb(Vector4::ZERO),
							  Lg(Vector4::ZERO),
							  Ld(Vector4::ZERO)
{}


void CMethodOld::averageNode(Node &n1, Node &n2, Node &nav)
{
	nav.ro = (n1.ro + n2.ro) / 2;
	
	//Скорость в "среднем" узле есть (i+1)-ая скорость, усреднять ее не нужно.
	//Старое:	nav.v  = (n1.ro*n1.v + n2.ro*n2.v) / 2 / nav.ro;
	//Ооо! Я понял эту тему. Сейчас сделаю в CMethodLagrange наследника функции 
	//averageNode c этим изменением, тогда все будет хорошо. Здесь оставляю то,
	//что было, коммент не стираю.

	nav.v  = nav.ro ? (n1.ro*n1.v + n2.ro*n2.v) / 2 / nav.ro : 0;
	nav.e  = nav.ro ? (n1.ro*n1.e + n2.ro*n2.e) / 2 / nav.ro : 0;
	nav.ei = nav.ro ? (n1.ro*n1.ei + n2.ro*n2.ei) / 2 / nav.ro : 0;
	
	updateNode(nav);

	//////////////////////DEBUG///////////////
	// Возможно, энтропия на границе падает потому, что мы неверно усредняем скорость C

	// if (n1.ro==0.) 
	//	nav.C = (n1.C + n2.C)/2.;

	//////////////////////////////////////////

}


void CMethodOld::updateNode(Node &n)
{
	if(eos.getType() == ideal
		/////////// DEBUG!!! ///////////////
		
		|| eos.getType() == table

		////////////////////////////////////
		) /* или если это 1T-задача -- временно добавляем этот пункт!*/
	{
		n.ee = 0.0;
		n.ti = eos.getti(n.ro, n.e);
	}
	else
	{
		n.ee = n.e - n.ei;
		n.ti = eos.getti(n.ro, n.ei); 
	}

	n.te = eos.gette(n.ro, n.ti, n.ee);

	n.p  = eos.getp (n.ro, n.ti, n.te);
	//n.pi = eos.getpi(n.ro, n.ti);
	//n.pe = eos.getpe(n.ro, n.ti, n.te);
	
	//n.ci = eos.getci(n.ro, n.te);
	//n.ce = eos.getce(n.ro, n.te);

	n.C		= eos.getC (n.ro, n.ti, n.te);
	//n.Alphaei     = eos.getAlpha(n.ro, n.ti, n.te);
	//n.kappa = eos.getkappa(n.ro, n.ti, n.te);
}


void CMethodOld::fillLambda(Vector4 &Fm1, Vector4 &F0, Vector4 &Fp1, Vector4 &Fp2,
						 Vector4 &L, double step)
{
	double curant;
	double criteria;
  
	for(int i=0; i<4; i++)
	{
		curant   = fabs(L[i]) * step;

//////////////// DEBUG ////////////////////
if(curant < 0 || curant > 1)
{
	printf("Curant achtung!! %lf\n", curant);
}
//////////////// DEBUG ////////////////////

		criteria = L[i] * (fabs(Fp2[i]-Fp1[i]) - fabs(F0[i]-Fm1[i]));
		fillLambdaComponent(i, L[i], criteria, curant);
	}
}


void CMethodOld::fillLambdaComponent(int i, double lambda, double criteria, double curant)
{
	double alpha, beta, gamma, delta;

	if( criteria >= 0 && lambda >= 0 )
	{
		alpha = -0.5 * (1.0-curant);
		 beta =  1.0 - alpha;
		gamma =  0.0; 
		delta =  0.0;
	}
	else if( criteria >= 0 && lambda < 0 )
	{
		alpha =  0.0;
		 beta =  0.0;
		delta = -0.5 * (1.0-curant);
		gamma =  1.0 - delta; 
	}
	else if( criteria < 0 && lambda >= 0 )
	{
		alpha =  0.0;
		gamma =  0.5 * (1.0-curant);
		 beta =  1.0 - gamma;
		delta =  0.0;
	}
	else
	{
		alpha =  0.0;
		 beta =  0.5 * (1.0-curant);
		gamma =  1.0 - beta;
		delta =  0.0;
	}

	La[i] = alpha*lambda;
	Lb[i] =  beta*lambda;
	Lg[i] = gamma*lambda;
	Ld[i] = delta*lambda;
}
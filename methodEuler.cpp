
#include "methodEuler.h"
#include <math.h>

void CMethodEuler::createGrid(CFieldOld &ms){
	vGrid = new double[ms.getSize()+1];
	X     = new double[ms.getSize()+1];
	for(int i=0; i<=ms.getSize(); i++){
		    X[i] = ms[i].x;
		vGrid[i] = 0.;
	}
}

void CMethodEuler::deleteGrid(){
	delete []X;
	delete []vGrid;
}

CMethodEuler::~CMethodEuler(){
	deleteGrid();
}

void CMethodEuler::matter2Flow(CFieldOld &ms) {
	double E=0.;
	for(int i=0; i<ms.getSize(); i++) {
		Node &n = ms[i];
		E = n.e + 0.5*n.v*n.v;
		n.W = Vector4(n.ro, n.ro*n.v, n.ro*E, 0);
		n.F = getOmegaInv(n) * n.W;
	}
}

void CMethodEuler::advanceFlow(CFieldOld &ms, double tau) {
	Matrix4 O, OI;
	Vector4 Fp, Fm;
	Vector4 L;
	Node	nav;
	int i=0, j=0;
	////////////////////DEBUG////////////////
    // Граничные условия уже заданы в функции CSolver::initVars() в самом конце с помощью функции setEdge().
	// Пока еще не разобрались, где лучше ставить граничные условия. Пусть будет тут.
	ms.setEdgeTransparent();
	/////////////////////////////////////////
	for(i=0; i<ms.getSize(); i++) {
		// F(i+1/2) flow
		j = i;
		averageNode(ms[j], ms[j+1], nav);
		O  = getOmega(nav);
		OI = getOmegaInv(nav);
		L  = getLambda(nav);
		fillLambda(ms[j-1].F, ms[j].F, ms[j+1].F, ms[j+2].F,
				   L, tau/(ms[j+1].x - ms[j].x));
		Fp = Vector4::ZERO;
		Fp += O * (La * (OI * ms[j-1].W));
		Fp += O * (Lb * (OI * ms[j  ].W));
		Fp += O * (Lg * (OI * ms[j+1].W));
		Fp += O * (Ld * (OI * ms[j+2].W));
		// F(i-1/2) flow
		j = i-1;
		averageNode(ms[j], ms[j+1], nav);
		O  = getOmega(nav);
		OI = getOmegaInv(nav);
		L  = getLambda(nav);
		fillLambda(ms[j-1].F, ms[j].F, ms[j+1].F, ms[j+2].F,
				   L, tau/(ms[j+1].x - ms[j].x));
		Fm = Vector4::ZERO;
		Fm += O * (La * (OI * ms[j-1].W));
		Fm += O * (Lb * (OI * ms[j  ].W));
		Fm += O * (Lg * (OI * ms[j+1].W));
		Fm += O * (Ld * (OI * ms[j+2].W));
		// Calculating W_temp 
		ms[i].W_temp = ms[i].W - (Fp-Fm) * tau/(ms[i+1].x - ms[i].x);
	}
	for(i=0; i<ms.getSize(); i++) {
		ms[i].W = ms[i].W_temp;
	}
}

void CMethodEuler::advanceFlowVacuum(CFieldOld &ms, double tau) {
	Matrix4 O, OI;
	Vector4 Fp, Fm;
	Vector4 L;
	Node	nav;
	int i=0, j=0, nSize=ms.getSize(); //N это количество узлов сетки (и количество промежутков плюс один)
	///DEBUG/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	//		(Дальше план такой
	//		+ 1. Делаем сетку новых сеточных скоростей. Координаты в конце каждого этапа сдвигаем на vGrid*tau.
	//		+ 2. В сеточных координатах пересчитываем значения всех функций (т.е. "пересаживаем" функции на новую сетку) (проигнорил)
	//		+ 3. Модернизируем схему (добавляем потоки, связанные с движением сетки).
	//		+ 4. Получаем результат, хорошо бы сравнить его с аналитическим.
	//		Не работает.)
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Алгоритм на подвижной сетке (какой сейчас есть):
	// 1. method.matter2Flow(ms); -- заполняет все n.W и n.F = OmegaInv * n.W (в ячейках)
	// 2. method.advanceFlowVacuum(ms, tau);
	//	  - говорит считать поток F_i-1/2 на левой границе нулевым
	//    - считает собственно схему (все матрицы, обратные матрицы и собственные значения в узлах, затем потоки и,
	//                                собственно, консервативные переменные на новом слое.) Собственные значения модернизированы!!!
	// 3. method.flow2MatterVacuum(ms, ms_temp, tau); 
	//    - из консервативных переменных получает ro, v, e и все остальное.
	//    - считает новые скорости сетки
	//	  - считает новое положение сетки (не сдвигает!) 
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Поехали
	for(i=0; i<ms.getSize(); i++)
	{
		// F(i+1/2) flow
		j = i;
		averageNode(ms[j], ms[j+1], nav);
		O  = getOmega(nav);
		OI = getOmegaInv(nav);
		L  = getLambda(nav);
		//////////////////////// DEBUG ///////////////////////////////
		/// Модифицируем собственные значения, вычитаем скорость сетки
		L[0]-=vGrid[i];
		L[1]-=vGrid[i];
		L[2]-=vGrid[i];
		L[3]-=vGrid[i];
		//////////////////////////////////////////////////////////////
		fillLambda(ms[j-1].F, ms[j].F, ms[j+1].F, ms[j+2].F,
				   L, tau/(ms[j+1].x - ms[j].x));
		Fp = Vector4::ZERO;
		Fp += O * (La * (OI * ms[j-1].W));
		Fp += O * (Lb * (OI * ms[j  ].W));
		Fp += O * (Lg * (OI * ms[j+1].W));
		Fp += O * (Ld * (OI * ms[j+2].W));
		// F(i-1/2) flow
		j = i-1;
		averageNode(ms[j], ms[j+1], nav);
		O  = getOmega(nav);
		OI = getOmegaInv(nav);
		L  = getLambda(nav);
		//////////////////////// DEBUG ///////////////////////////////
		/// Модифицируем собственные значения, вычитаем скорость сетки
		L[0]-=vGrid[i];
		L[1]-=vGrid[i];
		L[2]-=vGrid[i];
		L[3]-=vGrid[i];
		//////////////////////////////////////////////////////////////
		fillLambda(ms[j-1].F, ms[j].F, ms[j+1].F, ms[j+2].F,
				   L, tau/(ms[j+1].x - ms[j].x));
		Fm = Vector4::ZERO;
		Fm += O * (La * (OI * ms[j-1].W));
		Fm += O * (Lb * (OI * ms[j  ].W));
		Fm += O * (Lg * (OI * ms[j+1].W));
		Fm += O * (Ld * (OI * ms[j+2].W));
		// DEBUG ////////////////
		// Исправление номер два: если ячейка нулевая, то поток Fm -- нулевой
		if(i==0) {
			Fm[0] = 0.; Fm[1] = 0.; Fm[2] = 0.; Fm[3] = 0.;
		}
		/////////////////////////
		// Calculating W_temp 
		ms[i].W_temp = ms[i].W - (Fp-Fm) * tau/(ms[i+1].x - ms[i].x);
 		Node &n = ms[i];
		double W0 = n.W[0], W1 = n.W[1], W2 = n.W[2];
		double WT0 = n.W_temp[0], WT1 = n.W_temp[1], WT2 = n.W_temp[2];
		double xxx = 0.;
	}
	for(i=0; i<ms.getSize(); i++) {
		ms[i].W = ms[i].W_temp;
	}
}

void CMethodEuler::flow2Matter(CFieldOld &ms, double tau)
{
	double E=0.;
	for(int i=0; i<ms.getSize(); i++)
	{
		Node &n = ms[i];
		n.ro = n.W[0];
		n.v  = (n.ro!=0.0) ? n.W[1] / n.ro : 0.0;
		  E  = (n.ro!=0.0) ? n.W[2] / n.ro : 0.0;
		n.e  = E - 0.5*n.v*n.v;
		n.ei = 0;		
		updateNode(n);
		//////////////
		double q=0.;
		/////////////
	}
}

void CMethodEuler::flow2MatterVacuum(CFieldOld &ms, CFieldOld &ms_temp, double t, double tau)
{
	double E=0.;
	for(int i=0; i<ms.getSize(); i++) {
		Node &n = ms[i];
		n.ro = n.W[0];
		n.v  = (n.ro!=0.0) ? n.W[1] / n.ro : 0.0;
		  E  = (n.ro!=0.0) ? n.W[2] / n.ro : 0.0;
		n.e  = E - 0.5*n.v*n.v;
		n.ei = n.e;
		n.ee = 0.;
		updateNode(n);
		//////////////
		double q=0.;
		/////////////
	}
	///////////////////////    DEBUG ///////////////////////////////////////////////////////////////////////////////////////////////
	int nSize = ms.getSize();
	// Задаем скорости сетки на этапе. Условия: 1. Сетка всегда равномерная. 2. Скорость левого конца равна скорости v[0].
	for(int i=0; i < ms.getSize(); i++)
	{
		int j = nSize - 1 - i;
		//DEBUG//
		//Исправление номер один: скорость вычисляем из волны разрежения (вопрос -- всей волны или между последней ячейкой и вакуумом?) 
		//double _v0 = ms[0].v;
			double _v0 = -2.*ms[0].C/(5./3.-1.);
			////////////////
			// Здесь казус №1: расчетная скорость получается по модулю больше, чем аналитическая, а они должны быть одинаковы в силу того
			// способа, которым мы считаем скорость (скорость левого узла сетки равна аналитической) 
			// Казус №2 -- надо доразобраться с сеткой, что в целых точках, что в полуцелых. И будет, я думаю, всё ОК.
			////////////////
		////////////////////////
		double _vg = ms[0].v*j / (nSize-1);
		vGrid[i]= _vg;
	}
	double v0 = vGrid[0], v1=vGrid[1], v2=vGrid[2];
	// Сдвигаем сетку
	X[0] = ms[0].x + ms[0].v*tau;
	for(int i=1; i<ms.getSize(); i++){
		//X[i] = ms[i].x + vGrid[i]*tau;
		X[i] = ms[0].x + (ms[nSize].x - ms[0].x)/(nSize) * i;
	}
	// Значения сеточных функций меняются от сдвига сетки по профилю. Получаем новые значения интерполяцией.
	/*FILE *f = fopen("grid0.dat", "w");
	for(int i=0; i<nSize; i++)
		fprintf(f, "%e %e\n", ms[i].x, ms[i].ro);
	fclose(f);*/
				///////////////////////////////////
				/// Вот тут я что-то насиропил. Не работает ничего. Исправить.
				/// Прямо графически можно посмотреть, куда уезжает профиль.
		   	/*	ms_temp[0].ro = ms[0].ro; ms_temp[nSize-1].ro = ms[nSize-1].ro;
				ms_temp[0].v  = ms[0].v ; ms_temp[nSize-1].v  = ms[nSize-1].v;
				ms_temp[0].ei = ms[0].ei; ms_temp[nSize-1].ei = ms[nSize-1].ei;
				ms_temp[0].ee = ms[0].ee; ms_temp[nSize-1].ee = ms[nSize-1].ee;
				for(int i=1; i<nSize-1; i++){
					double alpha = (X[i]-ms[i-1].x)/(ms[i].x-ms[i-1].x);
					ms_temp[i].ro = ms[i-1].ro + (ms[i].ro - ms[i-1].ro)*alpha;
					ms_temp[i].v  = ms[i-1].v  + (ms[i].v  - ms[i-1].v )*alpha;
					ms_temp[i].ei = ms[i-1].ei + (ms[i].ei - ms[i-1].ei)*alpha;
					ms_temp[i].ee = ms[i-1].ee + (ms[i].ee - ms[i-1].ee)*alpha;
				}
				for(int i=0; i<nSize; i++){
					ms[i].ro = ms_temp[i].ro;
					ms[i].v  = ms_temp[i].v;
					ms[i].ei = ms_temp[i].ei;
					ms[i].ee = ms_temp[i].ee; 
				}*/
	///
	// В общем, если применить этот код, волна расползается по ширине примерно в два раза. Наверно, он лишний.
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	for(int i=0; i<ms.getSize(); i++){
		ms[i].x  = X[i];
				Node &n = ms[i];
				updateNode(n);
	}
	///////////////DEBUG////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/*f = fopen("grid1.dat", "w");
	for(int i=0; i<nSize; i++)
		fprintf(f, "%e %e\n", ms[i].x, ms[i].ro);
	fclose(f);*/
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}

Matrix4 CMethodEuler::getOmega(Node &n)
{
	double dpdroE = getdpdroE_ro_rov(n.ro, n.ti, n.te, n.v);
	double E = n.e + 0.5*n.v*n.v;
	double Epro = E + n.p/n.ro;

	return Matrix4(
		1.0,			1.0,				   1.0,			   0.0,
		n.v - n.C,		n.v,				   n.v + n.C,	   0.0,
		Epro - n.C*n.v, Epro - n.C*n.C/dpdroE, Epro + n.C*n.v, 0.0,
		0.0,			0.0,				   0.0,			   1.0
	);
}

Matrix4 CMethodEuler::getOmegaInv(Node &n)
{
	if(n.C!=0.){
		double dpdro  = getdpdro_rov_roE(n.ro, n.ti, n.te, n.v);
		double dpdrov = getdpdrov_ro_roE(n.ro, n.ti, n.te, n.v);
		double dpdroE = getdpdroE_ro_rov(n.ro, n.ti, n.te, n.v);
		double cc = 1.0/(2.0*n.C*n.C);
		return Matrix4(
			cc*(dpdro + n.v*n.C),	  cc*(dpdrov - n.C), cc*(dpdroE),	   0.0,
			cc*(2.0*(n.C*n.C-dpdro)), cc*(-2.0*dpdrov),  cc*(-2.0*dpdroE), 0.0,
			cc*(dpdro - n.v*n.C),	  cc*(dpdrov + n.C), cc*(dpdroE),	   0.0,
			0.0,					  0.0,				 0.0,			   1.0
		);
	} else // in case of vacuum
		return Matrix4(
			0., 0., 0., 0.,
			0., 0., 0., 0.,
			0., 0., 0., 0.,
			0., 0., 0., 0.
		);
}


Vector4 CMethodEuler::getLambda(Node &n)
{
	return Vector4(n.v-n.C, n.v, n.v+n.C, 0);
}

void CMethodEuler::averageNode(Node &n1, Node &n2, Node &nav)
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


void CMethodEuler::updateNode(Node &n)
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
	n.pi = eos.getpi(n.ro, n.ti);
	n.pe = eos.getpe(n.ro, n.ti, n.te);
	
	n.ci = eos.getci(n.ro, n.te);
	n.ce = eos.getce(n.ro, n.te);

	n.C		= eos.getC (n.ro, n.ti, n.te);
	n.Alphaei     = eos.getAlpha(n.ro, n.ti, n.te);
	n.kappa = eos.getkappa(n.ro, n.ti, n.te);
}


void CMethodEuler::fillLambda(Vector4 &Fm1, Vector4 &F0, Vector4 &Fp1, Vector4 &Fp2,
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


void CMethodEuler::fillLambdaComponent(int i, double lambda, double criteria, double curant)
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

Vector4 CMethodEuler::calcFlux(Matrix4 Omega, Matrix4 OmegaInv, Vector4 W_j_m, Vector4 W_j, Vector4 W_j_p, Vector4 W_j_p_p)
{
	Vector4 Fp = Vector4::ZERO;
	Fp += Omega * (La * (OmegaInv * W_j_m));
	Fp += Omega * (Lb * (OmegaInv * W_j));
	Fp += Omega * (Lg * (OmegaInv * W_j_p));		
	Fp += Omega * (Ld * (OmegaInv * W_j_p_p));
	return Fp;
}

Vector4 CMethodEuler::calcApprRPFlux(Vector4 W_j, Vector4 W_j_p, Node &n) {
	double pro  = getdpdro_rov_roE(n.ro,n.ti,n.te,n.v);
	double prov = getdpdrov_ro_roE(n.ro,n.ti,n.te,n.v);
	double proE = getdpdroE_ro_rov(n.ro,n.ti,n.te,n.v);
	double h	= n.e + 0.5*n.v*n.v + n.p/n.ro;
	Matrix4 A = Matrix4(0.0,		  1.,			0.0,			0.0,
						pro-n.v*n.v,  2*n.v+prov,	proE,		    0.0,
						n.v*(-h+pro), h+n.v*prov,   n.v*(1.+proE),  0.0,
						0.0, 0.0, 0.0, 1.0);
	Vector4 F_j   = A*W_j;
	Vector4 F_j_p = A*W_j_p;
	Matrix4  O = getOmega(n);
	Matrix4 OI = getOmegaInv(n);
	Vector4  L = getLambda(n);
	int j=0;
	L[0] -= vGrid[j+1];
	L[1] -= vGrid[j+1];
	L[2] -= vGrid[j+1];
	L[3] -= vGrid[j+1];
/*
	method.La[0] = fabs (nav.v - nav.C - method.vGrid[i+1]); 	
			method.La[1] = fabs (nav.v - method.vGrid[i+1]); 
			method.La[2] = fabs (nav.v + nav.C - method.vGrid[i+1]); 
			method.La[3] = 0.; 
			Vector4 F0 = Vector4(ms[0].ro*ms[0].v, 
								 ms[0].p + ms[0].ro*ms[0].v*ms[0].v, 
								 (ms[0].ro*(ms[0].e + 0.5*ms[0].v*ms[0].v) + ms[0].p) * ms[0].v,
								 0.);
			Vector4 F1 = Vector4(ms[1].ro*ms[1].v, 
								 ms[1].p + ms[1].ro*ms[1].v*ms[1].v, 
								 (ms[1].ro*(ms[1].e + 0.5*ms[1].v*ms[1].v) + ms[1].p) * ms[1].v,
								 0.);

*/
	Vector4 F1 = 0.5 * (F_j + F_j_p);
	Vector4 F2 = OI *(W_j - W_j_p);
	Vector4 F3 = L* (OI * (W_j - W_j_p));
	Vector4 F4 = 0.5 * ( O*(L* (OI * (W_j - W_j_p))));

	Vector4 F = 0.5 * (F_j + F_j_p) + 0.5 * ( O*(L* (OI * (W_j - W_j_p))));


	return F;
}

double CMethodEuler::getdpdro_rov_roE(double ro, double ti, double te, double v) 
{
	return 0.5*(gamma-1.0)*v*v;
}


double CMethodEuler::getdpdrov_ro_roE(double ro, double ti, double te, double v) 
{
	return -(gamma-1.0)*v;
}


double CMethodEuler::getdpdroE_ro_rov(double ro, double ti, double te, double v) 
{
	return (gamma-1.0);
}

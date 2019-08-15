#include "solver.h"


// Схема Годунова-Прокопова-Сафронова
void CSolver::calcHydroStageGPS(double t, double tau) {
	EOSOld &eos = task.getEOS(); 
	double gamma = eos.getGamma();
	double E=0.; 
	int i=0, nSize = ms.getSize();
	double h=ms[1].x-ms[0].x;
	Vector4 L = Vector4::ZERO, R = Vector4::ZERO, D = Vector4::ZERO, V = Vector4::ZERO;	
	// Значения на границе (в фиктивных ячейках)
	double roLB=0., vLB = 0., ELB=0., roRB=0., vRB=0., ERB=0.;
	// Задаем граничные условия
	// Transmissive left boundary
	roLB = ms[0].ro; vLB = ms[0].v;	ELB = ms[0].e + .5*ms[0].v*ms[0].v;
	// Transmissive right boundary
	roRB = ms[nSize-1].ro; vRB = ms[nSize-1].v; ERB = ms[nSize-1].e + .5*ms[nSize-1].v*ms[nSize-1].v;
	for(i=0; i<ms.getSize(); i++) {
		Node &n = ms[i];	
		if(i==0)
			n.F = calcGPSFlux(roLB, roLB*vLB, roLB*ELB, n.W[0], n.W[1], n.W[2]);
		else
			n.F = calcGPSFlux(ms[i-1].W[0], ms[i-1].W[1], ms[i-1].W[2], n.W[0], n.W[1], n.W[2]);		
	}
	ms[nSize].F = calcGPSFlux(ms[nSize-1].W[0], ms[nSize-1].W[1], ms[nSize-1].W[2], roRB, roRB*vRB, roRB*ERB);
	// Main cycle
	for(i=0; i<nSize; i++) {				
		ms[i].W_temp = ms[i].W - tau/h*(ms[i+1].F-ms[i].F);
	}
	for(i=0; i<nSize; i++) {
		// Обновляем консервативные переменные
		Node& n=ms[i]; 
		n.W[0] = n.W_temp[0];
		n.W[1] = n.W_temp[1];
		n.W[2] = n.W_temp[2];
		// Обновляем все переменные
		n.ro = n.W[0];
		n.v  = n.W[1]/n.ro;
		n.e  = n.W[2]/n.ro - 0.5*n.v*n.v;
		n.ti = eos.getti(n.ro, n.e);
		n.p  = eos.getpi(n.ro, n.ti);
		n.C  = eos.getC(n.ro, n.ti, n.te);
	}
}


// Первый вычислительный этап: гидродинамика
// Без вакуума и, желательно, с идеальным газом. То, что работет еще со времен Хаафа.
void CSolver::calcHydroStageGushchin(double t, double tau) {
	CMethodOld &method = task.getMethod();
	method.matter2Flow(ms);
	method.advanceFlow(ms, tau);
	method.flow2Matter(ms, tau);
}

void CSolver::calcHydroStageSimpleCIR(double t, double tau){
	int i=0;
	double p_plus = 0., p_minus = 0.;
	CFieldOld ms_temp;
    ms_temp.initData(&task);
	EOSOld &eos = task.getEOS();
	int nSize = ms.getSize();
	double h = ms[1].x - ms[0].x;
	for(i=0; i<nSize; i++) 	{
		ms_temp[i].x  = ms[i].x;
		ms_temp[i].v  = ms[i].v;
		ms_temp[i].ro = ms[i].ro;
		ms_temp[i].e  = ms[i].e;
		ms_temp[i].ti = ms[i].ti;
		ms_temp[i].te = ms[i].te;
		ms_temp[i].pi = ms[i].pi;
		ms_temp[i].pe = ms[i].pe;
		ms_temp[i].p  = ms[i].p;
		ms_temp[i].Z  = ms[i].Z;
	}
	i=nSize;
	ms_temp[i].x  = ms[i].x;
	ms_temp[i].v  = ms[i].v;
	// Поскольку индекс может выходить за границы, можно воспользоваться этим, чтобы 
	// делать цикл от 0 до nSize -- индексы 0 и -1 автоматически обработаются.
	// Задаем значения на границе расчетной области.
	ms.setEdgeTransparent();
	double E_i_plus = 0., E_i_minus = 0., E_i = 0., ro_i_E_i_temp=0.;
	for(int i=0; i<nSize+1; i++) {
		if(ms[i].v >= 0) {
			//Плотность
			ms_temp[i].ro = ms[i].ro - tau/h* (ms[i].ro*ms[i].v - ms[i-1].ro*ms[i-1].v); 
			//Скорость
			ms_temp[i].v  = 1./ms_temp[i].ro* (ms[i].ro*ms[i].v - 
				            tau/h *( 0.5*(ms[i+1].pi - ms[i-1].pi) + ms[i].ro*ms[i].v*ms[i].v - ms[i-1].ro*ms[i-1].v*ms[i-1].v ));
			//Энергия
			E_i           = ms[i].ei + 0.5*ms[i].v*ms[i].v;
			E_i_minus     = ms[i-1].ei + 0.5*ms[i-1].v*ms[i-1].v;
			ro_i_E_i_temp = ms[i].ro*E_i - 
								   tau/h*(ms[i].v*( ms[i].ro*E_i + ms[i].pi) - ms[i-1].v*(ms[i-1].ro*E_i_minus + ms[i-1].pi ));
			ms_temp[i].ei = ro_i_E_i_temp/ms_temp[i].ro - 0.5*ms_temp[i].v*ms_temp[i].v;
			////////////DEBUG//////
			double ro_i = ms[i].ro;
			double v_i  = ms[i].v;
			double ei_i = ms[i].ei;
			double pi_i = ms[i].pi;

			double ro_minus = ms[i-1].ro;
			double v_minus  = ms[i-1].v;
			double ei_minus = ms[i-1].ei;
			double pi_minus = ms[i-1].pi;

			double ro_plus = ms[i+1].ro;
			double v_plus  = ms[i+1].v;
			double ei_plus = ms[i+1].ei;
			double pi_plus = ms[i+1].pi;

			double ro_temp_i = ms_temp[i].ro;

			double xxx = 0.;
			///////////////////////
		} else {
			//Плотность
			ms_temp[i].ro = ms[i].ro - tau/h* (ms[i+1].ro*ms[i+1].v - ms[i].ro*ms[i].v);
			//Скорость
			ms_temp[i].v  = 1./ms_temp[i].ro* (ms[i].ro*ms[i].v - 
				            tau/h *( 0.5*(ms[i+1].pi - ms[i-1].pi) + ms[i+1].ro*ms[i+1].v*ms[i+1].v - ms[i].ro*ms[i].v*ms[i].v ));
			//Энергия
			E_i_plus  = ms[i+1].ei + 0.5*ms[i+1].v*ms[i+1].v;		
			E_i       = ms[i].ei + 0.5*ms[i].v*ms[i].v;			
			ro_i_E_i_temp = ms[i].ro*E_i - 
								   tau/h*(ms[i+1].v*( ms[i+1].ro*E_i_plus + ms[i+1].pi) - ms[i].v*(ms[i].ro*E_i + ms[i].pi ));
			ms_temp[i].ei = ro_i_E_i_temp/ms_temp[i].ro - 0.5*ms_temp[i].v*ms_temp[i].v;
			////////////DEBUG//////
			double ro_i = ms[i].ro;
			double v_i  = ms[i].v;
			double ei_i = ms[i].ei;
			double pi_i = ms[i].pi;

			double ro_minus = ms[i-1].ro;
			double v_minus  = ms[i-1].v;
			double ei_minus = ms[i-1].ei;
			double pi_minus = ms[i-1].pi;

			double ro_plus = ms[i+1].ro;
			double v_plus  = ms[i+1].v;
			double ei_plus = ms[i+1].ei;
			double pi_plus = ms[i+1].pi;

			double ro_temp_i = ms_temp[i].ro;

			double xxx = 0.;
			///////////////////////
		}
			////////////DEBUG//////
			double ro_i = ms[i].ro;
			double v_i  = ms[i].v;
			double ei_i = ms[i].ei;
			double pi_i = ms[i].pi;

			double ro_minus = ms[i-1].ro;
			double v_minus  = ms[i-1].v;
			double ei_minus = ms[i-1].ei;
			double pi_minus = ms[i-1].pi;

			double ro_plus = ms[i+1].ro;
			double v_plus  = ms[i+1].v;
			double ei_plus = ms[i+1].ei;
			double pi_plus = ms[i+1].pi;

			double ro_temp_i = ms_temp[i].ro;

			double xxx = 0.;
			///////////////////////
			//Температура
			ms_temp[i].ti = eos.getti(ms_temp[i].ro, ms_temp[i].ei);
			ms_temp[i].te = eos.gette(ms_temp[i].ro, ms_temp[i].ti, ms_temp[i].ee);
			//Давление 
			ms_temp[i].pi = eos.getpi(ms_temp[i].ro, ms_temp[i].ti);
			//Зануляем электронное давление и энергию
			ms_temp[i].ee = 0.; ms_temp[i].e = ms_temp[i].ei + ms_temp[i].ee;
			ms_temp[i].pe = 0.; ms_temp[i].p = ms_temp[i].pi + ms_temp[i].pe;
	}
	// Тест. Пока не приравняли значения на текущем слое значениям на новом слое, удобно всё проверить.
	double eps = .01; double epsAbsMax = 1.e-5;
	double l1=0., l2=0., l3=0.;
	double r1=0., r2=0., r3=0.;
	for(int i=0; i<nSize+1; i++) {
			//Левые части законов сохранения
			l1 = (ms_temp[i].ro - ms[i].ro)/tau; 
			if(fabs(ms_temp[i].ro-ms[i].ro)<epsAbsMax)
				l1 = 0.;
			l2 = (ms_temp[i].ro * ms_temp[i].v - ms[i].ro*ms[i].v)/tau;
			l3 = (ms_temp[i].ro*(ms_temp[i].ei + 0.5*ms_temp[i].v*ms_temp[i].v) - 
				  ms[i].ro*(ms[i].ei + 0.5*ms[i].v*ms[i].v) )/tau;
			if(fabs(ms_temp[i].ei-ms[i].ei)<epsAbsMax)
				l3 = 0.;

			//Правые части
			if(ms[i].v>=0) {
				r1 = - (ms[i].ro*ms[i].v - ms[i-1].ro*ms[i-1].v)/h;
				if(fabs(ms_temp[i].ro-ms[i].ro)<epsAbsMax)
					r1 = 0.;
				r2 = - ( 0.5*(ms[i+1].pi - ms[i-1].pi) + ms[i].ro*ms[i].v*ms[i].v - ms[i-1].ro*ms[i-1].v*ms[i-1].v )/h;
				r3 = - (ms[i].v*( ms[i].ro*(ms[i].ei + 0.5*ms[i].v*ms[i].v) + ms[i].pi) - 
				    ms[i-1].v*( ms[i-1].ro*(ms[i-1].ei + 0.5*ms[i-1].v*ms[i-1].v) + ms[i-1].pi))/h; 
				if(fabs(ms_temp[i].ei-ms[i].ei)<epsAbsMax)
					r3 = 0.;
			} else {
				r1 = - (ms[i+1].ro*ms[i+1].v - ms[i].ro*ms[i].v)/h;
				if(fabs(ms_temp[i].ro-ms[i].ro)<epsAbsMax)
					r1 = 0.;
				r2 = - ( 0.5*(ms[i+1].pi - ms[i-1].pi) + ms[i+1].ro*ms[i+1].v*ms[i+1].v - ms[i].ro*ms[i].v*ms[i].v )/h;
				r3 = - (ms[i+1].v*( ms[i+1].ro*(ms[i+1].ei + 0.5*ms[i+1].v*ms[i+1].v) + ms[i+1].pi) - 
				    ms[i].v*( ms[i].ro*(ms[i].ei + 0.5*ms[i].v*ms[i].v) + ms[i].pi))/h; 
				if(fabs(ms_temp[i].ei-ms[i].ei)<epsAbsMax)
					r3 = 0.;
			}
		if(fabs( (l1-r1)/l1 ) > eps)
			cout << "Cell number " << i << ": Conservation law 1 failed!" << endl;
		if(fabs( (l2-r2)/l2 ) > eps)
			cout << "Cell number " << i << ": Conservation law 2 failed!" << endl;
		if( (l3!=0) && fabs( (l3-r3)/l3 ) > eps)
			cout << "Cell number " << i << ": Conservation law 3 failed!" << endl;
	}

	for(int i=0; i<nSize; i++) {
		 ms[i].v = ms_temp[i].v;
		ms[i].ro = ms_temp[i].ro;
		 ms[i].e = ms_temp[i].e;
		ms[i].ee = ms_temp[i].ee;
		ms[i].ei = ms_temp[i].ei;
		ms[i].p  = ms_temp[i].p;
		ms[i].pi = ms_temp[i].pi;
		ms[i].pe = ms_temp[i].pe;
		ms[i].te = ms_temp[i].te;
		ms[i].ti = ms_temp[i].ti;
		ms[i].C     = eos.getC(ms[i].ro, ms[i].ti, ms[i].te);
		ms[i].Alphaei	    = eos.getAlpha(ms[i].ro, ms[i].ti, ms[i].te);
		ms[i].kappa = eos.getkappa(ms[i].ro, ms[i].ti, ms[i].te);
		ms[i].ci    = eos.getci(ms[i].ro, ms[i].ti);
		ms[i].ce    = eos.getce(ms[i].ro, ms[i].te);
	}
	

	// Тест для схемы. Хорошо бы написать.
	////////////////////// Проверка схемы (закона сохранения энергии) в одном узле.						///////
/*	for(i = 0; i<nSize; i++)
	double  x_old = ms[i].x,   x_new = ms_temp[i].x,   v_old = ms[i].v,   v_new = ms_temp[i].v;
	double ro_old = ms[i].ro, ro_new = ms_temp[i].ro, ei_old = ms[i].ei, ei_new = ms_temp[i].ei;
	double v_old_p = ms[i+1].v, v_new_p = ms_temp[i+1].v;
	double ee_old = ms_temp[i].ee, ee_new = ms_temp[i].ee, e_new = ms_temp[i].e, e_old = ms[i].e;
	double pi_new = ms_temp[i].pi, pi_old = ms_temp[i].pi; 
	double pe_new = ms_temp[i].pe, pe_old = ms_temp[i].pe;
	///////////Оххххххх! Граничные-то условия не очень понимаю, как выглядят, надо лопатить и разбираться./////
	if(i==0) {
		if(i==0) {
				// Эта аппроксимация граничных условий на давление, быть может, 
				// чуть менее точна на автомодельной ВР в идеальном газе, но гораздо 
				// монотоннее на алюминии (не порождает жутких осцилляций на границе).
				p_next_plus  = ms_prev[i].p;
			    p_next_minus = 0.0; //-ms_prev[i].p; 
			    p_plus       = ms[i].p;
			    p_minus      = 0.0; //-ms[i].p;  
			} else if(i==nSize) {
				p_next_plus  = 0.0; //-ms_prev[i-1].p;
			    p_next_minus = ms_prev[i-1].p; 
			    p_plus       = 0.0; //-ms[i-1].p;
			    p_minus      = ms[i-1].p;  
			} else {
				p_next_plus  = ms_prev[i].p + g[i];
				p_next_minus = ms_prev[i-1].p + g[i-1];
				p_plus       = ms[i].p + g[i];
			    p_minus      = ms[i-1].p + g[i-1];
			}
	} else if (i == nSize-1) {
	} else {
		p_plus = (ms[i].p + ms[i+1].p)/2.; p_minus = (ms[i].p + ms[i-1].p)/2.;
		p_next_plus  = (ms_temp[i].p + ms_temp[i+1].p)/2.; p_next_minus = (ms_temp[i].p + ms_temp[i-1].p)/2.;
	}
	double E1 = (v_new - v_old)/tau, E2 = - (p_next_plus - p_next_minus + p_plus - p_minus)/4.0/h;
	if(fabs(E1-E2/E2) > .01)
	{
		cout << "CSolver::calcHydroStage() test error!" << endl;
	}
	double E3 = (x_new - x_old)/tau, E4 = 0.5*(v_new + v_old);
	if(fabs(E3-E4/E4) > .01)
	{
		cout << "CSolver::calcHydroStage() test error!" << endl;
	}
	double E5 = (v_new_p + v_old_p - v_new - v_old)/2.0/h, E6 = (1/ro_new - 1.0/ro_old)/tau;
	if(fabs(E6-E5/E5) > .01)
	{
		cout << "CSolver::calcHydroStage() test error!" << endl;
	}
	double E7 = (ei_new - ei_old)/tau, E8  = - (pi_new + pi_old + pVisc)*(v_new_p + v_old_p - v_new - v_old)/4./h;
	if(fabs(E7-E8/E8) > .01)
	{
		cout << "CSolver::calcHydroStage() test error!" << endl;
	}
	double E9 = (ee_new - ee_old)/tau, E10 = - (pe_new + pe_old + pVisc)*(v_new_p + v_old_p - v_new - v_old)/4./h;
	if(fabs(E9-E10/E10) > .01)
	{
		cout << "CSolver::calcHydroStage() test error!" << endl;
	}
	double te_old		= ms[i].te;
	double te_minus_new = ms_temp[i-1].te;
	double te_new		= ms_temp[i].te;
	double te_plus_new  = ms_temp[i+1].te;
	kappa_minus  = (ms[i-1].kappa + ms[i].kappa)/2.;
	kappa_plus   = 0.;//(ms[i].kappa   + ms[i+1].kappa)/2.;	
	ro_minus	 = (ms[i-1].ro    + ms[i].ro)/2.;
	ro_plus		 = (ms[i].ro      + ms[i+1].ro)/2.;
	double ro_cell			= ms[i].ro;
	double dm_cell			= ms[i].dm;
	double ce_new			= ms[i].ce;
	double expX = exp(-(ms[i].x - ms[i_pulse_min].x)/deltaSkin); double expT = exp(-(t)*(t)/tauPulse/tauPulse);
	double src			= fluence/sqrt(3.14159)/deltaSkin/tauPulse*expX*expT;
	double E1			= ce_new*te_new;
	double E2			= ce_new*te_old;
	double E3			= tau * kappa_plus  / (dm_cell*dm_cell) * ro_cell * ro_plus * (te_plus_new - te_new);
	double E4			= tau * kappa_minus / (dm_cell*dm_cell) * ro_cell * ro_minus * (te_new - te_minus_new);	
	double E5			= src * tau;
	double left			= E1 - E2;
	double right		= E3 - E4 + E5;
	double delta		= left - right;
	//////////////////////////////////////////////////////////////////
*/
	cout << "calcHydroStageSimpleCIR() complete!" << endl;
}


// Упрощаем функцию по расчету одного шага схемы -- выкидываем внутренние зависимости и 
// обновление температур через УРС на каждом этапе. Обновляем только энергии.
// Нужно ли обновлять давления? -- это вопрос, надо разобраться.
void CSolver::calcHydroStageGushchinIdealSimple(double t, double tau) {
	EOSOld &eos = task.getEOS();
	CMethodOld &method = task.getMethod();
	ms_temp.initData(&task);
	double E=0.; int i=0;
	double test=0., h=0;
	double L0=0., L1=0., L2=0., L3=0., R0=0., R1=0., R2=0., R3=0., epsilon=0.; 
	double dx = ms[1].x-ms[0].x;
	for(i=0; i<ms.getSize(); i++) {
		Node &n = ms[i];
		E = n.e + 0.5*n.v*n.v;
		n.W = Vector4(n.ro, n.ro*n.v, n.ro*E, 0);
		n.F = method.getOmegaInv(n) * n.W;
	}
	// Внимание, граничные условия!!!
	// Сейчас просто сношу в фиктивную ячейку все значения в нулевой точке. Делаю, чтобы сетка не двигалась.
	
	/// Здесь нужно дать возможность для разных типов граничных условий -- стенка, прозрачное, периодическое
	// возможно -- прибавлять просто по две фиктивных ячейки без усложнений с индексированием.

	ms.setEdgeTransparent();

	///

	Matrix4 O, OI;
	Vector4 Fp, Fm;
	Vector4 L;
	Node	nav;
	int j=0, nSize=ms.getSize(); //nSize это количество узлов сетки (и количество промежутков плюс один)
	Vector4 V1 = Vector4::ZERO, V2 = Vector4::ZERO, V3 = Vector4::ZERO, V4 = Vector4::ZERO;
	// А вот это хорошо бы брать из УРС
	double gam = eos.getGamma();
	// Поехали
	for(i=0; i<ms.getSize(); i++)
	{
		// F(i+1/2) flow
		j = i;
		if(i!=nSize-1) 
			method.averageNode(ms[j], ms[j+1], nav);
		else
			method.averageNode(ms[j], ms[j], nav);
		O  = method.getOmega(nav);
		OI = method.getOmegaInv(nav);
		L  = method.getLambda(nav);
		if(j<nSize-2) {
			method.fillLambda(ms[j-1].F, ms[j].F, ms[j+1].F, ms[j+2].F,
					   L, tau/dx);
			Fp = method.calcFlux(O, OI, ms[j-1].W, ms[j].W, ms[j+1].W, ms[j+2].W);
		} else if (j==nSize-2) {
			method.fillLambda(ms[j-1].F, ms[j].F, ms[j+1].F, ms[j+1].F,
					   L, tau/dx);
			Fp = method.calcFlux(O, OI, ms[j-1].W, ms[j].W, ms[j+1].W, ms[j+1].W);
		} else if (j==nSize-1) {
			method.fillLambda(ms[j-1].F, ms[j].F, ms[j].F, ms[j].F,
					   L, tau/dx);
			Fp = method.calcFlux(O, OI, ms[j-1].W, ms[j].W, ms[j].W, ms[j].W);
		}
		// F(i-1/2) flow
		j = i-1;
		method.averageNode(ms[j], ms[j+1], nav);
		O  = method.getOmega(nav);
		OI = method.getOmegaInv(nav);
		L  = method.getLambda(nav);
		if(j>0) {
			method.fillLambda(ms[j-1].F, ms[j].F, ms[j+1].F, ms[j+2].F,
					   L, tau/dx);
			Fm = method.calcFlux(O, OI, ms[j-1].W, ms[j].W, ms[j+1].W, ms[j+2].W);
		} else {
			method.fillLambda(ms[j].F, ms[j].F, ms[j+1].F, ms[j+2].F,
					   L, tau/dx);
			Fm = method.calcFlux(O, OI, ms[j].W, ms[j].W, ms[j+1].W, ms[j+2].W);
		}
		//////////////////////////////
		// Вычисляем W_temp 
		ms[i].W_temp = ms[i].W - (Fp-Fm) * tau/dx;
	}
	/////////////////Здесь тестирование этапа////////////////////////////////////////////////
	/////////////////Начало тестирования/////////////////////////////////////////////////////
/*	Vector4 g = Vector4::ZERO;
	Vector4 F_p = Vector4::ZERO, F_m = Vector4::ZERO; 
	Vector4 alpha = Vector4::ZERO, beta = Vector4::ZERO, gamma = Vector4::ZERO, delta = Vector4::ZERO;
	Vector4 lambda = Vector4::ZERO;
	Matrix4 Omega = Matrix4::ZERO, OmegaInv = Matrix4::ZERO;
	Vector4 W_im2 = Vector4::ZERO, W_im1 = Vector4::ZERO, W_i = Vector4::ZERO,
		    W_ip1 = Vector4::ZERO, W_ip2 = Vector4::ZERO;
	Vector4 g_im2 = Vector4::ZERO, g_im1 = Vector4::ZERO, g_i = Vector4::ZERO,
		    g_ip1 = Vector4::ZERO, g_ip2 = Vector4::ZERO;
	int testFailed = 0;
	for(i=0; i<nSize; i++) {
		h = ms[i+1].x-ms[i].x;
	// Заполняем все W и g
    // g_im2
 		E = 0.;		
		W_im2 = Vector4::ZERO, W_im1 = Vector4::ZERO, W_i = Vector4::ZERO, W_ip1 = Vector4::ZERO, W_ip2 = Vector4::ZERO;
		g_im2 = Vector4::ZERO, g_im1 = Vector4::ZERO, g_i = Vector4::ZERO, g_ip1 = Vector4::ZERO, g_ip2 = Vector4::ZERO;
		if(i>=2) {
			E = ms[i-2].ei + 0.5*ms[i-2].v*ms[i-2].v;
			OmegaInv = method.getOmegaInv(ms[i-2]);
			W_im2 = Vector4(ms[i-2].ro, ms[i-2].ro*ms[i-2].v, ms[i-2].ro*E, 0.);
			g_im2 = OmegaInv*W_im2;			
		}
		if(i>=1) {
			E = ms[i-1].ei + 0.5*ms[i-1].v*ms[i-1].v;
			OmegaInv = method.getOmegaInv(ms[i-1]);
			W_im1 = Vector4(ms[i-1].ro, ms[i-1].ro*ms[i-1].v, ms[i-1].ro*E, 0.);
			g_im1 = OmegaInv*W_im1;			
		}
		E = ms[i].ei + 0.5*ms[i].v*ms[i].v;
		OmegaInv = method.getOmegaInv(ms[i]);
		W_i = Vector4(ms[i].ro, ms[i].ro*ms[i].v, ms[i].ro*E, 0.);
		g_i = OmegaInv*W_i;			
		if(i<=nSize-2) {
			E = ms[i+1].ei + 0.5*ms[i+1].v*ms[i+1].v;
			OmegaInv = method.getOmegaInv(ms[i+1]);
			W_ip1 = Vector4(ms[i+1].ro, ms[i+1].ro*ms[i+1].v, ms[i+1].ro*E, 0.);
			g_ip1 = OmegaInv*W_ip1;			
		} else {
			W_ip1 = W_i;
		}
		if(i<=nSize-3) {
			E = ms[i+2].ei + 0.5*ms[i+2].v*ms[i+2].v;
			OmegaInv = method.getOmegaInv(ms[i+2]);
			W_ip2 = Vector4(ms[i+2].ro, ms[i+2].ro*ms[i+2].v, ms[i+2].ro*E, 0.);
			g_ip2 = OmegaInv*W_ip2;			
		} else {
			W_ip2 = W_i;
		}
		// Теперь потоки
		// F_p
		Node nav;
		F_p = Vector4::ZERO;
		method.averageNode(ms[i], ms[i+1], nav);
		Omega    = method.getOmega(nav);
		OmegaInv = method.getOmegaInv(nav);
		lambda   = method.getLambda(nav);
		method.fillLambda(g_im1, g_i, g_ip1, g_ip2,
							 lambda, tau/(ms[i+1].x - ms[i].x));
		F_p = Vector4::ZERO;
		F_p += Omega * (method.La * (OmegaInv * W_im1));
		F_p += Omega * (method.Lb * (OmegaInv * W_i));
		F_p += Omega * (method.Lg * (OmegaInv * W_ip1));
		F_p += Omega * (method.Ld * (OmegaInv * W_ip2));
		//F_m
		// В нулевой ячейке ставим его равным 0: граница с вакуумом
		if(i!=0) {
			method.averageNode(ms[i-1], ms[i], nav);
			Omega    = method.getOmega(nav);
			OmegaInv = method.getOmegaInv(nav);
			lambda   = method.getLambda(nav);
			method.fillLambda(g_im2, g_im1, g_i, g_ip1,
								 lambda, tau/(ms[i+1].x - ms[i].x));
			F_m = Vector4::ZERO;
			F_m += Omega * (method.La * (OmegaInv * W_im2));
			F_m += Omega * (method.Lb * (OmegaInv * W_im1));
			F_m += Omega * (method.Lg * (OmegaInv * W_i));
			F_m += Omega * (method.Ld * (OmegaInv * W_ip1));
		} else {
			F_m = Vector4::ZERO;
		}
		// W_temp
		//ms[i].W_temp = ms[i].W - (Fp-Fm) * tau/(ms[i+1].x - ms[i].x);
	}	*/
	/////////////////Конец тестирования///////////////////////////////////////////////////////
	for(i=0; i<ms.getSize(); i++) {
		ms[i].W = ms[i].W_temp;
	}
	E=0.;
	for(int i=0; i<ms.getSize(); i++) {
		Node &n = ms[i];
		n.ro = n.W[0];
		n.v  = (n.ro!=0.0) ? n.W[1] / n.ro : 0.0;
		  E  = (n.ro!=0.0) ? n.W[2] / n.ro : 0.0;
		n.e  = E - 0.5*n.v*n.v;
		n.p  = (gam-1.) * n.ro * n.e;
		n.C  = sqrt(gam*n.p/n.ro);
	}
}




// Первый вычислительный этап: гидродинамика
// С вакуумом и подвижной сеткой. То, на основе чего мы собираемся считать реальные задачи с металлом.
void CSolver::calcHydroStageGushchinMovingGrid(double t, double tau) {
	EOSOld &eos = task.getEOS();
	CMethodOld &method = task.getMethod();
	ms_temp.initData(&task);
	//Сетка уже сгенерирована в функции go()
	double E=0.; int i=0;
	double test=0., h=0;
	double L0=0., L1=0., L2=0., L3=0., R0=0., R1=0., R2=0., R3=0., epsilon=0.; 
	for(i=0; i<ms.getSize(); i++) {
		Node &n = ms[i];
		E = n.e + 0.5*n.v*n.v;
		n.W = Vector4(n.ro, n.ro*n.v, n.ro*E, 0);
		n.F = method.getOmegaInv(n) * n.W;
	}
	// Сейчас просто сношу в фиктивную ячейку все значения в нулевой точке. Делаю, чтобы сетка не двигалась.
	ms.setEdgeTransparent();
	Matrix4 O, OI;
	Vector4 Fp, Fm;
	Vector4 L;
	Node	nav;
	int j=0, nSize=ms.getSize(); //nSize это количество узлов сетки (и количество промежутков плюс один)
	Vector4 V1 = Vector4::ZERO, V2 = Vector4::ZERO, V3 = Vector4::ZERO, V4 = Vector4::ZERO;
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
	// Сохраняем скорость левого конца для корректного расчета движения сетки
	// DEBUG 18.09.2014 -- делаю скорость сетки пока нулевой, тестирую сам солвер ////
	double gam = 1.4;
	Zone z = task.getZone(0);
	double vLeftEdge = z.v - 2.*sqrt(gam*eos.getpi(z.ro, z.ti)/z.ro)/(gam-1.);
	//////////////////////////////////////////////////////////////////////////////////
	// Теперь обновляем сеточные координаты и скорости для следующего этапа.
	// Задаем скорости сетки на этапе. Условия: 1. Сетка всегда равномерная. 2. Скорость левого конца равна скорости v[0].
	double _vg = 0.;
	for(i=0; i<nSize+1; i++)
	{
		j = nSize - i;
		//DEBUG//
		//Исправление номер один: скорость вычисляем из волны разрежения (вопрос -- всей волны или между последней ячейкой и вакуумом?) 
		//double _v0 = /*-2.*(ms[0].C)/(5./3.-1.) - ms[0].v;*/ms[0].v;
		double _v0 = vLeftEdge;
			////////////////
			// Здесь казус №1: расчетная скорость получается по модулю больше, чем аналитическая, а они должны быть одинаковы в силу того
			// способа, которым мы считаем скорость (скорость левого узла сетки равна аналитической) 
			// Казус №2 -- надо доразобраться с сеткой, что в целых точках, что в полуцелых. И будет, я думаю, всё ОК.
			// Номер 3 -- обращаю внимание на то, что со сменой парадигмы nSize (теперь в коде и в конфиге снаружи он означает
			// количество ячеек, а не количество узлов) я заменил расчеты по сетке (везде, где были разности nSize-1 я поставил nSize)
			// Номер 4 -- проблема может быть еще и в том, что вне вещества все W строго 0. А расчет потока опирается у нас на
			// W, выходящие за границы вещества, например, поток F_0-1/2 зависит от W_0-1. Разобраться, при необходимости спросить.
			////////////////
		////////////////////////
		//_vg = ms[0].v*j / (nSize);
		_vg = _v0*j / (nSize);
		method.vGrid[i]= _vg;
		method.X[i] += method.vGrid[i]*tau;
	}
	
	//double v0 = method.vGrid[0], v1=method.vGrid[1], v2=method.vGrid[2];
	// Сдвигаем сетку
	//method.X[0] = ms[0].x + v0*tau;
	//for(i=0; i<nSize+1; i++){
		//X[i] = ms[i].x + vGrid[i]*tau;
		//method.X[i] = ms[0].x + (ms[nSize].x - ms[0].x)/(nSize) * i;
		// Понятно, ms[nSize].x где-то недосчитывается, проверить, чтобы везде досчитывалось. Из-за этого и сдвиг. Проверить и пофиксить.
	//	method.X[i] = method.X[0] + (ms[nSize].x - method.X[0])/(nSize) * i;
		
	//}
	// Поехали
	for(i=0; i<ms.getSize(); i++)
	{
		// F(i+1/2) flow
		j = i;
		method.averageNode(ms[j], ms[j+1], nav);
		O  = method.getOmega(nav);
		OI = method.getOmegaInv(nav);
		L  = method.getLambda(nav);
		// Для подвижной сетки -- модифицируем собственные значения для расчета модифицированного потока
		L[0]-=method.vGrid[i+1];
		L[1]-=method.vGrid[i+1];
		L[2]-=method.vGrid[i+1];
		L[3]-=method.vGrid[i+1];
		method.fillLambda(ms[j-1].F, ms[j].F, ms[j+1].F, ms[j+2].F,
				   L, tau/(ms[j+1].x - ms[j].x));
		Fp = method.calcFlux(O, OI, ms[j-1].W, ms[j].W, ms[j+1].W, ms[j+2].W);
		// Альтернативные способы подсчета потока, потом посмотрим. Наш поток -- точный римановский. :) :) :)
		//Fp = Vector4(nav.ro*nav.v, nav.ro*nav.v*nav.v+nav.p, (nav.ro*(nav.e+0.5*nav.v*nav.v) + nav.p)*nav.v, 0.);
		//method.calcFlux(O, OI, -ms[j].W, ms[j].W, ms[j+1].W, ms[j+2].W);
		//Fp = method.calcApprRPFlux(ms[0].W, ms[1].W, nav);
		////////////DEBUG//////////
	/*	if(i==0) {	
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


			Fp = 0.5 * (ms[0].F + ms[1].F) + 0.5 * ( O*(method.La * (OI * (ms[0].W - ms[1].W))));
		} else {
		*/
		//////////////////////////
		//Fp += O * (method.La * (OI * ms[j-1].W));
		//Fp += O * (method.Lb * (OI * ms[j  ].W));
		//Fp += O * (method.Lg * (OI * ms[j+1].W));		
		//Fp += O * (method.Ld * (OI * ms[j+2].W));
		//}
		/////DEBUG////////////////////
		//if(i==0) {
		//	cout << "Fp_calc[" << i << "] = " << Fp << endl;
		//}
		//////////////////////////////
		// F(i-1/2) flow
		j = i-1;
		method.averageNode(ms[j], ms[j+1], nav);
		O  = method.getOmega(nav);
		OI = method.getOmegaInv(nav);
		L  = method.getLambda(nav);
		// Модифицируем собственные значения для учета подвижности сетки
		L[0]-=method.vGrid[i];
		L[1]-=method.vGrid[i];
		L[2]-=method.vGrid[i];
		L[3]-=method.vGrid[i];
		method.fillLambda(ms[j-1].F, ms[j].F, ms[j+1].F, ms[j+2].F,
				   L, tau/(ms[j+1].x - ms[j].x));
		
		//if(i==0) {
		//	Fm = Vector4::ZERO;
		////////// DEBUG ////////////////////////////
		// Пытаемся учесть то, что расчет потока в ячейке с номером 1 все-таки захватывает вакуум,
		// поэтому говорим схеме, что в вакууме решение как будто такое же, как и в первом потоке.
		//
		//} else if(i==1) {	
		//	Fm = method.calcApprRPFlux(ms[0].W, ms[1].W, nav);
		////////////////////////////////////////////
		//} else if(i==1) {
		//	Fp = Vector4(nav.ro*nav.v, nav.ro*nav.v*nav.v+nav.p, (nav.ro*(nav.e+0.5*nav.v*nav.v) + nav.p)*nav.v, 0.);
		//	//Fm = method.calcFlux(O, OI, -ms[j].W, ms[j].W, ms[j+1].W, ms[j+2].W);
		//} else {
		//////////////////////////
		Fm = method.calcFlux(O, OI, ms[j-1].W, ms[j].W, ms[j+1].W, ms[j+2].W);
		//}
		/////DEBUG////////////////////
		if(i==0) {
			Fm = Vector4(0., 0.001, 0., 0.);
			cout << "Fm_calc[" << i << "] = " << Fm << endl;
		}
		//////////////////////////////
		// Вычисляем W_temp 
		//ms[i].W_temp = ms[i].W - (Fp-Fm) * tau/(ms[i+1].x - ms[i].x);
	///DEBUG
		/*if(i!=0 && i!=1)
		{*/
	///
		ms[i].W_temp = ( (ms[i+1].x-ms[i].x)*ms[i].W - (Fp-Fm) * tau) / (method.X[i+1] - method.X[i]);
	
	///DEBUG
			/*	}else{
			
			//Плотность
			ms_temp[i].ro = (ms[i+1].x-ms[i].x)*ms[i].ro - 
				             tau /(method.X[i+1] - method.X[i]) * 
							 (ms[i+1].ro*(ms[i+1].v-method.vGrid[i+1]) - ms[i].ro*(ms[i].v-method.vGrid[i]));

			//Скорость
			ms_temp[i].v  = 1./ms_temp[i].ro* ((ms[i+1].x-ms[i].x)*ms[i].ro*ms[i].v - 
				            tau/ (method.X[i+1] - method.X[i]) *
							( 0.5*(ms[i+1].pi - ms[i-1].pi) + ms[i+1].ro*(ms[i+1].v-method.vGrid[i+1])*(ms[i+1].v-method.vGrid[i+1]) - ms[i].ro*(ms[i].v-method.vGrid[i])*(ms[i].v-method.vGrid[i]) ));
			//Энергия
			double E_i_plus  = ms[i+1].ei + 0.5*ms[i+1].v*ms[i+1].v;		
			double E_i       = ms[i].ei + 0.5*ms[i].v*ms[i].v;			
			double ro_i_E_i_temp = (ms[i+1].x-ms[i].x)*ms[i].ro*E_i - 
								   tau/(method.X[i+1] - method.X[i])*((ms[i+1].v-method.vGrid[i+1])*( ms[i+1].ro*E_i_plus + ms[i+1].pi) - (ms[i].v-method.vGrid[i+1])*(ms[i].ro*E_i + ms[i].pi ));
			ms_temp[i].ei = ro_i_E_i_temp/ms_temp[i].ro - 0.5*ms_temp[i].v*ms_temp[i].v;
			ms[i].W_temp = Vector4(ms_temp[i].ro, ms_temp[i].ro*ms_temp[i].v, ro_i_E_i_temp, 0.);
		}*/
	///
		
/*def getRPSolution(roL, uL, pL, roR, uR, pR, x, t):
    xi = x / t
    gamma = 1.4
    _p = 0.
    _u = 0.
    _ro = 0.
    _sigma = 0.
    cL = sqrt(gamma * pL / roL)
    cR = sqrt(gamma * pR / roR)
    sigmaL = pL / roL ** gamma
    sigmaR = pR / roR ** gamma
    if xi < uL - cL:
        return np.array([roL, uL, pL])
    if xi > uR + cR:
        return np.array([roR, uR, pR])
    _p = (pL / roL / cL + pR / roR / cR + uR - uL) / (1. / roL / cL + 1 / roR / cR)
    _u = (roL * cL * uL + roR * cR * uR + pR - pL) / (roL * cL + roR * cR)
    if _u <= xi:
        _sigma = sigmaR
    else:
        _sigma = sigmaL
    _ro = (_p / _sigma) ** 1. / gamma
    #    if uL>cL:
    #        return np.array([roL, uL, pL])
    #    elif uR<-cR:
    #        return np.array([roR, uR, pR])
    #    else:
    #        if _u > 0:
    #            _ro = roR + (_p - pR)/cR*cR
    #        else:
    #           _ro = roL + (_p - pL)/cL*cL
    return np.array([_ro, _u, _p])		*/
		
		
		
		
		
		///// DEBUG ///////////////////
/*		if (i == 1) {

			test = ms[i].W[1] - (Fp[1]-Fm[1]) * tau/h;
			cout << ms[i].W_temp << endl;
			cout << ms[i].W << endl;
		}*/
		/////////////////////////
	}

	
	
	


	/////////////////Здесь тестирование этапа////////////////////////////////////////////////
	/////////////////Начало тестирования/////////////////////////////////////////////////////
	Vector4 g = Vector4::ZERO;
	Vector4 F_p = Vector4::ZERO, F_m = Vector4::ZERO; 
	Vector4 alpha = Vector4::ZERO, beta = Vector4::ZERO, gamma = Vector4::ZERO, delta = Vector4::ZERO;
	Vector4 lambda = Vector4::ZERO;
	Matrix4 Omega = Matrix4::ZERO, OmegaInv = Matrix4::ZERO;
	Vector4 W_im2 = Vector4::ZERO, W_im1 = Vector4::ZERO, W_i = Vector4::ZERO,
		    W_ip1 = Vector4::ZERO, W_ip2 = Vector4::ZERO;
	Vector4 g_im2 = Vector4::ZERO, g_im1 = Vector4::ZERO, g_i = Vector4::ZERO,
		    g_ip1 = Vector4::ZERO, g_ip2 = Vector4::ZERO;
	int testFailed = 0;
	for(i=0; i<nSize; i++) {
		h = ms[i+1].x-ms[i].x;
	// Заполняем все W и g
    // g_im2
 		E = 0.;		
		W_im2 = Vector4::ZERO, W_im1 = Vector4::ZERO, W_i = Vector4::ZERO, W_ip1 = Vector4::ZERO, W_ip2 = Vector4::ZERO;
		g_im2 = Vector4::ZERO, g_im1 = Vector4::ZERO, g_i = Vector4::ZERO, g_ip1 = Vector4::ZERO, g_ip2 = Vector4::ZERO;
		if(i>=2) {
			E = ms[i-2].ei + 0.5*ms[i-2].v*ms[i-2].v;
			OmegaInv = method.getOmegaInv(ms[i-2]);
			W_im2 = Vector4(ms[i-2].ro, ms[i-2].ro*ms[i-2].v, ms[i-2].ro*E, 0.);
			g_im2 = OmegaInv*W_im2;			
		}
		if(i>=1) {
			E = ms[i-1].ei + 0.5*ms[i-1].v*ms[i-1].v;
			OmegaInv = method.getOmegaInv(ms[i-1]);
			W_im1 = Vector4(ms[i-1].ro, ms[i-1].ro*ms[i-1].v, ms[i-1].ro*E, 0.);
			g_im1 = OmegaInv*W_im1;			
		}
		E = ms[i].ei + 0.5*ms[i].v*ms[i].v;
		OmegaInv = method.getOmegaInv(ms[i]);
		W_i = Vector4(ms[i].ro, ms[i].ro*ms[i].v, ms[i].ro*E, 0.);
		g_i = OmegaInv*W_i;			
		if(i<=nSize-2) {
			E = ms[i+1].ei + 0.5*ms[i+1].v*ms[i+1].v;
			OmegaInv = method.getOmegaInv(ms[i+1]);
			W_ip1 = Vector4(ms[i+1].ro, ms[i+1].ro*ms[i+1].v, ms[i+1].ro*E, 0.);
			g_ip1 = OmegaInv*W_ip1;			
		} else {
			W_ip1 = W_i;
		}
		if(i<=nSize-3) {
			E = ms[i+2].ei + 0.5*ms[i+2].v*ms[i+2].v;
			OmegaInv = method.getOmegaInv(ms[i+2]);
			W_ip2 = Vector4(ms[i+2].ro, ms[i+2].ro*ms[i+2].v, ms[i+2].ro*E, 0.);
			g_ip2 = OmegaInv*W_ip2;			
		} else {
			W_ip2 = W_i;
		}
		// Теперь потоки
		// F_p
		Node nav;
		F_p = Vector4::ZERO;
		method.averageNode(ms[i], ms[i+1], nav);
		/////////////DEBUG////////////////////
		if(i == 2){
			double f = 0.;
		}
		/////////////////////////////////////
		Omega    = method.getOmega(nav);
		OmegaInv = method.getOmegaInv(nav);
		lambda   = method.getLambda(nav);
		/// Модифицируем собственные значения, вычитаем скорость сетки
		lambda[0]-=method.vGrid[i+1];
		lambda[1]-=method.vGrid[i+1];
		lambda[2]-=method.vGrid[i+1];
		lambda[3]-=method.vGrid[i+1];
		method.fillLambda(g_im1, g_i, g_ip1, g_ip2,
							 lambda, tau/(ms[i+1].x - ms[i].x));
		F_p = Vector4::ZERO;
		F_p += Omega * (method.La * (OmegaInv * W_im1));
		F_p += Omega * (method.Lb * (OmegaInv * W_i));
		F_p += Omega * (method.Lg * (OmegaInv * W_ip1));
		F_p += Omega * (method.Ld * (OmegaInv * W_ip2));
		/////////////DEBUG////////////////////
		//if(i == 3){
		//	V1 = OI * W_i;
		//	V2 = method.Lb * V1;
		//	V3 = O * V3;
		//}

		/////////////////////////////////////
		/////DEBUG////////////////////
		if(i==0) {
				cout << "Fp_test[" << i << "] = " << F_p << endl;
		}
		//////////////////////////////
		//F_m
		// В нулевой ячейке ставим его равным 0: граница с вакуумом
		if(i!=0) {
			method.averageNode(ms[i-1], ms[i], nav);
			Omega    = method.getOmega(nav);
			OmegaInv = method.getOmegaInv(nav);
			lambda   = method.getLambda(nav);
			/// Модифицируем собственные значения, вычитаем скорость сетки
			lambda[0]-=method.vGrid[i];
			lambda[1]-=method.vGrid[i];
			lambda[2]-=method.vGrid[i];
			lambda[3]-=method.vGrid[i];
			method.fillLambda(g_im2, g_im1, g_i, g_ip1,
								 lambda, tau/(ms[i+1].x - ms[i].x));
			F_m = Vector4::ZERO;
			/*if(i==1) {
				V1 = OmegaInv * W_im2;
				V2 = method.La * V1;
				V3 = Omega * V2;
			}*/
			F_m += Omega * (method.La * (OmegaInv * W_im2));
			F_m += Omega * (method.Lb * (OmegaInv * W_im1));
			F_m += Omega * (method.Lg * (OmegaInv * W_i));
			F_m += Omega * (method.Ld * (OmegaInv * W_ip1));
		} else {
			F_m = Vector4::ZERO;
		}
		/////DEBUG////////////////////
		if(i==0) {
			cout << "Fm_test[" << i << "] = " << F_m << endl;
		}
		//////////////////////////////
		//ms[i].W_temp = ( (ms[i+1].x-ms[i].x)*ms[i].W - (Fp-Fm) * tau) / (method.X[i+1] - method.X[i]);

		L0  = ( (method.X[i+1]-method.X[i])*ms[i].W_temp[0] - (ms[i+1].x-ms[i].x)*ms[i].W[0] ); R1 = -(F_p[0] - F_m[0])*tau;
		L1  = ( (method.X[i+1]-method.X[i])*ms[i].W_temp[1] - (ms[i+1].x-ms[i].x)*ms[i].W[1] ); R1 = -(F_p[1] - F_m[1])*tau;
		L2  = ( (method.X[i+1]-method.X[i])*ms[i].W_temp[2] - (ms[i+1].x-ms[i].x)*ms[i].W[2] ); R2 = -(F_p[2] - F_m[2])*tau;
		L3  = ( (method.X[i+1]-method.X[i])*ms[i].W_temp[3] - (ms[i+1].x-ms[i].x)*ms[i].W[3] ); R3 = -(F_p[3] - F_m[3])*tau;
		epsilon = .01; 
		testFailed = 0;
		// First component
		if(L0 == 0. && fabs(R0)>epsilon && fabs(ms[i].W_temp[0]-ms[i].W[0])> epsilon)
			testFailed++;
		if(L0!=0 && (fabs((L0-R0)/L0) > epsilon) && (fabs(ms[i].W_temp[0]-ms[0].W[0])> epsilon) )
			testFailed++;
		if(testFailed) {
			// Так считается W-temp: ms[i].W_temp = ms[i].W - (Fp-Fm) * tau/(ms[i+1].x - ms[i].x);

			L0 = ms[i].W_temp[0];
			R0 = ms[i].W[0]-(F_p[0]-F_m[0])*tau/(ms[i+1].x-ms[i].x);
			if(fabs((L0-R0)/L0)>epsilon) {
				cout << "node " << i << ": Conservation law violated in W[0]" << endl;	
			}
		}
		testFailed = 0;
		if(L1==0. && fabs(R1)>epsilon)
			testFailed++;
		if(L1!=0. && fabs((L1-R1)/L1) > epsilon && fabs(ms[i].W_temp[2]-ms[i].W[2])> epsilon) 
			testFailed++;
		if(testFailed) {
			L1 = ms[i].W_temp[1];
			R1 = ms[i].W[1]-(F_p[1]-F_m[1])*tau/(ms[i+1].x-ms[i].x);
			if( L1!=0. && fabs((L1-R1)/L1)>epsilon && fabs(ms_temp[i].v - ms[i].v)>epsilon ) {
				cout << "node " << i << ": Conservation law violated in W[1]" << endl;	
			}
		}
		testFailed = 0;
		if( (fabs((L2-R2)/L2) > epsilon && fabs(ms[i].W_temp[2]-ms[i].W[2])> epsilon) )
			cout << "node " << i << ": Conservation law violated in W[2]" << endl;	
		if(fabs((L3-R3)/L3) > epsilon)
			cout << "node " << i << ": Conservation law violated in W[3]" << endl; 
	}	
	/////////////////Конец тестирования///////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////






	


	////////// DEBUG ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// А здесь тест на целостность сетки.
	double LG = 0., RG = 0.;
	epsilon = .01;
	for(i=0; i<nSize; i++)
	{
		LG = (method.X[i+1]-method.X[i]) - (ms[i+1].x-ms[i].x);
		RG = (method.vGrid[i+1]-method.vGrid[i])*tau;
		if( LG!=0. && fabs(LG)>1.e-18 && fabs((LG-RG)/LG)>epsilon )
			cout << "Something wrong with grid in node " << i << endl; 
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	for(i=0; i<=ms.getSize(); i++){
		ms[i].x  = method.X[i];
	}

	for(i=0; i<ms.getSize(); i++) {
		ms[i].W = ms[i].W_temp;
	}
	////////////////////////DEBUG/////////////////
	//Третий, последний. Мы свободны от старых методов:
	//method.flow2MatterVacuum(ms, ms_temp, t, tau);
	//////////////////////////////////////////////
	E=0.;
	for(int i=0; i<ms.getSize(); i++) {
		Node &n = ms[i];
		n.ro = n.W[0];
		n.v  = (n.ro!=0.0) ? n.W[1] / n.ro : 0.0;
		  E  = (n.ro!=0.0) ? n.W[2] / n.ro : 0.0;
		n.e  = E - 0.5*n.v*n.v;
		n.ei = n.e;
		n.ee = 0;		
		method.updateNode(n);
		//////////////
		double q=0.;
		/////////////
	}
}

// Function finds resulting ro, v, p of Riemann task
RPSolutionPrimitive CSolver::solveRP(double roL, double vL, double pL, double roR, double vR, double pR) {
	// Решаем нелинейное уравнение относительно давления методом касательных Ньютона
	RPSolutionPrimitive res; res.roL = 0.; res.roR=0.; res.v = 0.; res.p = 0.;
	double p = 0., pPrev = 0.;
	double TOL = 1.e-6;
	double gamma = 1.4;
	int itCounter = 0;
	double cL = 0., cR = 0.;
	if(roL!=0.) cL = sqrt(gamma*pL/roL);
	if(roR!=0.) cR = sqrt(gamma*pR/roR);
	// Пытаюсь определить возможную конфигурацию решения, чтобы вернее выставить начальное приближение
	// Похоже, итерации нужны только в случаях "УВ+УВ" и "УВ + ВР", т.к. в случае ВР+ВР и ВР+вакуум есть 
	// аналитические решения для идеального газа
	//
	// Также вызывает вопрос последний тест Торо, где полученное решение отличается от его решения 
	// во втором знаке после запятой
	if(roL==roR && vL==vR && pL==pR) {
		res.type = RWRW;
		res.roL  = roL;
		res.roR  = roL;
		res.p    = pL;
		res.v	 = vL;
		return res;
	}
	if(roL==0.) {
		res.type = VacRW;
		res.roL  = 0.;
		res.roR  = roR;
		res.p	 = 0.;
		res.v    = vR - 2.*cR/(gamma-1.);
		return res;
	}
	if(roR==0.) {
		res.type = RWVac;
		res.roL  = roL;
		res.roR  = 0.;
		res.p	 = 0.;
		res.v    = vL + 2.*cL/(gamma-1.);
		return res;
	}
	if(2.*cL/(gamma-1) + 2*cR/(gamma-1.) < fabs(vL-vR)){
		res.type  = RWVacRW;
		res.roL	  = 0.;
		res.roR   = 0.;
		res.v     = 0.;
		res.p	  = 0.;
		return res;
	}

	double fLmin = fL(pL, roL, vL, pL) + fR(pL, roR, vR, pR) + vR-vL;
	double fRMax = fL(pR, roL, vL, pL) + fR(pR, roR, vR, pR) + vR-vL;
	// Начальное приближение
	//p = 0.5*(pL+pR);
	p=pL/2.;
	do {
		pPrev = p;
		p = pPrev - (fL(pPrev, roL, vL, pL) + fR(pPrev, roR, vR, pR) + vR - vL )/
			        (dfLdp(pPrev, roL, vL, pL) + dfRdp(pPrev, roR, vR, pR)); 
		if (p<=0.)
			p = TOL;
		itCounter++;
	} while (fabs(2*(p-pPrev)/(p+pPrev))>TOL);
	res.p   = p;
	res.v   = 0.5*(vL + vR) + 0.5*(fR(p, roR, vR, pR) - fL(p, roL, vL, pL));
	if( p<pL && p>pR) {
		res.type = RWSW;
		res.roL  = roL*pow(res.p/pL, 1./gamma); 
		res.roR  = roR*(res.p/pR + (gamma-1.)/(gamma+1.))/((gamma-1.)/(gamma+1.)*res.p/pR + 1.);
	} else if(p<=pL && p<=pR) {
		res.type = RWRW;
		res.roL  = roL*pow(res.p/pL, 1./gamma); 
		res.roR  = roR*pow(res.p/pR, 1./gamma); 
	} else if(p>pL && p<pR) {
		res.type = SWRW;
		res.roL  = roL*(res.p/pL + (gamma-1.)/(gamma+1.))/((gamma-1.)/(gamma+1.)*res.p/pL + 1.);
		res.roR  = roR*pow(res.p/pR, 1./gamma); 
	} else {
		res.type = SWSW;
		res.roL  = roL*(res.p/pL + (gamma-1.)/(gamma+1.))/((gamma-1.)/(gamma+1.)*res.p/pL + 1.);
		res.roR  = roR*(res.p/pR + (gamma-1.)/(gamma+1.))/((gamma-1.)/(gamma+1.)*res.p/pR + 1.);
	}
	//Тестирование (только без вакуума)
	/*double L =  - fL(p, roL, vL, pL);
	double R = fR(p, roR, vR, pR) + vR - vL;
	double delta = fabs((L-R)/0.5/(L+R));
	double LToro = - fL(res.p, roL, vL, pL);
	double RToro = fR(res.p, roR, vR, pR) + vR - vL;
	double deltaToro = fabs((LToro-RToro)/0.5/(LToro+RToro));*/
	return res;
}

double CSolver::fL(double p, double roL, double vL, double pL) {
	double gamma = 1.4;
	double f = 0.;
	if(p>pL) {
		double AL = 2./(gamma+1)/roL;
		double BL = (gamma-1.)/(gamma+1.)*pL;
		f = (p-pL) * sqrt(AL/(p+BL));
		return f;
	} else {
		double cL = sqrt(gamma*pL/roL);
		f = 2.*cL/(gamma-1.) * ( (pow(p/pL, (gamma-1.)/2./gamma)) - 1. );
		return f;	
	}
}

double CSolver::dfLdp(double p, double roL, double vL, double pL) {
	 double gamma = 1.4;
	double dfdp = 0.;
	if (p>pL) {
		double AL = 2./(gamma+1)/roL;
		double BL = (gamma-1.)/(gamma+1.)*pL;
		dfdp = sqrt(AL/(p+BL)) * (1. - (p-pL)/2./(p+BL));
		return dfdp;
	}
	else {
		double cL = sqrt(gamma*pL/roL);
		dfdp = cL/pL/gamma*pow(p/pL, -(gamma+1)/2./gamma); 
		return dfdp;
	}
}

double CSolver::fR(double p, double roR, double vR, double pR) {
	double gamma = 1.4;
	double f = 0.;
	if(p>pR) {
		double AR = 2./(gamma+1)/roR;
		double BR = (gamma-1.)/(gamma+1.)*pR;
		f = (p-pR) * sqrt(AR/(p+BR));
		return f;
	} else {
		double cR = sqrt(gamma*pR/roR);
		f = 2.*cR/(gamma-1.) * ( (pow(p/pR, (gamma-1.)/2./gamma)) - 1. );
		return f;
	}
}

double CSolver::dfRdp(double p, double roR, double vR, double pR) {
	double gamma = 1.4;
	double dfdp = 0.;
	if (p>pR) {
		double AR = 2./(gamma+1)/roR;
		double BR = (gamma-1.)/(gamma+1.)*pR;
		dfdp = sqrt(AR/(p+BR)) * (1. - (p-pR)/2./(p+BR));
		return dfdp;
	} else {
		double cR = sqrt(gamma*pR/roR);
		dfdp = cR/pR/gamma*pow(p/pR, -(gamma+1)/2./gamma); 
		return dfdp;
	}
;}

CVectorPrimitive CSolver::calcRPAnalyticalSolution(double roL, double vL, double pL, double roR, double vR, double pR, double x, double t){
	RPSolutionPrimitive res = solveRP(roL, vL, pL, roR, vR, pR);
	// V = (ro, v, p)T
	CVectorPrimitive V;
	double xi = x/t;
	double gamma = 1.4;
	double cL = 0., cR = 0.;
	if(roL!=0.) cL = sqrt(gamma*pL/roL);
	if(roR!=0.) cR = sqrt(gamma*pR/roR);
	double xiFront=0., xiHead=0., xiTail=0., xiHeadL=0., xiTailL=0., xiHeadR=0., xiTailR=0.;
	// Если вакуум
	if(res.type == VacRW) { 
		xiHead = vR + cR;
		xiTail = vR - 2.*cR/(gamma-1.);
		if(xi<=xiTail) {
			V.ro = 0.;
			V.v  = vR - 2.*cR/(gamma-1.);
			V.p  = 0.;
		} else if (xi<xiHead) {
			V.ro = roR*pow(2./(gamma+1.) - (gamma-1.)/(gamma+1.)/cR*(vR-xi), 2./(gamma-1.));
			V.v  = 2./(gamma+1)*(-cR + (gamma-1.)/2.*vR + xi);
			V.p  = pR*pow(2./(gamma+1.) - (gamma-1.)/(gamma+1.)/cR*(vR-xi), 2.*gamma/(gamma-1.)); 
		} else {
			V.ro = roR;
			V.v  = vR;
			V.p  = pR;
		}
		return V;
	}
	if(res.type == RWVac) {
		xiHead = vL - cL;
		xiTail = vL + 2.*cL/(gamma-1.);
		if(xi>=xiTail) {
			V.ro = 0.;
			V.v  = 0.;
			V.p  = 0.;
		} else if (xi>xiHead) {
			V.ro = roL*pow(2./(gamma+1.)+(gamma-1.)/(gamma+1.)/cL*(vL-xi), 2./(gamma-1.));
			V.v  = 2./(gamma+1)*(cL + (gamma-1.)/2.*vL + xi);
			V.p  = pL*pow(2./(gamma+1.)+(gamma-1.)/(gamma+1.)/cL*(vL-xi), 2.*gamma/(gamma-1.));
		} else {
			V.ro = roL;
			V.v = vL;
			V.p = pL;
		}
		return V;
	}
	if(res.type == RWVacRW) {
		xiHeadL = vL - cL;
		xiTailL = vL + 2.*cL/(gamma-1.);
		xiHeadR = vR + cR;
		xiTailR = vR - 2.*cR/(gamma-1.);
		if(xi<=xiHeadL) {
			V.ro = roL;
			V.v  = vL;
			V.p  = pL;
		} else if (xi<xiTailL) {
			V.ro = roL*pow(2./(gamma+1.)+(gamma-1.)/(gamma+1.)/cL*(vL-xi), 2./(gamma-1.));
			V.v  = 2./(gamma+1)*(cL + (gamma-1.)/2.*vL + xi);
			V.p  = pL*pow(2./(gamma+1.)+(gamma-1.)/(gamma+1.)/cL*(vL-xi), 2.*gamma/(gamma-1.));
		} else if (xi<=xiTailR) {
			V.ro = 0.;
			V.v  = 0.;
			V.p  = 0.;
		} else if (xi<xiHeadR) {
			V.ro = roR*pow(2./(gamma+1.) - (gamma-1.)/(gamma+1.)/cR*(vR-xi), 2./(gamma-1.));
			V.v  = 2./(gamma+1)*(-cR + (gamma-1.)/2.*vR + xi);
			V.p  = pR*pow(2./(gamma+1.) - (gamma-1.)/(gamma+1.)/cR*(vR-xi), 2.*gamma/(gamma-1.)); 
		} else {
			V.ro = roR;
			V.v  = vR;
			V.p  = pR;
		}
		return V;
	}
	double cLLocal = sqrt(gamma*res.p/res.roL), cRLocal = sqrt(gamma*res.p/res.roR);
	// Если не вакуум. Пусть точка слева от контактного разрыва (xiContact = res.v)
	if(xi<res.v) {
		if(res.type == SWSW || res.type == SWRW) { 
			xiFront = vL - cL*sqrt((gamma+1.)/2./gamma*res.p/pL + (gamma-1.)/2./gamma);
			if(xi<xiFront) {
				V.ro = roL;
				V.v  = vL;
				V.p  = pL;
			} else {
				V.ro = res.roL;
				V.v = res.v;
				V.p = res.p;
			}
		} else if (res.type == RWSW || res.type == RWRW) {
			xiHead = vL-cL;
			xiTail = res.v-cLLocal;
			if(xi<=xiHead) {
				V.ro = roL;
				V.v  = vL;
				V.p  = pL;
			} else if(xi>=xiTail) {
				V.ro = res.roL;
				V.v  = res.v;
				V.p  = res.p;
			} else {
				V.ro = roL*pow(2./(gamma+1.)+(gamma-1.)/(gamma+1.)/cL*(vL-xi), 2./(gamma-1.));
				V.v  = 2./(gamma+1)*(cL + (gamma-1.)/2.*vL + xi);
				V.p  = pL*pow(2./(gamma+1.)+(gamma-1.)/(gamma+1.)/cL*(vL-xi), 2.*gamma/(gamma-1.));
			}
		} 
	//Пусть точка справа от контактного разрыва (xiContact = res.v)
	} else {
		if(res.type == RWSW || res.type == SWSW) {
			xiFront = vR + cR*sqrt((gamma+1.)/2./gamma*res.p/pR + (gamma-1.)/2./gamma);
			if(xi>xiFront) {
				V.ro = roR;
				V.v  = vR;
				V.p  = pR;
			} else {
				V.ro = res.roR;
				V.v  = res.v;
				V.p  = res.p;
			}
		} else if(res.type == RWRW || res.type == SWRW) {
			xiHead = vR + cR;
			xiTail = res.v + cRLocal;
			if(xi >= xiHead) {
				V.ro = roR;
				V.v  = vR;
				V.p  = pR;
			} else if (xi <= xiTail) {
				V.ro = res.roR;
				V.v  = res.v;
				V.p  = res.p;
			} else {
				V.ro = roR*pow(2./(gamma+1.) - (gamma-1.)/(gamma+1.)/cR*(vR-xi), 2./(gamma-1.));
				V.v  = 2./(gamma+1)*(-cR + (gamma-1.)/2.*vR + xi);
				V.p  = pR*pow(2./(gamma+1.) - (gamma-1.)/(gamma+1.)/cR*(vR-xi), 2.*gamma/(gamma-1.));
			}
		}
	}
	return V;
}


// Точный Римановский солвер для (предположительно, произвольного) УРС в форме Ми-Грюнайзена [Miller, Puckett JoCP 1996]
CVectorPrimitive CSolver::calcRPExactMillerPuckett(CEOSMieGruneisen& eos, double roL, double uL, double pL, double roR, double uR, double pR, double x=0., double t=1.) {
	CVectorPrimitive V;
	const double gammaL = eos.getG(roL), gammaR = eos.getG(roR);
	const double c0 = eos.getc(eos.ro0, eos.getp0(eos.ro0));
	const double K0S = eos.ro0*c0*c0;
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




void CSolver::calcHydroStageMieGruneisen(CEOSMieGruneisen& eos, double t, double tau) {
	double roL = 0., uL = 0., eL = 0., pL = 0.,
		   roR = 0., uR = 0., eR = 0., pR = 0., 
		   E = 0.;
	Vector4 Fm, Fp;
	double _ro = 1000., _e = .637936e5, _p = 0.; 
	_p = eos.getp(_ro, _e); 
	cout << "_p=" << _p << endl;
	_e = eos.gete(_ro, _p);
	cout << "_e=" << _e << endl;

	double _x = 1., _G = eos.getG(998.2), _GPrime = eos.getGPrime(998.2), _L = _GPrime/_G + 1./_x + _G/_x;
	cout << "_L=" << _L << endl;


	// TODO: проверить на скорость выполнения операций, сравнить с реализацией через тип Vector4 -- если не медленнее, то в дальнейшем избавиться от Vector4 везде
	int nSize = ms.getSize();	
	double h = ms[1].x-ms[0].x;	
	int i=0;
	// Векторы W уже заполнены
	CVectorPrimitive res;
	// Потоки считаем по алгоритму решения задаче о распаде разрыва для УРС Ми-Грюнайзена из работы [Miller, Puckett]
	for(i=0; i<nSize; i++) {


		if(i==30) {
			double oo = 0;
		}




		Node& n=ms[i];			
		if(i!=0) {			
			roL = ms[i-1].W[0]; uL = ms[i-1].W[1]/roL; eL = ms[i-1].W[2]/roL - .5*uL*uL; pL = eos.getp(roL, eL);
			roR = ms[i].W[0];   uR = ms[i].W[1]/roR;   eR = ms[i].W[2]/roR - .5*uR*uR;   pR = eos.getp(roR, eR);						
		} else { 
		    // Левое граничное условие -- прозрачная граница
			roR = ms[i].W[0];   uR = ms[i].W[1]/roR;   eR = ms[i].W[2]/roR - .5*uR*uR;   pR = eos.getp(roR, eR);			
			roL = roR; uL = uR; pL = pR;			
		}
		res = calcRPExactMillerPuckett(eos, roL, uL, pL, roR, uR, pR);
		if(res.ro!=0.) E = eos.gete(res.ro, res.p) + .5*res.v*res.v; else E = 0.;
		Fm = Vector4(res.ro*res.v, res.ro*res.v*res.v + res.p, res.v*(res.ro*E+res.p), 0.);
		if(i!=nSize-1) {
			
			roL = ms[i].W[0];   uL = ms[i].W[1]/roL;   eL = ms[i].W[2]/roL - .5*uL*uL;   pL = eos.getp(roL, eL);						
			roR = ms[i+1].W[0]; uR = ms[i+1].W[1]/roR; eR = ms[i+1].W[2]/roR - .5*uR*uR; pR = eos.getp(roR, eR);			
		} else {
			// Правое граничное условие -- прозрачная граница
			roL = ms[i].W[0];   uL = ms[i].W[1]/roL;   eL = ms[i].W[2]/roL - .5*uL*uL;   pL = eos.getp(roL, eL);						
			roR = roL; uR = uL; pR = pL;					
		}
		res = calcRPExactMillerPuckett(eos, roL, uL, pL, roR, uR, pR);
		if(res.ro!=0.) E = eos.gete(res.ro, res.p) + .5*res.v*res.v; else E = 0.;		
		Fp = Vector4(res.ro*res.v, res.ro*res.v*res.v + res.p, res.v*(res.ro*E+res.p), 0.);
		n.W_temp = n.W - tau/h*(Fp-Fm);
	}
	for(i=0; i<nSize; i++) 
	{
		Node& n=ms[i];
		n.W[0] = n.W_temp[0];
		n.W[1] = n.W_temp[1];
		n.W[2] = n.W_temp[2];
		n.ro = n.W[0];
		n.v  = n.W[1]/n.ro;
		n.e  = n.W[2]/n.ro - 0.5*n.v*n.v;
		if(n.e < 0. ) {
			cout << "Warning: negative internal energy in cell i=" << i << ". Probably numerical instability." << endl;
			int qq = 0;
		}
		n.p = eos.getp(n.ro, n.e);
		n.C = eos.getc(n.ro, n.p);	
	}
	cout << "calcHydroStageMieGruneisen(): done!" << endl;
}
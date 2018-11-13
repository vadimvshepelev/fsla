
#include "solver.h"

void CSolver::calcHeatStage(double t, double tau) {
	int i = 0, itNum = 0, size = ms.getSize();;
	const double eps = 0.01; 
	double kappa_plus = 0., kappa_minus = 0., ro_plus = 0., ro_minus = 0.;
	EOS &eos = task.getEOS();
	CField ms_temp, ms_temp_temp;
	ms_temp.initData(&task);
	ms_temp_temp.initData(&task);
	for(i=0; i<size; i++)
	{
		     ms_temp[i].te = ms[i].te;
		ms_temp_temp[i].te = ms[i].te;
		
			 ms_temp[i].Z  = ms[i].Z;
		ms_temp_temp[i].Z  = ms[i].Z;

	}


	///////////////
//	printf("1 oh yeah!\n");
//	fflush(stdout);
	///////////////	

	// Столбцы коэффициентов системы. Система имеет вид:
	// A...-C...+B...=-F

	double* A  = new double[size];
	double* B  = new double[size];
	double* C  = new double[size];
	double* F  = new double[size];

	double *alpha = new double[size];
	double *beta  = new double[size];

	///////////////
//	printf("2 oh yeah!\n");
//	fflush(stdout);
	///////////////

	for(;;)
	{
		// Присваиваем значения коэффициентам 
		// (решается уравнение теплопроводности с помощью неявной схемы)

		for(i=0; i<size; i++)
		{
			if(i==0)
			{
				kappa_minus = ms[i].kappa;
				   ro_minus = ms[i].ro;
			}
			else
			{				
				kappa_minus = (ms[i].kappa + ms[i-1].kappa) / 2.0;
				   ro_minus = (ms[i].ro    + ms[i-1].ro   ) / 2.0;
			}
			
			if(i==size-1)
			{

				kappa_plus  = ms[i].kappa;
				   ro_plus  = ms[i].ro;
			}
			else
			{
				kappa_plus  = (ms[i].kappa + ms[i+1].kappa) / 2.0;
				   ro_plus  = (ms[i].ro    + ms[i+1].ro   ) / 2.0;
			}


//////////////  DEBUG ////////////////////////////////////////////////////
// Здесь делаем так, чтобы тепло не проходило через границу "стекло-металл"
// Возможно, это немного искусственно.
// Если наводить потом красоту, надо будет этот момент тоже сделать красивым
			
			if(i==i_pulse_min+1)
			{
				//////////
				double te_minus  = ms[i-1].te;
				double te_center = ms[i].te;
				double te_plus   = ms[i+1].te;
				//////////
        
				kappa_minus = 0.0;				   
			}
			
///////////////////////////////////////////////////////////////////////////



			A[i] = tau * kappa_minus / (ms[i].dm*ms[i].dm) * ms[i].ro * ro_minus;
			B[i] = tau * kappa_plus  / (ms[i].dm*ms[i].dm) * ms[i].ro * ro_plus;
			C[i] = ms[i].ce + A[i] + B[i];
// 14.06.2013
			double src = 0.0;
			double expX=0.0, expT=0.0;

			if(task.getSourceFlag() == 1)
			{
				expX = exp(-(ms[i].x)/deltaSkin);
				expT = exp(-(t)*(t)/tauPulse/tauPulse);
			}
			else if ((task.getSourceFlag() == 2) && (i>i_pulse_min))
			{
				expX = exp(-(ms[i].x - ms[i_pulse_min].x)/deltaSkin);
				expT = exp(-(t)*(t)/tauPulse/tauPulse);
			}
			else if ((task.getSourceFlag() == 3) && (i>i_pulse_min))
			{
				expX = exp(-(ms[i].x - ms[i_pulse_min].x)/deltaSkin)*sqrt(3.14159);
				if( (t>=0.0) && (t<tauPulse) )
					expT = 1.0;
				else
					expT = 0.0;
				
				// expT = exp(-(t)*(t)/tauPulse/tauPulse);
			}

			src=fluence/sqrt(3.14159)/deltaSkin/tauPulse*expX*expT;

			F[i] = ms[i].ce * ms[i].te + src * tau;
		}	

		sweepTe(ms_temp, ms_temp_temp, A, B, C, F, alpha, beta, size);

		for(i=0; i<size; i++)
		{
			Node &n = ms[i];

			n.ce    = eos.getce   (n.ro, ms_temp[i].te);
			n.kappa = eos.getkappa(n.ro, n.ti, ms_temp[i].te);
		}
	
		if(compTe(ms_temp, ms_temp_temp) < eps) break;

		itNum++;

		if(itNum == maxIt)
		{
			printf("Divergence=%f, iteration number %d\n", compTe(), itNum);
			printf("Error. Too many iterations for convergence.\n")	;
			cin.get();
			exit(1);
		}
	}

/*	if(t>5.e-12 && itNum <3)
	{
		double _Kur = getzKur();
		setzKur(1.01*_Kur);
	}
*/
	printf("calcHeatStage() complete: %d iterations of Te\n", itNum+1);
	for(i=0; i<size; i++) {
		Node &n = ms[i];
		n.te    = ms_temp[i].te;
		n.kappa = eos.getkappa(n.ro ,n.ti, n.te);
		n.ce	= eos.getce(n.ro, n.te); 
		n.C = eos.getC(n.ro, n.ti, n.te);
		n.pe = eos.getpe(n.ro, n.ti, n.te);
		n.p  = n.pi + n.pe;
		n.ee = eos.getee(n.ro, n.ti, n.te);
		n.e  = n.ei + n.ee;
	}
	delete[] A;
	delete[] B;
	delete[] C;
	delete[] F;
	delete[] alpha;
	delete[] beta;
}

void CSolver::calcHeatStage5LayersSi(double t, double tau) {
/*	int i = 0, itNum = 0, size = ms.getSize();
	double eps = 0.01, kappa_plus=0., kappa_minus=0., ro_plus=0., ro_minus=0., src = 0.;
	double x = 0., dx = 0.;
	EOS &eosCr = task.getEOSCr();
	EOS &eosAu = task.getEOSAu();
	EOS &eosSi = task.getEOSSi();
	MatterState ms_temp, ms_temp_temp;
	ms_temp.initData(&task); ms_temp_temp.initData(&task);
	for(i=0; i<size; i++) {
		     ms_temp[i].te = ms[i].te;
		ms_temp_temp[i].te = ms[i].te;
	}
	// Столбцы коэффициентов системы. Система имеет вид:
	// A...-C...+B...=-F
	double* A  = new double[size];
	double* B  = new double[size];
	double* C  = new double[size];
	double* F  = new double[size];
	double *alpha = new double[size];
	double *beta  = new double[size];
	for(;;) {
		// Присваиваем значения коэффициентам 
		// (решается уравнение теплопроводности с помощью неявной схемы)
		for(i=0; i<size; i++) {
			if(i==0) {
				kappa_minus = ms[i].kappa;
				   ro_minus = ms[i].ro;
			} else {				
				kappa_minus = (ms[i].kappa + ms[i-1].kappa) / 2.0;
				   ro_minus = (ms[i].ro    + ms[i-1].ro   ) / 2.0;
			}		
			if(i==size-1) {
				kappa_plus  = ms[i].kappa;
				   ro_plus  = ms[i].ro;
			} else {
				kappa_plus  = (ms[i].kappa + ms[i+1].kappa) / 2.0;
				   ro_plus  = (ms[i].ro    + ms[i+1].ro   ) / 2.0;
			}
			A[i] = tau * kappa_minus / (ms[i].dm*ms[i].dm) * ms[i].ro * ro_minus;
			B[i] = tau * kappa_plus  / (ms[i].dm*ms[i].dm) * ms[i].ro * ro_plus;
			C[i] = ms[i].ce + A[i] + B[i];		
			*/
		    //!!! Здесь самое интересное -- расчет флюенса
			/*if(task.getSourceFlag() == 1) {
				expX = exp(-(ms[i].x)/deltaSkin);
				expT = exp(-(t)*(t)/tauPulse/tauPulse);
			} else if ((task.getSourceFlag() == 2) && (i>=i_pulse_min)) { // '>' или '>=' ??? -- скорее всего, '>=', проверить 
				expX = exp(-(ms[i].x - ms[i_pulse_min].x)/deltaSkin);
				expT = exp(-(t)*(t)/tauPulse/tauPulse);
			} else if ((task.getSourceFlag() == 3) && (i>=i_pulse_min))	{
				expX = exp(-(ms[i].x - ms[i_pulse_min].x)/deltaSkin)*sqrt(3.14159);
				if( (t>=0.0) && (t<tauPulse) )
					expT = 1.0;
				else
					expT = 0.0;							
			}*/
		/*	if (t > 10.e-15) {
				src = 0.;
			} else {
				x = .5 * (ms[i].x + ms[i+1].x); dx = ms[i+1].x - ms[i].x; 
				if ( x > 250.e-9) 
					src = fluence/30.;///2330.;
				else if (x > 200.e-9)
					src = fluence/3.;///7190.;
				else if (x > 150.e-9)
					src = fluence;///19320.;
				else if (x > 100.e-9)
					src = fluence/3.;///7190.;
				else if (x > 50.e-9)
					src = fluence;///19320.;
				else
					src = fluence/3.;///7190.;
			}
			F[i] = ms[i].ce * ms[i].te + src * tau;
		}	
		sweepTe(ms_temp, ms_temp_temp, A, B, C, F, alpha, beta, size);
		for(i=0; i<size; i++) {
			Node &n = ms[i]; x = .5*(ms[i].x + ms[i+1].x);
			if ( x > 250.e-9) {
				n.ce    = eosSi.getce   (n.ro, ms_temp[i].te, n.Z);
				n.kappa = eosSi.getkappa(n.ro, n.ti, ms_temp[i].te, n.Z);
			} else if (x > 200.e-9) {
				n.ce    = eosCr.getce   (n.ro, ms_temp[i].te, n.Z);
				n.kappa = eosCr.getkappa(n.ro, n.ti, ms_temp[i].te, n.Z);
			} else if (x > 150.e-9) {
				n.ce    = eosAu.getce   (n.ro, ms_temp[i].te, n.Z);
				n.kappa = eosAu.getkappa(n.ro, n.ti, ms_temp[i].te, n.Z);
			} else if (x > 100.e-9) {
				n.ce    = eosCr.getce   (n.ro, ms_temp[i].te, n.Z);
				n.kappa = eosCr.getkappa(n.ro, n.ti, ms_temp[i].te, n.Z);
			} else if (x > 50.e-9) {
				n.ce    = eosAu.getce   (n.ro, ms_temp[i].te, n.Z);
				n.kappa = eosAu.getkappa(n.ro, n.ti, ms_temp[i].te, n.Z);
			} else {
				n.ce    = eosCr.getce   (n.ro, ms_temp[i].te, n.Z);
				n.kappa = eosCr.getkappa(n.ro, n.ti, ms_temp[i].te, n.Z);
			}
		}
		if(compTe(ms_temp, ms_temp_temp) < eps) break;
		itNum++;
		if(itNum == maxIt)
		{
			printf("Divergence=%f, iteration number %d\n", compTe(), itNum);
			printf("Error. Too many iterations for convergence.\n")	;
			cin.get();
			exit(1);
		}
	}
	printf("calcHeatStageGlass() complete: %d iterations of Te\n", itNum+1);
	for(i=0; i<size; i++)	{
		Node &n = ms[i];
		n.te    = ms_temp[i].te;
		x = .5*(ms[i].x + ms[i+1].x);
			if ( x > 250.e-9) {
				n.C  = eosSi.getC(n.ro, n.ti, n.te, n.Z);
				n.ce = eosSi.getce(n.ro, n.te, n.Z); 
				n.ci = eosSi.getci(n.ro, n.ti);
				n.pe = eosSi.getpe(n.ro, n.ti, n.te, n.Z);
				n.ee = eosSi.getee(n.ro, n.ti, n.te, n.Z);
				n.kappa   = eosSi.getkappa(n.ro ,n.ti, n.te, n.Z);
				n.Alphaei = eosSi.getAlpha(n.ro, n.ti, n.te);
			} else if (x > 200.e-9) {
				n.C  = eosCr.getC(n.ro, n.ti, n.te, n.Z);
				n.ce = eosCr.getce(n.ro, n.te, n.Z); 
				n.ci = eosCr.getci(n.ro, n.ti);
				n.pe = eosCr.getpe(n.ro, n.ti, n.te, n.Z);
				n.ee = eosCr.getee(n.ro, n.ti, n.te, n.Z);
				n.kappa   = eosCr.getkappa(n.ro ,n.ti, n.te, n.Z);
				n.Alphaei = eosCr.getAlpha(n.ro, n.ti, n.te);
			} else if (x > 150.e-9) {
				n.C  = eosAu.getC(n.ro, n.ti, n.te, n.Z);
				n.ce = eosAu.getce(n.ro, n.te, n.Z); 
				n.ci = eosAu.getci(n.ro, n.ti);
				n.pe = eosAu.getpe(n.ro, n.ti, n.te, n.Z);
				n.ee = eosAu.getee(n.ro, n.ti, n.te, n.Z);
				n.kappa   = eosAu.getkappa(n.ro ,n.ti, n.te, n.Z);
				n.Alphaei = eosAu.getAlpha(n.ro, n.ti, n.te);
			} else if (x > 100.e-9) {
				n.C  = eosCr.getC(n.ro, n.ti, n.te, n.Z);
				n.ce = eosCr.getce(n.ro, n.te, n.Z); 
				n.ci = eosCr.getci(n.ro, n.ti);
				n.pe = eosCr.getpe(n.ro, n.ti, n.te, n.Z);
				n.ee = eosCr.getee(n.ro, n.ti, n.te, n.Z);
				n.kappa   = eosCr.getkappa(n.ro ,n.ti, n.te, n.Z);
				n.Alphaei = eosCr.getAlpha(n.ro, n.ti, n.te);
			} else if (x > 50.e-9) {
				n.C  = eosAu.getC(n.ro, n.ti, n.te, n.Z);
				n.ce = eosAu.getce(n.ro, n.te, n.Z); 
				n.ci = eosAu.getci(n.ro, n.ti);
				n.pe = eosAu.getpe(n.ro, n.ti, n.te, n.Z);
				n.ee = eosAu.getee(n.ro, n.ti, n.te, n.Z);
				n.kappa   = eosAu.getkappa(n.ro ,n.ti, n.te, n.Z);
				n.Alphaei = eosAu.getAlpha(n.ro, n.ti, n.te);
			} else {
				n.C  = eosCr.getC(n.ro, n.ti, n.te, n.Z);
				n.ce = eosCr.getce(n.ro, n.te, n.Z); 
				n.ci = eosCr.getci(n.ro, n.ti);
				n.pe = eosCr.getpe(n.ro, n.ti, n.te, n.Z);
				n.ee = eosCr.getee(n.ro, n.ti, n.te, n.Z);
				n.kappa   = eosCr.getkappa(n.ro ,n.ti, n.te, n.Z);
				n.Alphaei = eosCr.getAlpha(n.ro, n.ti, n.te);
			}
		n.p  = n.pi + n.pe;
		n.e  = n.ei + n.ee;
	}
	delete[] A; delete[] B; delete[] C; delete[] F; delete[] alpha; delete[] beta;*/
}


int CSolver::calcHeatStageGlass(double t, double tau) {
	int i=0;
	int itNum = 0;
	const double eps = 0.01;
	double kappa_plus=0., kappa_minus=0.;
	double ro_plus=0., ro_minus=0.;
	int size = ms.getSize();
	EOS &eos = task.getEOS();
	EOS &eosGlass = task.getEOSGlass();
	int nBound = task.getZone(0).n;
	CField ms_temp, ms_temp_temp;
	ms_temp.initData(&task);
	ms_temp_temp.initData(&task);
	// Для поглощения
	const int N = task.getZone(0).n;
	const double L = task.getZone(0).l;
	double dx = L/N;
	double _x = 0.;
	for(i=0; i<size; i++) {
		     ms_temp[i].te = ms[i].te;
		ms_temp_temp[i].te = ms[i].te;
	}
	// Столбцы коэффициентов системы. Система имеет вид:
	// A...-C...+B...=-F

	double* A  = new double[size];
	double* B  = new double[size];
	double* C  = new double[size];
	double* F  = new double[size];

	double *alpha = new double[size];
	double *beta  = new double[size];

	//cout << "Heat:";
	for(;;)
	{
		// Присваиваем значения коэффициентам 
		// (решается уравнение теплопроводности с помощью неявной схемы)

		for(i=0; i<size; i++)
		{
			if(i==0)
			{
				kappa_minus = ms[i].kappa;
				   ro_minus = ms[i].ro;
			}
			else
			{				
				kappa_minus = (ms[i].kappa + ms[i-1].kappa) / 2.0;
				   ro_minus = (ms[i].ro    + ms[i-1].ro   ) / 2.0;
			}
			
			if(i==size-1)
			{

				kappa_plus  = ms[i].kappa;
				   ro_plus  = ms[i].ro;
			}
			else
			{
				kappa_plus  = (ms[i].kappa + ms[i+1].kappa) / 2.0;
				   ro_plus  = (ms[i].ro    + ms[i+1].ro   ) / 2.0;
			}
			A[i] = tau * kappa_minus / (ms[i].dm*ms[i].dm) * ms[i].ro * ro_minus;
			B[i] = tau * kappa_plus  / (ms[i].dm*ms[i].dm) * ms[i].ro * ro_plus;
			C[i] = ms[i].ce + A[i] + B[i];
			double src = 0.0;
			double expX=0.0, expT=0.0;

			if(task.getSourceFlag() == 1) 	{
				//expX = exp(-(ms[i].x)/deltaSkin);
				// Здесь для рентгеновского лазера не слишком точный учет поглощения -- на самом деле, нужно брать начальную координату
				// x с номером i. Также хочу сразу внести уточнение в интеграл -- брать не левую границу примоугольника, а среднее
				// арифметическое между двумя

				if(i >= N) 
					expX = 0.; 
				else {
					_x = (0.5 + i)*dx;
					expX = exp(-(_x/deltaSkin));
				}					
				expT = exp(-(t)*(t)/tauPulse/tauPulse);
			}
			else if ((task.getSourceFlag() == 2) && (i>=i_pulse_min)) // '>' или '>=' ??? -- скорее всего, '>=', проверить
			{
				expX = exp(-(ms[i].x - ms[i_pulse_min].x)/deltaSkin);
				expT = exp(-(t)*(t)/tauPulse/tauPulse);
			}
			else if ((task.getSourceFlag() == 3) && (i>=i_pulse_min))
			{
				expX = exp(-(ms[i].x - ms[i_pulse_min].x)/deltaSkin)*sqrt(3.14159);
				if( (t>=0.0) && (t<tauPulse) )
					expT = 1.0;
				else
					expT = 0.0;
				
				// expT = exp(-(t)*(t)/tauPulse/tauPulse);
			} 

			src=fluence/sqrt(3.14159)/deltaSkin/tauPulse*expX*expT;

			F[i] = ms[i].ce * ms[i].te + src * tau;
		}	

		sweepTe(ms_temp, ms_temp_temp, A, B, C, F, alpha, beta, size);

		for(i=0; i<size; i++) {
			Node &n = ms[i];
			if ( task.getSourceFlag() == 2 ) {
				if(i>=nBound) {
					n.ce    = eos.getce   (n.ro, ms_temp[i].te);
					n.kappa = eos.getkappa(n.ro, n.ti, ms_temp[i].te);
				} else {
					n.ce    = eosGlass.getce   (n.ro, ms_temp[i].te);
					n.kappa = eosGlass.getkappa(n.ro, n.ti, ms_temp[i].te);
				}
			} else if(i<nBound) {
				n.ce    = eos.getce   (n.ro, ms_temp[i].te);
				n.kappa = eos.getkappa(n.ro, n.ti, ms_temp[i].te);
			} else {
				n.ce    = eosGlass.getce   (n.ro, ms_temp[i].te);
				n.kappa = eosGlass.getkappa(n.ro, n.ti, ms_temp[i].te);
			}
		}
		if(compTe(ms_temp, ms_temp_temp) < eps) break;
		itNum++;
		if(itNum == maxIt)
		{
			printf("Divergence=%f, iteration number %d\n", compTe(), itNum);
			printf("Error. Too many iterations for convergence.\n")	;
			cin.get();
			exit(1);
		}
	}

/*	if(t>5.e-12 && itNum <3)
	{
		double _Kur = getzKur();
		setzKur(1.01*_Kur);
	}
*/
	//printf("calcHeatStageGlass() complete: %d iterations of Te\n", itNum+1);
//	cout << itNum+1 << "it ";
	for(i=0; i<size; i++)	{
		Node &n = ms[i];
		n.te    = ms_temp[i].te;
		if ( task.getSourceFlag() == 2 ) {
				if(i>=nBound) {
					n.kappa = eos.getkappa(n.ro ,n.ti, n.te);
					n.Alphaei = eos.getAlpha(n.ro, n.ti, n.te);
					n.C = eos.getC(n.ro, n.ti, n.te);
					n.ce	= eos.getce(n.ro, n.te); 
					n.pe = eos.getpe(n.ro, n.ti, n.te);
					n.ee = eos.getee(n.ro, n.ti, n.te);
				} else {
					n.kappa = eosGlass.getkappa(n.ro ,n.ti, n.te);
					n.Alphaei = eosGlass.getAlpha(n.ro, n.ti, n.te);
					n.C = eos.getC(n.ro, n.ti, n.te);
					n.ce	= eosGlass.getce(n.ro, n.te); 
					n.pe = eosGlass.getpe(n.ro, n.ti, n.te);
					n.ee = eosGlass.getee(n.ro, n.ti, n.te);
				}
		} else if(i<nBound) {
			n.kappa = eos.getkappa(n.ro ,n.ti, n.te);
			n.Alphaei = eos.getAlpha(n.ro, n.ti, n.te);
			n.C = eos.getC(n.ro, n.ti, n.te);
			n.ce	= eos.getce(n.ro, n.te); 
			n.pe = eos.getpe(n.ro, n.ti, n.te);
			n.ee = eos.getee(n.ro, n.ti, n.te);
		} else {
			n.kappa = eosGlass.getkappa(n.ro ,n.ti, n.te);
			n.Alphaei = eosGlass.getAlpha(n.ro, n.ti, n.te);
			n.C = eos.getC(n.ro, n.ti, n.te);
			n.ce	= eosGlass.getce(n.ro, n.te); 
			n.pe = eosGlass.getpe(n.ro, n.ti, n.te);
			n.ee = eosGlass.getee(n.ro, n.ti, n.te);
		}
		n.p  = n.pi + n.pe;
		n.e  = n.ei + n.ee;
	}
	delete[] A; delete[] B; delete[] C; delete[] F; delete[] alpha; delete[] beta;
	return itNum+1;
}

int CSolver::calcIonicHeatStageGlass(double t, double tau) {
	unsigned int i=0;
	unsigned int itNum = 0;
	const double eps = 0.01;
	double kappa_plus=0., kappa_minus=0.;
	double ro_plus=0., ro_minus=0.;
	unsigned int size = ms.getSize();
	EOS &eos = task.getEOS();
	EOS &eosGlass = task.getEOSGlass();
	unsigned int nBound = task.getZone(0).n;
	CField ms_temp, ms_temp_temp;
	ms_temp.initData(&task);
	ms_temp_temp.initData(&task);
	// Для поглощения
	const unsigned int N = task.getZone(0).n;
	const double L = task.getZone(0).l;
	double dx = L/N;
	double _x = 0.;
	for(i=0; i<size; i++) {
		     ms_temp[i].ti = ms[i].ti;
		ms_temp_temp[i].ti = ms[i].ti;
	}
	// Столбцы коэффициентов системы. Система имеет вид:
	// A...-C...+B...=-F
	double* A  = new double[size]; double* B  = new double[size]; 	double* C  = new double[size]; 	double* F  = new double[size];
	double *alpha = new double[size]; double *beta  = new double[size];	
	if(task.getSourceFlag() != SourceType::SrcGlass) 	{
				cerr << "CSolver::calcIonicHeatStage() reports error: only glass problems can have ionc heat stage now." << endl;
				exit(1);
			} 
	for(;;)
	{
		// Присваиваем значения коэффициентам 
		// (решается уравнение теплопроводности с помощью неявной схемы)
		/////////////////////////////////////////////////////////////////////////////////////////
		// !!! DEBUG
		// Пусть в стекле в kappa лежит у нас значение ионной теплопроводности, а не электронной
		// //////////////////////////////////////////////////////////////////////////////////////
		for(i=0; i<size; i++) {
			if(i==0) {
				kappa_minus = ms[i].kappa;
				   ro_minus = ms[i].ro;
			}
			else {				
				kappa_minus = (ms[i].kappa + ms[i-1].kappa) / 2.0;
				   ro_minus = (ms[i].ro    + ms[i-1].ro   ) / 2.0;
			} 			
			if(i==size-1) {

				kappa_plus  = ms[i].kappa;
				   ro_plus  = ms[i].ro;
			} else {
				kappa_plus  = (ms[i].kappa + ms[i+1].kappa) / 2.0;
				   ro_plus  = (ms[i].ro    + ms[i+1].ro   ) / 2.0;
			}
			A[i] = tau * kappa_minus / (ms[i].dm*ms[i].dm) * ms[i].ro * ro_minus;
			B[i] = tau * kappa_plus  / (ms[i].dm*ms[i].dm) * ms[i].ro * ro_plus;
			C[i] = ms[i].ci + A[i] + B[i];				
			F[i] = ms[i].ci * ms[i].ti;
		}	
		sweepTi(ms_temp, ms_temp_temp, A, B, C, F, alpha, beta, size);
		for(i=0; i<size; i++) {
			Node &n = ms[i];
			if(i>=nBound) {
				n.ci    = eos.getci(n.ro, ms_temp[i].ti);
				n.kappa = eos.getkappa(n.ro, ms_temp[i].ti, n.te);
			} else {
				n.ci    = eosGlass.getci   (n.ro, ms_temp[i].ti);
				n.kappa = eosGlass.getkappa(n.ro, ms_temp[i].ti, n.te);				
			}
		}
		if(compTi(ms_temp, ms_temp_temp) < eps) break;
		itNum++;
		if(itNum == maxIt) {
			cerr << "CSolver::calcIonicHeatStage() reports error: too many iterations for convergence." << endl;
			cerr << "Divergence=" << compTi() << ", iteration number " << itNum << endl;			
			exit(1);
		}		
	}
	for(i=0; i<size; i++) {
		Node &n = ms[i];
		n.ti = ms_temp[i].ti;
		if(i>=nBound) {
			n.kappa = eos.getkappa(n.ro ,n.ti, n.te);
			n.Alphaei = eos.getAlpha(n.ro, n.ti, n.te);
			n.C = eos.getC(n.ro, n.ti, n.te);
			n.ci	= eos.getci(n.ro, n.ti); 
			n.pi = eos.getpi(n.ro, n.ti);
			n.pe = eos.getpe(n.ro, n.ti, n.te);
			n.ei = eos.getei(n.ro, n.ti);
			n.ee = eos.getee(n.ro, n.ti, n.te);
		} else {
			n.kappa = eosGlass.getkappa(n.ro ,n.ti, n.te);
			n.Alphaei = eosGlass.getAlpha(n.ro, n.ti, n.te);
			n.C = eos.getC(n.ro, n.ti, n.te);
			n.ci	= eosGlass.getci(n.ro, n.ti); 
			n.pi = eosGlass.getpi(n.ro, n.ti);
			n.pe = eosGlass.getpe(n.ro, n.ti, n.te);
			n.ei = eosGlass.getei(n.ro, n.ti);
			n.ee = eosGlass.getee(n.ro, n.ti, n.te);
		}		
		n.p  = n.pi + n.pe;
		n.e  = n.ei + n.ee;
	}
	////////////////////////////////////////////////////////////////////////////////////////
	delete[] A; delete[] B; delete[] C; delete[] F; delete[] alpha; delete[] beta;
	return itNum+1;
}

void CSolver::calcHeatStageSpallation(double t, double tau) {
	int i = 0, itNum = 0, iSpall = getSpallCellNum(), size = ms.getSize();
	double eps = 0.01, kappa_plus = 0., kappa_minus = 0., ro_plus=0., ro_minus=0.;
	EOS &eos = task.getEOS();
	CField ms_temp, ms_temp_temp;
	ms_temp.initData(&task);
	ms_temp_temp.initData(&task);
	for(i=0; i<size; i++) {
		     ms_temp[i].te = ms[i].te;
		ms_temp_temp[i].te = ms[i].te;
			 ms_temp[i].Z  = ms[i].Z;
		ms_temp_temp[i].Z  = ms[i].Z;
	}
	// Столбцы коэффициентов системы. Система имеет вид:
	// A...-C...+B...=-F
	double* A  = new double[size]; double* B  = new double[size]; double* C  = new double[size]; double* F  = new double[size];
	double *alpha = new double[size]; double *beta  = new double[size];
	for(;;) {
		// Присваиваем значения коэффициентам (решается уравнение теплопроводности с помощью неявной схемы)
		for(i=0; i<size; i++) {
			if(i==0) {
				kappa_minus = ms[i].kappa; 
				   ro_minus = ms[i].ro;
			} else {				
				kappa_minus = (ms[i].kappa + ms[i-1].kappa) / 2.0;
				   ro_minus = (ms[i].ro    + ms[i-1].ro   ) / 2.0;
			}
			if(i==size-1) {
				kappa_plus  = ms[i].kappa;
				   ro_plus  = ms[i].ro;
			} else {
				kappa_plus  = (ms[i].kappa + ms[i+1].kappa) / 2.0;
				   ro_plus  = (ms[i].ro    + ms[i+1].ro   ) / 2.0;
			}
			if(i == iSpall) {
				kappa_plus  = 0.;
				kappa_minus = 0.;
			}
			if(i == iSpall-1) {
				kappa_plus  = 0.;
			}
			if(i == iSpall+1) {
				kappa_minus = 0.;
			}


//////////////  DEBUG ////////////////////////////////////////////////////
// Здесь делаем так, чтобы тепло не проходило через границу "стекло-металл"
// Возможно, это немного искусственно.
// Если наводить потом красоту, надо будет этот момент тоже сделать красивым
// 31.07.2013 -- я, честно, говоря, не понимаю, как этот фрагмент кода закрывает
// теплопередачу -- типа просто делает в формуле для неявной схемы температуру
// в точке i_pulse_min+1 независимой от температуры в точке i_pulse_min. Хм, 
// хотя может это и работает. Сейчас попробую сделать там тест без гидродинамики.
// Возможно этот небольшой прирост температуры возникает из-за гидродинамического 
// переноса энергии, а не теплового. -- Это, кстати, должно работать! -- Смешно, 
// но добавление этого пункта ничего особо не меняет, т.к. все равно прогонка 
//действует только начиная с i_pulse_min. Поэтому закомментим этот момент.
/*			if(i==i_pulse_min+1) {
				//////////
				double te_minus  = ms[i-1].te; double te_center = ms[i].te; double te_plus   = ms[i+1].te;
				//////////
				kappa_minus = 0.0;				   
			}*/
			
///////////////////////////////////////////////////////////////////////////



			A[i] = tau * kappa_minus / (ms[i].dm*ms[i].dm) * ms[i].ro * ro_minus;
			B[i] = tau * kappa_plus  / (ms[i].dm*ms[i].dm) * ms[i].ro * ro_plus;
			C[i] = ms[i].ce + A[i] + B[i];
// 14.06.2013
			double src = 0.0;
			double expX=0.0, expT=0.0;

			if(task.getSourceFlag() == 1)
			{
				expX = exp(-(ms[i].x)/deltaSkin);
				expT = exp(-(t)*(t)/tauPulse/tauPulse);
			}
			else if ((task.getSourceFlag() == 2) && (i>i_pulse_min))
			{
				expX = exp(-(ms[i].x - ms[i_pulse_min].x)/deltaSkin);
				expT = exp(-(t)*(t)/tauPulse/tauPulse);
			}
			else if ((task.getSourceFlag() == 3) && (i>i_pulse_min))
			{
				expX = exp(-(ms[i].x - ms[i_pulse_min].x)/deltaSkin)*sqrt(3.14159);
				if( (t>=0.0) && (t<tauPulse) )
					expT = 1.0;
				else
					expT = 0.0;
				
				// expT = exp(-(t)*(t)/tauPulse/tauPulse);
			}

			if (i!=iSpall) 
				src=fluence/sqrt(3.14159)/deltaSkin/tauPulse*expX*expT;
			else
				src=0.;
			F[i] = ms[i].ce * ms[i].te + src * tau;
		}	

		sweepTe(ms_temp, ms_temp_temp, A, B, C, F, alpha, beta, size);

		for(i=0; i<size; i++)
		{
			Node &n = ms[i];

			n.ce    = eos.getce   (n.ro, ms_temp[i].te);
			n.kappa = eos.getkappa(n.ro, n.ti, ms_temp[i].te);
		}
	
		if(compTe(ms_temp, ms_temp_temp) < eps) break;

		itNum++;

		if(itNum == maxIt)
		{
			printf("Divergence=%f, iteration number %d\n", compTe(), itNum);
			printf("Error. Too many iterations for convergence.\n")	;
			cin.get();
			exit(1);
		}
	}
	printf("calcHeatStage() complete: %d iterations of Te\n", itNum+1);
	for(i=0; i<size; i++)
	{
		Node &n = ms[i];

		n.te    = ms_temp[i].te;
		n.kappa = eos.getkappa(n.ro ,n.ti, n.te);
		n.ce	= eos.getce(n.ro, n.te); 
		n.pe = eos.getpe(n.ro, n.ti, n.te);
		n.p  = n.pi + n.pe;
		n.ee = eos.getee(n.ro, n.ti, n.te);
		n.e  = n.ei + n.ee;
	}
	delete[] A;
	delete[] B;
	delete[] C;
	delete[] F;

	delete[] alpha;
	delete[] beta;
}

// Второй вычислительный этап: электронная теплопроводность
// Имеет смысл только для EOSAnalytic и EOSTable

/*
1.На границах A,B действительно должны быть равны 0? Может можно усреднить ro,
  kappa с 0 или продлить их крайние значения, и считать границы как обычные точки?
2.Правда ли, что массивы коэффициентов A,B,C,F для каждой итерации одни и те же?
  Сейчас ce и kappa, по которым они считаются, изменяются на каждой итерации.
  Если правда, то ce и kappa можно считать один раз в конце. Если нет - нужно
  пересчитывать коэффициенты.
*/

// Решение системы линейных уравнений с трехдиагональной матрицей методом прогонки.

void CSolver::sweepTe(CField& ms_temp, CField& ms_temp_temp, double *A, double *B, double *C, double *F,
					  double *alpha, double *beta, int size) {
	int i = 0;
//	int iMin = 0;

	double teL=0.0, teR=0.0, alphaR=0.0, betaR=0.0;

	// Граничные условия

	double  AL = 0.0,
			BL = 1.0,
			CL = 1.0,
			FL = 0.0,

			AR = 1.0,
			BR = 0.0,
			CR = 1.0,
			FR = 0.0;

	// Прямая прогонка

		alpha[i_pulse_min] = BL/CL;
		 beta[i_pulse_min] = FL/CL;
		
	for(i=i_pulse_min; i<size-1; i++)
	{
		alpha[i+1] = B[i] / (C[i]-A[i]*alpha[i]);
		 beta[i+1] = (A[i]*beta[i]+F[i]) / (C[i]-A[i]*alpha[i]);
	}

	alphaR = B[size-1] / (C[size-1]-A[size-1]*alpha[size-1]);
	 betaR = (A[size-1]*beta[size-1]+F[size-1]) / 
		     (C[size-1]-A[size-1]*alpha[size-1]);
	
	// Обратная прогонка

	teR = (AR/CR*betaR + FR/CR) / (1.0 - AR/CR*alphaR);
	
	i = size-1;

	ms_temp_temp[i].te = ms_temp[i].te;
	     ms_temp[i].te = alphaR * teR + betaR;

	for(i=size-2; i>=(int)i_pulse_min; i--) 








		// Почему-то тут i может становиться меньше 0! Ошибка с типами signed/unsigned int???


	





	{


		if(i==1)
		{

			double qq=0.;
		}




		ms_temp_temp[i].te = ms_temp[i].te;
		     ms_temp[i].te = alpha[i+1] * ms_temp[i+1].te + beta[i+1];
	}

	teL = alpha[i_pulse_min] * ms_temp[i_pulse_min].te + beta[i_pulse_min];
}


void CSolver::sweepTi(CField& ms_temp, CField& ms_temp_temp, 
					  double *A, double *B, double *C, double *F,
					  double *alpha, double *beta, int size) {
	int i = 0;
	double tiL=0.0, tiR=0.0, alphaR=0.0, betaR=0.0;
	// Граничные условия
	double  AL = 0.0, BL = 1.0, CL = 1.0, FL = 0.0,
			AR = 1.0, BR = 0.0, CR = 1.0, FR = 0.0;
	// Прямая прогонка
		alpha[0] = BL/CL;
		 beta[0] = FL/CL;
		
	for(i=0; i<i_pulse_min; i++) {
		alpha[i+1] = B[i] / (C[i]-A[i]*alpha[i]);
		 beta[i+1] = (A[i]*beta[i]+F[i]) / (C[i]-A[i]*alpha[i]);
	}
	alphaR = B[size-1] / (C[size-1]-A[size-1]*alpha[size-1]);
	 betaR = (A[size-1]*beta[size-1]+F[size-1]) / 
		     (C[size-1]-A[size-1]*alpha[size-1]);	
	// Обратная прогонка
	tiR = (AR/CR*betaR + FR/CR) / (1.0 - AR/CR*alphaR);	
	i = i_pulse_min-1;
	ms_temp_temp[i].ti = ms_temp[i].ti;
	     ms_temp[i].ti = alphaR * tiR + betaR;
	for(i=i_pulse_min-2; i>=0; i--) {
		ms_temp_temp[i].ti = ms_temp[i].ti;
		     ms_temp[i].ti = alpha[i+1] * ms_temp[i+1].ti + beta[i+1];
	}
	tiL = alpha[0] * ms_temp[0].te + beta[0];
}

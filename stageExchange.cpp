#include <iostream>
#include "solver.h"
using namespace std;

// Третий вычислительный этап: электронно-ионный обмен
int CSolver::calcExchangeStage(double tau) {
	int i=0, itNum = 0, iMin=0;
	double eps = .01;	
	EOSOld &eos = task.getEOS();
	CFieldOld ms_temp, ms_temp_temp;
	ms_temp.initData(&task);
	ms_temp_temp.initData(&task);
	for(i=0; i<ms.getSize(); i++) {
		ms_temp[i].te = ms[i].te;
		ms_temp[i].ti = ms[i].ti;
		ms_temp[i].Z  = ms[i].Z;
	}
	for(;;) {
		for(i=iMin; i<ms.getSize(); i++) {
			Node &n = ms[i];
			ms_temp_temp[i].te=ms_temp[i].te;
			ms_temp_temp[i].ti=ms_temp[i].ti;	
			ms_temp[i].ti = (tau * n.Alphaei*n.te + (tau * n.Alphaei/n.ce - 1.0)*n.ci*n.ti) / 
				            (tau * n.Alphaei*(1.0 + n.ci/n.ce) - n.ci);
			ms_temp[i].te = n.te + n.ci/n.ce*n.ti - n.ci/n.ce * ms_temp[i].ti;
			/*
			te_temp[i] =te[i] + ci[i]/ce[i]*ti[i] - ci[i]/ce[i] *
					 (tau*Alphaei[i]*te[i] + (tau*Alphaei[i]/ce[i] - 1.0)*ci[i]*ti[i]) /
					 (tau*Alphaei[i]*(1.0 + ci[i]/ce[i]) - ci[i]);
	        */
			n.ce = eos.getce(n.ro, ms_temp[i].te);
			n.ci = eos.getci(n.ro, ms_temp[i].ti);
			n.Alphaei  = eos.getAlpha(n.ro, ms_temp[i].ti, ms_temp[i].te);
		}
		if(compTe(ms_temp, ms_temp_temp)<eps && compTi(ms_temp, ms_temp_temp)<eps) break;		
		itNum++;
		if(itNum > maxIt) {
			printf("Divergence=%f, iteration number=%d\n", compTe(), itNum);
			printf("Error. Too many iterations for convergence.\n")	;
			cin.get();
			exit(1);
		}
	}
	//cout << "calcExchangeStage() complete: " << itNum+1 << " iterations of Ti, Te" << endl;	
	for(i=0; i<ms.getSize(); i++) {
		Node &n = ms[i];
		n.te = ms_temp[i].te;
		n.ti = ms_temp[i].ti;
		n.pe = eos.getpe(n.ro, n.ti, n.te);
		n.pi = eos.getpi(n.ro, n.ti);
		n.p  = n.pe + n.pi;
		n.ee = eos.getee(n.ro, n.ti, n.te);
		n.ei = eos.getei(n.ro, n.ti);
		n.e  = n.ee + n.ei; 
		n.ci    = eos.getci(n.ro, n.ti);
		n.C     = eos.getC(n.ro, n.ti, n.te);				
		n.ce    = eos.getce(n.ro, n.te);
		n.kappa = eos.getkappa(n.ro, n.ti, n.te);		
	}
	return itNum+1;
}

int CSolver::calcExchangeStageGlass(double tau) {
	int i=0, itNum = 0;
	const int nBound = task.getZone(0).n;
	const double eps = 0.01;
	EOSOld &eos = task.getEOS();
	EOSOld &eosGlass = task.getEOSGlass();
	CFieldOld ms_temp, ms_temp_temp;
	ms_temp.initData(&task);
	ms_temp_temp.initData(&task);
	for(i=0; i<ms.getSize(); i++) {
		ms_temp[i].te = ms[i].te;
		ms_temp[i].ti = ms[i].ti;
	}
	int iInit = 0;
	if(task.getSourceFlag() == 2) 
		iInit = nBound;
	//cout << "Exchg:";
	for(;;) {
		for(i=0; i<ms.getSize(); i++) {
			Node &n = ms[i];
			ms_temp_temp[i].te=ms_temp[i].te;
			ms_temp_temp[i].ti=ms_temp[i].ti;
			ms_temp[i].ti = (tau * n.Alphaei*n.te + (tau * n.Alphaei/n.ce - 1.0)*n.ci*n.ti) / 
				            (tau * n.Alphaei*(1.0 + n.ci/n.ce) - n.ci);
			ms_temp[i].te = n.te + n.ci/n.ce*n.ti - n.ci/n.ce * ms_temp[i].ti;
			if(task.getSourceFlag() == 2) {
				if(i>=nBound) {
					n.ce = eos.getce(n.ro, ms_temp[i].te);
					n.ci = eos.getci(n.ro, ms_temp[i].ti);
					n.Alphaei  = eos.getAlpha(n.ro, ms_temp[i].ti, ms_temp[i].te);
				} else {
					n.ce = eosGlass.getce(n.ro, ms_temp[i].te);
					n.ci = eosGlass.getci(n.ro, ms_temp[i].ti);
					n.Alphaei  = eosGlass.getAlpha(n.ro, ms_temp[i].ti, ms_temp[i].te);
				}
			} else if(i<nBound) {
				n.ce = eos.getce(n.ro, ms_temp[i].te);
				n.ci = eos.getci(n.ro, ms_temp[i].ti);
				n.Alphaei  = eos.getAlpha(n.ro, ms_temp[i].ti, ms_temp[i].te);
			} else {
				n.ce = eosGlass.getce(n.ro, ms_temp[i].te);
				n.ci = eosGlass.getci(n.ro, ms_temp[i].ti);
				n.Alphaei  = eosGlass.getAlpha(n.ro, ms_temp[i].ti, ms_temp[i].te);
			}
		}

		if(compTe(ms_temp, ms_temp_temp)<eps && 
		   compTi(ms_temp, ms_temp_temp)<eps) break;
		itNum++;
		if(itNum >= maxIt) {
			printf("Divergence=%f, iteration number=%d\n", compTe(), itNum);
			printf("Error. Too many iterations for convergence.\n")	;
			cin.get();
			exit(1);
		}
	}
	//cout << itNum+1 << "it";
	for(i=0; i<ms.getSize(); i++) {
		Node &n = ms[i];
		n.te = ms_temp[i].te;
		n.ti = ms_temp[i].ti;
		if(task.getSourceFlag() == 2) {
			if(i>=nBound) {
				n.pe = eos.getpe(n.ro, n.ti, n.te);
				n.pi = eos.getpi(n.ro, n.ti);
				n.p  = n.pe + n.pi;
				n.ee = eos.getee(n.ro, n.ti, n.te);
				n.ei = eos.getei(n.ro, n.ti);
				n.e  = n.ee + n.ei;
				n.ci    = eos.getci(n.ro, n.ti);
				n.C     = eos.getC(n.ro, n.ti, n.te);
				n.ce    = eos.getce(n.ro, n.te);
				n.kappa = eos.getkappa(n.ro, n.ti, n.te);
			} else {
				n.pe = eosGlass.getpe(n.ro, n.ti, n.te);
				n.pi = eosGlass.getpi(n.ro, n.ti);
				n.p  = n.pe + n.pi;
				n.ee = eosGlass.getee(n.ro, n.ti, n.te);
				n.ei = eosGlass.getei(n.ro, n.ti);
				n.e  = n.ee + n.ei;
				n.ci    = eosGlass.getci(n.ro, n.ti);
				n.C     = eos.getC(n.ro, n.ti, n.te);
				n.ce    = eosGlass.getce(n.ro, n.te);	
				n.kappa = eosGlass.getkappa(n.ro, n.ti, n.te);
			}
		} else if(i<nBound) {
			n.pe = eos.getpe(n.ro, n.ti, n.te);
			n.pi = eos.getpi(n.ro, n.ti);
			n.p  = n.pe + n.pi;
			n.ee = eos.getee(n.ro, n.ti, n.te);
			n.ei = eos.getei(n.ro, n.ti);
			n.e  = n.ee + n.ei;
			n.ci    = eos.getci(n.ro, n.ti);
			n.C     = eos.getC(n.ro, n.ti, n.te);
			n.ce    = eos.getce(n.ro, n.te);
			n.kappa = eos.getkappa(n.ro, n.ti, n.te);
		} else {
			n.pe = eosGlass.getpe(n.ro, n.ti, n.te);
			n.pi = eosGlass.getpi(n.ro, n.ti);
			n.p  = n.pe + n.pi;
			n.ee = eosGlass.getee(n.ro, n.ti, n.te);
			n.ei = eosGlass.getei(n.ro, n.ti);
			n.e  = n.ee + n.ei;
			n.ci    = eosGlass.getci(n.ro, n.ti);
			n.C     = eos.getC(n.ro, n.ti, n.te);
			n.ce    = eosGlass.getce(n.ro, n.te);
			n.kappa = eosGlass.getkappa(n.ro, n.ti, n.te);
		}		
	}
	return itNum+1;
}

void CSolver::calcExchangeStage5LayersSi(double tau) {
/*	int i = 0, itNum = 0, size = ms.getSize();
	double eps = 0.01;
	double x = 0., dx = 0.;
	EOS &eosCr = task.getEOSCr(), &eosAu = task.getEOSAu(), &eosSi = task.getEOSSi();
	MatterState ms_temp, ms_temp_temp;
	ms_temp.initData(&task); ms_temp_temp.initData(&task);
	for(;;)
	{
		for(i=0; i<size; i++) {
			Node &n = ms[i]; x = .5*(ms[i].x + ms[i+1].x);
			ms_temp_temp[i].te=ms_temp[i].te;
			ms_temp_temp[i].ti=ms_temp[i].ti;
			ms_temp[i].ti = (tau * n.Alphaei*n.te + (tau * n.Alphaei/n.ce - 1.0)*n.ci*n.ti) / 
				            (tau * n.Alphaei*(1.0 + n.ci/n.ce) - n.ci);
			ms_temp[i].te = n.te + n.ci/n.ce*n.ti - n.ci/n.ce * ms_temp[i].ti;
			if ( x > 250.e-9) {
				n.ce = eosSi.getce(n.ro, ms_temp[i].te, n.Z);
				n.ci = eosSi.getci(n.ro, ms_temp[i].ti);
				n.Alphaei  = eosSi.getAlpha(n.ro, ms_temp[i].ti, ms_temp[i].te);
			} else if (x > 200.e-9) {
				n.ce = eosCr.getce(n.ro, ms_temp[i].te, n.Z);
				n.ci = eosCr.getci(n.ro, ms_temp[i].ti);
				n.Alphaei  = eosCr.getAlpha(n.ro, ms_temp[i].ti, ms_temp[i].te);
			} else if (x > 150.e-9) {
				n.ce = eosAu.getce(n.ro, ms_temp[i].te, n.Z);
				n.ci = eosAu.getci(n.ro, ms_temp[i].ti);
				n.Alphaei  = eosAu.getAlpha(n.ro, ms_temp[i].ti, ms_temp[i].te);
			} else if (x > 100.e-9) {
				n.ce = eosCr.getce(n.ro, ms_temp[i].te, n.Z);
				n.ci = eosCr.getci(n.ro, ms_temp[i].ti);
				n.Alphaei  = eosCr.getAlpha(n.ro, ms_temp[i].ti, ms_temp[i].te);
			} else if (x > 50.e-9) {
				n.ce = eosAu.getce(n.ro, ms_temp[i].te, n.Z);
				n.ci = eosAu.getci(n.ro, ms_temp[i].ti);
				n.Alphaei  = eosAu.getAlpha(n.ro, ms_temp[i].ti, ms_temp[i].te);
			} else {
				n.ce = eosCr.getce(n.ro, ms_temp[i].te, n.Z);
				n.ci = eosCr.getci(n.ro, ms_temp[i].ti);
				n.Alphaei  = eosCr.getAlpha(n.ro, ms_temp[i].ti, ms_temp[i].te);
			}
		}

		if(compTe(ms_temp, ms_temp_temp)<eps && 
		   compTi(ms_temp, ms_temp_temp)<eps) break;
		
		itNum++;

		if(itNum >= maxIt)
		{
			printf("Divergence=%f, iteration number=%d\n", compTe(), itNum);
			printf("Error. Too many iterations for convergence.\n")	;
			cin.get();
			exit(1);
		}
	}
	cout << "calcExchangeStage() complete: " << itNum+1 << " iterations of Ti, Te" << endl;
	for(i=0; i<size; i++)	{
		Node &n = ms[i]; x = .5*(ms[i].x + ms[i+1].x);
		n.te    = ms_temp[i].te;
		n.ti = ms_temp[i].ti;		
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
	}*/
}

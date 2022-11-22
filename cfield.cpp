#include "cfieldold.h"
#include "task.h"
#include "eosold.h"
#include "eosTable.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fstream>

using namespace std;

CFieldOld::CFieldOld() {
	nodes = 0;
	nSize = 0;
}


CFieldOld::~CFieldOld() {
	clearData();
}


void CFieldOld::initData(CTask *task) {
	clearData();	
	// Установка начальных условий
	double dx = 0.;
	double nextX   = 0.0;
	int i=0, counter=0;	
	EOSOld &eos = task->getEOS();
	EOSOld &eosGlass = task->getEOSGlass();
	EOSBin &eosBin = task->eosBin;
	int iMin=2, iMax=2+task->NX;
	double _e=0;
	if(task->type == TaskType::RP1D) {
		nSize = task->NX+2+2;
		nodes = new Node[nSize];
		dx = (task->xMax-task->xMin)/task->NX;
		for(i=iMin; i<iMax; i++) {
			Node &n = nodes[i];
			n.x = task->xMin + (.5 + i)*dx;
			if(n.x < task->qBound) {
				n.W[0] = task->roL;
				n.W[1] = task->roL*task->uL;
				_e = eosBin.gete(task->roL, task->pL); 
				n.W[2] = task->roL*(_e + .5*task->uL*task->uL);				
			} else {
				n.W[0] = task->roR;
				n.W[1] = task->roR*task->uR;
				_e = eosBin.gete(task->roR, task->pR); 
				n.W[2] = task->roR*(_e + .5*task->uR*task->uR);
			}
			n.W[3] = 0.;	
			n.W_temp = Vector4::ZERO;
			n.F = Vector4::ZERO;
		}		
		nodes[iMin-1].x = task->xMin + (.5-1)*dx; nodes[iMin-1].W = Vector4(nodes[iMin].W);   nodes[iMin-1].W_temp = Vector4::ZERO; nodes[iMin-1].F = Vector4::ZERO; 
		nodes[iMin-2].x = task->xMin + (.5-2)*dx; nodes[iMin-2].W = Vector4(nodes[iMin].W);   nodes[iMin-2].W_temp = Vector4::ZERO; nodes[iMin-2].F = Vector4::ZERO; 
		nodes[iMax].x = task->xMax + .5*dx;       nodes[iMax].W = Vector4(nodes[iMax-1].W);   nodes[iMax].W_temp = Vector4::ZERO;   nodes[iMax].F = Vector4::ZERO; 
		nodes[iMax+1].x = task->xMax + (.5+1)*dx; nodes[iMax+1].W = Vector4(nodes[iMax-1].W); nodes[iMax+1].W_temp = Vector4::ZERO; nodes[iMax+1].F = Vector4::ZERO; 
	} else {
		nSize = task->getTotalSize();
		nodes = new Node[nSize+1];	
		for(i=0; i<task->getNumZones(); i++) {
			Zone &zone = task->getZone(i);
			dx = zone.l / (zone.n);
			double xInit=nextX;
			for(int j=0; j<zone.n; j++)	{
				Node &n = nodes[counter];			
				n.x = nextX;
				n.v = zone.v;
				n.ro = zone.ro;
				n.dm = n.ro * dx;
				n.ti = zone.ti;
				n.te = zone.te;





				////////////// Patch for Toro test via BGK scheme
				n.C = eos.getC(n.ro, n.ti, n.te);
				///////////////////////////////////////////////////





				if( (task->getSourceFlag() == SourceType::SrcGlass) && (i==0) ) {
						n.pe = eosGlass.getpe(n.ro, n.ti, n.te);
						n.pi = eosGlass.getpi(n.ro, n.ti);
						n.p  = eosGlass.getp (n.ro, n.ti, n.te);
						n.ee = eosGlass.getee(n.ro, n.ti, n.te);
						n.ei = eosGlass.getei(n.ro, n.ti);
						n.e  = eosGlass.gete (n.ro, n.ti, n.te);
						n.ce = eosGlass.getce(n.ro, n.te);
						n.ci = eosGlass.getci(n.ro, n.ti);
						n.C  = eos.getC(n.ro, n.ti, n.te);
						n.Alphaei     = eosGlass.getAlpha(n.ro, n.ti, n.te);
						n.kappa = eosGlass.getkappa(n.ro, n.ti, n.te);
					} 
				if(task->getSourceFlag()==4) {} 
				if(task->type == TaskType::ruGlass) {  // Пока не навел порядок в инфраструктуре, проставляю этот флаг руками в CSolver::goGlass()			        
						if(i==0) {
							n.pe = eos.getpe(n.ro, n.ti, n.te);
							n.pi = eos.getpi(n.ro, n.ti);
							n.p  = eos.getp (n.ro, n.ti, n.te);
							n.ee = eos.getee(n.ro, n.ti, n.te);
							n.ei = eos.getei(n.ro, n.ti);
							n.e  = eos.gete (n.ro, n.ti, n.te);
							n.ce = eos.getce(n.ro, n.te);
							n.ci = eos.getci(n.ro, n.ti);
							n.C  = eos.getC(n.ro, n.ti, n.te);
							n.Alphaei     = eos.getAlpha(n.ro, n.ti, n.te);
							n.kappa = eos.getkappa(n.ro, n.ti, n.te);
						} else {
							n.pe = eosGlass.getpe(n.ro, n.ti, n.te);
							n.pi = eosGlass.getpi(n.ro, n.ti);
							n.p  = eosGlass.getp (n.ro, n.ti, n.te);
							n.ee = eosGlass.getee(n.ro, n.ti, n.te);
							n.ei = eosGlass.getei(n.ro, n.ti);
							n.e  = eosGlass.gete (n.ro, n.ti, n.te);
							n.ce = eosGlass.getce(n.ro, n.te);
							n.ci = eosGlass.getci(n.ro, n.ti);
							n.C  = eos.getC(n.ro, n.ti, n.te); // сознательно eos, а не EOSGlass, xчтобы не завысить шаг по времени
							n.Alphaei     = eosGlass.getAlpha(n.ro, n.ti, n.te);
							n.kappa = eosGlass.getkappa(n.ro, n.ti, n.te);
						}	
				}
				if(task->type == TaskType::MieGruneisenProblem) {
					n.ee = 0.;
					n.ei = 0.;
					n.e  = zone.e;
					n.pe = 0.;
					n.pi = 0.;
					n.p  = 0.;
					n.ce = 0.;
					n.ci = 0.;
					n.C  = 0.;
					n.Alphaei = 0.;
					n.kappa = 0.;
				} 				
				// flow
				double E = n.e + 0.5*n.v*n.v;
				n.W = Vector4(n.ro, n.ro*n.v, n.ro*E, 0);
				n.W_temp = Vector4::ZERO;
				n.F = Vector4::ZERO;
				counter++;
				nextX += dx;
			}
		}
	
		nodes[nSize].x = nodes[nSize-1].x + dx;
		nodes[nSize].v = nodes[nSize-1].v;
	}
}


void CFieldOld::clearData()
{
	delete[] nodes;
	nodes = 0;
	nSize = 0;
}

double CFieldOld::loadData(string fName, int nCut) {
	string buf = string("");
	string fullName = string(OUTPUT_FOLDER) + fName;
	ifstream fInput;
	fInput.open(fullName, ios::in);
	if(!fInput.is_open()){
		cout << "Error: cannot open data file!" << endl;
		exit(1);
	}
	int i=0, j=0;
	double _t=0.;
	for(i=0; i<nSize; i++) {
		if(i<nCut) {
			for(j=0; j<17; j++)
				fInput >> buf;
		} else {
			fInput >> nodes[i-nCut].x; 
			fInput >> nodes[i-nCut].v;
			fInput >> nodes[i-nCut].ro;
			fInput >> nodes[i-nCut].dm;
			fInput >> nodes[i-nCut].ti;
			fInput >> nodes[i-nCut].te;
			fInput >> nodes[i-nCut].pi;
			fInput >> nodes[i-nCut].pe;
			fInput >> nodes[i-nCut].p;
			fInput >> nodes[i-nCut].ei;
			fInput >> nodes[i-nCut].ee;
			fInput >> nodes[i-nCut].e;
			fInput >> nodes[i-nCut].ci;
			fInput >> nodes[i-nCut].ce;
			fInput >> nodes[i-nCut].C;
			fInput >> nodes[i-nCut].Alphaei;
			fInput >> nodes[i-nCut].kappa;
		}
	}
	// Меняем nSize
	setSize(nSize-nCut);
	i = nSize;
	fInput >> nodes[i].x; 
	fInput >> nodes[i].v;
	fInput >> nodes[i].dm;
	fInput >> _t;
	//Тесты
	cout << "Solution successfully read!" << endl;
	fInput.close();
	/*
	for(i=0; i<ms.getSize(); i++){
		Node &n = ms[i];
		fprintf(f, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
				    n.x,  n.v, n.ro, n.dm, 
					n.ti, n.te, n.pi, n.pe, n.p, n.ei, n.ee, n.e, 
					n.ci, n.ce, n.C,  n.Alphaei, n.kappa);
	}
	fLog << "Variables order: x,v,ro,dm,ti,te,pi,pe,p,ei,ee,e,ci,ce,C,Alphaei,kappa (17 variables)" << endl;
	double S1 = ms[0].x  + ms[0].v  + ms[0].ro + ms[0].dm + ms[0].ti + ms[0].te + ms[0].pi + ms[0].pe +
		        ms[0].p  + ms[0].ei + ms[0].ee + ms[0].e  + ms[0].ci + ms[0].ce + ms[0].C  + ms[0].Alphaei +
				ms[0].kappa;
	fLog << "First checksum: S1 = ms[0].x + ms[0].v + ... + ms[0].Alphaei + ms[0].kappa" << endl << "S1 = " << S1 << endl;
	i=ms.getSize();
	// Писать все сплошняком, как в цикле выше. В строчке с номером энСайз+1 давать во всех лишних переменных нули.
	//double x, v, ro, dm;
	//double te, ti, pe, pi, p, ee, ei, e;
	//double ce, ci, C, Alphaei, kappa, ne, Z, ti_temp, te_temp;
	Node &n = ms[i];
	fprintf(f, "%e %e %e\n", n.x, n.v, n.dm);
	fprintf(f, "%e\n", t);
	double S2 = ms[i].x + ms[i].v + ms[i].dm;
	fLog << "Second checksum: S1 = ms[nSize].x + ms[nSize].v + ms[nSize].dm" << endl << "S2 = " << S2 << endl;
	double S3 = 0;
	for(i=0; i<=ms.getSize();i++) {
		S3 += ms[i].x;
	}
	fLog << "Third checksum: S3 = ms[0].x + ms[1].x + ... + ms[nSize].x" << endl << "S3 = " << S3 << endl;
	fLog.close();
	fclose(f);
	*/


	return _t;
}



void CFieldOld::setEdge(Node &n, double x, double dm)
{
	n.x  = x;
	n.v  = 0;

	n.ro = 0;
	n.dm = dm;

	n.ti = 0;
	n.te = 0;

	n.pe = 0;
	n.pi = 0;
	n.p  = 0;

	n.ee = 0;
	n.ei = 0;
	n.e  = 0;

	n.ci = 0;
	n.ce = 0;

	n.C  = 0;
	n.Alphaei     = 0;
	n.kappa = 0;

	n.ne = 0.0;

	// flow

	n.F = Vector4::ZERO; //Vector4(1.0e-5);  //Vector4(0.00010188989189435); !!!!!!!!!
	n.W = Vector4::ZERO;

	// temporary

	n.W_temp = Vector4::ZERO;
	n.ti_temp = 0;
	n.te_temp = 0;
}


void CFieldOld::setEdgeTransparent()
{
	left_edge.x   = nodes[0].x - (nodes[1].x - nodes[0].x);
/*	left_edge.ro  = 0.; 
	left_edge.v   = 0.;//nodes[0].v;
	left_edge.pe  = 0.;
	left_edge.pi  = 0.;
	left_edge.p   = 0.;
	left_edge.ee  = 0.;
	left_edge.ei  = 0.;
	left_edge.e   = 0.;
	left_edge.C   = 0.;
	left_edge.W[0]   = 0.;
	left_edge.W[1]   = 0.;
	left_edge.W[2]   = 0.;
	left_edge.W[3]   = 0.;
	left_edge.F[0]   = 0.;
	left_edge.F[1]   = 0.;
	left_edge.F[2]   = 0.;
	left_edge.F[3]   = 0.;*/
	////////// DEBUG 17.09.2014 /////////
	left_edge.ro  = nodes[0].ro; 
	left_edge.v   = nodes[0].v;
	left_edge.pe  = nodes[0].pe;
	left_edge.pi  = nodes[0].pi;
	left_edge.p   = nodes[0].p;
	left_edge.ee  = nodes[0].ee;
	left_edge.ei  = nodes[0].ei;
	left_edge.e   = nodes[0].e;
	left_edge.C   = nodes[0].C;
	left_edge.W[0]   = nodes[0].W[0];
	left_edge.W[1]   = nodes[0].W[1];
	left_edge.W[2]   = nodes[0].W[2];
	left_edge.W[3]   = nodes[0].W[3];
	left_edge.F[0]   = nodes[0].F[0];
	left_edge.F[1]   = nodes[0].F[1];
	left_edge.F[2]   = nodes[0].F[2];
	left_edge.F[3]   = nodes[0].F[3];


	right_edge.x  = nodes[nSize-1].x + (nodes[nSize-1].x - nodes[nSize-2].x);
	right_edge.ro = nodes[nSize-1].ro; 
	right_edge.v  = nodes[nSize-1].v;
	right_edge.pe = nodes[nSize-1].pe;
	right_edge.pi = nodes[nSize-1].pi;
	right_edge.p  = nodes[nSize-1].p;
	right_edge.ee = nodes[nSize-1].ee;
	right_edge.ei = nodes[nSize-1].ei;
	right_edge.e  = nodes[nSize-1].e;
	right_edge.C  = nodes[nSize-1].C;
	right_edge.W[0]   = nodes[nSize-1].W[0];
	right_edge.W[1]   = nodes[nSize-1].W[1];
	right_edge.W[2]   = nodes[nSize-1].W[2];
	right_edge.W[3]   = nodes[nSize-1].W[3];
	right_edge.F[0]   = nodes[nSize-1].F[0];
	right_edge.F[1]   = nodes[nSize-1].F[1];
	right_edge.F[2]   = nodes[nSize-1].F[2];
	right_edge.F[3]   = nodes[nSize-1].F[3];

}


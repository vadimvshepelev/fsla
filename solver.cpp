#define _CRT_SECURE_NO_WARNINGS

#include "defines.h"
#include "solver.h"
#include "cfield.h"
#include "eosFigures.h"
#include "task.h"
#include "node.h"
#include "methodEuler.h"
#include "ceos.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <complex>
#include <cstring>
#include <fstream>
#include <vector>
#include <list>

#ifdef _WIN32
#include <conio.h>
#endif

using namespace std;

CSolver::CSolver() {
	maxIt = 40;
	epsE  = 0.1;
	CFL  =  0.0;
	tauPulse = 0.0;
	fluence  = 0.0;
	deltaSkin = 0.0;
	x_pulse_min = 0.0;
	i_pulse_min = 0;
	tInit = 0.0;
	iWeak = -1;
	spallFlag = 0;
	spallCellNum = -1;
	xVacBound = 0.;
}

void CSolver::goEuler(char* fName) {
	// Загрузка входного файла, инициализация объекта task типа CTask
	task.load(fName);
	CFL = task.getCFL();
	// Наполнение объектов ms, ms_temp данными на по начальным и граничным условиям из task
	ms.initData(&task);
	ms_temp.initData(&task);
	const int nSize = ms.getSize();	
	int i=0;
/*  // ENO-3 testing
	string _fName = string(OUTPUT_FOLDER) + "\\" + "test-reconstruction.dat";
	ofstream ofs; ofs.open(_fName, ios::out);
	ofs << "TITLE = \"ENO-3 reconstruction test, f(x) = x^2\"" << endl;
	ofs << "VARIABLES = \"x\", \"f(x)\", \"sq(x-100.)\"" << endl;
	for(i=0; i<nSize; i++) {
		double _x = 0.5*(ms[i].x + ms[i+1].x);
			ofs << _x << " " << ms[i].ro << " " << (_x)*(_x) << endl;
	}

	ofs.close();
	*/
	double t   = 0.;
	double tau = 0.0;
	int	counter = 0;
	//dumpToFileTestRP(t+100.e-6, counter);
    //	dumpToFileTestRP(t, counter);
	int nTimes = 20, nTimesCounter = 0;
	double *timesArray = new double[nTimes];
	const double tMax = task.getMaxTime();
	for(i=0; i<nTimes; i++) {
		timesArray[i] = (tMax-tInit)/(nTimes-1)*i;
	}



	// For Mie-Gruneisen Riemann solver tests
	CEOSMieGruneisen eos;





	for(;;)	{
		tau = calcTimeStepEuler(t);
		if(t+tau >tMax) tau = tMax-t;
		if(counter <=4) tau *=.2;
		cout << counter << ": " << "t=" << t+tau <<  " tau=" << tau << " CFL=" << getCFL() << endl;
		//if(task.getHydroStage()) calcHydroStageGodunov(t, tau);
		//if(task.getHydroStage()) calcHydroStageRoe(t, tau);		
		//if(task.getHydroStage()) calcHydroStageGPS(t, tau);	
		//if(task.getHydroStage()) calcHydroStageLaxFriedrichs(t, tau);	
		//if(task.getHydroStage()) calcHydroStageMccormack(t, tau);		
		//if(task.getHydroStage()) calcHydroStageMHM(t, tau);		
		//if(task.getHydroStage()) calcHydroStageGushchinIdealSimple(t, tau);
		//if(task.getHydroStage()) calcHydroStageG2(t, tau);	
		//if(task.getHydroStage()) calcHydroStageENO2G(t, tau);	
		//if(task.getHydroStage()) calcHydroStageENO3G(t, tau);	
		if(task.getHydroStage()) calcHydroStageMieGruneisen(eos, t, tau);	
		//if(task.getHydroStage()) calcHydroStageGodunovEOSBin(t, tau);		
		//if(task.getHydroStage()) calcHydroStageENO2G(t, tau);		
		if(handleKeys(t)) break;		
		// Regular file output
		t += tau;
		dumpToFileTestRP(t+tau, 100);		
		if (t>=timesArray[nTimesCounter]) {
			dumpToFileTestRP(t, nTimesCounter++);
			if(nTimesCounter==nTimes)
				break;
		}
		if( t>= tMax) {
			break;
		}		
		counter++;
	} 
	// Удаляем динамический массив времен выдачи
	delete[] timesArray;
}

void CSolver::goEulerMovingMesh(char* fName) {
	// Загрузка входного файла, инициализация объекта task типа CTask
	task.load(fName);
	// Наполнение объектов ms, ms_temp данными на по начальным и граничным условиям из task
	ms.initData(&task);
	ms_temp.initData(&task);
	initVars();
	double t   = 0.;
	double tau = 0.0;
	int	counter = 0;
	int i=0;
	int nSize = ms.getSize();
	dumpToFileTestRP(t+100.e-6, counter);
	int nTimes = 20, nTimesCounter = 0;
	CMethod& method = task.getMethod();
	method.vGrid = new double[nSize];
	method.X	 = new double[nSize];
	double *timesArray = new double[nTimes];
	double tMax = task.getMaxTime();
	for(i=0; i<nTimes; i++) {
		timesArray[i] = (tMax-tInit)/(nTimes-1)*i;
	}
	// Протестируем последовательность, на строгое возрастание (выход за границу мы не контролируем)
	for(i = 0; i<nTimes-1; i++)
		if(timesArray[i+1] <= timesArray[i]) {
			cout << "CSolver::go() error: Time points unordered in timesArray." << endl;
			exit(1);
		}
	// Начинаем счет
	for(;;)	{
		tau = calcTimeStepEuler(t);
		
		//if(counter <=4) tau *=.2;
		cout << counter << ": " << "t=" << t <<  " tau=" << tau << " courant=" << getCFL() << endl;
		//if(task.getHydroStage()) calcHydroStageGodunov(t, tau);
		
		if(task.getHydroStage()) calcHydroStageGodunovMovingMesh(t, tau)/*calcHydroStageGodunov(t, tau)*/;
		dumpToFileEuler(t+tau);

		//dumpToFileTestRP(t, 100);

		if(handleKeys(t)) break;
		// Regular file output
		if (t>=timesArray[nTimesCounter]) {
			dumpToFileTestRP(t+tau, nTimesCounter++);
			if(nTimesCounter==nTimes)
				break;
		}
		if( t>= task.getMaxTime()) {
			break;
		}
		t += tau;
		counter++;
	} 

	cout << endl << "Calculation finished!" << endl;
	// Удаляем сетку
	delete[] method.vGrid;
	delete[] method.X;
	// Удаляем динамический массив времен выдачи
	delete[] timesArray;

}



void CSolver::go(char* fName) {
	// Загрузка входного файла, инициализация объекта task типа CTask
	task.load(fName);
	// Наполнение объектов ms, ms_temp данными на по начальным и граничным условиям из task
	ms.initData(&task);
	ms_temp.initData(&task);
	initVars();	
	double t = 0., tau = 0.0;
	int	counter = 0, i=0;
	/*double timesArray[] = {      .0,   1.e-13,   3.e-13,   5.e-13,   1.e-12,   2.e-12,   3.e-12,    4.e-12,    5.e-12,   6.e-12,
		                    10.e-12,  20.e-12,  30.e-12,  40.e-12,  50.e-12,  60.e-12,  70.e-12,  80.1e-12,  90.2e-12, 100.e-12};
	int nTimes = 20;*/
	double timesArray[] = {0.,        2.e-12,  5.e-12,   10.e-12,   15.e-12, 16.e-12, 17.e-12,  18.e-12,  19.e-12,  20.e-12, 
		                   21.e-12,  22.e-12, 23.e-12, 23.87e-12, 23.88e-12, 25.e-12,  30.e-12, 50.e-12, 100.e-12, 200.e-12, 
						   500.e-12,   5.e-9,  50.e-9,   500.e-9};
	int nTimes = 24, timesCounter = 0;	
	dumpToFile(t);
	assert(task.getMethodFlag() == MethodType::samarskii);
	for(;;) {
		tau = calcTimeStep(t);
		if (t + tau > task.getMaxTime()) tau = task.getMaxTime() - t;
		cout << counter << ": " << "t=" << t <<  " tau=" << tau << " courant=" << getCFL() << endl;
		if(task.getHydroStage()) calcHydroStage(t, tau);
		// if(task.getHeatStage()) calcHeatStage5LayersSi(t, tau);
		if(task.getHeatStage()) calcHeatStage(t, tau);				
		//if(task.getExchangeStage()) calcExchangeStage5LayersSi(tau);
		if(task.getExchangeStage()) calcExchangeStage(tau);				
		// Regular file output
		if (t>=timesArray[timesCounter]) {
			if(task.getMethodFlag() == MethodType::samarskii) {
				dumpToFile(t); 
			} else {
				dumpToFileEuler(t);
				dumpToFileTestRP(t, counter);
			}
			timesCounter++;
			if(timesCounter==nTimes)
				break;
		}
		if( t>= task.getMaxTime()) {
			break;
		}
		t += tau;
		counter++;
		/* Cutting if(t > tCut) {
			dumpToFile(t);
			saveSolution("ms-alpha.dat", t);
			break;
		}*/
	} 
	
	/* Cutting
	/////////////////DEBUG////////////////////////////////////////////////
	// Здесь отрезаем nCut ячеек и считаем дальше с другим УРС
	if(! (task.getMethodFlag() == 0) ) {
		cout << endl << "Calculation finished!" << endl;
		// Если метод Эйлеров, то надо удалить сетку
		if(task.getMethodFlag() == 1)
			(task.getMethod()).deleteGrid();
		exit(1); 
	}
		
	cout << "Now cutting calculation region and changing EOS to fe_alpha" << endl;
	cout << "New left node is node number " << nCut << endl;
	// Устанавливаем новый УРС
	EOS* newEOS = new EOSTableFeAlpha("fe_alpha", 1, 7874.);
	task.setEOS(newEOS);
	// Устанавливаем новый nSize, новый ms
	ms.loadData("ms-alpha.dat", nCut);
	ms_temp.loadData("ms-alpha.dat", nCut);
	for(;;) {
		tau = calcTimeStep(t);
		cout << counter << ": " << "t=" << t <<  " tau=" << tau << " courant=" << getzKur() << endl;
		if(task.getHydroStage() && (!getSpallFlag())) 
			calcHydroStage(t, tau);
		else if (task.getHydroStage() && (getSpallFlag()))
			calcHydroStageSpallation(t, tau);
		if(task.getHeatStage() && (!getSpallFlag())) 
			calcHeatStage(t, tau);
		else if (task.getHeatStage() && (getSpallFlag()))
			calcHeatStageSpallation(t, tau);
		if(task.getExchangeStage()) calcExchangeStage(tau);
		if(task.getIonizationStage()) calcIonizationStage(t, tau);	
		if(handleKeys(t)) break;
		double xMeltL=1.0e-6;
		double xMeltR=1.0e-6;
		for(i=0; i<ms.getSize(); i++) {
			if( (task.getEOS().getphase(ms[i].ro, ms[i].ti) > 1.0)&&
			(task.getEOS().getphase(ms[i].ro, ms[i].ti) < 3.0) ) {
				xMeltL = ms[i].x;			
				break;
			}
		}
		for(i=0; i<ms.getSize(); i++) {
			if( (task.getEOS().getphase(ms[i].ro, ms[i].ti) > 1.0)&&
			(task.getEOS().getphase(ms[i].ro, ms[i].ti) < 3.0) )
				xMeltR = ms[i].x;			
		}
		//Energy conservation
		double eInner=0.0, eKinetic=0.0, eFull=0.0, deltae=0.0;
		for(int i1=0; i1<ms.getSize(); i1++) {
			eInner+=ms[i1].e*ms[0].dm;
			eKinetic+=ms[0].dm*ms[i1].v*ms[i1].v/2.0;
			eFull=eInner+eKinetic;
			deltae=eFull-eInit;
		}
		cout << "Energy is: inner " << eInner << " kinetic " << eKinetic << " full " << eFull << endl;
		cout << "de = " << deltae << endl;
		if(counter%10 == 0) {   
			if ((task.getSourceFlag()==1)||(task.getSourceFlag()==2)||(task.getSourceFlag()==3)) {
				f=fopen(rEdgeFileName, "a+");
				fprintf(f, "%e %e %e %e %e %f\n", t*1.0e12, ms[ms.getSize()].x*1.0e9, ms[ms.getSize()].v, xMeltL*1.0e9, xMeltR*1.0e9, deltae);
				fclose(f);
			}
		}
		// Regular file output
		if (t>=timesArray[timesCounter]) {
		if(task.getMethodFlag() == 0) {
			dumpToFile(t); 
		} else {
			dumpToFileEuler(t);
		}
		timesCounter++;
		if(timesCounter==nTimes)
			break;
		}
		if( t>= task.getMaxTime()) {
			break;
		}
		t += tau;
		counter++;
	}
	dumpToFile(t);
	*/

	cout << endl << "Calculation finished!" << endl;
	// Если метод Эйлеров, то надо удалить сетку
	if(task.getMethodFlag() == 1)
		(task.getMethod()).deleteGrid();
}

void CSolver::goAuSpall(char *fName) {
	// Сначала отсекаем все лишние конфигурации нам нужно только золото, только лагранжева сетка и схема Самарского
	if(task.getMethodFlag()!= 0) {
		cout << "goAuSpall() error: only lagrangian coordinates are avilable." << endl;
		cout << "Change \'method\' parameter to \'euler\' in task.txt file or use another go...() function for your simulation. " << endl;
		exit(1);
	}
	task.load(fName);
	// Если УРС не табличное, нам тоже не сюда.
	if(task.getEOS().getType()!=table) {
		cout << "goAuSpall() error: only table EOS is available.";
	}
	ms.initData(&task);
	ms_temp.initData(&task);
	initVars();
	int nSize = ms.getSize();
	double t = -3.*task.getTauPulse(), tau = 0.;
	int i = 0, counter = 0;
	dumpToFileTestRP(t, counter);
	/*double timesArray[] = {      .0,   1.e-13,   3.e-13,   5.e-13,   1.e-12,  2.e-12,   3.e-12,   4.e-12,   5.e-12,   6.e-12,
		                     7.e-12,  10.e-12,  20.e-12,  30.e-12,  40.e-12,  50.e-12,  60.e-12,  70.e-12,  80.e-12,  90.e-12,
						   100.e-12, 150.e-12, 2.e-12};*/
	double timesArray[] = {      .0,   1.e-9,   2.e-9,   5.e-9,   10.e-9,  20.e-9,   50.e-9,   100.e-9,   200.e-9,   500.e-9,
		                      1.e-6,   2.e-6,   3.e-6,   4.e-6};


	int nTimes = 14;
	int timesCounter = 0;
	// Протестируем последовательность, на строгое возрастание (выход за границу мы не контролируем)
	for(int i = 0; i<nTimes-1; i++)
		if(timesArray[i+1] <= timesArray[i]) {
			cout << "CSolver::go() error: Time points unordered in timesArray." << endl;
			exit(1);
		}
	// EOS control values testing
	//testEOSControlNumbers(19301., 1000., 3000.);
	// Preparing file for right end data writing
	char rEdgeFileName[255];
	FILE* f=0;
	strcpy(rEdgeFileName, OUTPUT_FOLDER);
	strcat(rEdgeFileName, "right_");
	strcat(rEdgeFileName, task.getTaskName().c_str());
	strcat(rEdgeFileName, ".dat");
	f=fopen(rEdgeFileName, "w");
	fprintf(f, "TITLE=\"Right edge trajectory\"\n");
	fprintf(f, "VARIABLES=\"t [ps]\",\"x_r [nm]\",\"v_r [m/s]\",\"xMeltL [nm]\",\"xMeltR [nm]\", \"deltae [J/m2]\" \n");
	fclose(f);
	// Left edge trajectory watching file
	string fTName = OUTPUT_FOLDER + string("t_fronts_") + task.getTaskName() + ".dat";
	ofstream fT;
	fT.open(fTName, fstream::out);
	fT << "TITLE=\"Temperature fronts trajectory\"\n";
	fT << "VARIABLES=\"t[ps]\",\"xL[nm]\",\"xR[nm]\"\n";
	fT.close();
	// Calculating total inner energy for testing
	double eInit=0.0;
	for(int i2=0; i2<ms.getSize(); i2++) {
		eInit+=ms[i2].e*ms[0].dm;
	}
	// Number of a cell and moment of time in which spallation takes place
	int nCut = 77;
	double tCut = 300.e-12;
	for(;;)
	{
		// Заметки на полях. При первой же генеральной уборке кода. Cписок вот тут в комментариях: CSolver::solveti()
		// Вот этот фрагмент кода включаем, если считаем лагранжевым методом.
		tau = calcTimeStep(t);
		cout << counter << ": " << "t=" << t <<  " tau=" << tau << " courant=" << getCFL() << endl;
		if(task.getHydroStage()) calcHydroStage(t, tau);
		if(task.getHeatStage()) calcHeatStage(t, tau);
		if(task.getExchangeStage()) calcExchangeStage(tau);
		// Keyboard input tracking
		if(handleKeys(t)) break;
		// Melting fronts tracking
		double xMeltL=1.0e-6;
		double xMeltR=1.0e-6;
		for(i=0; i<ms.getSize(); i++) {
			if( (task.getEOS().getphase(ms[i].ro, ms[i].ti) > 1.0)&&
 			(task.getEOS().getphase(ms[i].ro, ms[i].ti) < 3.0) ) {
				xMeltL = ms[i].x;			
				break;
			}
		}
		for(i=0; i<ms.getSize(); i++) {
			if( (task.getEOS().getphase(ms[i].ro, ms[i].ti) > 1.0)&&
			(task.getEOS().getphase(ms[i].ro, ms[i].ti) < 3.0) )
				xMeltR = ms[i].x;			
		}
		// Calculating energy for conservation test
		double eInner=0.0, eKinetic=0.0, eFull=0.0, deltae=0.0;
		for(int i1=0; i1<ms.getSize(); i1++) {
			/*eInner+=ms[i1].e*ms[0].dm;
			eKinetic+=ms[0].dm*ms[i1].v*ms[i1].v/2.0;*/
			eInner += (ms[i1+1].x-ms[i1].x)*ms[i1].ro*ms[i1].e;
			eKinetic += 0.5 * ms[i1].ro*ms[i1].v*ms[i1].v;
			eFull=eInner+eKinetic;
			deltae=eFull-eInit;
		}
		cout << "Energy is: inner " << eInner << " kinetic " << eKinetic << " full " << eFull << endl;
		cout << "de = " << deltae << endl;
		if(counter%10 == 0) {   
			f=fopen(rEdgeFileName, "a+");
			fprintf(f, "%e %e %e %e %e %f\n", t*1.0e12, ms[ms.getSize()].x*1.0e9, ms[ms.getSize()].v, xMeltL*1.0e9, xMeltR*1.0e9, deltae);
			fclose(f);
			double xL=0., xR=0.;
			for(i=0; i<nSize; i++) 
				if(ms[i].ti >=400.) {
					xL=ms[i].x;
					break;
				}			
			for(i=nSize*9/10; i>=0; i--)
				if(ms[i].ti >=400.) {
					xR=ms[i].x;
					break;
				}
			fT.open(fTName, fstream::app);
			fT << t*1.e12 << " " << xL*1.e9 << " " << xR*1.e9 << endl;
			fT.close();
		}
		// Regular file output
		if (t>=timesArray[timesCounter]) {
			if(task.getMethodFlag() == 0) {
				dumpToFile(t); 
			} else {
				dumpToFileEuler(t);
				dumpToFileTestRP(t, counter);
			}
			timesCounter++;
			//if(timesCounter==nTimes)
				//break;
		}
		if( t>= task.getMaxTime()) {
			break;
		}
		t += tau;
		counter++;
	}
	cout << endl << "Calculation finished!" << endl;
}

void CSolver::goGlass(char* fName) {
	// Загрузка входного файла, инициализация объекта task типа CTask	
	// task.type = TaskType::RuGlass; // Пока не навел порядок в инфраструктуре, проставляю этот флаг руками
	task.load(fName);
	// Наполнение объектов ms, ms_temp данными на по начальным и граничным условиям из task
	assert(task.getMethodFlag() == MethodType::samarskii);
	ms.initData(&task);
	ms_temp.initData(&task);
	initVars();
	double t   = -3*tauPulse,  tau = 0.0;
	int	counter = 0, i=0;
	double timesArray[] = {0.,        2.e-12,  5.e-12,   10.e-12,   15.e-12, 16.e-12, 17.e-12,  18.e-12,  19.e-12,  20.e-12, 
		                   21.e-12,  22.e-12, 23.e-12, 23.87e-12, 23.88e-12, 25.e-12,  30.e-12, 50.e-12, 100.e-12, 200.e-12, 
						   500.e-12,   5.e-9,  50.e-9,   500.e-9};
	int nTimes = 24, timesCounter = 0;
	dumpToFile(t);
	int itNumHydro = 0, itNumHeat = 0, itNumExchg = 0;
	for(;;)	{
		ostringstream oss;
		tau = calcTimeStep(t);			
		oss << counter << ": " << "t=" << setprecision(6) << t*1.e12 <<  "ps tau=" << tau*1.e12 << "ps CFL=" << getCFL() << " ";
		if(task.getHydroStage()) {
			if(task.type == TaskType::auWater) itNumHydro = calcHydroStageGlass(t, tau); else itNumHydro = calcHydroStage(t, tau);			
			oss << "Hydro:" << itNumHydro << " "; 
		}
		if(task.getHeatStage()) { 
			if(task.type == TaskType::auWater) itNumHeat = calcHeatStageGlass(t, tau); else itNumHeat = calcHeatStage(t, tau);		 
			oss << "Heat:" << itNumHeat << " "; 
		}
		if(task.getExchangeStage()) { 
			if(task.type == TaskType::auWater) itNumHeat = itNumExchg = calcExchangeStageGlass(tau); else itNumHeat = itNumExchg = calcExchangeStage(tau);		 			
			oss << "Exchg:" << itNumExchg << " "; 
		}
		oss << "(iters)";
		if(handleKeys(t)) break;
		// Regular file output
		if (t>=timesArray[timesCounter]) {
			dumpToFile(t); 
			timesCounter++;
			if(timesCounter==nTimes)
				break;
		}		
		if( t>= task.getMaxTime()) 
			break;		
		t += tau;
		counter++;
		cout << oss.str() << endl; 		
	} 
	cout << endl << "Calculation finished!" << endl;
}

/*
double CSolver::getEntropy(double ro, double ti, double te) {
	double ci = task.getEOS().getci(ro, ti);
	double  p = task.getEOS().getp(ro, ti, te);
	if(p<0) {
		cerr << "Error: CSolver::getEntropy(): negative pressure!" << endl;
		exit(1);
	}
	double ro_gamma = pow(ro, (task.getEOS()).getGamma());
	if(ro!=0.)
		return ci*log(p/ro_gamma);
	else 
		return 0.;
}*/

double CSolver::calcTimeStep(double t) {
	if( task.getSourceFlag()==2 && (tauPulse<1.e-12) && (t<5.0e-12)) {
		return 1.0e-15;
	} else if( task.getSourceFlag()==1 && (tauPulse<1.e-12) && (t<5.0e-12)) {
		if(task.type!=TaskType::ruGlass) 
		  return 1.0e-15;
		else
		  return .5e-15;
	}
	double tau_temp1 = (ms[1].x - ms[0].x) / ms[0].C;
	double tau_temp2;
	for(int i=1; i<ms.getSize()-1; i++)
	{
		tau_temp2 = (ms[i+1].x - ms[i].x) / ms[i].C;
		if(tau_temp1 > tau_temp2)
			tau_temp1 = tau_temp2;
	}
	return CFL * tau_temp1;
}


double CSolver::calcTimeStepEuler(double t) 
{
	double v_max = max(fabs(ms[0].v), max(fabs(ms[0].v - ms[0].C), fabs(ms[0].v + ms[0].C))); 
	
	double tau_temp1 = (ms[1].x - ms[0].x) / v_max;
	double tau_temp2 = 0.;

	CMethod& method = task.getMethod();

	for(int i=1; i<ms.getSize()-1; i++)
	{
		//v_max = max(fabs(ms[i].v), max(fabs(ms[i].v - ms[i].C), fabs(ms[i].v + ms[i].C))); 
		v_max = fabs(ms[i].v) + ms[i].C;   //, fabs(method.vGrid[i]) +  ms[i].C);
		tau_temp2 = (ms[i+1].x - ms[i].x) / v_max;

		if(tau_temp1 > tau_temp2)
			tau_temp1 = tau_temp2;
	}

	return CFL * tau_temp1;
}


bool CSolver::handleKeys(double t)
{
#ifdef _WIN32
	if(_kbhit())
		switch(toupper(_getch()))
		{
			case 'Q':
				return true;

			case 'D':
				if(task.getMethodFlag() == MethodType::samarskii) {
				    dumpToFile(t);
				} else {
					dumpToFileEuler(t);
				}
				break;
		}
#endif
	return false;
}


double CSolver::compTe()
{
	double rValue1 = fabs(ms[0].te_temp - ms[0].te);
	double rValue2;

	for(int i=1; i<ms.getSize(); i++)
	{
		rValue2 = fabs(ms[i].te_temp - ms[i].te);

		if(rValue1 < rValue2)
			rValue1 = rValue2;
	}

	return rValue1;
}


double CSolver::compTi()
{
	double rValue1 = fabs(ms[0].ti_temp - ms[0].ti);
	double rValue2 = 0;

	for(int i=1; i<ms.getSize(); i++)
	{
		rValue2 = fabs(ms[i].ti_temp - ms[i].ti);

		if(rValue1 < rValue2)
			rValue1 = rValue2;
	}

	return rValue1;
}

double CSolver::compTi(CFieldOld &ms1, CFieldOld &ms2)
{
	double rValue1 = fabs(ms1[0].ti - ms2[0].ti);
	double rValue2 = 0;
	for(int i=1; i<ms.getSize(); i++) {
		rValue2 = fabs(ms1[i].ti - ms2[i].ti);
		if(rValue1 < rValue2)
			rValue1 = rValue2;
	}
	return rValue1;
}


double CSolver::compTe(CFieldOld &ms1, CFieldOld &ms2)
{
	double rValue1 = fabs(ms1[0].te - ms2[0].te);
	double rValue2;

	for(int i=1; i<ms.getSize(); i++)
	{
		rValue2 = fabs(ms1[i].te - ms2[i].te);

		if(rValue1 < rValue2)
			rValue1 = rValue2;
	}

	return rValue1;
}

double CSolver::compE(CFieldOld &ms1, CFieldOld &ms2)
{
	Node* n1=ms1.getnodes();
	Node* n2=ms2.getnodes();

	double rValue1 = fabs(n1[0].e - n2[0].e);
	double rValue2;


	for(int i=1; i<ms.getSize(); i++)
	{
		rValue2 = fabs( n1[i].e - n2[i].e);

		if(rValue1 < rValue2)
			rValue1 = rValue2;
	}

	return rValue1;
}

double CSolver::compEi(CFieldOld &ms1, CFieldOld &ms2)
{
	double rValue1 = fabs(ms1[0].ei - ms2[0].ei);
	double rValue2;

	for(int i=1; i<ms.getSize(); i++)
	{
		rValue2 = fabs(ms1[i].ei - ms2[i].ei);

		if(rValue1 < rValue2)
			rValue1 = rValue2;
	}

	return rValue1;
}

void CSolver::dumpToFile(double t) {
	double p_an, ro_an, v_an;
	char buf[40];
	char frac_buf[40];
	char fName[256];
	int tConv = (int)(fabs(t*1.0e12));
    sprintf(buf, "%d", tConv);
	buf[7]='\0';
	if(tConv < 10) {
		buf[6]=buf[0];
		buf[5]='0';
		buf[4]='0';
		buf[3]='0';
		buf[2]='0';
		buf[1]='0';
		buf[0]='0';
	}
	else if(tConv < 100) {
		buf[6]=buf[1];
		buf[5]=buf[0];
		buf[4]='0';
		buf[3]='0';
		buf[2]='0';
		buf[1]='0';
		buf[0]='0';
	}
	else if(tConv < 1000) {
		buf[6]=buf[2];
		buf[5]=buf[1];
		buf[4]=buf[0];
		buf[3]='0';
		buf[2]='0';
		buf[1]='0';
		buf[0]='0';
	}
	else if(tConv < 10000) {
		buf[6]=buf[3];
		buf[5]=buf[2];
		buf[4]=buf[1];
		buf[3]=buf[0];
		buf[2]='0';
		buf[1]='0';
		buf[0]='0';
	}
	else if(tConv < 100000) {
		buf[6]=buf[4];
		buf[5]=buf[3];
		buf[4]=buf[2];
		buf[3]=buf[1];
		buf[2]=buf[0];
		buf[1]='0';
		buf[0]='0';
	}
	else if(tConv < 1000000) {
		buf[6]=buf[5];
		buf[5]=buf[4];
		buf[4]=buf[3];
		buf[3]=buf[2];
		buf[2]=buf[1];
		buf[1]=buf[0];
		buf[0]='0';
	}
	int tFrac = abs((int)(t*1.0e13) % 10);
    sprintf(frac_buf, "%d", tFrac);
	strcpy(fName, OUTPUT_FOLDER);
	if(t<0) {
		strcat(fName, "-");
	} else {
		strcat(fName, "");
	}
	strcat(fName, buf); strcat(fName, "."); strcat(fName, frac_buf); strcat(fName, "ps_"); strcat(fName, task.getTaskName().c_str()); strcat(fName, ".dat");
	printf("Dump matter to file: %s\n", fName);
	EOSType EType = EOSType::none;
	if(task.getSourceFlag()!=4) EOSType EType = task.getEOS().getType(); 
	FILE* f=fopen(fName, "w");
	if(!f) {
		cout << "Cannot open output file " << fName << ". Sorry." << endl;
		cin.get();
		exit(1);
	}
	if( (EType == ideal)||(EType == test) ) {
		fprintf(f, "TITLE=\"Physical quantities' profiles\"\n");
		fprintf(f, "VARIABLES=\"x\",\"RO\",\"v\",\"Pe\",");
		fprintf(f, "\"Pi\",\"P\",\"Ee\",\"Ei\",\"E\",\"Te\",");
		fprintf(f, "\"Ti\",\"C\",\"P_an\",\"v_an\",\"RO_an\"\n");
	}
	else {
		fprintf(f, "TITLE=\"Physical quantities' profiles\"\n");
		fprintf(f, "VARIABLES=\"x [nm]\",\"RO [kg/m3]\",\"v [km/s]\",\"Pe [GPa]\",");
		fprintf(f, "\"Pi [GPa]\",\"P [GPa]\",\"Ee [MJ/kg]\",\"Ei [MJ/kg]\",\"E [MJ/kg]\",\"Te [K]\",");
		fprintf(f, "\"Ti [K]\",\"C [km/s]\", \"Phase\", \"Z\", \"kappa[J/m/K]\" \n");
	}
	double mul_x=1.0e9, mul_p=1.0e-9, mul_e=1.0e-6, mul_v=1.0e-3;
	double x = 0., v = 0., phase = 0.;
	for(int i=0; i<ms.getSize(); i++) {
		Node &n = ms[i];		
		if( (EType == ideal)||(EType == test) ) {
			getLagrangeAnalyticApproximation(i, t, &p_an, &v_an, &ro_an);			
			x = 0.5 * (n.x + ms[i+1].x); v = 0.5 * (n.v + ms[i+1].v);
			fprintf(f, "%e %e %f %f %f %f %f %f %f %f %f %e %e %e %e\n",
					x *mul_x,  n.ro, v *  mul_v,
					n.pe*mul_p, n.pi*mul_p, n.p*mul_p,			
					n.ee, n.ei, n.e,
					n.te, n.ti, 
					n.C, p_an, v_an, ro_an);
		} else {
			x = 0.5 * (n.x + ms[i+1].x); v = 0.5 * (n.v + ms[i+1].v);
			if (EType == EOSType::none) phase = 1.; else phase = task.getEOS().getphase(n.ro, n.ti);
			fprintf(f, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", 
					x *  mul_x, n.ro,         v*mul_v,
					n.pe * mul_p, n.pi * mul_p, n.p * mul_p,			
					n.ee * mul_e, n.ei * mul_e, n.e * mul_e,
					n.te,  n.ti,  n.C  * mul_v, phase, /*1.0-mix*/ n.Z, n.kappa); 
		}
	}
	fclose(f);
}


void CSolver::dumpToFileTestRP(double t, int num) {
	char strNum[10];
	_itoa(num, strNum, 10);
	string fName = OUTPUT_FOLDER + string("out-RP-") + strNum + ".dat";
	printf("Dump matter to file: %s\n", fName.c_str());
	EOSType EType = EOSType::none;
	if(!task.eosBin)
		EType = task.getEOS().getType();
	else 
		EType = EOSType::bin;
	FILE* f=fopen(fName.c_str(), "w");
	if(!f) {
		cout << "Cannot open output file. Sorry." << endl;
		cin.get();
		exit(1);
	}
	if( (EType == ideal)||(EType == test) )	{		
		fprintf(f, "TITLE=\"Physical quantities' profiles\"\n");
		fprintf(f, "VARIABLES=\"x\",\"ro\",\"v\",\"p\",\"e\",\"s\",\"s_an\",\"IL\",\"IR\",\"IL_an\",\"IR_an\",\"ro_an\",\"v_an\",\"p_an\",\"e_an\",\"xi\",\"T\"\n");
	} else if (EType == EOSType::bin) {
		fprintf(f, "TITLE=\"RP test for binary eos, t=%f\"\n", t);
		fprintf(f, "VARIABLES=\"x\",\"ro\",\"v\",\"p\",\"e\",\"ro_an\",\"v_an\",\"p_an\",\"e_an\",\"xi\"\n");
	} else {
		cout << "Error CSolver::dumpToFileRP(): improper EOS type!" << endl;
		cin.get();
		exit(1);
	}
	CVectorPrimitive q;
	EOSOld* eos = &(task.getEOS());
	int nZones = task.getNumZones();
	double roL = 0., vL = 0., pL = 0., roR = 0., vR = 0., pR = 0., x0 = 0., e = 0., s = 0., IL = 0., IR = 0., IL_an = 0., IR_an = 0, T= 0.;
	double gamma = 0.; 
	if(EType != EOSType::bin)
		gamma = eos->getGamma();
	else
		gamma = task.eosBin->gamma;
	// If left zone is vacuum zone
	if(nZones == 1) { 
		roL = 0.; roR = task.getZone(0).ro;
		 vL = 0.;  vR = task.getZone(0).v;
		 pL = 0.;  
		if(EType != EOSType::bin) 
			pR = eos->getpi(roR, task.getZone(0).ti);
		else
		    pR = task.eosBin->getp(roR, task.getZone(0).e);		 
		 x0 = 0.;
	} else {
	// If two non-vacuum zones are present in IVP (initial value problem)
		roL = task.getZone(0).ro;				  roR = task.getZone(1).ro;
		 vL = task.getZone(0).v;				   vR = task.getZone(1).v;
		if(EType != EOSType::bin) {
			pL = eos->getpi(roR, task.getZone(0).ti);
			pR = eos->getpi(roR, task.getZone(1).ti);
		}
		else {
		    pL = task.eosBin->getp(roL, task.getZone(0).e);				 
			pR = task.eosBin->getp(roR, task.getZone(1).e);				 		 
		}		
	 	x0 = task.getZone(0).l;
	}
	double x = 0.0, xi = 0., sAn = 0.;
	for(int i=0; i<ms.getSize(); i++)	{
		Node &n = ms[i];		
		x = 0.5 * (n.x + ms[i+1].x);
		if(t!=0.) xi = x/t; else xi = 0.;
		if( (EType == ideal)||(EType == test) )	{		
			q = calcRPAnalyticalSolution(roL, vL, pL, roR, vR, pR, x-x0, t);
			if(q.ro!=0.) e = q.p/(gamma-1.)/q.ro; else e = 0.;			
			s=eos->getEntropy(ms[i].ro, ms[i].ti);			
			sAn=eos->getEntropy(q.ro, eos->getti(q.ro, q.p/(gamma-1.)/q.ro));			
			IL = ms[i].v + 2.0*ms[i].C/(gamma-1.);
			IR = ms[i].v - 2.0*ms[i].C/(gamma-1.);
			IL_an = q.v + 2.0*sqrt(gamma*q.p/q.ro)/(gamma-1);
			IR_an = q.v - 2.0*sqrt(gamma*q.p/q.ro)/(gamma-1);
			T = eos->getti(ms[i].ro, ms[i].p/(gamma-1.)/ms[i].ro);
			fprintf(f, "%e %f %f %f %f %e %e %f %f %f %f %f %f %f %f %f %e\n", x, n.ro, n.v, n.p, n.e, s, sAn, IL, IR, IL_an, IR_an, q.ro, q.v, q.p, e, xi, T);
		}
		else if (EType == EOSType::bin) {
			q = calcRPAnalyticalSolutionEOSBin(roL, vL, pL, roR, vR, pR, x-x0, t);
			e = task.eosBin->gete(q.ro, q.p);
			fprintf(f, "%e %f %f %f %f %f %f %f %f %e\n", x, n.ro, n.v, n.p, n.e, q.ro, q.v, q.p, e, xi);
		}
	}
	fclose(f);
}

void CSolver::dumpToFileEuler(double t)
{
	char buf[40];
	char frac_buf[40];
	char fName[256];
	double p_an=0., ro_an=0., v_an=0.;
	double orderMul = 1.;
	int tConv = (int)(fabs(t*orderMul));
    sprintf(buf, "%d", tConv);	
	buf[6]='\0';

	if(tConv < 10)
	{
		buf[5]=buf[0];
		buf[4]='0';
		buf[3]='0';
		buf[2]='0';
		buf[1]='0';
		buf[0]='0';
	}
	else if(tConv < 100)
	{
		buf[5]=buf[1];
		buf[4]=buf[0];
		buf[3]='0';
		buf[2]='0';
		buf[1]='0';
		buf[0]='0';
	}
	else if(tConv < 1000)
	{
		buf[5]=buf[2];
		buf[4]=buf[1];
		buf[3]=buf[0];
		buf[2]='0';
		buf[1]='0';
		buf[0]='0';
	}
	else if(tConv < 10000)
	{
		buf[5]=buf[3];
		buf[4]=buf[2];
		buf[3]=buf[1];
		buf[2]=buf[0];
		buf[1]='0';
		buf[0]='0';
	}
	else if(tConv < 100000)
	{
		buf[5]=buf[4];
		buf[4]=buf[3];
		buf[3]=buf[2];
		buf[2]=buf[1];
		buf[1]=buf[0];
		buf[0]='0';
	}

	int tFrac = abs((int)(t*orderMul*10.) % 10);
    sprintf(frac_buf, "%d", tFrac);
	strcpy(fName, OUTPUT_FOLDER);
	if(t<0)
		strcat(fName, "-");
	else
		strcat(fName, "");

	strcat(fName, buf);
	strcat(fName, ".");
	strcat(fName, frac_buf);
	strcat(fName, "_");
	strcat(fName, task.getTaskName().c_str());
	strcat(fName, ".dat");
	printf("Dump matter to file (euler): %s\n", fName);

	EOSType EType = task.getEOS().getType();
	FILE* f=fopen(fName, "w");
	fprintf(f, "TITLE=\"Physical quantities' profiles\"\n");
	fprintf(f, "VARIABLES=\"x [nm]\",\"RO [kg/m3]\",\"v [km/s]\",\"Pe [GPa]\",");
	fprintf(f, "\"Pi [GPa]\",\"P [GPa]\",\"Ee [MJ/kg]\",\"Ei [MJ/kg]\",\"E [MJ/kg]\",\"Te [K]\",");
	fprintf(f, "\"Ti [K]\",\"C [km/s]\", \"Phase\", \"liqMix\", \"Z\", \"kappa[J/m/K]\", \"ro_an[kg/m^3]\", \"v_an[m/s]\", \"P_an[GPa]\" \n");


	double mul_x=1.0e9, mul_p=1.0e-9, mul_e=1.0e-6, mul_v=1.0e-3;
	for(int i=0; i<ms.getSize(); i++) {
		Node &n = ms[i];
	 	double x = 0.5*(n.x + ms[i+1].x);
		double v = n.v;
		double phase = task.getEOS().getphase(n.ro, n.ti);
		double   mix = task.getEOS().getmix(n.ro, n.ti);
		getEulerAnalyticApproximationGrid(i, x, t, &p_an, &v_an, &ro_an);
		fprintf(f, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", 
					x *  mul_x, n.ro,         v*mul_v,
					n.pe * mul_p, n.pi * mul_p, n.p * mul_p,			
					n.ee * mul_e, n.ei * mul_e, n.e * mul_e,
					n.te,  n.ti,  n.C  * mul_v, phase, 1.0-mix, n.Z, n.kappa, ro_an, v_an * mul_v, p_an * mul_p); 
	}
	fclose(f);
	strcpy(fName, OUTPUT_FOLDER);
	strcat(fName, "progress.log");
	f=fopen(fName, "w");
	fprintf(f, "%e", t); 
	fclose(f);
}

void CSolver::setCFL(double _CFL) { CFL = _CFL; }

double CSolver::getCFL(void) { return CFL; }

void CSolver::saveSolution(const char* fName, double t) {
	char fullName[256];
	strcpy(fullName, OUTPUT_FOLDER);
	strcat(fullName, fName);
	unsigned int i=0, nSize = ms.getSize();
	printf("Saving solution at t=%e to file: %s\n", t, fullName);
	string fLogName = OUTPUT_FOLDER + string("saveSolution.log");
	FILE* f=fopen(fullName, "w");
	fprintf(f, "TITLE = \"Solution profile for t=%e\"\n", t);
	fprintf(f, "VARIABLES = \"x\", \"u\", \"ro\", \"T\", \"Te\", \"pi\", \"pe\", \"p\", \"ei\", \"ee\", \"e\"\n");
	for(i=0; i<nSize; i++){
		Node &n = ms[i];
		fprintf(f, "%e %e %e %e %e %e %e %e %e %e %e\n", n.x, n.v, n.ro, n.ti, n.te, n.pi, n.pe, n.p, n.ei, n.ee, n.e);
	}
	i=ms.getSize();
	Node &n = ms[i];
	fprintf(f, "%e %e %e %e %e %e %e %e %e %e %e\n", n.x, n.v, -1., -1., -1., -1., -1., -1., -1., -1., -1.);
	fclose(f);
}


// Эту функцию надо переписать под новый saveSolution()
double CSolver::loadSolution(char* fName)
{
	printf("Loading solution from file: %s\n", fName);
	double t = 0.0; 
	int i=0, nSize = 0;	
	char buf[256]; 
	buf[0] = 0;
	strcat(buf, "calc/output/");
	strcat(buf, fName);
	
	FILE* f=fopen(buf, "r");
	if(!f)  
		printf("Solution loading error!\n");
	
	if (fscanf(f, "%s", buf) != 1) 
		printf("Syntax error!\n");
	nSize = atoi(buf);
	char buf01[20], buf02[20], buf03[20], buf04[20], buf05[20], buf06[20], buf07[20], buf08[20], buf09[20], buf10[20],
		 buf11[20], buf12[20], buf13[20], buf14[20], buf15[20], buf16[20], buf17[20], buf18[20], buf19[20], buf20[20], buf21[20];


	for(i=0; i<ms.getSize(); i++)
	{
		Node &n = ms[i];
		if(fscanf(f, "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s",
				  buf01, buf02, buf03, buf04, buf05, buf06, buf07, buf08, buf09, buf10,
			      buf11, buf12, buf13, buf14, buf15, buf16, buf17, buf18, buf19, buf20, buf21
				   /* n.x,  n.v, n.ro, n.dm, 
					n.ti, n.te, n.pi, n.pe, n.p, n.ei, n.ee, n.e, 
					n.ci, n.ce, n.C,  n.Alphaei, n.kappa, n.ne, n.Z, n.ti_temp, n.te_temp*/)!=21)
			printf("Syntax error!\n");
		n.x=atof(buf01);  n.v=atof(buf02); n.ro=atof(buf03); n.dm=atof(buf04); n.ti=atof(buf05); n.te=atof(buf06); n.pi=atof(buf07); 
		n.pe=atof(buf08); n.p=atof(buf09); n.ei=atof(buf10); n.ee=atof(buf11); n.e=atof(buf12);  n.ci=atof(buf13); n.ce=atof(buf14);
		n.C=atof(buf15);  n.Alphaei=atof(buf16); n.kappa=atof(buf17); n.ne=atof(buf18); n.Z=atof(buf19); n.ti_temp=atof(buf20); 
		n.te_temp=atof(buf21);
	}
	i = ms.getSize();
	Node &n = ms[i];
	if(fscanf(f, "%s %s %s", buf01, buf02, buf03/*n.x, n.v, n.dm*/)!=3)
		printf("Syntax error!\n");
	n.x=atof(buf01); n.v=atof(buf02); n.dm=atof(buf03);
	if(fscanf(f, "%s\n", /*t*/buf01)!=1)
		printf("Syntax error!\n");
	t = atof(buf01);

	fclose(f);

	return t;
}


void CSolver::getEulerAnalyticApproximation(int i, double t, double *p_an, double *v_an, double *ro_an)
{
	double x0, c0, p0, ro0;
	Zone *zone;

	if(task.getNumZones() < 2 || t==0)
	{
		*v_an  = 0;
		*p_an  = 0;
		*ro_an = 0;
		//*F3_an = 0;
		return;
	}

	zone = &task.getZone(0);
	x0  = zone->l;

	zone = &task.getZone(1);
	c0  = task.getEOS().getC(zone->ro, zone->ti, zone->te);
	p0  = task.getEOS().getp(zone->ro, zone->ti, zone->te);
	ro0 = zone->ro;

	const double gamma = 5.0/3.0;

	if (ms[i].x <= x0 + c0 * t)
	{
		*v_an = -2.0 / (gamma+1.0) * (c0-(ms[i].x-x0)/t);

		double tmp = 1.0 - fabs(*v_an) / (3.0*c0);

		*p_an  = p0  * pow(tmp, 5);
		*ro_an = ro0 * pow(tmp, 3);

		//c_an = c0 - 1.0/3.0 * fabs(v_an);
		//*F3_an = p_an / (2.0*c_an);
	}
	else
	{
		*v_an  = 0;
		*p_an  = p0;
		*ro_an = ro0;
		//*F3_an = 0; 
	}
}


void CSolver::getLagrangeAnalyticApproximation(int i, double t, double *p_an, double *v_an, double *ro_an)
{
	double x0, c0, p0, ro0;
	Zone *zone;

	x0  = 0;

	zone = &task.getZone(0);
	c0  = task.getEOS().getC(zone->ro, zone->ti, zone->te);
	p0  = task.getEOS().getp(zone->ro, zone->ti, zone->te);
	ro0 = zone->ro;

	const double gamma = 5.0/3.0;

	double x = 0.5 * (ms[i].x + ms[i+1].x);

	if ( (x <= x0 + c0 * t) && (t!=0.0) )
	{
		*v_an = -2.0 / (gamma+1.0) * (c0-(x-x0)/t);

		double tmp = 1.0 - fabs(*v_an) / (3.0*c0);

		*p_an  = p0  * pow(tmp, 5);
		*ro_an = ro0 * pow(tmp, 3);

		//c_an = c0 - 1.0/3.0 * fabs(v_an);
	}
	else
	{
		*v_an  = 0;
		*p_an  = p0;
		*ro_an = ro0;
	}
}

void CSolver::getEulerAnalyticApproximationGrid(int i, double x, double t, double *p_an, double *v_an, double *ro_an)
{
	Zone *zone         = &task.getZone(0);
	double  x0=0., c0  =  task.getEOS().getC(zone->ro, zone->ti, zone->te), 
		           p0  =  task.getEOS().getp(zone->ro, zone->ti, zone->te), 
		  		   ro0 =  zone->ro;
	const double gamma =  5.0/3.0;
	if ( (x <= x0 + c0 * t) && (t!=0.0) )
	{
		*v_an = -2.0 / (gamma+1.0) * (c0-(x-x0)/t);
		double tmp = 1.0 - fabs(*v_an) / (3.0*c0);
		*p_an  = p0  * pow(tmp, 5);
		*ro_an = ro0 * pow(tmp, 3);
		//c_an = c0 - 1.0/3.0 * fabs(v_an);
	}
	else
	{
		*v_an  = 0;
		*p_an  = p0;
		*ro_an = ro0;
	}
}

void CSolver::initVars(void) {	
	CFL = task.getCFL();
	tauPulse  = task.getTauPulse();
	fluence   = task.getFluence();
	deltaSkin = task.getDeltaSkin(); 
	if(task.getSourceFlag()==SourceType::SrcUndef) {
		x_pulse_min = 0.;
		i_pulse_min = 0; 
		//tInit = 0.;
	}	 	
	if(task.getSourceFlag() == SourceType::SrcGlass) {
		x_pulse_min = task.getZone(0).l;
		i_pulse_min = task.getZone(0).n - 1;
		//tInit = -5.0*tauPulse;
	}	
	if(task.getSourceFlag() == SourceType::SrcMetal) {
		x_pulse_min = 0.; // task.getZone(0).l;
		i_pulse_min = 0; //task.getZone(0).n - 1; 
		//tInit = 0.0;
	}
}

int CSolver::calcHydroStage(double t, double tau) {	
	int i=0, counter=0, itNumFull=0, itNumIon=0;
	double e_prev=0., ei_prev=0., Q=0., dv=0.,
		   p_next_plus = 0., p_next_minus = 0., p_plus = 0., p_minus = 0.;
	CFieldOld ms_temp, ms_prev;
    ms_temp.initData(&task);
	ms_prev.initData(&task);
	EOSOld &eos = task.getEOS();
	int nSize = ms.getSize();
	double h = ms[0].dm;
	int itCounter = 0;
	double *g = new double[nSize];
	for(i=0; i<nSize; i++) g[i]=0.;
	if(task.getViscFlag()) {
		for(int i=1; i<nSize-1; i++) {
			double dv = ms[i+1].v - ms[i].v;
			 Q =6000.0*ms[i].ro*dv*dv; // Al
			//Q = 360.0*ms[i].ro*dv*dv; // Ni
			//Q = 3000.0*ms[i].ro*dv*dv; // Ta
			if (dv < 0) {
				g[i] += Q;
			}
		}
	}
	for(i=0; i<=nSize; i++) {
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
	do {
		itCounter++;
		if(itCounter > 30)
		{
			cout << endl << 
				 "CalcHydroStage() error: no convergence in 30 iterations" << endl;
			dumpToFile(t);
			exit(1);
		}
		for(i=0; i<nSize; i++) {
			ms_prev[i].x  = ms_temp[i].x;
			ms_prev[i].v  = ms_temp[i].v;
			ms_prev[i].ro = ms_temp[i].ro;
			ms_prev[i].e  = ms_temp[i].e;
			ms_prev[i].ti = ms_temp[i].ti;
			ms_prev[i].p  = ms_temp[i].p;
			ms_prev[i].te = ms_temp[i].te;
			ms_prev[i].pi = ms_temp[i].pi;
			ms_prev[i].pe = ms_temp[i].pe;
			ms_prev[i].p  = ms_temp[i].p;
		}	
		ms_prev[nSize].x  = ms_temp[nSize].x;
		ms_prev[nSize].v  = ms_temp[nSize].v;
		for(i=0; i<=nSize; i++) {
			if(i==0) {
				p_next_plus  = ms_prev[i].p;
			    p_next_minus = 0.0;
			    p_plus       = ms[i].p;
			    p_minus      = 0.0;
			} else if(i==nSize) {
				p_next_plus  = 0.0;
			    p_next_minus = ms_prev[i-1].p; 
			    p_plus       = 0.0;
			    p_minus      = ms[i-1].p;  
			} else {
				p_next_plus  = ms_prev[i].p + g[i];
				p_next_minus = ms_prev[i-1].p + g[i-1];
				p_plus       = ms[i].p + g[i];
			    p_minus      = ms[i-1].p + g[i-1];
			}
			ms_temp[i].v = ms[i].v - tau/2.0/h*(p_next_plus - p_next_minus + p_plus - p_minus);
			ms_temp[i].x = ms[i].x + 0.5*tau*(ms_temp[i].v + ms[i].v);
		}		

		for(int i=0; i<nSize; i++) {
			ms_temp[i].ro = 1.0/(1.0/ms[i].ro + tau/2.0/h*
				            (ms_temp[i+1].v + ms[i+1].v - ms_temp[i].v-ms[i].v));			
			ms_temp[i].ei = ms[i].ei - tau/4.0/h*(ms_prev[i].pi + ms[i].pi + g[i]) *
							(ms_temp[i+1].v + ms[i+1].v - ms_temp[i].v - ms[i].v);
			ms_temp[i].ee = ms[i].ee - tau/4.0/h*(ms_prev[i].pe + ms[i].pe + g[i]) *
							(ms_temp[i+1].v + ms[i+1].v - ms_temp[i].v - ms[i].v);
			ms_temp[i].e  = ms_temp[i].ei + ms_temp[i].ee;		
		}
		for(int i=0; i<nSize; i++)	{				
			ms_temp[i].ti = eos.getti(ms_temp[i].ro, ms_temp[i].ei);
			if(ms_temp[i].ti == -1.) {
				cout << "Ooops! No solution of non-linear equation of hydrodynamics in node " << i << endl 
					 << "ro = " << ms[i].ro << endl 
					 << "ti = " << ms[i].ti << endl
					 << "te = " << ms[i].te << endl
					 << "e = "  << ms[i].e  << endl
					 << "ei = " << ms[i].ei << endl
					 << "ee = " << ms[i].ee << endl
					 << "p = "  << ms[i].p  << endl
					 << "pi = " << ms[i].pi << endl
					 << "pe = " << ms[i].pi << endl;
				saveSolution("error.dat", t);
				cin.get();
			}
			ms_temp[i].te = eos.gette(ms_temp[i].ro, ms_temp[i].ti, ms_temp[i].ee);			
		}
		for(int i=0; i<nSize; i++) 	{			
			ms_temp[i].pi = eos.getpi(ms_temp[i].ro, ms_temp[i].ti);
			ms_temp[i].pe = eos.getpe(ms_temp[i].ro, ms_temp[i].ti, ms_temp[i].te);
			ms_temp[i].p  = ms_temp[i].pi + ms_temp[i].pe;
		}		
	} 	
	while (compTi(ms_temp, ms_prev) > .01);	
	for(int i=0; i<nSize; i++) {
		 ms[i].x = ms_temp[i].x;
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
		ms[i].ce    = eos.getce(ms[i].ro, ms[i].te);
	}
	ms[nSize].x = ms_temp[nSize].x;
	ms[nSize].v = ms_temp[nSize].v;	
	delete[] g;
	return 1;
}

int CSolver::calcHydroStageGlass(double t, double tau) {	
	int i=0, counter=0, itNumFull=0, itNumIon=0;
	double h = ms[0].dm, e_prev=0., ei_prev=0., Q=0, dv=0;
	double p_next_plus  = 0., p_next_minus = 0., p_plus = 0., p_minus = 0.;
	CFieldOld ms_temp, ms_prev;
    ms_temp.initData(&task);
	ms_prev.initData(&task);
	EOSOld &eos = task.getEOS();
	EOSOld &eosGlass = task.getEOSGlass();
	SourceType flag = task.getSourceFlag();
	int nSize = ms.getSize(), nBound = task.getZone(0).n, itCounter = 0;
	double *g = new double[nSize];
	for(i=0; i<nSize; i++) g[i]=0.;
	if(task.getViscFlag()) {
		for(i=1; i<nSize-1; i++) {
			double dv = ms[i+1].v - ms[i].v;
			Q = 8000.*ms[i].ro*dv*dv; // Al
			if (dv < 0) {
				g[i] += Q;
			}
		}
	}
	for(i=0; i<=nSize; i++) {
		ms_temp[i].x  = ms[i].x;
		ms_temp[i].v  = ms[i].v;
		ms_temp[i].ro = ms[i].ro;
		ms_temp[i].e  = ms[i].e;
		ms_temp[i].ti = ms[i].ti;
		ms_temp[i].te = ms[i].te;
		ms_temp[i].pi = ms[i].pi;
		ms_temp[i].pe = ms[i].pe;
		ms_temp[i].p  = ms[i].p;
	}
	i=nSize;
	ms_temp[i].x  = ms[i].x;
	ms_temp[i].v  = ms[i].v;
	do {
		itCounter++;	
		if(itCounter > 30)	{
			cout << endl << 
				 "CalcHydroStageGlass() error: no convergence in 30 iterations" << endl;
			dumpToFile(t);
			exit(1);
		}
		for(i=0; i<nSize; i++) { 
			ms_prev[i].x  = ms_temp[i].x;
			ms_prev[i].v  = ms_temp[i].v;
			ms_prev[i].ro = ms_temp[i].ro;
			ms_prev[i].e  = ms_temp[i].e;
			ms_prev[i].ti = ms_temp[i].ti;
			ms_prev[i].p  = ms_temp[i].p;
			ms_prev[i].te = ms_temp[i].te;
			ms_prev[i].pi = ms_temp[i].pi;
			ms_prev[i].pe = ms_temp[i].pe;
			ms_prev[i].p  = ms_temp[i].p;
		}	
		ms_prev[nSize].x  = ms_temp[nSize].x;
		ms_prev[nSize].v  = ms_temp[nSize].v;
		for(i=0; i<=nSize; i++) {
			if(i==0) {
				p_next_plus  = ms_prev[i].p;
			    p_next_minus = 0.0;
			    p_plus       = ms[i].p;
			    p_minus      = 0.0;
			} else if(i==nSize) {
				p_next_plus  = 0.0;
			    p_next_minus = ms_prev[i-1].p; 
			    p_plus       = 0.0;
			    p_minus      = ms[i-1].p;  
			} else {
				p_next_plus  = ms_prev[i].p + g[i];
				p_next_minus = ms_prev[i-1].p + g[i-1];
				p_plus       = ms[i].p + g[i];
			    p_minus      = ms[i-1].p + g[i-1];
			}
			ms_temp[i].v = ms[i].v - tau/2.0/h*(p_next_plus - p_next_minus + p_plus - p_minus);
			ms_temp[i].x = ms[i].x + 0.5*tau*(ms_temp[i].v + ms[i].v);
		}
		for(i=0; i<nSize; i++) {
			ms_temp[i].ro = 1.0/(1.0/ms[i].ro + tau/2.0/h*
				            (ms_temp[i+1].v + ms[i+1].v - ms_temp[i].v-ms[i].v));
			ms_temp[i].ei = ms[i].ei - tau/4.0/h*(ms_prev[i].pi + ms[i].pi + g[i]) *
							(ms_temp[i+1].v + ms[i+1].v - ms_temp[i].v - ms[i].v);
			ms_temp[i].ee = ms[i].ee - tau/4.0/h*(ms_prev[i].pe + ms[i].pe + g[i]) *
							(ms_temp[i+1].v + ms[i+1].v - ms_temp[i].v - ms[i].v);
			ms_temp[i].e  = ms_temp[i].ei + ms_temp[i].ee;		
		}
		for(i=0; i<nSize; i++)	{	
			if(flag == SourceType::SrcGlass) {
				if(i >= nBound) 
				    ms_temp[i].ti = eos.getti(ms_temp[i].ro, ms_temp[i].ei);
				else
				    ms_temp[i].ti = eosGlass.getti(ms_temp[i].ro, ms_temp[i].ei);
			} else if(i<nBound) 
				ms_temp[i].ti = eos.getti(ms_temp[i].ro, ms_temp[i].ei);
			else
				ms_temp[i].ti = eosGlass.getti(ms_temp[i].ro, ms_temp[i].ei);
			if(ms_temp[i].ti == -1.) {
				cout << "CSolver::calcHydroStageGlass() error: no ti equation root" << endl << "on the interval in node " << i << ":" << endl
					 << "i = " << i << endl
					 << "x = " << ms[i].x << endl 
					 << "ro = " << ms[i].ro << endl 
					 << "ti = " << ms[i].ti << endl
					 << "te = " << ms[i].te << endl
					 << "e = "  << ms[i].e  << endl
					 << "ei = " << ms[i].ei << endl
					 << "ee = " << ms[i].ee << endl
					 << "p = "  << ms[i].p  << endl
					 << "pi = " << ms[i].pi << endl
					 << "pe = " << ms[i].pe << endl;
				saveSolution("error-ti.dat", t);
				cin.get();
				exit(1);
			}
			if(flag == SourceType::SrcGlass) {
				if(i >= nBound) 
					ms_temp[i].te = eos.gette(ms_temp[i].ro, ms_temp[i].ti, ms_temp[i].ee);			
				else
					ms_temp[i].te = eosGlass.gette(ms_temp[i].ro, ms_temp[i].ti, ms_temp[i].ee);
			} else if(i<nBound) 
				ms_temp[i].te = eos.gette(ms_temp[i].ro, ms_temp[i].ti, ms_temp[i].ee);			
			else
				ms_temp[i].te = eosGlass.gette(ms_temp[i].ro, ms_temp[i].ti, ms_temp[i].ee);
			if(ms_temp[i].te == -1.) {
				cout << "CSolver::calcHydroStageGlass() error: no te equation root on the interval in node " << i << "." << endl 
					 << "ro = " << ms[i].ro << endl 
					 << "ti = " << ms[i].ti << endl
					 << "te = " << ms[i].te << endl
					 << "e = "  << ms[i].e  << endl
					 << "ei = " << ms[i].ei << endl
					 << "ee = " << ms[i].ee << endl
					 << "p = "  << ms[i].p  << endl
					 << "pi = " << ms[i].pi << endl
					 << "pe = " << ms[i].pe << endl;
				saveSolution("error-te.dat", t);
				cin.get();
				exit(1);
			}
		}
		for(i=0; i<nSize; i++) 	{			
			if(flag == SourceType::SrcGlass) {
				if(i >= nBound) {
					ms_temp[i].pi = eos.getpi(ms_temp[i].ro, ms_temp[i].ti);
					ms_temp[i].pe = eos.getpe(ms_temp[i].ro, ms_temp[i].ti, ms_temp[i].te);
				} else {
					ms_temp[i].pi = eosGlass.getpi(ms_temp[i].ro, ms_temp[i].ti);
					ms_temp[i].pe = eosGlass.getpe(ms_temp[i].ro, ms_temp[i].ti, ms_temp[i].te);
				}
			} else if(i<nBound) {
				ms_temp[i].pi = eos.getpi(ms_temp[i].ro, ms_temp[i].ti);
				ms_temp[i].pe = eos.getpe(ms_temp[i].ro, ms_temp[i].ti, ms_temp[i].te);
			} else {
				ms_temp[i].pi = eosGlass.getpi(ms_temp[i].ro, ms_temp[i].ti);
				ms_temp[i].pe = eosGlass.getpe(ms_temp[i].ro, ms_temp[i].ti, ms_temp[i].te);
			}
			ms_temp[i].p  = ms_temp[i].pi + ms_temp[i].pe;
		}		
	} while (compTi(ms_temp, ms_prev) > 0.01);
	//cout << itCounter << "it ";
	for(i=0; i<nSize; i++) {
		 ms[i].x = ms_temp[i].x;
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
		if(flag == SourceType::SrcGlass) {
				if(i >= nBound) {
					ms[i].C     = eos.getC(ms[i].ro, ms[i].ti, ms[i].te);
					ms[i].Alphaei	    = eos.getAlpha(ms[i].ro, ms[i].ti, ms[i].te);
					ms[i].kappa = eos.getkappa(ms[i].ro, ms[i].ti, ms[i].te);
					ms[i].ce    = eos.getce(ms[i].ro, ms[i].te);
				} else {
					ms[i].C     = eos.getC(ms[i].ro, ms[i].ti, ms[i].te); // 'eos', а не 'eosGlass' сознательно, чтобы не уменьшить шаг по времени случайно
					ms[i].Alphaei	    = eosGlass.getAlpha(ms[i].ro, ms[i].ti, ms[i].te);
					ms[i].kappa = eosGlass.getkappa(ms[i].ro, ms[i].ti, ms[i].te);
					ms[i].ce    = eosGlass.getce(ms[i].ro, ms[i].te);
				}
		} else if(i<nBound) {
			ms[i].C     = eos.getC(ms[i].ro, ms[i].ti, ms[i].te);
			ms[i].Alphaei	    = eos.getAlpha(ms[i].ro, ms[i].ti, ms[i].te);
			ms[i].kappa = eos.getkappa(ms[i].ro, ms[i].ti, ms[i].te);
			ms[i].ce    = eos.getce(ms[i].ro, ms[i].te);
		} else {
			ms[i].C     = eos.getC(ms[i].ro, ms[i].ti, ms[i].te); // 'eos', а не 'eosGlass' сознательно, чтобы не уменьшить шаг по времени случайно
			ms[i].Alphaei	    = eosGlass.getAlpha(ms[i].ro, ms[i].ti, ms[i].te);
			ms[i].kappa = eosGlass.getkappa(ms[i].ro, ms[i].ti, ms[i].te);
			ms[i].ce    = eosGlass.getce(ms[i].ro, ms[i].te);
		}
	}
	ms[nSize].x = ms_temp[nSize].x;
	ms[nSize].v = ms_temp[nSize].v;	
	delete[] g;
	return itCounter;
}



void CSolver::calcHydroStageSpallation(double t, double tau)
{	
	int i=0; int counter=0; int itCounter = 0; int iSpall = getSpallCellNum();
	double Q=0., dv=0.;
	double p_next_plus = 0.0, p_next_minus = 0.0, p_plus = 0.0, p_minus = 0.0;
	CFieldOld ms_temp, ms_prev;
    ms_temp.initData(&task); ms_prev.initData(&task);
	EOSOld &eos = task.getEOS(); int nSize = ms.getSize(); double h = ms[0].dm;
	///// Artificial viscosity
	double *g = new double[nSize]; for(int i=0; i<nSize; i++) { g[i] = 0.;}
	if(task.getViscFlag()) {
		for(int i=1; i<nSize-1; i++) {
			double dv = ms[i+1].v - ms[i].v;
			Q =160.0*ms[i].ro*dv*dv; 
			if (dv < 0)	
				g[i] = Q;
		}
	}

	for(i=0; i<nSize; i++)
	{
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


	do	{
		itCounter++;
		if(itCounter > 30)
		{
			cout << endl << 
				 "CalcHydroStageSpallation() error: no convergence in 30 iterations" << endl;
			dumpToFile(t);
			exit(1);
		}
		for(int i=0; i<nSize; i++)
		{
			ms_prev[i].x  = ms_temp[i].x;
			ms_prev[i].v  = ms_temp[i].v;
			ms_prev[i].ro = ms_temp[i].ro;
			ms_prev[i].e  = ms_temp[i].e;
			ms_prev[i].ti = ms_temp[i].ti;
			ms_prev[i].p  = ms_temp[i].p;
			ms_prev[i].te = ms_temp[i].te;
			ms_prev[i].pi = ms_temp[i].pi;
			ms_prev[i].pe = ms_temp[i].pe;
			ms_prev[i].p  = ms_temp[i].p;
		}	
		int j = nSize;
		ms_prev[j].x  = ms_temp[j].x;
		ms_prev[j].v  = ms_temp[j].v;

		for(int i=0; i<=nSize; i++)	{
			if((i==0)||(i==iSpall+1)) {
				p_next_plus  = ms_prev[i].p;
			    p_next_minus = 0.0;
			    p_plus       = ms[i].p;
			    p_minus      = 0.0;
				ms_temp[i].v = ms[i].v - tau/4.0/h*(p_next_plus - p_next_minus + p_plus - p_minus);
			}
			else if((i==nSize)||(i==iSpall)) {
				p_next_plus  = 0.0; 
			    p_next_minus = ms_prev[i-1].p; 
			    p_plus       = 0.0;
			    p_minus      = ms[i-1].p;  
				ms_temp[i].v = ms[i].v - tau/4.0/h*(p_next_plus - p_next_minus + p_plus - p_minus);
			}
			else {
				p_next_plus  = ms_prev[i].p + g[i];
				p_next_minus = ms_prev[i-1].p + g[i-1];
				p_plus       = ms[i].p + g[i];
			    p_minus      = ms[i-1].p + g[i-1];
				ms_temp[i].v = ms[i].v - tau/2.0/h*(p_next_plus - p_next_minus + p_plus - p_minus);
			}
			ms_temp[i].x = ms[i].x + 0.5*tau*(ms_temp[i].v + ms[i].v);
		}
		for(int i=0; i<nSize; i++) {
			if(i!=iSpall) {	
				ms_temp[i].ro = 1.0/(1.0/ms[i].ro + tau/2.0/h*
					            (ms_temp[i+1].v + ms[i+1].v - ms_temp[i].v-ms[i].v));
				ms_temp[i].ei = ms[i].ei - tau/4.0/h*(ms_prev[i].pi + ms[i].pi + g[i]) *
								(ms_temp[i+1].v + ms[i+1].v - ms_temp[i].v - ms[i].v);
				ms_temp[i].ee = ms[i].ee - tau/4.0/h*(ms_prev[i].pe + ms[i].pe + g[i]) *
								(ms_temp[i+1].v + ms[i+1].v - ms_temp[i].v - ms[i].v);
				ms_temp[i].e  = ms_temp[i].ei + ms_temp[i].ee;
			} else {
				ms_temp[i].ro = eos.getMIN_RO();
				ms_temp[i].ti = eos.getMIN_T();
				ms_temp[i].te = eos.getMIN_T();
				ms_temp[i].ei = eos.getei(ms_temp[i].ro, ms_temp[i].ti);
				ms_temp[i].ee = eos.getee(ms_temp[i].ro, ms_temp[i].ti, ms_temp[i].te);
				ms_temp[i].e  = ms_temp[i].ei + ms_temp[i].ee;
			}
		}

		for(int i=0; i<nSize; i++) 	{
			if(i != iSpall)	{
				ms_temp[i].ti = eos.getti(ms_temp[i].ro, ms_temp[i].ei);
				if(ms_temp[i].ti == -1.) {
					cout << "Ooops! No solution of non-linear equation of hydrodynamics in node " << i << endl 
						 << "ro = " << ms[i].ro << endl 
						 << "ti = " << ms[i].ti << endl
						 << "te = " << ms[i].te << endl
						 << "e = "  << ms[i].e  << endl
						 << "ei = " << ms[i].ei << endl
						 << "ee = " << ms[i].ee << endl
						 << "p = "  << ms[i].p  << endl
						 << "pi = " << ms[i].pi << endl
						 << "pe = " << ms[i].pi << endl;
					saveSolution("error.dat", t);
					cin.get();
				}
				ms_temp[i].te = eos.gette(ms_temp[i].ro, ms_temp[i].ti, ms_temp[i].ee);
			} 
		}
		for(int i=0; i<nSize; i++) {
			if(i!=iSpall) {	
				ms_temp[i].pi = eos.getpi(ms_temp[i].ro, ms_temp[i].ti);
				ms_temp[i].pe = eos.getpe(ms_temp[i].ro, ms_temp[i].ti, ms_temp[i].te);
			} else {
				ms_temp[i].pi = eos.getpi(ms_temp[i].ro, ms_temp[i].ti);
				ms_temp[i].pe = eos.getpe(ms_temp[i].ro, ms_temp[i].ti, ms_temp[i].te);
			}
			ms_temp[i].p  = ms_temp[i].pi + ms_temp[i].pe;
		}
	} while (compTi(ms_temp, ms_prev) > 0.01);

	printf("calcHydroStageSpallation() complete: %d iterations\n", itCounter);
	for(int i=0; i<nSize; i++) 	{
		ms[i].x  = ms_temp[i].x;
		ms[i].v  = ms_temp[i].v;
	 	ms[i].ro = ms_temp[i].ro;
		ms[i].e  = ms_temp[i].e;
		ms[i].ee = ms_temp[i].ee;
		ms[i].ei = ms_temp[i].ei;
		ms[i].p  = ms_temp[i].p;
		ms[i].pi = ms_temp[i].pi;
		ms[i].pe = ms_temp[i].pe;
		ms[i].te = ms_temp[i].te;
		ms[i].ti = ms_temp[i].ti;
		if(i!=iSpall) {
			ms[i].C     = eos.getC(ms[i].ro, ms[i].ti, ms[i].te);
			ms[i].Alphaei	    = eos.getAlpha(ms[i].ro, ms[i].ti, ms[i].te);
			ms[i].kappa = eos.getkappa(ms[i].ro, ms[i].ti, ms[i].te);
			ms[i].ce    = eos.getce(ms[i].ro, ms[i].te);
		} else {
			ms[i].C     = eos.getC(eos.getMIN_RO(), ms[i].ti, ms[i].te);
			ms[i].Alphaei	    = eos.getAlpha(eos.getMIN_RO(), ms[i].ti, ms[i].te);
			ms[i].kappa = eos.getkappa(eos.getMIN_RO(), ms[i].ti, ms[i].te);
			ms[i].ce    = eos.getce(eos.getMIN_RO(), ms[i].te);
		}
	}
	ms[nSize].x = ms_temp[nSize].x;
	ms[nSize].v = ms_temp[nSize].v;	
	delete[] g;
}

void CSolver::calcHydroStageNoElectron(double t, double tau) {	
	int i=0, counter=0, itNumFull=0, itNumIon=0;	
	double e_prev = 0., ei_prev = 0., Alphaei = 0., dv = 0.,
		   p_next_plus  = 0., p_next_minus = 0., p_plus = 0., p_minus = 0.;
	CFieldOld ms_temp, ms_prev;
    ms_temp.initData(&task);
	ms_prev.initData(&task);
	EOSOld &eos = task.getEOS();
	int nSize = ms.getSize();
	double h = ms[0].dm;
	int itCounter = 0;
	double *g = new double[nSize];
	for(i=0; i<nSize; i++) {
		g[i]=0.0;
	}
	
	if(task.getViscFlag())
	{
		for(i=1; i<nSize-1; i++)
		{
			double dv = ms[i+1].v - ms[i].v;
			// Alphaei=160.0*ms[i].ro*dv*dv; // Al
			
			 Alphaei=560.0*ms[i].ro*dv*dv; // Ni
			
			if (dv < 0)
			{
				g[i]+= Alphaei;
			}
		}
	}

	for(i=0; i<nSize; i++)
	{
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



	do
	{
		itCounter++;

		if(itCounter > 30)
		{
			cout << endl << 
				 "CalcHydroStage() error: no convergence in 30 iterations" << endl;
			dumpToFile(t);
			exit(1);
		}

		Node &n_temp = ms_temp[0];
		Node &n_prev = ms_prev[0];
		Node &n      = ms[0];


		for(i=0; i<nSize; i++)
		{
			ms_prev[i].x  = ms_temp[i].x;
			ms_prev[i].v  = ms_temp[i].v;
			ms_prev[i].ro = ms_temp[i].ro;
			ms_prev[i].e  = ms_temp[i].e;
			ms_prev[i].ti = ms_temp[i].ti;
			ms_prev[i].p  = ms_temp[i].p;
			ms_prev[i].te = ms_temp[i].te;
			ms_prev[i].pi = ms_temp[i].pi;
			ms_prev[i].pe = ms_temp[i].pe;
			ms_prev[i].p  = ms_temp[i].p;
		}
	
		int j = nSize;
		ms_temp[j].x  = ms[j].x;
		ms_temp[j].v  = ms[j].v;

	// О, ахтунг! В верхних трех строчках может крыться небольшая ошибка, влияющая
	// на точность.
	// Везде мы приравниваем мс_прев мс_темпу, а в конечной точке мс_темп мс'у.
	// Т.е., мс_прев в конечной точке отпущен на волю случая.


		for(i=0; i<=nSize; i++)
		{
			if(i==0)
			{
				
				// Эта аппроксимация граничных условий на давление, быть может, 
				// чуть менее точна на автомодельной ВР в идеальном газе, но гораздо 
				// монотоннее на алюминии (не порождает жутких осцилляций на границе).

				p_next_plus  = ms_prev[i].p;
			    p_next_minus = 0.0; //-ms_prev[i].p; 
			    p_plus       = ms[i].p;
			    p_minus      = 0.0; //-ms[i].p;  

				ms_temp[i].v = ms[i].v - tau/4.0/h*
				          (p_next_plus - p_next_minus + p_plus - p_minus);

			}

			else if(i==nSize)
			{
				p_next_plus  = 0.0; //-ms_prev[i-1].p;
			    p_next_minus = ms_prev[i-1].p; 
			    p_plus       = 0.0; //-ms[i-1].p;
			    p_minus      = ms[i-1].p;  
				
				ms_temp[i].v = ms[i].v - tau/4.0/h*
				          (p_next_plus - p_next_minus + p_plus - p_minus);
			}
			
			else
			{
				p_next_plus  = ms_prev[i].p + g[i];
				p_next_minus = ms_prev[i-1].p + g[i-1];
				p_plus       = ms[i].p + g[i];
			    p_minus      = ms[i-1].p + g[i-1];

				ms_temp[i].v = ms[i].v - tau/2.0/h*
				          (p_next_plus - p_next_minus + p_plus - p_minus);

			}

			ms_temp[i].x = ms[i].x + 0.5*tau*(ms_temp[i].v + ms[i].v);

		}
	
/*		i=nSize;

		ms_temp[i].v = 0.0;
		ms_temp[i].x = ms[i].x + 0.5*tau*(ms_temp[i].v + ms[i].v);
*/		

		for(i=0; i<nSize; i++)
		{
			


			ms_temp[i].ro = 1.0/(1.0/ms[i].ro + tau/2.0/h*
				            (ms_temp[i+1].v + ms[i+1].v - ms_temp[i].v-ms[i].v));

		
			ms_temp[i].ei = ms[i].ei - tau/4.0/h*(ms_prev[i].pi + ms[i].pi + g[i]) *
							(ms_temp[i+1].v + ms[i+1].v - ms_temp[i].v - ms[i].v);

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

			src=fluence/sqrt(3.14159)/deltaSkin/tauPulse/ms[i].ro*expX*expT;	

/*		// Здесь о том, что какая-то часть импульса все же поглащается 
		// в стекле у границы "стекло-алюминий". Спросить у Н.А. 

			if ((task.getSourceFlag() == 2) && (i<=i_pulse_min))
			{
				double F_glass = 1000.0;
				double d_glass = 100.0e-9; 

				expX = exp(- ((ms[i].x-x_pulse_min) * (ms[i].x-x_pulse_min)) 
				               /d_glass/d_glass);
				expT = exp(-(t)*(t)/tauPulse/tauPulse);

				src=F_glass/sqrt(3.14159)/d_glass/tauPulse/ms[i].ro*expX*expT;	
			}
*/

		/*	ms_temp[i].ee = ms[i].ee - tau/4.0/h*(ms_prev[i].pe + ms[i].pe + g[i]) *
							(ms_temp[i+1].v + ms[i+1].v - ms_temp[i].v - ms[i].v) + tau*src;*/

			ms_temp[i].ee = eos.getee(ms_temp[i].ro, ms_temp[i].ti, ms_temp[i].te);
			ms_temp[i].e  = ms_temp[i].ei + ms_temp[i].ee;

		
		}

		/////////////////////
		double _te = ms_temp[0].te;
		double _ti = ms_temp[0].ti;
		/////////////////////

		for(i=0; i<nSize; i++)
		{

			ms_temp[i].ti = eos.getti(ms_temp[i].ro, ms_temp[i].ei);
			// ms_temp[i].te = eos.gette(ms_temp[i].ro, ms_temp[i].ti, ms_temp[i].ee, ms_temp[i].Z);
			
		}

		for(i=0; i<nSize; i++)
		{
			
			ms_temp[i].pi = eos.getpi(ms_temp[i].ro, ms_temp[i].ti);
			ms_temp[i].pe = eos.getpe(ms_temp[i].ro, ms_temp[i].ti, ms_temp[i].te);
			ms_temp[i].p  = ms_temp[i].pi + ms_temp[i].pe;
		}
		
	} 	
	while (compTi(ms_temp, ms_prev) > 0.01); 
	printf("calcHydroStageNoElectron() complete: %d iterations\n",
		    itCounter);
	for(i=0; i<nSize; i++) 	{
		ms_temp[i].pi = eos.getpi(ms_temp[i].ro, ms_temp[i].ti);
		ms_temp[i].pe = eos.getpe(ms_temp[i].ro, ms_temp[i].ti, ms_temp[i].te);
		ms_temp[i].p  = ms_temp[i].pi + ms_temp[i].pe;
	}
	for(i=0; i<nSize; i++) 	{
		 ms[i].x = ms_temp[i].x;
		 ms[i].v = ms_temp[i].v;
		ms[i].ro = ms_temp[i].ro;
		 ms[i].e = ms_temp[i].e;
		ms[i].ee = ms_temp[i].ee;
		ms[i].ei = ms_temp[i].ei;
		ms[i].p  = ms_temp[i].p;
		ms[i].pi = ms_temp[i].pi;
		ms[i].pe = ms_temp[i].pe;
//		ms[i].te = ms_temp[i].te;
		ms[i].ti = ms_temp[i].ti;
		ms[i].C     = eos.getC(ms[i].ro, ms[i].ti, ms[i].te);
		ms[i].Alphaei	    = eos.getAlpha(ms[i].ro, ms[i].ti, ms[i].te);
		ms[i].kappa = eos.getkappa(ms[i].ro, ms[i].ti, ms[i].te);
		ms[i].ce    = eos.getce(ms[i].ro, ms[i].te);
	}
	ms[nSize].x = ms_temp[nSize].x;
	ms[nSize].v = ms_temp[nSize].v;	
	delete[] g;
}

double CSolver::solveti(double tau, int i, double ro_temp, double dv_temp) {
	double	ei_temp  = 0.0, 
			pi_temp  = 0.0,
			ti_temp  = 0.0,
			criteria = 0.0;
			
	double 	ei		 = ms[i].ei,
			dm		 = ms[i].dm;

	EOSOld& eos = task.getEOS();

	double  lowBorder = eos.getMIN_T() + 1.0,
		   highBorder = eos.getMAX_T() - 1.0;

	double lValue = eos.getei(ro_temp, lowBorder)  - ei + 
		            tau/dm*dv_temp*eos.getpi(ro_temp, lowBorder);

	double rValue = eos.getei(ro_temp, highBorder) - ei + 
		            tau/dm*dv_temp*eos.getpi(ro_temp, highBorder);

	if (( lValue > 0) || ( rValue < 0) )
	{
		cout << "Error solveei()No equation root on interval: i=" << i << endl;
		dumpToFile(0.0);
		exit(1);
	}
	
	do
	{
		ti_temp  = 0.5*(lowBorder + highBorder);
		ei_temp  = eos.getei(ro_temp, ti_temp);
		pi_temp  = eos.getpi(ro_temp, ti_temp);
		
		criteria = ei_temp - ei + tau/dm*dv_temp*pi_temp;

		if (criteria > 0.0)
		{
			highBorder = ti_temp;
		}
		else if (criteria < 0.0)
		{
			lowBorder  = ti_temp;
		}
		else
		{
			return ti_temp;
		}

		if (fabs(highBorder-lowBorder) < 1.e-8)
		{
			cout << "solveei(): Possibly incorrect convergence" << endl <<
				    "i="    << i       << ", criteria=" << criteria <<
					", ti=" << ti_temp << endl;		
			exit(1);
		}

	} while( fabs(criteria) > epsE );

	return ti_temp;
}

double CSolver::solvete(double t, double tau, int i, double ro_temp, double dv_temp, double ti_temp)
{
	double	ee_temp  = 0.0, 
			pe_temp  = 0.0,
			te_temp  = 0.0,
			criteria = 0.0;
			
	double 	ee		 = ms[i].ee,
			dm		 = ms[i].dm;

	double src = 0.0;

	EOSOld& eos = task.getEOS();
	
	if(eos.getType() == ideal )
		return 0.0;

	double  lowBorder = eos.getMIN_T() + 1.0,
		   highBorder = eos.getMAX_T() - 1.0;
	
	double lVal = eos.getee(ro_temp, ti_temp, lowBorder) - ee + 
			 	  tau/dm*dv_temp*eos.getpe(ro_temp, ti_temp,  lowBorder) 
				  - tau*src;
	double rVal = eos.getee(ro_temp, ti_temp, highBorder) - ee + 
			 	  tau/dm*dv_temp*eos.getpe(ro_temp, ti_temp,  highBorder) 
				  - tau*src;

	if ((lVal > 0) || (rVal < 0) )
	{
		cout << "Error solveee():No equation root on interval: i=" << i << endl;
		exit(1);
	}

	do
	{
		te_temp  = 0.5*(lowBorder + highBorder);
		ee_temp  = eos.getee(ro_temp, ti_temp, te_temp);
		pe_temp  = eos.getpe(ro_temp, ti_temp, te_temp);
		
		criteria = ee_temp - ee + tau/dm*dv_temp*pe_temp - tau*src;
		if (criteria > 0.0)
		{
			highBorder = te_temp;
		}
		else if (criteria < 0.0)
		{
			lowBorder  = te_temp;
		}
		else
		{
			return te_temp;
		}

	} while( fabs(criteria) > epsE );
	
	double val = ee_temp - ms[i].ee + tau/ms[i].dm*dv_temp*pe_temp - tau*src;

	return te_temp;
}

double CSolver::solveteConservative(double t, double tau, int i, double ro_temp, double dv_temp, double ti_temp)
{
	double	ee_temp  = 0.0, 
			pe_temp  = 0.0,
			te_temp  = 0.0,
			criteria = 0.0;
			
	double 	ee		 = ms[i].ee,
			dm		 = ms[i].dm;

	double src = 0.0;

	EOSOld& eos = task.getEOS();
	
	if(eos.getType() == ideal )
		return 0.0;

	double  lowBorder = eos.getMIN_T() + 1.0,
		   highBorder = eos.getMAX_T() - 1.0;
	
/////////////////////////DEBUG!!!////////////////////
	double tauLaser = 100.0e-15;

	if(task.getSourceFlag()==1)
	{

		double expX = exp(-(ms[i].x)/10.0e-9);
		double expT = exp(-(t-3.0*tauLaser)*(t-3.0*tauLaser)/tauLaser/tauLaser);
		src = 650.0/sqrt(3.14159)/10.0e-9/tauLaser/ro_temp * expX * expT;
	}
	else if (task.getSourceFlag()==2)
	{
		if(i>=600)
		{
			double expX = exp(-(ms[i].x-600.0e-9)/10.0e-9);
			double expT = exp(-(t-3.0*tauLaser)*(t-3.0*tauLaser)/tauLaser/tauLaser);
			src = 650.0/sqrt(3.14159)/10.0e-9/tauLaser/ro_temp * expX * expT;
		}
		else
			src =0.0;
	}
	else
		src = 0.0;
/////////////////////////////////////////////////////
	
	double lVal = eos.getee(ro_temp, ti_temp, lowBorder) - ee + 
			 	  tau/dm*dv_temp*eos.getpe(ro_temp, ti_temp,  lowBorder) 
				  - tau*src;
	double rVal = eos.getee(ro_temp, ti_temp, highBorder) - ee + 
			 	  tau/dm*dv_temp*eos.getpe(ro_temp, ti_temp,  highBorder) 
				  - tau*src;

	if ((lVal > 0) || (rVal < 0) )
	{
		cout << "Error solveee():No equation root on interval: i=" << i << endl;
		exit(1);
	}

	do
	{
		te_temp  = 0.5*(lowBorder + highBorder);
		ee_temp  = eos.getee(ro_temp, ti_temp, te_temp);
		pe_temp  = eos.getpe(ro_temp, ti_temp, te_temp);
		
		criteria = ee_temp - ee + tau/dm*dv_temp*pe_temp - tau*src;
		if (criteria > 0.0)
		{
			highBorder = te_temp;
		}
		else if (criteria < 0.0)
		{
			lowBorder  = te_temp;
		}
		else
		{
			return te_temp;
		}

	} while( fabs(criteria) > epsE );
	
	double val = ee_temp - ms[i].ee + tau/ms[i].dm*dv_temp*pe_temp - tau*src;

	return te_temp;
}

double CSolver::solveteSource(double t, double tau, int i, double ro_temp, double dv_temp, double ti_temp)
{
	double	ee_temp  = 0.0, 
			pe_temp  = 0.0,
			te_temp  = 0.0,
			criteria = 0.0,
			src		 = 0.0;

	double  tauLaser = 100.0e-15;
			
	double 	ee		 = ms[i].ee,
			dm		 = ms[i].dm;

	EOSOld& eos = task.getEOS();
	
	if(eos.getType() == ideal )
		return 0.0;

	double  lowBorder = eos.getMIN_T() + 1.0,
		   highBorder = eos.getMAX_T() - 1.0;
/////////////////////////DEBUG!!!////////////////////
	double expX = exp(-(ms[i].x-ms[600].x)/10.0e-9);
	double expT = exp(-t*t/tauLaser/tauLaser);
	
	if(i<600)
		src = 0.0;
	else
		src = 2360.0/sqrt(3.14159)/10.0e-9/tauLaser/ro_temp * expX * expT;
/////////////////////////////////////////////////////
	
	double lVal = eos.getee(ro_temp, ti_temp, lowBorder) - ee + 
			 	  tau/dm*dv_temp*eos.getpe(ro_temp, ti_temp,  lowBorder) 
				  - tau*src;
	double rVal = eos.getee(ro_temp, ti_temp, highBorder) - ee + 
			 	  tau/dm*dv_temp*eos.getpe(ro_temp, ti_temp,  highBorder) 
				  - tau*src;

	if ((lVal > 0) || (rVal < 0) )
	{
		cout << "Error solveee():No equation root on interval: i=" << i << endl;
		exit(1);
	}

	
	do
	{
		te_temp  = 0.5*(lowBorder + highBorder);
		ee_temp  = eos.getee(ro_temp, ti_temp, te_temp);
		pe_temp  = eos.getpe(ro_temp, ti_temp, te_temp);
		
		double deedte = eos.getdeedte_ro(ro_temp, ti_temp, te_temp);

		criteria = ee_temp - ee + tau/dm*dv_temp*pe_temp - tau*src;
		if (criteria > 0.0)
		{
			highBorder = te_temp;
		}
		else if (criteria < 0.0)
		{
			lowBorder  = te_temp;
		}
		else
		{
			return te_temp;
		}

	} while( fabs(criteria) > epsE );


	double val = ee_temp - ms[i].ee + tau/ms[i].dm*dv_temp*pe_temp - tau*src;


	return te_temp;
}


void CSolver::testHeatStage(char* inputFileName, char* outputFileName)
{
	/***************************************************************
		Тест для линейного модельного уравнения теплопроводности.
		Для осуществления теста требуется:

  1. Входной файл \input_files\input_lagrange_heat_test.txt.
  2. Функция getce() в табличном уравлении состояния 
     должны возвращать единицу.
  3. Функция getkappa() в табличном уравнении состояния должна возвращать 1e-8. 
  4. Должно быть задано гауссово (~exp(-x2)) начальное условие -- 
     задается в коде ниже.

	****************************************************************/

#define PI 3.14159

	task.load(inputFileName);
	ms.initData(&task);
	EOSOld &eos = task.getEOS();
	

	double t   = 1.0e-8;
	double tau = 1.0e-10;
	double kappa = 0.0;
	int i =0;

	// Задаем начальное условие (п.4 требований к тесту)

	for(i=0; i<ms.getSize(); i++)
	{
		kappa = eos.getkappa(ms[i].ro, ms[i].ti, ms[i].te);
		
		double te = 1.0 / sqrt( t * PI * kappa ) *
				          exp( - (ms[i].x-150.0e-9) * (ms[i].x-150.0e-9) / 4.0 / kappa / t);
		ms[i].te = te;
	}


	FILE *f1  = fopen(outputFileName, "w+");
    fprintf(f1, "TITLE=\"Heat test: Te(x), Te_an(x)\"\n");
	fprintf(f1, "VARIABLES=\"x\",\"Te\",\"Te_an\"\n");

	for(i=0; i<400; i++)
	{
		calcHeatStage(t, tau);
		t+=tau;
	}	

	for(i=0; i<ms.getSize(); i++)
	{
		kappa = eos.getkappa(ms[i].ro, ms[i].ti, ms[i].te);
		
		fprintf(f1, "%e %e %e\n", ms[i].x, ms[i].te, 
		    1.0 / sqrt( t * PI * kappa ) * exp( - (ms[i].x-150.0e-9) * (ms[i].x-150.0e-9) / 4.0 / kappa / t));
	}

	fclose(f1);
}


void CSolver::testExchangeStage(char* inputFileName, char* outputFileName)
{
	/***************************************************************
		Тест для линейной модельной системы уравнений обмена.
		Для осуществления теста требуется:

  1. Входной файл \input_files\input_lagrange_exchange_test.txt.
  2. Функции getce(), getci() в табличном уравлении состояния 
     должны возвращать единицу.
  3. Функция getAlpha() в табличном уравнении состояния должна возвращать -1. 
  
	****************************************************************/

	task.load(inputFileName);
	ms.initData(&task);
	EOSOld &eos = task.getEOS();

	double t   = 0;
	double tau = 5.0e-2;

	FILE *f1  = fopen(outputFileName, "w+");
    fprintf(f1, "TITLE=\"Exchange test: Ti(t), Te(t), Ti_an(t), Te_an(t)\"\n");
	fprintf(f1, "VARIABLES=\"t\",\"Ti\",\"Te\",\"Ti_an\",\"Te_an\"\n");

	int k = ms.getSize()/2;

	double te=ms[k].te,
		   ti=ms[k].ti,
		   ti_an = 3000.0 - 1000*exp(-2.0*t), 
		   te_an = 3000.0 + 1000*exp(-2.0*t);

	fprintf(f1, "%e %e %e %e %e\n", t, ms[k].ti, ms[k].te,
			    3000.0 - 1000*exp(-2.0*t), 
				3000.0 + 1000*exp(-2.0*t) );

	for(int i=0; i<100; i++)
	{
		calcExchangeStage(tau);
		t+=tau;
			
		fprintf(f1, "%e %e %e %e %e\n", t, ms[k].ti, ms[k].te,
			    3000.0 - 1000*exp(-2.0*t), 
				3000.0 + 1000*exp(-2.0*t) );
	}

	fclose(f1);
}


void CSolver::sweep(double *u, double *A, double *B, double *C, double *F, int size)
{
	double	*alpha = new double[size],
			 *beta = new double[size];
	int i=0;

	alpha[0] = 0.0;
	 beta[0] = 0.0;
		
	for(i=0; i<size-1; i++)
	{
		alpha[i+1] = B[i] / (C[i]-A[i]*alpha[i]);
		 beta[i+1] = (A[i]*beta[i]+F[i]) / (C[i]-A[i]*alpha[i]);
	}

	// Обратная прогонка

//	teR = (AR/CR*betaR + FR/CR) / (1.0 - AR/CR*alphaR);

	u[size-1] = (A[size-1]/C[size-1]*beta[size-1] + F[size-1]/C[size-1]) /
		        (1.0 - A[size-1]/C[size-1]*alpha[size-1]);

	for(i=size-2; i>=0; i--)
	{
		u[i] = alpha[i+1] * u[i+1] + beta[i+1];
	}

	delete []alpha;
	delete []beta;

}

void CSolver::calcHydroStageGodunov(double t, double tau) {
	double E = 0.;
	Vector4 Fm = Vector4::ZERO, Fp = Vector4::ZERO;
	int nSize = ms.getSize();
	EOSOld &eos=task.getEOS();
	double h = ms[1].x-ms[0].x;
	double gamma=eos.getGamma();
	int i=0;
	// Векторы W уже заполнены
	CVectorPrimitive res;
	// Потоки считаем по решению задачи Римана о распаде разрыва между ячейками
	for(i=0; i<nSize; i++) {
		Node& n=ms[i];
		if(i!=0) {
			Node& nm=ms[i-1];
			res = calcRPAnalyticalSolution(nm.ro, nm.v, nm.p, n.ro, n.v, n.p, 0., tau);
		} else
		    // Левое граничное условие -- прозрачная граница
			res = calcRPAnalyticalSolution(n.ro, n.v, n.p, n.ro, n.v, n.p, 0., tau);
		if(res.ro!=0.) E = res.p/(gamma-1.)/res.ro + 0.5*res.v*res.v; else E=0.;
		Fm = Vector4(res.ro*res.v, res.ro*res.v*res.v + res.p, res.v*(res.ro*E+res.p), 0.);
		if(i!=nSize-1) {
			Node& np=ms[i+1];
			res = calcRPAnalyticalSolution(n.ro, n.v, n.p, np.ro, np.v, np.p, 0., tau);
		} else 
			// Правое граничное условие -- прозрачная граница
			res = calcRPAnalyticalSolution(n.ro, n.v, n.p, n.ro, n.v, n.p, 0., tau);
		if(res.ro!=0.) E = res.p/(gamma-1.)/res.ro + 0.5*res.v*res.v; else E=0.;
		Fp = Vector4(res.ro*res.v, res.ro*res.v*res.v + res.p, res.v*(res.ro*E+res.p), 0.);
		n.W_temp = n.W - tau/h*(Fp-Fm);
	}
	// Тест решения
/*	CVectorPrimitive Um, Up;
	double Em=0., Ep=0.;
	double L0=0., R0=0., L1 = 0., L2=0., L3=0., R1 = 0., R2=0., R3=0.;
	double eps = 1.e-2;
	double epsAbs = 1.e-7;
	for(i=1; i<nSize-1; i++) {
		Um = calcRPAnalyticalSolution(ms[i-1].ro, ms[i-1].v, ms[i-1].p, ms[i].ro, ms[i].v, ms[i].p, 0., tau);
		Em = Um.p/(gamma-1)/Um.ro + 0.5*Um.v*Um.v;
		Fm = Vector4(Um.ro*Um.v, Um.ro*Um.v*Um.v + Um.p, Um.v*(Um.ro*Em+Um.p), 0.);
		Up = calcRPAnalyticalSolution(ms[i].ro, ms[i].v, ms[i].p, ms[i+1].ro, ms[i+1].v, ms[i+1].p, 0., tau);
		Ep = Up.p/(gamma-1)/Up.ro + 0.5*Up.v*Up.v;
		Fp = Vector4(Up.ro*Up.v, Up.ro*Up.v*Up.v + Up.p, Up.v*(Up.ro*Ep+Up.p), 0.);
		double L0 = (ms[i].W_temp[0] - ms[i].W[0])*h;
		double R0 = -(Fp[0]-Fm[0])*tau;
		double L1 = (ms[i].W_temp[1] - ms[i].W[1])*h;
		double R1 = -(Fp[1]-Fm[1])*tau;
		double L2 = (ms[i].W_temp[2] - ms[i].W[2])*h;
		double R2 = -(Fp[2]-Fm[2])*tau;
		if(fabs((L0-R0)/0.5/(L0+R0)) > eps && L0 != 0. && R0 != 0. && fabs(L0-R0)>epsAbs)
			cout << "calcHydroStageGodunov() test not passed: mesh conversion law #0 not satisfied in node " << i << endl;
		if(fabs((L1-R1)/0.5/(L1+R1))>eps && L1 != 0. && R1 != 0. && fabs(L1-R1)>epsAbs)
			cout << "calcHydroStageGodunov() test not passed: mesh conversion law #1 not satisfied in node " << i << endl;
		if(fabs((L2-R2)/0.5/(L2+R2))>eps && L2 != 0. && R2 != 0. && fabs(L2-R2)>epsAbs)
			cout << "calcHydroStageGodunov() test not passed: mesh conversion law #2 not satisfied in node " << i << endl;
	}
	// */
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
		n.ti = eos.getti(n.ro, n.e);
		n.p  = eos.getpi(n.ro, n.ti);
		n.C  = eos.getC(n.ro, n.ti, n.te);
	}
	cout << "calcHydroStageGodunov(): done!" << endl;
}

void CSolver::calcHydroStageGodunovMovingMesh(double t, double tau)
{
	// Версия солвера метода Годунова 1-го порядка на подвижной однородной сетке
	double E = 0.;
	Vector4 Fm = Vector4::ZERO, Fp = Vector4::ZERO;
	int nSize = ms.getSize();
	EOSOld &eos=task.getEOS();
	double h = 0.;
	double hNew = 0.;
	double gamma=eos.getGamma();
	int i=0;
	/// DEBUG ///
	double s=0., sNew=0., sTemp=0., sFluxM=0., sFluxP=0.;
	/////////////
	// Векторы W уже заполнены
	CVectorPrimitive res;
	res.ro = 0.; res.v = 0.; res.p = 0.;
	// Считаем скорости узлов сетки
	CMethod &method = task.getMethod();
	// Два варианта: либо константа из начальных условий (скорость левой границы, "хвоста" волны разрежения), либо обновлять на каждом шаге по времени из ВР между крайней ячейкой и вакуумом.
	// Первый вариант
	Zone z = task.getZone(0);
	double vLeft = z.v - 2.*sqrt(gamma*eos.getpi(z.ro, z.ti)/z.ro)/(gamma-1.);


	// Второй вариант
	//double vLeft = ms[0].v - 2.*sqrt(gamma*ms[0].p/ms[0].ro)/(gamma-1.);
	// Третий вариант: решать задачу о распаде разрыва и давать скорость левой границы
	// ???
	// Значения сеточных координат и скоростей на новом временном слое
	for(i=0; i<nSize+1; i++) {
		method.vGrid[i] = vLeft*(1.0-double(i)/double(nSize));
		method.X[i] = ms[i].x + method.vGrid[i]*tau;
	}

	/// DEBUG ///
	// Переинтерполяция узлов сетки
	// Авантюра?
	/*Vector4 lTerm=Vector4::ZERO, rTerm=Vector4::ZERO;
	for(i=0; i<nSize; i++){
		Node& n = ms[i];
		if(i==0) lTerm = 0.; else lTerm = ms[i-1].W*fabs(ms[i].x-method.X[i]);
		if(i==nSize-1) rTerm=0.; else rTerm = ms[i].W*fabs(method.X[i+1]-ms[i].x);
		n.W = (lTerm + rTerm)/(method.X[i+1]-method.X[i]);
		method.updateNode(n);
	}*/

	/////////////

	// Использую метод Munz
	xVacBound += vLeft*tau;
	double eps = 0.01;

	// Потоки считаем по решению задачи Римана о распаде разрыва между ячейками
	for(i=0; i<nSize; i++) {
		h    = ms[i+1].x - ms[i].x;
		hNew = method.X[i+1]-method.X[i];
		Node& n=ms[i];
		if(i!=0) {
			Node& nm=ms[i-1];
			res = calcRPAnalyticalSolution(nm.ro, nm.v, nm.p, n.ro, n.v, n.p, method.X[i]-ms[i].x, tau);
			if(res.ro!=0.) E = res.p/(gamma-1.)/res.ro + 0.5*res.v*res.v; else E=0.;
			Fm = Vector4(res.ro*(res.v-method.vGrid[i]), res.ro*res.v*(res.v-method.vGrid[i]) + res.p, res.ro*E*(res.v-method.vGrid[i])+res.p*res.v, 0.);
			/// DEBUG ///
			s=eos.getEntropy(res.ro, eos.getti(res.ro, res.p/(gamma-1.)/res.ro));
			sFluxM = res.ro * s * (res.v-method.vGrid[i]);
			/////////////
		} else {
		    // Левое граничное условие -- прозрачная граница
			// Один из двух вариантов -- либо нулевой поток всегда, либо решение задачи о распаде разрыва, во втором случае поток зависит от скорости левой границы сетки 
			// Первый вариант
			//Fm = Vector4::ZERO;
			/// DEBUG ///
			// Второй вариант
			res = calcRPAnalyticalSolution(0., vLeft, 0., n.ro, n.v, n.p, method.X[i]-ms[i].x, tau);
			//res = calcRPAnalyticalSolution(n.ro, /*vLeft*/n.v, n.p, n.ro, n.v, n.p, method.X[i]-ms[i].x, tau);
			if(res.ro!=0.) E = res.p/(gamma-1.)/res.ro + 0.5*res.v*res.v; else E=0.;
			//Fm = Vector4::ZERO; //Vector4(res.ro*(res.v-method.vGrid[i]), res.ro*res.v*(res.v-method.vGrid[i]) + res.p, res.ro*E*(res.v-method.vGrid[i])+res.p*res.v, 0.);
			if( fabs(xVacBound - ms[0].x) < eps*h ) {
				Fm = Vector4::ZERO;
			} else {
				if( fabs(ms[0].v) - ms[0].C > 0 ) {
					res = calcRPAnalyticalSolution(0., vLeft, 0., ms[0].ro, ms[0].v, ms[0].p, method.X[i]-ms[i].x, tau);
					Fm = Vector4(res.ro*(res.v-method.vGrid[0]), res.ro*res.v*(res.v-method.vGrid[i]) + res.p, res.ro*E*(res.v-method.vGrid[i])+res.p*res.v, 0.);
				}
				else { 
					res = calcRPAnalyticalSolution(0., vLeft, 0., 
						                           task.getZone(0).ro, task.getZone(0).v, eos.getp(task.getZone(0).ro, task.getZone(0).ti, 0.),
												   method.X[i]-ms[i].x, tau);
					Fm = Vector4(res.ro*(res.v-method.vGrid[0]), res.ro*res.v*(res.v-method.vGrid[i]) + res.p, res.ro*E*(res.v-method.vGrid[i])+res.p*res.v, 0.);
				}
			}
			
			
			/// DEBUG ///
			//Fm = Vector4(res.ro*(res.v-method.vGrid[i]), res.ro*res.v*(res.v-method.vGrid[i]) + res.p, res.ro*E*(res.v-method.vGrid[i])+res.p*res.v, 0.);
			// Энтропия и поток энтропии			
			//s=eos.getEntropy(res.ro, eos.getti(res.ro, res.p/(gamma-1.)/res.ro), 0.);
			//sFluxM = res.ro *s* (res.v-method.vGrid[i]);
			/////////////
			//Fm += Vector4(0., -0.00001, 0., 0.);
			////
			int qq = 0;
			////
		}
		if(i!=nSize-1) {
			Node& np=ms[i+1];
			res = calcRPAnalyticalSolution(n.ro, n.v, n.p, np.ro, np.v, np.p, method.X[i+1]-ms[i+1].x, tau);
		} else 
			// Правое граничное условие -- прозрачная граница
			res = calcRPAnalyticalSolution(n.ro, n.v, n.p, n.ro, n.v, n.p, method.X[i+1]-ms[i+1].x, tau);
		if(res.ro!=0.) E = res.p/(gamma-1.)/res.ro + 0.5*res.v*res.v; else E=0.;
		Fp = Vector4(res.ro*(res.v-method.vGrid[i+1]), res.ro*res.v*(res.v-method.vGrid[i+1]) + res.p, res.ro*E*(res.v-method.vGrid[i+1])+res.p*res.v, 0.);
		n.W_temp   = (h*n.W - tau*(Fp-Fm))/hNew;
	}
	// Тест решения
/*	CVectorPrimitive Um, Up;
	double Em=0., Ep=0.;
	double L0=0., R0=0., L1 = 0., L2=0., L3=0., R1 = 0., R2=0., R3=0.;
	double eps = 1.e-2;
	double epsAbs = 1.e-7;
	for(i=1; i<nSize-1; i++) {
		Um = calcRPAnalyticalSolution(ms[i-1].ro, ms[i-1].v, ms[i-1].p, ms[i].ro, ms[i].v, ms[i].p, 0., tau);
		Em = Um.p/(gamma-1)/Um.ro + 0.5*Um.v*Um.v;
		Fm = Vector4(Um.ro*Um.v, Um.ro*Um.v*Um.v + Um.p, Um.v*(Um.ro*Em+Um.p), 0.);
		Up = calcRPAnalyticalSolution(ms[i].ro, ms[i].v, ms[i].p, ms[i+1].ro, ms[i+1].v, ms[i+1].p, 0., tau);
		Ep = Up.p/(gamma-1)/Up.ro + 0.5*Up.v*Up.v;
		Fp = Vector4(Up.ro*Up.v, Up.ro*Up.v*Up.v + Up.p, Up.v*(Up.ro*Ep+Up.p), 0.);
		double L0 = (ms[i].W_temp[0] - ms[i].W[0])*h;
		double R0 = -(Fp[0]-Fm[0])*tau;
		double L1 = (ms[i].W_temp[1] - ms[i].W[1])*h;
		double R1 = -(Fp[1]-Fm[1])*tau;
		double L2 = (ms[i].W_temp[2] - ms[i].W[2])*h;
		double R2 = -(Fp[2]-Fm[2])*tau;
		if(fabs((L0-R0)/0.5/(L0+R0)) > eps && L0 != 0. && R0 != 0. && fabs(L0-R0)>epsAbs)
			cout << "calcHydroStageGodunov() test not passed: mesh conversion law #0 not satisfied in node " << i << endl;
		if(fabs((L1-R1)/0.5/(L1+R1))>eps && L1 != 0. && R1 != 0. && fabs(L1-R1)>epsAbs)
			cout << "calcHydroStageGodunov() test not passed: mesh conversion law #1 not satisfied in node " << i << endl;
		if(fabs((L2-R2)/0.5/(L2+R2))>eps && L2 != 0. && R2 != 0. && fabs(L2-R2)>epsAbs)
			cout << "calcHydroStageGodunov() test not passed: mesh conversion law #2 not satisfied in node " << i << endl;
	}*/
	// Тест закона сохранения энергии на сетке
/*	i=17;
	double G=0., G_plus=0.;
	double roL=0., vL=0., pL=0., roR=0., vR=0., pR=0.;
	double e = 0.;
	double FEp = 0., FEm=0.;
	double L=0, R=0;
	for(i=1; i<nSize-1; i++) {
		G      = ms[i+1].x-ms[i].x;
		G_plus = method.X[i+1]-method.X[i];
		roL=ms[i].ro; vL=ms[i].v; pL=ms[i].p; roR=ms[i+1].ro; vR=ms[i+1].v; pR=ms[i+1].p;
		res  = calcRPAnalyticalSolution(roL, vL, pL, roR, vR, pR, method.X[i+1]-ms[i+1].x, tau);
		e = res.p/(gamma-1.)/res.ro;
		FEp = res.v*(res.p + res.ro*e + 0.5*res.ro*res.v*res.v) - method.vGrid[i+1]*(res.ro*e + 0.5*res.ro*res.v*res.v);
		
		double qq = ms[0].v;
		roL=ms[i-1].ro; vL=ms[i-1].v; pL=ms[i-1].p; roR=ms[i].ro; vR=ms[i].v; pR=ms[i].p;
		res = calcRPAnalyticalSolution(roL, vL, pL, roR, vR, pR, method.X[i]-ms[i].x, tau);
		e = res.p/(gamma-1.)/res.ro;
		E = e + 0.5*res.v*res.v;
		FEm = res.v*(res.p + res.ro*e + res.ro*res.v*res.v/2.) - method.vGrid[i]*(res.ro*e + res.ro*res.v*res.v/2.);
		L = G_plus*ms[i].W_temp[2]-G*ms[i].W[2];
		R = -(FEp-FEm)*tau;
		if(fabs( (L-R)/(.5*(L+R)) ) > 0.002 && fabs(L)>1.e-6 && fabs(R)>1.e-6)
			cout << "I found it! I waited so long! -- i=" << i << endl;
	}*/

	// Значения на новом временнОм слое
	for(i=0; i<nSize; i++) 
	{
		Node& n=ms[i];
		n.W[0] = n.W_temp[0];
		n.W[1] = n.W_temp[1];
		n.W[2] = n.W_temp[2];
	/*	if(n.W[0]<TOL) {
			n.ro = 0.;
			n.v  = vLeftVac;
			n.p  = 0.;
			n.e  = 0.;
			n.ti = 0.;
			n.C	 = 0.;
		} else { */
			n.ro = n.W[0];
			n.v  = n.W[1]/n.ro;
			if(n.v < vLeft) 
			{
				int qq =0;
			}
			E = n.W[2]/n.ro; 
			n.e  = n.W[2]/n.ro - 0.5*n.v*n.v;
			if(n.e < 0. ) {
				cout << "Warning: negative internal energy in cell i=" << i << ". Probably numerical instability." << endl;
				int qq = 0;
			}
			n.ti = eos.getti(n.ro, n.e);
			n.p  = eos.getpi(n.ro, n.ti);
			n.C  = eos.getC(n.ro, n.ti, n.te);
	}

	// И проставляем новые координаты узлов сетки
	for(i=0; i<nSize+1; i++){
		ms[i].x = method.X[i];
	}
	cout << "calcHydroStageGodunovMovingMesh(): done!" << endl;
}

void CSolver::calcHydroStageMHM(double t, double tau) {
	EOSOld &eos = task.getEOS();
	double E=0.; 
	int i=0;
	int nSize = ms.getSize();
	double h=ms[1].x-ms[0].x;
	Vector4 L = Vector4::ZERO, R = Vector4::ZERO, D = Vector4::ZERO, V = Vector4::ZERO;
	double gamma = eos.getGamma();
	// Значения на границе (в фиктивных ячейках)
	double roLB=0., vLB = 0., ELB=0., roRB=0., vRB=0., ERB=0.;
	// Transmissive b.c.s
	Vector4 _W_LEFT= ms[0].W, _W_RIGHT = ms[nSize-1].W;	
	// Data reconstruction
	// Не забываем, что координатная таблица наклонов "сдвинута" на 2 ячейки (порядок схемы = количество фиктивных ячеек) 
	int nOrder = 2;
	double omega = .5;
	double r = 0.;
	double TOL = 1.e-5; Vector4 minDiff = Vector4(TOL, TOL, TOL, 0.);
	std::vector<Vector4>slope(nSize+nOrder+nOrder), _U_L(nSize+nOrder+nOrder), _U_R(nSize+nOrder+nOrder);
	slope[0] = minDiff; _U_L[0] = Vector4::ZERO; _U_R[0] = Vector4::ZERO;
	i=-1;
	//slope[i+nOrder] = .5*(1.+omega)*(_W_LEFT - _W_LEFT) + .5*(1.-omega)*(ms[i+1].W - _W_LEFT);
	//slope[i+nOrder] = calcMinmodSlope(Vector4::ZERO, ms[i+1].W - ms[i].W);
	//slope[i+nOrder] = calcSuperBEESlope(omega, Vector4::ZERO, ms[i+1].W -_W_LEFT);
	slope[i+nOrder] = calcMinmodSlope(Vector4::ZERO, ms[i+1].W -_W_LEFT);
	_U_L[i+nOrder] = _W_LEFT - .5*slope[i+nOrder]; _U_R[i+nOrder] = _W_LEFT + .5*slope[i+nOrder];
	i=0;
	//slope[i+nOrder] = .5*(1.+omega)*(ms[i].W - _W_LEFT) + .5*(1.-omega)*(ms[i+1].W - ms[i].W);
	//slope[i+nOrder] = calcMinmodSlope(ms[i].W - _W_LEFT, ms[i+1].W - ms[i].W);
	//slope[i+nOrder] = calcSuperBEESlope(omega, ms[i].W - _W_LEFT, ms[i+1].W - ms[i].W);
	slope[i+nOrder] = calcMinmodSlope(ms[i].W - _W_LEFT, ms[i+1].W - ms[i].W);
	_U_L[i+nOrder] = ms[i].W-.5*slope[i+nOrder];
	_U_R[i+nOrder] = ms[i].W+.5*slope[i+nOrder];
	i=nSize-1;
	//slope[i+nOrder] = .5*(1.+omega)*(ms[i].W - ms[i-1].W) + .5*(1.-omega)*(_W_RIGHT - ms[i].W);
	//slope[i+nOrder] = calcMinmodSlope(ms[i].W - ms[i-1].W, _W_RIGHT - ms[i].W);
	//slope[i+nOrder] = calcSuperBEESlope(omega, ms[i].W - ms[i-1].W, _W_RIGHT - ms[i].W);
	slope[i+nOrder] = calcMinmodSlope(ms[i].W - ms[i-1].W, _W_RIGHT - ms[i].W);
	_U_L[i+nOrder] = ms[i].W-.5*slope[i+nOrder];
	_U_R[i+nOrder] = ms[i].W+.5*slope[i+nOrder];
	i=nSize;
	//slope[i+nOrder] = .5*(1.+omega)*(_W_RIGHT - ms[i-1].W) + .5*(1.-omega)*(_W_RIGHT - _W_RIGHT);
	//slope[i+nOrder] = calcSuperBEESlope(omega, _W_RIGHT - ms[i-1].W, _W_RIGHT - _W_RIGHT);
	slope[i+nOrder] = calcMinmodSlope(_W_RIGHT - ms[i-1].W, _W_RIGHT - _W_RIGHT);
	
	//slope[i+nOrder] = calcMinmodSlope(_W_RIGHT - ms[i-1].W, _W_RIGHT - _W_RIGHT);
	_U_L[i+nOrder] = _W_RIGHT-.5*slope[i+nOrder];
	_U_R[i+nOrder] = _W_RIGHT+.5*slope[i+nOrder];
	slope[nSize+1+nOrder] = Vector4::ZERO; _U_L[nSize+1+nOrder] = Vector4::ZERO; _U_R[nSize+1+nOrder] = Vector4::ZERO;
	for(i=1; i<nSize-1; i++) {
		////
		if (i==79) {
			Vector4 V1 = ms[i].W, V2 = ms[i-1].W, V3 = ms[i+1].W, delta = V3-V1;
			double qq=0.;
		}
		////
		//slope[i+nOrder] = (.5*(1.+omega)*(ms[i].W - ms[i-1].W) + .5*(1.-omega)*(ms[i+1].W - ms[i].W));
		//slope[i+nOrder] = calcMinmodSlope(ms[i].W - ms[i-1].W, ms[i+1].W - ms[i].W);
		//slope[i+nOrder]=calcSuperBEESlope(omega, ms[i].W - ms[i-1].W, ms[i+1].W - ms[i].W);
		slope[i+nOrder]=calcMinmodSlope(ms[i].W - ms[i-1].W, ms[i+1].W - ms[i].W);
		_U_L[i+nOrder] = ms[i].W-.5*slope[i+nOrder];
		_U_R[i+nOrder] = ms[i].W+.5*slope[i+nOrder];
	} 
	///
	if(t==0.) {
		cout << "28:" << slope[28+nOrder]<<endl;
		cout << "29:" << slope[29+nOrder]<<endl;
		cout << "30:" << slope[30+nOrder]<<endl;
		cout << "31:" << slope[31+nOrder]<<endl;
	}

	// Vizualizing MUSCL linear data reconstruction
	/*string fName = string("reconstr.dat");
	ofstream ofs;
	ofs.open(fName, ios::out);
	double lTotal = h*nSize;
	int j=0, nZoom = 10;
	int nDetSize = nSize*nZoom;
	double hDet = lTotal/nDetSize;
	double _x_i = 0., _x=0., _ro=0., _rou=0., _roE=0.;
	for(i=0; i<nSize; i++) {
		_x_i = h*i + h/2.;
		for(j=0; j<nZoom; j++) {
			_x = _x_i - h/2. + hDet/2. + j*hDet;
			_ro = ms[i].ro + (_x-_x_i)/h*slope[i+2][0];
			_rou = ms[i].ro*ms[i].v + (_x-_x_i)/h*slope[i+2][1];
			_roE = ms[i].ro*(ms[i].e+0.5*ms[i].v*ms[i].v) + (_x-_x_i)/h*slope[i+2][2];
			ofs.setf(ios::scientific);
			ofs << _x << " " << _ro << " " << _rou << " " << _roE << endl;
		}
	}

	ofs.close();*/
	///


	// Data evolution
	for(i=0; i<nSize+nOrder+nOrder; i++) {
		Vector4 FL = calcPhysicalFlux(_U_L[i][0], _U_L[i][1], _U_L[i][2]),
			    FR = calcPhysicalFlux(_U_R[i][0], _U_R[i][1], _U_R[i][2]);
		_U_L[i] += .5*tau/h*(FL - FR);
		_U_R[i] += .5*tau/h*(FL - FR);		
	}
	// Godunov flux, based on exact RP solution
	for(i=0; i<nSize; i++) {
		////
		if(i==80) {
		    double qq=0.;	
		}
		////
		
		ms[i].F = calcGodunovFlux(_U_R[i-1+nOrder][0], _U_R[i-1+nOrder][1], _U_R[i-1+nOrder][2], 
			                      _U_L[i+nOrder][0], _U_L[i+nOrder][1], _U_L[i+nOrder][2]);

	}
	i = nSize;
	ms[i].F = calcGodunovFlux(_U_R[i-1+nOrder][0], _U_R[i-1+nOrder][1], _U_R[i-1+nOrder][2], 
			                      _U_L[i+nOrder][0], _U_L[i+nOrder][1], _U_L[i+nOrder][2]);
	// Riemann problem		
	// Main cycle
	for(i=0; i<nSize; i++) {				
		
		
		
		////
		if(i==80) {
		    double qq=0.;	
		}
		////



		
		
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

// "Godunov second order" method
// Riemann solver based on exact solution + simple LSQ reconstruction
void CSolver::calcHydroStageG2(double t, double tau) {
	EOSOld &eos = task.getEOS(); 
	double gamma = eos.getGamma();
	double E=0.; 
	int i=0, j=0, nSize = ms.getSize();
	double h=ms[1].x-ms[0].x;
	double _x = 0.;
	Vector4 L = Vector4::ZERO, R = Vector4::ZERO, D = Vector4::ZERO, V = Vector4::ZERO;	
	// Значения на границе (в фиктивных ячейках)
	double roim1=0., uim1=0., Eim1=0., roim2=0., uim2=0., Eim2=0.,
		   roip1=0., uip1=0., Eip1=0., roip2=0., uip2=0., Eip2=0.;	
	Vector4 slopeim1=Vector4::ZERO, slopeip1=Vector4::ZERO;
	// B.c.'s
	// Transmissive left boundary
	roim1 = ms[0].ro; uim1 = ms[0].v; Eim1 = ms[0].e + .5*ms[0].v*ms[0].v;
	roim2 = ms[0].ro; uim2 = ms[0].v; Eim2 = ms[0].e + .5*ms[0].v*ms[0].v;
	// Transmissive right boundary
	roip1= ms[nSize-1].ro; uip1 = ms[nSize-1].v; Eip1 = ms[nSize-1].e + .5*ms[nSize-1].v*ms[nSize-1].v;
	roip2= ms[nSize-1].ro; uip2 = ms[nSize-1].v; Eip2 = ms[nSize-1].e + .5*ms[nSize-1].v*ms[nSize-1].v;
	Vector4 Wim1 = Vector4(roim1, roim1*uim1, roim1*Eim1, 0.), Wim2 = Vector4(roim2, roim2*uim2, roim2*Eim2, 0.),
		    Wip1 = Vector4(roip1, roip1*uip1, roip1*Eip1, 0.), Wip2 = Vector4(roip2, roip2*uip2, roip2*Eip2, 0.);
	// Reconstruction
	// slope[i] = (f(x_i+1)-f(x_i-1))/(2*dx)
	std::vector<Vector4>slope(nSize);	
	/*slopeim1[0] = (ms[0].W[0] - roim2)/2./h;
	slopeim1[1] = (ms[0].W[1] - roim2*uim2)/2./h;
    slopeim1[2] = (ms[0].W[2] - roim2*Eim2)/2./h;
	slopeim1[3] = 0.;
	slopeip1[0] = (roip2      - ms[nSize-1].W[0])/2./h;
	slopeip1[1] = (roip2*uip2 - ms[nSize-1].W[1])/2./h;
	slopeip1[2] = (roip2*Eip2 - ms[nSize-1].W[2])/2./h;
	slopeip1[3] = 0.;
	for(i=0; i<nSize-1; i++) {
		if(i==0) {
			slope[i][0] = (ms[i+1].W[0] - roim1)/2./h;
			slope[i][1] = (ms[i+1].W[1] - roim1*uim1)/2./h;
			slope[i][2] = (ms[i+1].W[2] - roim1*Eim1)/2./h;
			slope[i][3] = 0.;

		} else if(i==nSize-1) {
			slope[i][0] = (roip1      - ms[i-1].W[0])/2./h;
			slope[i][1] = (roip1*uip1 - ms[i-1].W[1])/2./h;
			slope[i][2] = (roip1*Eip1 - ms[i-1].W[2])/2./h;
			slope[i][3] = 0.;
		} else 
			for(j=0; j<4; j++)
				slope[i][j] = (ms[i+1].W[j]-ms[i-1].W[j])/2./h;
	}*/
	// Minmod monotone reconstruction
	slopeim1 = calcMinmodSlope(Wim1-Wim2, ms[0].W-Wim1);
	slope[0] = calcMinmodSlope(ms[0].W-Wim1, ms[1].W-ms[0].W);
	slope[nSize-1] = calcMinmodSlope(ms[nSize-1].W-ms[nSize-2].W, Wip1-ms[nSize-1].W);
	slopeip1 = calcMinmodSlope(Wip1-ms[nSize-1].W, Wip2-Wip1);
	for(i=1; i<nSize-1; i++)
		slope[i] = calcMinmodSlope(ms[i].W-ms[i-1].W, ms[i+1].W-ms[i].W);	
	// Intercell fluxes
	for(i=0; i<nSize+1; i++) {
		Node &n = ms[i];	
		if(i==0)			
			n.F = calcGodunovFlux(roim1+slopeim1[0]*0.5*h, roim1*uim1+slopeim1[1]*0.5*h, roim1*Eim1+slopeim1[2]*0.5*h, 
			                     n.W[0]-slope[i][0]*0.5*h,     n.W[1]-slope[i][1]*0.5*h,     n.W[2]-slope[i][2]*0.5*h);
		else if(i==nSize)
			n.F = calcGodunovFlux(ms[i-1].W[0]+slope[i-1][0]*0.5*h, ms[i-1].W[1]+slope[i-1][1]*0.5*h, ms[i-1].W[2]+slope[i-1][2]*0.5*h, 
			                             Wip1[0]-slopeip1[0]*0.5*h,      Wip1[1]-slopeip1[1]*0.5*h,        Wip1[2]-slopeip1[2]*0.5*h);
		else
			n.F = calcGodunovFlux(ms[i-1].W[0]+slope[i-1][0]*0.5*h, ms[i-1].W[1]+slope[i-1][1]*0.5*h, ms[i-1].W[2]+slope[i-1][2]*0.5*h,
								          n.W[0]-slope[i][0]*0.5*h,         n.W[1]-slope[i][1]*0.5*h,         n.W[2]-slope[i][2]*0.5*h);		
	}
	
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

// Lax-Friedrichs 1st order method
void CSolver::calcHydroStageLaxFriedrichs(double t, double tau) {
	EOSOld &eos = task.getEOS();
	unsigned int i=0;
	const unsigned int nSize = ms.getSize();
	const double h=ms[1].x-ms[0].x;
	Vector4 UL = Vector4::ZERO, UR = Vector4::ZERO, FL = Vector4::ZERO, FR = Vector4::ZERO;
	const double gamma = eos.getGamma();
	// Значения на границе (в фиктивных ячейках)
	double roLB=0., vLB = 0., ELB=0., roRB=0., vRB=0., ERB=0.;
	// Задаем граничные условия
	// Transmissive left boundary
	roLB = ms[0].ro; vLB = ms[0].v;	ELB = ms[0].e + .5*ms[0].v*ms[0].v;
	// Transmissive right boundary
	roRB = ms[nSize-1].ro; vRB = ms[nSize-1].v; ERB = ms[nSize-1].e + .5*ms[nSize-1].v*ms[nSize-1].v;
	for(i=0; i<=nSize; i++) {
		Node &n = ms[i];	
		if(i==0) {
			UL = Vector4(roLB, roLB*vLB, roLB*ELB, 0.); UR = ms[i].W;			
		} else if(i==nSize) {
			UL = ms[i-1].W;	UR = Vector4(roRB, roRB*vRB, roRB*ERB, 0.); 
		} else {
			UL = ms[i-1].W;	UR = ms[i].W; 
		}
		FL = calcPhysicalFlux(UL[0], UL[1], UL[2]);
		FR = calcPhysicalFlux(UR[0], UR[1], UR[2]);
		n.F = .5*(FL + FR) - 0.5*h/tau*(UR - UL);
	}
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

Vector4 CSolver::calcMinmodSlope(Vector4 deltaMinus, Vector4 deltaPlus) {
	Vector4 res = Vector4::ZERO;
	for(int i=0; i<3; i++) {
		if(deltaPlus[i] > 0.) 
			res[i] = max(0., min(deltaMinus[i], deltaPlus[i])); 
		else
			res[i] = min(0., max(deltaMinus[i], deltaPlus[i])); 
	}
	return res;
}

Vector4 CSolver::calcMinmodSlopeModified(Vector4 deltaMinus, Vector4 deltaPlus) {
	Vector4 res = Vector4::ZERO;
	for(int i=0; i<3; i++) {
		if(fabs(deltaMinus[i]) < fabs(deltaPlus[i])) 
			res[i] = deltaMinus[i]; 
		else 
			res[i] = deltaPlus[i]; 
	}
	return res;
}

Vector4 CSolver::calcSuperBEESlope(double omega, Vector4 deltaMinus, Vector4 deltaPlus) {
	/*SUBROUTINE SBSLIC(R, OMEGA, DELTA)
	  REAL  DELTA, DENOR, OMEGA, PHI, PHIR, R
      PHI             = 0.0
      IF(R.GE.0.0)PHI = 2.0*R
      IF(R.GE.0.5)PHI = 1.0
      IF(R.GE.1.0)THEN
         DENOR = 1.0 - OMEGA + (1.0 + OMEGA)*R
         PHIR  = 2.0/DENOR
         PHI   = MIN(PHIR, R)
         PHI   = MIN(PHI, 2.0)
      ENDIF
      DELTA = PHI*DELTA
      END */
	int i = 0;
	double fi = 0., fiR=0., r = 0., TOL = 1.e-5;
	Vector4 res = .5*((1.+omega)*deltaMinus + (1-omega)*deltaPlus);
	for(i=0; i<4; i++) {
		fi = 0;
		fiR = 0;
		if(fabs(deltaPlus[i])<TOL && deltaPlus[i]>=0) deltaPlus[i] = TOL;
		if(fabs(deltaPlus[i])<TOL && deltaPlus[i]<0) deltaPlus[i] = -TOL;
		if(fabs(deltaMinus[i])<TOL && deltaMinus[i]>=0) deltaMinus[i] = TOL;
		if(fabs(deltaMinus[i])<TOL && deltaMinus[i]<0) deltaMinus[i] = -TOL;
		r = deltaMinus[i]/deltaPlus[i];
		if(r>0.) fi = 2.*r;
		if(r>.5) fi = 1.;
		if(r>1.) {
			fiR = 2./(1.-omega + (1+omega)*r);
			fi = min(fi, fiR);
			fi = min(fi, 2.);
		}
		res[i]*=fi;
	}
	return res;		
}

Vector4 CSolver::calcVanAlbadaSlope(double omega, Vector4 deltaMinus, Vector4 deltaPlus) {
	/*SUBROUTINE VASLIC(R, OMEGA, DELTA)
      IMPLICIT NONE
      REAL  DELTA, DENOR, OMEGA, PHI, PHIR, R
      PHI = 0.0
      IF(R.GE.0.0)THEN
         DENOR = 1.0 - OMEGA + (1.0 + OMEGA)*R
         PHIR  = 2.0/DENOR
         PHI   = R*(1.0 + R)/(1.0 + R*R)
         PHI   = MIN(PHI, PHIR)
      ENDIF
      DELTA    = PHI*DELTA
      END*/
	int i = 0;
	double fi = 0., fiR=0., r = 0., TOL = 1.e-5;
	Vector4 res = .5*((1.+omega)*deltaMinus + (1-omega)*deltaPlus);
	for(i=0; i<4; i++) {
		fi = 0;
		if(fabs(deltaPlus[i])<TOL && deltaPlus[i]>=0) deltaPlus[i] = TOL;
		if(fabs(deltaPlus[i])<TOL && deltaPlus[i]<0) deltaPlus[i] = -TOL;
		if(fabs(deltaMinus[i])<TOL && deltaMinus[i]>=0) deltaMinus[i] = TOL;
		if(fabs(deltaMinus[i])<TOL && deltaMinus[i]<0) deltaMinus[i] = -TOL;
		r = deltaMinus[i]/deltaPlus[i];
		if(r>0.) {
			fiR = 2./(1.-omega + (1.+omega)*r);
			fi = r*(1.+r)/(1+r*r);
			fi = min(fi, fiR);
		}
		res[i]*=fi;
	}
	return res;
}

Vector4 CSolver::calcVanLeerSlope(double omega, Vector4 deltaMinus, Vector4 deltaPlus) {
	/*SUBROUTINE VLSLIC(R, OMEGA, DELTA)
      REAL  DELTA, DENOR, OMEGA, PHI, PHIR, R
      PHI = 0.0
      IF(R.GE.0.0)THEN
         DENOR = 1.0 - OMEGA + (1.0 + OMEGA)*R
         PHIR  = 2.0/DENOR
         PHI   = 2.0*R/(1.0 + R)
         PHI   = MIN(PHI, PHIR)
      ENDIF
      DELTA    = PHI*DELTA
      END*/
	int i = 0;
	double fi = 0., fiR=0., r = 0., TOL = 1.e-5;
	Vector4 res = .5*((1.+omega)*deltaMinus + (1-omega)*deltaPlus);
	for(i=0; i<4; i++) {
		fi = 0;
		fiR = 0;
		if(fabs(deltaPlus[i])<TOL && deltaPlus[i]>=0) deltaPlus[i] = TOL;
		if(fabs(deltaPlus[i])<TOL && deltaPlus[i]<0) deltaPlus[i] = -TOL;
		if(fabs(deltaMinus[i])<TOL && deltaMinus[i]>=0) deltaMinus[i] = TOL;
		if(fabs(deltaMinus[i])<TOL && deltaMinus[i]<0) deltaMinus[i] = -TOL;
		r = deltaMinus[i]/deltaPlus[i];
		if(r>0.) {
			fi = 2.*r/(1.+r);
			fiR = 2./(1.-omega + (1+omega)*r);
			fi = min(fi, fiR);			
		}
		res[i]*=fi;
	}
	return res;		 


}



Vector4 CSolver::calcF(double ro, double v, double p, double e) {
	return Vector4(ro*v, p + ro*v*v, v*(p + ro*(e + 0.5*v*v)), 0.);
}

Vector4 CSolver::calcPhysicalFlux(double ro, double rou, double roE) {
	if (ro==0) return Vector4::ZERO; 
	EOSOld& eos = task.getEOS();
	double gamma = eos.getGamma();
	double u = rou/ro;
	double E = roE/ro;
	double p = (gamma-1.)*ro*(E-.5*u*u);
	return Vector4(rou, p + ro*u*u, u*(p + roE), 0.);
}

Vector4 CSolver::calcPhysicalFluxEOSBin(double ro, double rou, double roE) {
	if (ro==0) return Vector4::ZERO; 
	EOSBin* eos = task.eosBin;
	double gamma = eos->gamma;
	double u = rou/ro;
	double e = roE/ro-0.5*u*u;
	double p = eos->getp(ro, e);
	return Vector4(rou, p + ro*u*u, u*(p + roE), 0.);
}


void CSolver::calcHydroStageMccormack(double t, double tau) {
	EOSOld &eos = task.getEOS();
	ms_temp.initData(&task);
	ms_temp_temp.initData(&task);
	double E=0.; int i=0;
	double h=ms[1].x-ms[0].x;
	Vector4 L=Vector4::ZERO, R=Vector4::ZERO, D=Vector4::ZERO, V=Vector4::ZERO;
	double epsilon=0.01; 
	double gamma = eos.getGamma();
	for(i=0; i<ms.getSize(); i++) {
		Node &n = ms[i];
		Node &n_temp = ms_temp[i];
		Node &n_temp_temp = ms_temp_temp[i];
		E = n.e + 0.5*n.v*n.v;
		n.W = Vector4(n.ro, n.ro*n.v, n.ro*E, 0);
		n_temp.W = Vector4(n.ro, n.ro*n.v, n.ro*E, 0);
		n_temp_temp.W = Vector4(n.ro, n.ro*n.v, n.ro*E, 0);
	}
	Vector4 Fp=Vector4::ZERO, Fm=Vector4::ZERO, Fp_temp=Vector4::ZERO, Fm_temp=Vector4::ZERO;
	int nSize=ms.getSize(); //nSize это количество узлов сетки (и количество промежутков плюс один)

	// Предиктор
	for(i=0; i<nSize; i++) {
		Node &n = ms[i], &n_temp = ms_temp[i];
		// Прозрачная граница
		if(i==nSize-1)
			Fp = calcF(n.ro, n.v, n.p, n.e);
		else {
			Node &np = ms[i+1];
			Fp = calcF(np.ro, np.v, np.p, np.e);			
		}
		Fm = calcF(n.ro, n.v, n.p, n.e);
		n_temp.W = n.W - tau/h*(Fp-Fm);
		// Обновляем значения плотности, скорости, давления, энергии n_temp
		n_temp.ro = n_temp.W[0];
		n_temp.v  = n_temp.W[1]/n_temp.ro;
		n_temp.e  = n_temp.W[2]/n_temp.ro - 0.5*n_temp.v*n_temp.v;
		n_temp.p  = (gamma-1.) * n_temp.ro* n_temp.e;
	}
	// Корректор 
	for(i=0; i<nSize; i++) {
		Node &n = ms[i], &n_temp = ms_temp[i], &n_temp_temp = ms_temp_temp[i];
		Fp_temp = calcF(n_temp.ro, n_temp.v, n_temp.p, n_temp.e);			
		// Прозрачная граница
		if(i==0)
			Fm_temp = calcF(n_temp.ro, n_temp.v, n_temp.p, n_temp.e); 
		else {
			Node &nm_temp = ms_temp[i-1];
			Fm_temp = calcF(nm_temp.ro, nm_temp.v, nm_temp.p, nm_temp.e); 
		}
		n_temp_temp.W = 0.5*(n.W + n_temp.W) - tau/2.0/h*(Fp_temp-Fm_temp);
		// Обновляем значения плотности, скорости, давления, энергии n_temp
		n_temp_temp.ro = n_temp_temp.W[0];
		n_temp_temp.v  = n_temp_temp.W[1]/n_temp_temp.ro;
		n_temp_temp.e  = n_temp_temp.W[2]/n_temp_temp.ro - 0.5*n_temp_temp.v*n_temp_temp.v;
		n_temp_temp.p  = (gamma-1.) * n_temp_temp.ro* n_temp_temp.e;
	}
	// Искусственная вязкость
	Vector4 Qminus = Vector4::ZERO, Q=Vector4::ZERO, Qplus=Vector4::ZERO;
	double mu = 0.01;
	for(i=0; i<nSize; i++) {
		Node &n_temp_temp = ms_temp_temp[i];
		Q = n_temp_temp.W;
		// Прозрачная левая граница
		if(i==0)
			Qminus = ms_temp_temp[0].W;
		else 
			Qminus = ms_temp_temp[i-1].W;
		// Прозрачная правая граница
		if(i==nSize-1)
			Qplus = ms_temp_temp[nSize-1].W;
		else
			Qplus = ms_temp_temp[i+1].W;
		ms_temp_temp[i].W_temp = mu * (Qminus - 2.*Q + Qplus);
	}
	// Тест законов сохранения
	// Тестируем предиктор
	for(i=0; i<nSize; i++) {
		L = 1./tau * (ms_temp[i].W - ms[i].W);
		if(i!=nSize-1)
			Fp = Vector4(ms[i+1].ro*ms[i+1].v,
						 ms[i+1].p + ms[i+1].ro*ms[i+1].v*ms[i+1].v,
						 ms[i+1].v*(ms[i+1].p + ms[i+1].ro*(ms[i+1].e + 0.5*ms[i+1].v*ms[i+1].v)),
						 0.);
		else
			Fp = Vector4(ms[i].ro*ms[i].v,
						 ms[i].p + ms[i].ro*ms[i].v*ms[i].v,
						 ms[i].v*(ms[i].p + ms[i].ro*(ms[i].e + 0.5*ms[i].v*ms[i].v)),
						 0.);
		Fm = Vector4(ms[i].ro*ms[i].v,
					 ms[i].p + ms[i].ro*ms[i].v*ms[i].v,
					 ms[i].v*(ms[i].p+ms[i].ro*(ms[i].e+0.5*ms[i].v*ms[i].v)),
					 0.);
		R = -1./h * (Fp - Fm);
		D = L-R;
		if(fabs(D[0]/L[0]) > epsilon && fabs(D[0]) > epsilon)
			cout << "i=" << i << ": Conservation law on predictor step violated in eq.#1" << endl;
		if(fabs(D[1]/L[1]) > epsilon && fabs(D[1]) > epsilon)
			cout << "i=" << i << ": Conservation law on predictor step violated in eq.#2" << endl;
		if(fabs(D[2]/L[2]) > epsilon && fabs(D[2]) > epsilon)
			cout << "i=" << i << ": Conservation law on predictor step violated in eq.#3" << endl;
	}
	// Тестируем корректор и искусственную вязкость
	for(i=0; i<nSize; i++) {

		// Fp 				x	0.82565114466347900	double
		//y	1.6961706557372802	double
		//z	3.1537303898249780	double


		// Fm 		x	0.75000000000000000	double
		// y	1.5625000000000002	double
		// z	2.8359375000000013	double

		// n.W 		x	1.0000000000000000	double
		// y	0.75000000000000000	double
		// z	2.7812500000000013	double


		// n_temp.W x	1.0387954588017840	double
		// y	0.82565114466347900	double
		// z	2.9279453285942476	double

		// n_temp_temp.W 		x	1.0174411154901588	double
		// y	0.78436836271996413	double
		// z	2.8463783836539878	double


		// V  		x	0.00017441115490158811	double
		// y	0.00034368362719964130	double
		// z	0.00065128383653986429	double



		/////////////
		
		Fp_temp = Vector4(ms_temp[i].ro*ms_temp[i].v,
						  ms_temp[i].p + ms_temp[i].ro*ms_temp[i].v*ms_temp[i].v,
						  ms_temp[i].v*(ms_temp[i].p + ms_temp[i].ro*(ms_temp[i].e + 0.5*ms_temp[i].v*ms_temp[i].v)),
						  0.);
		if(i!=0)
			Fm_temp = Vector4(ms_temp[i-1].ro*ms_temp[i-1].v,
							  ms_temp[i-1].p + ms_temp[i-1].ro*ms_temp[i-1].v*ms_temp[i-1].v,
							  ms_temp[i-1].v*(ms_temp[i-1].p + ms_temp[i-1].ro*(ms_temp[i-1].e+0.5*ms_temp[i-1].v*ms_temp[i-1].v)),
							  0.);
		else
			Fm_temp = Vector4(ms_temp[i].ro*ms_temp[i].v,
							  ms_temp[i].p + ms_temp[i].ro*ms_temp[i].v*ms_temp[i].v,
							  ms_temp[i].v*(ms_temp[i].p + ms_temp[i].ro*(ms_temp[i].e + 0.5*ms_temp[i].v*ms_temp[i].v)),
							  0.);
		R = - 1./h/2. * (Fp_temp - Fm_temp);
		// Искусственная вязкость
		if(i==0) Qminus = ms_temp_temp[0].W; else Qminus = ms_temp_temp[i-1].W;
		if(i==nSize-1) Qplus = ms_temp_temp[nSize-1].W; else Qplus = ms_temp_temp[i+1].W;
		Q = ms_temp_temp[i].W;
		V = mu*(Qminus - 2.*Q + Qplus);
		L = 1./tau * ( ms_temp_temp[i].W + V - ms_temp_temp[i].W_temp - 0.5*(ms[i].W + ms_temp[i].W) );
		D = L-R;
		if(fabs(D[0]/L[0]) > epsilon && fabs(D[0]) > epsilon)
			cout << "i=" << i << ": Conservation law on corrector step violated in eq.#1" << endl;
		if(fabs(D[1]/L[1]) > epsilon && fabs(D[1]) > epsilon)
			cout << "i=" << i << ": Conservation law on corrector step violated in eq.#2" << endl;
		if(fabs(D[2]/L[2]) > epsilon && fabs(D[2]) > epsilon)
			cout << "i=" << i << ": Conservation law on corrector step violated in eq.#3" << endl;
	}
	// Присваиваем текущему временнОму слою рассчитанные значения
	for(i=0; i<nSize; i++) 
	{
		// Добавляем искусственную вязкость
		ms_temp_temp[i].W += ms_temp_temp[i].W_temp;
		// Обновляем консервативные переменные
		Node& n=ms[i], &n_temp_temp=ms_temp_temp[i]; 
		n.W[0] = n_temp_temp.W[0];
		n.W[1] = n_temp_temp.W[1];
		n.W[2] = n_temp_temp.W[2];
		// Обновляем все переменные
		n.ro = n.W[0];
		n.v  = n.W[1]/n.ro;
		n.e  = n.W[2]/n.ro - 0.5*n.v*n.v;
		n.ti = eos.getti(n.ro, n.e);
		n.p  = eos.getpi(n.ro, n.ti);
		n.C  = eos.getC(n.ro, n.ti, n.te);
	}
	cout << "calcHydroStageMccormack(): done!" << endl;
}

void CSolver::calcHydroStageRoe(double t, double tau) {
	EOSOld &eos = task.getEOS();
	double E=0.; 
	int i=0, nSize = ms.getSize();
	double h=ms[1].x-ms[0].x;
	Vector4 L = Vector4::ZERO, R = Vector4::ZERO, D = Vector4::ZERO, V = Vector4::ZERO;
	double gamma = eos.getGamma();
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
			n.F = calcRoeFlux(roLB, roLB*vLB, roLB*ELB, n.W[0], n.W[1], n.W[2]);
		else
			n.F = calcRoeFlux(ms[i-1].W[0], ms[i-1].W[1], ms[i-1].W[2], n.W[0], n.W[1], n.W[2]);
		double qqq =0.;
	}
	ms[nSize].F = calcRoeFlux(ms[nSize-1].W[0], ms[nSize-1].W[1], ms[nSize-1].W[2], roRB, roRB*vRB, roRB*ERB);
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

Vector4 CSolver::calcGPSFlux(double roL, double rouL, double roEL, double roR, double rouR, double roER) {
	EOSOld& eos = task.getEOS(); 
	double gamma = eos.getGamma(), _ro=0., _u=0., _p=0., _e=0., _E=0., _sigma = 0., sigmaL = 0., sigmaR = 0.;
	double uL = rouL/roL, uR = rouR/roR, EL = roEL/roL, ER = roER/roR, eL = EL - 0.5*uL*uL, eR = ER - 0.5*uR*uR, 
		   pL = (gamma-1.)*roL*eL, pR = (gamma-1.)*roR*eR, cL = sqrt(gamma*pL/roL), cR = sqrt(gamma*pR/roR);
	sigmaL = pL/pow(roL, gamma);
	sigmaR = pR/pow(roR, gamma);
	if(uL > cL) 
		return Vector4(roL*uL, pL+roL*uL*uL, uL*(pL+roL*EL), 0.);
	else if (uR < -cR)
		return Vector4(roR*uR, pR+roR*uR*uR, uR*(pR+roR*ER), 0.);
	_p = (pL / roL / cL + pR / roR / cR + uL - uR) / (1. / roL / cL + 1 / roR / cR);
	_u = (roL * cL * uL + roR * cR * uR + pL - pR) / (roL * cL + roR * cR);
	if(_u > 0.)
		_sigma = sigmaL; 
	else 
		_sigma = sigmaR;
	_ro = pow(_p/_sigma, 1./gamma);
	_e  = _p/(gamma-1.)/_ro;
	_E  = _e + 0.5*_u*_u;
	return Vector4(_ro*_u, _p+_ro*_u*_u, _u*(_p+_ro*_E), 0.);
}


Vector4 CSolver::calcRoeFlux(double roL, double rouL, double roEL, double roR, double rouR, double roER) {
	EOSOld &eos = task.getEOS(); double gamma = eos.getGamma();
	double uL = rouL/roL, eL = roEL/roL - 0.5*uL*uL, uR = rouR/roR, eR = roER/roR - 0.5*uR*uR;
	double pL = (gamma-1.)*roL*eL, pR = (gamma-1.)*roR*eR, cL = sqrt(gamma*pL/roL), cR = sqrt(gamma*pR/roR);
	double EL = roEL/roL, ER = roER/roR, HL = (roEL + pL)/roL, HR = (roER + pR)/roR;
	double sqroL = sqrt(roL), sqroR = sqrt(roR);
	double uAv = (sqroL*uL + sqroR*uR)/(sqroL+sqroR), HAv = (sqroL*HL + sqroR*HR)/(sqroL+sqroR), cAv = sqrt((HAv-.5*uAv*uAv)*(gamma-1.));
	Vector4 Lambda = Vector4(uAv - cAv, uAv, uAv + cAv, 0.),
		        K0 = Vector4(1., uAv-cAv, HAv-uAv*cAv, 0.),
			    K1 = Vector4(1., uAv,     uAv*uAv/2.,  0.),
				K2 = Vector4(1., uAv+cAv, HAv+uAv*cAv, 0.);
	Vector4 Alpha = Vector4::ZERO;
	Alpha[1] = (gamma-1.)/cAv/cAv * ((HAv-uAv*uAv)*(roR-roL) + uAv*(rouR-rouL) - (roER-roEL));
	Alpha[0] = 1./2./cAv * ((roR-roL)*(uAv+cAv) - (rouR-rouL) - cAv*Alpha[1]);
	Alpha[2] = (roR-roL) - Alpha[0] - Alpha[1];
	Vector4 FL = Vector4(roL*uL, pL+roL*uL*uL, uL*(pL+roEL), 0.);
	Vector4 FR = Vector4(roR*uR, pR+roR*uR*uR, uR*(pR+roER), 0.);
	// Hyman's entropy fix with Roe average v and c in the Star region (Toro)
	/*double lambdaL=0., lambdaR=0.;
	double roStar=0., uStar=0., pStar=0., EStar=0., cStar=0.;
	/// DEBUG ///
	if(uAv-cAv > 0.) {
		double qq=0.;
	}
	////////////
    roStar = roL+Alpha[0]; uStar = (roL*uL+Alpha[0]*(uAv-cAv))/(roL+Alpha[0]);
	pStar = (gamma-1.)*(EL+Alpha[0]*(HAv-uAv*cAv) - 0.5*roStar*uStar*uStar);
	cStar = sqrt(gamma*pStar/roStar);
	lambdaL = uL-cL; lambdaR = uStar-cStar;
	if( lambdaL < 0. && lambdaR > 0.) {		
		cout << "Entropy fix!" << endl;
		Lambda[0] = lambdaL*(lambdaR-Lambda[0])/(lambdaR-lambdaL);
	}
	roStar = roR-Alpha[2]; uStar = (roR*uR-Alpha[2]*(uAv+cAv))/(roR-Alpha[2]);
	pStar = (gamma-1.)*(ER-Alpha[2]*(HAv+uAv*cAv) - 0.5*roStar*uStar*uStar);
	cStar = sqrt(gamma*pStar/roStar);
	lambdaL = uStar+cStar; lambdaR = uR+cR;
	if( lambdaL < 0. && lambdaR > 0.) {
		cout << "Entropy fix!" << endl;
		Lambda[2] = lambdaR*(Lambda[2]-lambdaL)/(lambdaR-lambdaL);
	}*/
	// Hyman entropy fix, Katate Matatsuka's realization
	/*
	!Absolute values of the wave speeds (Eigenvalues)
   ws(1) = abs(v-a)
   ws(2) = abs(v  )
   ws(3) = abs(v+a)

!Modified wave speeds for nonlinear fields (the so-called entropy fix, which
!is often implemented to remove non-physical expansion shocks).
!There are various ways to implement the entropy fix. This is just one
!example. Try turn this off. The solution may be more accurate.
   Da = max(zero, four*((vR-aR)-(vL-aL)) )
   if (ws(1) < half*Da) ws(1) = ws(1)*ws(1)/Da + quarter*Da
   Da = max(zero, four*((vR+aR)-(vL+aL)) )
   if (ws(3) < half*Da) ws(3) = ws(3)*ws(3)/Da + quarter*Da
	*/
	double delta=0.;
	delta = max(0., 4.*((uR-cR)-(uL-cL)) );
    if (fabs(Lambda[0]) < 0.5*delta) Lambda[0] = Lambda[0]*Lambda[0]/delta + .25*delta;
    delta = max(0., 4.*((uR+cR)-(uL+cL)) );
    if (fabs(Lambda[2]) < 0.5*delta) Lambda[2] = Lambda[2]*Lambda[2]/delta + .25*delta;

	Vector4 FRoe = 0.5*(FL + FR) - 0.5*(Alpha[0]*fabs(Lambda[0])*K0 + 
		                                Alpha[1]*fabs(Lambda[1])*K1 + 
										Alpha[2]*fabs(Lambda[2])*K2);
	return FRoe;
}

Vector4 CSolver::calcGodunovFlux(double roL, double rouL, double roEL, double roR, double rouR, double roER) {
	EOSOld &eos=task.getEOS(); double gamma=eos.getGamma();	
	double uL = rouL/roL, uR = rouR/roR;
	double pL = (gamma-1.)*(roEL - .5*roL*uL*uL), pR = (gamma-1.)*(roER - .5*roR*uR*uR);
	CVectorPrimitive res = calcRPAnalyticalSolution(roL, uL, pL, roR, uR, pR, 0., .01);
	double E = res.p/(gamma-1.)/res.ro + .5*res.v*res.v;
	Vector4 FGodunov = Vector4(res.ro*res.v, res.ro*res.v*res.v + res.p, res.v*(res.ro*E+res.p), 0.);
	return FGodunov;
}

double CSolver::getdx() {
	return ms[1].x-ms[0].x;
}

// ТРЭШ -- убираем мусор пока под комментарии

/* Это не вполне верная реализация Годуновской схемы -- вместо решения
   задачи о распаде разрыва между ячейками для подсчета потоков, я 
   тупо беру полусуммы. Тем не менее, она почему-то как-то сходится 
   к решению. :) 
  

void CSolver::calcHydroStageGodunov(double t, double tau)
{
	MatterState ms_temp;
    ms_temp.initData(&task);
	EOS &eos = task.getEOS();
	int nSize = ms.getSize();

	double      E_temp = 0.0,
				 E_cur = 0.0,
			 v_av_temp = 0.0,
			  v_av_cur = 0.0,
		   pi_plus_cur = 0.0,
		  pi_minus_cur = 0.0;


	// W = ( 1/ro v E )T

	// ВНИМАНИЕ!!! Граничное условие: P=0, как справа, так и слева!!!

	// Первый компонент (1/ro)

	for(int i=0; i<nSize; i++)
	{
		ms_temp[i].ro = 1.0/
					   (1.0/ms[i].ro + (tau/ms[i].dm)*(ms[i+1].v-ms[i].v));
	}

	// Второй компонент (v) и сразу же x

	for(int i=0; i<nSize; i++)
	{	
		if(i==0) 	
		{
			ms_temp[i].v = (tau/ms[i].dm)*(0.0 - ms[i].p) + ms[i].v;
		}			
		else
		{
			ms_temp[i].v = (tau/ms[i].dm)*(ms[i-1].p - ms[i].p) + ms[i].v;
		}
	
		ms_temp[i].x = ms[i].x + ms_temp[i].v*tau;
	}

	i=nSize;

	ms_temp[i].v = (tau/ms[i].dm)*(ms[i-1].p - 0.0) + ms[i].v;
	ms_temp[i].x = ms[i].x + ms_temp[i].v*tau;

	//////////DEBUG/////////////
	
	double v0	   = ms_temp[0].v;
	double v0_plus = 0.5*(ms_temp[0].v + ms_temp[1].v);
	double v1	   = ms_temp[1].v;
	double v1_plus = 0.5*(ms_temp[1].v + ms_temp[2].v);
	double p0	   = ms[0].p;
	double p1	   = ms[1].p;
	////////////////////////////

	///////////////
		
	double v_bef_last = ms_temp[nSize-1].v;
	double v_last	  = ms_temp[nSize].v;

	///////////////

	// Третий компонент (E)
	
	for(int i=0; i<nSize; i++)
	{	
		v_av_cur = 0.5 * (ms[i].v + ms[i+1].v);
		E_cur = ms[i].e + v_av_cur*v_av_cur/2.0;
		
		if(i==0)
		{
			 pi_plus_cur = 0.5 * (ms[i].pi + ms[i+1].pi);
			pi_minus_cur = 0.5 * (     0.0 +   ms[i].pi);
		}
		else if(i==nSize-1)
		{
			 pi_plus_cur = 0.5 * (  ms[i].pi + 0.0);
			pi_minus_cur = 0.5 * (ms[i-1].pi + ms[i].pi);
		}
		else
		{
			 pi_plus_cur = 0.5 * (ms[i].pi + ms[i+1].pi);
			 pi_minus_cur = 0.5 * (ms[i-1].pi + ms[i].pi);
		}
	
		E_temp = (tau/ms[i].dm)*(pi_minus_cur*ms[i].v -
								pi_plus_cur*ms[i+1].v) + E_cur;

		v_av_temp = 0.5 * (ms_temp[i].v + ms_temp[i+1].v);		
		ms_temp[i].ei = E_temp - v_av_temp*v_av_temp/2.0;
	}

	// Остальные величины
	
	for(int i=0; i<nSize; i++)
	{
		ms_temp[i].ti = eos.getti(ms_temp[i].ro, ms_temp[i].ei);
		ms_temp[i].te = ms[i].te;	

		ms_temp[i].pi = eos.getpi(ms_temp[i].ro, ms_temp[i].ti);

		ms_temp[i].ee = eos.getee(ms_temp[i].ro, ms_temp[i].ti, ms_temp[i].te);
		ms_temp[i].pe = eos.getpe(ms_temp[i].ro, ms_temp[i].ti, ms_temp[i].te);

		ms_temp[i].e  = eos.gete(ms_temp[i].ro, ms_temp[i].ti, ms_temp[i].te);
		ms_temp[i].p  = eos.getp(ms_temp[i].ro, ms_temp[i].ti, ms_temp[i].te);
	}
	
	// Копируем темповые массивы в оригинальные

	for(int i=0; i<nSize; i++)
	{
		 ms[i].x = ms_temp[i].x;
		 ms[i].v = ms_temp[i].v;
		ms[i].ro = ms_temp[i].ro;
		
		 ms[i].e = ms_temp[i].e;
		ms[i].ee = ms_temp[i].ee;
		ms[i].ei = ms_temp[i].ei;

		 ms[i].p = ms_temp[i].p;
		ms[i].pi = ms_temp[i].pi;
		ms[i].pe = ms_temp[i].pe;

		ms[i].te = ms_temp[i].te;
		ms[i].ti = ms_temp[i].ti;

		 ms[i].C     = eos.getC(ms[i].ro, ms[i].ti, ms[i].te);
		 ms[i].kappa = eos.getkappa(ms[i].ro, ms[i].ti, ms[i].te);
		 ms[i].Alphaei		= eos.getAlpha(ms[i].ro, ms[i].ti, ms[i].te);
	}

	ms[nSize].x = ms_temp[nSize].x;
	ms[nSize].v = ms_temp[nSize].v;

}
*/

/* Это реализация "обычной" схемы Самарского, той, которая 
приведена для примера в книге Самарского и Попова. Плюс, 
искусственная вязкость почти наверняка работает неправильно,
потому что в процессе программирования этой схемы я с трудом
представлял себе, что это такое. :) 
  
	
void CSolver::calcHydroStageSamarskii(double t, double tau)
{	
	int i=0;
	int counter=0;
	int itNumFull=0, itNumIon=0;
	
	double e_prev=0.0, ei_prev=0.0;
	double e=0.0, dm=0.0, e_next=0.0, p_next=0.0, ro_next=0.0, 
		   ti_next=0.0, te_next=0.0;
	double Alphaei=0, dv=0;

	MatterState ms_temp;

    ms_temp.initData(&task);

	EOS &eos = task.getEOS();
	int nSize = ms.getSize();

// Первое уравнение ((2.7) из книжки Самарского и Попова)
// Граничное условие с обеих сторон -- нулевое давление (вакуум)

	for(int i=0; i<nSize-1; i++)
	{	
		if(i==0) 	
		{
			ms_temp[i].v = (tau/ms[i].dm)*(0.0 - ms[i].p) + ms[i].v;
		}				
		else
		{
			ms_temp[i].v = (tau/ms[i].dm)*(ms[i-1].p - ms[i].p) + ms[i].v;
		}
		ms_temp[i].x = ms[i].x + ms_temp[i].v*tau;
	}

	/////////////////DEBUG//////////////////
	//ms_temp[nSize-1].v = (tau/ms[nSize-1].dm)*(ms[nSize-2].p-0.0)+ms[nSize-1].v;	
	// Живет, пока не автоматизируем граничные условия
		ms_temp[nSize-1].v = 0.0;
	///////////////////////////////
	
	ms_temp[nSize-1].x = ms[nSize-1].x + ms_temp[nSize-1].v*tau;

//Уравнение (2.9) из книжки Самарского и Попова (массив v_temp уже вычислен полностью)

	for(int i=0; i<nSize-1; i++)
	{
		ms_temp[i].ro = 1.0/
					   (1.0/ms[i].ro + (tau/ms[i].dm)*(ms_temp[i+1].v-ms_temp[i].v));
		
		////////////DEBUG////////////////
		//Исключительно для борьбы с антисанитарией в схеме.
		
	//	if (ms_temp[i].ro > 3000.0)
	//		ms_temp[i].ro = 3000.0;


		//////////////////////////////////


	}

// Энергия	

	for(int i=0; i<nSize-1; i++)
	{

// Закон сохранения ионной энергии. В однотемпературном случае 
// является законом сохранения полной внутренней энергии (2.11)
// (без кинетической)


		///////////////////////DEBUG///////////////////

	//	if(i==595)
	//		dumpToFile(t);

		////////////////////////////////////////////////


		ms_temp[i].ti = solveti(tau, i, ms_temp[i].ro, ms_temp[i+1].v-ms_temp[i].v);
		ms_temp[i].ei = eos.getei(ms_temp[i].ro, ms_temp[i].ti); 	
		ms_temp[i].pi = eos.getpi(ms_temp[i].ro, ms_temp[i].ti);

// Закон сохранения для электронной энергии (равносильно закону 
// сохранения для полной энергии (2.11), но в двухтемпературном
// случае)

		ms_temp[i].te = solvete(t, tau, i, ms_temp[i].ro, ms_temp[i+1].v-ms_temp[i].v,
		  					            ms_temp[i].ti);
	
		ms_temp[i].ee = eos.getee(ms_temp[i].ro, ms_temp[i].ti, ms_temp[i].te);
		ms_temp[i].pe = eos.getpe(ms_temp[i].ro, ms_temp[i].ti, ms_temp[i].te);

		ms_temp[i].e  = eos.gete(ms_temp[i].ro, ms_temp[i].ti, ms_temp[i].te);
		ms_temp[i].p  = eos.getp(ms_temp[i].ro, ms_temp[i].ti, ms_temp[i].te);
	}

// Искусственная вязкость. Протестирована на ударной волне у А.М.

	if(task.getViscFlag())
	{
		for(int i=0; i<nSize-2; i++)
		{
			double dv = ms_temp[i+1].v - ms_temp[i].v;
			Alphaei=4.0*ms_temp[i+1].ro*dv*dv;
			if (dv < 0)
			{
				ms_temp[i].pi += Alphaei;
				ms_temp[i].p  += Alphaei;
			}
		}
	}

	for(int i=0; i<nSize-1; i++)
	{
		 ms[i].x = ms_temp[i].x;
		 ms[i].v = ms_temp[i].v;
		ms[i].ro = ms_temp[i].ro;
		
		 ms[i].e = ms_temp[i].e;
		ms[i].ee = ms_temp[i].ee;
		ms[i].ei = ms_temp[i].ei;

		 ms[i].p = ms_temp[i].p;
		ms[i].pi = ms_temp[i].pi;
		ms[i].pe = ms_temp[i].pe;

////////////!!!DEBUG!!!////////////////////////
		
		if(ms[i].pe<0)
		{
			double pe = ms[i].pe;
			double p  = ms[i].p;
			double ti = ms[i].ti;
			double te = ms[i].te;
			double ro = ms[i].ro;
			
			double pe2= eos.getpe(ro, ti, te);

		}

///////////////////////////////////////////////



		ms[i].te = ms_temp[i].te;
		ms[i].ti = ms_temp[i].ti;

		 ms[i].C     = eos.getC(ms[i].ro, ms[i].ti, ms[i].te);
		 ms[i].kappa = eos.getkappa(ms[i].ro, ms[i].ti, ms[i].te);
		 ms[i].Alphaei		= eos.getAlpha(ms[i].ro, ms[i].ti, ms[i].te);
	}

	ms[nSize-1].x = ms_temp[nSize-1].x;
	ms[nSize-1].v = ms_temp[nSize-1].v;	

}

*/


/* Это то же самое, что и calcHydroStageSamarskii(), но добавлен 
источниковый член для учета поглощения лазерного излучения. 
Искусственная вязкость тоже не работает. 



void CSolver::calcHydroStageSamarskiiSource(double t, double tau)
{
	int i=0;
	int counter=0;
	int itNumFull=0, itNumIon=0;
	
	double e_prev=0.0, ei_prev=0.0;
	double e=0.0, dm=0.0, e_next=0.0, p_next=0.0, ro_next=0.0, 
		   ti_next=0.0, te_next=0.0;
	double Alphaei=0, dv=0;

	MatterState ms_temp;

    ms_temp.initData(&task);

	EOS &eos = task.getEOS();
	int nSize = ms.getSize();

// Первое уравнение ((2.7) из книжки Самарского и Попова)
// Граничное условие с обеих сторон -- нулевое давление (вакуум)

	for(int i=0; i<nSize-1; i++)
	{	
		if(i==0) 	
		{
			ms_temp[i].v = (tau/ms[i].dm)*(0.0 - ms[i].p) + ms[i].v;
		}				
		else
		{
			ms_temp[i].v = (tau/ms[i].dm)*(ms[i-1].p - ms[i].p) + ms[i].v;
		}
		ms_temp[i].x = ms[i].x + ms_temp[i].v*tau;
	}


	/////////////////DEBUG//////////////////
	//ms_temp[nSize-1].v = (tau/ms[nSize-1].dm)*(ms[nSize-2].p-0.0)+ms[nSize-1].v;	
	// Живет, пока не автоматизируем граничные условия
		ms_temp[nSize-1].v = 0.0;
	///////////////////////////////
	
	ms_temp[nSize-1].x = ms[nSize-1].x + ms_temp[nSize-1].v*tau;

//Уравнение (2.9) из книжки Самарского и Попова (массив v_temp уже вычислен полностью)

	for(int i=0; i<nSize-1; i++)
	{
		ms_temp[i].ro = 1.0/
					   (1.0/ms[i].ro + (tau/ms[i].dm)*(ms_temp[i+1].v-ms_temp[i].v));

	}

// Энергия	

	for(int i=0; i<nSize-1; i++)
	{

// Закон сохранения ионной энергии. В однотемпературном случае 
// является законом сохранения полной внутренней энергии (2.11)
// (без кинетической)

		ms_temp[i].ti = solveti(tau, i, ms_temp[i].ro, ms_temp[i+1].v-ms_temp[i].v);
		ms_temp[i].ei = eos.getei(ms_temp[i].ro, ms_temp[i].ti); 	
		ms_temp[i].pi = eos.getpi(ms_temp[i].ro, ms_temp[i].ti);

// Закон сохранения для электронной энергии (равносильно закону 
// сохранения для полной энергии (2.11), но в двухтемпературном
// случае)

		ms_temp[i].te = solveteSource(t, tau, i, ms_temp[i].ro, ms_temp[i+1].v-ms_temp[i].v,
		  					            ms_temp[i].ti);
	
		ms_temp[i].ee = eos.getee(ms_temp[i].ro, ms_temp[i].ti, ms_temp[i].te);
		ms_temp[i].pe = eos.getpe(ms_temp[i].ro, ms_temp[i].ti, ms_temp[i].te);

		ms_temp[i].e  = eos.gete(ms_temp[i].ro, ms_temp[i].ti, ms_temp[i].te);
		ms_temp[i].p  = eos.getp(ms_temp[i].ro, ms_temp[i].ti, ms_temp[i].te);
	}

// Искусственная вязкость. Протестирована на ударной волне у А.М.

	if(task.getViscFlag())
	{
		for(int i=0; i<nSize-2; i++)
		{
			double dv = ms_temp[i+1].v - ms_temp[i].v;
			Alphaei=4.0*ms_temp[i+1].ro*dv*dv;
			if (dv < 0)
			{
				ms_temp[i].pi += Alphaei;
				ms_temp[i].p  += Alphaei;
			}
		}
	}


	for(int i=0; i<nSize-1; i++)
	{
		 ms[i].x = ms_temp[i].x;
		 ms[i].v = ms_temp[i].v;
		ms[i].ro = ms_temp[i].ro;
		
		 ms[i].e = ms_temp[i].e;
		ms[i].ee = ms_temp[i].ee;
		ms[i].ei = ms_temp[i].ei;

		 ms[i].p = ms_temp[i].p;
		ms[i].pi = ms_temp[i].pi;
		ms[i].pe = ms_temp[i].pe;

		ms[i].te = ms_temp[i].te;
		ms[i].ti = ms_temp[i].ti;

		 ms[i].C     = eos.getC(ms[i].ro, ms[i].ti, ms[i].te);
		 ms[i].kappa = eos.getkappa(ms[i].ro, ms[i].ti, ms[i].te);
		 ms[i].Alphaei		= eos.getAlpha(ms[i].ro, ms[i].ti, ms[i].te);
	}

	ms[nSize-1].x = ms_temp[nSize-1].x;
	ms[nSize-1].v = ms_temp[nSize-1].v;	
}

*/


/* void CSolver::testGodunovScheme(char* inputFileName)
{
	task.load(inputFileName);
	ms.initData(&task);
	EOS &eos = task.getEOS();

	double t = 0.0;
	double tau = 0.0;

	for(int i=0; i<10000; i++)
	{
		tau = calcTimeStep(t);
		tau = 0.1e-11;
		calcHydroStageGodunovQuasiVector(t, tau);
		//calcHydroStageGodunov(t, tau);
		t+=tau;
		dumpToFile(t);	
	}

}

void CSolver::testSaveSolution(char *inputFileName)
{
	task.load(inputFileName);
	ms.initData(&task);
	EOS &eos = task.getEOS();

	double t = 0.0;
	double tau = 0.0;

	for(int i=0; i<50; i++)
	{
		tau = calcTimeStep(t);
		calcHydroStageSamarskii(t, tau);
		calcHeatStage(tau);
		calcExchangeStage(tau);
		
		t+=tau;
		dumpToFile(t);	
	}

////////Дописать сюда сохранение и восстановление//////////

	for(int i=0; i<50; i++)
	{
		tau = calcTimeStep(t);
		calcHydroStageSamarskii(t, tau);
		calcHeatStage(tau);
		calcExchangeStage(tau);
		
		t+=tau;
		dumpToFile(t);	
	}

}

  */

void CSolver::testEOSControlNumbers(double _ro, double _ti, double _te) {
	//double _ro = 19301., _ti = 1000., _te = 1000.;
	EOSOld& eos = task.getEOS();
	cout << endl <<"Gold EOS test" << endl;
	cout << "ro = " << _ro << " kg/m3, ti = " << _ti <<" K, te = " << _te << " K" << endl;
	double _ei = eos.getei(_ro, _ti); double _eiToCGSMul = 1.e-6;
	double _pi = eos.getpi(_ro, _ti); double _piToGPaMul = 1.e-9;
	double _ee = eos.getee(_ro, _ti, _te); double _eeToCGSMul = _eiToCGSMul;
	double _pe = eos.getpe(_ro, _ti, _te); double _peToGPaMul = _piToGPaMul;
	double _ci = eos.getci(_ro, _ti)/_ro; double _ciToCGSMul = 1.e-3;
	double _ce = eos.getce(_ro, _te)/_ro; double _ceToCGSMul = 1.e-6;
	double _alpha = eos.getAlpha(_ro, _ti, _te); double _alphaToCGSMul = 1.e-9;
	double _kappa = eos.getkappa(_ro, _ti, _te); double _kappaToCGSMul = 1.e-5;
	cout << setw(20) << "ei(ro, ti) = "        << setw(12) << _ei     << setw(12) << " J/kg = " << _ei*_eiToCGSMul << " kJ/g" << endl;
	cout << setw(20) << "pi(ro, ti) = "        << setw(12) << _pi     << setw(12) << " Pa = " << _pi*_piToGPaMul << " GPa" << endl;
	cout << setw(20) << "ci(ro, ti, te) = "    << setw(12) << _ci     << setw(12) << " J/kg/K = " << _ci*_ciToCGSMul << " kJ/g/kK" << endl;
	cout << setw(20) << "ee(ro, ti, te) = "    << setw(12) << _ee     << setw(12) << " J/kg = " << _ee*_eeToCGSMul << " kJ/g" << endl;
	cout << setw(20) << "pe(ro, ti, te) = "    << setw(12) << _pe     << setw(12) << " Pa = " << _pe*_peToGPaMul << " GPa" << endl;
	cout << setw(20) << "ce(ro, ti, te) = "    << setw(12) << _ce     << setw(12) << " J/kg/K = " << _ce*_ceToCGSMul << " kJ/g/K" << endl;
	cout << setw(20) << "alpha(ro, ti, te) = " << setw(12) << _alpha  << setw(12) << " J/s/m3/K = " << _alpha*_alphaToCGSMul << " kW/cm3/K" << endl;
	cout << setw(20) << "kappa(ro, ti, te) = " << setw(12) << _kappa  << setw(12) << " J/s/m/K = " << _kappa*_kappaToCGSMul << " kW/sm/K" << endl;
	cout << "Test finished!" << endl;
	return;
}
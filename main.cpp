#define _CRT_SECURE_NO_WARNINGS

//#include <vld.h>

#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>

#include "defines.h"
#include "solver.h"
#include "C1DProblem.h"
#include "F1DSimulation.h"
#include "C1DBCs.h"
#include "F1DReconstruction.h"
#include "eos/FEOSMGLiF.h"

char* INPUT_FOLDER = new char[_MAX_PATH];  //"calc/";
char* OUTPUT_FOLDER = new char[_MAX_PATH]; //"calc/output/";


int main(int argc, char *argv[]) {
	if(argc > 1) {
		strcpy(INPUT_FOLDER, argv[1]);
		strcpy(OUTPUT_FOLDER, INPUT_FOLDER);
		strcat(OUTPUT_FOLDER, "output/");
	} else {
		strcpy(INPUT_FOLDER, "calc/");
		strcpy(OUTPUT_FOLDER, "calc/output/");
	}
	cout << "FSLA1D: hydrocode for numerical simulations in 1D-geometry v.0.1." << endl;
	cout << "Author: Vadim V. Shepelev, ICAD RAS, e-mail: vadim.v.shepelev@gmail.com" << endl;
	cout << "=======================================================================" << endl;
	string outputDir = string("output");
	// Uncomment for NB EOS test problem
	// CEOSMieGruneisenAl eos = CEOSMieGruneisenAl();
	//C1DProblem pr = prNBtest;
	// Uncomment for Toro #1 test problem with ideal EOS
	//CEOSIdeal eos = CEOSIdeal(3.9);
	//CEOSIdeal eosAl = CEOSIdeal(3.9);
	//C1DProblem pr = prLaserVTAlIdealTest1;
	//C1DField *fldptr = new C1DField(pr);
	// Uncomment for exact solver based Godunov-type method
	//CExactRiemannSolver ex;
	//C1DGodunovTypeMethod mtd = C1DGodunovTypeMethod(ex);

	/*
	CSolver *s = new CSolver;
	s->goEuler("task-toro-5.txt"); // для Эйлеровых задач
	delete s;
	*/

	// Uncomment for Lagrange 1D code in old architecture /*
	
	/*FEOSMGLiF eos = FEOSMGLiF();
	double e0 = eos.gete(2640., 0.);
	double c0 = eos.getc(2640., 0.);
	double e1 = eos.gete(6000., 1.e9); // -8.35e7 done		 e_cold = 6.81e7 done
	double c1 = eos.getc(6000., 1.e9); // 4.484e4
	double c2 = eos.getc(3000., 1.e9);*/
	// pCold = 6.46e11 -- done
	// p_e = 4.26e3 -- done 
	// pCold_tau = -1.37e17 = rho0*dpc_dx = -8254662065539,5137 * 6000 = -49 527 972 393 237 082,2	
	// eCold_tau = -5.45e16 = rho0*dec_dx = -991697679,97802985 * 6000 = -5 950 186 079 868,1791


	/*
	CSolver *s = new CSolver;
	s->go("task-LiF.txt"); // для Эйлеровых задач
	delete s;*/

	// Uncomment for Lagrange 1D code in new architecture
	/*FEOSMGLiF eos;
	C1DProblem pr =  prLiF; //prTestLagrange1D;
	C1DFieldPrimitive *fldptr = new C1DFieldPrimitive(pr);
	C1DMethodSamarskii mtd;
	double _dtt[] = {pr.tmin, 
		             .05e-9, 
		             .1e-9, .2e-9, .3e-9, .4e-9, .5e-9, .6e-9, .7e-9, .8e-9, .9e-9, 
		             1.e-9, 2.e-9, 3.e-9, 3.128e-9, 4.e-9, 5.e-9, 6.e-9, 7.e-9, 8.e-9, 9.e-9,
				 	 pr.tmax};
	vector<double> dtt = vector<double>(_dtt, _dtt+sizeof(_dtt)/sizeof(double));
	COutput outp = COutput(pr, outputDir, dtt);
	F1DSimulationLagrange sim = F1DSimulationLagrange(pr, eos, *fldptr, mtd, outp);
	sim.run();
	delete fldptr;*/

	FEOSMGLiF eos;
	C1DProblem pr = prLiFSliceEuler; //prLiF; //prTestLagrange1D;
	C1DField* fldptr = new C1DField(pr);
	CHLLCRiemannSolver hllc;
	F1DENO2Reconstruction rec(*fldptr);
	C1D2ndOrderMethod mtd(hllc, rec);
	double _dtt[] = { pr.tmin,
					 .05e-9,
					 .1e-9, .2e-9, .3e-9, .4e-9, .5e-9, .6e-9, .7e-9, .8e-9, .9e-9,
					 /*1.e-9, 2.e-9, 3.e-9, 3.128e-9, 4.e-9, 5.e-9, 6.e-9, 7.e-9, 8.e-9, 9.e-9,*/
					 pr.tmax };
	vector<double> dtt = vector<double>(_dtt, _dtt + sizeof(_dtt) / sizeof(double));
	COutput outp = COutput(pr, outputDir, dtt);
	F1DSimulation sim = F1DSimulation(pr, eos, *fldptr, mtd, outp);
	sim.run();
	delete fldptr;

	// Uncomment for Lagrange 1D code in new architecture
	/*FEOSIdeal eos(1.4);
	C1DProblem pr = prToro1Idealtest;
	C1DFieldPrimitive *fldptr = new C1DFieldPrimitive(pr);
	C1DMethodSamarskii mtd;
	double _dtt[] = {pr.tmin, pr.tmax};
	vector<double> dtt = vector<double>(_dtt, _dtt + sizeof(_dtt) / sizeof(double));
	COutput outp = COutput(pr, outputDir, dtt);
	F1DSimulationLagrange sim = F1DSimulationLagrange(pr, eos, *fldptr, mtd, outp);
	sim.run();
	delete fldptr; /


	// Uncomment for LaserVT test problem
    /*FEOSMGAlPrecise6 eos;
	FEOSMieGruneisenAl eos;
	FEOSIdeal eos = FEOSIdeal(3.9);
	C1DProblem pr = prVTAlMGTest2_2;
	FEOSIdeal eos = FEOSIdeal(1.4);

	double _rho = 2413.,_p = 1.e9;
	double _e = eos.gete(_rho, _p);
	double _pro = eos.getdpdrho(_rho, _p);
	double _pe = eos.getdpde(_rho, _e);
	double c1 = eos.getc(_rho, _p);

	double c2 = sqrt(_p*_pe/_rho/_rho + _pro);

	C1DProblem pr = prToro1Idealtest;

	C1DField *fldptr = new C1DField(pr);
	CHLLRiemannSolver hll;

	F1DENO2Reconstruction eno2rec=F1DENO2Reconstruction(*fldptr);
	C1D2ndOrderMethod mtd = C1D2ndOrderMethod(roe, eno2rec);
	//C1DGodunovTypeMethod mtd = C1DGodunovTypeMethod(ex);
	//C1DBGKMethod mtd = C1DBGKMethod(bgk);
	C1DGodunovTypeMethod mtd = C1DGodunovTypeMethod(hll);
	// C1DLFMethod mtd = C1DLFMethod(lf);
	double _dtt[] = {pr.tmin, pr.tmax};
	vector<double> dtt = vector<double>(_dtt, _dtt+sizeof(_dtt)/sizeof(double));
	COutput outp = COutput(pr, outputDir, dtt);
	F1DSimulation sim = F1DSimulation(pr, eos, *fldptr, mtd, outp);
	sim.run();
	delete fldptr; */

	// Uncomment for SW-induced mechanism of holes formation
//	C1DProblem pr = prToro1Idealtest;  // prFedorAl;
//	FEOSIdeal eos = FEOSIdeal(1.4);
//	C1DField *fldptr = new C1DField(pr);
//	CHLLRiemannSolver hll;
//	CHLLCRiemannSolver hllc;
//	CLFRiemannSolver lf;
//	CGPSRiemannSolver gps;
//	CRoeRiemannSolver roe;
//	CRoeGeneralRiemannSolver roegen;
//	CBGKRiemannSolver bgk;
//	CExactRiemannSolver ex;
//	//F1DENO2Reconstruction eno2rec=F1DENO2Reconstruction(*fldptr);
//	//C1D2ndOrderMethod mtd = C1D2ndOrderMethod(hll, eno2rec);
//	//C1DGodunovTypeMethod mtd = C1DGodunovTypeMethod(roegen);
//	//C1DBGKMethod mtd = C1DBGKMethod(bgk);
//	C1DGodunovTypeMethod mtd = C1DGodunovTypeMethod(ex);
//	// C1DLFMethod mtd = C1DLFMethod(lf);
//	double _dtt[] = {pr.tmin, pr.tmax};
//	vector<double> dtt = vector<double>(_dtt, _dtt+sizeof(_dtt)/sizeof(double));
//	COutput outp = COutput(pr, outputDir, dtt);
//	F1DSimulation sim = F1DSimulation(pr, eos, *fldptr, mtd, outp);
//	sim.run();
//	delete fldptr;
	// Uncomment for ideal gas vs vacuum test
	/* //CEOSIdeal eos = CEOSIdeal(3.9);
	CEOSIdeal eos = CEOSIdeal(2.5);
	C1DProblem pr =  prLaserVTAlIdealTest2; //prIdealVacTest;
	C1DField *fldptr = new C1DField(pr);
	CExactRiemannSolver exrslv;
	C1DGodunovTypeMethodVacuum gdn = C1DGodunovTypeMethodVacuum(exrslv, pr.x0);
	// C1DGodunovTypeMethod gdn = C1DGodunovTypeMethod(exrslv);
	double _dtt[] = {pr.tmin, pr.tmax};
	vector<double> dtt = vector<double>(_dtt, _dtt+sizeof(_dtt)/sizeof(double));
	COutput outp = COutput(pr, outputDir, dtt);
	F1DSimulation sim = F1DSimulation(pr, eos, *fldptr, gdn, outp);
	sim.run();
	delete fldptr;	*/

	// Uncomment for Toro #1 test problem with ideal EOS

//	FEOSIdeal eos = FEOSIdeal(1.4);
//	// FEOSMGAlPrecise6 eos = FEOSMGAlPrecise6();
//	// C1DProblem pr = prLaserVTAlMGTest1;
//	// C1DProblem pr = prVTAlMGTest1;
//	C1DProblem pr = prToro2Idealtest;
//	// C1DProblem pr = prVTAlBel;
//	// C1DProblem pr = prHoles;
//	std::unique_ptr<C1DField> fldptr = std::make_unique<C1DField>(pr);
//	CHLLRiemannSolver hll;
//	CHLLCRiemannSolver hllc;
//	CLFGlobalRiemannSolver lfgl;
//	CLFRiemannSolver lf_old;
//	CExactRiemannSolver ex;
//	CRoeRiemannSolver roe;
//	CRoeGeneralRiemannSolver roegen;
////	F1DWENO5Reconstruction eno5rec = F1DWENO5Reconstruction(*fldptr);
////	F1DCharWiseWENO5Reconstruction eno5rec_charwise
////			= F1DCharWiseWENO5Reconstruction(*fldptr, eos);
////	F1DENO2Reconstruction eno2rec = F1DENO2Reconstruction(*fldptr);
////	F1DENO3Reconstruction eno3rec = F1DENO3Reconstruction(*fldptr);
////	C1D2ndOrderLFGlobalMethod mtd = C1D2ndOrderLFGlobalMethod(
////				lfgl, eno5rec_charwise);
////	C1D2ndOrderMethod mtd = C1D2ndOrderMethod(hllc,eno2rec);
//	// C1DLFGlobalMethod mtd = C1DLFGlobalMethod(lfgl);
//	// C1DGodunovTypeMethodVacuum mtd = C1DGodunovTypeMethodVacuum(hllc);
//	// CExactRiemannSolver exrslv;
//	// C1DGodunovTypeMethodVacuum mtd = C1DGodunovTypeMethodVacuum(hllc, pr.x0);
//	C1DGodunovTypeMethod mtd = C1DGodunovTypeMethod(hllc);
//	// C1DGodunovTypeMethod mtd = C1DGodunovTypeMethod(ex);
//	double _dtt[] = {pr.tmin, pr.tmax};
//	vector<double> dtt = vector<double>(_dtt, _dtt+sizeof(_dtt)/sizeof(double));
//	COutput outp = COutput(pr, outputDir, dtt);
////	ForwardEuler default_ode_solver = ForwardEuler(pr, eos, *fldptr, mtd);
////	SSPERK3_3 ssprk = SSPERK3_3(pr, eos, *fldptr, mtd);
////	eBDF5 e_bdf5 = eBDF5(pr, eos, *fldptr, mtd, ssprk);
////	ERK6_5 erk = ERK6_5(pr, eos, *fldptr, mtd);
////	F1DSimulation sim = F1DSimulation(pr, eos, *fldptr, mtd, outp, ssprk);
//	F1DSimulation sim = F1DSimulation(pr, eos, *fldptr, mtd, outp);
//	sim.run();


	// Uncomment for metal problems
	// CSolver s = CSolver();
	// s->goGlass("task-Ru-glass.txt");
	//s->goGlass("task-Ru-glass-optic-1000.txt");
	//s->goGlass("task-Ru-glass-optic-2000.txt");
	//s->goGlass("task-Ru-glass-optic-3000.txt");
	//s->goGlass("task-Fe-glass.txt");

	//s->goGlass("task-Au-water-simple.txt");
	//s->goGlass("task-Au-water.txt");
	//s->goGlass("task-Au-water-vacuum.txt");

	//s->goGlass("task - Ni.txt");
	//s->goGlass("task.txt"); // для железа и вообще
	//s->go("task-5Si.txt");
	//s->goAuSpall("task.txt"); // для золота
	// s.goEuler("task-toro-1.txt"); // для Эйлеровых задач
	//s->goEulerMovingMesh("task.txt"); // для Эйлеровых задач c движущейся сеткой
	//s->goEuler("task-LH1D.txt");          // Задача LH (стекло-золото-вакуум), одномерное приближение
    //s->goEuler("task-LH1D-p=123GPA.txt"); // Задача LH (стекло-золото-вакуум), одномерное приближение, с более мелкой ступенькой

	// s->goEuler("task-eosbin-test-2.txt"); // Тест на двучленное УРС от Паши -- он же тест на лазерное облучение объемной мишеги с идеальным УРС
	//CSolver *s = new CSolver;
	//s->goEuler("task-laservt-ideal-test.txt"); // Тест на лазерное облучение идеальной мишени с УРС от Н.А.
	//s->goEuler("task-eosbin-toro-test-5.txt"); // Тест на двучленное УРС (при нулевых параметрах должно работать как идеальное)

	//s->goEuler("task-LH1D-aux-1.txt"); // Задача LH (стекло-золото-вакуум), левый разрыв стекло-золото
	//s->goEuler("task-eosbin-test-1.txt"); // Первый тест Глайстера по двучленному уравнению состояния
	// Uncomment for Euler problems
	//s->goEuler("task-eosbin-toro-test-1.txt");
	//s->goEuler("task-test-NB.txt");
	// Uncomment for HLL flux Toro test
	/*EOSBin eos = EOSBin(1.4, 0., 0.);
	CTask prHLLTest = CTask(TaskType::RP1D, eos, 1., .75, 1., .125, 0., 1., 0., 1., 0., .2, .3, 100, .9, MethodType::hll);
    CSolver* s = new CSolver(prHLLTest);
	s->goEuler();
	//
	delete s;*/

	///////////////////////////////////////
	// Uncomment for convergence rate calculation
	/*ofstream ofs, ofs1;
	ofs.open("conv-rate-runge.dat", ios::out);
	ofs << "TITLE = Convergence rates" << endl;
	ofs << "VARIABLES = \"h\",\"err_L1\",\"err_L2\",\"conv_Integral\"" << endl;
	ofs << "ZONE T=\"GPS\" I=6 F=POINT " << endl;
	ofs1.open("conv-rate-runge-all.dat", ios::out);
	ofs1 << "TITLE = Runge set of convergence rates for different schemes" << endl;
	ofs1 << "VARIABLES = \"h\",\"P\"" << endl;
	ofs1 << "ZONE T=\"BGK\" I=4 F=POINT" << endl;
	double f[] = {0., 0., 0., 0., 0., 0.};
	double h[] = {.01, .01/2., .01/4., .01/8., .01/16., .01/32.};
	unsigned counter = 0;
	CSolver *s = new CSolver;
	s->goEuler("task-toro-1.txt"); // для Эйлеровых задач
	ofs << s->getdx() << " " << s->calcProblemL1NormRo() << " " << s->calcProblemL2NormRo() << " " << s->convIntegral << endl;
	f[counter++] = s->convIntegral;
	delete s;
	s = new CSolver;
	s->goEuler("task-toro-1-2.txt"); // для Эйлеровых задач
	ofs << s->getdx() << " " << s->calcProblemL1NormRo() << " " << s->calcProblemL2NormRo() << " " << s->convIntegral << endl;
	f[counter++] = s->convIntegral;
	delete s;
	s = new CSolver;
	s->goEuler("task-toro-1-4.txt"); // для Эйлеровых задач
	ofs << s->getdx() << " " << s->calcProblemL1NormRo() << " " << s->calcProblemL2NormRo() << " " << s->convIntegral << endl;
	f[counter++] = s->convIntegral;
	delete s;
	s = new CSolver;
	s->goEuler("task-toro-1-8.txt"); // для Эйлеровых задач
	ofs << s->getdx() << " " << s->calcProblemL1NormRo() << " " << s->calcProblemL2NormRo() << " " << s->convIntegral << endl;
	f[counter++] = s->convIntegral;
	delete s;
	s = new CSolver;
	s->goEuler("task-toro-1-16.txt"); // для Эйлеровых задач
	ofs << s->getdx() << " " << s->calcProblemL1NormRo() << " " << s->calcProblemL2NormRo() << " " << s->convIntegral << endl;
	f[counter++] = s->convIntegral;
	delete s;
	s = new CSolver;
	s->goEuler("task-toro-1-32.txt"); // для Эйлеровых задач
	ofs << s->getdx() << " " << s->calcProblemL1NormRo() << " " << s->calcProblemL2NormRo() << " " << s->convIntegral << endl;
	f[counter++] = s->convIntegral;
	delete s;

	ofs1 << h[1] << " " << log10((f[0]-f[1])/(f[1]-f[2]))/log10(2.) << endl;
	ofs1 << h[2] << " " << log10((f[1]-f[2])/(f[2]-f[3]))/log10(2.) << endl;
	ofs1 << h[3] << " " << log10((f[2]-f[3])/(f[3]-f[4]))/log10(2.) << endl;
	ofs1 << h[4] << " " << log10((f[3]-f[4])/(f[4]-f[5]))/log10(2.) << endl;

	ofs.close();
	ofs1.close();*/
	///////////////////////////////////////

//	s = new CSolver;
//	s->go("al_glass_1000_400.txt");
//	delete s;


	//  s.go("test.txt");
	//  s.go("input6_euler.txt");
	//	s.testHeatStage("input_lagrange_heat_test.txt", "heat.dat");
	//	s.testExchangeStage("input_lagrange_exchange_test.txt", "exchange.dat");
	//  s.testGodunovScheme("input_godunov_test.txt");

/*	s->go("al_glass_testEOS_v200.txt");
	delete s;P

	s = new CSolver;
	s->go("al_glass_testEOS_v500.txt");
	delete s;

	s = new CSolver;
	s->go("al_glass_testEOS_v1000.txt");
	delete s;

	s = new CSolver;
	s->go("al_glass_testEOS_v2000.txt");
	delete s;

*/

	delete[] INPUT_FOLDER;
	delete[] OUTPUT_FOLDER;
	return 0;
}

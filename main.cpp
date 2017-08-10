#define _CRT_SECURE_NO_WARNINGS

#include<string.h>

#include "defines.h"
#include "solver.h"
char* INPUT_FOLDER = new char[_MAX_PATH];  //"calc/";
char* OUTPUT_FOLDER = new char[_MAX_PATH]; //"calc/output/";

//#include <vld.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

int main(int argc, char *argv[])
{
	if(argc > 1)
	{
		strcpy(INPUT_FOLDER, argv[1]);
		strcpy(OUTPUT_FOLDER, INPUT_FOLDER);
		strcat(OUTPUT_FOLDER, "output/");
	}
	else
	{
		strcpy(INPUT_FOLDER, "calc/");
		strcpy(OUTPUT_FOLDER, "calc/output/");
	}
		
    //	s.go("input6_ideal_lagrange.txt");
	//	s.go("input6_table_lagrange.txt");
	//	s.go("input6_table_lagrange_he.txt");
	//  s.go("input6_table_lagrange_inner.txt");
	//	s.go("input6_source_shock.txt");
	//  s.go("input6_source_shock_ideal.txt");
	//  s.go("input6_test_shockwave.txt");

	//	s->go("al_vac_350.txt");

//	s->go("al_vac_500_cold.txt");
//	delete s;
	
	// Uncomment for metal problems
	CSolver *s = new CSolver;
	//s->goGlass("task-Ru-glass-1000.txt");
	//s->goGlass("task-Ru-glass-optic-1000.txt");
	//s->goGlass("task-Ru-glass-optic-2000.txt");
	//s->goGlass("task-Ru-glass-optic-3000.txt");
	//s->goGlass("task - Ni.txt");
	//s->goGlass("task.txt"); // дл€ железа и вообще
	//s->go("task-5Si.txt");
	//s->goAuSpall("task.txt"); // дл€ золота
	//s->goEuler("task-toro-1.txt"); // дл€ Ёйлеровых задач
	//s->goEulerMovingMesh("task.txt"); // дл€ Ёйлеровых задач c движущейс€ сеткой
	s->goEuler("task-LH1D.txt"); // «адача LH (стекло-золото-вакуум), одномерное приближение
	//s->goEuler("task-LH1D-aux-1.txt"); // «адача LH (стекло-золото-вакуум), левый разрыв стекло-золото
	delete s;

	// Uncomment for Euler problems
/*	CSolver *s = new CSolver;
	s->goEuler("task-toro-1.txt");
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
	s->goEuler("task-toro-1.txt"); // дл€ Ёйлеровых задач
	ofs << s->getdx() << " " << s->calcProblemL1NormRo() << " " << s->calcProblemL2NormRo() << " " << s->convIntegral << endl;
	f[counter++] = s->convIntegral;
	delete s;
	s = new CSolver;
	s->goEuler("task-toro-1-2.txt"); // дл€ Ёйлеровых задач
	ofs << s->getdx() << " " << s->calcProblemL1NormRo() << " " << s->calcProblemL2NormRo() << " " << s->convIntegral << endl;
	f[counter++] = s->convIntegral;
	delete s;
	s = new CSolver;
	s->goEuler("task-toro-1-4.txt"); // дл€ Ёйлеровых задач
	ofs << s->getdx() << " " << s->calcProblemL1NormRo() << " " << s->calcProblemL2NormRo() << " " << s->convIntegral << endl;
	f[counter++] = s->convIntegral;
	delete s;
	s = new CSolver;
	s->goEuler("task-toro-1-8.txt"); // дл€ Ёйлеровых задач
	ofs << s->getdx() << " " << s->calcProblemL1NormRo() << " " << s->calcProblemL2NormRo() << " " << s->convIntegral << endl;
	f[counter++] = s->convIntegral;
	delete s;
	s = new CSolver;
	s->goEuler("task-toro-1-16.txt"); // дл€ Ёйлеровых задач
	ofs << s->getdx() << " " << s->calcProblemL1NormRo() << " " << s->calcProblemL2NormRo() << " " << s->convIntegral << endl;
	f[counter++] = s->convIntegral;
	delete s;
	s = new CSolver;
	s->goEuler("task-toro-1-32.txt"); // дл€ Ёйлеровых задач
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

	delete INPUT_FOLDER;
	delete OUTPUT_FOLDER;
	return 0;
}




#define _CRT_SECURE_NO_WARNINGS

#include "defines.h"
#include "eosFigures.h"
#include "task.h"

#include <stdio.h>
#include <iostream>
#include <math.h>
#include <string.h>

#include <fstream>

using namespace std;

EOSFigures::EOSFigures(EOSOld* _eos) : eos(_eos)
{}

void EOSFigures::writeIsoterms(char *outputFileName)
{


/* — кодом, приведенным ниже, разобратьс€ -- чтобы было достаточное 
   количество зон, сетка изотерм, ограничени€ по давлению/энергии 
   (пр€моугольник допустимых давлений/энергий), видны все эффекты 
   на изотермах (реальный газ, фазовые переходы, отрицательные 
   давлени€). ј там посмотрим. ;)
*/

	int i=0;
	double pi=0, ti=0, te=0, ei=0, ee=0, e=0, ro=0;

	char output_file[_MAX_PATH];
	
	int iMax = 1000;

	strcpy(output_file, OUTPUT_FOLDER);
	strcat(output_file, outputFileName);

	double    ro_0 =  eos->getMIN_RO() + 1.0;
/*	double    ro_N =  eos->getMAX_RO() - 1.0;
	double ro_step = (eos->getMAX_RO() - eos->getMIN_RO()) / iMax;*/
	
	double ro_step = 5.0;

	FILE *f1  = fopen(output_file, "w+");

	fprintf(f1, "TITLE=\"Isoterms: Pi(Ro), Ei(Ro)\"\n");
	fprintf(f1, "VARIABLES=\"Ro\",\"Pi\",\"Ei\"\n");

	ti= 300.0;
	fprintf(f1, "ZONE T=\"%dK\" I=%d F=POINT\n", (int)ti, iMax);
	for(i=0; i<iMax; i++)
	{
		ro=ro_0+i*ro_step;
		pi=eos->getpi(ro, ti);
		ei=eos->getei(ro, ti);
		fprintf(f1, "%f %e %e\n", ro/1000.0, pi/1.0e+9, ei);
	}

	ti= 400.0;
	fprintf(f1, "ZONE T=\"%dK\" I=%d F=POINT\n", (int)ti, iMax);
	for(i=0; i<iMax; i++)
	{
		ro=ro_0+i*ro_step;
		pi=eos->getpi(ro, ti);
		ei=eos->getei(ro, ti);
		fprintf(f1, "%f %e %e\n", ro/1000.0, pi/1.0e+9, ei);
	}

	ti= 500.0;	
	fprintf(f1, "ZONE T=\"%dK\" I=%d F=POINT\n", (int)ti, iMax);
	for(i=0; i<iMax; i++)
	{
		ro=ro_0+i*ro_step;
		pi=eos->getpi(ro, ti);
		ei=eos->getei(ro, ti);
		fprintf(f1, "%f %e %e\n", ro/1000.0, pi/1.0e+9, ei);
	}

	ti= 600.0;	
	fprintf(f1, "ZONE T=\"%dK\" I=%d F=POINT\n", (int)ti, iMax);
	for(i=0; i<iMax; i++)
	{
		ro=ro_0+i*ro_step;
		pi=eos->getpi(ro, ti);
		ei=eos->getei(ro, ti);
		fprintf(f1, "%f %e %e\n", ro/1000.0, pi/1.0e+9, ei);
	}

	ti= 700.0;	
	fprintf(f1, "ZONE T=\"%dK\" I=%d F=POINT\n", (int)ti, iMax);
	for(i=0; i<iMax; i++)
	{
		ro=ro_0+i*ro_step;
		pi=eos->getpi(ro, ti);
		ei=eos->getei(ro, ti);
		fprintf(f1, "%f %e %e\n", ro/1000.0, pi/1.0e+9, ei);
	}

	ti= 800.0;	
	fprintf(f1, "ZONE T=\"%dK\" I=%d F=POINT\n", (int)ti, iMax);
	for(i=0; i<iMax; i++)
	{
		ro=ro_0+i*ro_step;
		pi=eos->getpi(ro, ti);
		ei=eos->getei(ro, ti);
		fprintf(f1, "%f %e %e\n", ro/1000.0, pi/1.0e+9, ei);
	}


	ti= 900.0;
	fprintf(f1, "ZONE T=\"%dK\" I=%d F=POINT\n", (int)ti, iMax);
	for(i=0; i<iMax; i++)
	{
		ro=ro_0+i*ro_step;
		pi=eos->getpi(ro, ti);
		ei=eos->getei(ro, ti);
		fprintf(f1, "%f %e %e\n", ro/1000.0, pi/1.0e+9, ei);
	}

	ti= 950.0;
	fprintf(f1, "ZONE T=\"%dK\" I=%d F=POINT\n", (int)ti, iMax);
	for(i=0; i<iMax; i++)
	{
		ro=ro_0+i*ro_step;
		pi=eos->getpi(ro, ti);
		ei=eos->getei(ro, ti);
		fprintf(f1, "%f %e %e\n", ro/1000.0, pi/1.0e+9, ei);
	}

	ti= 928.0;
	fprintf(f1, "ZONE T=\"%dK\" I=%d F=POINT\n", (int)ti, iMax);
	for(i=0; i<iMax; i++)
	{
		ro=ro_0+i*ro_step;
		pi=eos->getpi(ro, ti);
		ei=eos->getei(ro, ti);
		fprintf(f1, "%f %e %e\n", ro/1000.0, pi/1.0e+9, ei);
	}

	ti= 1000.0;
	fprintf(f1, "ZONE T=\"%dK\" I=%d F=POINT\n", (int)ti, iMax);
	for(i=0; i<iMax; i++)
	{
		ro=ro_0+i*ro_step;
		pi=eos->getpi(ro, ti);
		ei=eos->getei(ro, ti);
		fprintf(f1, "%f %e %e\n", ro/1000.0, pi/1.0e+9, ei);
	}

	ti= 1100.0;
	fprintf(f1, "ZONE T=\"%dK\" I=%d F=POINT\n", (int)ti, iMax);
	for(i=0; i<iMax; i++)
	{
		ro=ro_0+i*ro_step;
		pi=eos->getpi(ro, ti);
		ei=eos->getei(ro, ti);
		fprintf(f1, "%f %e %e\n", ro/1000.0, pi/1.0e+9, ei);
	}

	ti= 1200.0;
	fprintf(f1, "ZONE T=\"%dK\" I=%d F=POINT\n", (int)ti, iMax);
	for(i=0; i<iMax; i++)
	{
		ro=ro_0+i*ro_step;
		pi=eos->getpi(ro, ti);
		ei=eos->getei(ro, ti);
		fprintf(f1, "%f %e %e\n", ro/1000.0, pi/1.0e+9, ei);
	}

	ti= 1300.0;
	fprintf(f1, "ZONE T=\"%dK\" I=%d F=POINT\n", (int)ti, iMax);
	for(i=0; i<iMax; i++)
	{
		ro=ro_0+i*ro_step;
		pi=eos->getpi(ro, ti);
		ei=eos->getei(ro, ti);
		fprintf(f1, "%f %e %e\n", ro/1000.0, pi/1.0e+9, ei);
	}

	ti= 1400.0;
	fprintf(f1, "ZONE T=\"%dK\" I=%d F=POINT\n", (int)ti, iMax);
	for(i=0; i<iMax; i++)
	{
		ro=ro_0+i*ro_step;
		pi=eos->getpi(ro, ti);
		ei=eos->getei(ro, ti);
		fprintf(f1, "%f %e %e\n", ro/1000.0, pi/1.0e+9, ei);
	}

	ti=1500.0;
	fprintf(f1, "ZONE T=\"%dK\" I=%d F=POINT\n", (int)ti, iMax);
	for(i=0; i<iMax; i++)
	{
		ro=ro_0+i*ro_step;
		pi=eos->getpi(ro, ti);
		ei=eos->getei(ro, ti);
		fprintf(f1, "%f %e %e\n", ro/1000.0, pi/1.0e+9, ei);
	}

	ti= 1600.0;
	fprintf(f1, "ZONE T=\"%dK\" I=%d F=POINT\n", (int)ti, iMax);
	for(i=0; i<iMax; i++)
	{
		ro=ro_0+i*ro_step;
		pi=eos->getpi(ro, ti);
		ei=eos->getei(ro, ti);
		fprintf(f1, "%f %e %e\n", ro/1000.0, pi/1.0e+9, ei);
	}
	
	ti= 1700.0;
	fprintf(f1, "ZONE T=\"%dK\" I=%d F=POINT\n", (int)ti, iMax);
	for(i=0; i<iMax; i++)
	{
		ro=ro_0+i*ro_step;
		pi=eos->getpi(ro, ti);
		ei=eos->getei(ro, ti);
		fprintf(f1, "%f %e %e\n", ro/1000.0, pi/1.0e+9, ei);
	}

	ti= 1800.0;
	fprintf(f1, "ZONE T=\"%dK\" I=%d F=POINT\n", (int)ti, iMax);
	for(i=0; i<iMax; i++)
	{
		ro=ro_0+i*ro_step;
		pi=eos->getpi(ro, ti);
		ei=eos->getei(ro, ti);
		fprintf(f1, "%f %e %e\n", ro/1000.0, pi/1.0e+9, ei);
	}

	ti= 1900.0;
	fprintf(f1, "ZONE T=\"%dK\" I=%d F=POINT\n", (int)ti, iMax);
	for(i=0; i<iMax; i++)
	{
		ro=ro_0+i*ro_step;
		pi=eos->getpi(ro, ti);
		ei=eos->getei(ro, ti);
		fprintf(f1, "%f %e %e\n", ro/1000.0, pi/1.0e+9, ei);
	}

	ti= 2000.0;
	fprintf(f1, "ZONE T=\"%dK\" I=%d F=POINT\n", (int)ti, iMax);
	for(i=0; i<iMax; i++)
	{
		ro=ro_0+i*ro_step;
		pi=eos->getpi(ro, ti);
		ei=eos->getei(ro, ti);
		fprintf(f1, "%f %e %e\n", ro/1000.0, pi/1.0e+9, ei);
	}
	
	ti= 3000.0;
	fprintf(f1, "ZONE T=\"%dK\" I=%d F=POINT\n", (int)ti, iMax);
	for(i=0; i<iMax; i++)
	{
		ro=ro_0+i*ro_step;
		pi=eos->getpi(ro, ti);
		ei=eos->getei(ro, ti);
		fprintf(f1, "%f %e %e\n", ro/1000.0, pi/1.0e+9, ei);
	}

	ti= 4000.0;
	fprintf(f1, "ZONE T=\"%dK\" I=%d F=POINT\n", (int)ti, iMax);
	for(i=0; i<iMax; i++)
	{
		ro=ro_0+i*ro_step;
		pi=eos->getpi(ro, ti);
		ei=eos->getei(ro, ti);
		fprintf(f1, "%f %e %e\n", ro/1000.0, pi/1.0e+9, ei);
	}

	
	ti= 5000.0;
	fprintf(f1, "ZONE T=\"%dK\" I=%d F=POINT\n", (int)ti, iMax);
	for(i=0; i<iMax; i++)
	{
		ro=ro_0+i*ro_step;
		pi=eos->getpi(ro, ti);
		ei=eos->getei(ro, ti);
		fprintf(f1, "%f %e %e\n", ro/1000.0, pi/1.0e+9, ei);
	}

	ti= 6000.0;
	fprintf(f1, "ZONE T=\"%dK\" I=%d F=POINT\n", (int)ti, iMax);
	for(i=0; i<iMax; i++)
	{
		ro=ro_0+i*ro_step;
		pi=eos->getpi(ro, ti);
		ei=eos->getei(ro, ti);
		fprintf(f1, "%f %e %e\n", ro/1000.0, pi/1.0e+9, ei);
	}

	ti= 7000.0;
	fprintf(f1, "ZONE T=\"%dK\" I=%d F=POINT\n", (int)ti, iMax);
	for(i=0; i<iMax; i++)
	{
		ro=ro_0+i*ro_step;
		pi=eos->getpi(ro, ti);
		ei=eos->getei(ro, ti);
		fprintf(f1, "%f %e %e\n", ro/1000.0, pi/1.0e+9, ei);
	}

	ti= 8000.0;
	fprintf(f1, "ZONE T=\"%dK\" I=%d F=POINT\n", (int)ti, iMax);
	for(i=0; i<iMax; i++)
	{
		ro=ro_0+i*ro_step;
		pi=eos->getpi(ro, ti);
		ei=eos->getei(ro, ti);
		fprintf(f1, "%f %e %e\n", ro/1000.0, pi/1.0e+9, ei);
	}

	ti= 9000.0;
	fprintf(f1, "ZONE T=\"%dK\" I=%d F=POINT\n", (int)ti, iMax);
	for(i=0; i<iMax; i++)
	{
		ro=ro_0+i*ro_step;
		pi=eos->getpi(ro, ti);
		ei=eos->getei(ro, ti);
		fprintf(f1, "%f %e %e\n", ro/1000.0, pi/1.0e+9, ei);
	}

	ti=10000.0;
	fprintf(f1, "ZONE T=\"%dK\" I=%d F=POINT\n", (int)ti, iMax);
	for(i=0; i<iMax; i++)
	{
		ro=ro_0+i*ro_step;
		pi=eos->getpi(ro, ti);
		ei=eos->getei(ro, ti);
		fprintf(f1, "%f %e %e\n", ro/1000.0, pi/1.0e+9, ei);
	}

	fclose(f1);

}

void EOSFigures::writePiEi(char* outputFileName)
{
	FILE *f1  = fopen(outputFileName, "w+");
	fprintf(f1, "TITLE=\"Graphics: Pi(Ei)\"\n");
	fprintf(f1, "VARIABLES=\"Ti\",\"Ei\",\"Pi\"\n");

	double ti = 300.0;
	double pi = 0.0;
	double ro = 15000.0;
	double ei = 0.0;

	for(int i=0; i<1000; i++)
	{
		ti+=90.0;
		ei  = eos->getei(ro, ti);
		pi  = eos->getpi(ro, ti);
		fprintf(f1, "%f %e %e\n", ti, ei, pi);
	}
	fclose(f1);
}

void EOSFigures::testPiT(double roMin, double roMax, double T, double nIntervals) {
	string fName = string(OUTPUT_FOLDER) + "\\" + "test-pi-ro.dat";
	ofstream ofs; 
	ofs.open(fName);
	ofs << "TITLE=\"EOS test curves: pi(ro) at T = " << T << "\"" << endl;
	ofs << "VARIABLES=\"ro [kg/m3]\",\"pi [Pa]\"" << endl;
	double ro = 0.;
	double delta = (roMax-roMin)/((double)(nIntervals+1));
	for(int i=0; i<=nIntervals; i++)	{
		ro = roMin + (double)i*delta;		
		ofs << ro << " " << eos->getpi(ro, T) << endl;
	}
	ofs.close();
}

void EOSFigures::testEiT(double roMin, double roMax, double T, double nIntervals) {
	string fName = string(OUTPUT_FOLDER) + "\\" + "test-ei-ro.dat";
	ofstream ofs; 
	ofs.open(fName);
	ofs << "TITLE=\"EOS test curves: ei(ro) at T = " << T << "\"" << endl;
	ofs << "VARIABLES=\"ro [kg/m3]\",\"ei [J/kg]\"" << endl;
	double ro = 0.;
	double delta = (roMax-roMin)/((double)(nIntervals+1));
	for(int i=0; i<=nIntervals; i++)	{
		ro = roMin + (double)i*delta;		
		ofs << ro << " " << eos->getei(ro, T) << endl;
	}
	ofs.close();
}

void EOSFigures::testeEOSTe(double ro, double T, double TeMin, double TeMax, double nIntervals) {
	string fName = string(OUTPUT_FOLDER) + "\\" + "test-pe-ee-ce-alpha-kappa-Te.dat";
	ofstream ofs; 
	ofs.open(fName);
	ofs << "TITLE=\"EOS test curves: pe(Te), ee(Te), ce(Te), alpha(Te), kappa(Te)\"" << endl;
	ofs << "VARIABLES=\"Te [K]\",\"pe [Pa]\",\"ee [Pa]\",\"ce [J/m3/K]\",\"alpha [W/m3/K]\",\"kappa [W/m/K]\"" << endl;
	double _te = 0.;
	double delta = (TeMax-TeMin)/((double)(nIntervals+1));
	for(int i=0; i<=nIntervals; i++) {
		_te = TeMin + (double)i*delta;		
		ofs << _te << " " << eos->getpe(ro, T, _te) << " " << eos->getee(ro, T, _te) << 
			          " " << eos->getce(ro, _te) << " " << eos->getAlpha(ro, T, _te) << 
					  " " << eos->getkappa(ro, T, _te) << endl;
	}
	ofs.close();
}





void EOSFigures::writeColdCurve(char* outputFileName)
{
	FILE *f1  = fopen(outputFileName, "w+");
	fprintf(f1, "TITLE=\"Cold curve: Pi(V/V0)\"\n");
	fprintf(f1, "VARIABLES=\"Ro, kg/m3\",\"V/V0\",\"Pi, GPa\"\n");

	double ti   = 1.0;
	double pi   = 0.0;
	double ro   = 499.0;
	double ro0  = 8900.0;
	double v_v0 = 0.0;
 
	for(int i=0; i<10000; i++)
	{
		ro += 1.0;
		pi  = eos->getpi(ro, ti)*1.0e-9;
		v_v0= 1.0/(ro/ro0);
		fprintf(f1, "%f %e %e\n", ro, v_v0, pi);
	}
	fclose(f1);
}



void EOSFigures::writePeEe(char* outputFileName)
{
	  FILE *f1  = fopen(outputFileName, "w+");
	fprintf(f1, "TITLE=\"Graphics: Pe(Ee)\"\n");
	fprintf(f1, "VARIABLES=\"Ee\",\"Pe\"\n");
	

	double ti = 300;
	double ro = 2700.0;
	double pe = 0.0;
	double te = 300;
	double ee = 0.0;

	for(int i=0; i<120; i++)
	{
		te += 100.0*i;
		ee  = eos->getee(ro, ti, te);
		pe  = eos->getee(ro, ti, te);

		fprintf(f1, "%e %e\n", ee, pe);
	}
	fclose(f1);
}


void EOSFigures::writeIsochor(char* outputFileName, double ro)
{
	double ti=0.0, ei=0.0, pi=0.0;
	char output_file[_MAX_PATH];
	double maxT  = eos->getMAX_T()-1.;
	double minT  = eos->getMIN_T()+1.;
	double maxRo = eos->getMAX_RO()-1.;
	double minRo = eos->getMIN_RO()+1.;
	int numPoints = 1000; 
	double tStep = (maxT-minT)/numPoints;
	strcpy(output_file, OUTPUT_FOLDER);
	strcat(output_file, outputFileName);
	FILE *f1  = fopen(output_file, "w+");
    fprintf(f1, "TITLE=\"Isochor: Pi(T), Ei(T)\"\n");
	fprintf(f1, "VARIABLES=\"Ti [K]\",\"Pi [Pa]\",\"Ei [J/kg]\"\n");
	for(int i=0; i<numPoints; i++)
	{
		ti=minT+i*tStep;
		pi=eos->getpi(ro, ti);
		ei=eos->getei(ro, ti);
		fprintf(f1, "%f %e %e\n", ti, pi, ei);
	}
	fclose(f1);
}



void EOSFigures::writeIsoterm(char* outputFileName, double ti)
{
	double ei=0.0, pi=0.0;

	char output_file[_MAX_PATH];

	double maxT  = eos->getMAX_T()-1.;
	double minT  = eos->getMIN_T()+1.;
	double maxRo = eos->getMAX_RO()-1.;
	double minRo = eos->getMIN_RO()+1.;
	int numPoints = 1000; 

	strcpy(output_file, OUTPUT_FOLDER);
	strcat(output_file, outputFileName);

	FILE *f1  = fopen(output_file, "w+");
    fprintf(f1, "TITLE=\"Isoterm: Pi(Ro), Ei(Ro)\"\n");
	fprintf(f1, "VARIABLES=\"Ro [kg/m3]\",\"Pi [Pa]\",\"Ei [J/kg]\"\n");

	double ro = minRo; double roStep = (maxRo-minRo)/numPoints;

	for(int i=0; i<numPoints; i++)
	{
		ro=minRo+roStep*i;
		pi=eos->getpi(ro, ti);
		ei=eos->getei(ro, ti);
		fprintf(f1, "%f %e %e\n", ro, pi, ei);
	}

	fclose(f1);
}



void EOSFigures::writeEiDiagram(char *outputFileName, double ro_temp, double ei,
		                double dv_temp, double tau, double dm)
{
	double ti_temp=0.0, ei_temp=0.0, pi_temp=0.0, F=0.0;
	double ei_temp_min=0.0, ei_temp_max=0.0, ei_step=0.0;

	FILE *f1  = fopen(outputFileName, "w+");
    fprintf(f1, "TITLE=\"F(Ei) diagram: F(Ei)\"\n");
	fprintf(f1, "VARIABLES=\"Ei\",\"F\"\n");

	ei_temp_min = eos->getei(ro_temp, 294);
	ei_temp_max = eos->getei(ro_temp, 99999);
	ei_step = (ei_temp_max - ei_temp_min)/100.0;
	ei_temp = ei_temp_min - ei_step;

	for(int i=0; i<100; i++)
	{
		ei_temp+=ei_step;
		ti_temp=eos->getti(ro_temp, ei_temp);
		pi_temp=eos->getpi(ro_temp, ti_temp);
		F = ei_temp - ei - tau/dm*dv_temp*pi_temp;

		fprintf(f1, "%e %e \n", ei_temp, F);
	}

	ei_temp = (ei_temp_max - ei_temp_min)/2.0;
	ti_temp=eos->getti(ro_temp, ei_temp);
	pi_temp=eos->getpi(ro_temp, ti_temp);
	F = ei_temp - ei - tau/dm*dv_temp*pi_temp;
	cout << "In the centre F = " << F << endl;

	fclose(f1);
}

void EOSFigures::writekappa(char *outputFileName, double ro)
{
	FILE *f1  = fopen(outputFileName, "w+");
    fprintf(f1, "TITLE=\"kappa(ti)\"\n");
	fprintf(f1, "VARIABLES=\"ti [eV]\",\"kappa [10^2 J/(m*s*K)\", \"v^2 [10e12 m2/s2]\",\"cev [J/(K*m3)]\"\n");

	double ti0 = 3500.0, dti = 200.0, ti = 0.0;

	for(int i=0; i<100; i++)
	{
		ti = ti0 + i*dti;
		
		double kB  = 1.38e-16;
		double TF  = 115000.0;
		double EF  = kB * TF;
		double me  = 9.109e-28;
		double vF  = sqrt(2.0*EF/me);

	double ro0 = 2700.0;
	
	double v  = sqrt(9.0/25.0*vF*vF*vF*vF + 9.0*kB*kB*ti*ti/me/me)*0.0001;
	double ce = eos->getce(ro, ti);

		fprintf(f1, "%e %e %e %e \n", ti/11605.0, eos->getkappa(ro, ti, ti)*1.0e-2, v/1.0e12, ce);
	}	

	fclose(f1);
}

void EOSFigures::writeEntropy(char *outputFileName)
{
/*	FILE *f1  = fopen(outputFileName, "w+");
    fprintf(f1, "TITLE=\"Entropy\"\n");
	fprintf(f1, "VARIABLES=\"Ro [kg/m3]\",\"Ti [K]\", \"S [J/K]\"\n");

	double ti0 = 300.0, dti = 3.0, ti = 0.0;

	double ro = 670.0;

	for(int i=0; i<10000; i++)
	{
		ti = ti0 + i*dti;
		
		double S = getEntropy(ro, ti, 300.0);

		fprintf(f1, "%f %f %f\n", ro, ti, S);
	}	

	fclose(f1);*/
}

void EOSFigures::writeAdiabatic(char *outputFileName)
{
	double ro0	   = 2700.0;
	double ti0     = 300.0;
	double step_ro = 50.0;

	double ro1=0.0, ro2=0.0, ti1=0.0, ti2=0.0, pi1=0.0, pi2=0.0, ei1=0.0, ei2=0.0;
	
	FILE *f1  = fopen(outputFileName, "w+");
    fprintf(f1, "TITLE=\"Adiabatic curve pi(ro)\"\n");
	fprintf(f1, "VARIABLES=\"Ro [kg/m3]\",\"Pi [GPa]\"\n");

	ro1 = ro0;
	ti1 = ti0;
	ei1 = eos->getei(ro1, ti1);
	pi1 = eos->getpi(ro1, ti1);
	fprintf(f1, "%e %e\n", ro1, pi1);

	for(int i=0; i<10; i++)
	{
		ei1 = eos->getei(ro1, ti1);
		pi1 = eos->getpi(ro1, ti1);
	
		ro2 = ro1 + step_ro;
		ei2 = ei1 - pi1*(1.0/ro2 - 1.0/ro1);
		ti2 = eos->getti(ro2, ei2);
		pi2  = eos->getpi(ro2, ti2);

		fprintf(f1, "%e %e\n", ro2, pi2);

		ro1 = ro2;
		ti1 = ti2;
	}
}

void EOSFigures::writeIsentrope(char *outputFileName)
{
	FILE *f1  = fopen(outputFileName, "w+");
    fprintf(f1, "TITLE=\"Isentrope Pi(Ro)\"\n");
	fprintf(f1, "VARIABLES=\"Ro [kg/m3]\",\"Ti [K]\",\"Pi [Pa]\"\n");

	double ti = 300.0, pi = 0.0;
	double s = 5.57278;
	double ro0 = 570.0, dro = 10.0, ro = 0.0;

	for(int i=0; i<250; i++)
	{
		ro = ro0 + i*dro;
		
		ti = getTiFromEntropy(ro, s, 300.0, 50000.0);
		pi = eos->getpi(ro, ti);

		fprintf(f1, "%f %f %f\n", ro, ti, pi);
	}	

	fclose(f1);


}

double EOSFigures::getTiFromEntropy(double ro, double s, double low_border, double high_border)
{
/*	double ti;
	double eps = 0.0001;
	double ti_dei, low_dei, high_dei;

	for(;;)
	{
		ti = (low_border + high_border)/2.0;

		if(fabs(low_border-high_border) < eps)
			return ti;

		ti_dei   = eos->getEntropy(ro, ti, 300.0) - s;
		low_dei  = eos->getEntropy(ro, low_border, 300.0) - s;
		high_dei = eos->getEntropy(ro, high_border, 300.0) - s;

		if(fabs(ti_dei) < eps)
			return ti;

		if(low_dei * high_dei > 0)
		{
			printf("\nEOSFigures::getTiFromEntropy(): Error in initial data -- no equation root on the interval\n");
			printf("ro=%e, s=%e, F(a)=%e, F(b)=%e\n", ro, s, low_dei, high_dei);
			exit(1);
		}

		if(ti_dei * low_dei > 0)
			low_border  += (high_border-low_border)/2.0;
		else
			high_border -= (high_border-low_border)/2.0;
	}*/
	return 0.;
}

void EOSFigures::writePiRoTDiagram(char *outputFileName)
{
    FILE *f1  = fopen(outputFileName, "w+");
    fprintf(f1, "TITLE=\"Pi(ro, Ti)\"\n");
    fprintf(f1, "VARIABLES=\"Ro, kg//m3\",\"Ti, K\", \"P, GPa\"\n");
    fprintf(f1, "ZONE I=100 J=100 F=POINT\n");

        double Tmin=1.0, Tmax=10000.0, dT=(Tmax-Tmin)/99.0;
        double Romin=500.0, Romax=5000.0, dRo=(Romax-Romin)/99.0;
        
    double ro=0.0, ti=0.0, pi=0.0;

        for(int i=0; i<100; i++)
        {               
                for(int j=0; j<100; j++)
                {
                        ro=Romin+i*dRo;
                        ti=Tmin+j*dT;
                        pi=eos->getpi(ro, ti)*1.0e-9;

                        fprintf(f1, "%f %f %f \n", ro, ti, pi);
                }
        }
        fclose(f1);
}

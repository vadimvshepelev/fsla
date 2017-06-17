#define _CRT_SECURE_NO_WARNINGS

#include "defines.h"
#include "eosTableTa.h"

#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include<iostream>

using namespace std;

double EOSTableTa::getpi(double ro, double ti)
{
   /* if(ro<MIN_RO){
		cout << "EOSTableTa::getpi(): MIN_RO reached! ro = " << ro << " ti = " << ti << endl;  
		return pi_table.interpolate(MIN_RO, ti);
	}*/
	double pi = pi_table.interpolate(ro, ti);	
	return pi;
}

double EOSTableTa::getei(double ro, double ti) 
{ 
//	if(ro<1500.) return 0.;
/*	if(ro<MIN_RO){
		cout << "EOSTableTa::getei(): MIN_RO reached! ro = " << ro << " ti = " << ti << endl;  
		return pi_table.interpolate(MIN_RO, ti);
	}*/
	double ei = ei_table.interpolate(ro, ti);
	return ei;
}

double EOSTableTa::getEntropy(double ro, double ti, double te)
{
	double s  = entropy_table.interpolate(ro, ti);
	return s;
}

double EOSTableTa::getti(double ro, double ei) 
{
	double ti = solve_ti(ro, ei, t_scale.getMin()*(1.+.01), t_scale.getMax()-1.0);
	return ti;
}

double EOSTableTa::getti_interp(double ro, double ei)
{
	
	return 0.0;
}

double EOSTableTa::getC(double ro, double ti, double te, double v)
{
	//double C = getC_an(ro, ti, te);
	//double C = phase_table.interpolate(ro, ti);
	return 3000.;
}

double EOSTableTa::getphase(double ro, double ti)
{
/*	if(ro<MIN_RO){
		cout << "EOSTableTa::getphase(): MIN_RO reached! ro = " << ro << " ti = " << ti << endl;  
		return pi_table.interpolate(MIN_RO, ti);
	}*/
	double phase = phase_table.interpolate(ro, ti);
	return phase;
}

double EOSTableTa::getmix(double ro, double ti)
{
	double mix = mix_table.interpolate(ro, ti);
	return mix;
}

double EOSTableTa::getnuWR(double ro, double ti, double te, double b, double Z)
{
	/* 2correct

	double kB  = 1.38e-16;
	double TF  = 115000.0;
	double EF  = kB * TF;
	double me  = 9.109e-28;
	double vF  = sqrt(2.0*EF/me);
	double h_  = 1.05e-27;
	double ro0 = 2700.0;
	double T3  = 933.0;
	double eV  = 11605.0;
	double kZ  = 0.45;
	double  e  = 4.80320e-10;
	
	double v  = sqrt(sqrt(9.0/25.0*vF*vF*vF*vF + 9.0*kB*kB*ti*ti/me/me));
	double ce = getce(ro, te, Z)* 10.0; // CGS!!!

	// Frequences

				double s  = 0.2474e15*te*te/eV/eV;
				double i  = 1.15e15*pow(te/eV, 0.28);
		double nu_ee_c = (1.0 / (1.0/s + 1.0/i)) * pow(2700.0/ro, 1.0/3.0);
	
				double nu_ei_sol = 3.2e11*ti;
				double nu_ei_liq = 0.72e14*ti/(130.0 + 0.0367*ti - 66700/ti);
				
				double     phase = 0.0;
				double       mix = 0.0;

				if (ti < getMAX_T())
				{
					phase = getphase(ro, ti);	
				      mix = getmix(ro, ti);
				}
				else
				{
					phase = 3.0;
					mix   = 0.0;
				}

				double  sol_frac = 0.0;
				double abs_phase = fabs(phase);

				if(abs_phase <= 1.0)
					sol_frac = 1.0;
				else if (abs_phase < 3.0)
					sol_frac = mix;
				else
					sol_frac = 0.0;
			
		double nu_ei_c = (sol_frac*nu_ei_sol + (1.0-sol_frac)*nu_ei_liq) * pow(2700.0/ro, 1.0/3.0);

				double kappa_plasma = 16.0 * sqrt(2.0)/pow(3.14159, 1.5) * kZ *kB *
					   pow(kB * te, 2.5)/Z/sqrt(me)/e/e/e/e;
		double nu_ee_pl = (1.0 - kZ) * v * v * ce / 3.0 / kappa_plasma;
		double nu_ei_pl =        kZ  * v * v * ce / 3.0 / kappa_plasma;
	
	double nu_ei  = 1.0/sqrt(1.0/nu_ei_c/nu_ei_c + 1.0/nu_ei_pl/nu_ei_pl); 
	double nu_ee  = 1.0/sqrt(1.0/nu_ee_c/nu_ee_c + 1.0/nu_ee_pl/nu_ee_pl); 

	double nu     = b*nu_ee + nu_ei;

	return nu;*/ return 0.0;
}

EOSTableTa::EOSTableTa(char* dirName, int EOSFlag, double _ro0)
{
	char buf[256];

	char filename[_MAX_PATH];
	strcpy(filename, TABLE_FOLDER);
	strcat(filename, dirName);
	strcat(filename, "/");
	strcat(filename, "data.man");

	FILE* f=fopen(filename, "r");
	
	printf("Processing %s:\n", filename);

	fgets(buf, sizeof(buf), f);
	fgets(buf, sizeof(buf), f);
	
	if(EOSFlag == 1)
	{
		fgets(buf, sizeof(buf), f);
		fgets(buf, sizeof(buf), f);
		fgets(buf, sizeof(buf), f);
		fgets(buf, sizeof(buf), f);
		fgets(buf, sizeof(buf), f);
		fgets(buf, sizeof(buf), f);
		fgets(buf, sizeof(buf), f);
		fgets(buf, sizeof(buf), f);
		fgets(buf, sizeof(buf), f);
	}

	fgets(buf, sizeof(buf), f);
	double tableTmin = atof(buf)*1.0e3;
	printf("tableTmin=%e\n", tableTmin);
	
	fgets(buf, sizeof(buf), f);
	double tableTmax = atof(buf)*1.0e3;
	printf("tableTmax=%e\n", tableTmax);
	
	fgets(buf, sizeof(buf), f);
	double tableVmin = atof(buf);
	printf("tableVmin=%e\n", tableVmin);

	fgets(buf, sizeof(buf), f);
	double tableVmax = atof(buf);
	printf("tableVmax=%e\n", tableVmax);
	
	fgets(buf, sizeof(buf), f);
	int nTScale = atoi(buf);
	printf("nTScale=%d\n", nTScale);

	fgets(buf, sizeof(buf), f);
	int nVScale = atoi(buf);
	printf("nVScale=%d\n", nVScale);

	fclose(f);

	if (EOSFlag == 0)
	{
		v_scale.create(nVScale, tableVmin, tableVmax, 1.0e-3);
		t_scale.create(nTScale, tableTmin, tableTmax, 1.0);
		
		entropy_table.create("entrop.tab", dirName, 1.0, &v_scale, &t_scale, EOSFlag);				
	}
	else if (EOSFlag == 1)
	{
		v_scale.getFromFile("volume.in", dirName, nVScale, tableVmin, tableVmax, 1.0e-3);
		t_scale.getFromFile("temper.in", dirName, nTScale, tableTmin, tableTmax, 1.0e3);		
	}
	
	   ei_table.create("energy.tab", dirName, 1.0e6, &v_scale, &t_scale, EOSFlag);
	   pi_table.create("pressu.tab", dirName, 1.0e9, &v_scale, &t_scale, EOSFlag);
	  //  C_table.create("sndvls.tab", dirName, 1.0e3, &v_scale, &t_scale, EOSFlag);

	if(EOSFlag == 1)
	{
		phase_table.create("physid.tab", dirName, 1.0, &v_scale, &t_scale, EOSFlag);
	    //  mix_table.create("PhsMix.tab", dirName, 1.0, &v_scale, &t_scale, EOSFlag); 
	}		

	MAX_T  = tableTmax; 
	MIN_T  = tableTmin;
	MAX_RO = 1.0/v_scale.getMin();
	MIN_RO = 1.0/v_scale.getMax();
	ro0 = _ro0;
	tableErrorFlag = 0;
	tableErrorCellNum = -1;

}

double EOSTableTa::getkappa(double ro, double ti, double te, double Z=3.0)
{
	double ro0		  = 16654.;
	double teKK		  = te/1000.;
	double tiKK		  = ti/1000.;
	double f_lin	  = 0.016 * tiKK;
	double f_sat	  = 0.063;
	double a          = pow( pow(f_lin, -2.5) + pow(f_sat, -2.5), -0.4 );
	double kappa_low  = teKK/a;
	double kappa_high = 520. + 0.12*(teKK + 7.)*(teKK + 7.);
	double kappa      = (ro*ro*ro/ro0/ro0/ro0)/(1./kappa_low + 1./kappa_high);
	return kappa; //[J/s/m/K]
}


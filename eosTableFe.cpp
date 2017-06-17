#define _CRT_SECURE_NO_WARNINGS
#include "defines.h"
#include "eosTableFe.h"

#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include<iostream>

using namespace std;


double EOSTableFe::getpi(double ro, double ti)
{
	double pi = pi_table.interpolate(ro, ti);

	//////////////////////////
//	if (pi < 0.0)
//		return 0.0;
	//////////////////////////
	
	return pi;
}


double EOSTableFe::getei(double ro, double ti) 
{ 
	double ei = ei_table.interpolate(ro, ti);
	return ei;
}

double EOSTableFe::getEntropy(double ro, double ti, double te)
{
	double s  = entropy_table.interpolate(ro, ti);
	return s;
}

double EOSTableFe::getti(double ro, double ei) 
{
	double ti = solve_ti(ro, ei, t_scale.getMin()+1.0, t_scale.getMax()-1.0);
	return ti;
}

double EOSTableFe::getti_interp(double ro, double ei)
{
	
	return 0.0;
}

double EOSTableFe::getC(double ro, double ti, double te) {	
	return 6000.;
}

double EOSTableFe::getphase(double ro, double ti)
{
	double phase = phase_table.interpolate(ro, ti);
	return phase;
}

double EOSTableFe::getmix(double ro, double ti)
{
	double mix = mix_table.interpolate(ro, ti);
	return mix;
}

double EOSTableFe::getnuWR(double ro, double ti, double te, double b)
{
	// Constants	
	/*
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

EOSTableFe::EOSTableFe(char* dirName, int EOSFlag, double _ro0)
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
	printf("tableTmin=%7.2eK\n", tableTmin);
	
	fgets(buf, sizeof(buf), f);
	double tableTmax = atof(buf)*1.0e3;
	printf("tableTmax=%7.2fK\n", tableTmax);
	
	fgets(buf, sizeof(buf), f);
	double tableVmin = atof(buf);            // 
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

	if(EOSFlag == 1)
	{
		phase_table.create("physid.tab", dirName, 1.0, &v_scale, &t_scale, EOSFlag);
	      mix_table.create("PhsMix.tab", dirName, 1.0, &v_scale, &t_scale, EOSFlag); 
	}		

	MAX_T  = tableTmax; 
	MIN_T  = tableTmin;
	MAX_RO = 1.0/tableVmin * 1.0e3;
	MIN_RO = 1.0/tableVmax * 1.0e3;

	ro0 = _ro0;

	tableErrorFlag = 0;
	tableErrorCellNum = -1;
}

double EOSTableFe::getkappa(double ro, double ti, double te)
{
/*	double        v2 = (2.52 + 0.405e-4*te)*1.e12; // [m2/s2]
	
	double nu_sol_ei = 4.9e15*pow( ((ti+520.0 )/1250.0), 1.7 );
	double nu_liq_ei = 1.0e15*pow( ((ti+1728.0)/1728.0), 0.7 );

	double A_ee   = 1.0;
	double kB  = 1.38e-23;
	double TF  = 115000.0;
	double EF  = kB * TF;
	double h_  = 1.05e-34;

	double nu_abr = (kB*kB*te*te)/h_/EF;
	double nu_pl  = EF/h_*pow(kB*te/EF, -1.5);
	double nu_ee  = A_ee * nu_abr*nu_pl/sqrt(nu_abr*nu_abr + nu_pl*nu_pl);
	
	double kappa_sol_int = 1.0/3.0 * getce(ro, te, Z) * v2 / (nu_sol_ei + nu_ee);
	double kappa_liq_int = 1.0/3.0 * getce(ro, te, Z) * v2 / (nu_liq_ei + nu_ee);
	double kappa_pl		 = 1.0e-5  * 0.31*pow(te, 2.5) /
						   log(1.0 + 1.74e-6*te*(1.34+1.38e-5*te)) * 1.0e-6;	

	double kappa_sol = sqrt(kappa_sol_int*kappa_sol_int +
							kappa_pl     *kappa_pl);
	double kappa_liq = sqrt(kappa_liq_int*kappa_liq_int +
							kappa_pl     *kappa_pl);


	double mix = getmix(ro, ti);
	double kappa     = mix*kappa_sol + (1.0-mix) * kappa_liq; */

	double ro0 = 7874.0;
	double K = 180.0 * te * te / 2.0;
	double S = 1.275e5 * pow(te, 1.3);
	double K_prime = 180.0 * te;
	double S_prime = 1.6575e5 * pow(te, 0.3);
	double cs = pow(ro/ro0, 0.6)*(K_prime*S*S + S_prime*K*K)/(K+S)/(K+S);

	double v2 = pow( (1.8 + 3.86e-6*pow(te, 1.25)), 0.8 )*1.e12; // [m2/s2]


			double	 r_l = -0.042 +  2.08e-4*ti + 5.126e-7*ti*ti;
			double	 r_h =  0.094 + 3.764e-4*ti - 5.644e-8*ti*ti; 
			double r_sol = 1.0e-6*pow(1.0/r_l/r_l/r_l/r_l/r_l/
										  r_l/r_l/r_l/r_l/r_l +
									  1.0/r_h/r_h/r_h/r_h/r_h/
										  r_h/r_h/r_h/r_h/r_h, -0.1);
	double nu_sol_si = 5.06e21 * r_sol;


					 r_l = 0.67 + 1.0e-4*ti;
					 r_h = 2.0;
			double r_liq = 1.0e-6*pow(1.0/r_l/r_l/r_l/r_l/r_l +
									  1.0/r_h/r_h/r_h/r_h/r_h, -0.2);							  
	double nu_liq_si = 5.06e21 * r_liq;

   		double   mix = getmix(ro, ti);
	double     nu_si = mix*nu_sol_si + (1.0-mix)*nu_liq_si;
    double     nu_se = 1.0/( 1.0/1.9e11/te + 1.0/14.7e14/pow(te, 0.1));
	
	double	    nu_s = nu_se + nu_si;

	double  kappa_s  = 1.0/3.0 * cs * v2 / nu_s;
	double  kappa_pl = 1.3e-10 * pow(te, 2.5);
	double  kappa    = pow(kappa_s*kappa_s + kappa_pl*kappa_pl, 0.5);
	
	return kappa;
}



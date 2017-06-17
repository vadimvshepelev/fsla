#define _CRT_SECURE_NO_WARNINGS

#include "defines.h"
#include "eosTableNi.h"

#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include<iostream>

using namespace std;


double EOSTableNi::getpi(double ro, double ti)
{
	double pi = pi_table.interpolate(ro, ti);

	//////////////////////////
//	if (pi < 0.0)
//		return 0.0;
	//////////////////////////
	
	return pi;
}


double EOSTableNi::getei(double ro, double ti) 
{ 
	double ei = ei_table.interpolate(ro, ti);
	return ei;
}

double EOSTableNi::getEntropy(double ro, double ti, double te)
{
	double s  = entropy_table.interpolate(ro, ti);
	return s;
}

double EOSTableNi::getti(double ro, double ei) 
{
	double ti = solve_ti(ro, ei, t_scale.getMin()+1.0, t_scale.getMax()-1.0);
	return ti;
}

double EOSTableNi::getti_interp(double ro, double ei)
{
	
	return 0.0;
}

double EOSTableNi::getC(double ro, double ti, double te) {
	return 3000.;
}

double EOSTableNi::getphase(double ro, double ti)
{
	double phase = phase_table.interpolate(ro, ti);
	return phase;
}

double EOSTableNi::getmix(double ro, double ti)
{
	double mix = mix_table.interpolate(ro, ti);
	return mix;
}

EOSTableNi::EOSTableNi(char* dirName, int EOSFlag, double _ro0)
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

double EOSTableNi::getkappa(double ro, double ti, double te)
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

	double ro0 = 8900.0;
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
	double  kappa    = pow(kappa_s*kappa_s + kappa_pl*kappa_pl, -0.5);
	
	return kappa;
}


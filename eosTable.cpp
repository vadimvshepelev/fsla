#include "defines.h"
#include "eosTable.h"

#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <iostream>

using namespace std;


double EOSTable::getpi(double ro, double ti)
{

/*	This is for 'normal', ordinary EOS.
	
	double pi = pi_table.interpolate(ro, ti);
*/


/*	This is for tables, divided by ro*T for precise derivatives calculations */

	double pi = pi_table.interpolate(ro, ti) * ro * ti;

/****************************************************************************/

	//////////// DEBUG!!!!!!!!!/////////////

//	if (pi<0) pi = 0.001;

	/////////////////////////////////////////

	return pi;
}


double EOSTable::getei(double ro, double ti) 
{ 

/*	This is for 'normal', ordinary EOS.*/
	
	double ei = ei_table.interpolate(ro, ti);
//////////////////////////////////

/*	This is for tables, divided by ro*T for precise derivatives calculations */

//	double ei = ei_table.interpolate(ro, ti) * ti;

/****************************************************************************/
	
	return ei;
}

double EOSTable::getEntropy(double ro, double ti)
{
	double s  = entropy_table.interpolate(ro, ti);
	return s;
}

double EOSTable::getti(double ro, double ei) {
/*	double ti = solve_ti(ro, ei, t_scale.getMin(), t_scale.getMax()-1.0);
    return ti;*/

	// Более точный способ. И, надеюсь, более быстрый.
	// Проверяем границы, есть ли решение у уравнения ei(ti)|ro=ro0 = ei
	double ei_min = getei(ro, getMIN_T());
	double ei_max = getei(ro, getMAX_T());
	if( (ei<ei_min)||(ei>ei_max) ) {
		printf("\nEOSTable::getti(): Error in initial data -- no equation root on the interval\n");
		printf("ro=%e, ei=%e, F(a)=%e, F(b)=%e\n", ro, ei, ei_min, ei_max);
		//return -1.;
		exit(1);
	} else if (ei == ei_min) {
		return getMIN_T();
	} else if (ei == ei_max) {
		return getMAX_T();
	}
	double ti = ei_table.calct(ro, ei);
	return ti;
}

double EOSTable::getti_interp(double ro, double ei)
{
	
	return 0.0;
}

double EOSTable::getC(double ro, double ti, double te)
{
	///////////DEBUG///////////////	
	//double C = C_table.interpolate(ro, ti);
	///////////////////////////////

	////// Это было раньше, пока обходились без Эйлера. Быть внимательным,
	////// поскольку эта скорость звука входит в рабочую Лагранжеву 
	////// конфигурацию для алюминия
	/////     double C = getC_an(ro, ti, te);
	////////////////////////////////////////////////////////////////////

//	double C = sqrt(getdpdro_rov_roE(ro, ti, te, v)+ v*getdpdrov_ro_roE(ro, ti, te, v)+(getei(ro, ti)+v*v/2.0+getpi(ro, ti)/ro)*getdpdroE_ro_rov(ro, ti, te, v));

	//////// DEBUG /////////
/*	double C2 = sqrt(5.0/3.0 * getpi(ro, ti) / ro);	

	if(fabs(C-C2) > 0.1*C)
	{
		printf("C achtung!!!\n");
	}
*/
	////////////////////////

	double C = 3000.;

	return C; // при замене этой строчки на return C2 аналитическое решение начинает сходиться с табличным

	//return sqrt(5.0/3.0*8.31*ti/27.0e-3);
}

double EOSTable::getphase(double ro, double ti)
{
	double phase = phase_table.interpolate(ro, ti);



	return phase;
}

double EOSTable::getmix(double ro, double ti)
{
	double mix = mix_table.interpolate(ro, ti);
	return mix;
}

double EOSTable::getnuWR(double ro, double ti, double te, double b, double Z)
{
	// Drude theory
/*Old Drude kappa
	// Constants	
	
	double kB  = 1.38e-16;
	double TF  = 115000.0;
	double EF  = kB * TF;
	double me  = 9.109e-28;
	double vF  = sqrt(2.0*EF/me);
	double h_  = 1.05e-27;
	double  b  = 1.0;
	double ro0 = 2700.0;
	double T3  = 933.0;
	
	double v  = sqrt(vF*vF + 3.0*kB*te/me);
	double ce = getce(ro, te);

	// Frequences

	double nu_ee = b * (EF/h_) * (te*te/TF/TF);

	double nu_ei_sol = 4.2e14*ti/T3;
	double nu_ei_liq = 1.1e14*ti/(130.0 + 0.0367*ti - 66700/ti);
	
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

	double nu_ei = sol_frac*nu_ei_sol + (1.0-sol_frac)*nu_ei_liq;

	double nu_deg = (nu_ei + nu_ee) * pow(ro/ro0, -1.0/3.0);
	double nu_pl  = (EF/h_) * pow(te/TF, -3.0/2.0) * pow(ro0/ro, 2.0/3.0);
	double nu     = 1.0/sqrt(1.0/nu_deg/nu_deg + 1.0/nu_pl/nu_pl); 
	
	double kappa  = 1.0/3.0*v*v*ce/nu;
 
	return kappa * 1.0e-5;*/
	
  

	// NEW DRUDE KAPPA !!!

	// Constants	
	
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
	double ce = getce(ro, te)* 10.0; // CGS!!!

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

	return nu;
}

complex<double> EOSTable::geteps(double ro, double ti, double te, double Z)
{	
	double NA			= 6.0e23;
	double Matom		= 27.0/NA;
	double nat			= ro/1000.0/Matom;
	double ne			= nat*Z;

	double lambda		= 400.0e-9;

	double omega_prob   = 2.0*3.14159/lambda*3.0e8;
	double omega_pl_2   = 4.0*3.14159*ne*4.8032e-10*4.8032e-10/1.7/9.109382e-28;
	double b			= 1.0; 
	double nuWR			= getnuWR(ro, ti, te, b);

	complex<double> eps0     = complex<double>(
							   
							   1.0 - omega_pl_2 / (omega_prob * omega_prob + nuWR * nuWR), 
							    
							   omega_pl_2 * nuWR / 
							          omega_prob / 
									( omega_prob * omega_prob + nuWR * nuWR ));

	//// DEBUG ////

	complex<double> epsPalik = complex<double>(-23.3795, 4.7628);

	complex<double> delta300 = epsPalik-eps0;

	////////////////


	// lambda 620 nm, b = 0.05
	// complex<double> delta_Palik = complex<double>(-18.70668, 19.04585); 

	// lambda 620 nm, b = 1
    // complex<double> delta_Palik = complex<double>(-18.70680, 19.04397); 
 
	// lambda 400 nm, b = 0.05
	// complex<double> delta_Palik = complex<double>(-9.21096, 4.45376);

	// lambda 400 nm, b = 1
	 complex<double> delta_Palik = complex<double>(-9.21098, 4.45326);


	complex<double> delta_bb = delta_Palik*
							   (0.7e15+getnuWR(ro, 300.0, 300.0, b))/
							   (0.7e15+getnuWR(ro, ti, te, b));
	
	complex<double> eps = eps0+delta_bb;

	return eps;
}

EOSTable::EOSTable()
{
	char buf[256];

	char filename[_MAX_PATH];
	strcpy(filename, TABLE_FOLDER);
	strcat(filename, "data.man");

	FILE* f=fopen(filename, "r");

	fgets(buf, sizeof(buf), f);
	fgets(buf, sizeof(buf), f);

	printf("Processing %s...\n", filename);

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

	v_scale.create(nVScale, tableVmin, tableVmax, 1.0e-3);
	t_scale.create(nTScale, tableTmin, tableTmax, 1.0);

	ei_table.create("energy.tab", "", 1.0e6, &v_scale, &t_scale, 0);
	pi_table.create("pressu.tab", "",1.0e9, &v_scale, &t_scale, 0);

////////////DEBUG!!!!!!!!!!!!!////////////////////////////////
	
//	C_table.create("sndvls.tab", "", 1.0e3, &v_scale, &t_scale);

//////////////////////////////////////////////////////////////

	////////DEBUG!!!!!!!!!!////////////

//	monotonizeTables();

	//////////////////////////////////

	MAX_T  = tableTmax; 
	MIN_T  = tableTmin;
	MAX_RO = 1.0/tableVmin * 1.0e3;
	MIN_RO = 1.0/tableVmax * 1.0e3;
}

EOSTable::EOSTable(string dirName, int EOSFlag, double _ro0) {
	char buf[256];
	char filename[_MAX_PATH];
	strcpy(filename, TABLE_FOLDER);
	strcat(filename, dirName.c_str());
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
	
	//////////////////////DEBUG!!!!!!!!!!!!!! Вставил идеальные, убрать за собой /////////

	   ei_table.create("energy.tab", dirName, 1.0e6, &v_scale, &t_scale, EOSFlag);
	   pi_table.create("pressu_roT.tab", dirName, 1.0e9, &v_scale, &t_scale, EOSFlag);
	
//	   dpdti_table.create("dpdti_ideal.tab", dirName, 1.0, &v_scale, &t_scale);
//	   dpdro_table.create("dpdro_ideal.tab", dirName, 1.0, &v_scale, &t_scale);
//	   dedti_table.create("dedti_ideal.tab", dirName, 1.0, &v_scale, &t_scale);

	
	//////////////////////////////////////////////////////////////////////////////////////

	if(EOSFlag == 1)
	{
		phase_table.create("physid.tab", dirName, 1.0, &v_scale, &t_scale, EOSFlag);
	      mix_table.create("PhsMix.tab", dirName, 1.0, &v_scale, &t_scale, EOSFlag); 
	}		
	
	if(EOSFlag == 1)
	{
		// Это вставка только для EOSFlag == 1. Исправляет следствия ошибки -- немонотонности файла с объемами.

	 	 ei_table.correctTable(nVScale, nTScale);
		 pi_table.correctTable(nVScale, nTScale);
	  phase_table.correctTable(nVScale, nTScale);
	    mix_table.correctTable(nVScale, nTScale);
	}


// 	C_table.create("sndvls.tab", dirName, 1.0e3, &v_scale, &t_scale, EOSFlag);
//	monotonizeTables();

	MAX_T  = tableTmax; 
	MIN_T  = tableTmin;
	MAX_RO = 1.0/tableVmin * 1.0e3;
	MIN_RO = 1.0/tableVmax * 1.0e3;
	ro0 = _ro0;
	tableErrorFlag = 0;
	tableErrorCellNum = -1;

}


void EOSTable::monotonizeTables(void)
{
	bool c1, c2;
	int nMax=0, nMin=0;
	
	for(int j=0; j<t_scale.getSize()-1; j++)
		for(int i=0; i<v_scale.getSize()-1; i++)
		{
			c1 = pi_table.isMonotone(i, j);
			c2 = ei_table.isMonotone(i, j);
			if( (!c1) && (!c2) )
			{
				pi_table.monotonize(i, j);
				ei_table.monotonize(i, j);
				 C_table.monotonize(i, j);
			}
		}
}


double EOSTable::getkappa(double ro, double ti, double te) {
	// Constants		
	double kB  = 1.38e-16;
	double TF  = 115000.0;
	double EF  = kB * TF;
	double me  = 9.109e-28;
	double vF  = sqrt(2.0*EF/me);
	double h_  = 1.05e-27;
	double  b  = 1.0;
	double ro0 = 2700.0;
	double T3  = 933.0;
	double eV  = 11605.0;
	double kZ  = 0.45;
	double  e  = 4.80320e-10;
	
	double v  = sqrt(sqrt(9.0/25.0*vF*vF*vF*vF + 9.0*kB*kB*ti*ti/me/me));
	double ce = getce(ro, te)* 10.0; // CGS!!!

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
					   pow(kB * te, 2.5)/3./sqrt(me)/e/e/e/e;
		double nu_ee_pl = (1.0 - kZ) * v * v * ce / 3.0 / kappa_plasma;
		double nu_ei_pl =        kZ  * v * v * ce / 3.0 / kappa_plasma;
	
	double nu_ei  = 1.0/sqrt(1.0/nu_ei_c/nu_ei_c + 1.0/nu_ei_pl/nu_ei_pl); 
	double nu_ee  = 1.0/sqrt(1.0/nu_ee_c/nu_ee_c + 1.0/nu_ee_pl/nu_ee_pl); 

	double nu     = nu_ee + nu_ei;

	double kappa  = 1.0/3.0*v*v*ce/nu;
 
	// return kappa * 1.0e-5;



	// Алюминий используется как стекло -- т.е. теплопроводность и обмен в нем отсутствуют
	
	return 0.;

}

////////////////// DEBUG!!!!!!!!!!!!!! ///////////////////////////////////

/////////////// Тут все от идеального газа, убрать за собой!!! //////////

double EOSTable::getdpdro(double ro, double ti, double te)
{
/*	double deriv1 = (getpi(1.001*ro, ti)-getpi(ro, ti))/(0.001*ro);
	double deriv2 = (getpi(1.002*ro, ti)-getpi(ro, ti))/(0.002*ro);
	double deriv3 = (getpi(1.005*ro, ti)-getpi(ro, ti))/(0.005*ro);
	double deriv4 = (getpi(1.01*ro, ti)-getpi(ro, ti))/(0.01*ro);
	double deriv5 = (getpi(1.02*ro, ti)-getpi(ro, ti))/(0.02*ro);
	double deriv6 = (getpi(1.05*ro, ti)-getpi(ro, ti))/(0.05*ro);	
	double deriv7 = (getpi(1.1*ro, ti)-getpi(ro, ti))/(0.1*ro);	*/
	double deriv8 =  8.31*ti/27.0e-3;	
	
//	return (getpi(1.005*ro, ti)-getpi(ro, ti))/(0.005*ro);
	
//	double dpdro = dpdro_table.interpolate(ro, ti);


	double dpdro = (pi_table.interpolate(1.005*ro, ti)*1.005*ro*ti -
					pi_table.interpolate(0.995*ro, ti)*0.995*ro*ti) / 
				   (0.01 * ro);	

	return dpdro;
	//return deriv8;
}

double EOSTable::getdpdt(double ro, double ti, double te)
{
/*	double deriv1 = (getpi(ro, 1.001*ti)-getpi(ro, ti))/(0.001*ti);
	double deriv2 = (getpi(ro, 1.002*ti)-getpi(ro, ti))/(0.002*ti);
	double deriv3 = (getpi(ro, 1.005*ti)-getpi(ro, ti))/(0.005*ti);
	double deriv4 = (getpi(ro, 1.01*ti)-getpi(ro, ti))/(0.01*ti);
	double deriv5 = (getpi(ro, 1.02*ti)-getpi(ro, ti))/(0.02*ti);
	double deriv6 = (getpi(ro, 1.05*ti)-getpi(ro, ti))/(0.05*ti);
	double deriv7 = (getpi(ro, 1.1*ti)-getpi(ro, ti))/(0.1*ti); */
	double deriv8 =  ro*8.31/27.0e-3;	
	
//	return (getpi(ro, 1.1*ti)-getpi(ro, ti))/(0.1*ti);

//	return deriv8;

//	double dpdti = dpdti_table.interpolate(ro, ti);

	double dpdti = (pi_table.interpolate(ro, 1.05*ti)*ro*1.05*ti -
					pi_table.interpolate(ro, 0.95*ti)*ro*0.95*ti) / 
				   (0.1 * ti);	


    return dpdti;
}

double EOSTable::getdedro(double ro, double ti, double te)
{
/*	double deriv1 = (getei(1.001*ro, ti)-getei(ro, ti))/(0.001*ro);
	double deriv2 = (getei(1.002*ro, ti)-getei(ro, ti))/(0.002*ro);
	double deriv3 = (getei(1.005*ro, ti)-getei(ro, ti))/(0.005*ro);
	double deriv4 = (getei(1.01*ro, ti)-getei(ro, ti))/(0.01*ro);
	double deriv5 = (getei(1.02*ro, ti)-getei(ro, ti))/(0.02*ro);
	double deriv6 = (getei(1.05*ro, ti)-getei(ro, ti))/(0.05*ro);
	double deriv7 = (getei(1.1*ro, ti)-getei(ro, ti))/(0.1*ro);
	double deriv8 = 0.0;*/

	return 0.0;
	//return (getei(1.005*ro, ti)-getei(ro, ti))/(0.005*ro);
}

double EOSTable::getdedt(double ro, double ti, double te)
{
/*	double deriv1 = (getei(ro, 1.001*ti)-getei(ro, ti))/(0.001*ti);
	double deriv2 = (getei(ro, 1.002*ti)-getei(ro, ti))/(0.002*ti);
	double deriv3 = (getei(ro, 1.005*ti)-getei(ro, ti))/(0.005*ti);
	double deriv4 = (getei(ro, 1.01*ti)-getei(ro, ti))/(0.01*ti);
	double deriv5 = (getei(ro, 1.02*ti)-getei(ro, ti))/(0.02*ti);
	double deriv6 = (getei(ro, 1.05*ti)-getei(ro, ti))/(0.05*ti);
	double deriv7 = (getei(ro, 1.1*ti)-getei(ro, ti))/(0.1*ti);*/
	double deriv8 =  8.31/(5.0/3.0-1)/27.0e-3;
	
//    return (getei(ro, 1.1*ti)-getei(ro, ti))/(0.1*ti);
//	return deriv8;
		
//	double dedti = dedti_table.interpolate(ro, ti);

	double dedti = (ei_table.interpolate(ro, 1.05*ti)*1.05*ti -
					ei_table.interpolate(ro, 0.95*ti)*0.95*ti) / 
				   (0.1 * ti);	

	return dedti;
}


double EOSTable::getdpdro_rov_roE(double ro, double ti, double te, double v)
{	

	double deriv8 = 0.5*(5.0/3.0-1.0)*v*v;

	double dpde = getdpdt(ro, ti, te)/getdedt(ro, ti, te); 
	double dpdro = (getdpdro(ro, ti, te) + dpde/ro*(v*v/2.0 - getei(ro, ti) - ro*getdedro(ro, ti, te)));		
	return dpdro;
	
//	return 0.5*(5.0/3.0-1.0)*v*v;


}

double EOSTable::getdpdrov_ro_roE(double ro, double ti, double te, double v) 
{

	double deriv8 = -(5.0/3.0-1.0)*v;

	double dpdrov = (- getdpdt(ro, ti, te)/getdedt(ro, ti, te)*v/ro);

	return  dpdrov;

//	return -(5.0/3.0-1.0)*v;
}

double EOSTable::getdpdroE_ro_rov(double ro, double ti, double te, double v) 
{ 
	double deriv8 = (5.0/3.0-1.0);

	double dpdroE = (getdpdt(ro, ti, te)/getdedt(ro, ti, te)/ro);
	return dpdroE;

//	return (5.0/3.0-1.0);
}

///////////// Конец DEBUGа!///////////////////////////////////////////////

#include "defines.h"
#include "eosTableAu.h"

#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<iostream>

using namespace std;

// Таблицы

double EOSTableAu::getpi(double ro, double ti) {
 	double pi = pi_table.interpolate(ro, ti);
	// Критерий откола	
	if(pi*1.e-9 < -5.04 + 0.0017*ti - 1.49e-7*ti*ti)
		pi = 0.;
	return pi; // [Pa]
}

double EOSTableAu::getei(double ro, double ti) { return ei_table.interpolate(ro, ti);}	      // [J/kg]
double EOSTableAu::getci(double ro, double ti) { return ci_table.interpolate(ro, ti)*ro; }	  // ci [J/kg/K], ci/ro [J/m3/K]
double EOSTableAu::getphase(double ro, double ti) { return phase_table.interpolate(ro, ti); } 
double EOSTableAu::getmix(double ro, double ti) { return mix_table.interpolate(ro, ti); }
double EOSTableAu::getC(double ro, double ti, double te) { return C_table.interpolate(ro, ti); }   // [m/s]
double EOSTableAu::getEntropy(double ro, double ti) { return entropy_table.interpolate(ro, ti); }

// Обратная функция (температура)

double EOSTableAu::getti(double ro, double ei) {
	double ti = solve_ti(ro, ei, t_scale.getMin()*(1.+.01), t_scale.getMax()-1.0);
	return ti; // [K]
}

// Теплопроводность 

double EOSTableAu::getkappa(double ro, double ti, double te) {
	double eF=8.849e-22*pow(ro/ro0, 2./3.); // [kJ]
	double kB=1.38041e-26;				    // [kJ/K] 
	double TF = eF/kB;						// [K]	
    double mEl = 9.1083E-28;				// [g]	
	double vF2 = 3./5.*2.0*eF/mEl*1.e10;	// [sm2/s2], [kJ/g] = [10^10 sm2/s2]
	double v2 = sqrt(vF2*vF2/3./3.+kB*kB/mEl/mEl*te*te*1.e10*1.e10);		// [sm2/s2]
	double h_ = 1.054887e-37;				// [kJ*s]
	double beta = 1.;
	double nu_sol = 0.37e14 * ti/300. * getAlpha(ro, ti, te)/getAlpha(ro0, ti, 300.) + 
		            beta * eF/h_ * te*te/TF/TF; 
	double frac = 0., phase = fabs(getphase(ro, ti));
	if( phase == 2. || phase == 6.) 
		frac = getmix(ro, ti); 
	else if ( phase >= 3. && phase <= 5.) 
		frac = 0.;
	else 
		frac = 1.;
	double nu_liq = (3.3e14 + 1.5e11*ti) * (1. - frac);
	double nu_deg = nu_sol + nu_liq;
	double nu = 1./sqrt(1./(nu_deg*nu_deg/pow(ro/ro0, 1.3)/pow(ro/ro0, 1.3)) +
                        1./(eF*eF/h_/h_/pow(te/TF, 3.)*pow (ro/ro0, 4/3))); 
	/////////
	double _ce = getce(ro, te);
	/////////
	double kappa = getce(ro, te)*v2/nu*1.e-9; // 1. Теплоемкость ce здесь в [J/m3/K] = [10^-9 kJ/sm3/K] 
	// 2. Квадрат скорости v2 в [sm2/s2] 3. Итог -- kappa в [kJ/s/sm/K] = [10^5 J/m/s/K]
    return kappa*1.e5;		//[W/m/K]
}

double getnuWR(double ro, double ti, double te, double b, double Z) {return 1.;}

// Служебное

EOSTableAu::EOSTableAu(string dirName, int EOSFlag, double _ro0) {
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
	if(EOSFlag == 1) {
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
	if(EOSFlag !=1) {
		cout << "EOS error: probably table EOS of unknown type. EOSFlag should be 1 for this EOS."<< endl;
		exit(1);
	}
    v_scale.getFromFile("volume.in", dirName, nVScale, 0., 0., 1.0e-3);
	t_scale.getFromFile("temper.in", dirName, nTScale, tableTmin, tableTmax, 1.0e3);		
       ei_table.create("energy.tab", dirName, 1.0e6, &v_scale, &t_scale, EOSFlag);
	   pi_table.create("pressu.tab", dirName, 1.0e9, &v_scale, &t_scale, EOSFlag);
	   ci_table.create("hcapv.tab",  dirName, 1.0e3,   &v_scale, &t_scale, EOSFlag);
	    C_table.create("sndvls.tab", dirName, 1.0e3, &v_scale, &t_scale, EOSFlag);
	phase_table.create("physid.tab", dirName, 1.0, &v_scale, &t_scale, EOSFlag);
	  mix_table.create("PhsMix.tab", dirName, 1.0, &v_scale, &t_scale, EOSFlag); 
	MAX_T  = tableTmax; 
	MIN_T  = tableTmin;
	MAX_RO = 1.0/v_scale.getMin();
	MIN_RO = 1.0/v_scale.getMax();
	ro0 = _ro0;
	tableErrorFlag = 0;
	tableErrorCellNum = -1;
}


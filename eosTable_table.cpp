#include "defines.h"
#include "eosTable.h"
#include "eosTable_table.h"

//////////////////////
#include "eosIdeal.h"
//////////////////////

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <iostream>

using namespace std;

EOSTTable::EOSTTable()
{
	v_scale = 0;
	t_scale = 0;
	tab = 0;
}


EOSTTable::~EOSTTable()
{
	clear();
}


void EOSTTable::create(const char fName[], std::string dirName, double scaler, EOSTScale *_v_scale, EOSTScale *_t_scale, int EOSFlag)
{
	clear();

	v_scale = _v_scale;
	t_scale = _t_scale;

	int vSize = v_scale->getSize();
	int tSize = t_scale->getSize();

	double tmp;
	int i,j;

	// allocate memory

	tab = new double *[vSize];
	
	for(i=0; i < vSize; i++)
		tab[i] = new double[tSize];

	// read table from file

	char filename[_MAX_PATH];
	strcpy(filename, TABLE_FOLDER);
	strcat(filename, dirName.c_str());
	strcat(filename, "/");
	strcat(filename, fName);

	FILE* f=fopen(filename, "r");

	
	//////////////////////DEBUG!!!////////////////
/*	FILE *f1=fopen("pressure.dat", "w");
	fprintf(f1, "TITLE=\"Cold profiles\"\n");
	fprintf(f1, "VARIABLES=\"v\",\"Ro\",\"Pi\",\"P_an\"\n");*/
	//////////////////////////////////////////////

	if(!f)
	{
		cout << "Cannot open file " << filename << endl;
		exit(1);
	}

	///// Creating EOS tables////////

/*	EOSIdeal *eosIdeal = new EOSIdeal;

	FILE *etF=fopen("dedti_ideal.dat","w+");
	FILE *ptF=fopen("dpdti_ideal.dat","w+");
	FILE *proF=fopen("dpdro_ideal.dat","w+");

	FILE *p_roT_F=fopen("pressu_ideal_roT.dat","w+");
	FILE *e_roT_F=fopen("energy_ideal_roT.dat","w+");
*/
	
	FILE *pF = fopen("pressu_roT.dat","w+");
	FILE *eF = fopen("energy_T.dat","w+");

	/////////////////////////////////


	for(i=0; i<tSize; i++)
	{
		for(j=0; j<vSize; j++)
		{
			if(!fscanf(f, "%lf", &tmp))
			{
				printf("Error reading table EOS: %s\n", filename);
				perror("Fail");
				printf("V node number = %d\n", j);
				printf("T node number = %d\n", i);

				exit(1);
			}

			tab[j][i] = tmp*scaler;




			///// Creating EOS tables////////

			if(j%5==0)
			{
/*				fprintf(etF, "\n");
				fprintf(ptF, "\n");
				fprintf(proF, "\n");

				fprintf(p_roT_F, "\n");
				fprintf(e_roT_F, "\n");*/
		 		
				fprintf(pF, "\n");
				fprintf(eF, "\n");

			}

			double ro = 1.0/(*v_scale)[j];
			double ti  = (*t_scale)[i];			
			
/*			double dedti = 8.31 / (5.0/3.0 - 1.0) / 27.0e-3;
			double dpdro = 8.31 * ti / 27.0e-3;
			double dpdti = ro * 8.31 / 27.0e-3;

			fprintf(etF, "%e ", dedti);
			fprintf(ptF, "%e ", dpdti);
			fprintf(proF, "%e ", dpdro);

			double p_roT = 8.31/27.0e-3*1.0e-9;
			double e_roT = 8.31/ro/(5.0/3.0-1.0)/27.0e-3*1.0e-6;

			fprintf(p_roT_F, "%e ", p_roT);
			fprintf(e_roT_F, "%e ", e_roT);
*/			

			double p_roT = tab[j][i]/ro/ti/scaler;
			double e_T = tab[j][i]/ti/scaler;

			fprintf(pF, "%e ", p_roT);
			fprintf(eF, "%e ", e_T);
			
			/////////////////////////////////

		}
	}

	//////////// Fixing t=0.0001K nodes -- for euler BGK method //////////////
// Temporarily block it 
	/*if(EOSFlag == 1)
	{
		for(j=0; j<vSize; j++)
		{
			tab[j][0] = tab[j][1];
		}
	} */

	//////////////////////////////////////////////////

		///// Creating EOS tables////////
		
/*		fclose(etF);
		fclose(ptF);
		fclose(proF);
		fclose(p_roT_F);
		fclose(e_roT_F);

		
		delete eosIdeal;
*/
		fclose(eF);
		fclose(pF);

		/////////////////////////////////	

	fclose(f);
}

void EOSTTable::clear()
{
	if(tab)
	{
		for(int i=0; i < v_scale->getSize(); i++)
			delete[] tab[i];
	
		delete[] tab;
		
		v_scale = 0;
		t_scale = 0;
		tab = 0;
	}
}


double EOSTTable::interpolate(double ro, double ti)
{
	EOSTScale &vscale = *v_scale;
	EOSTScale &tscale = *t_scale;

	double v=1.0/ro;
	
	int i = vscale.getLeftIdx(v);
	
	if (i<0)
	{
		cout << "EOSTTable::interpolate() error: value out of range of table EOS" << endl <<
			    "ro = "  << ro << endl <<
				"ti = "  << ti << endl;
		exit(1);
	}

	double log_v       = log(v);
	double log_v1      = log(vscale[i]);
	double log_v2      = log(vscale[i+1]);
	double log_v_coeff = (log_v-log_v1) / (log_v2-log_v1);

	int j = tscale.getLeftIdx(ti);
	if(ti == tscale[0]) 
		j = 0;
	
	if (j<0)
	{
		cout << "EOSTTable::interpolate() error: value out of range of table EOS" << endl <<
			    "ro = "  << ro << endl <<
				"ti = "  << ti << endl;
		exit(1);
	}

	double log_t       = log(ti);
	double log_t1      = log(tscale[j]);
	double log_t2	   = log(tscale[j+1]);
	double log_t_coeff = (log_t-log_t1) / (log_t2-log_t1);

	double l_t_appr = tab[i][j]   + log_t_coeff*(tab[i][j+1]  - tab[i][j]);
	double r_t_appr = tab[i+1][j] + log_t_coeff*(tab[i+1][j+1]- tab[i+1][j]);
	
	/////////////////////
	double t1 = tab[i][j];
	double t2 = tab[i+1][j];
	double t3 = tab[i][j+1];
	double t4 = tab[i+1][j+1];

	double v1 = vscale[i];
	double v2 = vscale[i+1];
	double ti1 = tscale[j];
	double ti2 = tscale[j+1];

	/////////////////////

	return l_t_appr + log_v_coeff*(r_t_appr-l_t_appr);	
}

void EOSTTable::monotonize(int i, int j)
{
	tab[i+1][j] = tab[i][j];
}

bool EOSTTable::isMonotone(int i, int j)
{
	if(tab[i+1][j]<tab[i][j]) 
		return true;
	else 
		return false;
}

void EOSTTable::xchg(double *a, double *b)
{
	double temp = 0.0;
	temp = *a;
	*a = *b;
	*b = temp;
}

/* Заморозим тут старую версию для опытов. Вот кому теперь сказать 
спасибо, что Сиплюсплюс не проверяет выходов за границы массива? 
И вообще, надо брать версии от разработчиков хоть на чем, потому что 
они РАБОТАЮТ!!!

Это можно оставить как пример. Ошибка в том, что просто перепутаны 
местами параметры в вызове.

void EOSTTable::correctTable(int nTScale, int nVScale)
{
	for(int i=0; i<nTScale; i++)
	{
		xchg(&(tab[117][i]), &(tab[118][i]));
		xchg(&(tab[233][i]), &(tab[234][i]));
		xchg(&(tab[232][i]), &(tab[233][i]));
	}

	// Вадим, внимание! В EOSTScale::getFromFile() уже также 
	// производится коррекция! Я думаю, это и есть источник 
	// необъяснимой ошибки.
}

*/

void EOSTTable::correctTable(int nVScale, int nTScale)
{
	for(int i=0; i<nTScale; i++)
	{
		xchg(&(tab[117][i]), &(tab[118][i]));
		xchg(&(tab[233][i]), &(tab[234][i]));
		xchg(&(tab[232][i]), &(tab[233][i]));
	}

	// Вадим, внимание! В EOSTScale::getFromFile() уже также 
	// производится коррекция! Я думаю, это и есть источник 
	// необъяснимой ошибки.
}



int cond(double, double, double, double, double);
int cond(double tabVal, double d1, double d2, double d3, double d4){
	int c_more=0, c_less=0;
	if(tabVal>=d1)
		c_more++;
	else
		c_less++;
	if(tabVal>=d2)
		c_more++;
	else
		c_less++;
	if(tabVal>=d3)
		c_more++;
	else
		c_less++;
	if(tabVal>=d4)
		c_more++;
	else
		c_less++;
	if(c_more == c_less)
		return 1;
	else
		return 0;
}

double EOSTTable::calct(double ro, double tabVal){
	int iV = v_scale->getLeftIdx(1./ro);
	int nSize = t_scale->getSize();
	int j = 0;
	EOSTScale &vscale = *v_scale;
	EOSTScale &tscale = *t_scale;
	double log_v = log(1/ro);
	double log_v1      = log(vscale[iV]);
	double log_v2      = log(vscale[iV+1]);
	double log_v_coeff = (log_v-log_v1) / (log_v2-log_v1);
	double l_v_appr = 0., r_v_appr=0.;
	for(;;) {
		double l_v_appr = tab[iV][j]   + log_v_coeff*(tab[iV+1][j]  - tab[iV][j]);
		double r_v_appr = tab[iV][j+1] + log_v_coeff*(tab[iV+1][j+1]- tab[iV][j+1]);
		if(tabVal >= l_v_appr && tabVal < r_v_appr || (j == nSize-1))
			break;
		j++;
	} 
	//////////////////////////////////////////////////////////////////////////////////////////////////
	int jT = j;
	// Теперь мы знаем, что точка со значением tabVal находится в ячейке с iV, jT.
	// Осталось, зная точный объем (по плотности), выяснить, какой температуре
	// соответствует значение энергии tabVal, полученное билинейной интерполяцией.
/*	EOSTScale &vscale = *v_scale;
	EOSTScale &tscale = *t_scale;
	double log_v       = log(1/ro);
	double log_v1      = log(vscale[iV]);
	double log_v2      = log(vscale[iV+1]);
	double log_v_coeff = (log_v-log_v1) / (log_v2-log_v1);*/

	/*double l_t_appr = tab[i][j]   + log_t_coeff*(tab[i][j+1]  - tab[i][j]);
	double r_t_appr = tab[i+1][j] + log_t_coeff*(tab[i+1][j+1]- tab[i+1][j]);
	return l_t_appr + log_v_coeff*(r_t_appr-l_t_appr);	
	=>

	ei =              tab[i][j]   + log_t_coeff*(tab[i][j+1]  - tab[i][j]) +
	   + log_v_coeff* (tab[i+1][j] + log_t_coeff*(tab[i+1][j+1]- tab[i+1][j])   -    tab[i][j]   - log_t_coeff*(tab[i][j+1]  - tab[i][j]) )		

	ei =  tab[i][j] + log_v_coeff*(tab[i+1][j] - tab[i][j]) +  
	                  log_t_coeff*(tab[i][j+1]- tab[i][j] + log_v_coeff*(tab[i+1][j+1]- tab[i+1][j]) - tab[i][j+1]  + tab[i][j])   )
	
	log_t_coeff = (ei - tab[i][j] - log_v_coeff*(tab[i+1][j] - tab[i][j])) / ( (tab[i][j+1]- tab[i][j] + log_v_coeff*(tab[i+1][j+1]- tab[i+1][j]) - tab[i][j+1]  + tab[i][j])  ) */

	// Вот эту штуку он считает неправильно! Перепроверить.
	double log_t_coeff = (tabVal - tab[iV][jT] - log_v_coeff*(tab[iV+1][jT] - tab[iV][jT])  )  /  ( tab[iV][jT+1]- tab[iV][jT] + log_v_coeff*(tab[iV+1][jT+1]- tab[iV+1][jT] - tab[iV][jT+1]  + tab[iV][jT]) );
	// После перепроверки все сойдется.
	/*
	double log_t       = log(ti);
	double log_t1      = log(tscale[j]);
	double log_t2	   = log(tscale[j+1]);
	double log_t_coeff = (log_t-log_t1) / (log_t2-log_t1);
	=>
	log_t - log_t1  = log_t_coeff * (log_t2-log_t1)
	log_t			= log_t1 + log_t_coeff * (log_t2-log_t1)
	t = exp (   log_t1 + log_t_coeff * (log_t2-log_t1)   )
	находим, в общем, логарифм t, а затем находим t.	
	*/
	double log_t1      = log(tscale[jT]);
	double log_t2	   = log(tscale[jT+1]);
	double log_t       = log_t1 + log_t_coeff * (log_t2-log_t1);
	double t		   = exp (log_t);
	return t;
}


#ifndef EOSTABLE_H
#define EOSTABLE_H

#include "defines.h"
#include "eosAnalytic.h"
#include "eosTable_scale.h"
#include "eosTable_table.h"


// Папка, где хранятся файлы таблиц

//0.184E+00  [cc/g]   initial volume
//0.184E+04  [cc/g]   final volume
//
//#define EOS_MIN_VOLUME 0.184e-3
//#define EOS_MAX_VOLUME 0.184e1
//
//0.293E+00  10**3[K] initial temperature(energy, pressure)
//0.100E+03  10**3[K] final temperature(...)
//
//#define EOS_MIN_T 0.293e3
//#define EOS_MAX_T 0.1e6

class EOSTable : public EOSAnalytic
{
public:
	EOSTable();
	EOSTable(string dirName, int EOSFlag, double _ro0);

	EOSType getType() { return table; }

	double    getpi(double ro, double ti);
	double    getei(double ro, double ti);
	double    getti(double ro, double ei);
	double	  getkappa(double ro, double ti, double te);
	double        getC(double ro, double ti, double te);
	double	  getphase(double ro, double ti);
	double	    getmix(double ro, double ti);	
	double     getnuWR(double ro, double ti, double te, double b, double Z=3.0);
	complex<double> geteps(double ro, double ti, double te, double Z=3.0);

	double getEntropy(double ro, double ti);

	double getdpdro  (double ro, double ti, double te);
	double getdpdt   (double ro, double ti, double te);
	double getdedro  (double ro, double ti, double te);
	double getdedt   (double ro, double ti, double te);

	double getdpdro_rov_roE(double ro, double ti, double te, double v);
	double getdpdrov_ro_roE(double ro, double ti, double te, double v);
	double getdpdroE_ro_rov(double ro, double ti, double te, double v);

	//////////////////////////////
	
	double	  getti_interp(double ro, double ei);
	void	  monotonizeTables(void);


private:

	EOSTScale v_scale;
	EOSTScale t_scale;

	EOSTTable pi_table;
	EOSTTable ei_table;
	EOSTTable C_table;
	EOSTTable entropy_table;
	EOSTTable phase_table;
	EOSTTable mix_table;

	int tableErrorFlag;
	int tableErrorCellNum;

//////////////////DEBUG ///////////////
// Пока только для идеального газа
	
	EOSTTable dpdti_table;
	EOSTTable dpdro_table;
	EOSTTable dedti_table;

///////////////////////////////////////

};


#endif

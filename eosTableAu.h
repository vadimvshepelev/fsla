#ifndef EOSTABLEAU_H
#define EOSTABLEAU_H

#include "defines.h"
#include "eosAnalyticAu.h"
#include "eosTable_scale.h"
#include "eosTable_table.h"

class EOSTableAu : public EOSAnalyticAu {
public:
	EOSTableAu();
	EOSTableAu(string dirName, int EOSFlag, double _ro0);
	EOSType getType() { return tableNi; }
	double    getpi(double ro, double ti);
	double    getei(double ro, double ti);
	double    getti(double ro, double ei);
	double	  getci(double ro, double ti);
	double	  getkappa(double ro, double ti, double te);
	double        getC(double ro, double ti, double te);
	double	  getphase(double ro, double ti);
	double	    getmix(double ro, double ti);	
	double getEntropy(double ro, double ti);
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
	EOSTTable ci_table;
	int tableErrorFlag;
	int tableErrorCellNum;
};

#endif
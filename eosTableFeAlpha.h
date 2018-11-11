#ifndef _EOS_TABLE_FE_ALPHA
#define _EOS_TABLE_FE_ALPHA

#include "defines.h"
#include "eosAnalyticFe.h"
#include "eosTable_scale.h"
#include "eosTable_table.h"

class EOSTableFeAlpha : public EOSAnalyticFe
{
public:
	EOSTableFeAlpha();
	EOSTableFeAlpha(string dirName, int EOSFlag, double _ro0);

	EOSType getType() { return tableFeAlpha; }

	double    getpi(double ro, double ti);
	double    getei(double ro, double ti);
	
	double    getti(double ro, double ei);
	
	double	  getkappa(double ro, double ti, double te, double Z);
	double        getC(double ro, double ti, double te, double v=0);

	double	  getphase(double ro, double ti);
	double	    getmix(double ro, double ti);	

	double     getnuWR(double ro, double ti, double te, double b, double Z);

	double getEntropy(double ro, double ti);

	////DEBUG//////////////
	
	double	  getti_interp(double ro, double ei);
	void	  monotonizeTables(void);

	///////////////////////

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
};



#endif
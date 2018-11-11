
#ifndef EOSTABLE_TABLE_H
#define EOSTABLE_TABLE_H


#include "eosTable_scale.h"


class EOSTTable
{
public:

	EOSTTable();
	~EOSTTable();

	void	create(char *fName, string dirName, double scaler, EOSTScale *_v_scale, EOSTScale *_t_scale, int EOSFlag);
	void	clear();
	double	interpolate(double ro, double ti);

	bool	isMonotone(int i, int j);
	void    monotonize(int i, int j);

	void	correctTable(int nVScale, int nVTscale); // Только для EOSFlag 1, исправляет следствия немонотонного файла объемов.
	//////////DEBUG/////////////
	double	calct(double ro, double tabVal);
	////////////////////////////

private:

	EOSTScale *v_scale;
	EOSTScale *t_scale;

	double **tab;

	void xchg(double *a, double *b);
};


#endif

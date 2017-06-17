#include "defines.h"
#include "eosTable_scale.h"
#include "eosTable.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>


EOSTScale::EOSTScale()
{
	scale = 0;
	nSize = 0;

	minV = 0.0;
	maxV = 0.0;
}


EOSTScale::~EOSTScale()
{
	clear();
}


void EOSTScale::create(int _nSize, double _minV, double _maxV, double scaler)
{
	clear();

	nSize = _nSize;
	minV  = _minV;
	maxV  = _maxV;

	scale = new double[nSize];

	double log_step = (log(maxV)-log(minV))/(nSize-1);
	for(int i=0; i<nSize; i++)
		scale[i] = exp(log(minV) + log_step*i) * scaler;
}

void EOSTScale::getFromFile(char* fName, char* dirName, int _nSize, double _minV, double _maxV, double scaler)
{
	clear();

	nSize = _nSize;
	

	scale = new double[nSize];
	
	double tmp=0.0;

	char filename[_MAX_PATH];
	strcpy(filename, TABLE_FOLDER);
	strcat(filename, dirName);	
	strcat(filename, "/");
	strcat(filename, fName);

	FILE* f=fopen(filename, "r");
	printf("Processing: %s\n", filename);
	
	for(int j=0; j<getSize(); j++) {
		if(!fscanf(f, "%lf", &tmp)) {
			printf("Error reading table EOS: %s\n", filename);
			exit(1);
		}
		scale[j] = tmp*scaler;
	}

	minV  = scale[0];
	maxV  = scale[getSize()-1];
}

void EOSTScale::clear()
{
	if(scale) delete[] scale;

	scale = 0;
	nSize = 0;

	minV = 0.0;
	maxV = 0.0;
}


int EOSTScale::getLeftIdx(double t)
{
	//////////////////Важный момент! Меняю второе неравенство на "<="/////////
	if(!(t >= scale[0] && t /*<*/ <= scale[nSize-1]))
		return -1;

	//////////////////Важный момент! Меняю второе неравенство на "<="/////////
	assert(t >= scale[0] && t <= scale[nSize-1]);
	int i = 0;
	while(t > scale[i] && i < nSize-1) i++;
	return i-1;
}

void EOSTScale::xchg(double *a, double *b)
{
	double temp = 0.0;
	temp = *a;
	*a = *b;
	*b = temp;
}


#ifndef NODE_H
#define NODE_H


#include "_matrix4.h"


class Node
{
public:

	double x;
	double v;
	    
	double ro;
	double dm;

	double te;
	double ti;

	double pe;
	double pi;
	double p;

	double ee;
	double ei;
	double e;

	double ce;
	double ci;

	double C;	
	double Alphaei;
	double kappa;

	double ne;
	double Z;

	// flow

	Vector4 W;
	Vector4 F;

	// temporary values for double buffering

	double ti_temp;
	double te_temp;

	Vector4 W_temp;
};


#endif

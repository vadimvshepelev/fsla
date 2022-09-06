#ifndef _C1DBCS_H_
#define _C1DBCS_H_

#include"C1DField.h"


class C1DBoundary {
public:
	virtual void set(C1DField &fld) = 0;
};


class CLeftTransmissiveBoundary : public C1DBoundary {
public:
	void set(C1DField &fld);
};


class CRightTransmissiveBoundary : public C1DBoundary {
public:
	void set(C1DField &fld);
};

// Жалко очень этот класс. Такой классный он! Не буду стирать, оставлю. Надеюсь, что по осмыслении органично встроится он в код.
class C1DBCs {
public:
	C1DBoundary& BCL;
	C1DBoundary& BCR;
	C1DBCs(C1DBoundary& _BCL, C1DBoundary& _BCR) : BCL(_BCL), BCR(_BCR) {}
	void set(C1DField &fld) { BCL.set(fld); BCR.set(fld); }
};


#endif
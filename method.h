
#ifndef METHOD_H
#define METHOD_H


#include "matterState.h"
#include "eos.h"


/*
   Базовый абстрактный класс для метода расчета на этапе
   гидродинамики. Реализует общий интерфейс и вспомогательные
   функции.
*/

class CMethod
{
public:

	CMethod(EOS *_eos);
	virtual void	matter2Flow(CField &ms) = 0;
	virtual void	advanceFlow(CField &ms, double tau) = 0;
	virtual void	advanceFlowVacuum(CField &ms, double tau) = 0;
	virtual void	flow2Matter(CField &ms, double tau) = 0;
	virtual void	flow2MatterVacuum(CField &ms, CField &ms_temp, double t, double tau) = 0;
	virtual	Matrix4 getOmega(Node &n) = 0;
	virtual Matrix4 getOmegaInv(Node &n) = 0;
	virtual Vector4 getLambda(Node &n) = 0;
	virtual void	createGrid(CField &ms) = 0;
	virtual void	deleteGrid() = 0;
	virtual void	averageNode(Node &n1, Node &n2, Node &nav) = 0;
	virtual void	updateNode(Node &n) = 0;
	virtual void	fillLambda(Vector4 &Fm1, Vector4 &F0, Vector4 &Fp1, Vector4 &Fp2,
					   Vector4 &L, double step) = 0;
	virtual void	fillLambdaComponent(int i, double lambda, double criteria, double curant) = 0;
	virtual Vector4 calcFlux(Matrix4 Omega, Matrix4 OmegaInv, Vector4 W_j_m, Vector4 W_j, Vector4 W_j_p, Vector4 W_j_p_p) = 0;
	virtual Vector4 calcApprRPFlux(Vector4 W_j, Vector4 W_j_p, Node &n) = 0;

	Vector4 La, Lb, Lg, Ld;		// Лямбды * веса для 4ех точек сетки

	double *vGrid;
	double *X;

protected:
	
	// Vector4 La, Lb, Lg, Ld;		// Лямбды * веса для 4ех точек сетки

	EOS &eos;
};


#endif
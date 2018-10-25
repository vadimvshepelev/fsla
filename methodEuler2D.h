#ifndef METHODEULER2D_H
#define METHODEULER2D_H

#include "method.h"
#include "_matrix4.h"

class CMethodEuler2D : public CMethod
{
public:
	CMethodEuler2D(EOS *_eos) : CMethod(_eos) { ///DEBUG////
											  
											  ro_min=2000.0;
											  ro_max=3000.0;
											  ti_min=10000.0;
											  ti_max=10000.0;

											  ////////////
											  /*assert(eos.getType() == ideal);*/ }

	void	matter2Flow(CField &ms);
	void	advanceFlow(CField &ms, double tau);
	void	advanceFlowVacuum(CField &ms, double tau);
	void	flow2Matter(CField &ms, double tau);
	void	flow2MatterVacuum(CField &ms, CField &ms_temp, double t, double tau);
	void	averageNode(Node &n1, Node &n2, Node &nav) {}
	void	updateNode(Node &n) {}
	void	fillLambda(Vector4 &Fm1, Vector4 &F0, Vector4 &Fp1, Vector4 &Fp2, Vector4 &L, double step) {}
	void	fillLambdaComponent(int i, double lambda, double criteria, double curant) {}
	void	createGrid(CField &ms);
	void	deleteGrid();
	Matrix4 getOmega(Node &n) {return Matrix4::ZERO;}
	Matrix4 getOmegaInv(Node &n) {return Matrix4::ZERO;}
	Vector4 getLambda(Node &n) {return Vector4::ZERO;}
	Vector4 calcFlux(Matrix4 Omega, Matrix4 OmegaInv, Vector4 W_j_m, Vector4 W_j, Vector4 W_j_p, Vector4 W_j_p_p);
	Vector4 calcApprRPFlux(Vector4 W_j, Vector4 W_j_p, Node &n);


	/////DEBUG////////
	double ro_min;
	double ro_max;
	double ti_min;
	double ti_max;
	//////////////////

private:


};



#endif

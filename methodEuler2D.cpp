
#include "methodEuler2D.h"
#include <math.h>

void CMethodEuler2D::matter2Flow(MatterState &ms)
{
	// We should take away this all in the nearest refactoring! Too inconvinient.
}

void CMethodEuler2D::advanceFlow(MatterState &ms, double tau)
{
	
}

void CMethodEuler2D::advanceFlowVacuum(MatterState &ms, double tau)
{
	
}

void CMethodEuler2D::flow2Matter(MatterState &ms, double tau) {}

void CMethodEuler2D::flow2MatterVacuum(MatterState &ms, MatterState &ms_temp, double t, double tau) {}

void CMethodEuler2D::createGrid(MatterState &ms){}
void CMethodEuler2D::deleteGrid(){}

Vector4 CMethodEuler2D::calcFlux(Matrix4 Omega, Matrix4 OmegaInv, Vector4 W_j_m, Vector4 W_j, Vector4 W_j_p, Vector4 W_j_p_p){
	return Vector4::ZERO;
}

Vector4 CMethodEuler2D::calcApprRPFlux(Vector4 W_j, Vector4 W_j_p, Node &n)
{
	return Vector4::ZERO;
}

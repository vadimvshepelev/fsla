#include "..\\solver.h"
#include <vector>
#include <fstream>

// ENO3-G - method
// Riemann solver based on exact solution + ENO reconstruction of 3rd order
void CSolver::calcHydroStageENO3G(double t, double tau) {
	cout << "calcHydroStageENO3G(): ";
	EOS &eos = task.getEOS(); 
	const double gamma = eos.getGamma();
	double E=0.; 
	unsigned int i=0, j=0;
	const unsigned int nSize = ms.getSize();
	const double h=ms[1].x-ms[0].x;
	Vector4 L = Vector4::ZERO, R = Vector4::ZERO, D = Vector4::ZERO, V = Vector4::ZERO;	
	// Специальные массивы для хранения переменных и потоков, с дополнительными фиктивными ячейками
	const unsigned int nGhostCells = 3;	
	const unsigned int mini = nGhostCells, maxi = nGhostCells+nSize;
	std::vector<Vector4> U, Up, Um;
	std::vector<Vector4> F;
	for(i=0; i<nSize+nGhostCells+nGhostCells; i++) {
		F.push_back(Vector4::ZERO); Up.push_back(Vector4::ZERO); Um.push_back(Vector4::ZERO);
		if(i>=nGhostCells && i<nGhostCells+nSize)
			U.push_back(ms[i-nGhostCells].W);
		else
			U.push_back(Vector4::ZERO);
	}
	// Transmissive b.c.s
	U[2] = U[mini]; U[1] = U[mini]; U[0] = U[mini];  
	U[maxi] = U[maxi-1]; U[maxi+1] = U[maxi-1]; U[maxi+2] = U[maxi-1];
	// Possible ENO-3 stencils
	double rm2r[] = { 1./3., -7./6.,  11./6.,     0.,     0.}, 
		   rm1r[] = {    0., -1./6.,  5./6.,   1./3.,     0.},
			r0r[] = {    0.,     0.,  1./3.,   5./6., -1./6.},
		   rm2l[] = {-1./6.,  5./6.,  1./3.,      0.,     0.},
		   rm1l[] = {    0.,  1./3.,  5./6.,  -1./6.,     0.},
		    r0l[] = {    0.,     0.,  11./6., -7./6.,  1./3.};
	ofstream ofs;
	ofs.open("reconstruction.dat", ios::out);
	double x_for_writing=0., step = h/20; 
	// ENO 3rd order reconstruction
	unsigned int stencil0=0, stencil1=0, stencil2=0;
	for(i=mini; i<maxi; i++) {		
		Vector4 diffPlus = U[i+1]-U[i], diffPlusPlus = U[i+2]-U[i+1], diffMinus = U[i]-U[i-1], diffMinusMinus = U[i-1]-U[i-2];
		stencil0 = i; stencil1 = i; stencil2 = i;
		double x0 = (double)(i-nGhostCells)*h;
		if( fabs(diffMinus[0])<=fabs(diffPlus[0]) ) {
			stencil0 = i-1;
			if( fabs(diffMinus[0] - diffMinusMinus[0])<=fabs(diffPlus[0] - diffMinus[0]) )
				stencil0 = i-2;		
		} else {
			if( fabs(diffPlusPlus[0] - diffPlus[0]) <= fabs(diffPlus[0] - diffMinus[0]) )
				stencil0 = i;
			else
				stencil0 = i-1;
		}
		if( fabs(diffMinus[1])<=fabs(diffPlus[1]) ) {
			stencil1 = i-1;
			if( fabs(diffMinus[1] - diffMinusMinus[1])<=fabs(diffPlus[1] - diffMinus[1]) )
				stencil1 = i-2;		
		} else {
			if( fabs(diffPlusPlus[1] - diffPlus[1]) <= fabs(diffPlus[1] - diffMinus[1]) )
				stencil1 = i;
			else
				stencil1 = i-1;
		}
		if( fabs(diffMinus[2])<=fabs(diffPlus[2]) ) {
			stencil2 = i-1;
			if( fabs(diffMinus[2] - diffMinusMinus[2])<=fabs(diffPlus[2] - diffMinus[2]) )
				stencil2 = i-2;	
		} else {
			if( fabs(diffPlusPlus[2] - diffPlus[2]) <= fabs(diffPlus[2] - diffMinus[2]) )
				stencil2 = i;
			else
				stencil2 = i-1;
		}
		// Interface values

		// TODO: Использовать "принцип отражения", как говорит А.В. Конюхов
		// Смысл -- чтобы подойти к точке с другой стороны, мысленно переворачиваем
		// профиль относительно границы i+1/2 -- имеем тот же набор средних, ту же задачу,
		// но в другом направлении. 

		double _up = 0., _um = 0., xim12 = (double)(i-mini)*h;



		if(i == 32) {

			double q =  0.;

		}



	/*	if(stencil0 == i-2) {
			Up[i][0] = rm2r[0]*U[i-2][0] + rm2r[1]*U[i-1][0] + rm2r[2]*U[i][0] + rm2r[3]*U[i+1][0] + rm2r[4]*U[i+2][0];
			Um[i][0] = rm2l[0]*U[i-2][0] + rm2l[1]*U[i-1][0] + rm2l[2]*U[i][0] + rm2l[3]*U[i+1][0] + rm2l[4]*U[i+2][0];
			_up = calcInterpolationPolynomialDerivative3(xim12-2.*h, xim12-h, xim12, xim12+h, 0., U[i-2][0]*h, (U[i-2][0]+U[i-1][0])*h, (U[i-2][0]+U[i-1][0]+U[i][0])*h, xim12+h);
			_um = calcInterpolationPolynomialDerivative3(xim12-2.*h, xim12-h, xim12, xim12+h, 0., U[i-2][0]*h, (U[i-2][0]+U[i-1][0])*h, (U[i-2][0]+U[i-1][0]+U[i][0])*h, xim12);
			if(fabs(_up-Up[i][0]) > 1.e-5 || fabs(_um-Um[i][0]) > 1.e-5){
				double qqqqqq=0.;
			}
		} else if (stencil0 == i-1) {
			Up[i][0] = rm1r[0]*U[i-2][0] + rm1r[1]*U[i-1][0] + rm1r[2]*U[i][0] + rm1r[3]*U[i+1][0] + rm1r[4]*U[i+2][0];
			Um[i][0] = rm1l[0]*U[i-2][0] + rm1l[1]*U[i-1][0] + rm1l[2]*U[i][0] + rm1l[3]*U[i+1][0] + rm1l[4]*U[i+2][0];
		} else if (stencil0 == i) {
			Up[i][0] = r0r[0]*U[i-2][0] + r0r[1]*U[i-1][0] + r0r[2]*U[i][0] + r0r[3]*U[i+1][0] + r0r[4]*U[i+2][0];
			Um[i][0] = r0l[0]*U[i-2][0] + r0l[1]*U[i-1][0] + r0l[2]*U[i][0] + r0l[3]*U[i+1][0] + r0l[4]*U[i+2][0];
		}
		if(stencil1 == i-2) {
			Up[i][1] = rm2r[0]*U[i-2][1] + rm2r[1]*U[i-1][1] + rm2r[2]*U[i][1] + rm2r[3]*U[i+1][1] + rm2r[4]*U[i+2][1];
			Um[i][1] = rm2l[0]*U[i-2][1] + rm2l[1]*U[i-1][1] + rm2l[2]*U[i][1] + rm2l[3]*U[i+1][1] + rm2l[4]*U[i+2][1];
		} else if (stencil1 == i-1) {
			Up[i][1] = rm1r[0]*U[i-2][1] + rm1r[1]*U[i-1][1] + rm1r[2]*U[i][1] + rm1r[3]*U[i+1][1] + rm1r[4]*U[i+2][1];
			Um[i][1] = rm1l[0]*U[i-2][1] + rm1l[1]*U[i-1][1] + rm1l[2]*U[i][1] + rm1l[3]*U[i+1][1] + rm1l[4]*U[i+2][1];
		} else if (stencil1 == i) {
			Up[i][1] = r0r[0]*U[i-2][1] + r0r[1]*U[i-1][1] + r0r[2]*U[i][1] + r0r[3]*U[i+1][1] + r0r[4]*U[i+2][1];
			Um[i][1] = r0l[0]*U[i-2][1] + r0l[1]*U[i-1][1] + r0l[2]*U[i][1] + r0l[3]*U[i+1][1] + r0l[4]*U[i+2][1];
		}
		if(stencil2 == i-2) {
			Up[i][2] = rm2r[0]*U[i-2][2] + rm2r[1]*U[i-1][2] + rm2r[2]*U[i][2] + rm2r[3]*U[i+1][2] + rm2r[4]*U[i+2][2];
			Um[i][2] = rm2l[0]*U[i-2][2] + rm2l[1]*U[i-1][2] + rm2l[2]*U[i][2] + rm2l[3]*U[i+1][2] + rm2l[4]*U[i+2][2];
		} else if (stencil2 == i-1) {
			Up[i][2] = rm1r[0]*U[i-2][2] + rm1r[1]*U[i-1][2] + rm1r[2]*U[i][2] + rm1r[3]*U[i+1][2] + rm1r[4]*U[i+2][2];
			Um[i][2] = rm1l[0]*U[i-2][2] + rm1l[1]*U[i-1][2] + rm1l[2]*U[i][2] + rm1l[3]*U[i+1][2] + rm1l[4]*U[i+2][2];
		} else if (stencil2 == i) {
			Up[i][2] = r0r[0]*U[i-2][2] + r0r[1]*U[i-1][2] + r0r[2]*U[i][2] + r0r[3]*U[i+1][2] + r0r[4]*U[i+2][2];
			Um[i][2] = r0l[0]*U[i-2][2] + r0l[1]*U[i-1][2] + r0l[2]*U[i][2] + r0l[3]*U[i+1][2] + r0l[4]*U[i+2][2];
		} */
		

		for(int j=0; j<20; j++) {
			x_for_writing = x0 + (double)(j)*step;
			ofs << x_for_writing << " ";
			if (stencil0 == i-2)
				ofs << calcInterpolationPolynomialDerivative3(x0-2.*h,x0-h,x0,x0+h,0.,U[i-2][0]*h,(U[i-2][0]+U[i-1][0])*h,(U[i-2][0]+U[i-1][0]+U[i][0])*h,x_for_writing) << endl;
			else if (stencil0 == i-1)
				ofs << calcInterpolationPolynomialDerivative3(x0-h,x0,x0+h,x0+2.*h,0.,U[i-1][0]*h,(U[i-1][0]+U[i][0])*h,(U[i-1][0]+U[i][0]+U[i+1][0])*h,x_for_writing) << endl;
			else if (stencil0 == i)
				ofs << calcInterpolationPolynomialDerivative3(x0,x0+h,x0+2.*h,x0+3.*h,0.,U[i][0]*h,(U[i][0]+U[i+1][0])*h,(U[i][0]+U[i+1][0]+U[i+2][0])*h,x_for_writing) << endl;
			else
				exit(1);
		}


		double _up_formula = 0., _um_formula = 0.;
		if(stencil0 == i-2) {
			Up[i][0] =  1./3.*U[i-2][0] - 7./6.*U[i-1][0] + 11./6.*U[i][0];
			Um[i][0] = -1./6.*U[i-2][0] + 5./6.*U[i-1][0] +  1./3.*U[i][0];
			//_up_formula = Up[i][0];
			//_um_formula = Um[i][0];
			//_up = calcInterpolationPolynomialDerivative3(xim12-2.*h, xim12-h, xim12, xim12+h, 0., U[i-2][0]*h, (U[i-2][0]+U[i-1][0])*h, (U[i-2][0]+U[i-1][0]+U[i][0])*h, xim12+h);
			//_um = calcInterpolationPolynomialDerivative3(xim12-2.*h, xim12-h, xim12, xim12+h, 0., U[i-2][0]*h, (U[i-2][0]+U[i-1][0])*h, (U[i-2][0]+U[i-1][0]+U[i][0])*h, xim12);
			//if(fabs(_up-_up_formula)>1.e-5)	{
			//	double qqqqqq=0.;
			//}
			//if(fabs(_um-_um_formula)>1.e-5)	{
			//	double qqqqqq=0.;
			//}			
		} else if (stencil0 == i-1) {
			Up[i][0] = -1./6.*U[i-1][0] + 5./6.*U[i][0] + 1./3*U[i+1][0];
			Um[i][0] =  1./3.*U[i-1][0] + 5./6.*U[i][0] - 1./6.*U[i+1][0];
			//_up_formula = Up[i][0];
			//_um_formula = Um[i][0];
			//_up = calcInterpolationPolynomialDerivative3(xim12-h, xim12, xim12+h, xim12+2.*h, 0., U[i-1][0]*h, (U[i-1][0]+U[i][0])*h, (U[i-1][0]+U[i][0]+U[i+1][0])*h, xim12+h);
			//_um = calcInterpolationPolynomialDerivative3(xim12-h, xim12, xim12+h, xim12+2.*h, 0., U[i-1][0]*h, (U[i-1][0]+U[i][0])*h, (U[i-1][0]+U[i][0]+U[i+1][0])*h, xim12);
			//if(fabs(_up-_up_formula)>1.e-5)	{
			//	double qqqqqq=0.;
			//}
			//if(fabs(_um-_um_formula)>1.e-5)	{
		    //		double qqqqqq=0.;
			//}
		} else if (stencil0 == i) {
			Up[i][0] =  1./3.*U[i][0] + 5./6.*U[i+1][0] - 1./6.*U[i+2][0];
			Um[i][0] = 11./6.*U[i][0] - 7./6.*U[i+1][0] + 1./3.*U[i+2][0];	

			//double up_i_minus_2 = 1./3.*U[i-2][0] - 7./6.*U[i-1][0] + 11./6.*U[i][0];
			//double um_i_minus_2 = -1./6.*U[i-2][0] + 5./6.*U[i-1][0] +  1./3.*U[i][0];
			//double up_i_minus_1 = -1./6.*U[i-1][0] + 5./6.*U[i][0] + 1./3*U[i+1][0];
			//double um_i_minus_1 = 1./3.*U[i-1][0] + 5./6.*U[i][0] - 1./6.*U[i+1][0];



			//_up_formula = Up[i][0];
			//_um_formula = Um[i][0];
			//_up = calcInterpolationPolynomialDerivative3(xim12, xim12+h, xim12+2.*h, xim12+3.*h, 0., U[i][0]*h, (U[i][0]+U[i+1][0])*h, (U[i][0]+U[i+1][0]+U[i+2][0])*h, xim12+h);
			//_um = calcInterpolationPolynomialDerivative3(xim12, xim12+h, xim12+2.*h, xim12+3.*h, 0., U[i][0]*h, (U[i][0]+U[i+1][0])*h, (U[i][0]+U[i+1][0]+U[i+2][0])*h, xim12);
			//if(fabs(_up-_up_formula)>1.e-5)	{
			//	double qqqqqq=0.;
			//}
			//if(fabs(_um-_um_formula)>1.e-5)	{
			//	double qqqqqq=0.;
			//}				
		}
		if(stencil1 == i-2) {
			Up[i][1] =  1./3.*U[i-2][1] - 7./6.*U[i-1][1] + 11./6.*U[i][1];
			Um[i][1] = -1./6.*U[i-2][1] + 5./6.*U[i-1][1] +  1./3.*U[i][1];
			//_up_formula = Up[i][1];
			//_um_formula = Um[i][1];
			//_up = calcInterpolationPolynomialDerivative3(xim12-2.*h, xim12-h, xim12, xim12+h, 0., U[i-2][1]*h, (U[i-2][1]+U[i-1][1])*h, (U[i-2][1]+U[i-1][1]+U[i][1])*h, xim12+h);
			//_um = calcInterpolationPolynomialDerivative3(xim12-2.*h, xim12-h, xim12, xim12+h, 0., U[i-2][1]*h, (U[i-2][1]+U[i-1][1])*h, (U[i-2][1]+U[i-1][1]+U[i][1])*h, xim12);
			//if(fabs(_up-_up_formula)>1.e-5)	{
			//	double qqqqqq=0.;
			//}
			//if(fabs(_um-_um_formula)>1.e-5)	{
			//	double qqqqqq=0.;
			//}
		} else if (stencil1 == i-1) {
			Up[i][1] = -1./6.*U[i-1][1] + 5./6.*U[i][1] + 1./3.*U[i+1][1];
			Um[i][1] =  1./3.*U[i-1][1] + 5./6.*U[i][1] - 1./6.*U[i+1][1];
			//_up_formula = Up[i][1];
			//_um_formula = Um[i][1];
			//_up = calcInterpolationPolynomialDerivative3(xim12-h, xim12, xim12+h, xim12+2.*h, 0., U[i-1][1]*h, (U[i-1][1]+U[i][1])*h, (U[i-1][1]+U[i][1]+U[i+1][1])*h, xim12+h);
			//_um = calcInterpolationPolynomialDerivative3(xim12-h, xim12, xim12+h, xim12+2.*h, 0., U[i-1][1]*h, (U[i-1][1]+U[i][1])*h, (U[i-1][1]+U[i][1]+U[i+1][1])*h, xim12);
			//if(fabs(_up-_up_formula)>1.e-5)	{
			//	double qqqqqq=0.;
			//}
			//if(fabs(_um-_um_formula)>1.e-5)	{
		    //		double qqqqqq=0.;
			//}
		} else if (stencil1 == i) {
			Up[i][1] =  1./3.*U[i][1] + 5./6.*U[i+1][1] - 1./6.*U[i+2][1];
			Um[i][1] = 11./6.*U[i][1] - 7./6.*U[i+1][1] + 1./3.*U[i+2][1];
			//_up_formula = Up[i][1];
			//_um_formula = Um[i][1];
			//_up = calcInterpolationPolynomialDerivative3(xim12, xim12+h, xim12+2.*h, xim12+3.*h, 0., U[i][1]*h, (U[i][1]+U[i+1][1])*h, (U[i][1]+U[i+1][1]+U[i+2][1])*h, xim12+h);
			//_um = calcInterpolationPolynomialDerivative3(xim12, xim12+h, xim12+2.*h, xim12+3.*h, 0., U[i][1]*h, (U[i][1]+U[i+1][1])*h, (U[i][1]+U[i+1][1]+U[i+2][1])*h, xim12);
			//if(fabs(_up-_up_formula)>1.e-5)	{
			//	double qqqqqq=0.;
			//}
			//if(fabs(_um-_um_formula)>1.e-5)	{
		    //		double qqqqqq=0.;
			//}		
		}
		if(stencil2 == i-2) {
			Up[i][2] =  1./3.*U[i-2][2] - 7./6.*U[i-1][2] + 11./6.*U[i][2];
			Um[i][2] = -1./6.*U[i-2][2] + 5./6.*U[i-1][2] +  1./3.*U[i][2];
			//_up_formula = Up[i][2];
			//_um_formula = Um[i][2];
			//_up = calcInterpolationPolynomialDerivative3(xim12-2.*h, xim12-h, xim12, xim12+h, 0., U[i-2][2]*h, (U[i-2][2]+U[i-1][2])*h, (U[i-2][2]+U[i-1][2]+U[i][2])*h, xim12+h);
			//_um = calcInterpolationPolynomialDerivative3(xim12-2.*h, xim12-h, xim12, xim12+h, 0., U[i-2][2]*h, (U[i-2][2]+U[i-1][2])*h, (U[i-2][2]+U[i-1][2]+U[i][2])*h, xim12);
			//if(fabs(_up-_up_formula)>1.e-5)	{
			//	double qqqqqq=0.;
			//}
			//if(fabs(_um-_um_formula)>1.e-5)	{
			//	double qqqqqq=0.;
			//}
		} else if (stencil2 == i-1) {
			Up[i][2] = -1./6.*U[i-1][2] + 5./6.*U[i][2] + 1./3.*U[i+1][2];
			Um[i][2] =  1./3.*U[i-1][2] + 5./6.*U[i][2] - 1./6.*U[i+1][2];
			//_up_formula = Up[i][2];
			//_um_formula = Um[i][2];
			//_up = calcInterpolationPolynomialDerivative3(xim12-h, xim12, xim12+h, xim12+2.*h, 0., U[i-1][2]*h, (U[i-1][2]+U[i][2])*h, (U[i-1][2]+U[i][2]+U[i+1][2])*h, xim12+h);
			//_um = calcInterpolationPolynomialDerivative3(xim12-h, xim12, xim12+h, xim12+2.*h, 0., U[i-1][2]*h, (U[i-1][2]+U[i][2])*h, (U[i-1][2]+U[i][2]+U[i+1][2])*h, xim12);
			//if(fabs(_up-_up_formula)>1.e-5)	{
		    //		double qqqqqq=0.;
			//}
			//if(fabs(_um-_um_formula)>1.e-5)	{
			//	double qqqqqq=0.;
			//}
		} else if (stencil2 == i) {
			Up[i][2] =  1./3.*U[i][2] + 5./6.*U[i+1][2] - 1./6.*U[i+2][2];
			Um[i][2] = 11./6.*U[i][2] - 7./6.*U[i+1][2] + 1./3.*U[i+2][2];
			//_up_formula = Up[i][2];
			//_um_formula = Um[i][2];
			//_up = calcInterpolationPolynomialDerivative3(xim12, xim12+h, xim12+2.*h, xim12+3.*h, 0., U[i][2]*h, (U[i][2]+U[i+1][2])*h, (U[i][2]+U[i+1][2]+U[i+2][2])*h, xim12+h);
			//_um = calcInterpolationPolynomialDerivative3(xim12, xim12+h, xim12+2.*h, xim12+3.*h, 0., U[i][2]*h, (U[i][2]+U[i+1][2])*h, (U[i][2]+U[i+1][2]+U[i+2][2])*h, xim12);
			//if(fabs(_up-_up_formula)>1.e-5)	{
			//	double qqqqqq=0.;
			//}
			//if(fabs(_um-_um_formula)>1.e-5)	{
			//	double qqqqqq=0.;
			//}		
		}
	}
	// Boundary nodes U+[0] and U-[nSize]
	Up[mini-1] = Up[mini]; Um[maxi] = Um[maxi-1];
	/////
/*	string _fName = string(OUTPUT_FOLDER) + "\\" + "test-reconstruction-2.dat";
	ofstream ofs; ofs.open(_fName, ios::out);
	ofs << "TITLE = \"ENO-3 reconstruction test on the cells' borders, f(x) = x^2\"" << endl;
	ofs << "VARIABLES = \"x\", \"f(x)\", \"sq(x)\"" << endl;
	for(i=mini; i<maxi; i++) {
		double _x = ms[i-nGhostCells].x;
		double _x_left = _x;
		double _x_right = ms[i-nGhostCells+1].x;
		double _left_rec = Um[i][0];
		double _right_rec = Up[i][0];
		double _left_abs = _x_left*_x_left;
		double _right_abs = _x_right*_x_right;
			ofs << _x << " " << _left_rec << " " << _left_abs << endl;
	}
*/



	for(i=0; i<=nSize; i++) {		
//		F[i] = calcGodunovFlux(Up[i-1][0], Up[i-1][1], Up[i-1][2], Um[i][0], Um[i][1], Um[i][2]);	
		
//		ms[i].F = calcGodunovFlux(Up[i-1+nGhostCells][0], Up[i-1+nGhostCells][1], Up[i-1+nGhostCells][2], 
//			                      Um[i+nGhostCells][0], Um[i+nGhostCells][1], Um[i+nGhostCells][2]);
		ms[i].F = calcRoeFlux(Up[i-1+nGhostCells][0], Up[i-1+nGhostCells][1], Up[i-1+nGhostCells][2], 
			                      Um[i+nGhostCells][0], Um[i+nGhostCells][1], Um[i+nGhostCells][2]);

//		ms[i].F = calcRoeFlux(U[i-1+nGhostCells][0], U[i-1+nGhostCells][1], U[i-1+nGhostCells][2], 
//			                      U[i+nGhostCells][0], U[i+nGhostCells][1], U[i+nGhostCells][2]);
		
	}
	ofs.close();
	// Main cycle
	for(i=0; i<nSize; i++) {				
		ms[i].W_temp = ms[i].W - tau/h*(ms[i+1].F-ms[i].F);     //(F[nGhostCells+i+1]-F[nGhostCells+i]);
	}
	for(i=0; i<nSize; i++) {
		// Обновляем консервативные переменные
		Node& n=ms[i]; 
		n.W[0] = n.W_temp[0];
		n.W[1] = n.W_temp[1];
		n.W[2] = n.W_temp[2];
		// Обновляем все переменные
		n.ro = n.W[0];
		n.v  = n.W[1]/n.ro;
		n.e  = n.W[2]/n.ro - 0.5*n.v*n.v;
		n.ti = eos.getti(n.ro, n.e);
		n.p  = eos.getpi(n.ro, n.ti);
		n.C  = eos.getC(n.ro, n.ti, n.te);
	}
	cout << " done!" << endl;
}

void CSolver::calcHydroStageENO2G(double t, double tau) {
	cout << "calcHydroStageENO2G(): ";
	EOS &eos = task.getEOS(); 
	const double gamma = eos.getGamma();
	double E=0.; 
	unsigned int i=0, j=0;
	const unsigned int nSize = ms.getSize();
	const double h=ms[1].x-ms[0].x;
	Vector4 L = Vector4::ZERO, R = Vector4::ZERO, D = Vector4::ZERO, V = Vector4::ZERO;	
	// Специальные массивы для хранения переменных и потоков, с дополнительными фиктивными ячейками
	const unsigned int nGhostCells = 3;	
	const unsigned int mini = nGhostCells, maxi = nGhostCells+nSize;
	std::vector<Vector4> U, Up, Um;
	std::vector<Vector4> F;
	for(i=0; i<nSize+nGhostCells+nGhostCells; i++) {
		F.push_back(Vector4::ZERO); Up.push_back(Vector4::ZERO); Um.push_back(Vector4::ZERO);
		if(i>=nGhostCells && i<nGhostCells+nSize)
			U.push_back(ms[i-nGhostCells].W);
		else
			U.push_back(Vector4::ZERO);
	}
	// Transmissive b.c.s
	U[2] = U[mini]; U[1] = U[mini]; U[0] = U[mini];  
	U[maxi] = U[maxi-1]; U[maxi+1] = U[maxi-1]; U[maxi+2] = U[maxi-1];
	// Possible ENO-2 stencils
	double  rm1r[] = {-1./2., 3./2.,  0.},       // {0.,     3./2.,  -1./2}, 
		    r0r[]  = {0.,     1./2.,  1./2.},
			rm1l[] = {1./2.,  1./2.,  0.},
			r0l[]  = {0.,     3./2., -1./2.};
	ofstream ofs;
	ofs.open("reconstruction.dat", ios::out);
	double x_for_writing=0., step = h/20; 
	// ENO 2nd order reconstruction
	unsigned int stencil0=0, stencil1=0, stencil2=0;
	for(i=mini; i<maxi; i++) {		
		Vector4 diffPlus = U[i+1]-U[i], diffMinus = U[i]-U[i-1];
		stencil0 = i; stencil1 = i; stencil2 = i;
		double x0 = (double)(i-nGhostCells)*h;
		// Choose stencil
		if( fabs(diffMinus[0])<=fabs(diffPlus[0]) ) stencil0 = i-1; else stencil0 = i;
		if( fabs(diffMinus[1])<=fabs(diffPlus[1]) ) stencil1 = i-1; else stencil1 = i;
		if( fabs(diffMinus[2])<=fabs(diffPlus[2]) ) stencil2 = i-1; else stencil2 = i;
		// Calculate ENO-2 approximations themselves


		if(i==32) {

			double qq = 0.;

		}


	   if (stencil0 == i-1) {
			Up[i][0] = rm1r[0]*U[i-1][0] + rm1r[1]*U[i][0] + rm1r[2]*U[i+1][0];
			Um[i][0] = rm1l[0]*U[i-1][0] + rm1l[1]*U[i][0] + rm1l[2]*U[i+1][0];
		} else if (stencil0 == i) {
			Up[i][0] = r0r[0]*U[i-1][0] + r0r[1]*U[i][0] + r0r[2]*U[i+1][0];
			Um[i][0] = r0l[0]*U[i-1][0] + r0l[1]*U[i][0] + r0l[2]*U[i+1][0];
		}
		if (stencil1 == i-1) {
			Up[i][1] = rm1r[0]*U[i-1][1] + rm1r[1]*U[i][1] + rm1r[2]*U[i+1][1];
			Um[i][1] = rm1l[0]*U[i-1][1] + rm1l[1]*U[i][1] + rm1l[2]*U[i+1][1];
		} else if (stencil1 == i) {
			Up[i][1] = r0r[0]*U[i-1][1] + r0r[1]*U[i][1] + r0r[2]*U[i+1][1];
			Um[i][1] = r0l[0]*U[i-1][1] + r0l[1]*U[i][1] + r0l[2]*U[i+1][1];
		}
	    if (stencil2 == i-1) {
			Up[i][2] = rm1r[0]*U[i-1][2] + rm1r[1]*U[i][2] + rm1r[2]*U[i+1][2];
			Um[i][2] = rm1l[0]*U[i-1][2] + rm1l[1]*U[i][2] + rm1l[2]*U[i+1][2];
		} else if (stencil2 == i) {
			Up[i][2] = r0r[0]*U[i-1][2] + r0r[1]*U[i][2] + r0r[2]*U[i+1][2];
			Um[i][2] = r0l[0]*U[i-1][2] + r0l[1]*U[i][2] + r0l[2]*U[i+1][2];
		}		
	
		 //Minmod slope limiter
/*		Up[i] = U[i]+0.5*calcMinmodSlopeModified(U[i]-U[i-1], U[i+1]-U[i]);
		Um[i] = U[i]-0.5*calcMinmodSlopeModified(U[i]-U[i-1], U[i+1]-U[i]);*/		
	}

	// Boundary nodes U+[0] and U-[nSize]
	Up[mini-1] = Up[mini]; Um[maxi] = Um[maxi-1];
	



	// Sign property test
	for(i=mini; i<=maxi; i++) {
		Vector4 deltaRec = Um[i+1]-Up[i], deltaVal = U[i+1]-U[i];
		for(int j=0; j<3; j++) {
			if (deltaRec[j]*deltaVal[j] < -1.e6) {
				cout << "Warning!!! ENO \'sign\' rule violated in node " << i << ", component " << j << endl;
			}
		}
	}

	// Calculating intercell fluxes
	for(i=0; i<=nSize; i++) {		
//		F[i] = calcGodunovFlux(Up[i-1][0], Up[i-1][1], Up[i-1][2], Um[i][0], Um[i][1], Um[i][2]);	
		
//		ms[i].F = calcGodunovFlux(Up[i-1+nGhostCells][0], Up[i-1+nGhostCells][1], Up[i-1+nGhostCells][2], 
//			                      Um[i+nGhostCells][0], Um[i+nGhostCells][1], Um[i+nGhostCells][2]);
		ms[i].F = calcRoeFlux(Up[i-1+nGhostCells][0], Up[i-1+nGhostCells][1], Up[i-1+nGhostCells][2], 
			                      Um[i+nGhostCells][0], Um[i+nGhostCells][1], Um[i+nGhostCells][2]);

//		ms[i].F = calcRoeFlux(U[i-1+nGhostCells][0], U[i-1+nGhostCells][1], U[i-1+nGhostCells][2], 
//			                      U[i+nGhostCells][0], U[i+nGhostCells][1], U[i+nGhostCells][2]);
		
	}
	ofs.close();
	// Main cycle
	for(i=0; i<nSize; i++) {				
		ms[i].W_temp = ms[i].W - tau/h*(ms[i+1].F-ms[i].F);     //(F[nGhostCells+i+1]-F[nGhostCells+i]);
	}
	for(i=0; i<nSize; i++) {
		// Обновляем консервативные переменные
		Node& n=ms[i]; 
		n.W[0] = n.W_temp[0];
		n.W[1] = n.W_temp[1];
		n.W[2] = n.W_temp[2];
		// Обновляем все переменные
		n.ro = n.W[0];
		n.v  = n.W[1]/n.ro;
		n.e  = n.W[2]/n.ro - 0.5*n.v*n.v;
		n.ti = eos.getti(n.ro, n.e);
		n.p  = eos.getpi(n.ro, n.ti);
		n.C  = eos.getC(n.ro, n.ti, n.te);
	}
	cout << " done!" << endl;
}


double CSolver::calcInterpolationPolynomialDerivative3(double xim32, double xim12, double xip12, double xip32, double fim32, double fim12, double fip12, double fip32, double x){
	// Calculating Lagrange interpolation polynomial of 3rd order for 4 control points	
	double pim32 = (3*x*x - 2*(xim12+xip12+xip32)*x + (xim12*xip12+xip12*xip32+xip32*xim12)) / ((xim32-xim12)*(xim32-xip12)*(xim32-xip32));
	double pim12 = (3*x*x - 2*(xim32+xip12+xip32)*x + (xim32*xip12+xip12*xip32+xip32*xim32)) / ((xim12-xim32)*(xim12-xip12)*(xim12-xip32));
	double pip12 = (3*x*x - 2*(xim32+xim12+xip32)*x + (xim32*xim12+xim12*xip32+xip32*xim32)) / ((xip12-xim32)*(xip12-xim12)*(xip12-xip32));
	double pip32 = (3*x*x - 2*(xim32+xim12+xip12)*x + (xim32*xim12+xim12*xip12+xip12*xim32)) / ((xip32-xim32)*(xip32-xim12)*(xip32-xip12));
	double p = fim32*pim32 + fim12*pim12 + fip12*pip12 + fip32*pip32;
	return p;
}
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
	double rm1[] = {11./6., -7./6.,  1./3}, 
		    r0[] = { 1./3.,  5./6., -1./6.},
			r1[] = {-1./6.,  5./6.,  1./3.},
			r2[] = { 1./3., -7./6., 11./6.};


	ofstream ofs;
	ofs.open("reconstruction.dat", ios::out);
	double x_for_writing=0., step = h/20; 
	double TV = 0.;
	for(i=mini+1; i<maxi; i++) {		
		TV += fabs(U[i][0]-U[i-1][0]);
	}







	// ENO 3rd order reconstruction
	unsigned int stencil0=0, stencil1=0, stencil2=0;
	for(i=mini; i<maxi; i++) {		
		Vector4 diffPlus = U[i+1]-U[i], diffPlusPlus = U[i+2]-U[i+1], diffMinus = U[i]-U[i-1], diffMinusMinus = U[i-1]-U[i-2];
		stencil0 = i; stencil1 = i; stencil2 = i;
		double x0 = (double)(i-nGhostCells)*h;



		if(i==34) {
		
			double qq = 0.;
			double qqq = calcInterpolationPolynomialDerivative3(.3, .4, .5, .6, .3, .4, .5, .55, 2.);
			double qqqq = calcInterpolationPolynomialDerivative3(.3, .4, .5, .6, 0., .1, .2, .25, 2.);
			double qqqqq = 0.;


		}
		


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

		// Поток Лакса-Фридрихса -- самый грубый и легкий в обращении, его попробовать

	/*	if(stencil0 == i-2) {
			Up[i][0] = U[i-2][0]*r2[0] + U[i-1][0]*r2[1] + U[i][0]*r2[2];
			Um[i][0] = U[i-2][0]*r1[0] + U[i-1][0]*r1[1] + U[i][0]*r1[2];
		} else if (stencil0 == i-1) {
			Up[i][0] = U[i-1][0]*r1[0] + U[i][0]*r1[1] + U[i+1][0]*r1[2];
			Um[i][0] = U[i-1][0]*r0[0] + U[i][0]*r0[1] + U[i+1][0]*r0[2];
		} else if (stencil0 == i) {
			Up[i][0] = U[i][0]*r0[0] + U[i+1][0]*r0[1] + U[i+2][0]*r0[2];
			Um[i][0] = U[i][0]*rm1[0] + U[i+1][0]*rm1[1] + U[i+2][0]*rm1[2];
		}
		if(stencil1 == i-2) {
			Up[i][1] = U[i-2][1]*r2[0] + U[i-1][1]*r2[1] + U[i][1]*r2[2];
			Um[i][1] = U[i-2][1]*r1[0] + U[i-1][1]*r1[1] + U[i][1]*r1[2];
		} else if (stencil1 == i-1) {
			Up[i][1] = U[i-1][1]*r1[0] + U[i][1]*r1[1] + U[i+1][1]*r1[2];
			Um[i][1] = U[i-1][1]*r0[0] + U[i][1]*r0[1] + U[i+1][1]*r0[2];
		} else if (stencil1 == i) {
			Up[i][1] = U[i][1]*r0[0] + U[i+1][1]*r0[1] + U[i+2][1]*r0[2];
			Um[i][1] = U[i][1]*rm1[0] + U[i+1][1]*rm1[1] + U[i+2][1]*rm1[2];
		}
		if(stencil2 == i-2) {
			Up[i][2] = U[i-2][2]*r2[0] + U[i-1][2]*r2[1] + U[i][2]*r2[2];
			Um[i][2] = U[i-2][2]*r1[0] + U[i-1][2]*r1[1] + U[i][2]*r1[2];
		} else if (stencil2 == i-1) {
			Up[i][2] = U[i-1][2]*r1[0] + U[i][2]*r1[1] + U[i+1][2]*r1[2];
			Um[i][2] = U[i-1][2]*r0[0] + U[i][2]*r0[1] + U[i+1][2]*r0[2];
		} else if (stencil2 == i) {
			Up[i][2] = U[i][2]*r0[0] + U[i+1][2]*r0[1] + U[i+2][2]*r0[2];
			Um[i][2] = U[i][2]*rm1[0] + U[i+1][2]*rm1[1] + U[i+2][2]*rm1[2];
		} */		
		if(i==34) {
			double q = 0.;
		}
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


		double _up = 0., _um = 0., _up_formula = 0., _um_formula = 0., xim12 = (double)(i-mini)*h;
		if(stencil0 == i-2) {
			Up[i][0] =  1./3.*U[i-2][0] - 7./6.*U[i-1][0] + 11./6.*U[i][0];
			Um[i][0] = -1./6.*U[i-2][0] + 5./6.*U[i-1][0] +  1./3.*U[i][0];
			_up_formula = Up[i][0];
			_um_formula = Um[i][0];
			_up = calcInterpolationPolynomialDerivative3(xim12-2.*h, xim12-h, xim12, xim12+h, 0., U[i-2][0]*h, (U[i-2][0]+U[i-1][0])*h, (U[i-2][0]+U[i-1][0]+U[i][0])*h, xim12+h);
			_um = calcInterpolationPolynomialDerivative3(xim12-2.*h, xim12-h, xim12, xim12+h, 0., U[i-2][0]*h, (U[i-2][0]+U[i-1][0])*h, (U[i-2][0]+U[i-1][0]+U[i][0])*h, xim12);
			if(fabs(_up-_up_formula)>1.e-5)	{
				double qqqqqq=0.;
			}
			if(fabs(_um-_um_formula)>1.e-5)	{
				double qqqqqq=0.;
			}			
		} else if (stencil0 == i-1) {
			Up[i][0] = -1./6.*U[i-1][0] + 5./6.*U[i][0] + 1./3*U[i+1][0];
			Um[i][0] =  1./3.*U[i-1][0] + 5./6.*U[i][0] - 1./6.*U[i+1][0];
			_up_formula = Up[i][0];
			_um_formula = Um[i][0];
			_up = calcInterpolationPolynomialDerivative3(xim12-h, xim12, xim12+h, xim12+2.*h, 0., U[i-1][0]*h, (U[i-1][0]+U[i][0])*h, (U[i-1][0]+U[i][0]+U[i+1][0])*h, xim12+h);
			_um = calcInterpolationPolynomialDerivative3(xim12-h, xim12, xim12+h, xim12+2.*h, 0., U[i-1][0]*h, (U[i-1][0]+U[i][0])*h, (U[i-1][0]+U[i][0]+U[i+1][0])*h, xim12);
			if(fabs(_up-_up_formula)>1.e-5)	{
				double qqqqqq=0.;
			}
			if(fabs(_um-_um_formula)>1.e-5)	{
				double qqqqqq=0.;
			}
		} else if (stencil0 == i) {
			Up[i][0] =  1./3.*U[i][0] + 5./6.*U[i+1][0] - 1./6.*U[i+2][0];
			Um[i][0] = 11./6.*U[i][0] - 7./6.*U[i+1][0] + 1./3.*U[i+2][0];	

			double up_i_minus_2 = 1./3.*U[i-2][0] - 7./6.*U[i-1][0] + 11./6.*U[i][0];
			double um_i_minus_2 = -1./6.*U[i-2][0] + 5./6.*U[i-1][0] +  1./3.*U[i][0];
			double up_i_minus_1 = -1./6.*U[i-1][0] + 5./6.*U[i][0] + 1./3*U[i+1][0];
			double um_i_minus_1 = 1./3.*U[i-1][0] + 5./6.*U[i][0] - 1./6.*U[i+1][0];



			_up_formula = Up[i][0];
			_um_formula = Um[i][0];
			_up = calcInterpolationPolynomialDerivative3(xim12, xim12+h, xim12+2.*h, xim12+3.*h, 0., U[i][0]*h, (U[i][0]+U[i+1][0])*h, (U[i][0]+U[i+1][0]+U[i+2][0])*h, xim12+h);
			_um = calcInterpolationPolynomialDerivative3(xim12, xim12+h, xim12+2.*h, xim12+3.*h, 0., U[i][0]*h, (U[i][0]+U[i+1][0])*h, (U[i][0]+U[i+1][0]+U[i+2][0])*h, xim12);
			if(fabs(_up-_up_formula)>1.e-5)	{
				double qqqqqq=0.;
			}
			if(fabs(_um-_um_formula)>1.e-5)	{
				double qqqqqq=0.;
			}				
		}
		if(stencil1 == i-2) {
			Up[i][1] =  1./3.*U[i-2][1] - 7./6.*U[i-1][1] + 11./6.*U[i][1];
			Um[i][1] = -1./6.*U[i-2][1] + 5./6.*U[i-1][1] +  1./3.*U[i][1];
			_up_formula = Up[i][1];
			_um_formula = Um[i][1];
			_up = calcInterpolationPolynomialDerivative3(xim12-2.*h, xim12-h, xim12, xim12+h, 0., U[i-2][1]*h, (U[i-2][1]+U[i-1][1])*h, (U[i-2][1]+U[i-1][1]+U[i][1])*h, xim12+h);
			_um = calcInterpolationPolynomialDerivative3(xim12-2.*h, xim12-h, xim12, xim12+h, 0., U[i-2][1]*h, (U[i-2][1]+U[i-1][1])*h, (U[i-2][1]+U[i-1][1]+U[i][1])*h, xim12);
			if(fabs(_up-_up_formula)>1.e-5)	{
				double qqqqqq=0.;
			}
			if(fabs(_um-_um_formula)>1.e-5)	{
				double qqqqqq=0.;
			}
		} else if (stencil1 == i-1) {
			Up[i][1] = -1./6.*U[i-1][1] + 5./6.*U[i][1] + 1./3.*U[i+1][1];
			Um[i][1] =  1./3.*U[i-1][1] + 5./6.*U[i][1] - 1./6.*U[i+1][1];
			_up_formula = Up[i][1];
			_um_formula = Um[i][1];
			_up = calcInterpolationPolynomialDerivative3(xim12-h, xim12, xim12+h, xim12+2.*h, 0., U[i-1][1]*h, (U[i-1][1]+U[i][1])*h, (U[i-1][1]+U[i][1]+U[i+1][1])*h, xim12+h);
			_um = calcInterpolationPolynomialDerivative3(xim12-h, xim12, xim12+h, xim12+2.*h, 0., U[i-1][1]*h, (U[i-1][1]+U[i][1])*h, (U[i-1][1]+U[i][1]+U[i+1][1])*h, xim12);
			if(fabs(_up-_up_formula)>1.e-5)	{
				double qqqqqq=0.;
			}
			if(fabs(_um-_um_formula)>1.e-5)	{
				double qqqqqq=0.;
			}
		} else if (stencil1 == i) {
			Up[i][1] =  1./3.*U[i][1] + 5./6.*U[i+1][1] - 1./6.*U[i+2][1];
			Um[i][1] = 11./6.*U[i][1] - 7./6.*U[i+1][1] + 1./3.*U[i+2][1];
			_up_formula = Up[i][1];
			_um_formula = Um[i][1];
			_up = calcInterpolationPolynomialDerivative3(xim12, xim12+h, xim12+2.*h, xim12+3.*h, 0., U[i][1]*h, (U[i][1]+U[i+1][1])*h, (U[i][1]+U[i+1][1]+U[i+2][1])*h, xim12+h);
			_um = calcInterpolationPolynomialDerivative3(xim12, xim12+h, xim12+2.*h, xim12+3.*h, 0., U[i][1]*h, (U[i][1]+U[i+1][1])*h, (U[i][1]+U[i+1][1]+U[i+2][1])*h, xim12);
			if(fabs(_up-_up_formula)>1.e-5)	{
				double qqqqqq=0.;
			}
			if(fabs(_um-_um_formula)>1.e-5)	{
				double qqqqqq=0.;
			}		
		}
		if(stencil2 == i-2) {
			Up[i][2] =  1./3.*U[i-2][2] - 7./6.*U[i-1][2] + 11./6.*U[i][2];
			Um[i][2] = -1./6.*U[i-2][2] + 5./6.*U[i-1][2] +  1./3.*U[i][2];
			_up_formula = Up[i][2];
			_um_formula = Um[i][2];
			_up = calcInterpolationPolynomialDerivative3(xim12-2.*h, xim12-h, xim12, xim12+h, 0., U[i-2][2]*h, (U[i-2][2]+U[i-1][2])*h, (U[i-2][2]+U[i-1][2]+U[i][2])*h, xim12+h);
			_um = calcInterpolationPolynomialDerivative3(xim12-2.*h, xim12-h, xim12, xim12+h, 0., U[i-2][2]*h, (U[i-2][2]+U[i-1][2])*h, (U[i-2][2]+U[i-1][2]+U[i][2])*h, xim12);
			if(fabs(_up-_up_formula)>1.e-5)	{
				double qqqqqq=0.;
			}
			if(fabs(_um-_um_formula)>1.e-5)	{
				double qqqqqq=0.;
			}
		} else if (stencil2 == i-1) {
			Up[i][2] = -1./6.*U[i-1][2] + 5./6.*U[i][2] + 1./3.*U[i+1][2];
			Um[i][2] =  1./3.*U[i-1][2] + 5./6.*U[i][2] - 1./6.*U[i+1][2];
			_up_formula = Up[i][2];
			_um_formula = Um[i][2];
			_up = calcInterpolationPolynomialDerivative3(xim12-h, xim12, xim12+h, xim12+2.*h, 0., U[i-1][2]*h, (U[i-1][2]+U[i][2])*h, (U[i-1][2]+U[i][2]+U[i+1][2])*h, xim12+h);
			_um = calcInterpolationPolynomialDerivative3(xim12-h, xim12, xim12+h, xim12+2.*h, 0., U[i-1][2]*h, (U[i-1][2]+U[i][2])*h, (U[i-1][2]+U[i][2]+U[i+1][2])*h, xim12);
			if(fabs(_up-_up_formula)>1.e-5)	{
				double qqqqqq=0.;
			}
			if(fabs(_um-_um_formula)>1.e-5)	{
				double qqqqqq=0.;
			}
		} else if (stencil2 == i) {
			Up[i][2] =  1./3.*U[i][2] + 5./6.*U[i+1][2] - 1./6.*U[i+2][2];
			Um[i][2] = 11./6.*U[i][2] - 7./6.*U[i+1][2] + 1./3.*U[i+2][2];
			_up_formula = Up[i][2];
			_um_formula = Um[i][2];
			_up = calcInterpolationPolynomialDerivative3(xim12, xim12+h, xim12+2.*h, xim12+3.*h, 0., U[i][2]*h, (U[i][2]+U[i+1][2])*h, (U[i][2]+U[i+1][2]+U[i+2][2])*h, xim12+h);
			_um = calcInterpolationPolynomialDerivative3(xim12, xim12+h, xim12+2.*h, xim12+3.*h, 0., U[i][2]*h, (U[i][2]+U[i+1][2])*h, (U[i][2]+U[i+1][2]+U[i+2][2])*h, xim12);
			if(fabs(_up-_up_formula)>1.e-5)	{
				double qqqqqq=0.;
			}
			if(fabs(_um-_um_formula)>1.e-5)	{
				double qqqqqq=0.;
			}		
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
	
	ofs.close();*/
	/////












	// Calculating intercell fluxes
	for(i=0; i<=nSize; i++) {
		
		if(i==32) {
			double gg = 9.;
		}
		
//		F[i] = calcGodunovFlux(Up[i-1][0], Up[i-1][1], Up[i-1][2], Um[i][0], Um[i][1], Um[i][2]);	

		ms[i].F = calcGodunovFlux(Up[i-1+nGhostCells][0], Up[i-1+nGhostCells][1], Up[i-1+nGhostCells][2], 
			                      Um[i+nGhostCells][0], Um[i+nGhostCells][1], Um[i+nGhostCells][2]);	

		// i=3 -- нормальный поток, распад разрыва
		// i=4 +		ms.nodes[4].F	{x=0.00000000000000000 y=0.10000000000000001 z=0.00000000000000000 ...}	Vector4
/*		roEL	0.25000000000000006	double ??? чем, интересно, roEL[4] отличается от roEL[5], [6], итд?
		
		roEL[4] = Up[4-1+3][2] = Up[6][2] -- stencil0==stencil1==stencil2==i==6
		Посмотреть также на Up[6][1] и Up[6][0]

		roER	0.25000000000000006	double
		roL	0.12500000000000000	double
		roR	0.12500000000000000	double
		rouL	0.00000000000000000	double
		rouR	0.00000000000000000	double
		*/
		// i=5 +		ms.nodes[5].F	{x=6.5566384441430873e-018 y=0.10000000000000002 z=1.8358587643600652e-017 ...}	Vector4
/*		roEL	0.25000000000000011	double
		roER	0.25000000000000006	double
		roL	0.12500000000000000	double
		roR	0.12500000000000000	double
		rouL	0.00000000000000000	double
		rouR	0.00000000000000000	double
		*/
		// i=6 +		ms.nodes[6].F	{x=6.5566384441430873e-018 y=0.10000000000000002 z=1.8358587643600652e-017 ...}	Vector4
/* 		roEL	0.25000000000000011	double
		roER	0.25000000000000006	double
		roL	0.12500000000000000	double
		roR	0.12500000000000000	double
		rouL	0.00000000000000000	double
		rouR	0.00000000000000000	double
		*/

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


	TV = 0.;
	for(i=1; i<nSize; i++) {		
		TV += fabs(ms[i].ro-ms[i-1].ro);
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
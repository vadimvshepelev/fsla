#include "..\\solver.h"
#include <vector>
#include <fstream>

// ENO3-G - method
// Riemann solver based on exact solution + ENO reconstruction of 3rd order
void CSolver::calcHydroStageENO3G(double t, double tau) {
	cout << "calcHydroStageENO3G(): ";
	EOSOld &eos = task.getEOS(); 
	const double gamma = eos.getGamma();
	double E=0.; 
	int i=0, nDimCounter=0;
	const int nSize = ms.getSize();
	const double h=ms[1].x-ms[0].x;
	Vector4 L = Vector4::ZERO, R = Vector4::ZERO, D = Vector4::ZERO, V = Vector4::ZERO;	
	// Специальные массивы для хранения переменных и потоков, с дополнительными фиктивными ячейками
	const int nGhostCells = 3;	
	const int mini = nGhostCells, maxi = nGhostCells+nSize;
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
	double _rm2r[] = { 1./3., -7./6.,  11./6.,     0.,     0.}, 
		   _rm1r[] = {    0., -1./6.,  5./6.,   1./3.,     0.},
			_r0r[] = {    0.,     0.,  1./3.,   5./6., -1./6.},
		   _rm2l[] = {-1./6.,  5./6.,  1./3.,      0.,     0.},
		   _rm1l[] = {    0.,  1./3.,  5./6.,  -1./6.,     0.},
		    _r0l[] = {    0.,     0.,  11./6., -7./6.,  1./3.};
	vector<double> rm2r(_rm2r, _rm2r + sizeof(_rm2r)/sizeof(_rm2r[0])),
				   rm1r(_rm1r, _rm1r + sizeof(_rm1r)/sizeof(_rm1r[0])),
				    r0r(_r0r, _r0r + sizeof(_r0r)/sizeof(_r0r[0])),
				   rm2l(_rm2l, _rm2l + sizeof(_rm2l)/sizeof(_rm2l[0])),
				   rm1l(_rm1l, _rm1l + sizeof(_rm1l)/sizeof(_rm1l[0])),
				   r0l(_r0l, _r0l + sizeof(_r0l)/sizeof(_r0l[0]));
	double  _rm1rENO2[] = {0., -1./2., 3./2.,     0., 0.},       // {0.,     3./2.,  -1./2}, 
		    _r0rENO2[]  = {0.,     0., 1./2.,  1./2., 0.},
			_rm1lENO2[] = {0.,  1./2., 1./2.,     0., 0.},
			_r0lENO2[]  = {0.,     0., 3./2., -1./2., 0.};
	vector<double> rm1rENO2(_rm1rENO2, _rm1rENO2 + sizeof(_rm1rENO2)/sizeof(_rm1rENO2[0])),
				   r0rENO2(_r0rENO2, _r0rENO2 + sizeof(_r0rENO2)/sizeof(_r0rENO2[0])),
				   rm1lENO2(_rm1lENO2, _rm1lENO2 + sizeof(_rm1lENO2)/sizeof(_rm1lENO2[0])),
				   r0lENO2(_r0lENO2, _r0lENO2 + sizeof(_r0lENO2)/sizeof(_r0lENO2[0]));

	vector<double> deltarm2rENO3(5), deltarm1rENO3(5), deltar0rENO3(5), deltarm2lENO3(5), deltarm1lENO3(5), deltar0lENO3(5);
    for(int j=0; j<5; j++) {
		deltarm2rENO3[j] = rm2r[j] - rm1rENO2[j],
		deltarm1rENO3[j] = rm1r[j] - rm1rENO2[j],
		deltar0rENO3[j]  = r0r[j] - r0rENO2[j],
		deltarm2lENO3[j] = rm2l[j] - rm1lENO2[j],
		deltarm1lENO3[j] = rm1l[j] - rm1lENO2[j],
		deltar0lENO3[j]  = r0l[j] - r0lENO2[j];
	}
	//Special arrays for finite differences for the density
	vector<double> diff1(nSize+nGhostCells+nGhostCells), diff2(nSize+nGhostCells+nGhostCells);
	for(i=mini; i<maxi; i++) 
		diff1[i]=U[i+1][0]-U[i][0];
	diff1[maxi]=0.;
	for(i=mini; i<maxi; i++) 
		diff2[i]=diff1[i+1]-diff1[i];
	diff2[maxi]=0.;
	


	// Just switching from ENO2 to ENO3 and back inside the code
	double flag=1., flag0 = 1., flag1 = 1., flag2 = 1.;
	int nimin2Counter = 0, nimin1Counter = 0, niCounter = 0;
	// DEBUG
	// В итоге здесь какие-то расхождения на графиках. Сделать четко, расписать подробно, что и где вычисляем,
	// при необходимости сделать вычисления руками в нескольких точках около вершины параболы, любая ошибка 
	// реконструкции на этих данных должны. 
	//
	// Разобраться с граничными условаями, почему они так жестко срезаются в ноль? Проверить алгоритм, корректные
	// ли значения на границах. Также есть несовпадение с точным решением (параболой) порядка процентов или десятых
	// долей процентов. Возможно, это и есть ошибка реконструкции. 
	//
	// Также нет совпадения между левыми и правыми приближениями в границах ячеек, а они как раз должны
	// быть одинаковыми, раз приближают один и тот же многочлен (параболу).

	// Testing recinstruction procedure
	// Let ro(x) = x^2-1.5*x. Then ENO-3 reconstruction will exactly represent the curve
	/*ofstream fcentres("test-ENO3-centres.dat", ios::out);
	ofstream fborders("test-ENO3-borders.dat", ios::out);
	for(i=0; i<maxi+nGhostCells; i++) {
		double xm = (double)(i-nGhostCells)*h;
		double xp = xm + h;
		double x = .5*(xm + xp);
		U[i][0] = 500.* 1./h * (xp*xp*xp/3. - 3.*xp*xp/4. - xm*xm*xm/3. + 3.*xm*xm/4.);		
	}*/
	// ENO 3rd order reconstruction
	for(i=mini; i<maxi; i++) {				
		Vector4 diffPlus = U[i+1]-U[i], diffPlusPlus = U[i+2]-U[i+1], diffMinus = U[i]-U[i-1], diffMinusMinus = U[i-1]-U[i-2];
		int stencil[] = {0, 0, 0};		
		Up[i] = 0.;
		Um[i] = 0.;




		if(i==3+32) {

			double qq = 0;
		}





		for(nDimCounter=0; nDimCounter<3; nDimCounter++) {
			if( fabs(diffMinus[nDimCounter])<=fabs(diffPlus[nDimCounter]) ) {
				stencil[nDimCounter] = i-1;
				if( fabs(diffMinus[nDimCounter] - diffMinusMinus[nDimCounter])<=fabs(diffPlus[nDimCounter] - diffMinus[nDimCounter]) ) {
					stencil[nDimCounter] = i-2;		
					if(nDimCounter==0) nimin2Counter++;
				} else
					if(nDimCounter==0) nimin1Counter++;
			} else {
				if( fabs(diffPlusPlus[nDimCounter] - diffPlus[nDimCounter]) <= fabs(diffPlus[nDimCounter] - diffMinus[nDimCounter]) || flag==0.) {
					stencil[nDimCounter] = i;
					if(nDimCounter==0) niCounter++;
				} else {
					stencil[nDimCounter] = i-1;
					if(nDimCounter==0) nimin1Counter++;
				}
			}
			// Interface values
			// TODO: Использовать "принцип отражения", как говорит А.В. Конюхов
			// Смысл -- чтобы подойти к точке с другой стороны, мысленно переворачиваем
			// профиль относительно границы i+1/2 -- имеем тот же набор средних, ту же задачу,
			// но в другом направлении.
			for(int r=0; r<nGhostCells+2; r++) 
				switch(i - stencil[nDimCounter]) {
				case 2:
					Up[i][nDimCounter] += (rm1rENO2[r] + flag0*deltarm2rENO3[r])*U[i-2+r][nDimCounter];
					Um[i][nDimCounter] += (rm1lENO2[r] + flag0*deltarm2lENO3[r])*U[i-2+r][nDimCounter];
					break;
				case 1:
					Up[i][nDimCounter] += (rm1rENO2[r] + flag1*deltarm1rENO3[r])*U[i-2+r][nDimCounter];
					Um[i][nDimCounter] += (rm1lENO2[r] + flag1*deltarm1lENO3[r])*U[i-2+r][nDimCounter];
					break;
				case 0:
					Up[i][nDimCounter] += (r0rENO2[r] + flag2*deltar0rENO3[r])*U[i-2+r][nDimCounter];
					Um[i][nDimCounter] += (r0lENO2[r] + flag2*deltar0lENO3[0])*U[i-2+r][nDimCounter];
					break;
				default:
					cerr << "CSolver::calcHydroStageENO3G() reports error: unknown ENO stencil type!" << endl; 
					exit(1);
				} 
			/*if(stencil[nDimCounter] == i-2) {				
				Up[i][nDimCounter] = (rm1rENO2[0] + flag0*deltarm2rENO3[0])*U[i-2][nDimCounter] + 
					                 (rm1rENO2[1] + flag0*deltarm2rENO3[1])*U[i-1][nDimCounter] + 
									 (rm1rENO2[2] + flag0*deltarm2rENO3[2])*U[i][nDimCounter]   + 
									 (rm1rENO2[3] + flag0*deltarm2rENO3[3])*U[i+1][nDimCounter] + 
									 (rm1rENO2[4] + flag0*deltarm2rENO3[4])*U[i+2][nDimCounter];
				Um[i][nDimCounter] = (rm1lENO2[0] + flag0*deltarm2lENO3[0])*U[i-2][nDimCounter] + 
									 (rm1lENO2[1] + flag0*deltarm2lENO3[1])*U[i-1][nDimCounter] + 
									 (rm1lENO2[2] + flag0*deltarm2lENO3[2])*U[i][nDimCounter] + 
									 (rm1lENO2[3] + flag0*deltarm2lENO3[3])*U[i+1][nDimCounter] + 
									 (rm1lENO2[4] + flag0*deltarm2lENO3[4])*U[i+2][nDimCounter];				
			} else if (stencil[nDimCounter] == i-1) {
				Up[i][nDimCounter] = (rm1rENO2[0] + flag1*deltarm1rENO3[0])*U[i-2][nDimCounter] + 
					                 (rm1rENO2[1] + flag1*deltarm1rENO3[1])*U[i-1][nDimCounter] + 
									 (rm1rENO2[2] + flag1*deltarm1rENO3[2])*U[i][nDimCounter]   + 
									 (rm1rENO2[3] + flag1*deltarm1rENO3[3])*U[i+1][nDimCounter] + 
									 (rm1rENO2[4] + flag1*deltarm1rENO3[4])*U[i+2][nDimCounter];
				Um[i][nDimCounter] = (rm1lENO2[0] + flag1*deltarm1lENO3[0])*U[i-2][nDimCounter] + 
									 (rm1lENO2[1] + flag1*deltarm1lENO3[1])*U[i-1][nDimCounter] + 
									 (rm1lENO2[2] + flag1*deltarm1lENO3[2])*U[i][nDimCounter] + 
									 (rm1lENO2[3] + flag1*deltarm1lENO3[3])*U[i+1][nDimCounter] + 
									 (rm1lENO2[4] + flag1*deltarm1lENO3[4])*U[i+2][nDimCounter];
			} else if (stencil[nDimCounter] == i) {				
				Up[i][nDimCounter] = (r0rENO2[0] + flag2*deltar0rENO3[0])*U[i-2][nDimCounter] + 
					                 (r0rENO2[1] + flag2*deltar0rENO3[1])*U[i-1][nDimCounter] + 
									 (r0rENO2[2] + flag2*deltar0rENO3[2])*U[i][nDimCounter]   + 
									 (r0rENO2[3] + flag2*deltar0rENO3[3])*U[i+1][nDimCounter] + 
									 (r0rENO2[4] + flag2*deltar0rENO3[4])*U[i+2][nDimCounter];
				Um[i][nDimCounter] = (r0lENO2[0] + flag2*deltar0lENO3[0])*U[i-2][nDimCounter] + 
									 (r0lENO2[1] + flag2*deltar0lENO3[1])*U[i-1][nDimCounter] + 
									 (r0lENO2[2] + flag2*deltar0lENO3[2])*U[i][nDimCounter] + 
									 (r0lENO2[3] + flag2*deltar0lENO3[3])*U[i+1][nDimCounter] + 
									 (r0lENO2[4] + flag2*deltar0lENO3[4])*U[i+2][nDimCounter];
			}*/
		}
	}
	cout << endl << "Debug: ENO3 reconstruction: " << nimin2Counter << " i-2-stencils, " << nimin1Counter << " i-1-stencils, " << niCounter << " i-stencils" << endl;
	// Boundary nodes U+[0] and U-[nSize]
	Up[mini-1] = Up[mini]; Um[maxi] = Um[maxi-1];
	/*	for(i=1; i<maxi+nGhostCells; i++) {
		double xm = (double)(i-nGhostCells)*h;
		double xp = xm + h;
		double x = .5*(xm + xp);
		U[i][0] = 500.*1./h * (xp*xp*xp/3. - 3.*xp*xp/4. - xm*xm*xm/3. + 3.*xm*xm/4.);
		fcentres << x  << " " << U[i][0] << endl;
		fborders << xm << " " << Up[i-1][0] << " " << Um[i][0] << " " << 500.*(xm*xm - 1.5*xm) << endl;
	}
	fcentres.close();
	fborders.close();	*/
 	for(i=0; i<=nSize; i++) {	
		ms[i].F = calcGodunovFlux(Up[i-1+nGhostCells][0], Up[i-1+nGhostCells][1], Up[i-1+nGhostCells][2], 
			                      Um[i+nGhostCells][0], Um[i+nGhostCells][1], Um[i+nGhostCells][2]);	
	}
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
	EOSOld &eos = task.getEOS(); 
	const double gamma = eos.getGamma();
	double E=0.; 
	int i=0, j=0;
	const int nSize = ms.getSize();
	const double h=ms[1].x-ms[0].x;
	Vector4 L = Vector4::ZERO, R = Vector4::ZERO, D = Vector4::ZERO, V = Vector4::ZERO;	
	// Специальные массивы для хранения переменных и потоков, с дополнительными фиктивными ячейками
	const int nGhostCells = 2;	
	const int mini = nGhostCells, maxi = nGhostCells+nSize;
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
	
	
	
	
	
	// Former:
	/*U[2] = U[mini]; U[1] = U[mini]; U[0] = U[mini];  
	U[maxi] = U[maxi-1]; U[maxi+1] = U[maxi-1]; U[maxi+2] = U[maxi-1];*/
	
	// Now:
	U[mini-1] = U[mini]; U[mini-2] = U[mini]; //U[0] = U[mini];  
	U[maxi] = U[maxi-1]; U[maxi+1] = U[maxi-1]; //U[maxi+2] = U[maxi-1];
	
	
	
	
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




		if(i==29) {

			double qq = 0.;

		}




		
		ms[i].F = calcGodunovFlux(Up[i-1+nGhostCells][0], Up[i-1+nGhostCells][1], Up[i-1+nGhostCells][2], 
//			                      Um[i+nGhostCells][0], Um[i+nGhostCells][1], Um[i+nGhostCells][2]);
//		ms[i].F = calcRoeFlux(Up[i-1+nGhostCells][0], Up[i-1+nGhostCells][1], Up[i-1+nGhostCells][2], 
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
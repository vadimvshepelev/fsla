#ifndef __CTEST_H_
#define __CTEST_H_

class CTestToro {
public:
	double gamma; 
	double roL, vL, pL, eL, EL, 
		   roR, vR, pR, eR, ER;
	double q0;   // Начальное значение координаты разрыва 
	double tMax; // Время, до которого ведется счет
	CTestToro(double _gamma, 
		      double _roL, double _uL, double _vL, double _wL, double _pL, 
			  double _roR, double _uR, double _vR, double _wR, double _pR, double _q0, double _tMax) :
	          gamma(_gamma), roL(_roL), vL(_vL), pL(_pL), roR(_roR), vR(_vR), pR(_pR), q0(_q0), tMax(_tMax) {
				  eL = pL/(gamma-1.)/roL;
				  EL = eL + (vL*vL)/2.;
				  eR = pR/(gamma-1.)/roR;
				  ER = eR + (vR*vR)/2.;
			  }
};


#endif
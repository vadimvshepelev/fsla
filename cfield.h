
#ifndef MATTERSTATE_H
#define MATTERSTATE_H


#include "node.h"
#include <cstring>

using namespace std;

class CTask;

class CFieldOld {
public:

	CFieldOld();
	~CFieldOld();

	void	initData(CTask *task);
	void	clearData();
	double	loadData(string fName, int nCut);
	int getSize() { return nSize; }
	void	setSize(int nSizeNew) { nSize = nSizeNew; }

	void	setEdgeMode(bool edge_m) { bEdgeMode = edge_m; }
	bool	getEdgeMode() { return bEdgeMode; }

	inline Node &operator [] (const int i)
    {
        if(i < 0)		return bEdgeMode ? nodes[0]		  : left_edge;
		if(i > (int)nSize)   return bEdgeMode ? nodes[nSize-1] : right_edge;
        return nodes[i];
    }

	// Добавил функцию для прямого обращения к узлам сетки

	Node* getnodes() { return nodes; }

	// Функции, которая ставят "прозрачные" граничные условия

	void setEdgeTransparent(void);

private:

	void	setEdge(Node &n, double x, double dm);

	Node	*nodes;		// Указатель на массив узлов сетки
	int nSize;		// Размер массива

	Node	left_edge;	// Левые граничные условия
	Node	right_edge;	// Правые граничные условия
	bool	bEdgeMode;  // Режим обработки граничных условий
};


#endif

#ifndef MATTERSTATE_H
#define MATTERSTATE_H


#include "node.h"
#include <cstring>
using namespace std;

/*
   Объявление класса CTask вместо включения task.h.
   Чтобы избежать закольцованных друг на друга хедеров.
*/

class CTask;


/*
   Класс MatterState отвечает за создание и инициализацию сеточных 
   функций (массивов) всех физических величин. Т.е., на  классе лежит:
   - хранение всех вычислительных сеток;
   - выделение и освобождение памяти для них;
   - постановка начальных условий, извлеченных из входного файла классом CTask.

   "Фишка" класса в обработке граничных условий.
   Предусмотрено два режима их обработки:
  
   1. setEdgeMode(false)  (по умолчанию)
      При обращении к любому элементу за пределами сеточной функции
	  будет подставлен специальный Node, являющийся левым или правым
	  граничным условием.

	  Граничные Node'ы для этого режима можно задать просто присвоив
	  значения элементам ms[-1] и ms[ms.getSize()]. Либо изменив
	  значения по умолчанию в private-методе setEdge(...).

   2. setEdgeMode(true)
      При обращении к любому элементу за пределами сеточной функции
	  будет подставлен первый или последний Node сетки.
*/

class MatterState
{
public:

	MatterState();
	~MatterState();

	void	initData(CTask *task);
	void	clearData();
	double	loadData(string fName, unsigned int nCut);
	unsigned int getSize() { return nSize; }
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
	unsigned int nSize;		// Размер массива

	Node	left_edge;	// Левые граничные условия
	Node	right_edge;	// Правые граничные условия
	bool	bEdgeMode;  // Режим обработки граничных условий
};


#endif
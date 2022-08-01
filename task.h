#ifndef _CTask_H_
#define _CTask_H_

#include <stdio.h>
#include <stdlib.h>

#include "defines.h"
#include "methodold.h"
#include "eosold.h"
#include "CTestToro.h"
#include "eos/EOSBin.h"
#include "cfieldold.h"

//#include "eosTableFeAlpha.h"

/* Структура, обозначающая зону. Вычислительная область -- набор зон. */

struct Zone  {
	int n;				// количество ячеек в зоне
	double l;			// ширина зоны
	double ro;			// плотность
	double te;			// электронная температура
	double ti;			// ионная температура
	double v;			// скорость
	double expProperty;	// показатель экспоненты. Если он нулевой, экспоненциального 
						// уменьшения температур нет, зона однородна.
	double e;
};

/*
   Класс CTask считывает постановку задачи из входного файла,
   управляет именами входных и выходных файлов, содержит
   уравнение состояния, метод вычисления гидродинамики.


   Спецификация входного файла:
   
   1. Файл читается в текстовом режиме. Каждая строка в файле -- присваивание
      значения какой-либо переменной. Знаки равенства отделяются пробелами.
	  "nZones=5", "nZones =5", "nZones= 5" -- примеры неправильных записей.
	  "nZones = 5" -- правильная запись.

   2. Порядок следования переменных важен.

   3. Список переменных в заголовке:
	  EOSType = ideal | analytic | table - уравнение состояния
	  Method = euler | lagrange			 - метод вычисления гидродинамики
	  HydroStage = 0 | 1				 - расчитывать гидродинамику?
	  HeatStage = 0 | 1					 - расчитывать электронную теплопроводность?
	  ExchangeStage = 0 | 1				 - расчитывать электронно-ионный обмен?

	  nZones = 2						 - количество зон, определенных в файле

   4. Список переменных для каждой зоны:
	  l = 100.0e-9  - длина зоны
	  nSize = 50    - количество точек в зоне
	  ro = 0.1e-8   - плотность в зоне
	  ti = 300.0    - ионная температура в зоне
	  te = 0.0      - электронная температура в зоне
	  v = 0.0       - скорость в зоне
	  exp = 0.0     - экспоненцийальное затухание в зоне или 0

   5. Данные разных зон удобно отделять пустой строкой.


   Ограничения на комбинации параметров:

   1. Method = euler поддерживает только EOSType = ideal.

   2. HeatStage = 1, ExchangeStage = 1 и te != 0 для 
      EOSType = ideal не имеют смысла.
*/

enum TaskType {undef, ruGlass, LH1D, auWater, MieGruneisenProblem, RP1D};

enum SourceType {SrcUndef, SrcNone, SrcGlass, SrcMetal, SrcSq, Src5Layers};

enum MethodType {nomtd, samarskii, godunov1, eno2g, muscl, bgk, eno3g, hll};

class CTask {
public:
	CTask() : type(TaskType::undef), bHydroStage(false), bHeatStage(false), bExchangeStage(false), 
		      eos(0), eosGlass(0), eosBin(EOSBin(-1., 0., 0.)), sourceFlag(SourceType::SrcUndef), tauPulse(0.), fluence(0.), deltaSkin(0.), 
			  zones(0), nZones(0), maxTime(0.), mtd(0), viscFlag(0), CFL(0.), totalSize(0), 
			  EOSFlag(0), methodFlag(MethodType::nomtd), 
			  roL(0.), uL(0.), pL(0.), roR(0.), uR(0.), pR(0.), xMin(0.), xMax(0.), tMin(0.), tMax(0.), qBound(0.), NX(0) {}
	CTask(TaskType _type, bool _bHydroStage=true, bool _bHeatStage=false, bool _bExchangeStage=false, EOSOld* _eos=0, EOSOld* _eosGlass=0, 
		  SourceType _sourceFlag=SourceType::SrcNone,
		  double _tauPulse=0., double _fluence=0., double _deltaSkin=0., Zone* _zones=0, int _nZones=2, double _maxTime=0., CMethodOld* _mtd=0,
		  int _viscFlag=0, double _CFL=0., int _totalSize=0, int _EOSFlag=0, MethodType _methodFlag=MethodType::godunov1) :
		  type(_type), 
		  bHydroStage(_bHydroStage), bHeatStage(_bHeatStage), bExchangeStage(_bExchangeStage), 
		  eos(_eos), eosGlass(_eosGlass),
		  sourceFlag(_sourceFlag),
		  tauPulse(_tauPulse), fluence(_fluence), deltaSkin(_deltaSkin), 
		  zones(_zones), nZones(_nZones), 
		  maxTime(_maxTime),
		  mtd(_mtd), viscFlag(_viscFlag), CFL(_CFL), totalSize(_totalSize), 
		  EOSFlag(_EOSFlag), methodFlag(_methodFlag), roL(0.), uL(0.), pL(0.), roR(0.), uR(0.), pR(0.), xMin(0.), xMax(0.), tMin(0.), tMax(0.), qBound(0.), NX(0) {}
    CTask(TaskType _type, EOSBin _eos, double _roL, double _uL, double _pL, double _roR, double _uR, double _pR, double _xMin, double _xMax, double _tMin, double _tMax, 
		  double _qBound, int _NX, double _CFL, MethodType _methodFlag) : 
		  type(_type), eosBin(_eos), roL(_roL), uL(_uL), pL(_pL), roR(_roR), uR(_uR), pR(_pR), 
		  xMin(_xMin), xMax(_xMax), tMin(_tMin), tMax(_tMax), qBound(_qBound), NX(_NX), CFL(_CFL), methodFlag(_methodFlag), 
		  bHydroStage(false), bHeatStage(false), bExchangeStage(false), 
		  eos(0), eosGlass(0), sourceFlag(SourceType::SrcUndef), tauPulse(0.), fluence(0.), deltaSkin(0.), 
		  zones(0), nZones(0), maxTime(0.), mtd(0), viscFlag(0), totalSize(0), EOSFlag(0) {} 
	~CTask();					
	void		load(char *fName);		// Считать условие задачи из файла	
	void		clear();				// Удалить задачу
	EOSOld			&getEOS() { return *eos; }	
	EOSOld			&getEOSGlass() { return *eosGlass; }	
	void		setEOS(EOSOld* newEOS) { eos = newEOS; } ;			
	int			getEOSFlag() { return EOSFlag; }
	CMethodOld		&getMethod() { return *mtd; }	
	Zone		&getZone(int i) { return zones[i]; }
	int getNumZones() { return nZones; }
	int getTotalSize() { return totalSize; }
	bool		getHydroStage() { return bHydroStage; }
	bool		getHeatStage() { return bHeatStage; }
	bool		getExchangeStage() { return bExchangeStage; }
	string getOutputFilename() { return outputFileName; }
	string getInputFilename() { return inputFileName; }
	string getFlowFilename() { return flowFileName; }
	string getTaskName() { return taskName; }	
	SourceType	getSourceFlag() { return sourceFlag; }
	int			getViscFlag() { return viscFlag; }	
	MethodType  getMethodFlag() { return methodFlag; }		
	double		getTauPulse() { return tauPulse;}
	double		getFluence() { return fluence; }
	double		getDeltaSkin() { return deltaSkin; }
	double		getCFL() { return CFL; }
	double		getMaxTime() { return maxTime; }
	TaskType type;
	EOSBin eosBin;			// Двучленное УРС, если понадобится
	const double roL, uL, pL, roR, uR, pR, xMin, xMax, tMin, tMax, qBound; 
    const int NX;
	double CFL;   
	MethodType methodFlag;	
private:
	void		buildFileNames(string inputName);
	const char *readStringParam(FILE *f, const char *name);
	double		readFloatParam(FILE *f, const char *name);
	int			readIntParam(FILE *f, const char *name);
	bool	bHydroStage;		// Расчитывать гидродинамику?
	bool	bHeatStage;			// Расчитывать электронную теплопроводность?
	bool	bExchangeStage;		// Расчитывать электронно-ионный обмен?
	EOSOld		*eos;				// Уравнение состояния
	EOSOld		*eosGlass;			// Уравнение состояния стекла	
	SourceType sourceFlag;			// Источник энергии в расчетной области 
	double  tauPulse;			// Длительность вспышки
	double  fluence;			// Флюенс -- энергия, поглощенная единицей поверхности
	double  deltaSkin;			// Толщина скин-слоя
	Zone	*zones;				// Указатель на массив зон	
	int nZones;				// Количество зон
	double  maxTime;			// Время, до которого ведется расчет
	CMethodOld	*mtd;				// Метод расчета - Эйлер/Лагранж
	int		viscFlag;			// Искусственная вязкость ( 0 - не использовать, 1 - использовать).
	int		totalSize;			// Суммарное количество точек во всех зонах
	int		EOSFlag;			// Для табличного УРС -- тип таблиц (0 - нормальные, 1 - расширенные)	
	string taskName, inputFileName, outputFileName, flowFileName, EOSDirName; // EOSDirName для табличного УРС -- каталог в table_data, в котором лежат таблицы
};

#endif

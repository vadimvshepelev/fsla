#ifndef _CTask_H_
#define _CTask_H_

#include <stdio.h>
#include <stdlib.h>

#include "defines.h"
#include "method.h"
#include "eos.h"
#include "CTestToro.h"
#include "eos\\EOSBin.h"
//#include "eosTableFeAlpha.h"

/* Структура, обозначающая зону. Вычислительная область -- набор зон. */

struct Zone 
{
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

enum TaskType { undef, RuGlass, LH1D };

class CTask {
public:
	CTask();
	CTask(TaskType _type, bool _bHydroStage, bool _bHeatStage, bool _bExchangeStage, EOS* _eos, EOS* _eosGlass, int _sourceFlag,
		  double _tauPulse, double _fluence, double _deltaSkin, Zone* _zones, unsigned int _nZones, double _maxTime, CMethod* _mtd,
		  int _viscFlag, double _courant, unsigned int _totalSize, int _EOSFlag, int _methodFlag) :
		  type(_type), 
		  bHydroStage(_bHydroStage), bHeatStage(_bHeatStage), bExchangeStage(_bExchangeStage), 
		  eos(_eos), eosGlass(_eosGlass),
		  sourceFlag(_sourceFlag),
		  tauPulse(_tauPulse), fluence(_fluence), deltaSkin(_deltaSkin), 
		  zones(_zones), nZones(_nZones), 
		  maxTime(_maxTime),
		  mtd(_mtd), viscFlag(_viscFlag), courant(_courant), totalSize(_totalSize), 
		  EOSFlag(_EOSFlag), methodFlag(_methodFlag) {}
	~CTask();					
	void		load(char *fName);		// Считать условие задачи из файла
	void		clear();				// Удалить задачу
	EOS			&getEOS() { return *eos; }	
	EOS			&getEOSGlass() { return *eosGlass; }	
	void		setEOS(EOS* newEOS) { eos = newEOS; } ;			
	int			getEOSFlag() { return EOSFlag; }
	CMethod		&getMethod() { return *mtd; }	
	Zone		&getZone(int i) { return zones[i]; }
	unsigned int getNumZones() { return nZones; }
	unsigned int getTotalSize() { return totalSize; }
	bool		getHydroStage() { return bHydroStage; }
	bool		getHeatStage() { return bHeatStage; }
	bool		getExchangeStage() { return bExchangeStage; }
	char	   *getOutputFilename() { return outputFileName; }
	char	   *getInputFilename() { return inputFileName; }
	const char *getFlowFilename() { return flowFileName; }
	char	   *getTaskName() { return taskName; }	
	int			getSourceFlag() { return sourceFlag; }
	int			getViscFlag() { return viscFlag; }	
	int			getMethodFlag() { return methodFlag; }	
	double		getTauPulse() { return tauPulse;}
	double		getFluence() { return fluence; }
	double		getDeltaSkin() { return deltaSkin; }
	double		getCourant() { return courant; }
	double		getMaxTime() { return maxTime; }
	TaskType type;

	EOSBin  *eosBin;			// Указатель на двучленное УРС, если понадобится
private:
	void		buildFileNames(const char *inputName);
	const char *readStringParam(FILE *f, const char *name);
	double		readFloatParam(FILE *f, const char *name);
	int			readIntParam(FILE *f, const char *name);
	bool	bHydroStage;		// Расчитывать гидродинамику?
	bool	bHeatStage;			// Расчитывать электронную теплопроводность?
	bool	bExchangeStage;		// Расчитывать электронно-ионный обмен?
	EOS		*eos;				// Уравнение состояния
	EOS		*eosGlass;			// Уравнение состояния стекла	
	int		sourceFlag;			// Источник энергии в расчетной области (0 - нет, 1 - алюминий, 2 - алюминий на стекле).
	double  tauPulse;			// Длительность вспышки
	double  fluence;			// Флюенс -- энергия, поглощенная единицей поверхности
	double  deltaSkin;			// Толщина скин-слоя
	Zone	*zones;				// Указатель на массив зон	
	unsigned int nZones;				// Количество зон
	double  maxTime;			// Время, до которого ведется расчет
	CMethod	*mtd;				// Метод расчета - Эйлер/Лагранж
	int		viscFlag;			// Искусственная вязкость ( 0 - не использовать, 1 - использовать).
	double  courant;			// Число Куранта (начальное)

	int		totalSize;			// Суммарное количество точек во всех зонах
	int		EOSFlag;			// Для табличного УРС -- тип таблиц (0 - нормальные, 1 - расширенные)
	int		methodFlag;			// Флаг, обозначающий выбранный метод. 0 - Лагранж, консервативная схема Самарского,
								// 1 - Эйлер, схема Белоцерковского-Гущина-Коньшина, 2 - Двумерный Эйлер (в разработке).

	char inputFileName[_MAX_PATH];
	char outputFileName[_MAX_PATH];
	char flowFileName[_MAX_PATH];
	char EOSDirName[_MAX_PATH];  // Для табличного УРС -- каталог в table_data, в котором лежат таблицы
	char taskName[_MAX_PATH];
};






#endif
#include "defines.h"
#include "task.h"
#include "eosAnalytic.h"
#include "eosIdeal.h"
#include "eosTable.h"
#include "eosTest.h"

#include "eosAnalyticNi.h"
#include "eosTableNi.h"
#include "eosAnalyticTa.h"
#include "eosTableTa.h"
#include "eosAnalyticFe.h"
#include "eosTableFe.h"
#include "eosTableFeAlpha.h"
#include "eosTableAu.h"

#include "methodEuler.h"
#include "methodEuler2D.h"

#include<stdlib.h>
#include<string.h>

#include<iostream>
using namespace std;

CTask::CTask()
{
	eos = 0;
	mtd = 0;
	bHydroStage = false;
	bHeatStage = false;
	bExchangeStage = false;

	    zones = 0;
	   nZones = 0;
	totalSize = 0;

	sourceFlag = 0;
	  viscFlag = 0;
	methodFlag = 0;

	  tauPulse = 0;	
   
	fluence = 0.0;
	courant = 0.0;
	maxTime = 0.0;

	type = TaskType::undef; // Пока не навел порядок в инфраструктуре, проставляю этот флаг руками в CSolver::goGlass()

}

//CTask::CTask(EOS *_eos, CMethod _mtd, ) 
//{}



CTask::~CTask()
{
	clear();
}


void CTask::load(char* fName) {
	clear();
	unsigned int j = 0;
	EOSFlag = 0;
	char task_id[_MAX_PATH];
	strcpy(task_id, fName);
	task_id[strlen(task_id)-4]='\0';
	strcpy(taskName, task_id); 
	buildFileNames(fName);
	FILE* f=fopen(inputFileName, "r");
	if(!f)	{
		printf("Can't open task file: %s\n", inputFileName);
		exit(1);
	}
	printf("Processing task file: %s\n", inputFileName);
	// parse EOSType
	const char *eos_type = readStringParam(f, "EOSType");
	
	if     (!strcmp(eos_type, "table"))
	{
		const char* tableDir  = readStringParam(f, "TableDir");
		
		// Вопрос -- почему, если эту строчку поставить сразу после следующей:
		// "const char* tableFlag = readStringParam(f, "TableFlag");", то
		// tableDir становится равным tableFlag и равным "normal" ???
		//
		// 13.01.2012: Сдается мне, он почему-то пишет в одну область памяти эти значения. Пичалька.
		// Тут, короче, какой-то баг, непонятка. Можно у Руди спросить.

		strcpy(EOSDirName, tableDir);
		
		//////////////////////////////////////////////////////////////////////

		const char* tableFlag = readStringParam(f, "TableFlag");

		if(!strcmp(tableFlag, "normal")) EOSFlag = 0;
		else if(!strcmp(tableFlag, "extended")) EOSFlag = 1;
		else
		{
			printf("Unknown TableFlag parameter: %s\n", tableFlag);
			exit(1);	
		}
	
		//////////////// Очень стремное условие -- проверка по каталогу. Переписать.
		if(!strcmp(EOSDirName, "ni")) {
			eos	= new EOSTableNi(EOSDirName, EOSFlag, 8900.);
		} 
		else if(!strcmp(EOSDirName, "au")) {
			eos	= new EOSTableAu(EOSDirName, EOSFlag, 19301.);
		}



/*		else if(!strcmp(EOSDir, "ta"))
		{
			eos	= new EOSTableTa(EOSDir, EOSFlag, 16654.);
		}*/
		else if(!strcmp(EOSDirName, "fe"))
		{
			eos	= new EOSTableFe(EOSDirName, EOSFlag, 7874.);
		}
		else if(!strcmp(EOSDirName, "fe_comb"))
		{
			eos	= new EOSTableFe(EOSDirName, EOSFlag, 7874.);
		}
		else if(!strcmp(EOSDirName, "fe_alpha"))
		{
			eos	= new EOSTableFeAlpha(EOSDirName, EOSFlag, 7874.);
		}
		else
		{	
			eos	= new EOSTable(EOSDirName, EOSFlag, 2700.);
		}
		////////////////////////////////////////////////////////////////////////////
	} else if(!strcmp(eos_type, "analytic"))	eos = new EOSAnalytic();
	else if(!strcmp(eos_type, "ideal"))		eos = new EOSIdeal(1.4);
	else if(!strcmp(eos_type, "test"))		eos = new EOSTest();
	else if(!strcmp(eos_type, "MieGruneisenRu")) {
		eos = new EOSMieGruneisenRu();
		eosGlass = new EOSTable("new", 1, 2700.);
	} else {
		printf("Unknown EOSType: %s\n", eos_type);
		exit(1);
	}
	printf("EOS: %s\n", eos_type);
	// parse Method
	const char *method = readStringParam(f, "Method");
	if(!strcmp(method, "lagrange")) {
		methodFlag = 0;
		printf("Method: Samarskii scheme, lagrangian variables\n");
	} else if(!strcmp(method, "euler"))	{
		mtd = new CMethodEuler(eos);
		methodFlag = 1;
		printf("Method: Belotserkovskii-Gushchin-Konshin scheme, eulerian variables\n");
	} else {
		printf("Unknown Method: %s\n", method);
		exit(1);
	}

	// parse Stages

	bHydroStage		 = readIntParam(f, "HydroStage")		? true : false;
	bHeatStage		 = readIntParam(f, "HeatStage")		    ? true : false;
	bExchangeStage	 = readIntParam(f, "ExchangeStage")	    ? true : false;
	// parse Source and Viscosity
	const char *sourceParam = readStringParam(f, "source");
	if(!strcmp(sourceParam, "no")) 	{	
		sourceFlag = 0;
		printf("Source type: none\n");
	} else if(!strcmp(sourceParam, "default")) {	
		sourceFlag = 1;
		printf("Source type: vacuum-metal-vacuum\n");
	} else if(!strcmp(sourceParam, "Al_glass")) {	
		sourceFlag = 2;
		printf("Source type: glass-metal-vacuum\n");
		eosGlass = new EOSPyrexGlass();
	} else if(!strcmp(sourceParam, "Al_glass_sq")) {	
		sourceFlag = 3;
		printf("Source type: Aluminium on the glass (square pulse)\n");
	} else if(!strcmp(sourceParam, "5_layers_Si")) {
		sourceFlag = 4;
		cout << "Source type: 5-layer Si source" << endl;
	} else {
		printf("Unknown Source type: %s\n", sourceParam);
		exit(1);
	}

	tauPulse = readFloatParam(f, "tauPulse");
	if (tauPulse>0)
		printf("Laser pulse duration: %e seconds\n", tauPulse);
	else 
	{
		printf("Incorrect \"tauPulse\" parameter (must be positive): %e\n", tauPulse);
		exit(1);
	}

	fluence = readFloatParam(f, "fluence");
	if (fluence>=0)
		printf("Fluence: %f J/m^2\n", fluence);
	else 
	{
		printf("Incorrect \"fluence\" parameter (must be positive): %f\n", fluence);
		exit(1);
	}

	deltaSkin = readFloatParam(f, "deltaSkin");
	if (deltaSkin>0)
		printf("Diameter of skin-layer: %e m\n", deltaSkin);
	else 
	{
		printf("Incorrect \"deltaSkin\" parameter (must be positive): %f\n", deltaSkin);
		exit(1);
	}

	courant = readFloatParam(f, "courant");
	if (courant>0)
		printf("Courant number: %e\n", courant);
	else 
	{
		printf("Incorrect courant parameter (must be positive): %e\n", courant);
		exit(1);
	}

	viscFlag = readIntParam(f, "viscosity");
	if (viscFlag==0)
		printf("Artifficial viscosity: off\n");
	else if (viscFlag==1)
		printf("Artifficial viscosity: on\n");
	else 
	{
		printf("Unknown Viscosity parameter (must be 0 or 1): %d\n", viscFlag);
		exit(1);
	}

	maxTime = readFloatParam(f, "maxTime");
	if (maxTime>0)
		printf("Calculating up to t: %e seconds\n", maxTime);
	else 
	{
		printf("Incorrect maxTime parameter (must be positive): %e\n", maxTime);
		exit(1);
	}


	// parse nZones

	nZones = readIntParam(f, "nZones");
	
	if(nZones<=0)
	{
		printf("Invalid number of zones (nZones): %d\n", nZones);
		exit(1);
	}

	zones = new Zone[nZones];

	if(!zones)
	{
		printf("Can't allocate zones. nZones = %d\n", nZones);
		exit(1);
	}
	// Parse zones data
	for(j=0; j<nZones; j++) {
		zones[j].l  = readFloatParam(f,	"l");
		zones[j].n  = readIntParam  (f,	"nSize");
		zones[j].ro = readFloatParam(f, "ro");
		if(methodFlag==1) {
			zones[j].v  = readFloatParam(f,	"v");
			double _p = readFloatParam(f, "p");
			double _e = _p/(eos->getGamma()-1.)/zones[j].ro;
			zones[j].ti = eos->getti(zones[j].ro, _e);
			zones[j].te = zones[j].ti;
		} else {
			zones[j].ti = readFloatParam(f, "ti");
			zones[j].te = readFloatParam(f, "te");
			zones[j].v  = readFloatParam(f,	"v");
		}
		zones[j].expProperty = readFloatParam(f, "exp");
		printf("Zone %d loaded!\n", j);
		totalSize += zones[j].n;
	}
	printf("Task loaded! Total nSize=%d\n", totalSize);
	fclose(f);
}

void CTask::clear()
{
	if(eos) delete eos;
	eos = 0;
	if(eosGlass) delete eosGlass;
	eosGlass = 0;
	if(mtd) delete mtd;
	mtd = 0;
	
	bHydroStage = false;
	bHeatStage = false;
	bExchangeStage = false;

	if(zones) delete[] zones;
	zones = 0;
	nZones = 0;
	totalSize = 0;
}


const char *CTask::readStringParam(FILE *f, const char *name)
{
	static char param_name[50];
	static char param_value[50];

	//DEBUG!!!/////////////////////////////
//	cout << "Debugging readStringParam()" << endl;
//	cout << "Before reading: param_name = " << param_name << " param_value = "  << param_value << endl;
	///////////////////////////////////////

	if(fscanf(f, "%s = %s", param_name, param_value) != 2 )
	{
		printf("Syntax error: \'%s\' variable initialization string!\n", param_name);
		exit(1);
	}

	//DEBUG!!!/////////////////////////////
//	cout << "After reading:  param_name = " << param_name << " param_value = "  << param_value << endl;
	///////////////////////////////////////

	if(strcmp(param_name, name))
	{
		printf("Unexpected variable: %s (instead of %s)\n", param_name, name);
		exit(1);
	}

	return param_value;
}


double CTask::readFloatParam(FILE *f, const char *name)
{
	const char *val = readStringParam(f, name);
	return atof(val);
}


int CTask::readIntParam(FILE *f, const char *name)
{
	const char *val = readStringParam(f, name);
	return atoi(val);
}


void CTask::buildFileNames(const char *inputName) {
	strcpy(inputFileName, INPUT_FOLDER);
	strcat(inputFileName, inputName);

	strcpy(outputFileName, OUTPUT_FOLDER);
	strcat(outputFileName, TIME_SUFFIX);
	strcat(outputFileName, inputName);
	strcpy(outputFileName + strlen(outputFileName)-3, "dat");

	// Следующие две строки портят значение переменной task.taskName, удаляю их
//	strcpy(outputFileName, OUTPUT_FOLDER);
//	strcat(flowFileName, "flow.dat");
}


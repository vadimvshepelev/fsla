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
#include "eos\\EOSBin.h"

#include "methodEuler.h"

#include "ceos.h"

#include<stdlib.h>
#include<string.h>
#include<iostream>

using namespace std;

CTask::~CTask() {
	clear();
}

 void CTask::load(char* fName) {
	clear();
	unsigned int j = 0;
	EOSFlag = 0;
	char task_id[_MAX_PATH];
	strcpy(task_id, fName);
	task_id[strlen(task_id)-4]='\0';
	//strcpy(taskName, task_id); 
	taskName += task_id; 
	buildFileNames(fName);
	FILE* f=fopen(inputFileName.c_str(), "r");
	if(!f)	{
		printf("Can't open task file: %s\n", inputFileName.c_str());
		exit(1);
	}
	printf("Processing task file: %s\n", inputFileName.c_str());
	// parse EOSType
	const char *eos_type = readStringParam(f, "EOSType");
	
	if (!strcmp(eos_type, "table")) {
		const char* tableDir  = readStringParam(f, "TableDir");
		//strcpy(EOSDirName.c_str(), tableDir);
		EOSDirName = string(tableDir);
		const char* tableFlag = readStringParam(f, "TableFlag");
		if(!strcmp(tableFlag, "normal")) EOSFlag = 0;
		else if(!strcmp(tableFlag, "extended")) EOSFlag = 1;
		else {
			printf("Unknown TableFlag parameter: %s\n", tableFlag);
			exit(1);	
		}	
		if(!strcmp(EOSDirName.c_str(), "ni")) {
			eos	= new EOSTableNi(EOSDirName, EOSFlag, 8900.);
		} 
		else if(!strcmp(EOSDirName.c_str(), "au")) {
			eos	= new EOSTableAu(EOSDirName, EOSFlag, 19301.);
		}
/*		else if(!strcmp(EOSDir, "ta"))
		{
			eos	= new EOSTableTa(EOSDir, EOSFlag, 16654.);
		}*/
		else if(!strcmp(EOSDirName.c_str(), "fe")) {
			eos	= new EOSTableFe(EOSDirName, EOSFlag, 7874.);
		}
		else if(!strcmp(EOSDirName.c_str(), "fe_comb")) {
			eos	= new EOSTableFe(EOSDirName, EOSFlag, 7874.);
		}
		else if(!strcmp(EOSDirName.c_str(), "fe_alpha")) {
			eos	= new EOSTableFeAlpha(EOSDirName, EOSFlag, 7874.);
		} else {	
			eos	= new EOSTable(EOSDirName, EOSFlag, 2700.);
		}		
	} else if(!strcmp(eos_type, "analytic"))	
		eos = new EOSAnalytic();
	else if(!strcmp(eos_type, "ideal"))		
		eos = new EOSIdeal(1.4);
	else if(!strcmp(eos_type, "test"))		
		eos = new EOSTest();
	else if(!strcmp(eos_type, "bin")) {
		const double gamma = readFloatParam(f, "gamma");
		const double ro0 = readFloatParam(f, "ro0");
		const double c0 = readFloatParam(f, "c0");		
		eos = 0;
		eosBin = new EOSBin(gamma, ro0, c0);
	} else if(!strcmp(eos_type, "MieGruneisenRu")) {
		eos = new EOSMieGruneisenRu();
		eosGlass = new EOSTable("new", 1, 2700.);
	} else if(!strcmp(eos_type, "MieGruneisen")) {
		type = TaskType::MieGruneisenProblem;
		eos = NULL;
		eosGlass = NULL;
	} else {
		printf("Unknown EOSType: %s\n", eos_type);
		exit(1);
	}
	printf("EOS: %s\n", eos_type);
	// parse Method
	const char *method = readStringParam(f, "Method");
	if(!strcmp(method, "samarskii")) {
		methodFlag = MethodType::samarskii;
		printf("Method: Samarskii scheme, lagrangian variables\n");
	} else if(!strcmp(method, "bgk"))	{
		mtd = new CMethodEuler(eos);
		methodFlag = MethodType::bgk;
		printf("Method: Belotserkovskii-Gushchin-Konshin scheme, eulerian variables\n");
	} else if(!strcmp(method, "godunov1"))	{
		mtd = new CMethodEuler(eos);
		methodFlag = MethodType::godunov1;
		printf("Method: Godunov 1st-order scheme, eulerian variables\n");
	} else if(!strcmp(method, "eno3g"))	{
		mtd = new CMethodEuler(eos);
		methodFlag = MethodType::eno3g;
		printf("Method: ENO-3 reconstruction, Godunov 1st-order scheme Riemann solver, eulerian variables\n");
	} else {
		printf("Unknown Method: %s\n", method);
		exit(1);
	}

	// parse Stages
		if(methodFlag == MethodType::samarskii)	{ 
		bHydroStage		 = readIntParam(f, "HydroStage")		? true : false;
		bHeatStage		 = readIntParam(f, "HeatStage")		    ? true : false;
		bExchangeStage	 = readIntParam(f, "ExchangeStage")	    ? true : false;
	} else {
		bHydroStage		 = true;
		bHeatStage		 = false;
		bExchangeStage   = false;
	}

	// parse Source 
	if (methodFlag == MethodType::samarskii) {
		const char *sourceParam = readStringParam(f, "source");
		if(!strcmp(sourceParam, "no")) 	{	
			sourceFlag = SrcNone;
			printf("Source type: none\n");
		} else if(!strcmp(sourceParam, "default")) {	
			sourceFlag = SourceType::SrcMetal;
			printf("Source type: vacuum-metal-vacuum\n");
		} else if(!strcmp(sourceParam, "glass-metal")) {	
			sourceFlag = SourceType::SrcGlass;
			printf("Source type: glass-metal-vacuum\n");
			//eosGlass = new EOSPyrexGlass();
			eosGlass = new EOSSimpleWater();
		} else if(!strcmp(sourceParam, "Al_glass_sq")) {	
			sourceFlag = SourceType::SrcSq;
			printf("Source type: Aluminium on the glass (square pulse)\n");	
		} else if(!strcmp(sourceParam, "5_layers_Si")) {
			sourceFlag = SourceType::Src5Layers;
			cout << "Source type: 5-layer Si source" << endl;
		} else {
			printf("Unknown Source type: %s\n", sourceParam);
			exit(1);
		}
	} else {
		sourceFlag = SourceType::SrcNone;
	}

	if(methodFlag == MethodType::samarskii)	{
		tauPulse = readFloatParam(f, "tauPulse");
		if (tauPulse>0)
			printf("Laser pulse duration: %e seconds\n", tauPulse);
		else {
			printf("Incorrect \"tauPulse\" parameter (must be positive): %e\n", tauPulse);
			exit(1);
		}
		fluence = readFloatParam(f, "fluence");
		if (fluence>=0)
			printf("Fluence: %f J/m^2\n", fluence);
		else {
			printf("Incorrect \"fluence\" parameter (must be positive): %f\n", fluence);
			exit(1);
		}
		deltaSkin = readFloatParam(f, "deltaSkin");
		if (deltaSkin>0)
			printf("Diameter of skin-layer: %e m\n", deltaSkin);
		else {
			printf("Incorrect \"deltaSkin\" parameter (must be positive): %f\n", deltaSkin);
			exit(1);
		}
	} else {
		tauPulse = 0.; fluence = 0.; deltaSkin = 0.;
	}
	CFL = readFloatParam(f, "CFL");
	if (CFL > 0)
		printf("Courant number: %e\n", CFL);
	else {
		printf("Incorrect courant parameter (must be positive): %e\n", CFL);
		exit(1);
	}
	if(methodFlag == MethodType::samarskii)	{
		viscFlag = readIntParam(f, "viscosity");
		if (viscFlag==0)
			printf("Artifficial viscosity: off\n");
		else if (viscFlag==1)
			printf("Artifficial viscosity: on\n");
		else {
			printf("Unknown Viscosity parameter (must be 0 or 1): %d\n", viscFlag);
			exit(1);
		}
	} else {
		viscFlag = 0;
	}
	maxTime = readFloatParam(f, "maxTime");
	if (maxTime>0)
		printf("Calculating up to t: %e seconds\n", maxTime);
	else {
		printf("Incorrect maxTime parameter (must be positive): %e\n", maxTime);
		exit(1);
	}
	
	// parse nZones
	nZones = readIntParam(f, "nZones");	
	if(nZones<=0) {
		printf("Invalid number of zones (nZones): %d\n", nZones);
		exit(1);
	}
	zones = new Zone[nZones];
	if(!zones) {
		printf("Can't allocate zones. nZones = %d\n", nZones);
		exit(1);
	}
	// Parse zones data
	for(j=0; j<nZones; j++) {
		zones[j].l  = readFloatParam(f,	"l");
		zones[j].n  = readIntParam  (f,	"nSize");
		zones[j].ro = readFloatParam(f, "ro");
		if(methodFlag!=MethodType::samarskii) {
			zones[j].v  = readFloatParam(f,	"v");
			double _p = readFloatParam(f, "p");
			if(eos) {   
				double _e = _p/(eos->getGamma()-1.)/zones[j].ro;
				zones[j].ti = eos->getti(zones[j].ro, _e);
				zones[j].te = zones[j].ti;
				zones[j].e = 0.;
			} else if (eosBin) {
				double _e = eosBin->gete(zones[j].ro, _p);
				zones[j].ti = 0.;
				zones[j].te = 0.;
				zones[j].e = _e;
			}   else if(type == TaskType::MieGruneisenProblem) {
				// —овершенно адска€ зќплатка, никогда так больше не делай!!!
				CEOSMieGruneisen eos;
				zones[j].ti = 0.;
				zones[j].te = 0.;
				zones[j].e = eos.gete(zones[j].ro, _p);
			}
		} else {
			zones[j].ti = readFloatParam(f, "ti");
			zones[j].te = readFloatParam(f, "te");
			zones[j].v  = readFloatParam(f,	"v");
		}		
		printf("Zone %d loaded!\n", j);
		totalSize += zones[j].n;
	}
	printf("Task loaded! Total nSize=%d\n", totalSize);

	// And now let's set CTask::taskType variable
	if(!strcmp(EOSDirName.c_str(), "au") && (nZones == 2)) 
		type = TaskType::auWater;
	else if(type == TaskType::MieGruneisenProblem)
		// ≈ще одна заплатка, порожденна€ неверо€тно кривой архитектурой, от которой надо избавл€тьс€
		type == TaskType::MieGruneisenProblem;
	else
		type = TaskType::undef;

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


void CTask::buildFileNames(string inputName) {
	// strcpy(inputFileName.c_str(), INPUT_FOLDER);
	inputName = string (INPUT_FOLDER) + inputName;
	// strcat(inputFileName.c_str(), inputName);
	inputFileName += inputName;
	outputFileName = string(OUTPUT_FOLDER);
	outputFileName += string(TIME_SUFFIX);
	outputFileName += inputName;
	outputFileName += string(".dat");
	// —ледующие две строки порт€т значение переменной task.taskName, удал€ю их
//	strcpy(outputFileName, OUTPUT_FOLDER);
//	strcat(flowFileName, "flow.dat");
}


#ifndef _DEFINES_H
#define _DEFINES_H

//#define INPUT_FOLDER  "input_files/"
//#define OUTPUT_FOLDER "output_files/"
#define TIME_SUFFIX   "______" 
#define TABLE_FOLDER  "table_data/"

// Disabling error messages like "Function strcpy is unsafe. Use "
#define _CRT_SECURE_NO_WARNINGS

#ifdef _WIN32
//#define PATH_SEPARATOR   "\\" 
#include<stdlib.h>
#else
//#define PATH_SEPARATOR   "/" 
#include <limits.h>
#define _MAX_PATH PATH_MAX
#endif

extern char* INPUT_FOLDER;
extern char* OUTPUT_FOLDER;

#endif
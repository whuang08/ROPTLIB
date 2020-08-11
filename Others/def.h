/*
This is the global head file. Every file in ROPTLIB will include this file.

---- WH
*/

#ifndef DEF_H
#define DEF_H

//#define MATLAB_MEX_FILE//For debug---
//#define DRIVERJULIAPROB//For debug---

#define DOUBLE_PRECISION  /*SINGLE_PRECISION DOUBLE_PRECISION*/
//#define ROPTLIB_WITH_FFTW//When FFTW library is needed

#undef real
#include <cmath>

/*std library*/
#include <cstdio>
#include <cstdlib>

/*For obtaining the lower bound, upper bound of numbers of realdp precision*/
#include <climits>
#include <limits>

#include <cassert>

#ifdef SINGLE_PRECISION
    #define realdp float
    #define realdpcomplex complex
#else
    #define realdp double
    #define realdpcomplex doublecomplex
#endif

#ifndef MATLAB_MEX_FILE
#include <f2c.h>
#endif /* end of ifndef MATLAB_MEX_FILE */

#ifdef _WIN64 /* The following code is compiled only when this library is compiled in Windows (64-bit only)
	If the code is compile under DEBUG mode, then test wheter there is memory leakage or not*/
	#ifdef _DEBUG
	#define DEBUG_CLIENTBLOCK   new( _CLIENT_BLOCK, __FILE__, __LINE__)
	/*Use my code to help checking memory leakage. One has to define a global variable:
		std::map<integer *, integer> *CheckMemoryDeleted;
	before running the code.
	*/
	#else
	#define DEBUG_CLIENTBLOCK
	#endif

	/*This is used for checking the memory leakage in windows system*/
	#define _CRTDBG_MAP_ALLOC
	#include <crtdbg.h>

	/*This is used for checking the memory leakage in windows system if the code is run in DEBUG mode*/
	#ifdef _DEBUG
	#define new DEBUG_CLIENTBLOCK
	#endif

#elif _WIN32 /* The following code is compiled only when this library is compiled in Windows (both 32-bit and 64-bit only)
   define something for Windows (32-bit and 64-bit, this part is common) */
#elif __APPLE__ /* The following code is compiled only when this library is compiled in MAC */
    #include "TargetConditionals.h"
    #if TARGET_IPHONE_SIMULATOR
         /* iOS Simulator */
    #elif TARGET_OS_IPHONE
        /* iOS device */
    #elif TARGET_OS_MAC
        /* Other kinds of Mac OS */
    #else
        /* Unsupported platform */
    #endif
#elif __linux /* The following code is compiled only when this library is compiled in Linux system */
    /* linux */
	#ifdef __GNUC__
/*
		const class {
		public:
			template<class T>
			operator T*(void) const
			{
				return 0;
			}
			template<class C, class T>
			operator T C::*(void) const
			{
				return 0;
			}
		private:
			void operator&(void) const;
		} nullptr = {};
*/
	#endif /* end of __GNUC__ */
#elif __unix /* all unices not caught above */
    /* Unix */
#elif __posix
    /* POSIX */
#endif /* end of checking platforms */

#include "Others/BlasLapackCppWrapper.h"

/*Help to debug the code*/
#include "Others/ForDebug.h"

/*For obtain the computational time*/
#include "Others/Timer.h"

#include <map>
#include <string.h>

typedef std::map<std::string, realdp> PARAMSMAP;

/*Define the number PI*/
#define PI 3.14159265358979323846264

#endif /* end of DEF_H */

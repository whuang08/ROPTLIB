/*
This file defines some functions to help debugging.
In the current version, it only contains on function which prints a realdp array.

-----WH
*/

#ifndef FORDEBUG_H
#define FORDEBUG_H

#include "Others/def.h"
#include <iostream>

#ifdef MATLAB_MEX_FILE
#include "blas.h"
#include "lapack.h"
#define integer ptrdiff_t
#endif


/*Define the namespace*/
namespace ROPTLIB{

    namespace GLOBAL{
        extern integer IZERO, IONE, ITWO;
        extern realdp DZERO, DONE, DTWO, DNONE, DNTWO;
        extern realdpcomplex ZZERO, ZONE, ZTWO, ZNONE;
        extern char *N, *T, *L, *R, *V, *C, *U, *A, *S, *O;
    };

	class ForDebug{
	public:
		static void Print(const char *name, const realdp *M, integer row, integer col = 1, integer num = 1);
		static realdp NormF(const realdp *V, integer length);
	};
}; /*end of ROPTLIB namespace*/
#endif /* end of FORDEBUG_H */

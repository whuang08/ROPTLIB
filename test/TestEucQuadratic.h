/*
This is the test file to run the problem defined in EucQuadratic.h and EucQuadratic.cpp.

---- WH
*/

#ifndef TESTEUCQUADRATIC_H
#define TESTEUCQUADRATIC_H


#include <iostream>
#include "Others/randgen.h"
#include "Manifolds/Manifold.h"
#include "Problems/Problem.h"
#include "Solvers/SolversSMLS.h"
#include <ctime>

/*If this test file is called from Matlab, then functions in DriverMexProb.h are used.*/
#include "test/DriverMexProb.h"

//#include "Manifolds/Euclidean/EucVariable.h"
#include "Problems/EucQuadratic.h"

//#include "Problems/StieBrockett/StieBrockett.h"
//#include "Manifolds/Stiefel/StieVector.h"
//#include "Manifolds/Stiefel/StieVariable.h"
//#include "Manifolds/Stiefel/Stiefel.h"

#include "Solvers/RSD.h"
#include "Solvers/RNewton.h"
#include "Solvers/RCG.h"
#include "Solvers/RBroydenFamily.h"
#include "Solvers/RWRBFGS.h"
#include "Solvers/RBFGS.h"
#include "Solvers/LRBFGS.h"

#include "Solvers/SolversSMTR.h"
#include "Solvers/RTRSD.h"
#include "Solvers/RTRNewton.h"
#include "Solvers/RTRSR1.h"
#include "Solvers/LRTRSR1.h"

#include "Others/def.h"

using namespace ROPTLIB;

void testEucQuadratic(void);

#endif // end of TESTEUCQUADRATIC_H

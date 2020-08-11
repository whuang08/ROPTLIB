/*
This is the test file for the Brocokett problem defined in StieBrockett.h and StieBrockett.cpp.

---- WH
*/

#ifndef TESTSTIEBROCKETT_H
#define TESTSTIEBROCKETT_H

/*Output to console*/
#include <iostream>

/*Generate random number*/
#include "Others/randgen.h"

/*Computational time*/
#include <ctime>

/*If this test file is called from Matlab, then functions in DriverMexProb.h are used.*/
#include "test/DriverMexProb.h"

/*Problem related classes*/
#include "Problems/Problem.h"
#include "Problems/StieBrockett.h"
#include "Problems/SphereTxRQ.h"

/*Manifold related classes*/
#include "Manifolds/Manifold.h"
//#include "Manifolds/Stiefel/StieVariable.h"
#include "Manifolds/Stiefel.h"
#include "Manifolds/SphereTx.h"

/*Linesearch based solvers*/
#include "Solvers/RSD.h"
#include "Solvers/RNewton.h"
#include "Solvers/RCG.h"
#include "Solvers/RBroydenFamily.h"
#include "Solvers/RWRBFGS.h"
#include "Solvers/RBFGS.h"
#include "Solvers/LRBFGS.h"
#include "Solvers/RGS.h"
#include "Solvers/LRBFGSSub.h"
//#include "Solvers/RBFGSLPSub.h"

/*Trust-region based solvers*/
#include "Solvers/SolversSMTR.h"
#include "Solvers/RTRSD.h"
#include "Solvers/RTRNewton.h"
#include "Solvers/RTRSR1.h"
#include "Solvers/LRTRSR1.h"

/*The global head file*/
#include "Others/def.h"

#include "test/DriverMexProb.h"

using namespace ROPTLIB;

double LinesearchInput(integer iter, const Variable &x1, const Vector &exeta1, realdp initialstepsize, realdp initialslope, const Problem *prob, const Solvers *solver);
bool MyStop(const Variable &x, const Vector &funSeries, integer lengthSeries, realdp finalval, realdp initval, const Problem *prob, const Solvers *solver);
/*The main test function*/
void testStieBrockett(void);

#endif // end of TESTSTIEBROCKETT_H

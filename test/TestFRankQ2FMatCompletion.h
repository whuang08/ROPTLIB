/*
This is the test file to run the problem defined in LRMatrixCompletion.h and LRMatrixCompletion.cpp.

---- WH
*/

#ifndef TESTFRANKQ2MATCOMPLETION_H
#define TESTFRANKQ2MATCOMPLETION_H


#include <iostream>
#include "Others/randgen.h"
#include "Manifolds/Manifold.h"
#include "Problems/Problem.h"
#include <ctime>

#include "test/DriverMexProb.h"

#include "Problems/FRankQ2FMatCompletion.h"
#include "Manifolds/FixedRankQ2F.h"

#include "Solvers/RSD.h"
#include "Solvers/RNewton.h"
#include "Solvers/RCG.h"
#include "Solvers/RBroydenFamily.h"
#include "Solvers/RWRBFGS.h"
#include "Solvers/RBFGS.h"
#include "Solvers/LRBFGS.h"

#include "Solvers/RTRSD.h"
#include "Solvers/RTRNewton.h"
#include "Solvers/RTRSR1.h"
#include "Solvers/LRTRSR1.h"

#include "Others/def.h"

using namespace ROPTLIB;

void testFRankQ2FMatCompletion(void);

#endif

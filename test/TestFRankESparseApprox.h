/*
This is the test file to run the problem defined in WeightedLowRank.h and WeightedLowRank.cpp.

---- WH
*/

#ifndef TESTFRANKESPARSEAPPROX_H
#define TESTFRANKESPARSEAPPROX_H


#include <iostream>
#include "Others/randgen.h"
#include "Manifolds/Manifold.h"
#include "Problems/Problem.h"
#include <ctime>

#include "Problems/FRankESparseApprox.h"
//#include "Manifolds/FixedRankE/FixedRankEVariable.h"
#include "Manifolds/FixedRankE.h"

#include "Solvers/RSD.h"
#include "Solvers/RNewton.h"
#include "Solvers/RCG.h"
#include "Solvers/RBroydenFamily.h"
#include "Solvers/RWRBFGS.h"
#include "Solvers/RBFGS.h"
#include "Solvers/LRBFGS.h"
#include "Solvers/ManPG.h"
#include "Solvers/AManPG.h"

#include "Solvers/RTRSD.h"
#include "Solvers/RTRNewton.h"
#include "Solvers/RTRSR1.h"
#include "Solvers/LRTRSR1.h"

#include "Others/def.h"

#include "test/DriverMexProb.h"

using namespace ROPTLIB;

void testFRankESparseApprox(void);

#endif

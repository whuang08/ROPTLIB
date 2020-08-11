/*
This is the test file to run the problem defined in NSOLyapunov.h and NSOLyapunov.cpp.

---- WH
*/

#ifndef TESTSFRQLYAPUNOV_H
#define TESTSFRQLYAPUNOV_H

#include <iostream>
#include "Others/randgen.h"
#include "Manifolds/Manifold.h"
#include "Problems/Problem.h"
#include <ctime>

#include "test/DriverMexProb.h"

#include "Problems/SFRQLyapunov.h"
#include "Manifolds/SphereTx.h"
#include "Manifolds/SymFixedRankQ.h"

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

void testSFRQLyapunov(void);

#endif //TESTSFRQLYAPUNOV_H

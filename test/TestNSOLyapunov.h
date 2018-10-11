/*
This is the test file to run the problem defined in NSOLyapunov.h and NSOLyapunov.cpp.

---- WH
*/

#ifndef TESTNSOLYAPUNOV_H
#define TESTNSOLYAPUNOV_H


#include <iostream>
#include "Others/randgen.h"
#include "Manifolds/Manifold.h"
#include "Problems/Problem.h"
#include "Solvers/SolversLS.h"
#include <ctime>

#include "test/DriverMexProb.h"

#include "Problems/NSOLyapunov/NSOLyapunov.h"
#include "Problems/SphereTxRQ/SphereTxRQ.h"
#include "Manifolds/NStQOrth/NStQOrth.h"
#include "Manifolds/SphereTx/SphereTx.h"

#include "Solvers/RSD.h"
#include "Solvers/RNewton.h"
#include "Solvers/RCG.h"
#include "Solvers/RBroydenFamily.h"
#include "Solvers/RWRBFGS.h"
#include "Solvers/RBFGS.h"
#include "Solvers/LRBFGS.h"

#include "Solvers/SolversTR.h"
#include "Solvers/RTRSD.h"
#include "Solvers/RTRNewton.h"
#include "Solvers/RTRSR1.h"
#include "Solvers/LRTRSR1.h"

#include "Others/def.h"

using namespace ROPTLIB;

//void testfft();
void testNSOLyapunov(void);

#endif //TESTNSOLYAPUNOV_H

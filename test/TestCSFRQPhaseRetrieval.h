/*
This is the test file to run the problem defined in LRBlindDeconvolution.h and LRBlindDeconvolution.cpp.

---- WH
*/

#ifndef TESTCSFRQPHASERETRIEVAL_H
#define TESTCSFRQPHASERETRIEVAL_H


#include <iostream>
#include "Others/randgen.h"
#include "Manifolds/Manifold.h"
#include "Problems/Problem.h"
#include <ctime>

#include "test/DriverMexProb.h"

#include "Problems/CSFRQPhaseRetrieval.h"
#include "Manifolds/SphereTx.h"
#include "Manifolds/CSymFixedRankQ.h"

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

#include "test/DriverMexProb.h"

using namespace ROPTLIB;

#ifdef ROPTLIB_WITH_FFTW

void testCSFRQPhaseRetrieval(void);

#endif
#endif

/*
This is the test file to run the problem defined in LRBlindDeconvolution.h and LRBlindDeconvolution.cpp.

---- WH
*/

#ifndef TESTCSOPHASERETRIEVAL_H
#define TESTCSOPHASERETRIEVAL_H


#include <iostream>
#include "Others/randgen.h"
#include "Manifolds/Manifold.h"
#include "Problems/Problem.h"
#include "Solvers/SolversLS.h"
#include <ctime>

#include "test/DriverMexProb.h"

#include "Problems/CSOPhaseRetrieval/CSOPhaseRetrieval.h"
#include "Problems/SphereTxRQ/SphereTxRQ.h"
#include "Manifolds/CpxNStQOrth/CpxNStQOrth.h"
#include "Manifolds/CpxNStQOrth/CSOVariable.h"
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

#include "Others/fftw/fftw3.h"

using namespace ROPTLIB;

#ifdef ROPTLIB_WITH_FFTW

//void testfft();
void testCSOPhaseRetrieval(void);
bool MyPRstop(Variable *x, Vector *gf, double f, double ngf, double ngf0, const Problem *prob, const Solvers *solver);
void testCSOPhaseRetrievalFixedRank(void);
void CSOPRRankReduce(double *initX, double *b, double *masks, integer n1, integer n2, integer l, integer r, double kappa, integer maxiter, double *outsoln, double &outtime, integer &outnfft, integer &outnn);
void WFegf(double *x, double *b, double *masks, double *egf, integer n1, integer n2, integer l, integer r);
void WFlow(double *initX, double *b, double *masks, integer n1, integer n2, integer l, integer r, integer maxiter, double *outsoln, double &outtime, integer &outnfft);

#endif
#endif

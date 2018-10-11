/*
This is the test file to check all problems.

---- WH
*/

#ifndef DRIVERCPP_H
#define DRIVERCPP_H

#include "test/TestCFR2BlindDecon2D.h"
#include "test/TestCFR2BlindDeconvolution.h"
#include "test/TestCSO.h"
#include "test/TestCSOPhaseRetrieval.h"
#include "test/TestElasticCurvesRO.h"
#include "test/TestEucBlindDeconvolution.h"
#include "test/TestEucFrechetMean.h"
#include "test/TestEucPosSpCd.h"
#include "test/TestEucQuadratic.h"
#include "test/TestGrassRQ.h"
#include "test/TestKarcherMean.h"
#include "test/TestLRBlindDeconvolution.h"
#include "test/TestLRMatrixCompletion.h"
#include "test/TestMyMatrix.h"
#include "test/TestNSOLyapunov.h"
#include "test/TestOrthBoundingBox.h"
#include "test/TestPreShapePathStraighten.h"
#include "test/TestProduct.h"
#include "test/TestShapePathStraighten.h"
#include "test/TestSparsePCA.h"
#include "test/TestSPDMean.h"
#include "test/TestSPDTensorDL.h"
#include "test/TestSphereRayQuo.h"
#include "test/TestSphereSparsestVector.h"
#include "test/TestStieBrockett.h"
#include "test/TestStieSoftICA.h"
#include "test/TestStieSparseBrockett.h"
#include "test/TestTestSparsePCA.h"
#include "test/TestWeightedLowRank.h"

using namespace ROPTLIB;

void testall(void);
int main(void);

#endif

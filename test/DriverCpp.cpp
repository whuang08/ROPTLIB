#include "test/DriverCpp.h"

using namespace ROPTLIB;

int main(void)
{
	//_CrtSetDbgFlag(_CRTDBG_LEAK_CHECK_DF); /*This can detect the memory leakage for global variables!!*/
	//_CrtSetBreakAlloc(781179);
	/*Set the random seed*/
	unsigned tt = (unsigned)time(NULL);
	tt = 0; /*The following test is only for random seed zero*/
	genrandseed(tt);
	//testDSYL();
	testNSOLyapunov();
	//testLRMatrixCompletionMore();
	//testStieBrockettMore();
	//testSparsePCA();
	//testSPDMeanMore();
	//testall();

#ifdef _WIN64
#ifdef _DEBUG
	_CrtDumpMemoryLeaks();
#endif
#endif
	return 0;
}

void testall(void)
{ /*The following test is only for random seed zero*/

  /*Set the random seed*/
	unsigned tt = (unsigned)time(NULL);
	tt = 0; /*The following test is only for random seed zero*/
	genrandseed(tt);
	
	//testDSYL();
	//testCSO();
	//testElasticCurvesRO();
	//testKarcherMean();            //TOCHECK
	//testPreShapePathStraighten(); //TOCHECK
	//testShapePathStraighten();    //TOCHECK
	//testEigenSymmetricM();
	//testExpSymmetricM();
	//testLogSymmetricM();
	//testProduct();

	printf("\n testStieBrockett\n");
	testStieBrockett(); /*Test all methods*/

	printf("\n testEucFrechetMean\n");
	testEucFrechetMean();

	printf("\n testEucPosSpCd\n");
	testEucPosSpCd();

	printf("\n testEucQuadratic\n");
	testEucQuadratic();

	printf("\n testGrassRQ\n");
	testGrassRQ();

	printf("\n testLRMatrixCompletion\n");
	testLRMatrixCompletion();

	printf("\n testOrthBoundingBox\n");
	testOrthBoundingBox();

	printf("\n testSparsePCA\n");
	testSparsePCA();

	printf("\n testSPDMean\n");
	testSPDMean(); /*for dist(x_i, x_*)*/

	printf("\n testSPDTensorDL\n");
	testSPDTensorDL();

	printf("\n testSphereRayQuo\n");
	testSphereRayQuo();

	printf("\n testSphereSparsestVector\n");
	testSphereSparsestVector();

	printf("\n testStieSoftICA\n");
	testStieSoftICA();

	printf("\n testStieSparseBrockett\n");
	testStieSparseBrockett();

	printf("\n testTestSparsePCA\n");
	testTestSparsePCA();

	printf("\n testWeightedLowRank\n");
	testWeightedLowRank();

#ifdef ROPTLIB_WITH_FFTW
	printf("\n testCFR2BlindDecon2D\n");
	testCFR2BlindDecon2D();

	printf("\n testCFR2BlindDeconvolution\n");
	testCFR2BlindDeconvolution();
	testCFR2BlindDeconvolutionSparse();

	printf("\n testCSOPhaseRetrievalFixedRank\n");
	testCSOPhaseRetrievalFixedRank();

	printf("\n testEucBlindDeconvolution\n");
	testEucBlindDeconvolution(); /*smallest/largest Eigenvalues at some point*/
	testEucBlindDeconvolutionSparse();

	printf("\n testLRBlindDeconvolution\n");
	testLRBlindDeconvolution();
	testLRBlindDeconvolutionSparse();
#endif

};

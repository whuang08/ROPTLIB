
#include "test/TestProduct.h"

using namespace ROPTLIB;

void testProduct(void)
{
	integer d = 3;
	integer n = 2;

	/*This is used to generate a Manifold. You can use this manifold to compute operations on this manifold, such as retraction, vector transport....
	This is not used to generate a point on the manifold.*/
	OrthGroup OrthMani(d);
	Euclidean LinearMani(d);
	ProductManifold CartanMani(2, OrthMani, 1, LinearMani, 1);
	ProductManifold ProdMani(1, CartanMani, n);

	/*This is used to generate a point on a manifold*/
	OrthGroupVariable OrthV(d);
	EucVariable EucV(d);
	ProductElement ProdX(2, &OrthV, 1, &EucV, 1); /*two kinds of manifolds. One for each.*/
	ProductElement OProdX(1, &ProdX, n);
	OProdX.RandInManifold();  /*Randomly generate a point on this manifold*/

	OProdX.Print("A point on the product manifold:");

	std::cout << "output numbers in the consecutive memory:" << std::endl;
	const double *Xptr = OProdX.ObtainReadData();
	// The command 
	//      ForDebug::Print("Xptr:", Xptr, d * d * n + n * d);
	// is a short-cut for the following for-loop
	for (integer i = 0; i < d * d * n + d * n; i++)
		std::cout << Xptr[i] << std::endl;


	//ProdMani.CheckIntrExtr(&OProdX);
	//ProdMani.CheckRetraction(&OProdX);
	//ProdMani.CheckDiffRetraction(&OProdX);
	//ProdMani.CheckLockingCondition(&OProdX);
	//ProdMani.CheckcoTangentVector(&OProdX);
	//ProdMani.CheckIsometryofVectorTransport(&OProdX);
	//ProdMani.CheckIsometryofInvVectorTransport(&OProdX);
	//ProdMani.CheckVecTranComposeInverseVecTran(&OProdX);
	//ProdMani.CheckTranHInvTran(&OProdX);
	//ProdMani.CheckHaddScaledRank1OPE(&OProdX);

	//
	//	//integer n = 4;
	//	//L2SphereVariable L2SphereX(n);
	//	//L2SphereX.RandInManifold();
	//	//L2SphereX.Print("X:");//---
	//
	//	//L2Sphere mani(n);
	//	//mani.CheckParams();
	//
	//	////mani.SetHasHHR(true);
	//	////mani.CheckIntrExtr(&L2SphereX);
	//	////mani.CheckRetraction(&L2SphereX);
	//	////mani.CheckDiffRetraction(&L2SphereX, true);
	//	//mani.CheckLockingCondition(&L2SphereX);
	//	////mani.CheckcoTangentVector(&L2SphereX);
	//	////mani.CheckIsometryofVectorTransport(&L2SphereX);
	//	////mani.CheckIsometryofInvVectorTransport(&L2SphereX);
	//	////mani.CheckVecTranComposeInverseVecTran(&L2SphereX);
	//	////mani.CheckTranHInvTran(&L2SphereX);
	//	////mani.CheckHaddScaledRank1OPE(&L2SphereX);
	//
	//
	//	integer n = 3, p = 2, m = 2;
	//	integer numofmanis = 2;
	//	integer numofmani1 = 1;
	//	integer numofmani2 = 1;
	//
	//	StieVariable StieX(n, p);
	////	StieX.RandInManifold();
	//	EucVariable EucX(m);
	////	EucX.RandInManifold();
	//	ProductElement ProdX(numofmanis, &StieX, numofmani1, &EucX, numofmani2);
	//	ProdX.RandInManifold();
	//
	//
	//	//ObliqueVariable ObliqueX(n, p);
	//	//Variable **vars = new Variable *[1];
	//	//vars[0] = &ObliqueX;
	//	//integer *powsinterval = new integer[2];
	//	//powsinterval[0] = 0;
	//	//powsinterval[1] = 1;
	//	//ProductElement ProdX(vars, 1, powsinterval, 1);
	//	//ProdX.Print("ProdX1:");//--
	//	//ProdX.RandInManifold();
	//	//ProdX.Print("ProdX2:");//---
	//	//delete[] vars;
	//	//delete[] powsinterval;
	//
	//	//ProductElement *ProdXptr = ProdX.ConstructEmpty();
	//	//ProdX.CopyTo(ProdXptr);
	//	//ProdXptr->Print("ptr:");
	//
	//	//delete ProdXptr;
	//
	//	Stiefel mani1(n, p);
	//	Euclidean mani2(m);
	//	ProductManifold ProdMani(numofmanis, &mani1, numofmani1, &mani2, numofmani2);
	//	ProdMani.SetHasHHR(true);
	//
	//	ProdMani.CheckIntrExtr(&ProdX);
	//	ProdMani.CheckRetraction(&ProdX);
	//	ProdMani.CheckDiffRetraction(&ProdX);
	//	ProdMani.CheckLockingCondition(&ProdX);
	//	ProdMani.CheckcoTangentVector(&ProdX);
	//	ProdMani.CheckIsometryofVectorTransport(&ProdX);
	//	ProdMani.CheckIsometryofInvVectorTransport(&ProdX);
	//	ProdMani.CheckVecTranComposeInverseVecTran(&ProdX);
	//	ProdMani.CheckTranHInvTran(&ProdX);
	//	ProdMani.CheckHaddScaledRank1OPE(&ProdX);
};

#ifdef MATLAB_MEX_FILE

std::map<integer *, integer> *CheckMemoryDeleted;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	genrandseed(0);

	CheckMemoryDeleted = new std::map<integer *, integer>;
	testProduct();
	std::map<integer *, integer>::iterator iter = CheckMemoryDeleted->begin();
	for (iter = CheckMemoryDeleted->begin(); iter != CheckMemoryDeleted->end(); iter++)
	{
		if (iter->second != 1)
			printf("Global address: %p, sharedtimes: %d\n", iter->first, iter->second);
	}
	delete CheckMemoryDeleted;
	return;
}

#endif

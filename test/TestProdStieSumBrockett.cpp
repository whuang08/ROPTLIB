
#include "test/TestProdStieSumBrockett.h"

using namespace ROPTLIB;

void testProdStieSumBrockett(void)
{
    // choose a random seed
    unsigned tt = (unsigned)time(NULL);
    tt = 0;
    genrandseed(tt);

    // size of the Stiefel manifold
    integer n = 4, p = 2, m = 3, q = 2;

    Vector B1(n, n), B2(n, n), B3(m, m);
    B1.RandGaussian(); B1 = B1 + B1.GetTranspose();
    B2.RandGaussian(); B2 = B2 + B2.GetTranspose();
    B3.RandGaussian(); B3 = B3 + B3.GetTranspose();
    Vector D1(p), D2(p), D3(q);
    realdp *D1ptr = D1.ObtainWriteEntireData();
    realdp *D2ptr = D2.ObtainWriteEntireData();
    realdp *D3ptr = D3.ObtainWriteEntireData();
    
    for (integer i = 0; i < p; i++)
    {
        D1ptr[i] = static_cast<realdp> (i + 1);
        D2ptr[i] = D1ptr[i];
    }
    for (integer i = 0; i < q; i++)
    {
        D3ptr[i] = static_cast<realdp> (i + 1);
    }
    
    // number of manifolds in product of manifold
    integer numoftypes = 2; // two kinds of manifolds
    integer numofmani1 = 2; // the first one has two
    integer numofmani2 = 1; // the second one has one

    // Define the Stiefel manifold
    Stiefel mani1(n, p);
    Stiefel mani2(m, q);
    mani1.ChooseParamsSet1();
    mani2.ChooseParamsSet1();
//    mani1.ChooseParamsSet5();
//    mani2.ChooseParamsSet4();

    ProductManifold Domain(numoftypes, &mani1, numofmani1, &mani2, numofmani2);
//    Domain.SetHasHHR(true);
    // Obtain an initial iterate
    Variable ProdX = Domain.RandominManifold();
//    ProdX.Print("ProdX:");
    
//    Domain.CheckIntrExtr(ProdX);
//    Domain.CheckRetraction(ProdX);
//    Domain.CheckDiffRetraction(ProdX, true);
//    Domain.CheckLockingCondition(ProdX);
//    Domain.CheckcoTangentVector(ProdX);
//    Domain.CheckIsometryofVectorTransport(ProdX);
//    Domain.CheckIsometryofInvVectorTransport(ProdX);
//    Domain.CheckVecTranComposeInverseVecTran(ProdX);
//    Domain.CheckTranHInvTran(ProdX);
//    Domain.CheckHaddScaledRank1OPE(ProdX);

    // Define the Brockett problem
    ProdStieSumBrockett Prob(B1, D1, B2, D2, B3, D3);

    // Set the domain of the problem to be the Stiefel manifold
    Prob.SetDomain(&Domain);

    // output the parameters of the manifold of domain
    Domain.CheckParams();
    
    Prob.SetNumGradHess(true);
    Prob.CheckGradHessian(ProdX);
//    Prob.MinMaxEigValHess(ProdX).Print("ProdX:");
//    return;
    // test LRBFGS
    printf("********************************Check all line search algorithm in LRBFGS*************************************\n");
    RTRNewton *LRBFGSsolver = new RTRNewton(&Prob, &ProdX);
//    LRBFGSsolver->LineSearch_LS = LSSM_ARMIJO;// LSSM_EXACT;
//    LRBFGSsolver->Initstepsize = LSSM_ONESTEP;
    LRBFGSsolver->Verbose = ITERRESULT; //ITERRESULT;//
	LRBFGSsolver->Max_Iteration = 2000;
    LRBFGSsolver->CheckParams();
    LRBFGSsolver->Run();
    Prob.CheckGradHessian(LRBFGSsolver->GetXopt());
    
    delete LRBFGSsolver;
    
//    // test RSD
//    printf("********************************Check all line search algorithm in RSD*****************************************\n");
//    for (integer i = 0; i < INPUTFUN; i++)
//    {
//        RSD *RSDsolver = new RSD(&Prob, &ProdX);
//        RSDsolver->LineSearch_LS = static_cast<LSAlgo> (i);
//        RSDsolver->Debug = FINALRESULT;// FINALRESULT;
//        RSDsolver->CheckParams();
//        RSDsolver->Run();
//        delete RSDsolver;
//    }
//
//    // test RNewton
//    printf("********************************Check all line search algorithm in RNewton*************************************\n");
//    for (integer i = 0; i < INPUTFUN; i++)
//    {
//        RNewton *RNewtonsolver = new RNewton(&Prob, &ProdX);
//        RNewtonsolver->LineSearch_LS = static_cast<LSAlgo> (i);
//        RNewtonsolver->Debug = FINALRESULT;
//        RNewtonsolver->CheckParams();
//        RNewtonsolver->Run();
//        delete RNewtonsolver;
//    }
//
//    // test RCG
//    printf("********************************Check all Formulas in RCG*************************************\n");
//    for (integer i = 0; i < RCGMETHODSLENGTH; i++)
//    {
//        RCG *RCGsolver = new RCG(&Prob, &ProdX);
//        RCGsolver->RCGmethod = static_cast<RCGmethods> (i);
//        RCGsolver->LineSearch_LS = STRONGWOLFE;
//        RCGsolver->LS_beta = 0.1;
//        RCGsolver->Debug = FINALRESULT;
//        RCGsolver->CheckParams();
//        RCGsolver->Run();
//        delete RCGsolver;
//    }
//
//    // test RBroydenFamily
//    printf("********************************Check all Formulas in RCG*************************************\n");
//    for (integer i = 0; i < INPUTFUN; i++)
//    {
//        RBroydenFamily *RBroydenFamilysolver = new RBroydenFamily(&Prob, &ProdX);
//        RBroydenFamilysolver->LineSearch_LS = static_cast<LSAlgo> (i);
//        RBroydenFamilysolver->Debug = FINALRESULT;
//        RBroydenFamilysolver->CheckParams();
//        RBroydenFamilysolver->Run();
//        delete RBroydenFamilysolver;
//    }
//
//    // test RWRBFGS
//    printf("********************************Check all line search algorithm in RWRBFGS*************************************\n");
//    for (integer i = 0; i < INPUTFUN; i++)
//    {
//        RWRBFGS *RWRBFGSsolver = new RWRBFGS(&Prob, &ProdX);
//        RWRBFGSsolver->LineSearch_LS = static_cast<LSAlgo> (i);
//        RWRBFGSsolver->Debug = FINALRESULT; //ITERRESULT;//
//        RWRBFGSsolver->CheckParams();
//        RWRBFGSsolver->Run();
//        delete RWRBFGSsolver;
//    }
//
//    // test RBFGS
//    printf("********************************Check all line search algorithm in RBFGS*************************************\n");
//    for (integer i = 0; i < INPUTFUN; i++)
//    {
//        RBFGS *RBFGSsolver = new RBFGS(&Prob, &ProdX);
//        RBFGSsolver->LineSearch_LS = static_cast<LSAlgo> (i);
//        RBFGSsolver->Debug = FINALRESULT;
//        RBFGSsolver->CheckParams();
//        RBFGSsolver->Run();
//        delete RBFGSsolver;
//    }
//
//    // test LRBFGS
//    printf("********************************Check all line search algorithm in LRBFGS*************************************\n");
//    for (integer i = 0; i < INPUTFUN; i++)//LSALGOLENGTH
//    {
//        LRBFGS *LRBFGSsolver = new LRBFGS(&Prob, &ProdX);
//        LRBFGSsolver->LineSearch_LS = static_cast<LSAlgo> (i);
//        LRBFGSsolver->Debug = FINALRESULT; //ITERRESULT;//
//        LRBFGSsolver->CheckParams();
//        LRBFGSsolver->Run();
//        delete LRBFGSsolver;
//    }
//
//    // test RTRSD
//    printf("********************************Check RTRSD*************************************\n");
//    RTRSD RTRSDsolver(&Prob, &ProdX);
//    RTRSDsolver.Debug = FINALRESULT;
//    RTRSDsolver.CheckParams();
//    RTRSDsolver.Run();
//
//    // test RTRNewton
//    printf("********************************Check RTRNewton*************************************\n");
//    RTRNewton RTRNewtonsolver(&Prob, &ProdX);
//    RTRNewtonsolver.Debug = FINALRESULT;
//    RTRNewtonsolver.CheckParams();
//    RTRNewtonsolver.Run();
//
//    // test RTRSR1
//    printf("********************************Check RTRSR1*************************************\n");
//    RTRSR1 RTRSR1solver(&Prob, &ProdX);
//    RTRSR1solver.Debug = FINALRESULT;
//    RTRSR1solver.CheckParams();
//    RTRSR1solver.Run();
//
//    // test LRTRSR1
//    printf("********************************Check LRTRSR1*************************************\n");
//    LRTRSR1 LRTRSR1solver(&Prob, &ProdX);
//    LRTRSR1solver.Debug = FINALRESULT;
//    LRTRSR1solver.CheckParams();
//    LRTRSR1solver.Run();
//
//    // Check gradient and Hessian
//    Prob.CheckGradHessian(&ProdX);
//    const Variable *xopt = RTRNewtonsolver.GetXopt();
//    Prob.CheckGradHessian(xopt);
//
//    return;
}

/*If it is compiled in Matlab, then the following "mexFunction" is used as the entrance.*/
#ifdef MATLAB_MEX_FILE

/*Help to check the memory leakage problem. No necesary any more.*/
std::map<integer *, integer> *CheckMemoryDeleted;

/*This function checks the number and formats of input parameters.
nlhs: the number of output in mxArray format
plhs: the output objects in mxArray format
nrhs: the number of input in mxArray format
prhs: the input objects in mxArray format */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if(nrhs < 10)
    {
        mexErrMsgTxt("The number of arguments should be at least ten.\n");
    }
    realdp *B1, *D1, *B2, *D2, *B3, *D3, *X;
    B1 = mxGetPr(prhs[0]);
    D1 = mxGetPr(prhs[1]);
    B2 = mxGetPr(prhs[2]);
    D2 = mxGetPr(prhs[3]);
    B3 = mxGetPr(prhs[4]);
    D3 = mxGetPr(prhs[5]);
    X = mxGetPr(prhs[6]);
    /* dimensions of input matrices */
    integer p, n, q, m, HasHHR, Paramset;
    n = mxGetM(prhs[0]);
    p = mxGetM(prhs[1]);
    m = mxGetM(prhs[4]);
    q = mxGetM(prhs[5]);

    /*Check the correctness of the inputs*/
    if(mxGetN(prhs[0]) != n)
    {
        mexErrMsgTxt("The size of matrix B1 is not correct.\n");
    }
    if(mxGetN(prhs[1]) != 1)
    {
        mexErrMsgTxt("The size of the D1 is not correct!\n");
    }
    if(mxGetN(prhs[2]) != n)
    {
        mexErrMsgTxt("The size of matrix B2 is not correct.\n");
    }
    if(mxGetN(prhs[3]) != 1)
    {
        mexErrMsgTxt("The size of the D2 is not correct!\n");
    }
    if(mxGetN(prhs[4]) != m)
    {
        mexErrMsgTxt("The size of matrix B3 is not correct.\n");
    }
    if(mxGetN(prhs[5]) != 1)
    {
        mexErrMsgTxt("The size of the D3 is not correct!\n");
    }
    if(mxGetM(prhs[6]) != 2 * n * p + m * q || mxGetN(prhs[6]) != 1)
    {
        mexErrMsgTxt("The size of the initial X is not correct!\n");
    }
    HasHHR = static_cast<integer> (mxGetScalar(prhs[7]));
    Paramset = static_cast<integer> (mxGetScalar(prhs[8]));

    genrandseed(0);

    CheckMemoryDeleted = new std::map<integer *, integer>;

    integer numoftypes = 2; // two kinds of manifolds
    integer numofmani1 = 2; // the first one has two
    integer numofmani2 = 1; // the second one has one

    Stiefel mani1(n, p);
    Stiefel mani2(m, q);
    if (Paramset == 1)
    {
        mani1.ChooseParamsSet1();
        mani2.ChooseParamsSet1();
    }
    else if (Paramset == 2)
    {
        mani1.ChooseParamsSet2();
        mani2.ChooseParamsSet2();
    }
    else if (Paramset == 3)
    {
        mani1.ChooseParamsSet3();
        mani2.ChooseParamsSet3();
    }
    

    ProductManifold Domain(numoftypes, &mani1, numofmani1, &mani2, numofmani2);
    Variable ProdX = Domain.RandominManifold();
    realdp *ProdXptr = ProdX.ObtainWriteEntireData();
    for(integer i = 0; i < ProdX.Getlength(); i++)
        ProdXptr[i] = X[i];
    
    Vector BB1(n, n), BB2(n, n), BB3(m, m), DD1(p), DD2(p), DD3(q);
    realdp *BB1ptr = BB1.ObtainWriteEntireData();
    realdp *BB2ptr = BB2.ObtainWriteEntireData();
    realdp *BB3ptr = BB3.ObtainWriteEntireData();
    realdp *DD1ptr = DD1.ObtainWriteEntireData();
    realdp *DD2ptr = DD2.ObtainWriteEntireData();
    realdp *DD3ptr = DD3.ObtainWriteEntireData();
    for(integer i = 0; i < n * n; i++)
    {
        BB1ptr[i] = B1[i];
        BB2ptr[i] = B2[i];
    }
    for(integer i = 0; i < m * m; i++)
        BB3ptr[i] = B3[i];
    for(integer i = 0; i < p; i++)
    {
        DD1ptr[i] = D1[i];
        DD2ptr[i] = D2[i];
    }
    for(integer i = 0; i < q; i++)
        DD3ptr[i] = D3[i];

    ProdStieSumBrockett Prob(BB1, DD1, BB2, DD2, BB3, DD3);

    Prob.SetDomain(&Domain);
    
    Domain.SetHasHHR(HasHHR != 0);

    // Call the function defined in DriverMexProb.h
    ParseSolverParamsAndOptimizing(prhs[9], &Prob, &ProdX, plhs);

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

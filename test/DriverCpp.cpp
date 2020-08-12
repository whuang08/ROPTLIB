#include "test/DriverCpp.h"

using namespace ROPTLIB;

int main(void)
{
//	_CrtSetDbgFlag(_CRTDBG_LEAK_CHECK_DF); /*This can detect the memory leakage for global variables!!*///
//	_CrtSetBreakAlloc(160);
	/*Set the random seed*/
	unsigned tt = (unsigned)time(NULL);
    tt = 0; /*The following test is only for random seed zero*/
	std::cout << "seed:" << tt << std::endl;
	genrandseed(tt);
    
//    testElement();
//    testEucQuadratic();
//    testCFRankQ2FBlindDecon2D();
//    testCStieBrockett();
//    testFRankQ2FMatCompletion();
//    testFRankESparseApprox();
//    testFRankEWeightApprox();
//    testGrassRQ();
//    testProdStieSumBrockett();
//    testSPDKarcherMean();
//    testSphereSparsestVector();
//    testStieBrockett();
//    testStieSPCA();
//    testSFRQLyapunov();
//    testCSFRQPhaseRetrieval();

//    integer length = 2;
//    realdp *vals = new realdp [2];
//    vals[0] = 1; vals[1] = 1;
//    integer inc = 1;
//    std::cout << ":::" << sdot_(&length, vals, &inc, vals, &inc) << std::endl;//---
    
	testall();

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
	unsigned seed = (unsigned)time(NULL);
	seed = 1; /*The following test is only for random seed zero*/

#ifdef ROPTLIB_WITH_FFTW
    {
        genrandseed(seed);
        printf("\n testCSFRQPhaseRetrieval\n");

        integer n1 = 2, n2 = 2, r = 1, l = 2;
        integer n = n1 * n2, m = n * l;
        realdp kappa = 0.000;

        CSymFixedRankQ Domain(n, r);
        Vector InitialX = 3 * (Domain.RandominManifold());
        Vector masks(n, l, "complex"); masks.RandGaussian();
        Vector xtrue(n, r, "complex"); xtrue.RandGaussian();
        Vector ZY(m, r, "complex");
        realdp sqn = sqrt(static_cast<realdp> (n1 * n2));
        Vector tmpv;

        for(integer i = 0; i < l; i++)
        {
            for(integer j = 0; j < r; j++)
            {
                tmpv = (masks.GetSubmatrix(0, n - 1, i, i).GetHadamardProduct(xtrue.GetSubmatrix(0, n - 1, j, j))) / sqn;
                ZY.SubmatrixAssignment(i * n, (i + 1) * n - 1, j, j, tmpv.Reshape(n1, n2).GetFFT2D(FFTW_FORWARD).Reshape(n, 1));
            }
        }
        Vector b = ZY.GetTranspose().GetColNormsSquare().GetTranspose();

        CSFRQPhaseRetrieval Prob(b, masks, kappa, n1, n2, l, r);
        Prob.SetDomain(&Domain);

        std::vector<std::string> names = {"RSD", "RNewton", "RCG", "RBroydenFamily", "RBFGS", "LRBFGS", "RTRNewton", "RTRSR1", "LRTRSR1"}; /*"RTRSD", */
        testSmoothProblem(&Prob, &InitialX, "testCSFRQPhaseRetrieval", names);
    }

    {
        genrandseed(seed);
        printf("\n testCFRankQ2FBlindDecon2D\n");

        integer n1 = 2, n2 = 2, r = 1, L = n1 * n2;
        CFixedRankQ2F Domain(L, L, r);
        Domain.SetHasHHR(false);
        Variable InitialX = Domain.RandominManifold();

        Vector y(L, "complex");
        y.RandGaussian();
        // Generate the matrices in the Low rank approximation problem.
        realdp *B = new realdp[L * L * 2 + L * L * 2];
        realdp *C = B + L * L * 2;
        for (integer i = 0; i < L * L * 2 + L * L * 2; i++)
            B[i] = genrandnormal();
        integer nzmaxB = L * L;
        integer nzmaxC = L * L;
        integer *irB = new integer[2 * L * L + 2 * L * L];
        integer *jcB = irB + L * L;
        integer *irC = jcB + L * L;
        integer *jcC = irC + L * L;
        for (integer i = 0; i < L; i++)
        {
            for (integer j = 0; j < L; j++)
            {
                irB[j + i * L] = j;
                jcB[j + i * L] = i;
            }
        }
        for (integer i = 0; i < L; i++)
        {
            for (integer j = 0; j < L; j++)
            {
                irC[j + i * L] = j;
                jcC[j + i * L] = i;
            }
        }

        SparseMatrix sB(L, L, irB, jcB, (realdpcomplex *) B, nzmaxB);
        SparseMatrix sC(L, L, irC, jcC, (realdpcomplex *) C, nzmaxC);

        CFRankQ2FBlindDecon2D Prob(y, sB, sC, n1, n2, r, 0, 1, 1);
        Prob.SetDomain(&Domain);

        std::vector<std::string> names = {"RSD", "RNewton", "RCG", "RBroydenFamily", "RBFGS", "LRBFGS", "RTRSD", "RTRNewton", "RTRSR1", "LRTRSR1"};

        testSmoothProblem(&Prob, &InitialX, "testCFRankQ2FBlindDecon2D", names);

        delete[] B;
        delete[] irB;
    }
#endif

    {
        genrandseed(seed);
        printf("\n testSPDKarcherMean\n");

        /*Randomly generate a point on the SPD manifold*/
        integer n = 3, num = 3;

        // Define the manifold
        SPDManifold Domain(n);
        Domain.ChooseParamsSet1();
        Variable InitialX = Domain.RandominManifold();

        Vector EE(n, n), tmp(n, n);
        Vector Ls(1, &EE, num);
        Ls.NewMemoryOnWrite();
        for(integer i = 0; i < num; i++)
        {
            tmp.RandGaussian(); tmp.QRDecom();
            Ls.GetElement(i) = tmp.Field("_R").GetTranspose();
        }
        
        // Define the problem
        SPDKarcherMean Prob(Ls, n, num);
        /*The domain of the problem is a SPD manifold*/
        Prob.SetDomain(&Domain);
        
        std::vector<std::string> names = {"RSD", "RNewton", "RCG", "RBroydenFamily", "RBFGS", "LRBFGS", "RTRSD", "RTRNewton", "RTRSR1", "LRTRSR1"};
        testSmoothProblem(&Prob, &InitialX, "testSPDKarcherMean", names);
    }

    {
        genrandseed(seed);
        printf("\n testGrassRQ\n");

        /* size of the Grassmann manifold */
        integer n = 5, p = 2;
        /* Generate the matrices in the Rayleigh Quotient problem. */
        Vector B(n, n);
        B.RandGaussian();
        B = B + B.GetTranspose();

        // Define the manifold
        Grassmann Domain(n, p);
        Variable InitialX = Domain.RandominManifold();

        /* Define the  problem*/
        GrassRQ Prob(B, n, p);
        /*The domain of the problem is a Stiefel manifold*/
        Prob.SetDomain(&Domain);

        std::vector<std::string> names = {"RSD", "RNewton", "RCG", "RBroydenFamily", "RWRBFGS", "RBFGS", "LRBFGS", "RTRSD", "RTRNewton", "RTRSR1", "LRTRSR1"};
        testSmoothProblem(&Prob, &InitialX, "testGrassRQ", names);
    }

    {
        genrandseed(seed);
        printf("\n testSFRQLyapunov\n");

        integer n = 5, p = 2, pC = 1;
        integer nzmaxA = n + 2 * (n - 1);
        realdp *A = new realdp[n + 2 * (n - 1)];
        integer *inirA = new integer[2 * nzmaxA];
        integer *injcA = inirA + nzmaxA;
        for (integer i = 0; i < n; i++)
        {
            A[i] = 2;
            inirA[i] = i;
            injcA[i] = i;
        }
        for (integer i = 0; i < n - 1; i++)
        {
            A[n + i] = -1;
            inirA[n + i] = i;
            injcA[n + i] = i + 1;
            A[2 * n - 1 + i] = -1;
            inirA[2 * n - 1 + i] = i + 1;
            injcA[2 * n - 1 + i] = i;
        }
        SparseMatrix sA(n, n, inirA, injcA, A, nzmaxA);

        integer nzmaxM = n;
        realdp *M = new realdp[n];
        integer *inirM = new integer[2 * nzmaxM];
        integer *injcM = inirM + nzmaxM;
        for(integer i = 0; i < n; i++)
        {
            M[i] = 1;
            inirM[i] = i;
            injcM[i] = i;
        }
        SparseMatrix sM(n, n, inirM, injcM, M, nzmaxM);

        Vector C(n, pC);
        C.RandGaussian();

        SymFixedRankQ Domain(n, p);
        Vector InitialX = Domain.RandominManifold();
        SFRQLyapunov Prob(sA, sM, C, p);
        Prob.SetDomain(&Domain);

        std::vector<std::string> names = {"RSD", "RNewton", "RCG", "RBroydenFamily", "RBFGS", "LRBFGS", "RTRSD", "RTRNewton", "RTRSR1", "LRTRSR1"};
        testSmoothProblem(&Prob, &InitialX, "testSFRQLyapunov", names);
		delete[] inirM;
		delete[] M;
		delete[] inirA;
		delete[] A;
    }

    {
        genrandseed(seed);
        printf("\n testFRankQ2FMatCompletion\n");

        integer m = 8, n = 7, r = 2;
        FixedRankQ2F Domain(m, n, r);
        Variable InitialX = Domain.RandominManifold();
        Vector G(m, r); G.RandGaussian(); G.QRDecom(); InitialX.GetElement(0) = G.Field("_Q");
        Vector H(n, r); H.RandGaussian(); H.QRDecom(); InitialX.GetElement(1) = H.Field("_Q");
        // Generate the matrices in the matrix completion approximation problem.
        integer dim = (m + n - r) * r;
        integer nz = 2 * dim;
        integer *ir = new integer[nz * 2];
        integer *jc = ir + nz;
        integer *tmpforidx = new integer[m * n];
        for (integer i = 0; i < m * n; i++)
            tmpforidx[i] = i;
        /*nz number of indices*/
        integer idx = 0, itmp;
        for (integer i = 0; i < nz; i++)
        {
            /*idx is an integer in [0, m - i - 1]*/
            idx = static_cast<integer> ((m * n - i) * genrandreal());
            while (idx >= m * n - i)
                idx = static_cast<integer> ((m * n - i) * genrandreal());
            /*the chosen idx is put at the end of the array*/
            itmp = tmpforidx[m * n - i - 1];
            tmpforidx[m * n - i - 1] = tmpforidx[idx];
            tmpforidx[idx] = itmp;
        }
        for (integer i = 0; i < nz; i++)
        {
            /*tmpforidx[nz - 1 - i]*/
            ir[i] = static_cast<integer> (tmpforidx[nz - 1 - i] / n);
            jc[i] = tmpforidx[nz - 1 - i] - n * ir[i];
        }
        delete[] tmpforidx;
        integer mr = m * r, nr = n * r;
        realdp *A_G = new realdp[mr];
        realdp *A_H = new realdp[nr];
        for (integer i = 0; i < m * r; i++)
        {
            A_G[i] = genrandnormal();
        }
        for (integer i = 0; i < n * r; i++)
        {
            A_H[i] = genrandnormal();
        }
        realdp *A_vals = new realdp[nz];
        for (integer i = 0; i < nz; i++)
        {
            A_vals[i] = 0;
            for (integer j = 0; j < r; j++)
            {
                A_vals[i] += A_G[ir[i] + j * m] * A_H[jc[i] + j * n];
            }
        }
        delete[] A_G;
        delete[] A_H;
        FRankQ2FMatCompletion Prob(ir, jc, A_vals, nz, m, n, r);
        Prob.SetDomain(&Domain);

        std::vector<std::string> names = {"RSD", "RNewton", "RCG", "RBroydenFamily", "RBFGS", "LRBFGS", "RTRSD", "RTRNewton", "RTRSR1", "LRTRSR1"};
        testSmoothProblem(&Prob, &InitialX, "testFRankQ2FMatCompletion", names);
		delete[] A_vals;
		delete[] ir;
    }

    {
        genrandseed(seed);
        printf("\n testStieBrockett\n");
        // size of the Stiefel manifold
        integer n = 5, p = 3;
        // Generate the matrices in the Brockett problem.
        Vector B(n, n), D(p);
        B.RandGaussian();
        B = B + B.GetTranspose();
        realdp *Dptr = D.ObtainWriteEntireData();
        /*D is a diagonal matrix.*/
        for (integer i = 0; i < p; i++)
            Dptr[i] = static_cast<realdp> (i + 1);
        // Define the manifold
        Stiefel Domain(n, p);
//        Domain.ChooseParamsSet2();
        //Grassmann Domain(n, p);
//        Domain.SetHasHHR(true); /*set whether the manifold uses the idea in [HGA2015, Section 4.3] or not*/
        Variable StieX = Domain.RandominManifold();
        // Define the Brockett problem
        StieBrockett Prob(B, D);
        /*The domain of the problem is a Stiefel manifold*/
        Prob.SetDomain(&Domain);
        std::vector<std::string> names = {"RSD", "RNewton", "RCG", "RBroydenFamily", "RWRBFGS", "RBFGS", "LRBFGS", "RTRSD", "RTRNewton", "RTRSR1", "LRTRSR1"};
        testSmoothProblem(&Prob, &StieX, "testStieBrockett", names);
    }

    {
        genrandseed(seed);
        printf("\n testProdStieSumBrockett\n");
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
        ProductManifold Domain(numoftypes, &mani1, numofmani1, &mani2, numofmani2);
        // Obtain an initial iterate
        Variable ProdX = Domain.RandominManifold();
        // Define the Brockett problem
        ProdStieSumBrockett Prob(B1, D1, B2, D2, B3, D3);
        // Set the domain of the problem to be the Stiefel manifold
        Prob.SetDomain(&Domain);
        std::vector<std::string> names = {"RSD", "RNewton", "RCG", "RBroydenFamily", "RWRBFGS", "RBFGS", "LRBFGS", "RTRSD", "RTRNewton", "RTRSR1", "LRTRSR1"};
        testSmoothProblem(&Prob, &ProdX, "testProdStieSumBrockett", names);
    }

    {
        genrandseed(seed);
        printf("\n testEucQuadratic\n");
        // size of the domain
        integer dim = 10;
        Vector O(dim, dim);
        O.RandGaussian();
        O = O.GetOrth();
        Vector D(dim);
        D.RandUnform();
        D = D + 0.1;
        Vector A = O.GetTranspose() * D.GetDiagTimesM(O);
        // Obtain an initial iterate
        Euclidean EucDomain(dim);
        Variable EucX = EucDomain.RandominManifold();
        // Define the problem
        EucQuadratic Prob(A);
        Prob.SetDomain(&EucDomain);
        std::vector<std::string> names = {"RSD", "RNewton", "RCG", "RBroydenFamily", "RWRBFGS", "RBFGS", "LRBFGS", "RTRSD", "RTRNewton", "RTRSR1", "LRTRSR1"};
        testSmoothProblem(&Prob, &EucX, "testEucQuadratic", names);
    }

    {
        genrandseed(seed);
        printf("\n testFRankEWeightApprox\n");
        integer m = 5, n = 4, r = 2;
        //integer m = 100, n = 15, r = 5;
        FixedRankE Domain(m, n, r);
        Domain.SetHasHHR(true);
        Variable InitialX = Domain.RandominManifold();
        // Generate the matrices in the Low rank approximation problem.
        Vector A(m, n); A.RandGaussian();
        Vector O(m * n, m * n);
        O.RandGaussian();
        O = O.GetOrth();
        Vector D(m * n);
        D.RandUnform();
        D = D + 0.1;
        Vector W = O.GetTranspose() * D.GetDiagTimesM(O);

        FRankEWeightApprox Prob(A, W, m, n, r);
        Prob.SetDomain(&Domain);
        /* Extrinsic approach is used here and transporting linear opeartor has not been done. Therefore, all algorithms
        that need transporting operators are not tested here. */
        std::vector<std::string> names = {"RSD", "RNewton", "RCG", "LRBFGS", "RTRSD", "RTRNewton", "LRTRSR1"};
        testSmoothProblem(&Prob, &InitialX, "testFRankEWeightApprox", names);
    }

    {
        genrandseed(seed);
        printf("\n testFRankESparseApprox\n");
        integer m = 5, n = 5, r = 2;
        //integer m = 100, n = 15, r = 5;
        FixedRankE Domain(m, n, r);
        Domain.SetHasHHR(false);
        Variable InitialX = Domain.RandominManifold();

        Vector A(m, n); A.RandGaussian();
        realdp lambda = 0.1;
        integer lengthW = 1;

        FRankESparseApprox Prob(A, lambda, m, n, r, lengthW);
        Prob.SetDomain(&Domain);

        std::vector<std::string> names = {"AManPG", "ManPG"};
        testProxGradProblem(&Prob, &InitialX, "testFRankESparseApprox", names);
    }

    {
        genrandseed(seed);
        printf("\n testStieSPCA\n");
        // size of the Stiefel manifold
        integer n = 5, m = 4, p = 3;
        realdp lambda = 1;
        Vector B(m, n);
        realdp *Bptr = B.ObtainWriteEntireData();
        /*B is an n by n matrix*/
        for (integer i = 0; i < n * m; i++)
        {
            Bptr[i] = genrandnormal();
        }
        for(integer i = 0; i < n; i++)
        {
            realdp s = 0;
            for(integer j = 0; j < m; j++)
            {
                s += Bptr[j + i * m];
            }
            s /= m;
            for(integer j = 0; j < m; j++)
                Bptr[j + i * m] -= s;
            s = 0;
            for(integer j = 0; j < m; j++)
            {
                s += Bptr[j + i * m] * Bptr[j + i * m];
            }
            s = std::sqrt(s);
            for(integer j = 0; j < m; j++)
                Bptr[j + i * m] /= s;
        }
        Stiefel Domain(n, p);
        Domain.ChooseParamsSet4();
        Variable StieX = Domain.RandominManifold();
        integer lengthW = 1;
        StieSPCA Prob(B, lambda, n, m, p, lengthW);
        /*The domain of the problem is a Stiefel manifold*/
        Prob.SetDomain(&Domain);
        std::vector<std::string> names = {"AManPG", "ManPG"};
        testProxGradProblem(&Prob, &StieX, "testStieSPCA", names);
    }

    {
        genrandseed(seed);
        printf("\n testSphereSparsestVector\n");
        // size of the matrix Q
        integer m = 10, n = 3;
        // Generate the matrix
        Vector Q(m, n);
        Q.RandGaussian();
        // Define the manifold
        Sphere Domain(n);
        Domain.SetHasHHR(true);
        Variable SphereX = Domain.RandominManifold();
        //Domain.SetHasHHR(true); /*set whether the manifold uses the idea in [HGA2015, Section 4.3] or not*/
        // Define the SparestVector problem
        SphereSparsestVector Prob(Q);
        /*The domain of the problem is a Stiefel manifold*/
        Prob.SetDomain(&Domain);
        std::vector<std::string> names = {"LRBFGSSub", "RBFGSSub", "RGS"};
        testSubGradProblem(&Prob, &SphereX, "testSphereSparsestVector", names);
    }
};

void testSubGradProblem(Problem *prob, Variable *initx, const char *probname, std::vector<std::string> Methodnames)
{
    #ifdef SINGLE_PRECISION
        realdp tol = 5e-2;
    #else
        realdp tol = 1e-6;
    #endif
    if(stringinclude(Methodnames, "LRBFGSSub"))
    {
        printf("\n********************************Check LRBFGSSub in %s*************************************\n", probname);
        LRBFGSSub *LRBFGSSubsolver = new LRBFGSSub(prob, initx);
        LRBFGSSubsolver->Verbose = FINALRESULT;
        LRBFGSSubsolver->Tolerance = tol;
        LRBFGSSubsolver->Max_Iteration = 200;
    //    LRBFGSSubsolver->CheckParams();
        LRBFGSSubsolver->Run();
        if (LRBFGSSubsolver->Getnormndnd0() < tol)
            printf("SUCCESS!\n");
        else
            printf("FAIL!\n");
        delete LRBFGSSubsolver;
    }
    
    if(stringinclude(Methodnames, "RBFGSSub"))
    {
        printf("\n********************************Check RBFGSSub in %s*************************************\n", probname);
        RBFGSSub *RBFGSSubsolver = new RBFGSSub(prob, initx);
        RBFGSSubsolver->Verbose = FINALRESULT;
        RBFGSSubsolver->Tolerance = tol;
        RBFGSSubsolver->Max_Iteration = 200;
    //    RBFGSSubsolver->CheckParams();
        RBFGSSubsolver->Run();
        if (RBFGSSubsolver->Getnormndnd0() < tol)
            printf("SUCCESS!\n");
        else
            printf("FAIL!\n");
        delete RBFGSSubsolver;
    }
    
    if(stringinclude(Methodnames, "RGS"))
    {
        printf("\n********************************Check RGS in %s*************************************\n", probname);
        RGS *RGSsolver = new RGS(prob, initx);
        RGSsolver->Verbose = FINALRESULT;
        RGSsolver->Tolerance = tol;
        RGSsolver->Max_Iteration = 200;
    //    RGSsolver->CheckParams();
        RGSsolver->Run();
        if (RGSsolver->Getnormndnd0() < tol)
            printf("SUCCESS!\n");
        else
            printf("FAIL!\n");
        delete RGSsolver;
    }
};

void testProxGradProblem(Problem *prob, Variable *initx, const char *probname, std::vector<std::string> Methodnames)
{
    #ifdef SINGLE_PRECISION
        realdp tol = 1e-1;
    #else
        realdp tol = 1e-4;
    #endif
    if(stringinclude(Methodnames, "AManPG"))
    {
        printf("\n********************************Check AManPG in %s*************************************\n", probname);
        AManPG *AManPGsolver = new AManPG(prob, initx);
        AManPGsolver->Max_Iteration = 500;
        AManPGsolver->Verbose = FINALRESULT;
        AManPGsolver->Tolerance = tol;
    //    AManPGsolver->Variant = LSPG_REGULAR; //-- LSPG_REGULAR; //-- LSPG_ADALIPSCHITZ;
    //    AManPGsolver->CheckParams();
        AManPGsolver->Run();
        
        if (AManPGsolver->Getnormndnd0() < tol)
            printf("SUCCESS!\n");
        else
            printf("FAIL!\n");

        delete AManPGsolver;
    }
    
    if(stringinclude(Methodnames, "ManPG"))
    {
        printf("\n********************************Check ManPG in %s*************************************\n", probname);
        ManPG *ManPGsolver = new ManPG(prob, initx);
        ManPGsolver->Max_Iteration = 500;
        ManPGsolver->Verbose = FINALRESULT;
        ManPGsolver->Tolerance = tol;
    //    ManPGsolver->Variant = LSPG_ADALIPSCHITZ; //-- LSPG_REGULAR; //-- LSPG_ADALIPSCHITZ;
    //    ManPGsolver->CheckParams();
        ManPGsolver->Run();
        
        if (ManPGsolver->Getnormndnd0() < tol)
            printf("SUCCESS!\n");
        else
            printf("FAIL!\n");

        delete ManPGsolver;
    }
};

void testSmoothProblem(Problem *Prob, Variable *initx, const char *probname, std::vector<std::string> Methodnames)
{
#ifdef SINGLE_PRECISION
    realdp tol = 5e-2;
#else
    realdp tol = 1e-6;
#endif
    // test RSD
    if(stringinclude(Methodnames, "RSD"))
    {
        integer RSDMI[5] = { 2000, 2000, 2000, 2000, 2000 };
        printf("********************************Check all line search algorithm in RSD in %s*****************************************\n", probname);
        for (integer i = 0; i < LSSM_INPUTFUN; i++) //LSSM_INPUTFUN
        {
            Prob->SetNumGradHess(false);
            RSD *RSDsolver = new RSD(Prob, initx);
            RSDsolver->LineSearch_LS = static_cast<LSAlgoSM> (i);
            RSDsolver->Verbose = FINALRESULT;
            RSDsolver->Max_Iteration = RSDMI[i];
            RSDsolver->Run();
            if (RSDsolver->Getnormgfgf0() < tol)
                printf("SUCCESS!\n");
            else
                printf("FAIL!\n");
            delete RSDsolver;
        }
        printf("********************************Check all line search algorithm in RSD in %s with numerical gradient and Hessian*****************************************\n", probname);
        for (integer i = 0; i < LSSM_INPUTFUN; i++)
        {
            Prob->SetNumGradHess(true);
            RSD *RSDsolver = new RSD(Prob, initx);
            RSDsolver->LineSearch_LS = static_cast<LSAlgoSM> (i);
            RSDsolver->Verbose = FINALRESULT;
            RSDsolver->Max_Iteration = RSDMI[i];
            RSDsolver->Run();
            if (RSDsolver->Getnormgfgf0() < tol)
                printf("SUCCESS!\n");
            else
                printf("FAIL!\n");
            delete RSDsolver;
        }
    }
    
    // test RNewton
    if(stringinclude(Methodnames, "RNewton"))
    {
        integer RNewtonMI[5] = { 80, 80, 80, 80, 80 };
        printf("\n********************************Check all line search algorithm in RNewton in %s*************************************\n", probname);
        for (integer i = 0; i < LSSM_INPUTFUN; i++) // LSALGOLENGTH
        {
            Prob->SetNumGradHess(false);
            RNewton *RNewtonsolver = new RNewton(Prob, initx);
            RNewtonsolver->LineSearch_LS = static_cast<LSAlgoSM> (i);
            RNewtonsolver->Verbose = FINALRESULT;
            /*Uncomment following two lines to use the linesearch algorithm defined by the function "LinesearchInput".*/
            //RNewtonsolver->LineSearch_LS = INPUTFUN;
            //RNewtonsolver->LinesearchInput = &LinesearchInput;
			RNewtonsolver->Max_Iteration = RNewtonMI[i];
            //RNewtonsolver->CheckParams();
            RNewtonsolver->Run();

            if (RNewtonsolver->Getnormgfgf0() < tol)
                printf("SUCCESS!\n");
            else
                printf("FAIL!\n");

            delete RNewtonsolver;
        }
        printf("\n********************************Check all line search algorithm in RNewton in %s with numerical gradient and Hessian*************************************\n", probname);
        for (integer i = 0; i < LSSM_INPUTFUN; i++) //LSSM_INPUTFUN  LSALGOLENGTH
        {
            Prob->SetNumGradHess(true);
            RNewton *RNewtonsolver = new RNewton(Prob, initx);
            RNewtonsolver->LineSearch_LS = static_cast<LSAlgoSM> (i);
            RNewtonsolver->Verbose = FINALRESULT;
            /*Uncomment following two lines to use the linesearch algorithm defined by the function "LinesearchInput".*/
            //RNewtonsolver->LineSearch_LS = INPUTFUN;
            //RNewtonsolver->LinesearchInput = &LinesearchInput;
            RNewtonsolver->Max_Iteration = RNewtonMI[i];
            //RNewtonsolver->CheckParams();
            RNewtonsolver->Run();

            if (RNewtonsolver->Getnormgfgf0() < tol)
                printf("SUCCESS!\n");
            else
                printf("FAIL!\n");

            delete RNewtonsolver;
        }
    }

    // test RCG
    if(stringinclude(Methodnames, "RCG"))
    {
        integer RCGMI[6] = { 2000, 2000, 2000, 2000, 2000, 2000 };
        printf("\n********************************Check all Formulas in RCG in %s*************************************\n", probname);
        for (integer i = 0; i < RCGMETHODSLENGTH; i++)
        {
            Prob->SetNumGradHess(false);
            RCG *RCGsolver = new RCG(Prob, initx);
            RCGsolver->RCGmethod = static_cast<RCGmethods> (i);
            //RCGsolver->LineSearch_LS = LSSM_ARMIJO;
            //RCGsolver->InitSteptype = QUADINTMOD;
            RCGsolver->Verbose = FINALRESULT;
            RCGsolver->Max_Iteration = RCGMI[i];
            //RCGsolver->CheckParams();
            RCGsolver->Run();
            if (RCGsolver->Getnormgfgf0() < tol)
                printf("SUCCESS!\n");
            else
                printf("FAIL!\n");
            delete RCGsolver;
        }
        printf("\n********************************Check all Formulas in RCG in %s with numerical gradient and Hessian*************************************\n", probname);
        for (integer i = 0; i < RCGMETHODSLENGTH; i++)
        {
            Prob->SetNumGradHess(true);
            RCG *RCGsolver = new RCG(Prob, initx);
            RCGsolver->RCGmethod = static_cast<RCGmethods> (i);
            //RCGsolver->LineSearch_LS = LSSM_ARMIJO;
            //RCGsolver->InitSteptype = QUADINTMOD;
            RCGsolver->Verbose = FINALRESULT;
            RCGsolver->Max_Iteration = RCGMI[i];
            //RCGsolver->CheckParams();
            RCGsolver->Run();
            if (RCGsolver->Getnormgfgf0() < tol)
                printf("SUCCESS!\n");
            else
                printf("FAIL!\n");
            delete RCGsolver;
        }
    }
    
    // test RBroydenFamily
    if(stringinclude(Methodnames, "RBroydenFamily"))
    {
        integer RBroydenFamilyMI[5] = { 200, 200, 200, 200, 200 };
        printf("\n********************************Check all Formulas in RBroydenFamily in %s*************************************\n", probname);
        for (integer i = 0; i < LSSM_INPUTFUN; i++)
        {
            Prob->SetNumGradHess(false);
            RBroydenFamily *RBroydenFamilysolver = new RBroydenFamily(Prob, initx);
            RBroydenFamilysolver->LineSearch_LS = static_cast<LSAlgoSM> (i);
            RBroydenFamilysolver->Verbose = FINALRESULT;
            RBroydenFamilysolver->Max_Iteration = RBroydenFamilyMI[i];
//            RBroydenFamilysolver->CheckParams();
            RBroydenFamilysolver->Run();
            if (RBroydenFamilysolver->Getnormgfgf0() < tol)
                printf("SUCCESS!\n");
            else
                printf("FAIL!\n");
            delete RBroydenFamilysolver;
        }
        printf("\n********************************Check all Formulas in RBroydenFamily in %s with numerical gradient and Hessian*************************************\n", probname);
        for (integer i = 0; i < LSSM_INPUTFUN; i++)
        {
            Prob->SetNumGradHess(true);
            RBroydenFamily *RBroydenFamilysolver = new RBroydenFamily(Prob, initx);
            RBroydenFamilysolver->LineSearch_LS = static_cast<LSAlgoSM> (i);
            RBroydenFamilysolver->Verbose = FINALRESULT;
            RBroydenFamilysolver->Max_Iteration = RBroydenFamilyMI[i];
            //RBroydenFamilysolver->CheckParams();
            RBroydenFamilysolver->Run();
            if (RBroydenFamilysolver->Getnormgfgf0() < tol)
                printf("SUCCESS!\n");
            else
                printf("FAIL!\n");
            delete RBroydenFamilysolver;
        }
    }
    
    // test RWRBFGS
    if(stringinclude(Methodnames, "RWRBFGS"))
    {
        integer RWRBFGSMI[5] = { 200, 200, 200, 200, 200 };
        printf("\n********************************Check all line search algorithm in RWRBFGS in %s*************************************\n", probname);
        for (integer i = 0; i < LSSM_INPUTFUN; i++)
        {
            Prob->SetNumGradHess(false);
            RWRBFGS *RWRBFGSsolver = new RWRBFGS(Prob, initx);
            RWRBFGSsolver->LineSearch_LS = static_cast<LSAlgoSM> (i);
            RWRBFGSsolver->Verbose = FINALRESULT; //ITERRESULT;//
            RWRBFGSsolver->Max_Iteration = RWRBFGSMI[i];
            //RWRBFGSsolver->CheckParams();
            RWRBFGSsolver->Run();
            if (RWRBFGSsolver->Getnormgfgf0() < tol)
                printf("SUCCESS!\n");
            else
                printf("FAIL!\n");
            delete RWRBFGSsolver;
        }
        printf("\n********************************Check all line search algorithm in RWRBFGS in %s with numerical gradient and Hessian*************************************\n", probname);
        for (integer i = 0; i < LSSM_INPUTFUN; i++)
        {
            Prob->SetNumGradHess(true);
            RWRBFGS *RWRBFGSsolver = new RWRBFGS(Prob, initx);
            RWRBFGSsolver->LineSearch_LS = static_cast<LSAlgoSM> (i);
            RWRBFGSsolver->Verbose = FINALRESULT; //ITERRESULT;//
            RWRBFGSsolver->Max_Iteration = RWRBFGSMI[i];
            //RWRBFGSsolver->CheckParams();
            RWRBFGSsolver->Run();
            if (RWRBFGSsolver->Getnormgfgf0() < tol)
                printf("SUCCESS!\n");
            else
                printf("FAIL!\n");
            delete RWRBFGSsolver;
        }
    }
    
    // test RBFGS
    if(stringinclude(Methodnames, "RBFGS"))
    {
        integer RBFGSMI[5] = { 200, 200, 200, 200, 200 };
        printf("\n********************************Check all line search algorithm in RBFGS in %s*************************************\n", probname);
        for (integer i = 0; i < LSSM_INPUTFUN; i++)
        {
            Prob->SetNumGradHess(false);
            RBFGS *RBFGSsolver = new RBFGS(Prob, initx);
            RBFGSsolver->LineSearch_LS = static_cast<LSAlgoSM> (i);
            RBFGSsolver->Verbose = FINALRESULT;
            RBFGSsolver->Max_Iteration = RBFGSMI[i];
            //RBFGSsolver->CheckParams();
            RBFGSsolver->Run();
            if (RBFGSsolver->Getnormgfgf0() < tol)
                printf("SUCCESS!\n");
            else
                printf("FAIL!\n");
            delete RBFGSsolver;
        }
        printf("\n********************************Check all line search algorithm in RBFGS in %s with numerical gradient and Hessian*************************************\n", probname);
        for (integer i = 0; i < LSSM_INPUTFUN; i++)
        {
            Prob->SetNumGradHess(true);
            RBFGS *RBFGSsolver = new RBFGS(Prob, initx);
            RBFGSsolver->LineSearch_LS = static_cast<LSAlgoSM> (i);
            RBFGSsolver->Verbose = FINALRESULT;
            RBFGSsolver->Max_Iteration = RBFGSMI[i];
            //RBFGSsolver->CheckParams();
            RBFGSsolver->Run();
            if (RBFGSsolver->Getnormgfgf0() < tol)
                printf("SUCCESS!\n");
            else
                printf("FAIL!\n");
            delete RBFGSsolver;
        }
    }
    
    // test LRBFGS
    if(stringinclude(Methodnames, "LRBFGS"))
    {
        integer LRBFGSMI[5] = { 200, 200, 200, 200, 200 };
        printf("\n********************************Check all line search algorithm in LRBFGS in %s*************************************\n", probname);
        for (integer i = 0; i < LSSM_INPUTFUN; i++)//
        {
            Prob->SetNumGradHess(false);
            LRBFGS *LRBFGSsolver = new LRBFGS(Prob, initx);
            LRBFGSsolver->LineSearch_LS = static_cast<LSAlgoSM> (i);
            LRBFGSsolver->Verbose = FINALRESULT; //ITERRESULT;//
            LRBFGSsolver->Max_Iteration = LRBFGSMI[i];
            //LRBFGSsolver->CheckParams();
            LRBFGSsolver->Run();

            if (LRBFGSsolver->Getnormgfgf0() < tol)
                printf("SUCCESS!\n");
            else
                printf("FAIL!\n");

            delete LRBFGSsolver;
        }
        printf("\n********************************Check all line search algorithm in LRBFGS in %s with numerical gradient and Hessian*************************************\n", probname);
        for (integer i = 0; i < LSSM_INPUTFUN; i++)//
        {
            Prob->SetNumGradHess(true);
            LRBFGS *LRBFGSsolver = new LRBFGS(Prob, initx);
            LRBFGSsolver->LineSearch_LS = static_cast<LSAlgoSM> (i);
            LRBFGSsolver->Verbose = FINALRESULT; //ITERRESULT;//
            LRBFGSsolver->Max_Iteration = LRBFGSMI[i];
            //LRBFGSsolver->CheckParams();
            LRBFGSsolver->Run();

            if (LRBFGSsolver->Getnormgfgf0() < tol)
                printf("SUCCESS!\n");
            else
                printf("FAIL!\n");

            delete LRBFGSsolver;
        }
    }
    
    // test RTRSD
    if(stringinclude(Methodnames, "RTRSD"))
    {
        {
            printf("\n********************************Check RTRSD in %s*************************************\n", probname);
            Prob->SetNumGradHess(false);
            RTRSD RTRSDsolver(Prob, initx);
            RTRSDsolver.Verbose = FINALRESULT;
            RTRSDsolver.kappa = 0.001;
            RTRSDsolver.Max_Iteration = 5000;
            //RTRSDsolver.CheckParams();
            RTRSDsolver.Run();
            if (RTRSDsolver.Getnormgfgf0() < tol)
                printf("SUCCESS!\n");
            else
                printf("FAIL!\n");
        }
        {
            printf("\n********************************Check RTRSD in %s with numerical gradient and Hessian*************************************\n", probname);
            Prob->SetNumGradHess(true);
            RTRSD RTRSDsolver(Prob, initx);
            RTRSDsolver.Verbose = FINALRESULT;
            RTRSDsolver.kappa = 0.001;
            RTRSDsolver.Max_Iteration = 5000;
            //RTRSDsolver.CheckParams();
            RTRSDsolver.Run();
            if (RTRSDsolver.Getnormgfgf0() < tol)
                printf("SUCCESS!\n");
            else
                printf("FAIL!\n");
        }
    }
    
    // test RTRNewton
    if(stringinclude(Methodnames, "RTRNewton"))
    {
        {
            printf("\n********************************Check RTRNewton in %s*************************************\n", probname);
            Prob->SetNumGradHess(false);
            RTRNewton RTRNewtonsolver(Prob, initx);
            RTRNewtonsolver.Verbose = FINALRESULT;
            RTRNewtonsolver.Max_Iteration = 50;
            //RTRNewtonsolver.CheckParams();
            RTRNewtonsolver.Run();
            if (RTRNewtonsolver.Getnormgfgf0() < tol)
                printf("SUCCESS!\n");
            else
                printf("FAIL!\n");
            Prob->MinMaxEigValHess(RTRNewtonsolver.GetXopt()).Print("Min and Max of the eigenvalues of Hessian at x_k");
        }
        {
            printf("\n********************************Check RTRNewton in %s with numerical gradient and Hessian*************************************\n", probname);
            Prob->SetNumGradHess(true);
            RTRNewton RTRNewtonsolver(Prob, initx);
            RTRNewtonsolver.Verbose = FINALRESULT;
            RTRNewtonsolver.Max_Iteration = 50;
            //RTRNewtonsolver.CheckParams();
            RTRNewtonsolver.Run();
            if (RTRNewtonsolver.Getnormgfgf0() < tol)
                printf("SUCCESS!\n");
            else
                printf("FAIL!\n");
        }
    }
    
    // test RTRSR1
    if(stringinclude(Methodnames, "RTRSR1"))
    {
        {
            printf("\n********************************Check RTRSR1 in %s*************************************\n", probname);
            Prob->SetNumGradHess(false);
            RTRSR1 RTRSR1solver(Prob, initx);
            RTRSR1solver.Verbose = FINALRESULT;
            RTRSR1solver.Max_Iteration = 200;
            //RTRSR1solver.CheckParams();
            RTRSR1solver.Run();
            if (RTRSR1solver.Getnormgfgf0() < tol)
                printf("SUCCESS!\n");
            else
                printf("FAIL!\n");
        }
        {
            printf("\n********************************Check RTRSR1 in %s with numerical gradient and Hessian*************************************\n", probname);
            Prob->SetNumGradHess(true);
            RTRSR1 RTRSR1solver(Prob, initx);
            RTRSR1solver.Verbose = FINALRESULT;
            RTRSR1solver.Max_Iteration = 200;
            //RTRSR1solver.CheckParams();
            RTRSR1solver.Run();
            if (RTRSR1solver.Getnormgfgf0() < tol)
                printf("SUCCESS!\n");
            else
                printf("FAIL!\n");
        }
    }
    
    // test LRTRSR1
    if(stringinclude(Methodnames, "LRTRSR1"))
    {
        {
            printf("\n********************************Check LRTRSR1 in %s*************************************\n", probname);
            Prob->SetNumGradHess(false);
            LRTRSR1 LRTRSR1solver(Prob, initx);
            LRTRSR1solver.OutputGap = 1;
            LRTRSR1solver.Max_Iteration = 500;
            //LRTRSR1solver.Shrinked_tau = 0.1;
            //LRTRSR1solver.LengthSY = 4;
            LRTRSR1solver.Verbose = FINALRESULT;//-- FINALRESULT;
            //LRTRSR1solver.CheckParams();
            LRTRSR1solver.Run();
            if (LRTRSR1solver.Getnormgfgf0() < tol)
                printf("SUCCESS!\n");
            else
                printf("FAIL!\n");
        }
        {
            printf("\n********************************Check LRTRSR1 in %s with numerical gradient and Hessian*************************************\n", probname);
            Prob->SetNumGradHess(true);
            LRTRSR1 LRTRSR1solver(Prob, initx);
            LRTRSR1solver.OutputGap = 1;
            LRTRSR1solver.Max_Iteration = 500;
            //LRTRSR1solver.Shrinked_tau = 0.1;
            //LRTRSR1solver.LengthSY = 4;
            LRTRSR1solver.Verbose = FINALRESULT;//-- FINALRESULT;
            //LRTRSR1solver.CheckParams();
            LRTRSR1solver.Run();
            if (LRTRSR1solver.Getnormgfgf0() < tol)
                printf("SUCCESS!\n");
            else
                printf("FAIL!\n");
        }
    }
}

bool stringinclude(std::vector<std::string> names, std::string name)
{
    for(integer i = 0; i < names.size(); i++)
        if(strcmp(names[i].c_str(), name.c_str()) == 0)
            return true;
    return false;
}

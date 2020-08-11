
#include "Solvers/LRTRSR1.h"

/*Define the namespace*/
namespace ROPTLIB{

	LRTRSR1::LRTRSR1(const Problem *prob, const Variable *initialx)
	{
		Initialization(prob, initialx);
	};

	void LRTRSR1::SetProbX(const Problem *prob, const Variable *initialx)
	{
		SolversSMTR::SetProbX(prob, initialx);
		prob->SetUseGrad(true);
		prob->SetUseHess(false);
        s = Prob->GetDomain()->GetEMPTY();
        y = Prob->GetDomain()->GetEMPTY();
	};

	void LRTRSR1::SetDefaultParams(void)
	{
		SolversSMTR::SetDefaultParams();
		theta = static_cast<realdp> (0.1);
		kappa = static_cast<realdp> (0.1);
		LengthSY = 4;
        LMrestart = true;
		S = nullptr;
		Y = nullptr;
		YMGS = nullptr;
		inpss = 0;
		inpsy = 0;
		inpyy = 0;
		Currentlength = 0;
		beginidx = 0;
		SS = nullptr;
		SY = nullptr;
		gamma = 1;
		SolverName.assign("LRTRSR1");
	};

    void LRTRSR1::SetParams(PARAMSMAP params)
    {
        SolversSMTR::SetParams(params);
        PARAMSMAP::iterator iter;
        for (iter = params.begin(); iter != params.end(); iter++)
        {
            if (iter->first == static_cast<std::string> ("LengthSY"))
            {
                LengthSY = static_cast<integer> (iter->second);
            }
            else
            if (iter->first == static_cast<std::string> ("LMrestart"))
            {
                LMrestart = static_cast<integer> (iter->second);
            }
        }
    };

	LRTRSR1::~LRTRSR1(void)
	{
		DeleteVectors(S, LengthSY);
		DeleteVectors(Y, LengthSY);
		DeleteVectors(YMGS, LengthSY);
		if (SS != nullptr)
			delete[] SS;
        SS = nullptr;
		if (SY != nullptr)
			delete[] SY;
        SY = nullptr;
	};

	void LRTRSR1::Run(void)
	{
		DeleteVectors(S, LengthSY);
		NewVectors(S, LengthSY);
		DeleteVectors(Y, LengthSY);
		NewVectors(Y, LengthSY);
		DeleteVectors(YMGS, LengthSY);
		NewVectors(YMGS, LengthSY);
        
		if (SS != nullptr)
			delete[] SS;
		SS = new realdp[LengthSY * LengthSY];
		if (SY != nullptr)
			delete[] SY;
		SY = new realdp[LengthSY * LengthSY];
		SolversSMTR::Run();
	};

    void LRTRSR1::tCG_TR(void)
    {/*This method is only used when the Euclidean metric or an intrinsic approach is used.*/
        /*note that if the manifold is a complex, then we still view it as a real manifold with real reparameterizations
        therefore, we first set the types of vectors to be real and then recover the types of them at the end of this function*/
        bool gf1iscomplex = gf1.Getiscomplex(); gf1.Setiscomplex(false);
        
        integer n = gf1.Getlength(); /*number of rows*/
        integer k = Currentlength;    /*number of columns*/
        integer idx;
#ifdef SINGLE_PRECISION
        realdp tolLocalNewton = 1e-5;
#else
        realdp tolLocalNewton = 1e-10;
#endif
        Vector Psi(n, k);
        PMGQ = Vector (k, k);
        realdp *Psiptr = Psi.ObtainWriteEntireData();
        realdp *PMGQptr = PMGQ.ObtainWriteEntireData();
        
        /*Compute Psi = Y - gamma S, PMGS = SY - gamma SS */
        for (integer i = 0; i < k; i++)
        {
            idx = (i + beginidx) % LengthSY;
            Mani->ScalarVectorAddVector(x1, -gamma, S[idx], Y[idx], &YMGS[i]);
            copy_(&n, const_cast<realdp *> (YMGS[i].ObtainReadData()), &GLOBAL::IONE, Psiptr + i * n, &GLOBAL::IONE);
        }
        
        for (integer i = 0; i < k; i++)
        {
            for (integer j = 0; j < k; j++)
            {
                PMGQptr[i + j * Currentlength] = SY[i + j * LengthSY] - gamma * SS[i + j * LengthSY];
            }
        }
        
        /* QR decomposion:           Psi = Q RT^T
        eigenvalue decomposition: U D U^T = eig(RT^T * PMGQ^{-1} RT)
        P_parallel = Q * U
        Therefore: Psi * (PMGQ)^{-1} * Psi = Q * U * D * U^T * Q^T = P_parallel * D * P_parallel^T */
        Vector RinvPMGQRT, P_parallel(Psi), D, Psig(Psi.Getcol(), 1), g_parallel(Psi.Getcol(), 1);
        if(k > 0)
        {
            Psi.QRDecom();
            RinvPMGQRT = Psi.Field("_R") * (PMGQ % (Psi.Field("_R").GetTranspose()));
            RinvPMGQRT.EigenDecomSym();
            
            P_parallel.AlphaABaddBetaThis(1, Psi.Field("_Q"), GLOBAL::N, RinvPMGQRT.Field("_EigVec"), GLOBAL::N, 0); /* P_parallel = Psi.Field("_Q") * RinvPMGQRT.Field("_EigVec"); */
            D = RinvPMGQRT.Field("_EigVal");
            Vector gf1reshape(gf1); gf1reshape.Reshape(n);
            Psig.AlphaABaddBetaThis(1, Psi, GLOBAL::T, gf1reshape, GLOBAL::N, 0); /* Psig = Psi.GetTranspose() * gf1.GetReshape(n); */
            g_parallel.AlphaABaddBetaThis(1, P_parallel, GLOBAL::T, gf1reshape, GLOBAL::N, 0); /* g_parallel = P_parallel.GetTranspose() * gf1.GetReshape(n); */
        }
        
        Vector Lambda(k + 1);
        realdp *Lambdaptr = Lambda.ObtainWriteEntireData();
        const realdp *Dptr = D.ObtainReadData();
        for (integer i = 0; i < k; i++)
            Lambdaptr[i] = Dptr[i] + gamma;
        Lambdaptr[k] = gamma;
        
        for (integer i = 0; i < k + 1; i++)
            if (fabs(Lambdaptr[i]) <= tolLocalNewton)
                Lambdaptr[i] = 0;
        
        realdp lambda_min = Lambdaptr[0];
        for(integer i = 0; i < k + 1; i++)
            if(lambda_min > Lambdaptr[i])
                lambda_min = Lambdaptr[i];
        realdp a_kp2 = std::sqrt(Mani->Metric(x1, gf1, gf1) - ((k == 0) ? 0 : g_parallel.DotProduct(g_parallel)));
        
        if (a_kp2 * a_kp2 < tolLocalNewton)
            a_kp2 = 0;
        
        Vector a_j(k + 1);
        realdp *a_jptr = a_j.ObtainWriteEntireData();
        const realdp *g_parallelptr = g_parallel.ObtainReadData();
        
        for (integer i = 0; i < k; i++)
            a_jptr[i] = g_parallelptr[i];
        a_jptr[k] = a_kp2;
        realdp tmpv = 0;
        for (integer i = 0; i < k + 1; i++)
            tmpv += a_jptr[i] * a_jptr[i] / Lambdaptr[i] / Lambdaptr[i];
        tmpv = sqrt(tmpv);

        realdp *pStar = eta2.ObtainWriteEntireData();

        realdp sigmaStar = 0;
        tCGstatusSM = TRSM_MIN;

        if (lambda_min > 0 && tmpv <= Delta)
        {
            if(k == 0)
                ComputeSBySMW(gamma, Psig, Vector (), PMGQ, Psi);
            else
                ComputeSBySMW(gamma, Psig, Psi.Field("_R"), PMGQ, Psi); /*gf1, Psig, YMGS, PMGQ, RT*/
        }
        else
            if (lambda_min <= 0 && PhiBar_fg(-lambda_min, Delta, Lambda, a_j) > 0)
            {
                sigmaStar = -lambda_min;
                realdp *v = new realdp[k + 1];
                for (integer i = 0; i < k + 1; i++)
                {
                    if (fabs(Lambdaptr[i] + sigmaStar) > tolLocalNewton)
                    {
                        v[i] = a_jptr[i] / (Lambdaptr[i] + sigmaStar);
                    }
                    else
                    {
                        v[i] = 0;
                    }
                }
                const realdp *P_parallelptr = P_parallel.ObtainReadData();
                gemm_(GLOBAL::N, GLOBAL::N, &n, &GLOBAL::IONE, &k, &GLOBAL::DNONE, const_cast<realdp *> (P_parallelptr), &n, v, &k, &GLOBAL::DZERO, pStar, &n);
                delete[] v;
                if (fabs(gamma + sigmaStar) > tolLocalNewton)
                {
                    eta2 = (P_parallel * g_parallel - gf1) / (gamma + sigmaStar) + eta2;
                }

                if (lambda_min < 0)
                {
                    realdp alpha = sqrt(Delta * Delta - dot_(&n, pStar, &GLOBAL::IONE, pStar, &GLOBAL::IONE));

                    Vector zstar(n);
                    if (fabs(lambda_min - Lambdaptr[0]) < tolLocalNewton)
                    {
                        zstar = (alpha / std::sqrt(P_parallel.DotProduct(P_parallel))) * P_parallel;
                    }
                    else
                    {
                        realdp norm_umin = 1;
                        for (integer i = 0; i < k; i++)
                        {
                            Vector subM = P_parallel.GetSubmatrix(i, i, 0, k - 1);
                            zstar.AlphaABaddBetaThis(1, P_parallel, GLOBAL::N, subM, GLOBAL::T, 0); /*zstar = P_parallel * P_parallel.GetSubmatrix(i, i, 0, k - 1).GetTranspose();*/
                            realdp *zstarptr = zstar.ObtainWritePartialData();
                            zstarptr[i] += 1;
                            norm_umin = std::sqrt(zstar.DotProduct(zstar));
                            if (norm_umin > tolLocalNewton)
                                break;
                        }
                        zstar = (alpha / norm_umin) * zstar;
                    }
                    eta2 = zstar + eta2;
                }
                tCGstatusSM = TRSM_NEGCURVTURE;
            }
            else
            {
                if (lambda_min > 0)
                {
                    sigmaStar = LocalNewton(0, 10, tolLocalNewton, Lambda, a_j);
                }
                else
                {
                    realdp sigmaHat = 0;
                    for (integer i = 0; i < k + 1; i++)
                    {
                        if (sigmaHat < a_jptr[i] / Delta - Lambdaptr[i])
                            sigmaHat = a_jptr[i] / Delta - Lambdaptr[i];
                    }
                    if (sigmaHat > -lambda_min)
                        sigmaStar = LocalNewton(sigmaHat, 10, tolLocalNewton, Lambda, a_j);
                    else
                        sigmaStar = LocalNewton(-lambda_min, 10, tolLocalNewton, Lambda, a_j);
                }
                if(k == 0)
                    ComputeSBySMW(gamma + sigmaStar, Psig, Vector (), PMGQ, Psi);
                else
                    ComputeSBySMW(gamma + sigmaStar, Psig, Psi.Field("_R"), PMGQ, Psi); /*gf1, Psig, YMGS, PMGQ, RT*/
                tCGstatusSM = TRSM_EXCREGION;
            }
        Mani->Projection(x1, eta2, &eta2);
        /*recover the types of data*/
        gf1.Setiscomplex(gf1iscomplex);
        eta2.Setiscomplex(gf1iscomplex);
    };

	void LRTRSR1::ComputeSBySMW(realdp tauStar, const Vector &Psig, const Vector &R, const Vector &PMGQ, const Vector &Psi)
	{/*gamma, gf1, Psig, YMGS, PMGQ, R*/
		if (Currentlength != 0)
		{
            realdp eps = 100 * std::numeric_limits<realdp>::epsilon();
            Vector vw(PMGQ);
            vw.AlphaABaddBetaThis(tauStar, R, GLOBAL::T, R, GLOBAL::N, tauStar * tauStar);/*Vector vw = tauStar * (tauStar * PMGQ + (R.GetTranspose() * R));*/
            
            /*if vw is a rank deficient matrix, then eta2 = (static_cast<realdp> (-1) / tauStar) * gf1; */
            vw.HHRDecom();
            const realdp *HHRptr = vw.Field("_HHR").ObtainReadData();
            realdp minv = std::abs(HHRptr[0]), maxv = std::abs(HHRptr[0]);
            for(integer i = 0; i < vw.Getrow(); i++)
            {
                if(minv > std::abs(HHRptr[i + vw.Getrow() * i]))
                    minv = std::abs(HHRptr[i + vw.Getrow() * i]);
                if(maxv < std::abs(HHRptr[i + vw.Getrow() * i]))
                    maxv = std::abs(HHRptr[i + vw.Getrow() * i]);
            }
            
            if(minv / maxv < eps || maxv == static_cast<realdp> (0)) /*rank deficient*/
            {
                eta2 = gf1; eta2.ScalarTimesThis(static_cast<realdp> (-1) / tauStar); /* eta2 = (static_cast<realdp> (-1) / tauStar) * gf1; */
                return;
            }

            /*otherwise, use below eta2*/
            /* eta2 = Psi * (vw % Psig) + (static_cast<realdp> (-1) / tauStar) * gf1; */
            Vector tmp = vw % Psig;
            eta2 = gf1;
            eta2.AlphaABaddBetaThis(1, Psi, GLOBAL::N, tmp, GLOBAL::N, static_cast<realdp> (-1) / tauStar);
            
			return;
		}
        /*eta2 = gf1 * (static_cast<realdp> (-1) / tauStar);*/
        eta2 = gf1;
        eta2.ScalarTimesThis(static_cast<realdp> (-1) / tauStar);
	};

	realdp LRTRSR1::PhiBar_fg(realdp nlambda_min, realdp Delta, const Vector &Lambda, const Vector &a_j, realdp *gf)
	{
		integer k = Currentlength;    /*number of columns*/
        Vector tmp = Lambda + nlambda_min;
        
		realdp eps_tol = static_cast<realdp> (1e-10);
		realdp gradf = 0;
        const realdp *a_jptr = a_j.ObtainReadData();
        realdp *tmpptr = tmp.ObtainWritePartialData();

		bool flag = false;
		for (integer i = 0; i < k + 1; i++)
		{
			if (fabs(a_jptr[i]) < eps_tol || fabs(tmpptr[i]) < eps_tol)
			{
				flag = true;
				break;
			}
		}
        
		if (flag)
		{
			realdp pnorm2 = 0;
			for (integer i = 0; i < k + 1; i++)
			{
				if (fabs(a_jptr[i]) > eps_tol && fabs(tmpptr[i]) < eps_tol)
				{
					gradf = static_cast<realdp> (1) / eps_tol;
                    *gf = gradf;
					return static_cast<realdp> (-1) / Delta;
				}
				else
				if (fabs(a_jptr[i]) > eps_tol && fabs(tmpptr[i]) > eps_tol)
				{
					pnorm2 += a_jptr[i] * a_jptr[i] / tmpptr[i] / tmpptr[i];
					gradf += a_jptr[i] * a_jptr[i] / tmpptr[i] / tmpptr[i] / tmpptr[i];
				}
			}
			realdp norm = sqrt(pnorm2);
			gradf = gradf / norm / norm / norm;
            *gf = gradf;
			return static_cast<realdp> (1) / norm - static_cast<realdp> (1) / Delta;
		}
        
		for (integer i = 0; i < k + 1; i++)
			gradf += a_jptr[i] * a_jptr[i] / tmpptr[i] / tmpptr[i] / tmpptr[i];
        
		for (integer i = 0; i < k + 1; i++)
			tmpptr[i] = a_jptr[i] / tmpptr[i];
        
        realdp norm = tmp.Fnorm();
        
		gradf = gradf / norm / norm / norm;
        *gf = gradf;
		return static_cast<realdp> (1) / norm - static_cast<realdp> (1) / Delta;
	};

    realdp LRTRSR1::PhiBar_fg(realdp nlambda_min, realdp Delta, const Vector &Lambda, const Vector &a_j)
    {
        integer k = Currentlength;    /*number of columns*/
        Vector tmp = Lambda + nlambda_min;
        
        realdp eps_tol = static_cast<realdp> (1e-10);
        const realdp *a_jptr = a_j.ObtainReadData();
        realdp *tmpptr = tmp.ObtainWritePartialData();

        bool flag = false;
        for (integer i = 0; i < k + 1; i++)
        {
            if (fabs(a_jptr[i]) < eps_tol || fabs(tmpptr[i]) < eps_tol)
            {
                flag = true;
                break;
            }
        }
        
        if (flag)
        {
            realdp pnorm2 = 0;
            for (integer i = 0; i < k + 1; i++)
            {
                if (fabs(a_jptr[i]) > eps_tol && fabs(tmpptr[i]) < eps_tol)
                {
                    return static_cast<realdp> (-1) / Delta;
                }
                else
                if (fabs(a_jptr[i]) > eps_tol && fabs(tmpptr[i]) > eps_tol)
                {
                    pnorm2 += a_jptr[i] * a_jptr[i] / tmpptr[i] / tmpptr[i];
                }
            }
            realdp norm = sqrt(pnorm2);
            return static_cast<realdp> (1) / norm - static_cast<realdp> (1) / Delta;
        }
        for (integer i = 0; i < k + 1; i++)
            tmpptr[i] = a_jptr[i] / tmpptr[i];
        realdp norm = tmp.Fnorm();
        return static_cast<realdp> (1) / norm - static_cast<realdp> (1) / Delta;
    };

	realdp LRTRSR1::LocalNewton(realdp x0, integer maxIter, realdp tol, const Vector &Lambda, const Vector &a_j)
	{
		integer k = 0;
		realdp gf = 0;
		realdp fv = PhiBar_fg(x0, Delta, Lambda, a_j, &gf);
		while (fabs(fv) > std::numeric_limits<realdp>::epsilon() && k < maxIter)
		{
			x0 -= fv / gf;
			fv = PhiBar_fg(x0, Delta, Lambda, a_j, &gf);
			k++;
		}
		return x0;
	};

	void LRTRSR1::CheckParams(void)
	{
		SolversSMTR::CheckParams();
		char YES[] = "YES";
		char NO[] = "NO";
		char *status;

		printf("LRTRSR1 METHOD PARAMETERS:\n");
        status = (LengthSY >= 0) ? YES : NO;
        printf("LengthSY      :%15d[%s],\t", LengthSY, status);
        status = YES;
        printf("LMrestart     :%15d[%s],\n", LMrestart, status);
	};

	Vector &LRTRSR1::HessianEta(const Vector &Eta, Vector *result)
	{
		return HvLRTRSR1(Eta, result);
	};

	void LRTRSR1::UpdateData(void)
	{
		UpdateDataLRTRSR1();
	};

	void LRTRSR1::Acceptence(void)
	{
		for (integer i = 0; i < Currentlength; i++)
		{
			Mani->VectorTransport(x1, eta2, x2, S[i], &S[i]);
			Mani->VectorTransport(x1, eta2, x2, Y[i], &Y[i]);
		}
	};

	void LRTRSR1::PrintInfo(void)
	{
        printf("i:%d,f:%.3e,df/f:%.3e,", iter, f2, ((f1 - f2) / std::fabs(f2)));

        printf("|gf|:%.3e,time:%.2g,", ngf2, static_cast<realdp>(getTickCount() - starttime) / CLK_PS);

        printf("rho:%.2e,radius:%.3e,tCGstatus:%s,innerIter:%d,", rho, Delta, tCGstatusSetSMnames[tCGstatusSM].c_str(), innerIter);
        
        printf("gamma:%.3e,inpss:%.3e,inpsy:%.3e,inpyy:%.3e,IsUpdateHessian:%d,", gamma, inpss, inpsy, inpyy, isupdated);
        
        printf("nf:%d,ng:%d,", nf, ng);
        
        if (nH != 0)
            printf("nH:%d,", nH);
        
        printf("nR:%d,", nR);
        
        if (nV != 0)
            printf("nV(nVp):%d(%d),", nV, nVp);
        
        printf("\n");
	};

    Vector &LRTRSR1::HvLRTRSR1(const Vector &Eta, Vector *result)
    {
        /* This function makes use of SS, SY and gamma to evaluate the action of Hessian approximation [HAG2014, (64)].
        [HAG2014]: W. Huang, P.-A. Absil, and K. A. Gallivan. A Riemannian symmetric rank-one trustregion method.
        Mathematical Programming, 150(2):179?16, February 2015.

        SS is the Q in (46), SY is the P in (46), PMGQ is the P - gamma Q in (46).
        */
        Vector v(Currentlength), v2;
        realdp *vptr = v.ObtainWriteEntireData();

        for (integer i = 0; i < Currentlength; i++)
            vptr[i] = Mani->Metric(x1, YMGS[i], Eta);
        
        if (Currentlength > 0)
        {
            v2 = PMGQ % v;
        }
        const realdp *v2ptr = v2.ObtainReadData();

        Mani->ScalarTimesVector(x1, gamma, Eta, result);
        for (integer i = 0; i < Currentlength; i++)
        {
            Mani->ScalarVectorAddVector(x1, v2ptr[i], YMGS[i], *result, result);
        }
        return *result;
    };

    void LRTRSR1::UpdateDataLRTRSR1(void)
    {
        realdp denorminator, norm2ymBs;
        realdp mintolsq = std::numeric_limits<realdp>::epsilon();
        realdp mintol = sqrt(mintolsq);
        Prob->Grad(x2, &gf2); ng++;
        s = eta2;
        Mani->InverseVectorTransport(x1, eta2, x2, gf2, &eta1); nV++;
        Mani->VectorLinearCombination(x1, 1, eta1, -1, gf1, &y);
        
        if (iter == 0) /* This is for the robustness when the cost function is quadratic and its Hessian is identity everywhere.*/
        {
            inpsy = Mani->Metric(x1, s, y);
            inpyy = Mani->Metric(x1, y, y);
            gamma = inpyy / inpsy;
            Mani->ScalarTimesVector(x1, gamma, s, &Heta2);
        }
        Vector zeta(y); Mani->ScalarVectorAddVector(x1, -1, Heta2, y, &zeta);
        denorminator = Mani->Metric(x1, s, zeta);
        norm2ymBs = Mani->Metric(x1, zeta, zeta);
        inpss = Mani->Metric(x1, s, s);
        inpsy = Mani->Metric(x1, s, y);
        
        if (denorminator * denorminator >= mintolsq * inpss * norm2ymBs && (norm2ymBs >= mintolsq || ngf2 / ngf0 < 1e-3)
            && (iter != 0 || fabs(gamma - inpsy / inpss) > mintol)) /* This is for the robustness when the cost
             function is quadratic and its Hessian is identity everywhere. */
        {
            inpyy = Mani->Metric(x1, y, y);
            if(isupdated == true)
                gamma = inpyy / inpsy;
            
            if(LMrestart && Currentlength >= LengthSY)
                Currentlength = 0;

            /*if s and y are accepted, then S and Y need to be updated. It follows that the matrices SY and SS need to be update.*/
            if (Currentlength < LengthSY)
            {
                S[Currentlength] = s;
                Y[Currentlength] = y;
                SS[Currentlength + Currentlength * LengthSY] = Mani->Metric(x1, S[Currentlength], S[Currentlength]);
                SY[Currentlength + Currentlength * LengthSY] = Mani->Metric(x1, S[Currentlength], Y[Currentlength]);
                for (integer i = 0; i < Currentlength; i++)
                {
                    SS[Currentlength + i * LengthSY] = Mani->Metric(x1, S[Currentlength], S[i]);
                    SS[i + Currentlength * LengthSY] = SS[Currentlength + i * LengthSY];
                    SY[Currentlength + i * LengthSY] = Mani->Metric(x1, S[Currentlength], Y[i]);
                    SY[i + Currentlength * LengthSY] = SY[Currentlength + i * LengthSY];
                }
                Currentlength++;
            }
            else
            {
                if(LengthSY > 0)
                {
                    S[beginidx] = s;
                    Y[beginidx] = y;
                    for (integer i = 0; i < LengthSY - 1; i++)
                    {
                        for (integer j = 0; j < LengthSY - 1; j++)
                        {
                            SS[i + j * LengthSY] = SS[i + 1 + (j + 1) * LengthSY];
                            SY[i + j * LengthSY] = SY[i + 1 + (j + 1) * LengthSY];
                        }
                    }
                    SS[LengthSY * LengthSY - 1] = Mani->Metric(x1, S[beginidx], S[beginidx]);
                    SY[LengthSY * LengthSY - 1] = Mani->Metric(x1, S[beginidx], Y[beginidx]);
                    integer idx = 0;
                    for (integer i = 0; i < LengthSY - 1; i++)
                    {
                        idx = (i + beginidx + 1) % LengthSY;
                        SS[i + (LengthSY - 1) * LengthSY] = Mani->Metric(x1, S[idx], S[beginidx]);
                        SS[LengthSY - 1 + i * LengthSY] = SS[i + (LengthSY - 1) * LengthSY];
                        SY[i + (LengthSY - 1) * LengthSY] = Mani->Metric(x1, Y[idx], S[beginidx]);
                        SY[LengthSY - 1 + i * LengthSY] = SY[i + (LengthSY - 1) * LengthSY];
                    }
                    beginidx = (++beginidx) % LengthSY;
                }
            }

            isupdated = true;
        }
        else
        {
            isupdated = false;
        }
    };

}; /*end of ROPTLIB namespace*/

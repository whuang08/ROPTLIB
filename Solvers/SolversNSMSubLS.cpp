
#include "Solvers/SolversNSMSubLS.h"

/*Define the namespace*/
namespace ROPTLIB{
	void SolversNSMSubLS::Run(void)
	{
        Variable xTemp(x1);
        Vector gfTemp = Prob->GetDomain()->GetEMPTY();
		SolversNSMSub::Run();

		LSstatus = LSNSM_SUCCESS;
		f1 = Prob->f(x1); nf++;
		f2 = f1 + 1;
        Prob->Grad(x1, &gf1); ng++;
		ndir0 = sqrt(Mani->Metric(x1, gf1, gf1));
		newslope = 0;
		iter = 0;
        Vector initialstepsize(1 + Max_Iteration), acceptedstepsize(1 + Max_Iteration);
        realdp *timeSeriesptr = timeSeries.ObtainWritePartialData();
        realdp *funSeriesptr = funSeries.ObtainWritePartialData();
        realdp *dirSeriesptr = dirSeries.ObtainWritePartialData();
        realdp *initialstepsizeptr = initialstepsize.ObtainWritePartialData();
        realdp *acceptedstepsizeptr = acceptedstepsize.ObtainWritePartialData();
        
        if (Verbose >= FINALRESULT)
            printf("i:%d,f:%.3e,\n", iter, f1);
		if (Verbose >= ITERRESULT)
		{
			timeSeriesptr[iter] = static_cast<realdp > (getTickCount() - starttime) / CLK_PS;
			funSeriesptr[iter] = f1;
			dirSeriesptr[iter] = 0;
            initialstepsizeptr[iter] = 0;
            acceptedstepsizeptr[iter] = 0;
		}
		bool isstop = false;
		realdp ftmp = 0;
        LSstatus = LSNSM_SUCCESS;

		/*Start the loop*/
		while ( iter < Max_Iteration && ((!isstop) || iter < Min_Iteration || Eps > Min_Eps) && (LSstatus == LSNSM_SUCCESS || LSstatus == LSNSM_MAXSTEPSIZE) )
		{
			GetSearchDir(); /* Obtain search direction eta1, the minimum length vector minPv, and ndir1 = \|minPv\|_P */

			/*Call the function to check whether the stopping criterion is satisfied or not.
			The default function is written in Solvers.h and Solvers.cpp*/
			isstop = IsStopped();
			if (fabs(ndir1) <= Del || fabs(ndir1) <= Tolerance)
			{
				Eps *= Theta_eps;
				Eps = (Eps > Min_Eps) ? Eps : Min_Eps;
				if (Eps > Min_Eps)
				{
					Del *= Theta_del;
				}
				else
				{
					while (fabs(ndir1) < Del)
					{
						Del *= Theta_del;
					}
				}
				if (Verbose >= ITERRESULT)
					printf("Shinking Epsilon and Delta to %g and %g respectively, ndir1:%g.\n", Eps, Del, ndir1);
				continue;
			}

			initialslope = Mani->Metric(x1, minPv, eta1);
			/*Compute initial step size for the next iteration*/
            stepsize = 1; /*Initial step size is fixed to be one*/
            initiallength = stepsize;

			/* Call the specified linesearch algorithm.
			Note that in the linesearch algorithm, we need to obtain
			accepted stepsize, eta2=stepsize*eta1, x2 = R_{x_1}(eta_2), f2 = f(x2), and gf2 = grad f(x_2) */
            LinesearchWolfeLipschitz();
            
            /*Output debug information if necessary.*/
            if (LSstatus < LSNSM_SUCCESS && Verbose > FINALRESULT)
            {
                printf("Linesearch fails! LSstatus:%s\n", LSstatusSetnames[LSstatus].c_str());
            }

			iter++;

			/*Update some data, such as Hessian approximation in quasi-Newton methods, etc*/
			UpdateData();

			if (Verbose >= ITERRESULT)
			{
				/*Output information*/
				if (iter % OutputGap == 0)
				{
					PrintInfo();
				}
				/*Store debug information in the arrays*/
				timeSeriesptr[iter] = static_cast<realdp>(getTickCount() - starttime) / CLK_PS;
				funSeriesptr[iter] = f2; dirSeriesptr[iter - 1] = ndir1;
                initialstepsizeptr[iter] = initiallength;
                acceptedstepsizeptr[iter] = stepsize;
			}

			/*Switch information at x1 and x2*/
			xTemp = x1; x1 = x2; x2 = xTemp;
			gfTemp = gf1; gf1 = gf2; gf2 = gfTemp;
			ftmp = f1;
			f1 = f2;
			f2 = ftmp;
		}
		ComTime = static_cast<realdp>(getTickCount() - starttime) / CLK_PS;
		if (Verbose >= ITERRESULT)
			lengthSeries = iter + 1;
        PrintFinalInfo();
	};

    void SolversNSMSubLS::PrintInfo(void)
    {
        printf("i:%d,f:%.3e,df/f:%.3e,", iter, f2, ((f1 - f2) / std::fabs(f2)));

        printf("|nd|:%.3e,t0:%.2e,t:%.2e,s0:%.2e,s:%.2e,time:%.2e,", ndir1, initiallength, stepsize, initialslope, newslope,  static_cast<realdp>(getTickCount() - starttime) / CLK_PS);

        if (subprobtimes != 0)
            printf("nsubprob:%d,", subprobtimes);

        printf("nf:%d,ng:%d,", nf, ng);
        
        if (nH != 0)
            printf("nH:%d,", nH);
        
        printf("nR:%d,", nR);
        
        if (nV != 0)
            printf("nV(nVp):%d(%d),", nV, nVp);

        printf("Eps:%.3e,", Eps);
        printf("Del:%.3e,", Del);
        printf("\n");
    };

    void SolversNSMSubLS::PrintFinalInfo(void)
    {
        if (Verbose >= FINALRESULT)
        {
            printf("Iter:%d,f:%.3e,ndir:%.3e,|nd|/|nd0|:%.3e,time:%.2e,nsubprob:%d,nf:%d,ng:%d,nR:%d,", iter, f1, ndir1, ndir1/ndir0, ComTime, subprobtimes, nf, ng, nR);
            if (nH != 0)
            {
                printf("nH:%d,", nH);
            }
            if (nV != 0)
            {
                printf("nV(nVp):%d(%d),", nV, nVp);
            }
            printf("Eps:%.3e,", Eps);
            printf("Del:%.3e,", Del);
            printf("\n");
        }
    };

	void SolversNSMSubLS::SetProbX(const Problem *prob, const Variable *initialx)
	{
		SolversNSMSub::SetProbX(prob, initialx);
		prob->SetUseGrad(true);
		prob->SetUseHess(false);
        minPv = Prob->GetDomain()->GetEMPTY();
        eta1 = Prob->GetDomain()->GetEMPTY();
        eta2 = Prob->GetDomain()->GetEMPTY();
	};

	void SolversNSMSubLS::SetDefaultParams(void)
	{
		SolversNSMSub::SetDefaultParams();
		Theta_eps = static_cast<realdp> (0.01);
		Theta_del = static_cast<realdp> (0.01);
		Eps = 1;
		Min_Eps = static_cast<realdp> (1e-6);
		Del = 1;
        LS_alpha = 0.0001;
        LS_beta = 0.999;
        Minstepsize = 2e-16;
        Maxstepsize = 1000;
        lambdaLower = 1e-7;
        lambdaUpper = 1e7;

        LSstatusSetnames = new std::string[LSSTATUSSETNSMLENGTH];
        LSstatusSetnames[LSNSM_NOCURVATURE].assign("LSNSM_NOCURVATURE");
        LSstatusSetnames[LSNSM_MINSTEPSIZE].assign("LSNSM_MINSTEPSIZE");
        LSstatusSetnames[LSNSM_MAXSTEPSIZE].assign("LSNSM_MAXSTEPSIZE");
        LSstatusSetnames[LSNSM_LSERROR].assign("LSNSM_LSERROR");
        LSstatusSetnames[LSNSM_SUCCESS].assign("LSNSM_SUCCESS");
	};

	SolversNSMSubLS::~SolversNSMSubLS(void)
	{
		delete[] LSstatusSetnames;
	};

	void SolversNSMSubLS::CheckParams(void)
	{
		SolversNSMSub::CheckParams();
		char YES[] = "YES";
		char NO[] = "NO";
		char *status;

		printf("LINE SEARCH SUBGRADIENT-BASED METHODs PARAMETERS:\n");
		status = (Eps > 0) ? YES : NO;
		printf("Eps           :%15g[%s],\t", Eps, status);
		status = (Theta_eps > 0 && Theta_eps < 1) ? YES : NO;
		printf("Theta_eps     :%15g[%s],\n", Theta_eps, status);
		status = (Min_Eps > 0 && Min_Eps < 1) ? YES : NO;
		printf("Min_Eps       :%15g[%s],\t", Min_Eps, status);
		status = (Del > 0) ? YES : NO;
		printf("Del           :%15g[%s],\n", Del, status);
        status = (Theta_del > 0 && Theta_del < 1) ? YES : NO;
        printf("Theta_del     :%15g[%s],\t", Theta_del, status);
        status = (LS_alpha > 0 && LS_alpha < 1) ? YES : NO;
        printf("LS_alpha      :%15g[%s],\n", LS_alpha, status);
        status = (LS_beta > 0 && LS_beta < 1) ? YES : NO;
        printf("LS_beta       :%15g[%s],\t", LS_beta, status);
        status = (Minstepsize > 0 && Minstepsize <= 1) ? YES : NO;
        printf("Minstepsize   :%15g[%s],\n", Minstepsize, status);
        status = (Maxstepsize >= 1) ? YES : NO;
        printf("Maxstepsize   :%15g[%s],\t", Maxstepsize, status);
        status = (lambdaLower > 0 && lambdaLower < 1) ? YES : NO;
        printf("lambdaLower   :%15g[%s],\n", lambdaLower, status);
        status = (lambdaUpper > 1) ? YES : NO;
        printf("lambdaUpper   :%15g[%s],\n", lambdaUpper, status);
	};

	void SolversNSMSubLS::GetSearchDir(void)
	{
        realdp nminPv;
        gfs[0] = gf1;
		Currentlengthgfs = 1;
		realdp Pngfsq = 0, hb;
		for (integer i = 0; i < Lengthgfs - 1; i++)
		{
			/*minPv = argmin \| v \|_P^2, v in convex hall of W, let eta1 = - P minPv.*/
			Pngfsq = MinPNormConHull(Mani, x1, gfs, Currentlengthgfs, this, Hv, minPv);
			(this->*Hv)(minPv, &eta1);
            Mani->ScalarTimesVector(x1, -1.0, eta1, &eta1);
			if (Currentlengthgfs > 1)
				subprobtimes++;

            /*if nminPv is sufficient small, then stop*/
			if (sqrt(Pngfsq) > Del * lambdaLower)
				nminPv = sqrt(Mani->Metric(x1, minPv, minPv));
			else
			{
				nminPv = sqrt(MinPNormConHull(Mani, x1, gfs, Currentlengthgfs, nullptr, nullptr, minPv));
				subprobtimes++;
			}
			if (nminPv <= Del)
			{
				return;
			}
            /*end the stopping criterion*/
            
			ndir1 = sqrt(Mani->Metric(x1, eta1, eta1));
			stepsize = Eps / ndir1;
			f2 = h();
			hb = f2 - f1 + LS_alpha * stepsize * Pngfsq; /* - Pngfsq is the initial slope*/
			if (hb <= 0)
			{
				return;
			}
			/*get gf2 and stepsize*/
			Increasing(ndir1, Pngfsq, hb);
            Vector invTgf2(gf2);
            Mani->InverseVectorTransport(x1, eta2, x2, gf2, &invTgf2);
			realdp betay = Mani->Beta(x1, eta2);
			Mani->ScalarTimesVector(x1, static_cast<realdp> (1) / betay, invTgf2, &invTgf2);
            gfs[Currentlengthgfs] = invTgf2;
			Currentlengthgfs++;
			if (Currentlengthgfs == Lengthgfs && Verbose >= ITERRESULT)
			{
				printf("Warning: the number of W reaches its upper-bound: %d!\n", Currentlengthgfs);
			}
		}
	};

	void SolversNSMSubLS::Increasing(realdp neta1, realdp Pngfsq, realdp hb)
	{
		realdp a = 0, b = stepsize, ht;
		integer times = 0;
		while (times < 10)
		{
			if (dh() + LS_alpha * Pngfsq < 0)
			{
				stepsize = (a + b) / 2;
				f2 = h();
				ht = f2 - f1 + LS_alpha * Eps / neta1 * Pngfsq;
				if (hb > ht)
				{
					a = stepsize;
				}
				else
				{
					b = stepsize;
				}
			}
			else
			{
				break;
			}
			times++;
		}
		if (times == 10 && Verbose >= ITERRESULT)
		{
			printf("warning: the loop in SolversNSMSubLS::Increasing reaches the upperbound!\n");
		}
	};

    Vector &SolversNSMSubLS::HvSub(const Vector &v, Vector *result)
    {
        *result = v;
        return *result;
    };

    void SolversNSMSubLS::LinesearchWolfeLipschitz(void)
    {
        LSstatus = LSNSM_SUCCESS;
        realdp prestepsize = 0;
        while (1)
        {
            f2 = h();
            if (f2 - f1 - LS_alpha * initialslope * stepsize > 0) /*-initialslope = \|g\|_P^2*/
            {
                ZoomLP(prestepsize, stepsize);
                break;
            }
            if (dh() >= LS_beta * initialslope)
            {
                break;
            } else
            {
                prestepsize = stepsize;
                stepsize = (2 * stepsize < Maxstepsize) ? 2 * stepsize : Maxstepsize;

                if (stepsize == Maxstepsize)
                    break;
            }
        }
        if (stepsize <= Minstepsize)
        {
            LSstatus = LSNSM_MINSTEPSIZE;
        }
        if (stepsize >= Maxstepsize)
        {
            LSstatus = LSNSM_MAXSTEPSIZE;
        }
    };

    void SolversNSMSubLS::ZoomLP(realdp a, realdp b)
    {
        integer times = 0;
        while (1)
        {
            stepsize = (a + b) / 2;

            if (stepsize < Minstepsize)
            {
                if (Verbose > FINALRESULT)
                {
                    printf("Warning: step size reaches the minimum: %3.2e!\n", Minstepsize);
                }
                LSstatus = LSNSM_MINSTEPSIZE;
                dh();
                break;
            }

            f2 = h();
            if (f2 - f1 - LS_alpha * initialslope * stepsize > 0)
            {
                b = stepsize;
            }
            else
            {
                if (dh() >= LS_beta * initialslope)
                {
                    break;
                }
                else
                {
                    a = stepsize;
                    if (times >= 30)
                    {
                        if (Verbose > FINALRESULT)
                        {
                            printf("warning: line search algorithm reaches upper bound iterations when finding curvature condition.\n");
                        }
                        LSstatus = LSNSM_NOCURVATURE;
                        break;
                    }
                }
            }
            times++;
        };
    };

    realdp SolversNSMSubLS::h(void)
    {
        Mani->ScalarTimesVector(x1, stepsize, eta1, &eta2);
        Mani->Retraction(x1, eta2, &x2); nR++;
        nf++;
        return Prob->f(x2);
    };

    realdp SolversNSMSubLS::dh(void)
    {
        Prob->Grad(x2, &gf2); ng++;
        nV++;
        Vector diffeta1(eta1); Mani->DiffRetraction(x1, eta2, x2, eta1, &diffeta1, true);
        newslope = Mani->Metric(x2, gf2, diffeta1);
        return newslope;
    };

	void SolversNSMSubLS::SetParams(PARAMSMAP params)
	{
		SolversNSMSub::SetParams(params);
		PARAMSMAP::iterator iter;
		for (iter = params.begin(); iter != params.end(); iter++)
		{
			if (iter->first == static_cast<std::string> ("Eps"))
			{
				Eps = static_cast<realdp> (iter->second);
			}
			else
			if (iter->first == static_cast<std::string> ("Theta_eps"))
			{
				Theta_eps = static_cast<realdp> (iter->second);
			}
			else
			if (iter->first == static_cast<std::string> ("Min_Eps"))
			{
				Min_Eps = static_cast<realdp> (iter->second);
			}
			else
			if (iter->first == static_cast<std::string> ("Del"))
			{
				Del = static_cast<realdp> (iter->second);
			}
            else
            if (iter->first == static_cast<std::string> ("Theta_del"))
            {
                Theta_del = static_cast<realdp> (iter->second);
            }
            else
            if (iter->first == static_cast<std::string> ("LS_alpha"))
            {
                LS_alpha = static_cast<realdp> (iter->second);
            }
            else
            if (iter->first == static_cast<std::string> ("LS_beta"))
            {
                LS_beta = static_cast<realdp> (iter->second);
            }
            else
            if (iter->first == static_cast<std::string> ("Minstepsize"))
            {
                Minstepsize = static_cast<realdp> (iter->second);
            }
            else
            if (iter->first == static_cast<std::string> ("Maxstepsize"))
            {
                Maxstepsize = static_cast<realdp> (iter->second);
            }
            else
            if (iter->first == static_cast<std::string> ("lambdaLower"))
            {
                lambdaLower = static_cast<realdp> (iter->second);
            }
            if (iter->first == static_cast<std::string> ("lambdaUpper"))
            {
                lambdaUpper = static_cast<realdp> (iter->second);
            }
		}
	};
}; /*end of ROPTLIB namespace*/


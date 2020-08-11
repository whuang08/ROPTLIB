
#include "Solvers/SolversSMTR.h"

/*Define the namespace*/
namespace ROPTLIB{

	void SolversSMTR::Run(void)
	{
		Variable xTemp(x1);
		Vector gfTemp = Prob->GetDomain()->GetEMPTY();
		SolversSM::Run();
		starttime = getTickCount();
		realdp sqeps = sqrt(std::numeric_limits<realdp>::epsilon());

		f1 = Prob->f(x1); nf++;
		f2 = f1;
		Prob->Grad(x1, &gf1); ng++;

		ngf0 = sqrt(Mani->Metric(x1, gf1, gf1));
		ngf1 = ngf0; ngf2 = ngf1;

		iter = 0;
        realdp *timeSeriesptr = timeSeries.ObtainWritePartialData();
        realdp *funSeriesptr = funSeries.ObtainWritePartialData();
        realdp *gradSeriesptr = gradSeries.ObtainWritePartialData();
        if (Verbose >= FINALRESULT)
            printf("i:%d,f:%.3e,|gf|:%.3e,\n", iter, f1, ngf1);
		if (Verbose >= ITERRESULT)
		{
			timeSeriesptr[iter] = static_cast<realdp>(getTickCount() - starttime) / CLK_PS;
			funSeriesptr[iter] = f1;
            gradSeriesptr[iter] = ngf1;
		}
		Delta = initial_Delta;
		bool isstop = IsStopped();
		while (((!isstop) && iter < Max_Iteration) || iter < Min_Iteration)
		{
			InitialVector(); /* Obtain initial guess, eta1, for local model */
			tCG_TR(); /* obtain eta2 */
			Mani->Retraction(x1, eta2, &x2);	nR++;
			f2 = Prob->f(x2); nf++;
			if (std::isnan(f2) || std::isinf(f2)) /*Stop when got a nan or inf*/
			{
				printf("New function value is either nan or inf. Stop!\n");
				break;
			}
            HessianEta(eta2, &Heta2); nH++;
            Mani->ScalarVectorAddVector(x1, 0.5, Heta2, gf1, &eta1);
            rho = (f1 - f2) / (-Mani->Metric(x1, eta2, eta1));
			UpdateData(); /* Update S&Y or H or B */

			if (rho > 0.75)
			{
				if (tCGstatusSM == TRSM_EXCREGION || tCGstatusSM == TRSM_NEGCURVTURE || tCGstatusSM == TRSM_MAXITER)
					Delta *= Magnified_tau;
				if (Delta > maximum_Delta)
				{
					if (Verbose > ITERRESULT)
					{
						printf("reach the maximum of radius\n");
					}
					Delta = maximum_Delta;
				}
			}
			else
				if (rho < 0.25 && (ngf1 / (ngf0 + Tolerance) >= Accuracy)) /*If accurate enough, then always not shrink the region.*/
				{
					Delta *= Shrinked_tau;
					if (Delta < minimum_Delta)
					{
						if (Verbose > FINALRESULT)
						{
							printf("reach the minimum of radius. Stop!\n");
						}
						break;
					}
				}
			
			if (rho > Acceptence_Rho ||
				(fabs(f1 - f2) / (fabs(f1) + 1) < sqeps && f2 < f1) ||
				(ngf1 / (ngf0 + Tolerance) < Accuracy && (f2 < f1 || (tCGstatusSM == TRSM_MIN || tCGstatusSM == TRSM_LCON || tCGstatusSM == TRSM_SCON))) /*If accurate enough and the solution is in the trust region, then always accept the new iterate.*/
				)
			{
				Acceptence(); /* Algorithm specific operations */
				ngf2 = sqrt(Mani->Metric(x2, gf2, gf2));
				isstop = IsStopped(); /*This is done when the candidate is accepted. This is necessary for partly smooth stopping criterion*/
                xTemp = x1;
                x1 = x2;
                x2 = xTemp; xTemp.Delete();
                gfTemp = gf1; gf1 = gf2; gf2 = gfTemp;
				iter++;
				if (Verbose >= ITERRESULT && iter % OutputGap == 0)
				{
					PrintInfo(); /* Output information specific to Algorithms */
				}
				f1 = f2; ngf1 = ngf2;
			}
			else
			{
				iter++;
				if (Verbose >= ITERRESULT && iter % OutputGap == 0)
				{
					printf("X_{%d} WAS REJECTED.\n", iter);
					PrintInfo();
				}
			}

			if (Verbose >= ITERRESULT)
			{
				timeSeriesptr[iter] = static_cast<realdp>(getTickCount() - starttime) / CLK_PS;
				funSeriesptr[iter] = f2;
                gradSeriesptr[iter] = ngf2;
			}
		}
		ComTime = static_cast<realdp>(getTickCount() - starttime) / CLK_PS;
		if (Verbose >= ITERRESULT)
			lengthSeries = iter + 1;
        if (Verbose >= FINALRESULT)
            PrintFinalInfo();
	};

    void SolversSMTR::PrintFinalInfo(void)
    {
        printf("Iter:%d,f:%.3e,", iter, f1);
        
        printf("|gf|:%.3e,|gf|/|gf0|:%.3e,time:%.2e,nf:%d,ng:%d,nR:%d,", ngf1, ngf1 / ngf0, ComTime, nf, ng, nR);
        if (nH != 0)
        {
            printf("nH:%d,", nH);
        }
        if (nV != 0)
        {
            printf("nV(nVp):%d(%d),", nV, nVp);
        }
        printf("\n");
    };

	void SolversSMTR::InitialVector(void)
	{
		Mani->ScalarTimesVector(x1, 0, gf1, &eta1);
	};

	void SolversSMTR::tCG_TR(void)
	{
        Vector r(Prob->GetDomain()->GetEMPTY()), z(r), delta(r), Hd(r), Heta(r); /*Used for solving the local model*/
		realdp e_Pe, r_r, norm_r, norm_r0, d_Pd, z_r, e_Pd, d_Hd, alphatemp, e_Pe_new, tempnum, zold_rold, betatemp, tautemp;
		integer j;

		if (useRand)
		{
			HessianEta(eta1, &Heta); nH++;
            Mani->ScalarVectorAddVector(x1, 1, gf1, Heta, &r);
			e_Pe = Mani->Metric(x1, eta1, eta1);
		}
		else
		{
            Mani->ScalarTimesVector(x1, 1, eta1, &Heta);
            r = gf1;
			e_Pe = 0;
		}

		r_r = Mani->Metric(x1, r, r);
		norm_r = sqrt(r_r);
		norm_r0 = norm_r;

		Prob->PreConditioner(x1, r, &z);

		z_r = Mani->Metric(x1, z, r);
		d_Pd = z_r;
		Mani->ScalarTimesVector(x1, -1.0, z, &delta);
        
		if (useRand)
			e_Pd = Mani->Metric(x1, eta1, delta);
		else
			e_Pd = 0;

		tCGstatusSM = TRSM_MAXITER; /* pre-assume termination j == max_inner*/

        eta2 = eta1;

		realdp new_modelv = 0, modelv = 0;
		if (useRand)
			modelv = Mani->Metric(x1, eta1, gf1) + 0.5 * Mani->Metric(x1, eta1, Heta);
        
		for (j = 0; j < Max_Inner_Iter; j++)
		{
			HessianEta(delta, &Hd); nH++;
			d_Hd = Mani->Metric(x1, delta, Hd);
			alphatemp = z_r / d_Hd;
			e_Pe_new = e_Pe + static_cast<realdp> (2.0) * alphatemp * e_Pd + alphatemp * alphatemp * d_Pd;

			if (d_Hd <= 0 || e_Pe_new >= (Delta * Delta))
			{
				tautemp = (-e_Pd + sqrt(e_Pd * e_Pd + d_Pd * (Delta * Delta - e_Pe))) / d_Pd;
                Mani->ScalarVectorAddVector(x1, tautemp, delta, eta2, &eta2);
                
				if (d_Hd < 0)
					tCGstatusSM = TRSM_NEGCURVTURE; /* negative curvature*/
				else
					tCGstatusSM = TRSM_EXCREGION; /* exceeded trust region*/
				break;
			}
			e_Pe = e_Pe_new;
            
            Mani->ScalarVectorAddVector(x1, alphatemp, delta, eta2, &eta2);
            Mani->ScalarVectorAddVector(x1, alphatemp, Hd, Heta, &Heta);
            
			new_modelv = Mani->Metric(x1, eta2, gf1) + 0.5 * Mani->Metric(x1, eta2, Heta);
			if (new_modelv >= modelv)
			{
				tCGstatusSM = TRSM_ERROR;
				break;
			}
			modelv = new_modelv;
            Mani->ScalarVectorAddVector(x1, alphatemp, Hd, r, &r);
            Mani->Projection(x1, r, &r);

			r_r = Mani->Metric(x1, r, r);
			norm_r = sqrt(r_r);

			tempnum = pow(norm_r0, theta);

			if (j >= Min_Inner_Iter && norm_r <= norm_r0 * ((tempnum < kappa) ? tempnum : kappa))
			{
				if (kappa < tempnum)
					tCGstatusSM = TRSM_LCON; /* linear convergence*/
				else
					tCGstatusSM = TRSM_SCON; /* superlinear convergence*/
				break;
			}

			Prob->PreConditioner(x1, r, &z);

			zold_rold = z_r;
			z_r = Mani->Metric(x1, z, r);
			betatemp = z_r / zold_rold;
            Mani->ScalarTimesVector(x1, betatemp, delta, &delta);
            Mani->ScalarVectorAddVector(x1, -1, z, delta, &delta);
			e_Pd = betatemp * (e_Pd + alphatemp * d_Pd);
			d_Pd = z_r + betatemp * betatemp * d_Pd;
		}
		innerIter = j;
	};

	void SolversSMTR::PrintInfo(void)
	{
        printf("i:%d,f:%.3e,df/f:%.3e,", iter, f2, ((f1 - f2) / std::fabs(f2)));

        printf("|gf|:%.3e,time:%.2g,", ngf2, static_cast<realdp>(getTickCount() - starttime) / CLK_PS);

        printf("rho:%.2e,radius:%.3e,tCGstatus:%s,innerIter:%d,", rho, Delta, tCGstatusSetSMnames[tCGstatusSM].c_str(), innerIter);
        
        printf("nf:%d,ng:%d,", nf, ng);
        
        if (nH != 0)
            printf("nH:%d,", nH);
        
        printf("nR:%d,", nR);
        
        if (nV != 0)
            printf("nV(nVp):%d(%d),", nV, nVp);
        
        printf("\n");
	};

	void SolversSMTR::CheckParams(void)
	{
		SolversSM::CheckParams();

		char YES[] = "YES";
		char NO[] = "NO";
		char *status;

		printf("TRUST REGION TYPE METHODS PARAMETERS:\n");
		status = (initial_Delta > 0) ? YES : NO;
		printf("initial_Delta :%15g[%s],\t", initial_Delta, status);
		status = (Acceptence_Rho > 0 && Acceptence_Rho < 0.25) ? YES : NO;
		printf("Acceptence_Rho:%15g[%s],\n", Acceptence_Rho, status);
		status = (Shrinked_tau > 0 && Shrinked_tau < 1) ? YES : NO;
		printf("Shrinked_tau  :%15g[%s],\t", Shrinked_tau, status);
		status = (Magnified_tau > 1) ? YES : NO;
		printf("Magnified tau :%15g[%s],\n", Magnified_tau, status);
		status = (minimum_Delta > 0 && minimum_Delta <= maximum_Delta) ? YES : NO;
		printf("minimum_Delta :%15g[%s],\t", minimum_Delta, status);
		status = (maximum_Delta > 0 && maximum_Delta >= minimum_Delta) ? YES : NO;
		printf("maximum_Delta :%15g[%s],\n", maximum_Delta, status);
		status = (Min_Inner_Iter >= 0 && Min_Inner_Iter <= Max_Inner_Iter) ? YES : NO;
		printf("Min_Inner_Iter:%15d[%s],\t", Min_Inner_Iter, status);
		status = (Max_Inner_Iter >= 0 && Max_Inner_Iter >= Min_Inner_Iter) ? YES : NO;
		printf("Max_Inner_Iter:%15d[%s],\n", Max_Inner_Iter, status);
		status = (theta >= 0) ? YES : NO;
		printf("theta         :%15g[%s],\t", theta, status);
		status = (kappa > 0 && kappa < 1) ? YES : NO;
		printf("kappa         :%15g[%s],\n", kappa, status);
		status = YES;
		printf("useRand       :%15d[%s]\n", useRand, status);
	};

	void SolversSMTR::UpdateData(void)
	{
	};

	void SolversSMTR::Acceptence(void)
	{
		Prob->Grad(x2, &gf2); ng++;
	};

	void SolversSMTR::SetProbX(const Problem *prob, const Variable *initialx)
	{
		Solvers::SetProbX(prob, initialx);
        Heta2 = Prob->GetDomain()->GetEMPTY();
	};

	void SolversSMTR::SetDefaultParams()
	{
		SolversSM::SetDefaultParams();
		nH = 0;
		Acceptence_Rho = static_cast<realdp> (0.1);
		Shrinked_tau = static_cast<realdp> (0.25);
		Magnified_tau = static_cast<realdp> (2);
		minimum_Delta = std::numeric_limits<realdp>::epsilon();
		maximum_Delta = static_cast<realdp> (10000);
		useRand = false;
		Max_Inner_Iter = 1000;
		Min_Inner_Iter = 0;
		theta = static_cast<realdp> (1);
		kappa = static_cast<realdp> (0.1);
		initial_Delta = static_cast<realdp> (1);
		tCGstatusSetSMnames = new std::string[TCGSTATUSSETSMLENGTH];
		tCGstatusSetSMnames[TRSM_NEGCURVTURE].assign("TRSM_NEGCURVTURE");
		tCGstatusSetSMnames[TRSM_EXCREGION].assign("TRSM_EXCREGION");
		tCGstatusSetSMnames[TRSM_MIN].assign("TRSM__MIN");
		tCGstatusSetSMnames[TRSM_LCON].assign("TRSM_LCON");
		tCGstatusSetSMnames[TRSM_SCON].assign("TRSM_SCON");
		tCGstatusSetSMnames[TRSM_MAXITER].assign("TRSM_MAXITER");
		tCGstatusSetSMnames[TRSM_ERROR].assign("TRSM_ERROR");
	};

	SolversSMTR::~SolversSMTR(void)
	{
		delete[] tCGstatusSetSMnames;
	};

	void SolversSMTR::SetParams(PARAMSMAP params)
	{
		SolversSM::SetParams(params);
		PARAMSMAP::iterator iter;
		for (iter = params.begin(); iter != params.end(); iter++)
		{
			if (iter->first == static_cast<std::string> ("Acceptence_Rho"))
			{
				Acceptence_Rho = iter->second;
			}
			else
			if (iter->first == static_cast<std::string> ("Shrinked_tau"))
			{
				Shrinked_tau = iter->second;
			}
			else
			if (iter->first == static_cast<std::string> ("Magnified_tau"))
			{
				Magnified_tau = iter->second;
			}
			else
			if (iter->first == static_cast<std::string> ("minimum_Delta"))
			{
				minimum_Delta = iter->second;
			}
			else
			if (iter->first == static_cast<std::string> ("maximum_Delta"))
			{
				maximum_Delta = iter->second;
			}
			else
			if (iter->first == static_cast<std::string> ("useRand"))
			{
				useRand = ((static_cast<integer> (iter->second)) != 0);
			}
			else
			if (iter->first == static_cast<std::string> ("Max_Inner_Iter"))
			{
				Max_Inner_Iter = static_cast<integer> (iter->second);
			}
			else
			if (iter->first == static_cast<std::string> ("Min_Inner_Iter"))
			{
				Min_Inner_Iter = static_cast<integer> (iter->second);
			}
			else
			if (iter->first == static_cast<std::string> ("theta"))
			{
				theta = iter->second;
			}
			else
			if (iter->first == static_cast<std::string> ("kappa"))
			{
				kappa = iter->second;
			}
			else
			if (iter->first == static_cast<std::string> ("initial_Delta"))
			{
				initial_Delta = iter->second;
			}
		}
	};
}; /*end of ROPTLIB namespace*/

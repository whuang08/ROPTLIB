#include "Solvers/RNewton.h"

/*Define the namespace*/
namespace ROPTLIB{

	RNewton::RNewton(const Problem *prob, const Variable *initialx)
	{
		Initialization(prob, initialx);
	};

	void RNewton::SetProbX(const Problem *prob, const Variable *initialx)
	{
		SolversSMLS::SetProbX(prob, initialx);
		prob->SetUseGrad(true);
		prob->SetUseHess(true);
	};

	void RNewton::SetDefaultParams(void)
	{
		SolversSMLS::SetDefaultParams();
		theta = 1;
		kappa = static_cast<realdp> (0.1);
		Min_Inner_Iter = 0;
		Max_Inner_Iter = 1000;
		useRand = false;
		InitSteptype = LSSM_QUADINTMOD;

		SolverName.assign("RNewton");

		tCGLSSMstatusSetnames = new std::string[TCGLSSMSTATUSSETLENGTH];
		tCGLSSMstatusSetnames[LSSM_NEGCURVTURE].assign("LSSM_NEGCURVTURE");
		tCGLSSMstatusSetnames[LSSM_LCON].assign("LSSM_LCON");
		tCGLSSMstatusSetnames[LSSM_SCON].assign("LSSM_SCON");
		tCGLSSMstatusSetnames[LSSM_MAXITER].assign("LSSM_MAXITER");
		tCGLSSMstatusSetnames[LSSM_ERROR].assign("LSSM_ERROR");
	};

	void RNewton::CheckParams(void)
	{
		SolversSMLS::CheckParams();

		char YES[] = "YES";
		char NO[] = "NO";
		char *status;

		printf("RNEWTON METHOD PARAMETERS:\n");
		status = (Min_Inner_Iter >= 0 && Min_Inner_Iter <= Max_Inner_Iter) ? YES : NO;
		printf("Min_Inner_Iter:%15d[%s],\t", Min_Inner_Iter, status);
		status = (Max_Inner_Iter >= 0 && Max_Inner_Iter >= Min_Inner_Iter) ? YES : NO;
		printf("Max_Inner_Iter:%15d[%s],\n", Max_Inner_Iter, status);
		status = (theta >= 1) ? YES : NO;
		printf("theta         :%15g[%s],\t", theta, status);
		status = (kappa > 0 && kappa < 1) ? YES : NO;
		printf("kappa         :%15g[%s],\n", kappa, status);
		printf("useRand       :%15d[%s],\n", useRand, status);
	};

	RNewton::~RNewton(void)
	{
		delete[] tCGLSSMstatusSetnames;
	};

	void RNewton::GetSearchDir(void)
	{
		Mani->ScalarTimesVector(x1, 0, gf1, &eta1); /*initial iterate is a zero tangent vector*/
		tCG_LS();
	};

	void RNewton::PrintInfo(void)
	{
        printf("i:%d,f:%.3e,df/f:%.3e,", iter, f2, ((f1 - f2) / std::fabs(f2)));

        printf("|gf|:%.3e,t0:%.2e,t:%.2e,s0:%.2e,s:%.2e,time:%.2g,", ngf2, initiallength, stepsize, initialslope, newslope, static_cast<realdp>(getTickCount() - starttime) / CLK_PS);

        printf("nf:%d,ng:%d,", nf, ng);
        
        if (nH != 0)
            printf("nH:%d,", nH);
        
        printf("tCGstatus:%s,innerIter:%d,", tCGLSSMstatusSetnames[tCGLSSMstatus].c_str(), innerIter);
        
        printf("nR:%d,", nR);
        
        if (nV != 0)
            printf("nV(nVp):%d(%d),", nV, nVp);
        
        printf("\n");
	};

	void RNewton::tCG_LS(void)
	{
        Vector Heta(eta1), r(Prob->GetDomain()->GetEMPTY()), z(r), delta(r), Hd(r);
		realdp r_r, norm_r, norm_r0, z_r, d_Hd, alphatemp, tempnum, zold_rold, betatemp;
		integer j;

		if (useRand)
		{
			Prob->HessianEta(x1, eta1, &Heta); nH++;
            Mani->ScalarVectorAddVector(x1, 1, gf1, Heta, &r); /*r = gf1 + Heta*/
		}
		else
		{
            r = gf1;
		}

		r_r = Mani->Metric(x1, r, r);
		norm_r = sqrt(r_r);
		norm_r0 = norm_r;

		Prob->PreConditioner(x1, r, &z);

		z_r = Mani->Metric(x1, z, r);

		Mani->ScalarTimesVector(x1, -1.0, z, &delta);

		tCGLSSMstatus = LSSM_MAXITER; /* pre-assume termination j == max_inner*/

		realdp new_modelv = 0, modelv = 0;
		if (useRand)
			modelv = Mani->Metric(x1, eta1, gf1) + 0.5 * Mani->Metric(x1, eta1, Heta);

		for (j = 0; j < Max_Inner_Iter; j++)
		{
			Prob->HessianEta(x1, delta, &Hd); nH++;
			d_Hd = Mani->Metric(x1, delta, Hd);
			alphatemp = z_r / d_Hd;

			if (d_Hd <= 0)
			{
				tCGLSSMstatus = LSSM_NEGCURVTURE; /* negative curvature*/
				break;
			}
            Mani->ScalarVectorAddVector(x1, alphatemp, delta, eta1, &eta1); /*eta1 = eta1 + alphatemp * delta*/
            
			new_modelv = Mani->Metric(x1, eta1, gf1) + 0.5 * Mani->Metric(x1, eta1, Heta);
			if (new_modelv >= modelv || Mani->Metric(x1, eta1, gf1) >= -std::numeric_limits<realdp>::epsilon())
			{
				tCGLSSMstatus = LSSM_ERROR;
				break;
			}

			modelv = new_modelv;
            
            Mani->ScalarVectorAddVector(x1, alphatemp, Hd, r, &r); /*r = r + alphatemp * Hd*/
            Mani->Projection(x1, r, &r);

			r_r = Mani->Metric(x1, r, r);
			norm_r = sqrt(r_r);

			tempnum = pow(norm_r0, theta);

			if (j >= Min_Inner_Iter && norm_r <= norm_r0 * ((tempnum < kappa) ? tempnum : kappa))
			{
				if (kappa < tempnum)
					tCGLSSMstatus = LSSM_LCON; /* linear convergence*/
				else
					tCGLSSMstatus = LSSM_SCON; /* superlinear convergence*/
				break;
			}
            Prob->PreConditioner(x1, r, &z);
			zold_rold = z_r;
			z_r = Mani->Metric(x1, z, r);
			betatemp = z_r / zold_rold;
            Mani->ScalarTimesVector(x2, betatemp, delta, &delta);
            Mani->ScalarVectorAddVector(x2, -1, z, delta, &delta); /*delta =  - z + betatemp * delta*/
		}
		innerIter = j;
		if (j == 0)
            eta1 = delta;
	};

	void RNewton::SetParams(PARAMSMAP params)
	{
		SolversSMLS::SetParams(params);
		PARAMSMAP::iterator iter;
		for (iter = params.begin(); iter != params.end(); iter++)
		{
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
								kappa = static_cast<realdp> (iter->second);
							}
		}
	};
}; /*end of ROPTLIB namespace*/

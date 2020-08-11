
#include "Solvers/RCG.h"

/*Define the namespace*/
namespace ROPTLIB{

	RCG::RCG(const Problem *prob, const Variable *initialx)
	{
		Initialization(prob, initialx);
	};

	void RCG::SetProbX(const Problem *prob, const Variable *initialx)
	{
		SolversSMLS::SetProbX(prob, initialx);
		prob->SetUseGrad(true);
		prob->SetUseHess(false);
	};

	void RCG::SetDefaultParams(void)
	{
		SolversSMLS::SetDefaultParams();
		sigma = 0;
		RCGmethod = HESTENES_STIEFEL;
		InitSteptype = LSSM_BBSTEP;

		SolverName.assign("RCG");

		RCGmethodSetnames = new std::string[RCGMETHODSLENGTH];
		RCGmethodSetnames[FLETCHER_REEVES].assign("FLETCHER_REEVES");
		RCGmethodSetnames[POLAK_RIBIERE_MOD].assign("POLAK_RIBIERE_MOD");
		RCGmethodSetnames[HESTENES_STIEFEL].assign("HESTENES_STIEFEL");
		RCGmethodSetnames[FR_PR].assign("FR_PR");
		RCGmethodSetnames[DAI_YUAN].assign("DAI_YUAN");
		RCGmethodSetnames[HAGER_ZHANG].assign("HAGER_ZHANG");
	};

	RCG::~RCG(void)
	{
		delete[] RCGmethodSetnames;
	};

	void RCG::CheckParams(void)
	{
		SolversSMLS::CheckParams();
		char YES[] = "YES";
		char NO[] = "NO";
		char *status;

		printf("RCG METHOD PARAMETERS:\n");
		status = (RCGmethod >= 0 && RCGmethod <= RCGMETHODSLENGTH) ? YES : NO;
		printf("RCGmethod     :%15s[%s]\n", RCGmethodSetnames[RCGmethod].c_str(), status);
	};

	void RCG::PrintInfo(void)
	{
        printf("i:%d,f:%.3e,df/f:%.3e,", iter, f2, ((f1 - f2) / std::fabs(f2)));

        printf("|gf|:%.3e,sigma:%.2e,t0:%.2e,t:%.2e,s0:%.2e,s:%.2e,time:%.2g,", ngf2, sigma, initiallength, stepsize, initialslope, newslope, static_cast<realdp>(getTickCount() - starttime) / CLK_PS);

        printf("nf:%d,ng:%d,", nf, ng);
        
        if (nH != 0)
            printf("nH:%d,", nH);
        
        printf("nR:%d,", nR);
        
        if (nV != 0)
            printf("nV(nVp):%d(%d),", nV, nVp);
        
        printf("\n");
	};

	void RCG::GetSearchDir(void)
	{
		if (iter == 0 || Mani->Metric(x1, eta1, gf1) / ngf1 / ngf1 >= -std::sqrt(std::numeric_limits<realdp>::epsilon())) /* restart and safeguard */
		{
			Prob->PreConditioner(x1, gf1, &Pgf1); /* Use preconditioner */
            Mani->ScalarTimesVector(x1, -1.0, Pgf1, &eta1);

            if (Mani->Metric(x1, eta1, gf1) / ngf1 / ngf1 >= -std::sqrt(std::numeric_limits<realdp>::epsilon()))
            { /* Without preconditioner */
                Mani->ScalarTimesVector(x1, -1, gf1, &eta1);
            }
		}
	};

	void RCG::UpdateData(void)
	{
        Vector Tgf1(gf1); Mani->VectorTransport(x1, eta2, x2, gf1, &Tgf1); nV++;
        Vector Teta1 = Prob->GetDomain()->GetEMPTY();
        Prob->PreConditioner(x2, gf2, &Pgf2);
        if(std::fabs(Mani->Metric(x2, gf2, Tgf1)) / Mani->Metric(x1, gf1, gf1) > 0.1)
        { /*Restart using the safeguard (5.52) in [NW06] */
            sigma = 0;
        }
        else
        {
            switch (RCGmethod)
            {
                case FLETCHER_REEVES:
                    sigma = Mani->Metric(x2, gf2, Pgf2) / Mani->Metric(x1, gf1, Pgf1);
                    Mani->VectorTransport(x1, eta2, x2, eta1, &Teta1); nVp++;
                    break;
                case POLAK_RIBIERE_MOD:
                    sigma = Mani->Metric(x2, Pgf2, Mani->VectorLinearCombination(x2, 1, gf2, -1, Tgf1, &Teta1)) / Mani->Metric(x1, gf1, Pgf1); /*The Teat1 here is used as a temporary space*/
                     Mani->VectorTransport(x1, eta2, x2, eta1, &Teta1); nVp++;
                    sigma = (sigma < 0) ? 0 : sigma;
                    break;
                case HESTENES_STIEFEL:
                    Mani->VectorTransport(x1, eta2, x2, eta1, &Teta1); nVp++;
                    Mani->VectorLinearCombination(x2, 1, gf2, -1, Tgf1, &Tgf1);
                    sigma = Mani->Metric(x2, Pgf2, Tgf1) / Mani->Metric(x2, Tgf1, Teta1);
                    break;
                case FR_PR:
                {
                    Mani->VectorTransport(x1, eta2, x2, eta1, &Teta1); nVp++;
                    realdp sigmaFR = Mani->Metric(x2, gf2, Pgf2) / Mani->Metric(x1, gf1, Pgf1);
                    Mani->VectorLinearCombination(x2, 1, gf2, -1, Tgf1, &Tgf1);
                    realdp sigmaPR = Mani->Metric(x2, Pgf2, Tgf1) / Mani->Metric(x1, gf1, Pgf1);
                    sigma = (sigmaPR < -sigmaFR) ? - sigmaFR : ((sigmaPR > sigmaFR) ? sigmaFR : sigmaPR);
                }
                    break;
                case DAI_YUAN:
                    Mani->VectorTransport(x1, eta2, x2, eta1, &Teta1); nVp++;
                    Mani->VectorLinearCombination(x2, 1, gf2, -1, Tgf1, &Tgf1);
                    sigma = Mani->Metric(x2, gf2, Pgf2) / Mani->Metric(x2, Tgf1, Teta1);
                    break;
                case HAGER_ZHANG:
                {
                    Vector y(Tgf1);Mani->VectorLinearCombination(x2, 1, gf2, -1, Tgf1, &y);
                    Mani->VectorTransport(x1, eta2, x2, eta1, &Teta1); nVp++;
                    realdp yp = Mani->Metric(x2, y, Teta1);
                    Mani->VectorLinearCombination(x2, 1, y, - 2.0 * Mani->Metric(x2, y, y) / yp, Teta1, &Tgf1);
                    sigma = Mani->Metric(x2, Tgf1, Pgf2) / yp;
                }
                    break;
                default:
                    printf("Warning: RCG scheme does not exists!");
                    sigma = 0;
            }
        }
        Pgf1 = Pgf2;
        if(sigma == 0)
            Mani->ScalarTimesVector(x2, -1, gf2, &eta1);
        else
            Mani->VectorLinearCombination(x2, -1, gf2, sigma, Teta1, &eta1);
	};

	void RCG::SetParams(PARAMSMAP params)
	{
		SolversSMLS::SetParams(params);
		PARAMSMAP::iterator iter;
		for (iter = params.begin(); iter != params.end(); iter++)
		{
            if (iter->first == static_cast<std::string> ("RCGmethod"))
            {
                RCGmethod = static_cast<RCGmethods> (static_cast<integer> (iter->second));
            }
		}
	};
}; /*end of ROPTLIB namespace*/

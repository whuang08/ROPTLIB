
#include "Solvers/SolversNSMPGLS.h"

/*Define the namespace*/
namespace ROPTLIB{

    SolversNSMPGLS::~SolversNSMPGLS()
    {
		delete[] LSstatusSetnames;
    };

	void SolversNSMPGLS::SetProbX(const Problem *prob, const Variable *initialx)
	{
		SolversNSM::SetProbX(prob, initialx);
		prob->SetUseGrad(true);
		prob->SetUseHess(false);
        dimNorVec = Prob->GetDomain()->GetExtrDim() - Prob->GetDomain()->GetIntrDim();
        eta1 = Prob->GetDomain()->GetEMPTY();
        eta2 = Prob->GetDomain()->GetEMPTY();
	};

	void SolversNSMPGLS::SetDefaultParams()
	{
		SolversNSM::SetDefaultParams();
        Minstepsize = 2e-2;
        LS_ratio = 0.5;
        LS_alpha = 0.0001;
        SMlambda = 0.2;
        SMtol = 1e-10;
        Tolerance = 1e-4;
        adavalue = 1;
        Variant = LSPG_ADALIPSCHITZ;
        totalSMiter = 0;
        SMCGiter = 0;
        totalSMCGiter = 0;
        LSstatusSetnames = new std::string[LSSTATUSSETPGLENGTH];
        LSstatusSetnames[LSPG_MINSTEPSIZE].assign("LSPG_MINSTEPSIZE");
        LSstatusSetnames[LSPG_SUCCESS].assign("LSPG_SUCCESS");
	};

    void SolversNSMPGLS::CheckParams(void)
    {
        SolversNSM::CheckParams();

        std::string RPGVariantnames[RPGVARIANTLENGTH] = { "REGULAR", "ADALIPSCHITZ" };

        char YES[] = "YES";
        char NO[] = "NO";
        char *status;

        printf("RPG LINESEARCH-BASED METHODS PARAMETERS:\n");
        status = (Variant >= 0 && Variant < RPGVARIANTLENGTH) ? YES : NO;
        printf("RPGVariant    :%15s[%s],\t", RPGVariantnames[Variant].c_str(), status);
        status = (Minstepsize > 0 && Minstepsize <= 1) ? YES : NO;
        printf("Minstepsize   :%15g[%s],\n", Minstepsize, status);
        status = (LS_ratio > 0 && LS_ratio < 1) ? YES : NO;
        printf("LS_ratio      :%15g[%s],\t", LS_ratio, status);
        status = (LS_alpha > 0 && LS_alpha < 1) ? YES : NO;
        printf("LS_alpha      :%15g[%s],\n", LS_alpha, status);
    };

    void SolversNSMPGLS::LinesearchArmijo(void)
    {
        LSstatus = LSPG_SUCCESS;
        f2 = h();
        realdp maxpref = f1;
        /* simple backtracking */
        while (maxpref - f2 < - LS_alpha * initialslope * stepsize)
        {
            stepsize *= LS_ratio;
            if (stepsize < Minstepsize)
            {
                if (Verbose > FINALRESULT)
                {
                    printf("Warning: step size reaches the minimum: %3.2e!\n", Minstepsize);
                }
                LSstatus = LSPG_MINSTEPSIZE;
                break;
            }
            f2 = h();
        }
        newslope = dh();
    };

    realdp SolversNSMPGLS::h(void)
    {
        Mani->ScalarTimesVector(x1, stepsize, eta1, &eta2);
        Mani->Retraction(x1, eta2, &x2); nR++;
        nf++;
        return Prob->f(x2);
    };

    realdp SolversNSMPGLS::dh(void)
    {
        Prob->Grad(x2, &gf2); ng++;
        nV++;
        Vector diffeta1(eta1); Mani->DiffRetraction(x1, eta2, x2, eta1, &diffeta1, true);
        return Mani->Metric(x2, gf2, diffeta1);
    };

    void SolversNSMPGLS::UpdateData(void)
    {
    };

    void SolversNSMPGLS::SetParams(PARAMSMAP params)
    {
        SolversNSM::SetParams(params);
        PARAMSMAP::iterator iter;
        for (iter = params.begin(); iter != params.end(); iter++)
        {
            if (iter->first == static_cast<std::string> ("RPGLSVariant"))
            {
                Variant = static_cast<RPGLSVariant> (static_cast<integer> (iter->second));
            }
            else
            if (iter->first == static_cast<std::string> ("LS_ratio"))
            {
                LS_ratio = static_cast<integer> (iter->second);
            }
            else
            if (iter->first == static_cast<std::string> ("LS_alpha"))
            {
                LS_alpha = static_cast<integer> (iter->second);
            }
            else
            if (iter->first == static_cast<std::string> ("Minstepsize"))
            {
                Minstepsize = static_cast<integer> (iter->second);
            }
        }
    };
}; /*end of ROPTLIB namespace*/

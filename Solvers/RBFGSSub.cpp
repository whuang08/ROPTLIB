
#include "Solvers/RBFGSSub.h"

/*Define the namespace*/
namespace ROPTLIB{

	RBFGSSub::RBFGSSub(const Problem *prob, const Variable *initialx, LinearOPE *initialH)
	{
		Initialization(prob, initialx, initialH);
	};

	void RBFGSSub::Initialization(const Problem *prob, const Variable *initialx, LinearOPE *initialH)
	{
		SetProbX(prob, initialx, initialH);
		SetDefaultParams();
	};

	void RBFGSSub::SetProbX(const Problem *prob, const Variable *initialx, LinearOPE *initialH)
	{
        SolversNSMSubLS::SetProbX(prob, initialx);
        
        bool initHisnull = (initialH == nullptr);
        if (initHisnull)
        {
            if (prob->GetDomain()->GetIsIntrinsic())
            {
                initialH = new LinearOPE(prob->GetDomain()->GetEMPTYINTR().Getlength(), prob->GetDomain()->GetEMPTYINTR().Getlength());
            }
            else
            {
                initialH = new LinearOPE(prob->GetDomain()->GetEMPTYEXTR().Getlength(), prob->GetDomain()->GetEMPTYEXTR().Getlength());
            }
            initialH->ScaledIdOPE();
        }
        H = *initialH;
        if (initHisnull)
            delete initialH;
        prob->SetUseGrad(true);
        prob->SetUseHess(false);
        s = Prob->GetDomain()->GetEMPTY();
        y = Prob->GetDomain()->GetEMPTY();
        Py = Prob->GetDomain()->GetEMPTY();
	};

	void RBFGSSub::SetDefaultParams()
	{
		SolversNSMSubLS::SetDefaultParams();
		SolversNSMSubLS::SolverName.assign("RBFGSSub");
	};

    Vector &RBFGSSub::HvSub(const Vector &v, Vector *result)
    {
        return HvRBFGSSub(v, H, result);
    };

	RBFGSSub::~RBFGSSub(void)
	{
	};

    Vector &RBFGSSub::HvRBFGSSub(const Vector &v, const LinearOPE &H, Vector *result)
    {
        nH++;
        return Mani->LinearOPEEta(x1, H, v, result);
    };

    void RBFGSSub::UpdateDataRBFGSSub(void)
    {
        Mani->VectorTransport(x1, eta2, x2, eta2, &s); nV++;
        Vector zeta(gf1); Mani->VectorTransport(x1, eta2, x2, gf1, &zeta); nVp++;
        betay = Mani->Beta(x1, eta2);
        Mani->VectorLinearCombination(x2, static_cast<realdp> (1) / betay, gf2, -1, zeta, &y);
        inpyy = Mani->Metric(x2, y, y);
        inpsy = Mani->Metric(x2, s, y);
        realdp tmp = static_cast<realdp> (1) / lambdaUpper - inpsy / inpyy;
        tmp = (tmp > 0) ? tmp : 0;
        Mani->ScalarVectorAddVector(x2, tmp, y, s, &s);
        inpsy = Mani->Metric(x2, s, y);

        if (iter == 1 && inpsy > 0)
            H.ScaledIdOPE(inpsy / inpyy);

        Mani->TranHInvTran(x1, eta2, x2, H, &H);
        inpss = Mani->Metric(x2, s, s);
        if (inpsy / inpss > lambdaLower)
        {
            Mani->LinearOPEEta(x2, H, y, &zeta);
            Mani->HaddScaledRank1OPE(x2, H, static_cast<realdp> (-1) / inpsy, s, zeta, &H);
            Mani->LinearOPEEta(x2, H, y, &zeta);

            Mani->HaddScaledRank1OPE(x2, H, static_cast<realdp> (-1) / inpsy, zeta, s, &H);
            Mani->HaddScaledRank1OPE(x2, H, static_cast<realdp> (1) / inpsy, s, s, &H);
            isupdated = true;
        }
        else
        {
            isupdated = false;
            H.ScaledIdOPE(1);
        }
    };

	void RBFGSSub::PrintInfo(void)
	{
        printf("i:%d,f:%.3e,df/f:%.3e,", iter, f2, ((f1 - f2) / std::fabs(f2)));

        printf("|nd|:%.3e,t0:%.2e,t:%.2e,s0:%.2e,s:%.2e,time:%.2e,", ndir1, initiallength, stepsize, initialslope, newslope,  static_cast<realdp>(getTickCount() - starttime) / CLK_PS);
        
        printf("betay:%.3e,inpss:%.3e,inpsy:%.3e,inpyy:%.3e,IsUpdateHessian:%d,", betay, inpss, inpsy, inpyy, isupdated);
        
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

	void RBFGSSub::UpdateData(void)
	{
        UpdateDataRBFGSSub();
	};
}; /*end of ROPTLIB namespace*/

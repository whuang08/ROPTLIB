
#include "Solvers/RTRSR1.h"

/*Define the namespace*/
namespace ROPTLIB{

	RTRSR1::RTRSR1(const Problem *prob, const Variable *initialx, LinearOPE *initialB)
	{
		Initialization(prob, initialx, initialB);
	};

	void RTRSR1::Initialization(const Problem *prob, const Variable *initialx, LinearOPE *initialB)
	{
		SetProbX(prob, initialx, initialB);
		SetDefaultParams();
	};

	void RTRSR1::SetProbX(const Problem *prob, const Variable *initialx, LinearOPE *initialB)
	{
		SolversSMTR::SetProbX(prob, initialx);
        
        bool initBisnull = (initialB == nullptr);
        if (initBisnull)
        {
            if (prob->GetDomain()->GetIsIntrinsic())
            {
                initialB = new LinearOPE(prob->GetDomain()->GetEMPTYINTR().Getlength(), prob->GetDomain()->GetEMPTYINTR().Getlength());
            }
            else
            {
                initialB = new LinearOPE(prob->GetDomain()->GetEMPTYEXTR().Getlength(), prob->GetDomain()->GetEMPTYEXTR().Getlength());
            }
            initialB->ScaledIdOPE();
        }
        B = *initialB;
        if (initBisnull)
            delete initialB;
        prob->SetUseGrad(true);
        prob->SetUseHess(false);
        s = Prob->GetDomain()->GetEMPTY();
        y = Prob->GetDomain()->GetEMPTY();
	};

	void RTRSR1::SetDefaultParams(void)
	{
		SolversSMTR::SetDefaultParams();
		theta = static_cast<realdp> (0.1);
		kappa = static_cast<realdp> (0.1);
		SolverName.assign("RTRSR1");
	};

	RTRSR1::~RTRSR1(void)
	{
	};

	Vector &RTRSR1::HessianEta(const Vector &Eta, Vector *result)
	{
		return HvRTRSR1(Eta, B, result);
	};

	void RTRSR1::UpdateData(void)
	{
		UpdateDataRTRSR1();
	};

	void RTRSR1::Acceptence(void)
	{
		Mani->TranHInvTran(x1, eta2, x2, B, &B);
	};

    Vector &RTRSR1::HvRTRSR1(const Vector &v, const LinearOPE &B, Vector *result)
    {
        return Mani->LinearOPEEta(x1, B, v, result);
    };

    void RTRSR1::UpdateDataRTRSR1(void)
    {
        realdp denorminator, norm2ymBs;
        realdp mintolsq = std::numeric_limits<realdp>::epsilon();
        Prob->Grad(x2, &gf2); ng++;
        Mani->InverseVectorTransport(x1, eta2, x2, gf2, &y);  nV++;
        Mani->ScalarVectorAddVector(x1, -1, gf1, y, &y);
        s = eta2;
        Vector zeta(y); Mani->ScalarVectorAddVector(x1, -1, Heta2, y, &zeta);
        denorminator = Mani->Metric(x1, s, zeta);
        
        if (iter == 1)
            B.ScaledIdOPE(Mani->Metric(x1, y, y) / Mani->Metric(x1, s, y));
        
        inpss = Mani->Metric(x1, s, s);
        norm2ymBs = Mani->Metric(x1, zeta, zeta);
        if (denorminator * denorminator >= mintolsq * inpss * norm2ymBs && (norm2ymBs >= mintolsq || ngf2 / ngf0 < 1e-3))
        {
            Mani->HaddScaledRank1OPE(x1, B, static_cast<realdp> (1) / denorminator, zeta, zeta, &B);
            isupdated = true;
        }
        else
        {
            isupdated = false;
        }
    };

	void RTRSR1::PrintInfo(void)
	{
        printf("i:%d,f:%.3e,df/f:%.3e,", iter, f2, ((f1 - f2) / std::fabs(f2)));
        
        printf("|gf|:%.3e,time:%.2g,", ngf2, static_cast<realdp>(getTickCount() - starttime) / CLK_PS);

        printf("rho:%.2e,radius:%.3e,tCGstatus:%s,innerIter:%d,", rho, Delta, tCGstatusSetSMnames[tCGstatusSM].c_str(), innerIter);
        
        printf("\n\tinpss:%.3e,IsUpdateHessian:%d,", inpss, isupdated);
        
        printf("nf:%d,ng:%d,", nf, ng);
        
        if (nH != 0)
            printf("nH:%d,", nH);
        
        printf("nR:%d,", nR);
        
        if (nV != 0)
            printf("nV(nVp):%d(%d),", nV, nVp);
        
        printf("\n");
	};
}; /*end of ROPTLIB namespace*/

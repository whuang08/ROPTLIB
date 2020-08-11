
#include "Solvers/LRBFGSSub.h"

/*Define the namespace*/
namespace ROPTLIB{

	void LRBFGSSub::Run(void)
	{
        DeleteVectors(S, LengthSY);
        NewVectors(S, LengthSY);
        DeleteVectors(Y, LengthSY);
        NewVectors(Y, LengthSY);
        if (RHO != nullptr)
            delete[] RHO;
        RHO = new realdp[LengthSY];
        SolversNSMSubLS::Run();
	};

	LRBFGSSub::LRBFGSSub(const Problem *prob, const Variable *initialx)
	{
		Initialization(prob, initialx);
	};

	void LRBFGSSub::Initialization(const Problem *prob, const Variable *initialx)
	{
		SetProbX(prob, initialx);
		SetDefaultParams();
	};

	void LRBFGSSub::SetProbX(const Problem *prob, const Variable *initialx)
	{
		SolversNSMSubLS::SetProbX(prob, initialx);
		prob->SetUseGrad(true);
		prob->SetUseHess(false);
        s = Prob->GetDomain()->GetEMPTY();
        y = Prob->GetDomain()->GetEMPTY();
	};

	void LRBFGSSub::SetDefaultParams()
	{
		SolversNSMSubLS::SetDefaultParams();
		LengthSY = 4;
		S = nullptr;
		Y = nullptr;
		RHO = nullptr;
		Currentlength = 0;
		beginidx = 0;
		gamma = 1;
		SolversNSMSubLS::SolverName.assign("LRBFGSLPSub");
	};

    void LRBFGSSub::SetParams(PARAMSMAP params)
    {
        SolversNSMSubLS::SetParams(params);
        PARAMSMAP::iterator iter;
        for (iter = params.begin(); iter != params.end(); iter++)
        {
            if (iter->first == static_cast<std::string> ("LengthSY"))
            {
                LengthSY = static_cast<integer> (iter->second);
            }
        }
    };

	LRBFGSSub::~LRBFGSSub(void)
	{
		DeleteVectors(S, LengthSY);
		DeleteVectors(Y, LengthSY);
		if (RHO != nullptr)
			delete[] RHO;
	};

    Vector &LRBFGSSub::HvSub(const Vector &v, Vector *result)
    {
        return HvLRBFGSSub(v, result);
    };

	void LRBFGSSub::PrintInfo(void)
	{
        printf("i:%d,f:%.3e,df/f:%.3e,", iter, f2, ((f1 - f2) / std::fabs(f2)));

        printf("|nd|:%.3e,t0:%.2e,t:%.2e,s0:%.2e,s:%.2e,time:%.2e,", ndir1, initiallength, stepsize, initialslope, newslope,  static_cast<realdp>(getTickCount() - starttime) / CLK_PS);
        
        printf("betay:%.3e,rho:%.3e,gamma:%.3e,inpss:%.3e,inpsy:%.3e,inpyy:%.3e,IsUpdateHessian:%d,", betay, rho, gamma, inpss, inpsy, inpyy, isupdated);
        
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

	void LRBFGSSub::CheckParams(void)
	{
		SolversNSMSubLS::CheckParams();
		char YES[] = "YES";
		char NO[] = "NO";
		char *status;

		printf("LRBFGSLPSub METHOD PARAMETERS:\n");
		status = (LengthSY >= 0) ? YES : NO;
		printf("LengthSY      :%15d[%s]\n", LengthSY, status);
	};

	void LRBFGSSub::UpdateData(void)
	{
        UpdateDataLRBFGSSub();
	};

    Vector &LRBFGSSub::HvLRBFGSSub(const Vector &v, Vector *result)
    {
        realdp *xi = new realdp[Currentlength];
        realdp omega;
        integer idx;
        *result = v;
        for (integer i = Currentlength - 1; i >= 0; i--)
        {
            idx = (beginidx + i) % LengthSY;
            xi[idx] = RHO[idx] * Mani->Metric(x1, S[idx], *result);
            Mani->ScalarVectorAddVector(x1, -xi[idx], Y[idx], *result, result);
        }
        Mani->ScalarTimesVector(x1, gamma, *result, result);
        for (integer i = 0; i < Currentlength; i++)
        {
            idx = (beginidx + i) % LengthSY;
            omega = RHO[idx] * Mani->Metric(x1, Y[idx], *result);
            Mani->ScalarVectorAddVector(x1, xi[idx] - omega, S[idx], *result, result);
        }

        delete[] xi;
        return *result;
    };

    void LRBFGSSub::UpdateDataLRBFGSSub(void)
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
        inpss = Mani->Metric(x2, s, s);
        rho = static_cast<realdp> (1) / inpsy;
        if (inpsy / inpss >= lambdaLower)
        {
            gamma = inpsy / inpyy; /*Suggested in NW2006*/
            if (Currentlength < LengthSY)
            {
                Y[Currentlength] = y;
                S[Currentlength] = s;
                RHO[Currentlength] = rho;
                for (integer i = 0; i < Currentlength; i++)
                {
                    Mani->VectorTransport(x1, eta2, x2, Y[i], &Y[i]); nVp++;
                    Mani->VectorTransport(x1, eta2, x2, S[i], &S[i]); nVp++;
                }
                Currentlength++;
            }
            else
                if (LengthSY > 0)
                {
                    integer idx;
                    Y[beginidx] = y;
                    S[beginidx] = s;
                    RHO[beginidx] = rho;
                    beginidx = (++beginidx) % LengthSY;
                    for (integer i = beginidx; i < beginidx + LengthSY - 1; i++)
                    {
                        idx = i % LengthSY;
                        Mani->VectorTransport(x1, eta2, x2, Y[idx], &Y[idx]); nVp++;
                        Mani->VectorTransport(x1, eta2, x2, S[idx], &S[idx]); nVp++;
                    }
                }
            isupdated = true;
        }
        else
        {
            Currentlength = 0;
            beginidx = 0;
            for (integer i = 0; i < Currentlength; i++)
            {
                Mani->VectorTransport(x1, eta2, x2, Y[i], &Y[i]); nVp++;
                Mani->VectorTransport(x1, eta2, x2, S[i], &S[i]); nVp++;
            }
            isupdated = false;
        }
    };

}; /*end of ROPTLIB namespace*/


#include "Solvers/LRBFGS.h"

/*Define the namespace*/
namespace ROPTLIB{

	LRBFGS::LRBFGS(const Problem *prob, const Variable *initialx)
	{
		Initialization(prob, initialx);
	};

	void LRBFGS::SetProbX(const Problem *prob, const Variable *initialx)
	{
		SolversSMLS::SetProbX(prob, initialx);
		prob->SetUseGrad(true);
		prob->SetUseHess(false);
        s = Prob->GetDomain()->GetEMPTY();
        y = Prob->GetDomain()->GetEMPTY();
        Py = Prob->GetDomain()->GetEMPTY();
	};

	void LRBFGS::SetDefaultParams(void)
	{
		SolversSMLS::SetDefaultParams();
		isconvex = false;
		nu = static_cast<realdp> (1e-4);
		mu = 1;
		LengthSY = 4;
        LMrestart = false;
		S = nullptr;
		Y = nullptr;
        RHO = nullptr;
		Currentlength = 0;
		beginidx = 0;
		gamma = 1;
		InitSteptype = LSSM_QUADINTMOD;
		SolverName.assign("LRBFGS");
	};

    void LRBFGS::SetParams(PARAMSMAP params)
    {
        SolversSMLS::SetParams(params);
        PARAMSMAP::iterator iter;
        for (iter = params.begin(); iter != params.end(); iter++)
        {
            if (iter->first == static_cast<std::string> ("isconvex"))
            {
                isconvex = ((static_cast<integer> (iter->second)) != 0);
            }
            else
            if (iter->first == static_cast<std::string> ("LengthSY"))
            {
                LengthSY = static_cast<integer> (iter->second);
            }
            else
            if (iter->first == static_cast<std::string> ("LMrestart"))
            {
                LMrestart = static_cast<integer> (iter->second);
            }
            else
            if (iter->first == static_cast<std::string> ("nu"))
            {
                nu = iter->second;
            }
            else
            if (iter->first == static_cast<std::string> ("mu"))
            {
                mu = iter->second;
            }
        }
    };

	LRBFGS::~LRBFGS(void)
	{
		DeleteVectors(S, LengthSY);
		DeleteVectors(Y, LengthSY);
		if (RHO != nullptr)
			delete[] RHO;
	};

	void LRBFGS::Run(void)
	{
		DeleteVectors(S, LengthSY);
		NewVectors(S, LengthSY);
		DeleteVectors(Y, LengthSY);
		NewVectors(Y, LengthSY);
		if (RHO != nullptr)
			delete[] RHO;
		RHO = new realdp[LengthSY];
		SolversSMLS::Run();
	};

	void LRBFGS::CheckParams(void)
	{
		SolversSMLS::CheckParams();
		char YES[] = "YES";
		char NO[] = "NO";
		char *status;

		printf("LRBFGS METHOD PARAMETERS:\n");
		status = (nu >= 0 && nu < 1) ? YES : NO;
		printf("nu            :%15g[%s],\t", nu, status);
		status = (mu >= 0) ? YES : NO;
		printf("mu            :%15g[%s],\n", mu, status);
		status = YES;
		printf("isconvex      :%15d[%s],\t", isconvex, status);
		status = (LengthSY >= 0) ? YES : NO;
		printf("LengthSY      :%15d[%s],\n", LengthSY, status);
        status = YES;
        printf("LMrestart     :%15d[%s],\n", LMrestart, status);
	};

	void LRBFGS::GetSearchDir(void)
	{
        HvLRBFGS(gf1, &eta1);
		Mani->ScalarTimesVector(x1, -1.0, eta1, &eta1);
	};

	void LRBFGS::UpdateData(void)
	{
		UpdateDataLRBFGS();
	};

	void LRBFGS::PrintInfo(void)
	{
        printf("i:%d,f:%.3e,df/f:%.3e,", iter, f2, ((f1 - f2) / std::fabs(f2)));

        printf("|gf|:%.3e,t0:%.2e,t:%.2e,s0:%.2e,s:%.2e,time:%.2g,", ngf2, initiallength, stepsize, initialslope, newslope, static_cast<realdp>(getTickCount() - starttime) / CLK_PS);

        printf("\n\tbetay:%.3e,rho:%.3e,gamma:%.3e,inpss:%.3e,inpsy:%.3e,inpyy:%.3e,IsUpdateHessian:%d,", betay, rho, gamma, inpss, inpsy, inpyy, isupdated);
        
        printf("nf:%d,ng:%d,", nf, ng);
        
        if (nH != 0)
            printf("nH:%d,", nH);
        
        printf("nR:%d,", nR);
        
        if (nV != 0)
            printf("nV(nVp):%d(%d),", nV, nVp);
        
        printf("\n");
	};

    Vector &LRBFGS::HvLRBFGS(const Vector &v, Vector *result)
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
        Prob->PreConditioner(x1, *result, &Py);
        Mani->ScalarTimesVector(x1, gamma, Py, result);
        
        for (integer i = 0; i < Currentlength; i++)
        {
            idx = (beginidx + i) % LengthSY;
            omega = RHO[idx] * Mani->Metric(x1, Y[idx], *result);
            Mani->ScalarVectorAddVector(x1, xi[idx] - omega, S[idx], *result, result);
        }
        delete[] xi;
        return *result;
    };

    realdp LRBFGS::InitialHessian(realdp inpss, realdp inpsy, realdp inpyy)
    { /*Suggested in NW2006*/
        return inpsy / inpyy;
    };

    void LRBFGS::UpdateDataLRBFGS(void)
    {
        Mani->VectorTransport(x1, eta2, x2, eta2, &s); nV++;
        Vector Tgf1(gf1); Mani->VectorTransport(x1, eta2, x2, gf1, &Tgf1); nVp++;
        betay = Mani->Beta(x1, eta2);
        Mani->VectorLinearCombination(x2, static_cast<realdp> (1) / betay, gf2, -1, Tgf1, &y);
        Prob->PreConditioner(x2, y, &Py);
        
        inpsy = Mani->Metric(x2, s, y);
        inpss = Mani->Metric(x2, s, s);
        inpyy = Mani->Metric(x2, y, Py);
        rho = static_cast<realdp> (1) / inpsy;
        if (inpsy / inpss >= nu * pow(ngf2, mu) && (ngf2 / ngf0 < 1e-3 ||
            (inpss > std::numeric_limits<realdp>::epsilon() && inpsy > std::numeric_limits<realdp>::epsilon())))
        {
            gamma = InitialHessian(inpss, inpsy, inpyy);
            
            if(LMrestart && Currentlength >= LengthSY)
                Currentlength = 0;
            
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
            for (integer i = 0; i < Currentlength; i++)
            {
                Mani->VectorTransport(x1, eta2, x2, Y[i], &Y[i]); nVp++;
                Mani->VectorTransport(x1, eta2, x2, S[i], &S[i]); nVp++;
            }
            isupdated = false;
        }
    };

}; /*end of ROPTLIB namespace*/

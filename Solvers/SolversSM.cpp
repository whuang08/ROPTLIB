
#include "Solvers/SolversSM.h"

/*Define the namespace*/
namespace ROPTLIB{

	void SolversSM::Run(void)
	{
		Solvers::Run();
		if (Verbose >= ITERRESULT)
		{
			gradSeries = Vector(1 + Max_Iteration);
		}
	};

	void SolversSM::PrintInfo(void)
	{
        printf("i:%d,f:%.3e,df/f:%.3e,", iter, f2, ((f1 - f2) / std::fabs(f2)));

        printf("|gf|:%.3e,time:%.2g,", ngf2, static_cast<realdp>(getTickCount() - starttime) / CLK_PS);

        printf("nf:%d,ng:%d,", nf, ng);
        
        if (nH != 0)
            printf("nH:%d,", nH);
        
        printf("nR:%d,", nR);
        
        if (nV != 0)
            printf("nV(nVp):%d(%d),", nV, nVp);
        
		printf("\n");
	};

    void SolversSM::PrintFinalInfo(void)
    {
        printf("i:%d,f:%.3e,", iter, f2);

        printf("|gf|:%.3e,|gf|/|gf0|:%.3e,time:%.2g,", ngf2, ngf2/ngf0, static_cast<realdp>(getTickCount() - starttime) / CLK_PS);
        
        printf("nf:%d,ng:%d,", nf, ng);
        
        if (nH != 0)
            printf("nH:%d,", nH);
        
        printf("nR:%d,", nR);
        
        if (nV != 0)
            printf("nV(nVp):%d(%d),", nV, nVp);
        
        printf("\n");
    };

	bool SolversSM::IsStopped(void)
	{
		if (static_cast<realdp>(getTickCount() - starttime) / CLK_PS > TimeBound)
			return true;
        
		if (StopPtr != nullptr)
            return StopPtr(x2, funSeries, lengthSeries, ngf2, ngf0, Prob, this);

		if (Stop_Criterion == SM_FUN_REL)
			return ((fabs((f1 - f2) / (fabs(f1) + 1)) < Tolerance) && iter > 0);
        
		if (Stop_Criterion == SM_GRAD_F)
			return ngf2 < Tolerance;
        
		if (Stop_Criterion == SM_GRAD_F_0)
			return (ngf2 / ngf0) < Tolerance;
        
		printf("Error: Stopping Criterion is not specefic!\n");
		return true;
	};

	void SolversSM::CheckParams(void)
	{
        Solvers::CheckParams();
		std::string STOPCRITnames[STOPCRITSMLENGTH] = { "SM_FUN_REL", "SM_GRAD_F", "SM_GRAD_F_0" };
		char YES[] = "YES";
		char NO[] = "NO";
		char *status;
		printf("SMOOTH OPTIMIZATION PARAMETERS:\n");
		status = (Stop_Criterion >= 0 && Stop_Criterion < STOPCRITSMLENGTH) ? YES : NO;
		printf("Stop_Criterion:%15s[%s],\t", STOPCRITnames[Stop_Criterion].c_str(), status);
		status = (Accuracy >= 0 && Accuracy <= 1) ? YES : NO;
		printf("Accuracy      :%15g[%s],\n", Accuracy, status);
	};

	void SolversSM::Initialization(const Problem *prob, const Variable *initialx)
	{
        SetDefaultParams();
        /*Some problem setting requires the default parameters*/
		SetProbX(prob, initialx);
	};

    void SolversSM::SetDefaultParams(void)
    {
        Solvers::SetDefaultParams();
        Stop_Criterion = SM_GRAD_F_0;
        Accuracy = 1e-6;
    };

    void SolversSM::SetProbX(const Problem *prob, const Variable *initialx)
    {
        Solvers::SetProbX(prob, initialx);
        Pgf1 = Prob->GetDomain()->GetEMPTY();
        Pgf2 = Prob->GetDomain()->GetEMPTY();
        eta1 = Prob->GetDomain()->GetEMPTY();
        eta2 = Prob->GetDomain()->GetEMPTY();
    };

	SolversSM::~SolversSM(void)
	{
	};

	void SolversSM::SetParams(PARAMSMAP params)
	{
        Solvers::SetParams(params);
		PARAMSMAP::iterator iter;
		for (iter = params.begin(); iter != params.end(); iter++)
		{
			if (iter->first == static_cast<std::string> ("Stop_Criterion"))
			{
				Stop_Criterion = static_cast<StopCritSM> (static_cast<integer> (iter->second));
			}
			else
			if (iter->first == static_cast<std::string> ("Accuracy"))
			{
				Accuracy = iter->second;
			}
		}
	};
}; /*end of ROPTLIB namespace*/

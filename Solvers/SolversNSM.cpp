
#include "Solvers/SolversNSM.h"

/*Define the namespace*/
namespace ROPTLIB{

    void SolversNSM::Run(void)
    {
        /*For partly smooth functions*/
        Solvers::Run();
        if (Verbose >= ITERRESULT)
        {
            dirSeries = Vector(1 + Max_Iteration);
        }
    };

    void SolversNSM::CheckParams(void)
    {
        Solvers::CheckParams();
        std::string STOPCRITnames[STOPCRITNNSMLENGTH] = { "NSM_FUN_REL", "NSM_DIR_F", "NSM_DIR_F_0" };
        char YES[] = "YES";
        char NO[] = "NO";
        char *status;
        printf("NONSMOOTH OPTIMIZATION PARAMETERS:\n");
        status = (Stop_Criterion >= 0 && Stop_Criterion < STOPCRITNNSMLENGTH) ? YES : NO;
        printf("Stop_Criterion:%15s[%s],\n", STOPCRITnames[Stop_Criterion].c_str(), status);
    };

    void SolversNSM::SetDefaultParams(void)
    {
        Solvers::SetDefaultParams();
        Stop_Criterion = NSM_DIR_F_0;
    };

    void SolversNSM::SetProbX(const Problem *prob, const Variable *initialx)
    {
        Solvers::SetProbX(prob, initialx);
        dir1 = Prob->GetDomain()->GetEMPTY();
    };

    void SolversNSM::SetParams(PARAMSMAP params)
    {
        Solvers::SetParams(params);
        PARAMSMAP::iterator iter;
        for (iter = params.begin(); iter != params.end(); iter++)
        {
            if (iter->first == static_cast<std::string> ("Stop_Criterion"))
            {
                Stop_Criterion = static_cast<StopCritNSM> (static_cast<integer> (iter->second));
            }
        }
    };
	void SolversNSM::PrintInfo(void)
	{
        printf("i:%d,f:%.3e,df/f:%.3e,", iter, f2, ((f1 - f2) / std::fabs(f2)));

        printf("|nd|:%.3e,time:%.2e,", ndir1, static_cast<realdp>(getTickCount() - starttime) / CLK_PS);

        printf("nf:%d,ng:%d,", nf, ng);
        
        if (nH != 0)
            printf("nH:%d,", nH);
        
        printf("nR:%d,", nR);
        
        if (nV != 0)
            printf("nV(nVp):%d(%d),", nV, nVp);
		printf("\n");
	};

    void SolversNSM::PrintFinalInfo(void)
    {
        printf("i:%d,f:%.3e,df/f:%.3e,", iter, f2, ((f1 - f2) / std::fabs(f2)));

        printf("|nd|:%.3e,|nd|/|nd0|:%.3e,time:%.2e,", ndir1, ndir1 / ndir0, static_cast<realdp>(getTickCount() - starttime) / CLK_PS);

        printf("nf:%d,ng:%d,", nf, ng);
        
        if (nH != 0)
            printf("nH:%d,", nH);
        
        printf("nR:%d,", nR);
        
        if (nV != 0)
            printf("nV(nVp):%d(%d),", nV, nVp);
        printf("\n");
        
        printf("\n");
    };

	bool SolversNSM::IsStopped(void)
	{
		if (static_cast<realdp>(getTickCount() - starttime) / CLK_PS > TimeBound)
			return true;
        
        if (StopPtr != nullptr)
            return StopPtr(x2, funSeries, lengthSeries, ndir1, ndir0, Prob, this);
        
        if (Stop_Criterion == NSM_FUN_REL)
            return ((fabs((f1 - f2) / (fabs(f1) + 1)) < Tolerance) && iter > 0);
        
        if (Stop_Criterion == NSM_DIR_F)
            return ndir1 < Tolerance;
        if (Stop_Criterion == NSM_DIR_F_0)
            return (ndir1 / ndir0) < Tolerance;
        
        printf("Error: Stopping Criterion is not specefic!\n");
        return true;
	};

	void SolversNSM::Initialization(const Problem *prob, const Variable *initialx)
	{
        SetDefaultParams();
        /*Some problem setting requires the default parameters*/
		SetProbX(prob, initialx);
	};

	SolversNSM::~SolversNSM(void)
	{
	};

}; /*end of ROPTLIB namespace*/

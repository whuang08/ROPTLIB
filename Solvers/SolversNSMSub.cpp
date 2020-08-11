
#include "Solvers/SolversNSMSub.h"

/*Define the namespace*/
namespace ROPTLIB{

    void SolversNSMSub::Run(void)
    {
        SolversNSM::Run();
        /*For subgradient-based functions*/
        Lengthgfs = Mani->GetIntrDim() + NumExtraGF;

        DeleteVectors(gfs, Lengthgfs);
        NewVectors(gfs, Lengthgfs);
        DeleteVariables(Xs, Lengthgfs);
        NewVariables(Xs, Lengthgfs);
    };

	void SolversNSMSub::PrintInfo(void)
	{
        printf("i:%d,f:%.3e,df/f:%.3e,", iter, f2, ((f1 - f2) / std::fabs(f2)));

        printf("|nd|:%.3e,time:%.2e,", ndir1, static_cast<realdp>(getTickCount() - starttime) / CLK_PS);

        if (subprobtimes != 0)
            printf("nsubprob:%d,", subprobtimes);

        printf("nf:%d,ng:%d,", nf, ng);
        
        if (nH != 0)
            printf("nH:%d,", nH);
        
        printf("nR:%d,", nR);
        
        if (nV != 0)
            printf("nV(nVp):%d(%d),", nV, nVp);
		printf("\n");
	};

    void SolversNSMSub::PrintFinalInfo(void)
    {
        printf("i:%d,f:%.3e,df/f:%.3e,", iter, f2, ((f1 - f2) / std::fabs(f2)));

        printf("|nd|:%.3e,|nd|/|nd0|:%.3e,time:%.2e,", ndir1, ndir1 / ndir0, static_cast<realdp>(getTickCount() - starttime) / CLK_PS);

        if (subprobtimes != 0)
            printf("nsubprob:%d,", subprobtimes);

        printf("nf:%d,ng:%d,", nf, ng);
        
        if (nH != 0)
            printf("nH:%d,", nH);
        
        printf("nR:%d,", nR);
        
        if (nV != 0)
            printf("nV(nVp):%d(%d),", nV, nVp);
        printf("\n");
        
        printf("\n");
    };

	void SolversNSMSub::CheckParams(void)
	{
        SolversNSM::CheckParams();
		char YES[] = "YES";
		char NO[] = "NO";
		char *status;
		printf("SUBGRADIENT-BASED METHODS PARAMETERS:\n");
		status = (Diffx > 0) ? YES : NO;
		printf("Diffx         :%15g[%s],\t", Diffx, status);
		status = (NumExtraGF > 0) ? YES : NO;
		printf("NumExtraGF    :%15d[%s]\n", NumExtraGF, status);
	};

	void SolversNSMSub::Initialization(const Problem *prob, const Variable *initialx)
	{
        SetDefaultParams();
        /*Some problem setting requires the default parameters*/
		SetProbX(prob, initialx);
	};

    void SolversNSMSub::SetDefaultParams()
    {
        SolversNSM::SetDefaultParams();
        Currentlengthgfs = 0;
        idxgfs = 0;
        Diffx = static_cast<realdp> (1e-6);
        subprobtimes = 0;
        Hv = &SolversNSMSub::HvSub;
    };

	void SolversNSMSub::SetProbX(const Problem *prob, const Variable *initialx)
	{
        SolversNSM::SetProbX(prob, initialx);
        Xs = nullptr;
        gfs = nullptr;
        NumExtraGF = Mani->GetIntrDim();
	};

    Vector &SolversNSMSub::HvSub(const Vector &v, Vector *result)
    {
        *result = v;
        return *result;
    };

	SolversNSMSub::~SolversNSMSub(void)
	{
        DeleteVectors(gfs, Lengthgfs);
        DeleteVariables(Xs, Lengthgfs);
	};

	void SolversNSMSub::SetParams(PARAMSMAP params)
	{
        SolversNSM::SetParams(params);
		PARAMSMAP::iterator iter;
		for (iter = params.begin(); iter != params.end(); iter++)
		{
			if (iter->first == static_cast<std::string> ("Diffx"))
			{
				Diffx = static_cast<realdp> (iter->second);
			}
			else
			if (iter->first == static_cast<std::string> ("NumExtraGF"))
			{
				NumExtraGF = static_cast<integer> (iter->second);
			}
		}
	};
}; /*end of ROPTLIB namespace*/

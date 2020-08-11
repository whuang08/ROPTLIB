
#include "Solvers/Solvers.h"

/*Define the namespace*/
namespace ROPTLIB{

    unsigned long starttime = 0;
	void Solvers::CheckParams(void)
	{
		std::string VERBOSEnames[VERBOSELENGTH] = { "NOOUTPUT", "FINALRESULT", "ITERRESULT", "DETAILED" };
		char YES[] = "YES";
		char NO[] = "NO";
		char *status;
		printf("GENERAL PARAMETERS:\n");
		status = (Tolerance > 0) ? YES : NO;
		printf("Tolerance     :%15g[%s],\t", Tolerance, status);
        status = (TimeBound > 0) ? YES : NO;
        printf("TimeBound     :%15g[%s],\n", TimeBound, status);
        status = (OutputGap > 0) ? YES : NO;
        printf("OutputGap     :%15d[%s],\t", OutputGap, status);
		status = (Max_Iteration > 0 && Max_Iteration >= Min_Iteration) ? YES : NO;
		printf("Max_Iteration :%15d[%s],\n", Max_Iteration, status);
		status = (Min_Iteration >= 0 && Min_Iteration <= Max_Iteration) ? YES : NO;
		printf("Min_Iteration :%15d[%s],\t", Min_Iteration, status);
        status = (Verbose >= 0 && Verbose < VERBOSELENGTH) ? YES : NO;
        printf("Verbose       :%15s[%s],\n", VERBOSEnames[Verbose].c_str(), status);
	};

	void Solvers::Run(void)
	{
#ifdef MATLAB_MEX_FILE /*If ROPTLIB is used in Matlab, then start the Matlab timer.*/
		mxArray *lhs[1], *rhs[1];
		rhs[0] = mxCreateString("tic");
		mexCallMATLAB(0, lhs, 1, rhs, "feval");
#endif
		starttime = getTickCount();
		if (Verbose >= ITERRESULT)
		{
            timeSeries = Vector(1 + Max_Iteration);
            funSeries = Vector(1 + Max_Iteration);
		}
		if (Verbose >= FINALRESULT)
			printf("=========================%s=========================\n", SolverName.c_str());
	};

	void Solvers::Initialization(const Problem *prob, const Variable *initialx)
	{
        SetDefaultParams();
        /*Some problem setting requires the default parameters*/
		SetProbX(prob, initialx);
	};

    void Solvers::SetDefaultParams(void)
    {
        nf = 0; ng = 0; nV = 0; nVp = 0; nR = 0; nH = 0; lengthSeries = 0;
        StopPtr = nullptr;

        TimeBound = 60 * 60 * 24 * 365; /*one year; */
        Tolerance = static_cast<realdp> (1e-6);
        Max_Iteration = 500;
        Min_Iteration = 0;
        OutputGap = 1;
        Verbose = ITERRESULT;
    };

	void Solvers::SetProbX(const Problem *prob, const Variable *initialx)
	{
		Mani = prob->GetDomain();
		Prob = prob;
		x1 = *initialx;
		x2 = *initialx;
        gf1 = Prob->GetDomain()->GetEMPTY();
        gf2 = Prob->GetDomain()->GetEMPTY();
	};

	Solvers::~Solvers(void)
	{
	};

	void Solvers::NewVectors(Vector * &Vs, integer l)
	{
		Vs = new Vector [l];
		for (integer i = 0; i < l; i++)
            Vs[i] = Prob->GetDomain()->GetEMPTY();
	};

	void Solvers::DeleteVectors(Vector * &Vs, integer l)
	{
		if (Vs != nullptr)
			delete[] Vs;
        Vs = nullptr;
	};

	void Solvers::NewVariables(Vector * &Xs, integer l)
	{
        Xs = new Vector [l];
        for (integer i = 0; i < l; i++)
            Xs[i] = Prob->GetDomain()->GetEMPTYEXTR();
	};

	void Solvers::DeleteVariables(Vector * &Xs, integer l)
	{
		if (Xs != nullptr)
			delete[] Xs;
        Xs = nullptr;
	};

	void Solvers::SetParams(PARAMSMAP params)
	{
		PARAMSMAP::iterator iter;
		for (iter = params.begin(); iter != params.end(); iter++)
		{
			if (iter->first == static_cast<std::string> ("Tolerance"))
			{
				Tolerance = iter->second;
			}
			else
			if (iter->first == static_cast<std::string> ("TimeBound"))
			{
				TimeBound = iter->second;
			}
			else
			if (iter->first == static_cast<std::string> ("Max_Iteration"))
			{
				Max_Iteration = static_cast<integer> (iter->second);
			}
			else
			if (iter->first == static_cast<std::string> ("Min_Iteration"))
			{
				Min_Iteration = static_cast<integer> (iter->second);
			}
			else
			if (iter->first == static_cast<std::string> ("OutputGap"))
			{
				OutputGap = static_cast<integer> (iter->second);
			}
			else
			if (iter->first == static_cast<std::string> ("Verbose"))
			{
				Verbose = static_cast<VERBOSEINFO> (static_cast<integer> (iter->second));
			}
		}
	};
}; /*end of ROPTLIB namespace*/

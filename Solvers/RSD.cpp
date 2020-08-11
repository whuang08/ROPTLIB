
#include "Solvers/RSD.h"

/*Define the namespace*/
namespace ROPTLIB{

	RSD::RSD(const Problem *prob, const Variable *initialx)
	{
		Initialization(prob, initialx);
	};

	void RSD::SetProbX(const Problem *prob, const Variable *initialx)
	{
		SolversSMLS::SetProbX(prob, initialx);
		prob->SetUseGrad(true);
		prob->SetUseHess(false);
	};

	void RSD::SetDefaultParams(void)
	{
		SolversSMLS::SetDefaultParams();
		InitSteptype = LSSM_BBSTEP;
		SolverName.assign("RSD");
	};

	void RSD::GetSearchDir(void)
	{
		Prob->PreConditioner(x1, gf1, &Pgf1);
		Mani->ScalarTimesVector(x1, -1, Pgf1, &eta1);
		
		if (Mani->Metric(x1, eta1, gf1) / ngf1 / ngf1 >= -std::sqrt(std::numeric_limits<realdp>::epsilon()))
		{
			printf("Warning:Preconditioner gives an non-decreasing direction!Reset the search direction to the negative gradient\n");
			Mani->ScalarTimesVector(x1, -1, gf1, &eta1);
		}
	};
}; /*end of ROPTLIB namespace*/

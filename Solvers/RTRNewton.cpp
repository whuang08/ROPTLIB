
#include "Solvers/RTRNewton.h"

/*Define the namespace*/
namespace ROPTLIB{

	RTRNewton::RTRNewton(const Problem *prob, const Variable *initialx)
	{
		Initialization(prob, initialx);
	};

	void RTRNewton::SetProbX(const Problem *prob, const Variable *initialx)
	{
		SolversSMTR::SetProbX(prob, initialx);
		prob->SetUseGrad(true);
		prob->SetUseHess(true);
	};

	void RTRNewton::SetDefaultParams()
	{
		SolversSMTR::SetDefaultParams();
		SolverName.assign("RTRNewton");
	};

	Vector &RTRNewton::HessianEta(const Vector &Eta, Vector *result)
	{
		return Prob->HessianEta(x1, Eta, result);
	};
}; /*end of ROPTLIB namespace*/

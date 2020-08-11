
#include "Solvers/RTRSD.h"

/*Define the namespace*/
namespace ROPTLIB{

	RTRSD::RTRSD(const Problem *prob, const Variable *initialx)
	{
		Initialization(prob, initialx);
	};

	void RTRSD::SetProbX(const Problem *prob, const Variable *initialx)
	{
		SolversSMTR::SetProbX(prob, initialx);
		prob->SetUseGrad(true);
		prob->SetUseHess(false);
	};

	void RTRSD::SetDefaultParams(void)
	{
		SolversSMTR::SetDefaultParams();
		theta = static_cast<realdp> (0.1);
		kappa = static_cast<realdp> (0.9);
		SolverName.assign("RTRSD");
	};

	Vector &RTRSD::HessianEta(const Vector &Eta, Vector *result)
	{
        *result = Eta;
        return *result;
	};
}; /*end of ROPTLIB namespace*/

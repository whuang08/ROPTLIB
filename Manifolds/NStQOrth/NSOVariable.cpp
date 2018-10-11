
#include "Manifolds/NStQOrth/NSOVariable.h"

/*Define the namespace*/
namespace ROPTLIB{

	NSOVariable::NSOVariable(integer r, integer c)
	{
		Element::Initialization(2, r, c);
	};

	NSOVariable *NSOVariable::ConstructEmpty(void) const
	{
		return new NSOVariable(size[0], size[1]);
	};

	void NSOVariable::RandInManifold(void)
	{
		Element::RandGaussian();
	};
}; /*end of ROPTLIB namespace*/

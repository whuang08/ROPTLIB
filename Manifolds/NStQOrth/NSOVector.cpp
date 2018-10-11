
#include "Manifolds/NStQOrth/NSOVector.h"

/*Define the namespace*/
namespace ROPTLIB{

	NSOVector::NSOVector(integer r, integer c)
	{
		Element::Initialization(2, r, c);
	};

	NSOVector *NSOVector::ConstructEmpty(void) const
	{
		return new NSOVector(size[0], size[1]);
	};
}; /*end of ROPTLIB namespace*/

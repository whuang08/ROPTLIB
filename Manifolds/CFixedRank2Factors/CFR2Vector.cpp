
#include "Manifolds/CFixedRank2Factors/CFR2Vector.h"

/*Define the namespace*/
namespace ROPTLIB{

	CFR2Vector::CFR2Vector(integer m, integer n, integer r)
	{
		EucVector dG(2 * m, r);
		EucVector dH(2 * n, r);

		Element **Elems = new Element *[2];
		Elems[0] = &dG;
		Elems[1] = &dH;
		integer *powsintev = new integer[3];
		powsintev[0] = 0;
		powsintev[1] = 1;
		powsintev[2] = 2;

		ProductElementInitialization(Elems, 2, powsintev, 2);

		delete[] powsintev;
		delete[] Elems;
	};

	CFR2Vector::~CFR2Vector(void)
	{
	};

	CFR2Vector *CFR2Vector::ConstructEmpty(void) const
	{
		integer m = elements[0]->Getsize()[0] / 2;
		integer r = elements[0]->Getsize()[1];
		integer n = elements[1]->Getsize()[0] / 2;
		return new CFR2Vector(m, n, r);
	};
}; /*end of ROPTLIB namespace*/

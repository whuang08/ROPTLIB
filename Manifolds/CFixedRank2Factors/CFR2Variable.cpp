
#include "Manifolds/CFixedRank2Factors/CFR2Variable.h"

/*Define the namespace*/
namespace ROPTLIB{

	CFR2Variable::CFR2Variable(integer m, integer n, integer r)
	{
		EucVariable G(2 * m, r);
		EucVariable H(2 * n, r);

		Element **Elems = new Element *[2];
		Elems[0] = &G;
		Elems[1] = &H;
		integer *powsintev = new integer[3];
		powsintev[0] = 0;
		powsintev[1] = 1;
		powsintev[2] = 2;

		ProductElementInitialization(Elems, 2, powsintev, 2);

		delete[] powsintev;
		delete[] Elems;
	};

	CFR2Variable::~CFR2Variable(void)
	{
	};

	CFR2Variable *CFR2Variable::ConstructEmpty(void) const
	{
		integer m = elements[0]->Getsize()[0] / 2;
		integer r = elements[0]->Getsize()[1];
		integer n = elements[1]->Getsize()[0] / 2;
		return new CFR2Variable(m, n, r);
	};
}; /*end of ROPTLIB namespace*/

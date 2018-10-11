/*
This file defines the class of a point on the tangent space of R_*^{n \times p} / O_p, where R_*^{n \times p} is a n by p full column rank
matrix and O_p is a p-by-p orthogonal group.

SmartSpace --> Element --> NSOVector

---- WH
*/

#ifndef NSOVECTOR_H
#define NSOVECTOR_H

#include "Manifolds/Element.h"
#include <new>
#include <iostream>
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class NSOVector : public Element{
	public:
		/*Construct an empty variable on the tangent space of R_*^{n \times p} / O_p with only size information. */
		NSOVector(integer n, integer p = 1);

		/*Create an object of NSOVector with same size as this NSOVector.*/
		virtual NSOVector *ConstructEmpty(void) const;
	};
}; /*end of ROPTLIB namespace*/

#endif // end of NSOVECTOR_H

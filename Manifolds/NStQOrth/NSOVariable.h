/*
This file defines the class of a point on the manifold R_*^{n \times p} / O_p, where R_*^{n \times p} is a n by p full column rank
matrix and O_p is a p-by-p orthogonal group.

SmartSpace --> Element --> NSOVariable

---- WH
*/

#ifndef NSOVARIABLE_H
#define NSOVARIABLE_H

#include "Manifolds/Element.h"
#include <new>
#include <iostream>
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

class NSOVariable : public Element{
	public:
		/*Construct an empty variable on the manifold C_*^{n \times p} / O_p with only size information. */
		NSOVariable(integer n, integer p = 1);

		/*Create an object of CSOVariable with same size as this CSOVariable.*/
		virtual NSOVariable *ConstructEmpty(void) const;

		/*This function randomly generates a point on the manifold.*/
		virtual void RandInManifold();
	};
}; /*end of ROPTLIB namespace*/

#endif // end of NSOVARIABLE_H

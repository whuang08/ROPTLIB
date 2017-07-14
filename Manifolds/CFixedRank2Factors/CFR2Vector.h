/*
This file defines the class of a point on the tangent space of the complex quotient manifold C_*^{m \times r} \times C_*^{n \times r} / GL(r), 
where C_*^{s \times t} denotes the complex noncompact Stiefel manifold and GL(r) denotes the complex general group. 

SmartSpace --> ProductElement --> CFR2Vector

---- WH
*/

#ifndef CFR2VECTOR_H
#define CFR2VECTOR_H

#include "Manifolds/ProductElement.h"
#include "Manifolds/Euclidean/EucVector.h"

/*Define the namespace*/
namespace ROPTLIB{

	class CFR2Vector : public ProductElement{
	public:
		/*Construct an empty vector on the C_*^{m \times r} \times C_*^{n \times r} / GL(r) with only size information.*/
		CFR2Vector(integer m, integer n, integer r);

		/*Destruct by deleting variables*/
		virtual ~CFR2Vector(void);

		/*Create an object of CFR2Vector with same size as this CFR2Vector.*/
		virtual CFR2Vector *ConstructEmpty(void) const;
	};
}; /*end of ROPTLIB namespace*/
#endif // end of OBLIQUEVECTOR_H

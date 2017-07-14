/*
This file defines the class of a point on the complex quotient manifold C_*^{m \times r} \times C_*^{n \times r} / GL(r), 
where C_*^{s \times t} denotes the complex noncompact Stiefel manifold and GL(r) denotes the complex general group. 

SmartSpace --> ProductElement --> CFR2Variable

---- WH
*/

#ifndef CFR2VARIABLE_H
#define CFR2VARIABLE_H

#include "Manifolds/ProductElement.h"
#include "Manifolds/Euclidean/EucVariable.h"

/*Define the namespace*/
namespace ROPTLIB{

	class CFR2Variable : public ProductElement{
	public:
		/*Construct an empty variable on the complex quotient manifold C_*^{m \times r} \times C_*^{n \times r} / GL(r) with only size information.*/
		CFR2Variable(integer m, integer n, integer r);

		/*Destruct by deleting all variables*/
		virtual ~CFR2Variable(void);

		/*Create an object of CFR2Variable with same size as this CFR2Variable.*/
		virtual CFR2Variable *ConstructEmpty(void) const;
	};
}; /*end of ROPTLIB namespace*/
#endif // end of LOWRANKVARIABLE_H

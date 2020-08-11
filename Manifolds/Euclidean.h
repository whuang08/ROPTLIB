/*
This file defines the class for the Eucldean space.

Manifold --> Euclidean

---- WH
*/
#ifndef EUCLIDEAN_H
#define EUCLIDEAN_H


#include "Manifolds/Manifold.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class Euclidean : public Manifold{
	public:
        /*Construct the Euclidean space*/
        Euclidean(integer r, integer c = 1, integer n = 1, const char * = "real");
        
        /*Construct the Euclidean space*/
        Euclidean(integer r, integer c, const char * type);
        
        /*Construct the Euclidean space*/
        Euclidean(integer r, const char * type);
        
		/*Delete EMPTYINTR and EMPTYEXTR*/
		virtual ~Euclidean(void);
        
        /*Randomly generate a point in the manifold.*/
        virtual Variable RandominManifold() const;

		/*Check whether all the parameters are legal or not.*/
		virtual void CheckParams(void) const;

		/*gf <-- egf*/
		virtual Vector &EucGradToGrad(const Variable &x, const Vector &egf, const Problem *prob, Vector *result) const;

		/*xix <-- exix*/
		virtual Vector &EucHvToHv(const Variable &x, const Vector &etax, const Vector &exix, const Problem *prob, Vector *result) const;

		integer row; /*The first dimension of the space, i.e., the number of rows */
		integer col; /*The second dimension of the space, i.e., the number of columns */
		integer num; /*The third dimension of the space*/
        
        bool iscomplex;
	};
}; /*end of ROPTLIB namespace*/
#endif /* end of EUCLIDEAN_H */

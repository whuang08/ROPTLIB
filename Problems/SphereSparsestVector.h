/*
This file defines the class for 
min_{X \in S^{n-1}} \|Q x\|_1, where Q is a m by n matrix.
This cost function is partly smooth. The task is to find the
sparest vector in the space spanning by the columns of Q.

Problem --> SphereSparsestVector

---- WH
*/

#ifndef SPHERESPARSESTVECTOR_H
#define SPHERESPARSESTVECTOR_H

#include "Manifolds/Sphere.h"
#include "Problems/Problem.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class SphereSparsestVector : public Problem{
	public:
		SphereSparsestVector(Vector Q);
		virtual ~SphereSparsestVector(void);
		virtual realdp f(const Variable &x) const;

		virtual Vector &EucGrad(const Variable &x, Vector *result) const;
		virtual Vector &EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const;

		Vector Q;
		mutable integer m;
		mutable integer n;
	};

}; /*end of ROPTLIB namespace*/
#endif /* end of SPHERESPARSESTVECTOR_H */

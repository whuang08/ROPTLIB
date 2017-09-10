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

#include "Manifolds/Sphere/Sphere.h"
#include "Manifolds/Sphere/SphereVariable.h"
#include "Manifolds/Sphere/SphereVector.h"
#include "Problems/Problem.h"
#include "Manifolds/SharedSpace.h"
#include "Others/def.h"
#include "Others/MyMatrix.h"

/*Define the namespace*/
namespace ROPTLIB{

	class SphereSparsestVector : public Problem{
	public:
		SphereSparsestVector(double *Q, integer m, integer n);
		virtual ~SphereSparsestVector();
		virtual double f(Variable *x) const;

		virtual void EucGrad(Variable *x, Vector *egf) const;
		virtual void EucHessianEta(Variable *x, Vector *etax, Vector *exix) const;

		double *Q;
		mutable integer m;
		mutable integer n;
	};

}; /*end of ROPTLIB namespace*/
#endif // end of SPHERESPARSESTVECTOR_H

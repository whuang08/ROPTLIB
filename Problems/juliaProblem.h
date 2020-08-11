/*
This file defines the class for problem from julia function handles. This is the interface for julia users.

Problem

---- WH
*/

#ifndef JULIAPROBLEM_H
#define JULIAPROBLEM_H

#include "Problems/Problem.h"
#include <cstring>
#include "Others/def.h"

#ifdef DRIVERJULIAPROB

#include "julia.h"

/*Define the namespace*/
namespace ROPTLIB{

    class juliaProblem : public Problem{
	public:
        /*Construct a julia Problem with function handles of cost function, Euclidean gradient and action of
        Euclidean of Hessian from julia*/
        juliaProblem(jl_function_t *inf, jl_function_t *ingf, jl_function_t *inHess);

		/*Destructor*/
        virtual ~juliaProblem();

        /*call the Julia function handle of the cost function*/
		virtual double f(const Variable &x) const;

        /*call the Julia function handle of the Euclidean gradient*/
		virtual Vector &EucGrad(const Variable &x, Vector *result) const;

        /*call the Julia function handle of the action of the Euclidean Hessian*/
		virtual Vector &EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const;

	protected:
        jl_function_t *jl_f; /*Julia function handle of the cost function*/
        jl_function_t *jl_gf; /*Julia function handle of the Euclidean gradient*/
        jl_function_t *jl_Hess; /*Julia function handle of the action of the Euclidean Hessian.*/
	};
}; /*end of ROPTLIB namespace*/


#endif // end of DRIVERJULIAPROB

#endif // end of MEXPROBLEM_H

/*
This file defines the class for problem from matlab handles. This is the interface for Matlab users.

Problem

---- WH
*/

#ifndef MEXPROBLEM_H
#define MEXPROBLEM_H

#include "Problems/Problem.h"
#include <cstring>
#include "Others/def.h"

#ifdef MATLAB_MEX_FILE
#include "matrix.h"

/*Define the namespace*/
namespace ROPTLIB{

	class mexProblem : public Problem{
	public:
		/*Construct a mex Problem with function handles of cost function, Euclidean gradient and action of
		Euclidean of Hessian from Matlab*/
		mexProblem(const mxArray *inf, const mxArray *ingf, const mxArray *inHess, const mxArray *inPreCon);

		/*Destructor*/
		virtual ~mexProblem(void);

		/*call the Matlab function handle of the cost function*/
		virtual realdp f(const Variable &x) const;

		/*call the Matlab function handle of the Euclidean gradient*/
		virtual Vector &EucGrad(const Variable &x, Vector *result) const;

		/*call the Matlab function handle of the action of the Euclidean Hessian*/
		virtual Vector &EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const;

		/*The preconditioner for this problem*/
		virtual Vector &PreConditioner(const Variable &x, const Vector &eta, Vector *result) const;

		/*This function converts the format of storage in ROPTLIB to the format
		of storage in Matlab*/
		static void ObtainMxArrayFromElement(mxArray *&Xmx, const Element *X);

		/*This function converts the format of storage in Matlab to the format
		of storage in ROPTLIB*/
		static void ObtainElementFromMxArray(Element *X, const mxArray *Xmx);

		/*This function add the values of fields in the format of storage in Matlab to the format
		of storage in ROPTLIB*/
		static void AddToElementFromMxArray(Element *X, const mxArray *Xmx);

		/*S is a Matlab structure. This function obtain its field by key = name */
		static mxArray *GetFieldbyName(const mxArray *S, integer idxstruct, const char *name);


	protected:
		const mxArray *mxf; /*Matlab function handle of the cost function*/
		const mxArray *mxgf; /*Matlab function handle of the Euclidean gradient*/
		const mxArray *mxHess; /*Matlab function handle of the action of the Euclidean Hessian.*/
		const mxArray *mxPreCon; /*Matlab function handle of the Preconditioner.*/
	};
}; /*end of ROPTLIB namespace*/
#endif /* end of MATLAB_MEX_FILE */

#endif /* end of MEXPROBLEM_H */

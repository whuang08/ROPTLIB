/*
This file defines the abstract base class for all the subgradient-based solvers for nonsmooth objectives

Solvers --> SolversNSM --> SolversNSMSub

---- WH
*/

#ifndef SOLVERSNSMSUB_H
#define SOLVERSNSMSUB_H

#include <iostream>
#include <iomanip>
#include <ctime>
#include "Manifolds/Manifold.h"
#include "Problems/Problem.h"
#include "Solvers/SolversNSM.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class SolversNSMSub : public SolversNSM{
	public:
        /*Run the algorithm. This function initializes the variables for all the subgradient-based methods*/
        virtual void Run(void);

		/*Check whether the general parameters are legal or not.*/
		virtual void CheckParams(void);

		/*PARAMSMAP is defined in "def.h" and it is a map from string to realdp, i.e., std::map<std::string, realdp> .
		This function is used to set the parameters by the mapping*/
		virtual void SetParams(PARAMSMAP params);

		/*Destructor. It is a pure virtual function*/
		virtual ~SolversNSMSub(void) = 0;
        
        /*Compute result = H v in eps-subgradient-based methods*/
        virtual Vector &HvSub(const Vector &v, Vector *result);
        
		/*This parameter is used in the stopping criterion for optimizing partly smooth functions.
		If ||x_{k+1} - x_k|| / (||x_{k+1}|| + 1) is less than Diffx, then store the gradient of the latter one.
		Compute the minimum norm vector in the convex hull of the gradients stored. If the norm of the minimum norm
		vector is less than Tolerance, then algorithm terminates.
		Default: 1e-6*/
		realdp Diffx;

		/*In eps-subgradient-based algorithms, one needs to compute a minimum length in a convex hull of a few vectors. The number
        of those vectors need be greater than the dimenion of the domain manifold. This parameter specifies the number greater than
        the dimension.
		Default: Intrinsic dimension of the domain manifold*/
		integer NumExtraGF;

	protected:
		/*Initialize the solvers by calling the "SetProbX" and "SetDefultParams" functions.
			INPUT:	prob is the problem which defines the cost function, gradient and possible the action of Hessian
			and specifies the manifold of domain.
			initialx is the initial iterate. */
		virtual void Initialization(const Problem *prob, const Variable *initialx);

		/*Initialize the type of iterates x1, x2 and tangent vectors gf1, gf2 and obtian the problem and manifold information
			INPUT:	prob is the problem which defines the cost function, gradient and possible the action of Hessian
			and specifies the manifold of domain.
			initialx is the initial iterate. */
		virtual void SetProbX(const Problem *prob, const Variable *initialx);

		/*Setting parameters (member variables) to be default values */
		virtual void SetDefaultParams(void);

        /*Print information in every few iterations specific to an algorithm*/
        virtual void PrintInfo(void);
        
        /*Print last information in an algorithm*/
        virtual void PrintFinalInfo(void);
        
        /*A function pointer to an action of the Hessian*/
        Vector &(SolversNSMSub::*Hv)(const Vector &v, Vector *result);

		/* the next six parameters are for stopping criterion of optimizing partly smooth cost functions */
		Vector *gfs; /*The gradients of previous few steps*/
		integer Lengthgfs; /*The maximum length of gfs*/
		integer Currentlengthgfs; /*the current of gfs*/
		integer idxgfs; /*current idx of gfs*/
		integer subprobtimes; /*the number of runs for solving the sub convex problem*/

		/*points in a neighborhood of current iterate. It is used in RGS.*/
		Variable *Xs;
	};

}; /*end of ROPTLIB namespace*/

#endif /* end of SOLVERSNSM_H */

/*
This file defines the abstract base class for all the solvers for nonsmooth objectives
It defines the common properties and features of all the solvers for nonsmooth objectives

Solvers --> SolversNSM

---- WH
*/

#ifndef SOLVERSNSM_H
#define SOLVERSNSM_H

#include <iostream>
#include <iomanip>
#include <ctime>
#include "Manifolds/Manifold.h"
#include "Problems/Problem.h"
#include "Solvers/Solvers.h"
#include "Others/def.h"

/*The algorithm is stopped when a value (specified by ther parameter) is less than the "Tolerance" (a member variable)
The value should be assigned to the member variable: "Stop_Criterion" and the applicable values are
NSM_FUN_REL: |f_k - f_{k+1}| / max(|f_k|, 1)
NSM_DIR_F: \|gf_k\|
NSM_DIR_F_0: \|gf_k\| / \|gf_0\|*/
enum StopCritNSM{ NSM_FUN_REL, NSM_DIR_F, NSM_DIR_F_0, STOPCRITNNSMLENGTH };

/*Define the namespace*/
namespace ROPTLIB{

	/*Compute min_{y in convex hull of gfs and prefgs are tangent vectors at the tangent space at x} ||y||
	It is defined in MinPNormConHull.h and MinPNormConHull.cpp */
	extern realdp MinPNormConHull(const Manifold *Mani, Variable *x, Vector **Ys, integer LYs, Vector *Soln, realdp *YtY, integer inc);

	class SolversNSM : public Solvers{
	public:
        /*Run the algorithm. This function initializes the variables for all the non-smooth optimization methods*/
        virtual void Run(void);

        /*Check whether the general parameters are legal or not.*/
        virtual void CheckParams(void);

		/*Get the norm of the final search directionn*/
		inline realdp Getnormnd(void) const { return ndir1; };

		/*Get the norm of the search directionn at final iterate over the norm of the search directionn at initiate iterate*/
		inline realdp Getnormndnd0(void) const { return ndir1 / ndir0; };

		/* gradSeries is an array to store the norm of search directionn after each iteration*/
		inline Vector GetdirSeries(void) const { return dirSeries; };

		/*Destructor. It is a pure virtual function*/
		virtual ~SolversNSM(void) = 0;
        
        /*member variable: specify the stopping criterion. The applicable values are defined in enumerate "StopCrit"
        Default: NSM_DIR_F_0*/
        StopCritNSM Stop_Criterion;
        
        /*PARAMSMAP is defined in "def.h" and it is a map from string to realdp, i.e., std::map<std::string, realdp> .
        This function is used to set the parameters by the mapping*/
        virtual void SetParams(PARAMSMAP params);

	protected:
		/*Initialize the solvers by calling the "SetProbX" and "SetDefultParams" functions.
			INPUT:	prob is the problem which defines the cost function, gradient and possible the action of Hessian
			and specifies the manifold of domain.
			initialx is the initial iterate. */
		virtual void Initialization(const Problem *prob, const Variable *initialx);
        
        /*Setting parameters (member variables) to be default values */
        virtual void SetDefaultParams(void);
        
        /*Initialize the type of iterates x1, x2 and tangent vectors gf1, gf2 and obtian the problem and manifold information
            INPUT:    prob is the problem which defines the cost function, gradient and possible the action of Hessian
            and specifies the manifold of domain.
            initialx is the initial iterate. */
        virtual void SetProbX(const Problem *prob, const Variable *initialx);

		/*When one iteration, some algorithms need to update some information. For example,
			quasi-Newton methods need to update the Hessian approximation. They are done in the following function*/
		virtual void UpdateData(void) = 0;

		/*Check whether a stopping criterion is satisfied or not*/
		virtual bool IsStopped(void);
        
        /*Print information in every few iterations specific to an algorithm*/
        virtual void PrintInfo(void);
        
        /*Print last information in an algorithm*/
        virtual void PrintFinalInfo(void);

		/* algorithm-related variables: */
		Vector dir1;	/*dir1: search directionn at x1;*/
		realdp ndir0, ndir1; /*ndir0: the norm of the search direction at the initial iterate, ndir1: the norm of the search direction at x1;*/

		/* For debug information */
		Vector dirSeries; /*an array to store the norm of search direction after each iteration*/
	};

}; /*end of ROPTLIB namespace*/

#endif /* end of SOLVERSNSM_H */

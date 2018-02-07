/*
This file defines the class of the limited-memory Riemannian trust-region symmetric rank one update method in [HAG2014]
[HAG2014]: W. Huang, P.-A. Absil, and K. A. Gallivan. A Riemannian symmetric rank-one trustregion method. 
		Mathematical Programming, 150(2):179?16, February 2015.

Solvers --> SolversTR --> LRTRSR1

---- WH
*/

#ifndef LRTRSR1_H
#define LRTRSR1_H

#include "Solvers/SolversTR.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class LRTRSR1 : public SolversTR{
	public:
		/*The contructor of LRTRSR1 method. It calls the function Solvers::Initialization.
		INPUT : prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.
		insoln is the true solution. It is not required and only used for research.*/
		LRTRSR1(const Problem *prob, const Variable *initialx, const Variable *insoln = nullptr);

		/*Destructor. Delete the arrays and vectors used in LRTRSR1,
		i.e., vector series S and Y and YMGS, double series RHO, matrices PMGQ, SS, and YY, and permutation indices P (see below for details of the series)*/
		virtual ~LRTRSR1();

		/*Check whether the parameters about LRBFGS are legal or not.*/
		virtual void CheckParams();

		/*Run the algorithm. New memory for vector series S and Y and YMGS, double series RHO,
			matrices PMGQ, SS, and YY, and permutation indices P (see below for details of the series). Then call SolversTR::Run*/
		virtual void Run();

		/*Computes the solution to a trust-region subproblem when the quadratic model is defined by 
		limited-memory symmetric rank-one (L-SR1) quasi-Newton matrix. The details can be found in [JJR17]
		[JJR17]: On solving L-SR1 trust-region subproblems */
		virtual void tCG_TR(void);

		/*Call Solvers::SetProbX function and set up the temporary objects for LRTRSR1 algorithm.
		INPUT:	prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.
		insoln is the true solution. It is not required and only used for research.*/
		virtual void SetProbX(const Problem *prob, const Variable *initialx, const Variable *insoln);

		/*Setting parameters (member variables) to be default values */
		virtual void SetDefaultParams();

	protected:
		/*Compute result = H[Eta], where H is the Hessian or the Hessian approximation
			See details in [HAG2014, Section 4]
			[HAG2014]: W. Huang, P.-A. Absil, and K. A. Gallivan. A Riemannian symmetric rank-one trustregion method.
			Mathematical Programming, 150(2):179?16, February 2015.*/
		virtual void HessianEta(Vector *Eta, Vector *result);

		/*update the pairs of s and y. Add the latest one and remove the oldest one if necessary.*/
		virtual void UpdateData();

		/*This function is called when the candidate is accepted. It transports tangent vectors paris of s and y
		from the tangent space at x1 to the tangent space at x2.*/
		virtual void Acceptence();

		/*Print information*/
		virtual void PrintGenInfo(void);

		/*Print information specific to LRTRSR1*/
		virtual void PrintInfo();

		void ComputeSBySMW(double tauStar, double *Psig, double *RT);

		double PhiBar_fg(double nlambda_min, double Delta, double *Lambda, double *a_j, double *gf = nullptr);

		double LocalNewton(double x0, integer maxIter, double tol, double *Lambda, double *a_j);
	};
}; /*end of ROPTLIB namespace*/
#endif // end of LRTRSR1_H

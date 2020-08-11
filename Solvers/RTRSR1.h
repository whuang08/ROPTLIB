/*
This file defines the class of the Riemannian trust-region Newton method in [HAG2015]
	[HAG2015]: W. Huang, P.-A. Absil, and K. A. Gallivan. A Riemannian symmetric rank-one trust-region method. 
	Mathematical Programming, 150(2):179?16, 2015

Solvers --> SolversSM --> SolversSMTR --> RTRSR1

---- WH
*/

#ifndef RTRSR1_H
#define RTRSR1_H

#include "Solvers/SolversSMTR.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class RTRSR1 : public SolversSMTR{
	public:
		/*The contructor of RTRSR1 method. It calls the function Solvers::Initialization.
		INPUT : prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.
		initialB is the initial Hessian approximation. If the input is nullptr, then the identity is used as the initial approximation.*/
		RTRSR1(const Problem *prob, const Variable *initialx, LinearOPE *initialB = nullptr);

		/*Destructor.*/
		virtual ~RTRSR1(void);

	protected:
		/*Compute the action of the Hessian approximation. result = B [Eta]*/
		virtual Vector &HessianEta(const Vector &Eta, Vector *result);

		/*Update the Hessian approximation if necessary*/
		virtual void UpdateData(void);

		/*This function is called when the candidate is accepted. It transports the Hessian approximation B
		from the tangent space of x1 to x2*/
		virtual void Acceptence(void);

		/*Print information specific to RTRSR1*/
		virtual void PrintInfo(void);
        
        /*Initialize the solvers by calling the "SetProbX" and "SetDefultParams" functions.
        INPUT:    prob is the problem which defines the cost function, gradient and possible the action of Hessian
        and specifies the manifold of domain.
        initialx is the initial iterate.
        initialH is the initial inverse Hessian approximation. If the input is nullptr, then the identity is used as the initial approximation.
        insoln is the true solution. It is not required and only used for research.*/
        virtual void Initialization(const Problem *prob, const Variable *initialx, LinearOPE *initialB = nullptr);

        /*Initialize the type of iterates x1, x2 and tangent vectors gf1, gf2, s, y, H and tildeH and obtian the problem and manifold information
        INPUT:    prob is the problem which defines the cost function, gradient and possible the action of Hessian
        and specifies the manifold of domain.
        initialx is the initial iterate.
        initialB is the initial Hessian approximation. If the input is nullptr, then the identity is used as the initial approximation.
        insoln is the true solution. It is not required and only used for research.*/
        virtual void SetProbX(const Problem *prob, const Variable *initialx, LinearOPE *initialB = nullptr);

        /*Setting parameters (member variables) to be default values */
        virtual void SetDefaultParams(void);
        
        /*Compute result = H v in RTRSR1*/
        virtual Vector &HvRTRSR1(const Vector &v, const LinearOPE &B, Vector *result);
        
        /*Update the Hessian approximation for RTRSR1 if necessary*/
        virtual void UpdateDataRTRSR1(void);

        bool isupdated; /*Mark whether the (inverse) Hessian approximation is updated*/
        realdp inpss;  /*inpsy: g(s, y); inpss: g(s, s); inpyy: g(y, y); */

        Vector s, y;/*the s, y, and u of current step*/
        LinearOPE B; /*The Hessian approximations for current and next iterations respectively*/
	};
}; /*end of ROPTLIB namespace*/

#endif /* end of RTRSR1_H */

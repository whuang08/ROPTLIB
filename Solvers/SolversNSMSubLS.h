/*
This file defines the class of the line search algorithm for locally lipschitz functions on Riemannian manifolds

Solvers --> SolversNSM --> SolversNSMSub --> SolversNSMSubLS

---- WH
*/

#ifndef SOLVERSNNSMSUBLS_H
#define SOLVERSNNSMSUBLS_H

#include "Others/MinPNormConHull.h"
#include "Solvers/SolversNSMSub.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

    /* Linesearch status. It is an output argument and users don't need to assign this enumerate to any member variable.
    NOCURVATURE: the second Wolfe condition is not satisfied
    MINSTEPSIZE: line search algorithm reaches the minimum stepsize
    MAXSTEPSIZE: line search algorithm reaches the maximum stepsize
    SUCCESS: line search algorithm succeeds in finding a point satisfying the line search condition
    */
    enum LSstatusSetNSM{ LSNSM_NOCURVATURE, LSNSM_MINSTEPSIZE, LSNSM_MAXSTEPSIZE, LSNSM_LSERROR, LSNSM_SUCCESS, LSSTATUSSETNSMLENGTH };

	class SolversNSMSubLS : public SolversNSMSub{
	public:
		/*Run the algorithm. This function gives the framework for the linesearch method*/
		virtual void Run(void);

		/*Destructor. Delete the vectors and Hessian approximation used in RBFGSLPSub, i.e., s and y, H and tildeH*/
		virtual ~SolversNSMSubLS(void);

		/*PARAMSMAP is defined in "def.h" and it is a map from string to realdp, i.e., std::map<std::string, realdp> .
		This function is used to set the parameters by the mapping*/
		virtual void SetParams(PARAMSMAP params);

		/*Check whether the parameters about RBFGSLPSub are legal or not.*/
		virtual void CheckParams(void);

		/*Initialize the type of iterates x1, x2 and tangent vectors gf and obtian the problem and manifold information
		INPUT:	prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.*/
		virtual void SetProbX(const Problem *prob, const Variable *initialx);

		/*Setting parameters (member variables) to be default values */
		virtual void SetDefaultParams(void);

		/*Epsilon-subdifferential
		Default: 1*/
		realdp Eps;
		/*Reduce the Epsilon by Theta_eps each time, i.e., Eps *= Theta_eps.
		Default: 0.01*/
		realdp Theta_eps;
		/*The minimum value of Eps. Too small epsilon have non-negligible numerical error.
		Default: 1e-6*/
		realdp Min_Eps;

		/*The bound on the square of P-norm of gradient
		Default: 1*/
		realdp Del;
		/*Reduce Del by Theta_del each time.
		Default: 0.01*/
		realdp Theta_del;
        
        /*the coefficient of the Wolfe first condition
        Default: 0.0001 */
        realdp LS_alpha;

        /*the coefficient of the Wolfe second condition
        Default: 0.999*/
        realdp LS_beta;
        
        /*the minimum stepsize allowed in the linesearch algorithm
        Default: machine eps*/
        realdp Minstepsize;

        /*the maximum stepsize allowed in the linesearch algorithm
        Default: 1000 */
        realdp Maxstepsize;
        
        /*lower bound of the smallest eigenvalue of the hessian approximation
        Default: 1e-7*/
        realdp lambdaLower;
        
        /*lower bound of the largest eigenvalue of the hessian approximation
        Default: 1e7*/
        realdp lambdaUpper;
        
        /*Compute result = H v in eps-subgradient-based methods*/
        virtual Vector &HvSub(const Vector &v, Vector *result);
        
	protected:
        /*Print information in every few iterations specific to an algorithm*/
        virtual void PrintInfo(void);
        
        /*Print last information in an algorithm*/
        virtual void PrintFinalInfo(void);
        
        /*Evaluate the cost function h(stepsize) = f(R_{x_1}(stepsize * eta1))*/
        virtual realdp h(void);
        
        /*Evaluate the derivative of cost function h, i.e., h'(stepsize) = \frac{d}{d stepsize} f(R_{x_1}(stepsize * eta1))*/
        virtual realdp dh(void);

		void Increasing(realdp neta1, realdp Pngfsq, realdp hb);

		/*Compute the search direction */
		virtual void GetSearchDir(void);
        
        /*The weak Wolfe condition for Lipschitz functions*/
        virtual void LinesearchWolfeLipschitz(void);
     
		/*Algorithm-related variables*/
		Vector minPv; /*the subgradient with minimum P norm in eps-subgradient of f*/
        Vector eta1;
        Vector eta2;
        
        LSstatusSetNSM LSstatus;
        realdp initiallength;
        realdp stepsize;        /*The step size*/
        realdp initialslope;
        realdp newslope;
        std::string *LSstatusSetnames;    /*This string array is to store the line search status names*/
    private:
        void ZoomLP(realdp a, realdp b);

	};
}; /*end of ROPTLIB namespace*/
#endif /* end of SOLVERSNNSMSUBLS_H */


#include "Problems/Problem.h"

/*Define the namespace*/
namespace ROPTLIB{

	void Problem::CheckGradHessian(Variable xin) const
	{
		UseGrad = true;
		UseHess = true;
		integer length;
		realdp normxi;
		realdp t, fx, fy;
		realdp *X, *Y;
        
        Variable x = xin, y = x;
        fx = f(x);
        Vector gfx = Domain->GetEMPTY(); Grad(x, &gfx);
        Vector etax = gfx;
        Vector xi = Domain->GetEMPTY(), Hv(xi);
        printf("f:%f\n", fx);
        Domain->Projection(x, etax, &xi);
        normxi = sqrt(Domain->Metric(x, xi, xi));
        Domain->ScalarTimesVector(x, 100 / normxi, xi, &xi); /* initial length of xi is 100 */
        /* the length of xi variances from 100 to 100*2^(-35) approx 6e-9 */
        t = 1;
        length = 35;
        X = new realdp[length * 2];
        Y = X + length;

        for (integer i = 0; i < length; i++)
        {
            Domain->Retraction(x, xi, &y);
            fy = f(y);
            HessianEta(x, xi, &Hv);
            Y[i] = log(fabs(fy - fx - Domain->Metric(x, gfx, xi) - static_cast<realdp> (0.5) * Domain->Metric(x, xi, Hv)));
            X[i] = static_cast<realdp> (0.5) * log(Domain->Metric(x, xi, xi));
            printf("i:%d,|eta|:%.3e,(fy-fx)/<gfx,eta>:%.3e,(fy-fx-<gfx,eta>)/<0.5 eta, Hessian eta>:%.3e\n", i,
                sqrt(Domain->Metric(x, xi, xi)), (fy - fx) / Domain->Metric(x, gfx, xi),
                (fy - fx - Domain->Metric(x, gfx, xi)) / (0.5 * Domain->Metric(x, xi, Hv)));
            Domain->ScalarTimesVector(x, 0.5, xi, &xi);
        }
        delete[] X;

        printf("CHECK GRADIENT:\n");
        printf("\tSuppose the point is not a critical point.\n");
        printf("\tIf there exists an interval of |eta| such that (fy - fx) / <gfx, eta>\n");
        printf("\tapproximates ONE, then the gradient is probably correct!\n");

        printf("CHECK THE ACTION OF THE HESSIAN (PRESUME GRADIENT IS CORRECT):\n");
        printf("\tSuppose the retraction is second order or the point is a critical point.\n");
        printf("\tIf there exists an interval of |eta| such that (fy-fx-<gfx,eta>)/<0.5 eta, Hessian eta>\n");
        printf("\tapproximates ONE, then the action of Hessian is probably correct.\n");
	};

    Vector Problem::MinMaxEigValHess(Variable x) const
    {
        return MinMaxEigValHessian(&x, Domain, this);
    };

	Vector &Problem::Grad(const Variable &x, Vector *result) const
	{
		if (!Domain->GetIsIntrinsic())
			return RieGrad(x, result);
        
        Vector ExGrad(Domain->GetEMPTYEXTR());
        RieGrad(x, &ExGrad);
        Domain->ObtainIntr(x, ExGrad, result);
        return *result;
	};

	Vector &Problem::HessianEta(const Variable &x, const Vector &etax, Vector *result) const
	{
		if (!Domain->GetIsIntrinsic())
			return RieHessianEta(x, etax, result);
        
        Vector ExHeta(Domain->GetEMPTYEXTR()), Exetax(Domain->GetEMPTYEXTR());
        Domain->ObtainExtr(x, etax, &Exetax);
        RieHessianEta(x, Exetax, &ExHeta);
        return Domain->ObtainIntr(x, ExHeta, result);
	};

	Vector &Problem::RieGrad(const Variable &x, Vector *result) const
	{
        Vector EGrad(*result);
        
        if(NumGradHess)
        {
            Domain->EucGradToGrad(x, Problem::EucGrad(x, &EGrad), this, result);
            return *result;
        }
        
        Domain->EucGradToGrad(x, EucGrad(x, &EGrad), this, result);
        return *result;
	};

	Vector &Problem::RieHessianEta(const Variable &x, const Vector &etax, Vector *result) const
	{
        Vector EHeta(etax);
        
        if(NumGradHess)
            return Domain->EucHvToHv(x, etax, Problem::EucHessianEta(x, etax, &EHeta), this, result);
        
        return Domain->EucHvToHv(x, etax, EucHessianEta(x, etax, &EHeta), this, result);
	};

	Vector &Problem::EucGrad(const Variable &x, Vector *result) const
    {/*Compute Numerical Gradient: by Sean Martin, modified by WH*/
        
#ifdef SINGLE_PRECISION
        realdp _eps = 1e-4;
#else
        realdp _eps = 1e-8;
#endif
        
        realdp fx = f(x);
        size_t nn = x.Getlength();
        Variable x_eps = x;
        
        const realdp *x_ptr = x.ObtainReadData();
        realdp *x_eps_ptr = x_eps.ObtainWriteEntireData();
        Vector egf(x);
        realdp *egf_ptr = result->ObtainWriteEntireData();
        
        for (size_t i = 0; i < nn; ++i) {
            x_eps_ptr[i] = x_ptr[i];
        }
        
        for (size_t i = 0; i < nn; ++i) {
            x_eps_ptr[i] += _eps;
            x_eps.RemoveAllFromFields();
            double fp = f(x_eps);
            egf_ptr[i] = (fp - fx) / _eps;
            x_eps_ptr[i] = x_ptr[i];
        }
        return *result;
	};

	Vector &Problem::EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const
	{ /*By finite difference*/
		realdp normetax = sqrt(Domain->Metric(x, etax, etax));
#ifdef SINGLE_PRECISION
        realdp factor = static_cast<realdp> (5e-2) / normetax;
#else
        realdp factor = static_cast<realdp> (1e-5) / normetax;
#endif
        Vector exix = Domain->GetEMPTYEXTR(), y(x);
		Domain->ScalarTimesVector(x, factor, etax, &exix); /*exix = alpha * etax */
        y.AlphaXaddThis(1, exix); /*y = x + exix*/
		/*
		f(y) is evaluated before EucGrad since some computations, which are needed in EucGrad, are done in f(y).
		The Euclidean gradient uses extrinsic approach.
		*/
        f(y); Vector gfy(Domain->GetEMPTYEXTR()); EucGrad(y, &gfy);
        
        /*In the function EucGradToGrad, The Euclidean gradient has been added to x with field name: "EGrad" */
        Vector gfx = x.Field("EGrad");
		Domain->VectorLinearCombination(x, static_cast<realdp> (1) / factor, gfy, static_cast<realdp> (-1) / factor, gfx, result);
        return *result;
	};

	Vector &Problem::ProxW(const Vector &x, const Vector &Weight, Vector *result) const
	{
		/* default one without nonsmooth term, i.e., lambda = 0 */
        printf("Warning: Problem::ProxW has not been overridden! It may not be correct!\n");
        *result = x;
        return *result;
	};

	Vector &Problem::CalJW(const Vector &x, const Vector &eta, const Vector &Weight, Vector *result) const
	{
		/* default one without nonsmooth term, i.e., lambda = 0 */
        printf("Warning: Problem::CalJW has not been overridden! It may not be correct!\n");
        *result = eta;
        return *result;
	};

	Vector &Problem::PreConditioner(const Variable &x, const Vector &eta, Vector *result) const
	{
		/* default one means no preconditioner. */
        *result = eta;
        return *result;
	};

	void Problem::SetDomain(Manifold *inDomain)
	{
		Domain = inDomain;
	};

	Problem::~Problem(void)
	{
	};

	void Problem::SetUseGrad(bool usegrad) const
	{
		UseGrad = usegrad;
	};

    void Problem::SetUseHess(bool usehess) const
    {
        UseHess = usehess;
    };

    void Problem::SetNumGradHess(bool inNumGradHess) const
    {
        NumGradHess = inNumGradHess;
    };
}; /*end of ROPTLIB namespace*/

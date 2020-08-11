
#include "Manifolds/Euclidean.h"

/*Define the namespace*/
namespace ROPTLIB{

	Euclidean::Euclidean(integer r, integer c, integer n, const char *type)
	{
        if(strcmp(type, "real") == 0) /*type == 'real'*/
        {
            iscomplex = false;
            row = r;
            col = c;
            num = n;
            IsIntrApproach = false;
            HasHHR = false;
            name.assign("Euclidean");
            IntrinsicDim = n * r * c;
            ExtrinsicDim = n * r * c;
            Vector tmp(r, c, n);
            EMPTYEXTR = tmp;
            EMPTYINTR = tmp;
            return;
        }

        iscomplex = true;
        row = r;
        col = c;
        num = n;
        IsIntrApproach = false;
        HasHHR = false;
        name.assign("Euclidean");
        IntrinsicDim = 2 * n * r * c;
        ExtrinsicDim = 2 * n * r * c;
        Vector tmp(r, c, n, "complex");
        EMPTYEXTR = tmp;
        EMPTYINTR = tmp;
	};

    Euclidean::Euclidean(integer r, integer c, const char * type)
    {
        integer n = 1;
        if(strcmp(type, "real") == 0) /*type == 'real'*/
        {
            iscomplex = false;
            row = r;
            col = c;
            num = n;
            IsIntrApproach = false;
            HasHHR = false;
            name.assign("Euclidean");
            IntrinsicDim = n * r * c;
            ExtrinsicDim = n * r * c;
            Vector tmp(r, c, n);
            EMPTYEXTR = tmp;
            EMPTYINTR = tmp;
            return;
        }

        iscomplex = true;
        row = r;
        col = c;
        num = n;
        IsIntrApproach = false;
        HasHHR = false;
        name.assign("Euclidean");
        IntrinsicDim = 2 * n * r * c;
        ExtrinsicDim = 2 * n * r * c;
        Vector tmp(r, c, n, "complex");
        EMPTYEXTR = tmp;
        EMPTYINTR = tmp;
    };

    Euclidean::Euclidean(integer r, const char * type)
    {
        integer n = 1, c = 1;
        if(strcmp(type, "real") == 0) /*type == 'real'*/
        {
            iscomplex = false;
            row = r;
            col = c;
            num = n;
            IsIntrApproach = false;
            HasHHR = false;
            name.assign("Euclidean");
            IntrinsicDim = n * r * c;
            ExtrinsicDim = n * r * c;
            Vector tmp(r, c, n);
            EMPTYEXTR = tmp;
            EMPTYINTR = tmp;
            return;
        }

        iscomplex = true;
        row = r;
        col = c;
        num = n;
        IsIntrApproach = false;
        HasHHR = false;
        name.assign("Euclidean");
        IntrinsicDim = 2 * n * r * c;
        ExtrinsicDim = 2 * n * r * c;
        Vector tmp(r, c, n, "complex");
        EMPTYEXTR = tmp;
        EMPTYINTR = tmp;
    };

	Euclidean::~Euclidean(void)
	{
	};

    Variable Euclidean::RandominManifold(void) const
    {
        Vector result(row, col, num, iscomplex);
        result.RandGaussian();
        return result;
    };

	void Euclidean::CheckParams(void) const
	{
		Manifold::CheckParams();
		printf("%s PARAMETERS:\n", name.c_str());
		if (col == 1 && num == 1)
			printf("row           :%15d\n", row);
		else
			if (num == 1)
			{
				printf("row           :%15d,\t", row);
				printf("col           :%15d\n", col);
			}
			else
			{
				printf("row           :%15d,\t", row);
				printf("col           :%15d\n", col);
				printf("num           :%15d\n", num);
			}
	};

	Vector &Euclidean::EucGradToGrad(const Variable &x, const Vector &egf, const Problem *prob, Vector *result) const
	{
        if(prob->GetUseHess())
        {
            x.AddToFields("EGrad", egf);
        }
        *result = egf;
        return *result;
	};

	Vector &Euclidean::EucHvToHv(const Variable &x, const Vector &etax, const Vector &exix, const Problem *prob, Vector *result) const
	{
        *result = exix;
        return *result;
	};
}; /*end of ROPTLIB namespace*/

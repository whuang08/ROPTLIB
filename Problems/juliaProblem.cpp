
#include "Problems/juliaProblem.h"

#ifdef DRIVERJULIAPROB

/*Define the namespace*/
namespace ROPTLIB{

    juliaProblem::juliaProblem(jl_function_t *inf, jl_function_t *ingf, jl_function_t *inHess)
    {
        jl_f = inf;
        jl_gf = ingf;
        jl_Hess = inHess;
	};

    juliaProblem::~juliaProblem()
	{
	};

    double juliaProblem::f(const Variable &x) const
    {
//        x.Print("x in f:");//---
		jl_value_t* array_type = jl_apply_array_type((jl_value_t *) jl_float64_type, 1);
		jl_array_t *arrtmp = nullptr;
		if (x.FieldsExist("Tmp"))
		{
            const double *tmpptr = x.Field("Tmp").ObtainReadData();
			arrtmp = jl_ptr_to_array_1d(array_type, const_cast<double *> (tmpptr), x.Field("Tmp").Getlength(), 0);
		}
		else
		{
			arrtmp = jl_ptr_to_array_1d(array_type, nullptr, 0, 0);
		}
//        std::cout << "t1" << std::endl;//---
        const double *xptr = x.ObtainReadData();
//        std::cout << "t2" << std::endl;//---
        jl_array_t *arrx = jl_ptr_to_array_1d(array_type, const_cast<double *> (xptr), x.Getlength(), 0);
        
//        std::cout << "t3" << std::endl;//---

//        std::cout << (long) jl_f << std::endl;//---
        
        jl_value_t *retresult = jl_call2(jl_f, (jl_value_t *) arrx, (jl_value_t *) arrtmp);
//        std::cout << "t4" << std::endl;//---
//        std::cout << (long) retresult << std::endl;//---
        
//        jl_get_nth_field(retresult, 0);
//        std::cout << "t41" << std::endl;//---
        jl_value_t *fx = jl_get_nth_field(retresult, 0);
//        std::cout << "t5" << std::endl;//---
        jl_array_t *outtmp = (jl_array_t *) jl_get_nth_field(retresult, 1);
//        std::cout << "t6" << std::endl;//---

        integer outtmplen = jl_array_len(outtmp);
        
//        std::cout << "t7" << std::endl;//---
        Vector sharedouttmp(outtmplen);
        double *sharedouttmpptr = sharedouttmp.ObtainWriteEntireData();
//        std::cout << "t8" << std::endl;//---
        dcopy_(&outtmplen, (double*)jl_array_data(outtmp), &GLOBAL::IONE, sharedouttmpptr, &GLOBAL::IONE);
        x.AddToFields("Tmp", sharedouttmp);
        
//        std::cout << "t9" << std::endl;//---
        double result = jl_unbox_float64(fx);
//        std::cout << "t10" << std::endl;//---
        return result;
	};

    Vector &juliaProblem::EucGrad(const Variable &x, Vector *result) const
    {
		jl_value_t* array_type = jl_apply_array_type((jl_value_t *) jl_float64_type, 1);
		jl_array_t *arrtmp = nullptr;
		if (x.FieldsExist("Tmp"))
		{
            const double *tmpptr = x.Field("Tmp").ObtainReadData();
            arrtmp = jl_ptr_to_array_1d(array_type, const_cast<double *> (tmpptr), x.Field("Tmp").Getlength(), 0);
		}
		else
		{
			arrtmp = jl_ptr_to_array_1d(array_type, nullptr, 0, 0);
		}

        const double *xptr = x.ObtainReadData();
        jl_array_t *arrx = jl_ptr_to_array_1d(array_type, const_cast<double *> (xptr), x.Getlength(), 0);

        if(jl_gf == nullptr) /*use numerical gradient*/
        {
            return Problem::EucGrad(x, result);
        }
        
        jl_value_t *retresult = jl_call2(jl_gf, (jl_value_t *) arrx, (jl_value_t *) arrtmp);
        jl_array_t *jl_egf = (jl_array_t *) jl_get_nth_field(retresult, 0);
        jl_array_t *outtmp = (jl_array_t *) jl_get_nth_field(retresult, 1);

        if(jl_array_len(jl_egf) != result->Getlength())
        {
            std::cout << "error: the size of the Euclidean gradient is not correct!" << std::endl;
            exit(EXIT_FAILURE);
        }

        integer egflen = result->Getlength();
        double *egfptr = result->ObtainWriteEntireData();
        dcopy_(&egflen, (double*)jl_array_data(jl_egf), &GLOBAL::IONE, egfptr, &GLOBAL::IONE);

        integer outtmplen = jl_array_len(outtmp);
        if(outtmplen != 0)
        {
            Vector sharedouttmp(outtmplen);
            double *sharedouttmpptr = sharedouttmp.ObtainWriteEntireData();
            dcopy_(&outtmplen, (double*)jl_array_data(outtmp), &GLOBAL::IONE, sharedouttmpptr, &GLOBAL::IONE);
            x.AddToFields("Tmp", sharedouttmp);
        }
        return *result;
	};

    Vector &juliaProblem::EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const
    {
		jl_value_t* array_type = jl_apply_array_type((jl_value_t *) jl_float64_type, 1);
		jl_array_t *arrtmp = nullptr;
        if (x.FieldsExist("Tmp"))
        {
            const double *tmpptr = x.Field("Tmp").ObtainReadData();
            arrtmp = jl_ptr_to_array_1d(array_type, const_cast<double *> (tmpptr), x.Field("Tmp").Getlength(), 0);
        }
        else
        {
            arrtmp = jl_ptr_to_array_1d(array_type, nullptr, 0, 0);
        }
        
        const double *xptr = x.ObtainReadData();
        jl_array_t *arrx = jl_ptr_to_array_1d(array_type, const_cast<double *> (xptr), x.Getlength(), 0);
        const double *etaxptr = etax.ObtainReadData();
        jl_array_t *arretax = jl_ptr_to_array_1d(array_type, const_cast<double *> (etaxptr), etax.Getlength(), 0);

        if(jl_Hess == nullptr)  /*use numerical Hessian*/
        {
            return Problem::EucHessianEta(x, etax, result);
        }
        
        jl_value_t *retresult = jl_call3(jl_Hess, (jl_value_t *) arrx, (jl_value_t *) arrtmp, (jl_value_t *) arretax);
        jl_array_t *jl_exix = (jl_array_t *) jl_get_nth_field(retresult, 0);
        jl_array_t *outtmp = (jl_array_t *) jl_get_nth_field(retresult, 1);

        if(jl_array_len(jl_exix) != etax.Getlength())
        {
            std::cout << "error: the size of the action of the Hessian is not correct!" << std::endl;
            exit(EXIT_FAILURE);
        }

        integer exixlen = result->Getlength();
        double *exixptr = result->ObtainWriteEntireData();
        dcopy_(&exixlen, (double*)jl_array_data(jl_exix), &GLOBAL::IONE, exixptr, &GLOBAL::IONE);

        integer outtmplen = jl_array_len(outtmp);
        if(outtmplen != 0)
        {
            Vector sharedouttmp(outtmplen);
            double *sharedouttmpptr = sharedouttmp.ObtainWriteEntireData();
            dcopy_(&outtmplen, (double*)jl_array_data(outtmp), &GLOBAL::IONE, sharedouttmpptr, &GLOBAL::IONE);
            x.AddToFields("Tmp", sharedouttmp);
        }
        return *result;
	};

}; /*end of ROPTLIB namespace*/

#endif // end of DRIVERJULIAPROB

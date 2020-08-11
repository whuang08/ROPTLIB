
#include "Others/SparseMatrix.h"

/*Define the namespace*/
namespace ROPTLIB{

    SparseMatrix::SparseMatrix(integer inm, integer inn, integer *ir, integer *jc, realdp *vals, integer nz)
    {
        m = inm;
        n = inn;
        iscomplex = false;
#ifdef SINGLE_PRECISION
        SparseM = BLAS_suscr_begin(m, n);
#else
        SparseM = BLAS_duscr_begin(m, n);
#endif
        BLAS_uscr_insert_entries(SparseM, nz, vals, ir, jc);
#ifdef SINGLE_PRECISION
        BLAS_suscr_end(SparseM);
#else
        BLAS_duscr_end(SparseM);
#endif
    };

    SparseMatrix::SparseMatrix(integer inm, integer inn, integer *ir, integer *jc, realdpcomplex *vals, integer nz)
    {
        m = inm;
        n = inn;
        iscomplex = true;
#ifdef SINGLE_PRECISION
        SparseM = BLAS_cuscr_begin(m, n);
#else
        SparseM = BLAS_zuscr_begin(m, n);
#endif
        BLAS_uscr_insert_entries(SparseM, nz, vals, ir, jc);
#ifdef SINGLE_PRECISION
        BLAS_cuscr_end(SparseM);
#else
        BLAS_zuscr_end(SparseM);
#endif
    };

    void SparseMatrix::Print(const char *name) const
    {
        printf("%s:\n", name);
        NIST_SPBLAS::print(SparseM);
    };

    SparseMatrix::~SparseMatrix(void)
    {
        BLAS_usds(SparseM);
    };

}; /*end of ROPTLIB namespace*/

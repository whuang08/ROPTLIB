/*
This file defines the class of sparse matrix. It uses the sparse BLAS
to do matrix operations:
math.nist.gov/spblas/

SparseMatrix

---- WH
*/

#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

#include "Others/randgen.h"
#include <cstdarg>
#include <map>
#include <string>
#include "Others/SparseBLAS/blas_sparse.h"
#include "Others/def.h"
#include <sstream>

/*Define the namespace*/
namespace ROPTLIB{

    class SparseMatrix {
    public:
        
        /*construct real sparse matrix. ir and jc are arrays of indices of nonzero entries.
        vals is an array of nonzero values. nz is the number of nonzero entries. */
        SparseMatrix(integer inm, integer inn, integer *ir, integer *jc, realdp *vals, integer nz);
        
        /*construct complex sparse matrix. ir and jc are arrays of indices of nonzero entries.
        vals is an array of nonzero values. nz is the number of nonzero entries. */
        SparseMatrix(integer inm, integer inn, integer *ir, integer *jc, realdpcomplex *vals, integer nz);
        
        void Print(const char *name = "") const;
        
        inline blas_sparse_matrix GetSparseM(void) const { return SparseM; };
        
        inline bool Getiscomplex(void) const { return iscomplex; };
        
        inline integer Getrow(void) const { return m; };
        
        inline integer Getcol(void) const { return n; };
        
        ~SparseMatrix(void);
        
    protected:
        integer m; /*number of rows*/
        integer n; /*number of columns*/
        
        blas_sparse_matrix SparseM;
        bool iscomplex;
    };
}; /*end of ROPTLIB namespace*/

#endif

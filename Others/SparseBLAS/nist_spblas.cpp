/*
This is the C++ version of Sparse BLAS downloaded from
math.nist.gov/spblas/
I slightly modified it such that it is compatible with
ROPTLIB.


--Wen Huang
*/

/*
*
* Sparse BLAS (Basic Linear Algebra Subprograms) Library
*
* A C++ implementation of the routines specified by the ANSI C 
* interface specification of the Sparse BLAS in the BLAS Technical 
* Forum Standard[1].   For details, see [2].
*
* Mathematical and Computational Sciences Division
* National Institute of Technology,
* Gaithersburg, MD USA
*
*
* [1] BLAS Technical Forum: www.netlib.org/blas/blast-forum/
* [2] I. S. Duff, M. A. Heroux, R. Pozo, "An Overview of the Sparse Basic
*     Linear Algebra Subprograms: The new standard of the BLAS Techincal
*     Forum,"  Vol. 28, No. 2, pp. 239-267,ACM Transactions on Mathematical 
*     Software (TOMS), 2002.
*
*
* DISCLAIMER:
*
* This software was developed at the National Institute of Standards and
* Technology (NIST) by employees of the Federal Government in the course
* of their official duties. Pursuant to title 17 Section 105 of the
* United States Code, this software is not subject to copyright protection
* and is in the public domain. NIST assumes no responsibility whatsoever for
* its use by other parties, and makes no guarantees, expressed or implied,
* about its quality, reliability, or any other characteristic.
*
*
*/

#include "Others/SparseBLAS/nist_spblas.h"


namespace NIST_SPBLAS{

    /*
      finds an empty slot in global sparse marix table, and fills it with
      given entry.  Returns -1 if no spot found, or table is corrupt.
    */
    int Table_insert(Sp_mat* S)
    {
      if (Table_active_matrices <= Table.size())
      {
        Table.push_back(S);
        Table_active_matrices++;
        return static_cast<int> (Table.size() - 1);
      }
      else
      {
        /* there is an available slot; find it. */
        for (unsigned int i=0; i<Table.size(); i++)
        {
          if (Table[i] == NULL)
          {
            Table[i] = S;
            Table_active_matrices++;
            return i;
          }
        }
      }

      return -1;
    }

    /*
      removes an exisiting sparse matrix from global table.  Returns 0, if
      successfull, 1 otherwise.
    */
    int Table_remove(unsigned int i)
    {
      if (i < Table.size() && Table[i] != NULL)
      {
        Table[i] = NULL;
        Table_active_matrices--;
        return 0;
      }
      else
        return -1;
    }

    Sp_mat *GetSpMFromTable(unsigned int i)
    {
        if (i < Table.size() && Table[i] != nullptr)
        {
          return Table[i];
        }
        else
          return nullptr;
    };

    void Sp_mat::print() const
    {
      

      std::cout << "State : " <<
        (is_void() ? "void" :
         is_new()  ? "new" :
         is_open() ? "open" :
         is_valid() ? "valid" : "unknown") << "\n";

      std::cout << "M = " <<  num_rows() <<  "  N = " << num_cols() <<
            "  nz = " << num_nonzeros() << "\n";

    #define yesno(exp) ( (exp) ? "yes" : "no" )

      std::cout << "real: "     << yesno(is_real()) << "\n";
      std::cout << "complex: "  << yesno(is_complex()) << "\n";
      std::cout << "double "    << yesno(is_double_precision()) << "\n";
      std::cout << "single "    << yesno(is_single_precision()) << "\n";

      std::cout << "upper_triangular: " << yesno(is_upper_triangular()) << "\n";
      std::cout << "lower_triangular: " << yesno(is_lower_triangular()) << "\n";

      std::cout << "regular:    " << yesno(is_opt_regular()) << "\n";
      std::cout << "irregular:  " << yesno(is_opt_irregular()) << "\n";
      std::cout << "block:      " << yesno(is_opt_block()) << "\n";
      std::cout << "unassembled:" << yesno(is_opt_unassembled()) << "\n";
          
    #undef yesno
    }

    void table_print()
    {
      std::cout << "Table has " << Table.size() << " element(s). \n";
      for (unsigned int i=0; i< Table.size(); i++)
      {
        if (Table[i] != 0)
        {
          std::cout << "***** Table[" << i << "]: \n";
          Table[i]->print();
          std::cout << "\n\n";
        }
      }
    }

    void print(int A)
    {
      std::cout << "\n";
      Table[A]->print();
      std::cout << "\n";
    }
};

using namespace NIST_SPBLAS;

/* Level 1 */

/* these macros are useful for creating some consistency between the
   various precisions and floating point types.
*/

typedef float    FLOAT;
typedef double   DOUBLE;
typedef std::complex<float> COMPLEX_FLOAT;
typedef std::complex<double> COMPLEX_DOUBLE;



typedef float          SPBLAS_FLOAT_IN;
typedef double         SPBLAS_DOUBLE_IN;
typedef const void *   SPBLAS_COMPLEX_FLOAT_IN;
typedef const void *   SPBLAS_COMPLEX_DOUBLE_IN;

typedef float *  SPBLAS_FLOAT_OUT;
typedef double * SPBLAS_DOUBLE_OUT;
typedef void *   SPBLAS_COMPLEX_FLOAT_OUT;
typedef void *   SPBLAS_COMPLEX_DOUBLE_OUT;



typedef float *  SPBLAS_FLOAT_IN_OUT;
typedef double * SPBLAS_DOUBLE_IN_OUT;
typedef void *   SPBLAS_COMPLEX_FLOAT_IN_OUT;
typedef void *   SPBLAS_COMPLEX_DOUBLE_IN_OUT;

typedef const float *  SPBLAS_VECTOR_FLOAT_IN;
typedef const double * SPBLAS_VECTOR_DOUBLE_IN;
typedef const void *   SPBLAS_VECTOR_COMPLEX_FLOAT_IN;
typedef const void *   SPBLAS_VECTOR_COMPLEX_DOUBLE_IN;

typedef float *  SPBLAS_VECTOR_FLOAT_OUT;
typedef double * SPBLAS_VECTOR_DOUBLE_OUT;
typedef void *   SPBLAS_VECTOR_COMPLEX_FLOAT_OUT;
typedef void *   SPBLAS_VECTOR_COMPLEX_DOUBLE_OUT;

typedef float *  SPBLAS_VECTOR_FLOAT_IN_OUT;
typedef double * SPBLAS_VECTOR_DOUBLE_IN_OUT;
typedef void *   SPBLAS_VECTOR_COMPLEX_FLOAT_IN_OUT;
typedef void *   SPBLAS_VECTOR_COMPLEX_DOUBLE_IN_OUT;



#define SPBLAS_TO_FLOAT_IN(x)   x
#define SPBLAS_TO_DOUBLE_IN(x)  x
#define SPBLAS_TO_COMPLEX_FLOAT_IN(x) \
        (* reinterpret_cast<const std::complex<float> *>(x))
#define SPBLAS_TO_COMPLEX_DOUBLE_IN(x)  \
        (* reinterpret_cast<const std::complex<double> *>(x))

#define SPBLAS_TO_FLOAT_OUT(x)  x
#define SPBLAS_TO_DOUBLE_OUT(x) x
#define SPBLAS_TO_COMPLEX_FLOAT_OUT(x)  reinterpret_cast<std::complex<float> *>(x)
#define SPBLAS_TO_COMPLEX_DOUBLE_OUT(x) reinterpret_cast<std::complex<double> *>(x)

#define SPBLAS_TO_FLOAT_IN_OUT(x)   x
#define SPBLAS_TO_DOUBLE_IN_OUT(x)  x
#define SPBLAS_TO_COMPLEX_FLOAT_IN_OUT(x)  reinterpret_cast<std::complex<float> *>(x)
#define SPBLAS_TO_COMPLEX_DOUBLE_IN_OUT(x) reinterpret_cast<std::complex<double>*>(x)


#define SPBLAS_TO_VECTOR_DOUBLE_IN(x)   x
#define SPBLAS_TO_VECTOR_FLOAT_IN(x)  x
#define SPBLAS_TO_VECTOR_COMPLEX_FLOAT_IN(x) \
                          reinterpret_cast<const std::complex<float>*>(x)
#define SPBLAS_TO_VECTOR_COMPLEX_DOUBLE_IN(x) \
                          reinterpret_cast<const std::complex<double>*>(x)

#define SPBLAS_TO_VECTOR_DOUBLE_OUT(x)  x
#define SPBLAS_TO_VECTOR_FLOAT_OUT(x)   x
#define SPBLAS_TO_VECTOR_COMPLEX_FLOAT_OUT(x) \
                          reinterpret_cast<std::complex<float>*>(x)
#define SPBLAS_TO_VECTOR_COMPLEX_DOUBLE_OUT(x) \
                          reinterpret_cast<std::complex<double>*>(x)

#define SPBLAS_TO_VECTOR_DOUBLE_IN_OUT(x)   x
#define SPBLAS_TO_VECTOR_FLOAT_IN_OUT(x)  x
#define SPBLAS_TO_VECTOR_COMPLEX_FLOAT_IN_OUT(x) \
                          reinterpret_cast<std::complex<float>*>(x)
#define SPBLAS_TO_VECTOR_COMPLEX_DOUBLE_IN_OUT(x) \
                          reinterpret_cast<std::complex<double>*>(x)






#define BLAS_FLOAT_NAME(routine_name) BLAS_s##routine_name
#define BLAS_DOUBLE_NAME(routine_name) BLAS_d##routine_name
#define BLAS_COMPLEX_FLOAT_NAME(routine_name) BLAS_c##routine_name
#define BLAS_COMPLEX_DOUBLE_NAME(routine_name) BLAS_z##routine_name


#define TSp_MAT_SET_FLOAT(A) {A->set_single_precision(); A->set_real();}
#define TSp_MAT_SET_DOUBLE(A) {A->set_double_precision(); A->set_real();}
#define TSp_MAT_SET_COMPLEX_FLOAT(A) {A->set_single_precision(); A->set_complex();}
#define TSp_MAT_SET_COMPLEX_DOUBLE(A) {A->set_double_precision(); A->set_complex();}


        /*------------------------------------*/
        /* Non-precision Sparse BLAS routines */
        /*------------------------------------*/

/* -------- */
/*  USSP()  */
/* -------- */
int BLAS_ussp(blas_sparse_matrix A, int pname)
{
    Sp_mat *S = Table[A];


                                /* Note: these are returns, in the case */
                                /* statement, so "break" is not needed.  */
    switch (pname)
    {
        case (blas_zero_base) : S->set_zero_base(); break;
        case (blas_one_base)  : S->set_one_base(); break;

        case (blas_unit_diag) : S->set_unit_diag(); break;
        case (blas_complex)   : S->set_complex(); break;
        case (blas_real)      : S->set_real(); break;

        case (blas_single_precision) : S->set_single_precision(); break;
        case (blas_double_precision) : S->set_double_precision(); break;


        case (blas_lower_triangular) : S->set_lower_triangular(); break;
        case (blas_upper_triangular) : S->set_upper_triangular(); break;

        case (blas_lower_symmetric) : S->set_lower_symmetric(); break;
        case (blas_upper_symmetric) : S->set_upper_symmetric(); break;

        case (blas_lower_hermitian) : S->set_lower_hermitian(); break;
        case (blas_upper_hermitian) : S->set_upper_hermitian(); break;


                                        /* optimizations not used */
        case (blas_regular )        :
        case (blas_irregular)       :
        case (blas_block)           :
        case (blas_unassembled)     :    return 0;

        default:                return -1;  /* invalid property */
    }

    return 0;

}

/* -------- */
/*  USGP()  */
/* -------- */

int BLAS_usgp(blas_sparse_matrix A, int pname)
{
  Sp_mat *S = Table[A];

    switch (pname)
    {
        case (blas_num_rows)          : return S->num_rows();
        case (blas_num_cols)          : return S->num_cols();
        case (blas_num_nonzeros)      : return S->num_nonzeros();

        case (blas_complex)           : return S->is_complex();
        case (blas_real)              : return S->is_real();
        case (blas_single_precision)  : return S->is_single_precision();
        case (blas_double_precision)  : return S->is_double_precision();


        case (blas_lower_triangular) : return S->is_lower_triangular(); break;
        case (blas_upper_triangular) : return S->is_upper_triangular(); break;

        case (blas_general)           : return S->is_general();
        case (blas_symmetric)         : return S->is_symmetric();
        case (blas_hermitian)         : return S->is_hermitian();

        case (blas_zero_base) : return S->is_zero_base();
        case (blas_one_base) : return  S->is_one_base();


        case (blas_rowmajor)          : return S->is_rowmajor();
        case (blas_colmajor)          : return S->is_colmajor();
        case (blas_new_handle)        : return S->is_new();
        case (blas_valid_handle)      : return S->is_valid();
        case (blas_open_handle)       : return S->is_open();
        case (blas_invalid_handle)    : return S->is_void();

        case (blas_regular)           : return S->is_opt_regular();
        case (blas_irregular)         : return S->is_opt_irregular();
        case (blas_block)             : return S->is_opt_block();
        case (blas_unassembled)       : return S->is_opt_unassembled();

        default:                return -1;  /* invalid property */
    }
}

/* -------- */
/*  USDS()  */
/* -------- */
int BLAS_usds(int A)
{
  Sp_mat *S  = Table[A];
  S->destroy();
  delete S; /* Added to avoid memory leaks --- by WH */
  Table_remove(A);

  return 0;
}




/* --------------------------- */
/*  Level 1 generic routines   */
/* --------------------------- */

template <class T>
void BLAS_xusdot( enum blas_conj_type conj_flag, int nz,
    const T *x,  const int *index,  const T *y, int incy,
    T *r, enum blas_base_type index_base)
{

  T t(0);

  if (index_base == blas_one_base)
    y -= incy;

  if (conj_flag == blas_no_conj)
  {
    for (int i=0; i<nz; i++)
      t += x[i] * y[index[i]*incy];
  }
  else
    for (int i=0; i<nz; i++)
      t += conj(x[i]) * y[index[i]*incy];


  *r = t;
}


template <class T>
void BLAS_xusaxpy(int nz, T alpha, const T *x,
    const int *index, T *y, int incy,
    enum blas_base_type index_base)
{

  if (index_base == blas_one_base)
    y -= incy;

  for (int i=0; i<nz; i++)
  {
/*     y[index[i]*incy] +=  (alpha * x[i]); */

  }
}

template <class T>
void BLAS_xusga( int nz, const T *y, int incy, T *x, const int *indx,
              enum blas_base_type index_base )
{
  if (index_base == blas_one_base)
    y -= incy;

  for (int i=0; i<nz; i++)
    x[i] = y[indx[i]*incy];

}

template <class T>
void BLAS_xusgz( int nz, T *y, int incy, T *x, const int *indx,
              enum blas_base_type index_base )
{
  if (index_base == blas_one_base)
    y -= incy;

  for (int i=0; i<nz; i++)
  {
    x[i] = y[indx[i]*incy];
    y[indx[i]*incy] = (T) 0.0;
  }

}


template <class T>
void BLAS_xussc(int nz, const T *x, T *y, int incy, const int *index,
    enum blas_base_type index_base)
{
  if (index_base == blas_one_base)
    y -= incy;

  for (int i=0; i<nz; i++)
    y[index[i]*incy] = x[i];

}

/* --------------------------- */
/* Level 2&3 generic precision */
/* --------------------------- */


template <class T>
int BLAS_xuscr_insert_entry(blas_sparse_matrix A,  const T& val, int i, int j)
{
  return ((TSp_mat<T> *)Table[A])->insert_entry(val, i, j);

}

template <class T>
int BLAS_xuscr_insert_entries(blas_sparse_matrix A, int nz, const T* Val,
      const int* I, const int *J)
{
  return ((TSp_mat<T>*) Table[A])->insert_entries(nz, Val, I, J);
}


template <class T>
int BLAS_xuscr_insert_col(blas_sparse_matrix A, int j, int nz, const T* Val,
      const int* indx)
{
  return ((TSp_mat<T>*) Table[A])->insert_col(j, nz, Val, indx);
}

template <class T>
int BLAS_xuscr_insert_row(blas_sparse_matrix A, int i, int nz, const T* Val,
      const int* indx)
{
  return ((TSp_mat<T>*) Table[A])->insert_row(i, nz, Val, indx);
}

template <class T>
int BLAS_xuscr_insert_clique(blas_sparse_matrix A, int k, int l, const T* Val,
      const int row_stride, const int col_stride, const int *indx,
      const int *jndx)
{
  return ((TSp_mat<T>*) Table[A])->insert_clique(k, l, Val, row_stride,
            col_stride, indx, jndx);
}

template <class T>
int BLAS_xuscr_insert_block(blas_sparse_matrix A, const T* Val,
      const int row_stride, const int col_stride, int bi, int bj )
{
  return ((TSp_mat<T>*) Table[A])->insert_block(Val,
        row_stride, col_stride, bi, bj);
}

inline int BLAS_xuscr_end(blas_sparse_matrix A)
{
  return (Table[A])->end_construction();
}


template <class T>
int BLAS_xusmv(enum blas_trans_type transa, const T& alpha,
    blas_sparse_matrix A, const T *x, int incx, T *y, int incy )
{
  TSp_mat<T> *M = (TSp_mat<T> *) Table[A];

  ASSERT_RETURN(M->is_valid(), 1);

  return M->usmv(transa, alpha, x, incx, y, incy);
}


template <class T>
int BLAS_xusmm(enum blas_order_type ordera, enum blas_trans_type transa,
    int nrhs, const T& alpha, blas_sparse_matrix A,
    const T *B, int ldB, T* C, int ldC)
{
  TSp_mat<T> *M = (TSp_mat<T> *) Table[A];

  ASSERT_RETURN(M->is_valid(), 1);

  return M->usmm(ordera, transa, nrhs, alpha, B, ldB, C, ldC);
}


template <class T>
int BLAS_xussv(enum blas_trans_type transa, const T& alpha,
    blas_sparse_matrix A, T *x, int incx)
{
  TSp_mat<T> *M =
        (TSp_mat<T> *) Table[A];

  ASSERT_RETURN(M->is_valid(), 1);

  return M->ussv(transa, alpha, x, incx);
}

template <class T>
int BLAS_xussm(enum blas_order_type orderA, enum blas_trans_type transa,
    int nrhs, const T& alpha, blas_sparse_matrix A, T *C, int ldC)
{
  TSp_mat<T> *M =
        (TSp_mat<T> *) Table[A];

  ASSERT_RETURN(M->is_valid(), 1);

  return M->ussm(orderA, transa, nrhs, alpha, C, ldC);
}

/*   --- end of generic rouintes ---- */


/*********/
/*  ----  double  Level 1 rouintes ----- */
/*********/

 void BLAS_DOUBLE_NAME(usdot)(
    enum blas_conj_type conj_flag,
    int  nz,
    SPBLAS_VECTOR_DOUBLE_IN x,
    const int *index,
    SPBLAS_VECTOR_DOUBLE_IN y,
    int incy,
    SPBLAS_DOUBLE_OUT r,
    enum blas_base_type index_base)
{
   BLAS_xusdot(conj_flag, nz,
          SPBLAS_TO_VECTOR_DOUBLE_IN( x ), index,
          SPBLAS_TO_VECTOR_DOUBLE_IN( y ), incy,
          SPBLAS_TO_DOUBLE_OUT( r ), index_base);
}

 void BLAS_DOUBLE_NAME(usaxpy)(
    int nz,
    SPBLAS_DOUBLE_IN alpha,
    SPBLAS_VECTOR_DOUBLE_IN x,
    const int *index,
    SPBLAS_VECTOR_DOUBLE_IN_OUT y,
    int incy,
    enum blas_base_type index_base)
{
  BLAS_xusaxpy(nz, SPBLAS_TO_DOUBLE_IN( alpha ),
      SPBLAS_TO_VECTOR_DOUBLE_IN( x ), index,
      SPBLAS_TO_VECTOR_DOUBLE_IN_OUT( y ),
      incy, index_base);
}

 void BLAS_DOUBLE_NAME(usga)(
    int nz,
    SPBLAS_VECTOR_DOUBLE_IN y,
    int incy,
    SPBLAS_VECTOR_DOUBLE_IN_OUT x,
    const int *indx,
    enum blas_base_type index_base )
{
  BLAS_xusga( nz, SPBLAS_TO_VECTOR_DOUBLE_IN( y ), incy,
        SPBLAS_TO_VECTOR_DOUBLE_IN_OUT( x ), indx, index_base);

}

 void BLAS_DOUBLE_NAME(usgz)(
    int nz,
    SPBLAS_VECTOR_DOUBLE_IN_OUT y,
    int incy,
    SPBLAS_VECTOR_DOUBLE_OUT x,
    const int *indx,
    enum blas_base_type index_base )
{
  BLAS_xusgz(nz, SPBLAS_TO_DOUBLE_IN_OUT(y ), incy, SPBLAS_TO_DOUBLE_OUT(  x ),
    indx, index_base);
}


 void BLAS_DOUBLE_NAME(ussc)(
    int nz,
    SPBLAS_VECTOR_DOUBLE_IN x,
    SPBLAS_VECTOR_DOUBLE_IN_OUT y,
    int incy,
    const int *index,
    enum blas_base_type index_base)
{
  BLAS_xussc(nz, SPBLAS_TO_VECTOR_DOUBLE_IN( x ),
  SPBLAS_TO_DOUBLE_IN_OUT( y ), incy, index,
  index_base);
}

/*  DOUBLE Level 2/3 creation routines */

int BLAS_DOUBLE_NAME(uscr_begin)(int M, int N)
{
  TSp_mat<DOUBLE> *A = new TSp_mat<DOUBLE>(M, N);
  TSp_MAT_SET_DOUBLE(A);
  return Table_insert(A);;
}

blas_sparse_matrix BLAS_DOUBLE_NAME(uscr_block_begin)(
    int Mb, int Nb, int k, int l )
{
  TSp_mat<DOUBLE> *A = new TSp_mat<DOUBLE>(Mb*k, Nb*l);

  TSp_MAT_SET_DOUBLE(A);
  A->set_const_block_parameters(Mb, Nb, k, l);

  return Table_insert(A);
}

blas_sparse_matrix BLAS_DOUBLE_NAME(uscr_variable_block_begin)(
    int Mb, int Nb, const int *k, const int *l )
{
  TSp_mat<DOUBLE> *A = new TSp_mat<DOUBLE>(
                std::accumulate(k, k+Mb, 0), std::accumulate(l, l+Nb, 0) );

  TSp_MAT_SET_DOUBLE(A);
  A->set_var_block_parameters(Mb, Nb, k, l);

  return Table_insert(A);

}

/*  DOUBLE Level 2/3 insertion routines */

int BLAS_DOUBLE_NAME(uscr_insert_entry)(
    blas_sparse_matrix A, SPBLAS_DOUBLE_IN val, int i, int j )
{
  return BLAS_xuscr_insert_entry(A, SPBLAS_TO_DOUBLE_IN( val ), i, j);
}

int BLAS_DOUBLE_NAME(uscr_insert_entries)(
    blas_sparse_matrix A, int nz,
    SPBLAS_VECTOR_DOUBLE_IN val,
    const int *indx, const int *jndx )
{
  return BLAS_xuscr_insert_entries(A, nz, SPBLAS_TO_VECTOR_DOUBLE_IN( val ), indx, jndx);
}

int BLAS_DOUBLE_NAME(uscr_insert_col)(
    blas_sparse_matrix A, int j, int nz,
    SPBLAS_VECTOR_DOUBLE_IN val, const int *indx )
{
  return BLAS_xuscr_insert_col(A, j, nz, SPBLAS_TO_VECTOR_DOUBLE_IN( val ), indx);
}

int BLAS_DOUBLE_NAME(uscr_insert_row)(
  blas_sparse_matrix A, int i, int nz,
  SPBLAS_VECTOR_DOUBLE_IN val, const int *indx );

int BLAS_DOUBLE_NAME(uscr_insert_clique)(
    blas_sparse_matrix A,
    const int k,
    const int l,
    SPBLAS_VECTOR_DOUBLE_IN val,
    const int row_stride,
    const int col_stride,
    const int *indx,
    const int *jndx );

int BLAS_DOUBLE_NAME(uscr_insert_block)(
    blas_sparse_matrix A,
    SPBLAS_VECTOR_DOUBLE_IN val,
    int row_stride,
    int col_stride,
    int i, int j )
{
  return BLAS_xuscr_insert_block(
        A, SPBLAS_TO_VECTOR_DOUBLE_IN( val ),
        row_stride, col_stride, i, j);
}

int BLAS_DOUBLE_NAME(uscr_end)(blas_sparse_matrix A)
{
  return BLAS_xuscr_end(A);
}

/*  DOUBLE Level 2/3 computational routines */

 int BLAS_DOUBLE_NAME(usmv)(enum
    blas_trans_type transa,
    SPBLAS_DOUBLE_IN alpha,
    blas_sparse_matrix A,
    SPBLAS_VECTOR_DOUBLE_IN x,
    int incx,
    SPBLAS_VECTOR_DOUBLE_IN_OUT y,
    int incy )
{
  return BLAS_xusmv(
      transa,   SPBLAS_TO_DOUBLE_IN( alpha ), A,
      SPBLAS_TO_VECTOR_DOUBLE_IN( x ), incx,
      SPBLAS_TO_VECTOR_DOUBLE_IN_OUT( y ), incy);
}

int BLAS_DOUBLE_NAME(usmm)(
    enum blas_order_type order,
    enum blas_trans_type transa,
    int nrhs,
    SPBLAS_DOUBLE_IN alpha,
    blas_sparse_matrix A,
    SPBLAS_VECTOR_DOUBLE_IN b,
    int ldb,
    SPBLAS_VECTOR_DOUBLE_IN_OUT c,
    int ldc )
{
  return BLAS_xusmm(
      order, transa, nrhs,
      SPBLAS_TO_DOUBLE_IN( alpha), A,
      SPBLAS_TO_VECTOR_DOUBLE_IN(b), ldb,
      SPBLAS_TO_VECTOR_DOUBLE_IN_OUT( c ), ldc);
}

int BLAS_DOUBLE_NAME(ussv)(
    enum blas_trans_type transa,
    SPBLAS_DOUBLE_IN alpha,
    blas_sparse_matrix A,
    SPBLAS_VECTOR_DOUBLE_IN_OUT x,
    int incx )
{
  return BLAS_xussv( transa,
        SPBLAS_TO_DOUBLE_IN( alpha ), A,
        SPBLAS_TO_VECTOR_DOUBLE_IN_OUT( x ),
        incx);
}


int BLAS_DOUBLE_NAME(ussm)(
    enum blas_order_type order,
    enum blas_trans_type transt,
    int nrhs,
    SPBLAS_DOUBLE_IN alpha,
    blas_sparse_matrix A,
    SPBLAS_VECTOR_DOUBLE_IN_OUT b,
    int ldb )
{
  return BLAS_xussm(order, transt, nrhs,
      SPBLAS_TO_DOUBLE_IN( alpha ), A,
      SPBLAS_TO_VECTOR_DOUBLE_IN_OUT( b ), ldb);
}


/*  ----   end of DOUBLE routines -------  */




 void BLAS_COMPLEX_DOUBLE_NAME(usdot)(
    enum blas_conj_type conj_flag,
    int  nz,
    SPBLAS_VECTOR_COMPLEX_DOUBLE_IN x,
    const int *index,
    SPBLAS_VECTOR_COMPLEX_DOUBLE_IN y,
    int incy,
    SPBLAS_COMPLEX_DOUBLE_OUT r,
    enum blas_base_type index_base)
{
   BLAS_xusdot(conj_flag, nz,
          SPBLAS_TO_VECTOR_COMPLEX_DOUBLE_IN( x ), index,
          SPBLAS_TO_VECTOR_COMPLEX_DOUBLE_IN( y ), incy,
          SPBLAS_TO_COMPLEX_DOUBLE_OUT( r ), index_base);
}

 void BLAS_COMPLEX_DOUBLE_NAME(usaxpy)(
    int nz,
    SPBLAS_COMPLEX_DOUBLE_IN alpha,
    SPBLAS_VECTOR_COMPLEX_DOUBLE_IN x,
    const int *index,
    SPBLAS_VECTOR_COMPLEX_DOUBLE_IN_OUT y,
    int incy,
    enum blas_base_type index_base)
{
  BLAS_xusaxpy(nz, SPBLAS_TO_COMPLEX_DOUBLE_IN( alpha ),
      SPBLAS_TO_VECTOR_COMPLEX_DOUBLE_IN( x ), index,
      SPBLAS_TO_VECTOR_COMPLEX_DOUBLE_IN_OUT( y ),
      incy, index_base);
}

 void BLAS_COMPLEX_DOUBLE_NAME(usga)(
    int nz,
    SPBLAS_VECTOR_COMPLEX_DOUBLE_IN y,
    int incy,
    SPBLAS_VECTOR_COMPLEX_DOUBLE_IN_OUT x,
    const int *indx,
    enum blas_base_type index_base )
{
  BLAS_xusga( nz, SPBLAS_TO_VECTOR_COMPLEX_DOUBLE_IN( y ), incy,
        SPBLAS_TO_VECTOR_COMPLEX_DOUBLE_IN_OUT( x ), indx, index_base);

}

 void BLAS_COMPLEX_DOUBLE_NAME(usgz)(
    int nz,
    SPBLAS_VECTOR_COMPLEX_DOUBLE_IN_OUT y,
    int incy,
    SPBLAS_VECTOR_COMPLEX_DOUBLE_OUT x,
    const int *indx,
    enum blas_base_type index_base )
{
  BLAS_xusgz(nz, SPBLAS_TO_COMPLEX_DOUBLE_IN_OUT(y ), incy, SPBLAS_TO_COMPLEX_DOUBLE_OUT(  x ),
    indx, index_base);
}


 void BLAS_COMPLEX_DOUBLE_NAME(ussc)(
    int nz,
    SPBLAS_VECTOR_COMPLEX_DOUBLE_IN x,
    SPBLAS_VECTOR_COMPLEX_DOUBLE_IN_OUT y,
    int incy,
    const int *index,
    enum blas_base_type index_base)
{
  BLAS_xussc(nz, SPBLAS_TO_VECTOR_COMPLEX_DOUBLE_IN( x ),
  SPBLAS_TO_COMPLEX_DOUBLE_IN_OUT( y ), incy, index,
  index_base);
}

/*  COMPLEX_DOUBLE Level 2/3 creation routines */

int BLAS_COMPLEX_DOUBLE_NAME(uscr_begin)(int M, int N)
{
  TSp_mat<COMPLEX_DOUBLE> *A = new TSp_mat<COMPLEX_DOUBLE>(M, N);
  TSp_MAT_SET_COMPLEX_DOUBLE(A);

  return Table_insert(A);
}

blas_sparse_matrix BLAS_COMPLEX_DOUBLE_NAME(uscr_block_begin)(
    int Mb, int Nb, int k, int l )
{
  TSp_mat<COMPLEX_DOUBLE> *A = new TSp_mat<COMPLEX_DOUBLE>(Mb*k, Nb*l);

  TSp_MAT_SET_COMPLEX_DOUBLE(A);
  A->set_const_block_parameters(Mb, Nb, k, l);

  return Table_insert(A);
}

blas_sparse_matrix BLAS_COMPLEX_DOUBLE_NAME(uscr_variable_block_begin)(
    int Mb, int Nb, const int *k, const int *l )
{
  TSp_mat<COMPLEX_DOUBLE> *A = new TSp_mat<COMPLEX_DOUBLE>(
                std::accumulate(k, k+Mb, 0), std::accumulate(l, l+Nb, 0) );

  TSp_MAT_SET_COMPLEX_DOUBLE(A);
  A->set_var_block_parameters(Mb, Nb, k, l);

  return Table_insert(A);

}

/*  COMPLEX_DOUBLE Level 2/3 insertion routines */

int BLAS_COMPLEX_DOUBLE_NAME(uscr_insert_entry)(
    blas_sparse_matrix A, SPBLAS_COMPLEX_DOUBLE_IN val, int i, int j )
{
  return BLAS_xuscr_insert_entry(A, SPBLAS_TO_COMPLEX_DOUBLE_IN( val ), i, j);
}

int BLAS_COMPLEX_DOUBLE_NAME(uscr_insert_entries)(
    blas_sparse_matrix A, int nz,
    SPBLAS_VECTOR_COMPLEX_DOUBLE_IN val,
    const int *indx, const int *jndx )
{
  return BLAS_xuscr_insert_entries(A, nz, SPBLAS_TO_VECTOR_COMPLEX_DOUBLE_IN( val ), indx, jndx);
}

int BLAS_COMPLEX_DOUBLE_NAME(uscr_insert_col)(
    blas_sparse_matrix A, int j, int nz,
    SPBLAS_VECTOR_COMPLEX_DOUBLE_IN val, const int *indx )
{
  return BLAS_xuscr_insert_col(A, j, nz, SPBLAS_TO_VECTOR_COMPLEX_DOUBLE_IN( val ), indx);
}

int BLAS_COMPLEX_DOUBLE_NAME(uscr_insert_row)(
  blas_sparse_matrix A, int i, int nz,
  SPBLAS_VECTOR_COMPLEX_DOUBLE_IN val, const int *indx );

int BLAS_COMPLEX_DOUBLE_NAME(uscr_insert_clique)(
    blas_sparse_matrix A,
    const int k,
    const int l,
    SPBLAS_VECTOR_COMPLEX_DOUBLE_IN val,
    const int row_stride,
    const int col_stride,
    const int *indx,
    const int *jndx );

int BLAS_COMPLEX_DOUBLE_NAME(uscr_insert_block)(
    blas_sparse_matrix A,
    SPBLAS_VECTOR_COMPLEX_DOUBLE_IN val,
    int row_stride,
    int col_stride,
    int i, int j )
{
  return BLAS_xuscr_insert_block(
        A, SPBLAS_TO_VECTOR_COMPLEX_DOUBLE_IN( val ),
        row_stride, col_stride, i, j);
}

int BLAS_COMPLEX_DOUBLE_NAME(uscr_end)(blas_sparse_matrix A)
{
  return BLAS_xuscr_end(A);
}

/*  COMPLEX_DOUBLE Level 2/3 computational routines */

 int BLAS_COMPLEX_DOUBLE_NAME(usmv)(enum
    blas_trans_type transa,
    SPBLAS_COMPLEX_DOUBLE_IN alpha,
    blas_sparse_matrix A,
    SPBLAS_VECTOR_COMPLEX_DOUBLE_IN x,
    int incx,
    SPBLAS_VECTOR_COMPLEX_DOUBLE_IN_OUT y,
    int incy )
{
  return BLAS_xusmv(
      transa,   SPBLAS_TO_COMPLEX_DOUBLE_IN( alpha ), A,
      SPBLAS_TO_VECTOR_COMPLEX_DOUBLE_IN( x ), incx,
      SPBLAS_TO_VECTOR_COMPLEX_DOUBLE_IN_OUT( y ), incy);
}

int BLAS_COMPLEX_DOUBLE_NAME(usmm)(
    enum blas_order_type order,
    enum blas_trans_type transa,
    int nrhs,
    SPBLAS_COMPLEX_DOUBLE_IN alpha,
    blas_sparse_matrix A,
    SPBLAS_VECTOR_COMPLEX_DOUBLE_IN b,
    int ldb,
    SPBLAS_VECTOR_COMPLEX_DOUBLE_IN_OUT c,
    int ldc )
{
  return BLAS_xusmm(
      order, transa, nrhs,
      SPBLAS_TO_COMPLEX_DOUBLE_IN( alpha), A,
      SPBLAS_TO_VECTOR_COMPLEX_DOUBLE_IN(b), ldb,
      SPBLAS_TO_VECTOR_COMPLEX_DOUBLE_IN_OUT( c ), ldc);
}

int BLAS_COMPLEX_DOUBLE_NAME(ussv)(
    enum blas_trans_type transa,
    SPBLAS_COMPLEX_DOUBLE_IN alpha,
    blas_sparse_matrix A,
    SPBLAS_VECTOR_COMPLEX_DOUBLE_IN_OUT x,
    int incx )
{
  return BLAS_xussv( transa,
        SPBLAS_TO_COMPLEX_DOUBLE_IN( alpha ), A,
        SPBLAS_TO_VECTOR_COMPLEX_DOUBLE_IN_OUT( x ),
        incx);
}


int BLAS_COMPLEX_DOUBLE_NAME(ussm)(
    enum blas_order_type order,
    enum blas_trans_type transt,
    int nrhs,
    SPBLAS_COMPLEX_DOUBLE_IN alpha,
    blas_sparse_matrix A,
    SPBLAS_VECTOR_COMPLEX_DOUBLE_IN_OUT b,
    int ldb )
{
  return BLAS_xussm(order, transt, nrhs,
      SPBLAS_TO_COMPLEX_DOUBLE_IN( alpha ), A,
      SPBLAS_TO_VECTOR_COMPLEX_DOUBLE_IN_OUT( b ), ldb);
}




/*  ----   end of COMPLEX_COMPLEX_COMPLEX_DOUBLE routines -------  */


/*********/
/*  ----  double  Level 1 rouintes ----- */
/*********/

 void BLAS_FLOAT_NAME(usdot)(
    enum blas_conj_type conj_flag,
    int  nz,
    SPBLAS_VECTOR_FLOAT_IN x,
    const int *index,
    SPBLAS_VECTOR_FLOAT_IN y,
    int incy,
    SPBLAS_FLOAT_OUT r,
    enum blas_base_type index_base)
{
   BLAS_xusdot(conj_flag, nz,
          SPBLAS_TO_VECTOR_FLOAT_IN( x ), index,
          SPBLAS_TO_VECTOR_FLOAT_IN( y ), incy,
          SPBLAS_TO_FLOAT_OUT( r ), index_base);
}

 void BLAS_FLOAT_NAME(usaxpy)(
    int nz,
    SPBLAS_FLOAT_IN alpha,
    SPBLAS_VECTOR_FLOAT_IN x,
    const int *index,
    SPBLAS_VECTOR_FLOAT_IN_OUT y,
    int incy,
    enum blas_base_type index_base)
{
  BLAS_xusaxpy(nz, SPBLAS_TO_FLOAT_IN( alpha ),
      SPBLAS_TO_VECTOR_FLOAT_IN( x ), index,
      SPBLAS_TO_VECTOR_FLOAT_IN_OUT( y ),
      incy, index_base);
}

 void BLAS_FLOAT_NAME(usga)(
    int nz,
    SPBLAS_VECTOR_FLOAT_IN y,
    int incy,
    SPBLAS_VECTOR_FLOAT_IN_OUT x,
    const int *indx,
    enum blas_base_type index_base )
{
  BLAS_xusga( nz, SPBLAS_TO_VECTOR_FLOAT_IN( y ), incy,
        SPBLAS_TO_VECTOR_FLOAT_IN_OUT( x ), indx, index_base);

}

 void BLAS_FLOAT_NAME(usgz)(
    int nz,
    SPBLAS_VECTOR_FLOAT_IN_OUT y,
    int incy,
    SPBLAS_VECTOR_FLOAT_OUT x,
    const int *indx,
    enum blas_base_type index_base )
{
  BLAS_xusgz(nz, SPBLAS_TO_FLOAT_IN_OUT(y ), incy, SPBLAS_TO_FLOAT_OUT(  x ),
    indx, index_base);
}


 void BLAS_FLOAT_NAME(ussc)(
    int nz,
    SPBLAS_VECTOR_FLOAT_IN x,
    SPBLAS_VECTOR_FLOAT_IN_OUT y,
    int incy,
    const int *index,
    enum blas_base_type index_base)
{
  BLAS_xussc(nz, SPBLAS_TO_VECTOR_FLOAT_IN( x ),
  SPBLAS_TO_FLOAT_IN_OUT( y ), incy, index,
  index_base);
}

/*  FLOAT Level 2/3 creation routines */

int BLAS_FLOAT_NAME(uscr_begin)(int M, int N)
{
  TSp_mat<FLOAT> *A = new TSp_mat<FLOAT>(M, N);
  TSp_MAT_SET_FLOAT(A);

  return Table_insert(A);
}

blas_sparse_matrix BLAS_FLOAT_NAME(uscr_block_begin)(
    int Mb, int Nb, int k, int l )
{
  TSp_mat<FLOAT> *A = new TSp_mat<FLOAT>(Mb*k, Nb*l);

  TSp_MAT_SET_FLOAT(A);
  A->set_const_block_parameters(Mb, Nb, k, l);

  return Table_insert(A);
}

blas_sparse_matrix BLAS_FLOAT_NAME(uscr_variable_block_begin)(
    int Mb, int Nb, const int *k, const int *l )
{
  TSp_mat<FLOAT> *A = new TSp_mat<FLOAT>(
                std::accumulate(k, k+Mb, 0), std::accumulate(l, l+Nb, 0) );

  TSp_MAT_SET_FLOAT(A);
  A->set_var_block_parameters(Mb, Nb, k, l);

  return Table_insert(A);

}

/*  FLOAT Level 2/3 insertion routines */

int BLAS_FLOAT_NAME(uscr_insert_entry)(
    blas_sparse_matrix A, SPBLAS_FLOAT_IN val, int i, int j )
{
  return BLAS_xuscr_insert_entry(A, SPBLAS_TO_FLOAT_IN( val ), i, j);
}

int BLAS_FLOAT_NAME(uscr_insert_entries)(
    blas_sparse_matrix A, int nz,
    SPBLAS_VECTOR_FLOAT_IN val,
    const int *indx, const int *jndx )
{
  return BLAS_xuscr_insert_entries(A, nz, SPBLAS_TO_VECTOR_FLOAT_IN( val ), indx, jndx);
}

int BLAS_FLOAT_NAME(uscr_insert_col)(
    blas_sparse_matrix A, int j, int nz,
    SPBLAS_VECTOR_FLOAT_IN val, const int *indx )
{
  return BLAS_xuscr_insert_col(A, j, nz, SPBLAS_TO_VECTOR_FLOAT_IN( val ), indx);
}

int BLAS_FLOAT_NAME(uscr_insert_row)(
  blas_sparse_matrix A, int i, int nz,
  SPBLAS_VECTOR_FLOAT_IN val, const int *indx );

int BLAS_FLOAT_NAME(uscr_insert_clique)(
    blas_sparse_matrix A,
    const int k,
    const int l,
    SPBLAS_VECTOR_FLOAT_IN val,
    const int row_stride,
    const int col_stride,
    const int *indx,
    const int *jndx );

int BLAS_FLOAT_NAME(uscr_insert_block)(
    blas_sparse_matrix A,
    SPBLAS_VECTOR_FLOAT_IN val,
    int row_stride,
    int col_stride,
    int i, int j )
{
  return BLAS_xuscr_insert_block(
        A, SPBLAS_TO_VECTOR_FLOAT_IN( val ),
        row_stride, col_stride, i, j);
}

int BLAS_FLOAT_NAME(uscr_end)(blas_sparse_matrix A)
{
  return BLAS_xuscr_end(A);
}

/*  FLOAT Level 2/3 computational routines */

 int BLAS_FLOAT_NAME(usmv)(enum
    blas_trans_type transa,
    SPBLAS_FLOAT_IN alpha,
    blas_sparse_matrix A,
    SPBLAS_VECTOR_FLOAT_IN x,
    int incx,
    SPBLAS_VECTOR_FLOAT_IN_OUT y,
    int incy )
{
  return BLAS_xusmv(
      transa,   SPBLAS_TO_FLOAT_IN( alpha ), A,
      SPBLAS_TO_VECTOR_FLOAT_IN( x ), incx,
      SPBLAS_TO_VECTOR_FLOAT_IN_OUT( y ), incy);
}

int BLAS_FLOAT_NAME(usmm)(
    enum blas_order_type order,
    enum blas_trans_type transa,
    int nrhs,
    SPBLAS_FLOAT_IN alpha,
    blas_sparse_matrix A,
    SPBLAS_VECTOR_FLOAT_IN b,
    int ldb,
    SPBLAS_VECTOR_FLOAT_IN_OUT c,
    int ldc )
{
  return BLAS_xusmm(
      order, transa, nrhs,
      SPBLAS_TO_FLOAT_IN( alpha), A,
      SPBLAS_TO_VECTOR_FLOAT_IN(b), ldb,
      SPBLAS_TO_VECTOR_FLOAT_IN_OUT( c ), ldc);
}

int BLAS_FLOAT_NAME(ussv)(
    enum blas_trans_type transa,
    SPBLAS_FLOAT_IN alpha,
    blas_sparse_matrix A,
    SPBLAS_VECTOR_FLOAT_IN_OUT x,
    int incx )
{
  return BLAS_xussv( transa,
        SPBLAS_TO_FLOAT_IN( alpha ), A,
        SPBLAS_TO_VECTOR_FLOAT_IN_OUT( x ),
        incx);
}


int BLAS_FLOAT_NAME(ussm)(
    enum blas_order_type order,
    enum blas_trans_type transt,
    int nrhs,
    SPBLAS_FLOAT_IN alpha,
    blas_sparse_matrix A,
    SPBLAS_VECTOR_FLOAT_IN_OUT b,
    int ldb )
{
  return BLAS_xussm(order, transt, nrhs,
      SPBLAS_TO_FLOAT_IN( alpha ), A,
      SPBLAS_TO_VECTOR_FLOAT_IN_OUT( b ), ldb);
}


/*  ----   end of FLOAT routines -------  */




 void BLAS_COMPLEX_FLOAT_NAME(usdot)(
    enum blas_conj_type conj_flag,
    int  nz,
    SPBLAS_VECTOR_COMPLEX_FLOAT_IN x,
    const int *index,
    SPBLAS_VECTOR_COMPLEX_FLOAT_IN y,
    int incy,
    SPBLAS_COMPLEX_FLOAT_OUT r,
    enum blas_base_type index_base)
{
   BLAS_xusdot(conj_flag, nz,
          SPBLAS_TO_VECTOR_COMPLEX_FLOAT_IN( x ), index,
          SPBLAS_TO_VECTOR_COMPLEX_FLOAT_IN( y ), incy,
          SPBLAS_TO_COMPLEX_FLOAT_OUT( r ), index_base);
}

 void BLAS_COMPLEX_FLOAT_NAME(usaxpy)(
    int nz,
    SPBLAS_COMPLEX_FLOAT_IN alpha,
    SPBLAS_VECTOR_COMPLEX_FLOAT_IN x,
    const int *index,
    SPBLAS_VECTOR_COMPLEX_FLOAT_IN_OUT y,
    int incy,
    enum blas_base_type index_base)
{
  BLAS_xusaxpy(nz, SPBLAS_TO_COMPLEX_FLOAT_IN( alpha ),
      SPBLAS_TO_VECTOR_COMPLEX_FLOAT_IN( x ), index,
      SPBLAS_TO_VECTOR_COMPLEX_FLOAT_IN_OUT( y ),
      incy, index_base);
}

 void BLAS_COMPLEX_FLOAT_NAME(usga)(
    int nz,
    SPBLAS_VECTOR_COMPLEX_FLOAT_IN y,
    int incy,
    SPBLAS_VECTOR_COMPLEX_FLOAT_IN_OUT x,
    const int *indx,
    enum blas_base_type index_base )
{
  BLAS_xusga( nz, SPBLAS_TO_VECTOR_COMPLEX_FLOAT_IN( y ), incy,
        SPBLAS_TO_VECTOR_COMPLEX_FLOAT_IN_OUT( x ), indx, index_base);

}

 void BLAS_COMPLEX_FLOAT_NAME(usgz)(
    int nz,
    SPBLAS_VECTOR_COMPLEX_FLOAT_IN_OUT y,
    int incy,
    SPBLAS_VECTOR_COMPLEX_FLOAT_OUT x,
    const int *indx,
    enum blas_base_type index_base )
{
  BLAS_xusgz(nz, SPBLAS_TO_COMPLEX_FLOAT_IN_OUT(y ), incy, SPBLAS_TO_COMPLEX_FLOAT_OUT(  x ),
    indx, index_base);
}


 void BLAS_COMPLEX_FLOAT_NAME(ussc)(
    int nz,
    SPBLAS_VECTOR_COMPLEX_FLOAT_IN x,
    SPBLAS_VECTOR_COMPLEX_FLOAT_IN_OUT y,
    int incy,
    const int *index,
    enum blas_base_type index_base)
{
  BLAS_xussc(nz, SPBLAS_TO_VECTOR_COMPLEX_FLOAT_IN( x ),
  SPBLAS_TO_COMPLEX_FLOAT_IN_OUT( y ), incy, index,
  index_base);
}

/*  COMPLEX_FLOAT Level 2/3 creation routines */

int BLAS_COMPLEX_FLOAT_NAME(uscr_begin)(int M, int N)
{
  TSp_mat<COMPLEX_FLOAT> *A = new TSp_mat<COMPLEX_FLOAT>(M, N);
  TSp_MAT_SET_COMPLEX_FLOAT(A);

  return Table_insert(A);
}

blas_sparse_matrix BLAS_COMPLEX_FLOAT_NAME(uscr_block_begin)(
    int Mb, int Nb, int k, int l )
{
  TSp_mat<COMPLEX_FLOAT> *A = new TSp_mat<COMPLEX_FLOAT>(Mb*k, Nb*l);

  TSp_MAT_SET_COMPLEX_FLOAT(A);
  A->set_const_block_parameters(Mb, Nb, k, l);

  return Table_insert(A);
}

blas_sparse_matrix BLAS_COMPLEX_FLOAT_NAME(uscr_variable_block_begin)(
    int Mb, int Nb, const int *k, const int *l )
{
  TSp_mat<COMPLEX_FLOAT> *A = new TSp_mat<COMPLEX_FLOAT>(
                std::accumulate(k, k+Mb, 0), std::accumulate(l, l+Nb, 0) );

  TSp_MAT_SET_COMPLEX_FLOAT(A);
  A->set_var_block_parameters(Mb, Nb, k, l);

  return Table_insert(A);

}

/*  COMPLEX_FLOAT Level 2/3 insertion routines */

int BLAS_COMPLEX_FLOAT_NAME(uscr_insert_entry)(
    blas_sparse_matrix A, SPBLAS_COMPLEX_FLOAT_IN val, int i, int j )
{
  return BLAS_xuscr_insert_entry(A, SPBLAS_TO_COMPLEX_FLOAT_IN( val ), i, j);
}

int BLAS_COMPLEX_FLOAT_NAME(uscr_insert_entries)(
    blas_sparse_matrix A, int nz,
    SPBLAS_VECTOR_COMPLEX_FLOAT_IN val,
    const int *indx, const int *jndx )
{
  return BLAS_xuscr_insert_entries(A, nz, SPBLAS_TO_VECTOR_COMPLEX_FLOAT_IN( val ), indx, jndx);
}

int BLAS_COMPLEX_FLOAT_NAME(uscr_insert_col)(
    blas_sparse_matrix A, int j, int nz,
    SPBLAS_VECTOR_COMPLEX_FLOAT_IN val, const int *indx )
{
  return BLAS_xuscr_insert_col(A, j, nz, SPBLAS_TO_VECTOR_COMPLEX_FLOAT_IN( val ), indx);
}

int BLAS_COMPLEX_FLOAT_NAME(uscr_insert_row)(
  blas_sparse_matrix A, int i, int nz,
  SPBLAS_VECTOR_COMPLEX_FLOAT_IN val, const int *indx );

int BLAS_COMPLEX_FLOAT_NAME(uscr_insert_clique)(
    blas_sparse_matrix A,
    const int k,
    const int l,
    SPBLAS_VECTOR_COMPLEX_FLOAT_IN val,
    const int row_stride,
    const int col_stride,
    const int *indx,
    const int *jndx );

int BLAS_COMPLEX_FLOAT_NAME(uscr_insert_block)(
    blas_sparse_matrix A,
    SPBLAS_VECTOR_COMPLEX_FLOAT_IN val,
    int row_stride,
    int col_stride,
    int i, int j )
{
  return BLAS_xuscr_insert_block(
        A, SPBLAS_TO_VECTOR_COMPLEX_FLOAT_IN( val ),
        row_stride, col_stride, i, j);
}

int BLAS_COMPLEX_FLOAT_NAME(uscr_end)(blas_sparse_matrix A)
{
  return BLAS_xuscr_end(A);
}

/*  COMPLEX_FLOAT Level 2/3 computational routines */

 int BLAS_COMPLEX_FLOAT_NAME(usmv)(enum
    blas_trans_type transa,
    SPBLAS_COMPLEX_FLOAT_IN alpha,
    blas_sparse_matrix A,
    SPBLAS_VECTOR_COMPLEX_FLOAT_IN x,
    int incx,
    SPBLAS_VECTOR_COMPLEX_FLOAT_IN_OUT y,
    int incy )
{
  return BLAS_xusmv(
      transa,   SPBLAS_TO_COMPLEX_FLOAT_IN( alpha ), A,
      SPBLAS_TO_VECTOR_COMPLEX_FLOAT_IN( x ), incx,
      SPBLAS_TO_VECTOR_COMPLEX_FLOAT_IN_OUT( y ), incy);
}

int BLAS_COMPLEX_FLOAT_NAME(usmm)(
    enum blas_order_type order,
    enum blas_trans_type transa,
    int nrhs,
    SPBLAS_COMPLEX_FLOAT_IN alpha,
    blas_sparse_matrix A,
    SPBLAS_VECTOR_COMPLEX_FLOAT_IN b,
    int ldb,
    SPBLAS_VECTOR_COMPLEX_FLOAT_IN_OUT c,
    int ldc )
{
  return BLAS_xusmm(
      order, transa, nrhs,
      SPBLAS_TO_COMPLEX_FLOAT_IN( alpha), A,
      SPBLAS_TO_VECTOR_COMPLEX_FLOAT_IN(b), ldb,
      SPBLAS_TO_VECTOR_COMPLEX_FLOAT_IN_OUT( c ), ldc);
}

int BLAS_COMPLEX_FLOAT_NAME(ussv)(
    enum blas_trans_type transa,
    SPBLAS_COMPLEX_FLOAT_IN alpha,
    blas_sparse_matrix A,
    SPBLAS_VECTOR_COMPLEX_FLOAT_IN_OUT x,
    int incx )
{
  return BLAS_xussv( transa,
        SPBLAS_TO_COMPLEX_FLOAT_IN( alpha ), A,
        SPBLAS_TO_VECTOR_COMPLEX_FLOAT_IN_OUT( x ),
        incx);
}


int BLAS_COMPLEX_FLOAT_NAME(ussm)(
    enum blas_order_type order,
    enum blas_trans_type transt,
    int nrhs,
    SPBLAS_COMPLEX_FLOAT_IN alpha,
    blas_sparse_matrix A,
    SPBLAS_VECTOR_COMPLEX_FLOAT_IN_OUT b,
    int ldb )
{
  return BLAS_xussm(order, transt, nrhs,
      SPBLAS_TO_COMPLEX_FLOAT_IN( alpha ), A,
      SPBLAS_TO_VECTOR_COMPLEX_FLOAT_IN_OUT( b ), ldb);
}

/*  ----   end of COMPLEX_COMPLEX_COMPLEX_FLOAT routines -------  */

/*  ----   end of COMPLEX_COMPLEX_COMPLEX_FLOAT routines -------  */





















#include "Manifolds/Element.h"

/*Define the namespace*/
namespace ROPTLIB{

    namespace GLOBAL{
        integer IZERO = 0, IONE = 1, ITWO = 2;
        realdp DZERO = 0, DONE = 1, DTWO = 2, DNONE = -1, DNTWO = -2;
        realdpcomplex ZZERO = { 0, 0 }, ZONE = { 1, 0 }, ZTWO = { 2, 0 }, ZNONE = { -1, 0 };
        char *N = const_cast<char *> ("N"), *T = const_cast<char *> ("T");
        char *L = const_cast<char *> ("L"), *R = const_cast<char *> ("R");
        char *V = const_cast<char *> ("V"), *C = const_cast<char *> ("C");
        char *U = const_cast<char *> ("U"), *A = const_cast<char *> ("A");
        char *S = const_cast<char *> ("S"), *O = const_cast<char *> ("O");
    };

    std::ostream &operator<<(std::ostream &output, const Element &x)
    {
        integer ls = x.Getls();
        const integer *size = x.Getsize();
        const realdp *Space = x.GetSpace();
        const integer length = x.Getlength();
        const integer *sharedtimes = x.GetSharedTimes();
        const bool iscomplex = x.Getiscomplex();
        
        integer product = 1;
        for (integer i = 2; i < ls; i++)
        {
            product *= size[i];
        }
        
        std::string type;
        if(x.Getiscomplex())
            type = "complex";
        else
            type = "real";
            
        if (Space == nullptr)
        {
            if (size == nullptr)
            {
                printf("it is an empty %s data with size 0", type.c_str());
            }
            else
            {
                if(iscomplex)
                    printf("it is an empty %s data with size %d", type.c_str(), size[0] / 2);
                else
                    printf("it is an empty %s data with size %d", type.c_str(), size[0]);
            }
            for (integer i = 1; i < ls; i++)
                printf(" x %d", size[i]);
            printf("\n");
        }
        else
            if (ls == 1 || (ls > 1 && size[1] * product == 1))
            {
                printf("Type:%s , shared times:%d, shared times address:%p\n", type.c_str(), *sharedtimes, sharedtimes);
                if(iscomplex)
                {
                    for (integer i = 0; i < length; i++, i++)
                        printf("%.10e + %.10e*i\n", Space[i], Space[i+1]);
                } else
                {
                    for (integer i = 0; i < length; i++)
                        printf("%.10e\n", Space[i]);
                }
            }
            else
                if (ls == 2 || product == 1)
                {
                    printf("Type:%s , shared times:%d, shared times address:%p\n", type.c_str(), *sharedtimes, sharedtimes);
                    if(iscomplex)
                    {
                        for (integer j = 0; j < size[0]; j++, j++)
                        {
                            for (integer k = 0; k < size[1]; k++)
                            {
                                printf("%.10e + %.10e*i\t", Space[j + size[0] * k], Space[j + 1 + size[0] * k]);
                            }
                            printf("\n");
                        }
                    } else
                    {
                        for (integer j = 0; j < size[0]; j++)
                        {
                            for (integer k = 0; k < size[1]; k++)
                            {
                                printf("%.10e\t", Space[j + size[0] * k]);
                            }
                            printf("\n");
                        }
                    }
                }
                else
                {
                    integer row = size[0], col = size[1];
                    integer *idices = new integer[ls + 1];
                    const realdp *ptr = Space;
                    for (integer i = 2; i < ls + 1; i++)
                        idices[i] = 0;
                    while (1)
                    {
                        printf("Type: %s (:,:", type.c_str());
                        for (integer i = 2; i < ls; i++)
                            printf(",%d", idices[i]);
                        printf("), shared times:%d\n", *sharedtimes);
                        
                        if(iscomplex)
                        {
                            for (integer j = 0; j < row; j++, j++)
                            {
                                for (integer k = 0; k < col; k++)
                                {
                                    printf("%.10e + %.10e*i\t", ptr[j + row * k], ptr[j + 1 + row * k]);
                                }
                                printf("\n");
                            }
                        }
                        else
                        {
                            for (integer j = 0; j < row; j++)
                            {
                                for (integer k = 0; k < col; k++)
                                {
                                    printf("%.10e\t", ptr[j + row * k]);
                                }
                                printf("\n");
                            }
                        }
                        ptr += row * col;
                        idices[2]++;
                        for (integer i = 2; i < ls; i++)
                        {
                            if (idices[i] == size[i])
                            {
                                idices[i] = 0;
                                idices[i + 1]++;
                            }
                        }
                        if (idices[ls] == 1)
                            break;
                    }
                    delete[] idices;
                }
        
        return output;
    };

    Element operator+(Element left, Element right)
    {
        assert(left.Getlength() == right.Getlength());
        Element result(right);
        realdp *resultptr = result.ObtainWritePartialData();
        integer length = left.Getlength();
        const realdp *leftptr = left.ObtainReadData();
        axpy_(&length, &GLOBAL::DONE, const_cast<realdp *> (leftptr), &GLOBAL::IONE, resultptr, &GLOBAL::IONE);
        
        return result;
    };

    Element operator+(Element left, realdp right)
    {
        assert(!left.Getiscomplex());
        Element result(left);
        realdp *resultptr = result.ObtainWriteEntireData();
        const realdp *leftptr = left.ObtainReadData();
        for(integer i = 0; i < left.Getlength(); i++)
            resultptr[i] = leftptr[i] + right;

        return result;
    };

    Element operator+(Element left, realdpcomplex right)
    {
        assert(left.Getiscomplex());
        Element result(left);
        realdp *resultptr = result.ObtainWriteEntireData();
        const realdp *leftptr = left.ObtainReadData();
        for(integer i = 0; i < left.Getlength(); i++, i++)
        {
            resultptr[i] = leftptr[i] + right.r;
            resultptr[i + 1] = leftptr[i + 1] + right.i;
        }

        return result;
    };

    Element operator+(realdp left, Element right)
    {
        assert(!right.Getiscomplex());
        Element result(right);
        realdp *resultptr = result.ObtainWriteEntireData();
        const realdp *rightptr = right.ObtainReadData();
        for(integer i = 0; i < right.Getlength(); i++)
            resultptr[i] = rightptr[i] + left;

        return result;
    };

    Element operator+(realdpcomplex left, Element right)
    {
        assert(right.Getiscomplex());
        Element result(right);
        realdp *resultptr = result.ObtainWriteEntireData();
        const realdp *rightptr = right.ObtainReadData();
        for(integer i = 0; i < right.Getlength(); i++, i++)
        {
            resultptr[i] = rightptr[i] + left.r;
            resultptr[i + 1] = rightptr[i + 1] + left.i;
        }

        return result;
    };

    Element operator+=(Element left, Element right)
    {
        assert(left.Getlength() == right.Getlength());
        realdp *leftptr = left.ObtainWritePartialData();
        const realdp *rightptr = right.ObtainReadData();
        integer length = left.Getlength();
        axpy_(&length, &GLOBAL::DONE, const_cast<realdp *> (rightptr), &GLOBAL::IONE, leftptr, &GLOBAL::IONE);
        return left;
    };

    Element operator-(Element left, Element right)
    {
        assert(left.Getlength() == right.Getlength());
        Element result(left);
        realdp *resultptr = result.ObtainWritePartialData();
        integer length = left.Getlength();
        const realdp *rightptr = right.ObtainReadData();
        axpy_(&length, &GLOBAL::DNONE, const_cast<realdp *> (rightptr), &GLOBAL::IONE, resultptr, &GLOBAL::IONE);
        
        return result;
    };

    Element operator-(Element left, realdp right)
    {
        assert(!left.Getiscomplex());
        Element result(left);
        realdp *resultptr = result.ObtainWriteEntireData();
        const realdp *leftptr = left.ObtainReadData();
        for(integer i = 0; i < left.Getlength(); i++)
            resultptr[i] = leftptr[i] - right;

        return result;
    };

    Element operator-(Element left, realdpcomplex right)
    {
        assert(left.Getiscomplex());
        Element result(left);
        realdp *resultptr = result.ObtainWriteEntireData();
        const realdp *leftptr = left.ObtainReadData();
        for(integer i = 0; i < left.Getlength(); i++, i++)
        {
            resultptr[i] = leftptr[i] - right.r;
            resultptr[i + 1] = leftptr[i + 1] - right.i;
        }

        return result;
    };

    Element operator-(realdp left, Element right)
    {
        assert(!right.Getiscomplex());
        Element result(right);
        realdp *resultptr = result.ObtainWriteEntireData();
        const realdp *rightptr = right.ObtainReadData();
        for(integer i = 0; i < right.Getlength(); i++)
            resultptr[i] = left - rightptr[i];

        return result;
    };

    Element operator-(realdpcomplex left, Element right)
    {
        assert(right.Getiscomplex());
        Element result(right);
        realdp *resultptr = result.ObtainWriteEntireData();
        const realdp *rightptr = right.ObtainReadData();
        for(integer i = 0; i < right.Getlength(); i++, i++)
        {
            resultptr[i] = left.r - rightptr[i];
            resultptr[i + 1] = left.i - rightptr[i + 1];
        }

        return result;
    };

    Element operator-=(Element left, Element right)
    {
        assert(left.Getlength() == right.Getlength());
        realdp *leftptr = left.ObtainWritePartialData();
        const realdp *rightptr = right.ObtainReadData();
        integer length = left.Getlength();
        axpy_(&length, &GLOBAL::DNONE, const_cast<realdp *> (rightptr), &GLOBAL::IONE, leftptr, &GLOBAL::IONE);
        return left;
    };

    Element operator*(Element left, Element right)
    {
        assert(left.Getiscomplex() == right.Getiscomplex());
        
        if(left.Getiscomplex())
        {
            integer m1 = left.Getsize()[0] / 2, n1 = left.Getsize()[1];
            integer m2 = right.Getsize()[0] / 2, n2 = right.Getsize()[1];
            assert(n1 == m2);
            Element result(m1, n2, 1, "complex");
            const realdpcomplex *leftptr = (realdpcomplex *) left.ObtainReadData();
            const realdpcomplex *rightptr = (realdpcomplex *) right.ObtainReadData();
            realdpcomplex *resultptr = (realdpcomplex *) result.ObtainWriteEntireData();
            
            gemm_(GLOBAL::N, GLOBAL::N, &m1, &n2, &n1, &GLOBAL::ZONE, const_cast<realdpcomplex *> (leftptr), &m1, const_cast<realdpcomplex *> (rightptr), &m2, &GLOBAL::ZZERO, resultptr, &m1);
            
            return result;
        }
        
        integer m1 = left.Getsize()[0], n1 = left.Getsize()[1];
        integer m2 = right.Getsize()[0], n2 = right.Getsize()[1];
        assert(n1 == m2);
        Element result(m1, n2);
        const realdp *leftptr = left.ObtainReadData();
        const realdp *rightptr = right.ObtainReadData();
        realdp *resultptr = result.ObtainWriteEntireData();
        
        gemm_(GLOBAL::N, GLOBAL::N, &m1, &n2, &n1, &GLOBAL::DONE, const_cast<realdp *> (leftptr), &m1, const_cast<realdp *> (rightptr), &m2, &GLOBAL::DZERO, resultptr, &m1);
        
        return result;
    };

    Element operator*(realdp value, Element mat)
    {
        Element result(mat);
        const realdp *matptr = mat.ObtainReadData();
        realdp *resultptr = result.ObtainWriteEntireData();
        
        for(integer i = 0; i < mat.Getlength(); i++)
            resultptr[i] = value * matptr[i];

        return result;
    };

    Element operator*(realdpcomplex value, Element mat)
    {
        assert(mat.Getiscomplex());
        Element result(mat);
        const realdp *matptr = mat.ObtainReadData();
        realdp *resultptr = result.ObtainWriteEntireData();
        
        for(integer i = 0; i < mat.Getlength(); i++, i++)
        {
            resultptr[i] = value.r * matptr[i] - value.i * matptr[i + 1];
            resultptr[i + 1] = value.r * matptr[i + 1] + value.i * matptr[i];
        }

        return result;
    };

    Element operator*(Element mat, realdp value)
    {
        Element result(mat);
        const realdp *matptr = mat.ObtainReadData();
        realdp *resultptr = result.ObtainWriteEntireData();
        
        for(integer i = 0; i < mat.Getlength(); i++)
            resultptr[i] = value * matptr[i];

        return result;
    };

    Element operator*(Element mat, realdpcomplex value)
    {
        assert(mat.Getiscomplex());
        Element result(mat);
        const realdp *matptr = mat.ObtainReadData();
        realdp *resultptr = result.ObtainWriteEntireData();
        
        for(integer i = 0; i < mat.Getlength(); i++, i++)
        {
            resultptr[i] = value.r * matptr[i] - value.i * matptr[i + 1];
            resultptr[i + 1] = value.r * matptr[i + 1] + value.i * matptr[i];
        }

        return result;
    };

    Element operator/(Element left, Element right)
    {
        assert(left.Getiscomplex() == right.Getiscomplex());
        
        right.LUdecom();
        if(!left.Getiscomplex())
        {
            assert(left.Getsize()[1] == right.Getsize()[0] && right.Getsize()[0] == right.Getsize()[1]);
            integer m = left.Getsize()[0];
            integer n = left.Getsize()[1];
            
            Element LU = right.Field("_LU");
            Element P = right.Field("_P");
            const realdp *LUptr = LU.ObtainReadData();
            const integer *Pptr = (integer *) P.ObtainReadData();
            
            Element lefttran;
            lefttran = left.GetTranspose();
            realdp *lefttranptr = lefttran.ObtainWritePartialData();
            integer info = 0;
            
            /* solve linear system: PMGQ * X = v using the LU decomposition results from getrf, then solution is stored in v.
             details: www.netlib.org/lapack/explore-html/d6/d49/getrs_8f.html */
            getrs_(GLOBAL::T, &n, &m, const_cast<realdp *>(LUptr), &n, const_cast<integer *> (Pptr), lefttranptr, &n, &info);
            
            lefttran.Transpose();
            return lefttran;
        }

        assert(left.Getsize()[1] == right.Getsize()[0]/2 && right.Getsize()[0]/2 == right.Getsize()[1]);
        integer m = left.Getsize()[0] / 2;
        integer n = left.Getsize()[1];
        Element LU = right.Field("_LU");
        Element P = right.Field("_P");
        const realdpcomplex *LUptr = (realdpcomplex *) LU.ObtainReadData();
        const integer *Pptr = (integer *) P.ObtainReadData();
        
        Element lefttran;
        lefttran = left.GetTranspose();
        realdpcomplex *lefttranptr = (realdpcomplex *) lefttran.ObtainWritePartialData();
        integer info = 0;
        
        /* solve linear system: PMGQ * X = v using the LU decomposition results from getrf, then solution is stored in v.
        details: www.netlib.org/lapack/explore-html/d6/d49/getrs_8f.html */
        getrs_(GLOBAL::C, &n, &m, const_cast<realdpcomplex *>(LUptr), &n, const_cast<integer *> (Pptr), lefttranptr, &n, &info);
        
        lefttran.Transpose();
        
        return lefttran;
    };

    Element operator%(Element left, Element right)
    {
        assert(left.Getiscomplex() == right.Getiscomplex());
        
        left.LUdecom();
        if(!left.Getiscomplex())
        {
            assert(left.Getsize()[0] == left.Getsize()[1] && left.Getsize()[1] == right.Getsize()[0]);
            integer m = left.Getsize()[0];
            integer n = right.Getsize()[1];
            
            Element LU = left.Field("_LU");
            Element P = left.Field("_P");
            const realdp *LUptr = LU.ObtainReadData();
            const integer *Pptr = (integer *) P.ObtainReadData();
            
            Element result(right);
            realdp *resultptr = result.ObtainWritePartialData();
            integer info = 0;
            
            /* solve linear system: PMGQ * X = v using the LU decomposition results from getrf, then solution is stored in v.
            details: www.netlib.org/lapack/explore-html/d6/d49/getrs_8f.html */
            getrs_(GLOBAL::N, &m, &n, const_cast<realdp *>(LUptr), &m, const_cast<integer *> (Pptr), resultptr, &m, &info);
            
            return result;
        }

        assert(left.Getsize()[0]/2 == left.Getsize()[1] && left.Getsize()[1] == right.Getsize()[0]/2);
        
        integer m = left.Getsize()[0] / 2;
        integer n = right.Getsize()[1];
        Element LU = left.Field("_LU");
        Element P = left.Field("_P");
        const realdpcomplex *LUptr = (realdpcomplex *) LU.ObtainReadData();
        const integer *Pptr = (integer *) P.ObtainReadData();
        
        Element result(right);
        realdpcomplex *resultptr = (realdpcomplex *) result.ObtainWritePartialData();
        integer info = 0;
        
        /* solve linear system: PMGQ * X = v using the LU decomposition results from getrf, then solution is stored in v.
         details: www.netlib.org/lapack/explore-html/d6/d49/getrs_8f.html */
        getrs_(GLOBAL::N, &m, &n, const_cast<realdpcomplex *>(LUptr), &m, const_cast<integer *> (Pptr), resultptr, &m, &info);
        
        return result;
    };

    Element operator/(Element mat, realdp value)
    {
        Element result(mat);
        const realdp *matptr = mat.ObtainReadData();
        realdp *resultptr = result.ObtainWriteEntireData();
        
        for(integer i = 0; i < mat.Getlength(); i++)
            resultptr[i] = matptr[i] / value;

        return result;
    };

    Element operator/(Element mat, realdpcomplex value)
    {
        assert(mat.Getiscomplex());
        Element result(mat);
        const realdp *matptr = mat.ObtainReadData();
        realdp *resultptr = result.ObtainWriteEntireData();
        realdp tmp = value.r * value.r + value.i * value.i;
        
        for(integer i = 0; i < mat.Getlength(); i++, i++)
        {
            resultptr[i] = (value.r * matptr[i] + value.i * matptr[i + 1]) / tmp;
            resultptr[i + 1] = (value.r * matptr[i + 1] - value.i * matptr[i]) / tmp;
        }

        return result;
    };

    Element operator/(realdp value, Element mat)
    {
        assert(!mat.Getiscomplex());
        Element result(mat);
        const realdp *matptr = mat.ObtainReadData();
        realdp *resultptr = result.ObtainWriteEntireData();
        
        for(integer i = 0; i < mat.Getlength(); i++)
            resultptr[i] = value / matptr[i];

        return result;
    };

    Element operator/(realdpcomplex value, Element mat)
    {
        assert(mat.Getiscomplex());
        Element result(mat);
        const realdp *matptr = mat.ObtainReadData();
        realdp *resultptr = result.ObtainWriteEntireData();
        realdp tmp = 0;
        
        for(integer i = 0; i < mat.Getlength(); i++, i++)
        {
            tmp = matptr[i] * matptr[i] + matptr[i + 1] * matptr[i + 1];
            resultptr[i] = (value.r * matptr[i] + value.i * matptr[i + 1]) / tmp;
            resultptr[i + 1] = (-value.r * matptr[i + 1] + value.i * matptr[i]) / tmp;
        }

        return result;
    };

    Element operator==(realdp value, Element mat)
    {
        Element result(mat);
        const realdp *matptr = mat.ObtainReadData();
        realdp *resultptr = result.ObtainWriteEntireData();
        for(integer i = 0; i < mat.Getlength(); i++)
            resultptr[i] = (value == matptr[i]);
        return result;
    };

    Element operator==(Element mat, realdp value)
    {
        Element result(mat);
        const realdp *matptr = mat.ObtainReadData();
        realdp *resultptr = result.ObtainWriteEntireData();
        for(integer i = 0; i < mat.Getlength(); i++)
            resultptr[i] = (value == matptr[i]);
        return result;
    };

    Element operator==(Element left, Element right)
    {
        assert(left.Getlength() == right.Getlength());
        
        const realdp *leftptr = left.ObtainReadData();
        const realdp *rightptr = right.ObtainReadData();
        
        Element result(left);
        realdp *resultptr = result.ObtainWriteEntireData();
        for(integer i = 0; i < left.Getlength(); i++)
            resultptr[i] = (leftptr[i] == rightptr[i]);
        return result;
    };

    Element operator>(realdp value, Element mat)
    {
        Element result(mat);
        const realdp *matptr = mat.ObtainReadData();
        realdp *resultptr = result.ObtainWriteEntireData();
        for(integer i = 0; i < mat.Getlength(); i++)
            resultptr[i] = (value > matptr[i]);
        return result;
    };

    Element operator>(Element mat, realdp value)
    {
        Element result(mat);
        const realdp *matptr = mat.ObtainReadData();
        realdp *resultptr = result.ObtainWriteEntireData();
        for(integer i = 0; i < mat.Getlength(); i++)
            resultptr[i] = (value < matptr[i]);
        return result;
    };

    Element operator>(Element left, Element right)
    {
        assert(left.Getlength() == right.Getlength());
        
        const realdp *leftptr = left.ObtainReadData();
        const realdp *rightptr = right.ObtainReadData();
        
        Element result(left);
        realdp *resultptr = result.ObtainWriteEntireData();
        for(integer i = 0; i < left.Getlength(); i++)
            resultptr[i] = (leftptr[i] > rightptr[i]);
        return result;
    };

    Element operator>=(realdp value, Element mat)
    {
        Element result(mat);
        const realdp *matptr = mat.ObtainReadData();
        realdp *resultptr = result.ObtainWriteEntireData();
        for(integer i = 0; i < mat.Getlength(); i++)
            resultptr[i] = (value >= matptr[i]);
        return result;
    };

    Element operator>=(Element mat, realdp value)
    {
        Element result(mat);
        const realdp *matptr = mat.ObtainReadData();
        realdp *resultptr = result.ObtainWriteEntireData();
        for(integer i = 0; i < mat.Getlength(); i++)
            resultptr[i] = (value <= matptr[i]);
        return result;
    };

    Element operator>=(Element left, Element right)
    {
        assert(left.Getlength() == right.Getlength());
        
        const realdp *leftptr = left.ObtainReadData();
        const realdp *rightptr = right.ObtainReadData();
        
        Element result(left);
        realdp *resultptr = result.ObtainWriteEntireData();
        for(integer i = 0; i < left.Getlength(); i++)
            resultptr[i] = (leftptr[i] >= rightptr[i]);
        return result;
    };

    Element operator<(realdp value, Element mat)
    {
        Element result(mat);
        const realdp *matptr = mat.ObtainReadData();
        realdp *resultptr = result.ObtainWriteEntireData();
        for(integer i = 0; i < mat.Getlength(); i++)
            resultptr[i] = (value < matptr[i]);
        return result;
    };

    Element operator<(Element mat, realdp value)
    {
        Element result(mat);
        const realdp *matptr = mat.ObtainReadData();
        realdp *resultptr = result.ObtainWriteEntireData();
        for(integer i = 0; i < mat.Getlength(); i++)
            resultptr[i] = (value > matptr[i]);
        return result;
    };

    Element operator<(Element left, Element right)
    {
        assert(left.Getlength() == right.Getlength());
        
        const realdp *leftptr = left.ObtainReadData();
        const realdp *rightptr = right.ObtainReadData();
        
        Element result(left);
        realdp *resultptr = result.ObtainWriteEntireData();
        for(integer i = 0; i < left.Getlength(); i++)
            resultptr[i] = (leftptr[i] < rightptr[i]);
        return result;
    };

    Element operator<=(realdp value, Element mat)
    {
        Element result(mat);
        const realdp *matptr = mat.ObtainReadData();
        realdp *resultptr = result.ObtainWriteEntireData();
        for(integer i = 0; i < mat.Getlength(); i++)
            resultptr[i] = (value <= matptr[i]);
        return result;
    };

    Element operator<=(Element mat, realdp value)
    {
        Element result(mat);
        const realdp *matptr = mat.ObtainReadData();
        realdp *resultptr = result.ObtainWriteEntireData();
        for(integer i = 0; i < mat.Getlength(); i++)
            resultptr[i] = (value >= matptr[i]);
        return result;
    };

    Element operator<=(Element left, Element right)
    {
        assert(left.Getlength() == right.Getlength());
        
        const realdp *leftptr = left.ObtainReadData();
        const realdp *rightptr = right.ObtainReadData();
        
        Element result(left);
        realdp *resultptr = result.ObtainWriteEntireData();
        for(integer i = 0; i < left.Getlength(); i++)
            resultptr[i] = (leftptr[i] <= rightptr[i]);
        return result;
    };

    Element operator*(const SparseMatrix &left, Element right)
    {
        NIST_SPBLAS::Sp_mat *SM = NIST_SPBLAS::GetSpMFromTable(left.GetSparseM());
        integer leftrow = SM->num_rows();
        integer leftcol = SM->num_cols();
        bool leftiscomplex = SM->is_complex();
        
        assert(leftcol == right.Getrow() && leftiscomplex == right.Getiscomplex());
        if(!leftiscomplex)
        {
            Element result(leftrow, right.Getcol());
            result.SetToZeros();
            realdp *resultptr = result.ObtainWritePartialData();
            const realdp *rightptr = right.ObtainReadData();
            BLAS_usmm(blas_colmajor, blas_no_trans, right.Getcol(), 1, left.GetSparseM(), rightptr, right.Getrow(), resultptr, leftrow);
            return result;
        }
        
        Element result(leftrow, right.Getcol(), "complex");
        result.SetToZeros();
        realdp *resultptr = result.ObtainWritePartialData();
        const realdp *rightptr = right.ObtainReadData();
        BLAS_usmm(blas_colmajor, blas_no_trans, right.Getcol(), &GLOBAL::ZONE, left.GetSparseM(), (realdpcomplex *) rightptr, right.Getrow(), (realdpcomplex *) resultptr, leftrow);
        return result;
    };

    Element operator*(Element left, const SparseMatrix &right)
    {
        NIST_SPBLAS::Sp_mat *SM = NIST_SPBLAS::GetSpMFromTable(right.GetSparseM());
        
        integer rightrow = SM->num_rows();
        integer rightcol = SM->num_cols();
        bool rightiscomplex = SM->is_complex();
        assert(left.Getcol() == rightrow && left.Getiscomplex() == rightiscomplex);
        Element lefttrans = left.GetTranspose();
        if(!rightiscomplex)
        { /*right^T * lefttrans*/
            Element result(rightcol, lefttrans.Getcol());
            result.SetToZeros();
            realdp *resultptr = result.ObtainWritePartialData();
            const realdp *lefttransptr = lefttrans.ObtainReadData();
            BLAS_usmm(blas_colmajor, blas_trans, lefttrans.Getcol(), 1, right.GetSparseM(), lefttransptr, lefttrans.Getrow(), resultptr, rightcol);
            return result.GetTranspose();
        }

        Element result(rightcol, lefttrans.Getcol(), "complex");
        result.SetToZeros();
        realdp *resultptr = result.ObtainWritePartialData();
        const realdp *lefttransptr = lefttrans.ObtainReadData();
        BLAS_usmm(blas_colmajor, blas_conj_trans, lefttrans.Getcol(), &GLOBAL::ZONE, right.GetSparseM(), (realdpcomplex *) lefttransptr, lefttrans.Getrow(), (realdpcomplex *) resultptr, rightcol);
        return result.GetTranspose();
    };

    void Element::LUdecom() const
    {
        if(FieldsExist("_LU") && FieldsExist("_P"))
            return;
        
        if(!iscomplex)
        {
            assert(ls >= 2 && (size[0] == size[1]));
            integer k = size[0];
            Element LU(*this);
            realdp *LUptr = LU.ObtainWritePartialData();
            integer info = 0;
            Element P(k);
            integer *Pptr = (integer *) P.ObtainWriteEntireData();
            
            /*Lapack function for LU decomposition*/
            getrf_(&k, &k, LUptr, &k, Pptr, &info);
            
            this->AddToFields("_LU", LU);
            this->AddToFields("_P", P);
            return;
        }

        assert(ls >= 2 && (size[0] / 2 == size[1]));
        integer k = size[0] / 2;
        Element LU(*this);
        realdpcomplex *LUptr = (realdpcomplex *) LU.ObtainWritePartialData();
        integer info = 0;
        Element P(k);
        integer *Pptr = (integer *) P.ObtainWriteEntireData();
        
        /*Lapack function for LU decomposition*/
        getrf_(&k, &k, LUptr, &k, Pptr, &info);
        
        this->AddToFields("_LU", LU);
        this->AddToFields("_P", P);
    };

    void Element::CholDecom(void) const
    {
        if(FieldsExist("_L"))
            return;
        
        if(!iscomplex)
        {
            assert(ls >= 2 && (size[0] == size[1]));
            integer k = size[0];
            Element L(*this);
            realdp *Lptr = L.ObtainWriteEntireData();
            const realdp *ptr = this->ObtainReadData();
            for(integer i = 0; i < k; i++)
            {
                for(integer j = 0; j < k; j++)
                {
                    Lptr[j + i * k] = (j >= i) ? ptr[j + i * k] : 0;
                }
            }
            integer info = 0;
            potrf_(GLOBAL::L, &k, Lptr, &k, &info);
            if(info != 0)
            {
                printf("The Cholesky decomposition in Element::CholDecom fails with info: %d\n", info);
            }
            this->AddToFields("_L", L);
            return;
        }

        assert(ls >= 2 && (size[0] / 2 == size[1]));
        integer k = size[0] / 2;
        Element L(*this);
        realdp *Lptr = L.ObtainWriteEntireData();
        const realdp *ptr = this->ObtainReadData();
        for(integer i = 0; i < k; i++)
        {
            for(integer j = 0; j < k; j++)
            {
                Lptr[2 * j + i * 2 * k] = (j >= i) ? ptr[2 * j + i * 2 * k] : 0;
                Lptr[2 * j + 1 + i * 2 * k] = (j >= i) ? ptr[2 * j + 1 + i * 2 * k] : 0;
            }
        }
        integer info = 0;
        potrf_(GLOBAL::L, &k, (realdpcomplex *) Lptr, &k, &info);
        if(info != 0)
        {
            printf("The Cholesky decomposition in Element::CholDecom fails with info: %d\n", info);
        }
        this->AddToFields("_L", L);
    };
            
    Element Element::TriangleLinSol(Element L, char *trans) const
    {
        /* Solver a linear system L^* X = B, where L is a lower triangle matrix.
        trans = "N": solve L X = B
        trans = "T": solve L^T X = B
        trans = "C": solve L^H X = B */

        assert(L.Getrow() == L.Getcol() && Getrow() == L.Getcol() && iscomplex == L.Getiscomplex());
        integer k = Getrow(), n = Getcol();
        
        if(!iscomplex)
        {
            const realdp *Lptr = L.ObtainReadData();
            Element result(*this);
            realdp *resultptr = result.ObtainWritePartialData();
            integer info = 0;
            trtrs_(GLOBAL::L, trans, GLOBAL::N, &k, &n, const_cast<realdp *> (Lptr), &k, resultptr, &k, &info);
            if(info != 0)
            {
                printf("Solving linear system in Element::TriangleLinSol fails with info: %d\n", info);
            }
            return result;
        }

        const realdp *Lptr = L.ObtainReadData();
        Element result(*this);
        realdp *resultptr = result.ObtainWritePartialData();
        integer info = 0;
        trtrs_(GLOBAL::L, trans, GLOBAL::N, &k, &n, const_cast<realdpcomplex *> ((realdpcomplex *) Lptr), &k, (realdpcomplex *) resultptr, &k, &info);
        if(info != 0)
        {
            printf("Solving linear system in Element::TriangleLinSol fails with info: %d\n", info);
        }
        return result;
    };

    void Element::EigenDecomSym(char *JobZ, char *UorL)
    {
        if(FieldsExist("_EigVal") && (strcmp(JobZ, "N") == 0))
            return;
        
        if(FieldsExist("_EigVal") && FieldsExist("_EigVec"))
            return;
        
        if(!iscomplex)
        {
            assert(ls >= 2 && size[0] == size[1]);
            /* eigenvalue decomposition: syev approach. */
            integer N = size[0];
            
            Element EigVec(*this);
            realdp *EigVecptr = EigVec.ObtainWritePartialData();
            Element EigVal(N);
            realdp *EigValptr = EigVal.ObtainWriteEntireData();
            
            integer lwork = -1, info;
            realdp lworkopt;
            /*Obtain the size of memory required for eigenvalue decomposition*/
            syev_(JobZ, UorL, &N, EigVecptr, &N, EigValptr, &lworkopt, &lwork, &info);
            /*allocate the desired memory*/
            lwork = static_cast<integer> (lworkopt);
            realdp *work = new realdp[lwork];
            syev_(JobZ, UorL, &N, EigVecptr, &N, EigValptr, work, &lwork, &info);
            delete[] work;
            this->AddToFields("_EigVal", EigVal);
            if(strcmp(JobZ, "V") == 0) /*if JobZ == "V"*/
                this->AddToFields("_EigVec", EigVec);
            
            return;
        }
        assert(ls >= 2 && size[0] / 2 == size[1]);
        /* eigenvalue decomposition: syev approach. */
        integer N = size[0] / 2;
        
        Element EigVec(*this);
        realdpcomplex *EigVecptr = (realdpcomplex *) EigVec.ObtainWritePartialData();
        Element EigVal(N, 1, 1, "complex");
        realdp *EigValptr = EigVal.ObtainWriteEntireData();
        
        integer lwork = -1, info;
        realdpcomplex lworkopt;
        realdp *rwork = new realdp[(1 < 3 * N - 2) ? 3 * N - 2 : 1];
        /*Obtain the size of memory required for eigenvalue decomposition*/
        heev_(JobZ, UorL, &N, EigVecptr, &N, EigValptr, &lworkopt, &lwork, rwork, &info);
        /*allocate the desired memory*/
        lwork = static_cast<integer> (lworkopt.r);
        realdpcomplex *work = new realdpcomplex[lwork];
        heev_(JobZ, UorL, &N, EigVecptr, &N, EigValptr, work, &lwork, rwork, &info);
        delete[] work;
        delete[] rwork;
        
        for(integer i = N-1; i >= 0; i--)
        {
            EigValptr[2 * i + 1] = 0;
            EigValptr[2 * i] = EigValptr[i];
        }
        
        this->AddToFields("_EigVal", EigVal);
        if(strcmp(JobZ, "V") == 0) /*if JobZ == "V"*/
            this->AddToFields("_EigVec", EigVec);

        return;
    };

    Element Element::ExpSym(char *UorL)
    {
        if(!iscomplex)
        {
            assert(ls >= 2 && size[0] == size[1]);
            integer N = size[0];
            EigenDecomSym(GLOBAL::V, UorL);
            
            Element EigVec = this->Field("_EigVec");
            Element EigVal = this->Field("_EigVal");
            const realdp *EigValptr = EigVal.ObtainReadData();
            Element Tmp(EigVec);
            realdp *Tmpptr = Tmp.ObtainWritePartialData();
            
            for(integer i = 0; i < N; i++)
            {
                realdp a = exp(EigValptr[i]);
                scal_(&N, &a, Tmpptr + i * N, &GLOBAL::IONE);
            }
            Element result = Tmp * EigVec.GetTranspose();
            return result;
        }
        assert(ls >= 2 && size[0] / 2 == size[1]);
        integer N = size[0] / 2;
        EigenDecomSym(GLOBAL::V, UorL);
        Element EigVec = this->Field("_EigVec");
        Element EigVal = this->Field("_EigVal");
        const realdpcomplex *EigValptr = (realdpcomplex *) EigVal.ObtainReadData();
        
        Element Tmp(EigVec);
        realdp *Tmpptr = Tmp.ObtainWritePartialData();
        
        for(integer i = 0; i < N; i++)
        {
            realdp a = exp(EigValptr[i].r);
            integer N2 = N * 2;
            scal_(&N2, &a, Tmpptr + i * 2 * N, &GLOBAL::IONE);
        }
        
        Element result = Tmp * EigVec.GetTranspose();
        return result;
    };


    Element Element::LogSym(char *UorL)
    {
        if(!iscomplex)
        {
            assert(ls >= 2 && size[0] == size[1]);
            integer N = size[0];
            EigenDecomSym(GLOBAL::V, UorL);
            
            Element EigVec = this->Field("_EigVec");
            Element EigVal = this->Field("_EigVal");
            const realdp *EigValptr = EigVal.ObtainReadData();
            Element Tmp(EigVec);
            realdp *Tmpptr = Tmp.ObtainWritePartialData();
            realdp a = 0;
            for(integer i = 0; i < N; i++)
            {
                if(EigValptr[i] <= 0)
                {
                    printf("Error: Eigenvalues have negative values. Log can not be computed!\n");
                    return (*this);
                }
                a = log(EigValptr[i]);
                scal_(&N, &a, Tmpptr + i * N, &GLOBAL::IONE);
            }
            Element result = Tmp * EigVec.GetTranspose();
            return result;
        }
        assert(ls >= 2 && size[0] / 2 == size[1]);
        integer N = size[0] / 2;
        EigenDecomSym(GLOBAL::V, UorL);
        Element EigVec = this->Field("_EigVec");
        Element EigVal = this->Field("_EigVal");
        const realdpcomplex *EigValptr = (realdpcomplex *) EigVal.ObtainReadData();
        
        Element Tmp(EigVec);
        realdp *Tmpptr = Tmp.ObtainWritePartialData();
        realdp a = 0;
        
        for(integer i = 0; i < N; i++)
        {
            if(EigValptr[i].r <= 0)
            {
                printf("Error: Eigenvalues have negative values. Log can not be computed!\n");
                return (*this);
            }
            a = log(EigValptr[i].r);
            integer N2 = N * 2;
            scal_(&N2, &a, Tmpptr + i * 2 * N, &GLOBAL::IONE);
        }
        
        Element result = Tmp * EigVec.GetTranspose();
        return result;
    };

    void Element::HHRDecom(void) const
    {
        if(FieldsExist("_HHR") && FieldsExist("_tau"))
            return;
        
        if(!iscomplex)
        {
            integer m = size[0], n = size[1], minmn = (m < n) ? m : n;
            Element HHR(*this), tau(minmn);
            realdp *HHRptr = HHR.ObtainWritePartialData();
            realdp *tauptr = tau.ObtainWriteEntireData();
            integer *jpvt = new integer[n];
            integer info;
            integer lwork = -1;
            realdp lworkopt;
            for (integer i = 0; i < n; i++)
                jpvt[i] = i + 1;
            /*  compute the space required in the geqp3 */
            geqp3_(&m, &n, HHRptr, &m, jpvt, tauptr, &lworkopt, &lwork, &info);
            lwork = static_cast<integer> (lworkopt);
            realdp *work = new realdp[lwork];
            /* QR decomposition for ptrHHR using Householder reflections. Householder reflectors and R are stored in ptrHHR.
            details: www.netlib.org/lapack/explore-html/db/de5/geqp3_8f.html */
            geqp3_(&m, &n, HHRptr, &m, jpvt, tauptr, work, &lwork, &info);
            if (info < 0)
                printf("Error in qr decomposition!\n");

            for (integer i = 0; i < n; i++)
            {
                if (jpvt[i] != (i + 1))
                    printf("Error in qr decomposition!\n");
            }
            AddToFields("_HHR", HHR);
            AddToFields("_tau", tau);
            delete[] jpvt;
            delete[] work;
            return;
        }
        
        integer m = size[0] / 2, n = size[1], minmn = (m < n) ? m : n;
        Element HHR(*this), tau(minmn, 1, 1, "complex");
        realdpcomplex *HHRptr = (realdpcomplex *) HHR.ObtainWritePartialData();
        realdpcomplex *tauptr = (realdpcomplex *) tau.ObtainWriteEntireData();
        integer *jpvt = new integer[n];
        integer info;
        integer lwork = -1;
        realdpcomplex lworkopt;
        for (integer i = 0; i < n; i++)
            jpvt[i] = i + 1;
        realdp *rwork = new realdp[2 * n];
        /* compute the space required in the geqp3 */
        geqp3_(&m, &n, HHRptr, &m, jpvt, tauptr, &lworkopt, &lwork, rwork, &info);
        lwork = static_cast<integer> (lworkopt.r);
        realdpcomplex *work = new realdpcomplex[lwork];
        /* QR decomposition for ptrHHR using Householder reflections. Householder reflectors and R are stored in ptrHHR.
        details: www.netlib.org/lapack/explore-html/db/de5/geqp3_8f.html */
        geqp3_(&m, &n, HHRptr, &m, jpvt, tauptr, work, &lwork, rwork, &info);
        if (info < 0)
            printf("Error in qr decomposition!\n");

        for (integer i = 0; i < n; i++)
        {
            if (jpvt[i] != (i + 1))
                printf("Error in qr decomposition!\n");
        }
        AddToFields("_HHR", HHR);
        AddToFields("_tau", tau);
        delete[] jpvt;
        delete[] work;
        delete[] rwork;
    };

    void Element::QRDecom(void)
    {
        if(FieldsExist("_Q") && FieldsExist("_R"))
            return;
        
        HHRDecom();
        if(!iscomplex)
        {
            integer m = size[0], n = size[1], minmn = (m < n) ? m : n;
            Element HHR = Field("_HHR"), tau = Field("_tau");
            Element Q(m, minmn), R(minmn, n);
            const realdp *HHRptr = HHR.ObtainReadData();
            realdp *Rptr = R.ObtainWriteEntireData();
            realdp *Qptr = Q.ObtainWriteEntireData();
            for(integer i = 0; i < Q.Getlength(); i++)
                Qptr[i] = HHRptr[i];
            
            const realdp *tauptr = tau.ObtainReadData();
            for(integer i = 0; i < n; i++)
            {
                for(integer j = 0; j < minmn; j++)
                {
                    if(j <= i)
                        Rptr[j + i * minmn] = HHRptr[j + i * m];
                    else
                        Rptr[j + i * minmn] = 0;
                }
            }

            integer info;
            integer lwork = -1;
            realdp lworkopt;
            /* compute the space required in the orgqr */
            orgqr_(&m, &minmn, &minmn, Qptr, &m, const_cast<realdp *> (tauptr), &lworkopt, &lwork, &info);
            lwork = static_cast<integer> (lworkopt);
            realdp *work = new realdp[lwork];
            orgqr_(&m, &minmn, &minmn, Qptr, &m, const_cast<realdp *> (tauptr), work, &lwork, &info);
            if (info < 0)
                printf("Error in forming Q matrix!\n");
            
            AddToFields("_Q", Q);
            AddToFields("_R", R);
			delete[] work;
            return;
        }

        integer m = size[0] / 2, n = size[1], minmn = (m < n) ? m : n;
        Element HHR = Field("_HHR"), tau = Field("_tau");
        Element Q(m, minmn, 1, "complex"), R(minmn, n, 1, "complex");
        const realdp *HHRptr = HHR.ObtainReadData();
        realdp *Rptr = R.ObtainWriteEntireData();
        realdp *Qptr = Q.ObtainWritePartialData();
        for(integer i = 0; i < Q.Getlength(); i++)
            Qptr[i] = HHRptr[i];
        const realdpcomplex *tauptr = (realdpcomplex *) tau.ObtainReadData();
        
        for(integer i = 0; i < n; i++)
        {
            for(integer j = 0; j < minmn; j++)
            {
                if(j <= i)
                {
                    Rptr[2 * j + i * 2 * minmn] = HHRptr[2 * j + i * 2 * m];
                    Rptr[2 * j + 1 + i * 2 * minmn] = HHRptr[2 * j + 1 + i * 2 * m];
                }
                else
                {
                    Rptr[2 * j + 2 * i * minmn] = 0;
                    Rptr[2 * j + 1 + 2 * i * minmn] = 0;
                }
            }
        }

        integer info;
        integer lwork = -1;
        realdpcomplex lworkopt;
        /* compute the space required in the ungqr */
        ungqr_(&m, &minmn, &minmn, (realdpcomplex *) Qptr, &m, const_cast<realdpcomplex *> (tauptr), &lworkopt, &lwork, &info);
        lwork = static_cast<integer> (lworkopt.r);
        realdpcomplex *work = new realdpcomplex[lwork];
        ungqr_(&m, &minmn, &minmn, (realdpcomplex *) Qptr, &m, const_cast<realdpcomplex *> (tauptr), work, &lwork, &info);
        if (info < 0)
            printf("Error in forming Q matrix!\n");
        
        AddToFields("_Q", Q);
        AddToFields("_R", R);
		delete[] work;
        return;
    };

    Element Element::HHRMtp(Element HHR, Element tau, char *Trans, char *Side) const
    {
        assert(iscomplex == HHR.Getiscomplex() && iscomplex == tau.Getiscomplex());
        
        if(!iscomplex)
        {
            integer m = size[0], n = size[1], r = HHR.Getsize()[0], k = HHR.Getsize()[1], minrk = (r < k) ? r : k;
            Element result(*this);
            realdp *resultptr = result.ObtainWritePartialData();
            const realdp *HHRptr = HHR.ObtainReadData();
            const realdp *tauptr = tau.ObtainReadData();
            
            integer info;
            integer lwork = -1;
            realdp lworkopt;
            /* compute the size of space required in the ormqr */
            ormqr_(Side, Trans, &m, &n, &minrk, const_cast<realdp *> (HHRptr), &r, const_cast<realdp *> (tauptr), resultptr, &m, &lworkopt, &lwork, &info);
            lwork = static_cast<integer> (lworkopt);
            realdp *work = new realdp[lwork];
            
            ormqr_(Side, Trans, &m, &n, &minrk, const_cast<realdp *> (HHRptr), &r, const_cast<realdp *> (tauptr), resultptr, &m, work, &lwork, &info);
            delete [] work;
            
            return result;
        }
        
        integer m = size[0] / 2, n = size[1], r = HHR.Getsize()[0] / 2, k = HHR.Getsize()[1], minrk = (r < k) ? r : k;
        Element result(*this);
        realdpcomplex *resultptr = (realdpcomplex *) result.ObtainWritePartialData();
        const realdpcomplex *HHRptr = (realdpcomplex *) HHR.ObtainReadData();
        const realdpcomplex *tauptr = (realdpcomplex *) tau.ObtainReadData();
        
        integer info;
        integer lwork = -1;
        realdpcomplex lworkopt;

        unmqr_(Side, Trans, &m, &n, &minrk, const_cast<realdpcomplex *> (HHRptr), &r, const_cast<realdpcomplex *> (tauptr), resultptr, &m, &lworkopt, &lwork, &info);
        lwork = static_cast<integer> (lworkopt.r);
        realdpcomplex *work = new realdpcomplex[lwork];

        unmqr_(Side, Trans, &m, &n, &minrk, const_cast<realdpcomplex *> (HHRptr), &r, const_cast<realdpcomplex *> (tauptr), resultptr, &m, work, &lwork, &info);
        delete [] work;
        
        return result;
    };

    void Element::SVDDecom(void)
    {
        if(FieldsExist("_U") && FieldsExist("_S") && FieldsExist("_Vt"))
            return;
        
        if(!iscomplex)
        {
            integer m = size[0], n = size[1], minmn = (m < n) ? m : n;
            Element A(*this);
            Element S(minmn), U(m, minmn), Vt(minmn, n);
            realdp *Aptr = A.ObtainWritePartialData();
            realdp *Sptr = S.ObtainWriteEntireData();
            realdp *Uptr = U.ObtainWriteEntireData();
            realdp *Vtptr = Vt.ObtainWriteEntireData();

            realdp workoptsize;
            integer lwork = -1, info;
            /* compute the space required in the SVD computation. */
            gesvd_(GLOBAL::S, GLOBAL::S, &m, &n, Aptr, &m, Sptr, Uptr, &m, Vtptr, &minmn, &workoptsize, &lwork, &info);
            lwork = static_cast<integer> (workoptsize);
            realdp *work = new realdp[lwork];
            /* SVD: U * S * Vt = M, details: www.netlib.org/lapack/explore-html/d8/d2d/gesvd_8f.html */
            gesvd_(GLOBAL::S, GLOBAL::S, &m, &n, Aptr, &m, Sptr, Uptr, &m, Vtptr, &minmn, work, &lwork, &info);
            if (info != 0)
            {
                printf("Error:singular value decomposition failed!\n");
            }
            delete[] work;
            AddToFields("_U", U);
            AddToFields("_S", S);
            AddToFields("_Vt", Vt);
            return;
        }

        integer m = size[0] / 2, n = size[1], minmn = (m < n) ? m : n;
        Element A(*this);
        Element S(minmn, "complex"), U(m, minmn, "complex"), Vt(minmn, n, "complex");
        realdp *Aptr = A.ObtainWritePartialData();
        realdp *Sptr = S.ObtainWriteEntireData();
        realdp *Uptr = U.ObtainWriteEntireData();
        realdp *Vtptr = Vt.ObtainWriteEntireData();

        realdpcomplex workoptsize;
        integer lwork = -1, info;
        realdp *rwork = new realdp[5 * minmn];
        /* compute the space required in the SVD computation. */
        gesvd_(GLOBAL::S, GLOBAL::S, &m, &n, (realdpcomplex *) Aptr, &m, Sptr, (realdpcomplex *) Uptr, &m, (realdpcomplex *) Vtptr, &minmn, &workoptsize, &lwork, rwork, &info);
        lwork = static_cast<integer> (workoptsize.r);
        realdpcomplex *work = new realdpcomplex[lwork];
        /* SVD: U * S * Vt = M, details: www.netlib.org/lapack/explore-html/d8/d2d/gesvd_8f.html */
        gesvd_(GLOBAL::S, GLOBAL::S, &m, &n, (realdpcomplex *) Aptr, &m, Sptr, (realdpcomplex *) Uptr, &m, (realdpcomplex *) Vtptr, &minmn, work, &lwork, rwork, &info);
        if (info != 0)
        {
            printf("Error:singular value decomposition failed!\n");
        }
        delete[] work;
        delete[] rwork;
        
        for(integer i = minmn - 1; i >= 0; i--)
        {
            Sptr[2 * i + 1] = 0;
            Sptr[2 * i] = Sptr[i];
        }
        
        AddToFields("_U", U);
        AddToFields("_S", S);
        AddToFields("_Vt", Vt);
    };

    void Element::SchurForm(char *jobvs)
    {
        assert(Getrow() == Getcol());
        if(FieldsExist("_SchurForm") && (strcmp(jobvs, "N") == 0))
            return;
        if(FieldsExist("_SchurForm") && FieldsExist("_SchurVec"))
            return;
        integer n = Getrow(), sdim;
        if(!iscomplex)
        {
            Element SchurForm(*this);
            Element SchurVec(n, n);
            realdp *SchurFormptr = SchurForm.ObtainWritePartialData();
            realdp *eigsr = new realdp[2 * n];
            realdp *eigsi = eigsr + n;
            realdp lworkopt;
            integer lwork = -1;
            integer info;
            if(strcmp(jobvs, "N") == 0)
            {
                /* compute the size of space required in the gees */
                gees_(jobvs, GLOBAL::N, nullptr, &n, SchurFormptr, &n, &sdim, eigsr, eigsi, nullptr, &n, &lworkopt, &lwork, nullptr, &info);
                lwork = static_cast<integer> (lworkopt);
                realdp *work = new realdp[lwork];
                /* generalized schur decomposition for matrices A.
                details: www.netlib.org/lapack/explore-html/d8/d7e/gees_8f.html */
                gees_(jobvs, GLOBAL::N, nullptr, &n, SchurFormptr, &n, &sdim, eigsr, eigsi, nullptr, &n, work, &lwork, nullptr, &info);
                delete[] work;
            } else
            {
                realdp *SchurVecptr = SchurVec.ObtainWriteEntireData();
                /* compute the size of space required in the gees */
                gees_(jobvs, GLOBAL::N, nullptr, &n, SchurFormptr, &n, &sdim, eigsr, eigsi, SchurVecptr, &n, &lworkopt, &lwork, nullptr, &info);
                lwork = static_cast<integer> (lworkopt);
                realdp *work = new realdp[lwork];
                /* generalized schur decomposition for matrices A.
                details: www.netlib.org/lapack/explore-html/d8/d7e/gees_8f.html */
                gees_(jobvs, GLOBAL::N, nullptr, &n, SchurFormptr, &n, &sdim, eigsr, eigsi, SchurVecptr, &n, work, &lwork, nullptr, &info);
                delete[] work;
                AddToFields("_SchurVec", SchurVec);
            }
            delete[] eigsr;
            AddToFields("_SchurForm", SchurForm);
            if(info != 0)
                printf("Error: Schur decomposition failed!\n");
            return;
        }

        Element SchurForm(*this);
        Element SchurVec(n, n, "complex");
        realdpcomplex *SchurFormptr = (realdpcomplex *) SchurForm.ObtainWritePartialData();
        realdpcomplex *eigs = new realdpcomplex[n];
        realdp *rwork = new realdp[n];
        realdpcomplex lworkopt;
        integer lwork = -1;
        integer info;
        if(strcmp(jobvs, "N") == 0)
        {
            /* compute the size of space required in the gees */
            gees_(jobvs, GLOBAL::N, nullptr, &n, SchurFormptr, &n, &sdim, eigs, nullptr, &n, &lworkopt, &lwork, rwork, nullptr, &info);
            lwork = static_cast<integer> (lworkopt.r);
            realdpcomplex *work = new realdpcomplex[lwork];
            /* generalized schur decomposition for matrices A.
            details: www.netlib.org/lapack/explore-html/d8/d7e/gees_8f.html */
            gees_(jobvs, GLOBAL::N, nullptr, &n, SchurFormptr, &n, &sdim, eigs, nullptr, &n, work, &lwork, rwork, nullptr, &info);
            delete[] work;
        } else
        {
            realdpcomplex *SchurVecptr = (realdpcomplex *) SchurVec.ObtainWriteEntireData();
            /* compute the size of space required in the gees */
            gees_(jobvs, GLOBAL::N, nullptr, &n, SchurFormptr, &n, &sdim, eigs, SchurVecptr, &n, &lworkopt, &lwork, rwork, nullptr, &info);
            lwork = static_cast<integer> (lworkopt.r);
            realdpcomplex *work = new realdpcomplex[lwork];
            /* generalized schur decomposition for matrices A.
            details: www.netlib.org/lapack/explore-html/d8/d7e/gees_8f.html */
            gees_(jobvs, GLOBAL::N, nullptr, &n, SchurFormptr, &n, &sdim, eigs, SchurVecptr, &n, work, &lwork, rwork, nullptr, &info);
            delete[] work;
            AddToFields("_SchurVec", SchurVec);
        }
        delete[] eigs;
        delete[] rwork;
        AddToFields("_SchurForm", SchurForm);
        if(info != 0)
            printf("Error: Schur decomposition failed!\n");
        return;
    };

    Element Element::SYL(Element A, Element B)
    {
        assert(iscomplex == A.Getiscomplex() && iscomplex == B.Getiscomplex());
        assert(A.Getrow() == A.Getcol() && B.Getrow() == B.Getcol() && A.Getcol() == Getrow() && B.Getrow() == Getcol());
        integer n = A.Getrow(), m = B.Getrow(); /*A is n by n, B is m by m, C is n by m, X is n by m*/
        A.SchurForm();
        B.SchurForm();
        Element ASF = A.Field("_SchurForm"), ASV = A.Field("_SchurVec"), BSF = B.Field("_SchurForm"), BSV = B.Field("_SchurVec");
        Element C(*this);
        
        if(!iscomplex)
        {
            C = ASV.GetTranspose() * C * BSV;
            Element DIdentity(n, n), EIdentity(m, m), FZeros(n, m);
            DIdentity.SetToIdentity(); EIdentity.SetToIdentity(); FZeros.SetToZeros();
            realdp *DIdentityptr = DIdentity.ObtainWritePartialData();
            realdp *EIdentityptr = EIdentity.ObtainWritePartialData();
            realdp *FZerosptr = FZeros.ObtainWritePartialData();
            realdp *Cptr = C.ObtainWritePartialData();
            BSF = -1 * BSF;
            const realdp *ASFptr = ASF.ObtainReadData();
            const realdp *BSFptr = BSF.ObtainReadData();
            
            realdp scalar, dif;
            integer lwork = -1, info;
            integer *iwork = new integer[n + m + 6];
            realdp lworkopt;
            /* compute the size of space required in the tgsyl */
            tgsyl_(GLOBAL::N, &GLOBAL::IZERO, &n, &m, const_cast<realdp *> (ASFptr), &n, const_cast<realdp *> (BSFptr), &m, Cptr, &n,
                DIdentityptr, &n, EIdentityptr, &m, FZerosptr, &n, &scalar, &dif, &lworkopt, &lwork, iwork, &info);
            lwork = static_cast<integer> (lworkopt);
            realdp *work = new realdp[lwork];
            /* generalized Sylvester equation.
            details: www.netlib.org/lapack/explore-html/d8/d7e/gees_8f.html */
            tgsyl_(GLOBAL::N, &GLOBAL::IZERO, &n, &m, const_cast<realdp *> (ASFptr), &n, const_cast<realdp *> (BSFptr), &m, Cptr, &n,
                DIdentityptr, &n, EIdentityptr, &m, FZerosptr, &n, &scalar, &dif, work, &lwork, iwork, &info);
            
            delete[] work;
            delete[] iwork;
            if(info != 0)
                printf("warning: Matrix::DSYL may not be solved correctly!\n");
            
            return (ASV * C * BSV.GetTranspose()) / scalar;
        }

        C = ASV.GetTranspose() * C * BSV;
        Element DIdentity(n, n, "complex"), EIdentity(m, m, "complex"), FZeros(n, m, "complex");
        DIdentity.SetToIdentity(); EIdentity.SetToIdentity(); FZeros.SetToZeros();
        realdpcomplex *DIdentityptr = (realdpcomplex *) DIdentity.ObtainWritePartialData();
        realdpcomplex *EIdentityptr = (realdpcomplex *) EIdentity.ObtainWritePartialData();
        realdpcomplex *FZerosptr = (realdpcomplex *) FZeros.ObtainWritePartialData();
        realdpcomplex *Cptr = (realdpcomplex *) C.ObtainWritePartialData();
        realdpcomplex a={-1, 0};
        BSF = a * BSF;
        const realdpcomplex *ASFptr = (realdpcomplex *) ASF.ObtainReadData();
        const realdpcomplex *BSFptr = (realdpcomplex *) BSF.ObtainReadData();
        
        realdp scalar, dif;
        integer lwork = -1, info;
        integer *iwork = new integer[n + m + 2];
        realdpcomplex lworkopt;
        /* compute the size of space required in the tgsyl */
        tgsyl_(GLOBAL::N, &GLOBAL::IZERO, &n, &m, const_cast<realdpcomplex *> (ASFptr), &n, const_cast<realdpcomplex *> (BSFptr), &m, Cptr, &n,
            DIdentityptr, &n, EIdentityptr, &m, FZerosptr, &n, &scalar, &dif, &lworkopt, &lwork, iwork, &info);
        lwork = static_cast<integer> (lworkopt.r);
        realdpcomplex *work = new realdpcomplex[lwork];
        /* generalized Sylvester equation.
        details: www.netlib.org/lapack/explore-html/d8/d7e/gees_8f.html */
        tgsyl_(GLOBAL::N, &GLOBAL::IZERO, &n, &m, const_cast<realdpcomplex *> (ASFptr), &n, const_cast<realdpcomplex *> (BSFptr), &m, Cptr, &n,
            DIdentityptr, &n, EIdentityptr, &m, FZerosptr, &n, &scalar, &dif, work, &lwork, iwork, &info);
        
        delete[] work;
        delete[] iwork;
        if(info != 0)
            printf("warning: Matrix::DSYL may not be solved correctly!\n");
        a.r = scalar; a.i = 0;
        return (ASV * C * BSV.GetTranspose()) / a;
    };

    Element Element::GetDiagTimesM(Element M, char *side) const
    {
        assert(Space != nullptr);
        assert(M.Getiscomplex() == iscomplex);
        Element result(M);
        realdp *resultptr = result.ObtainWritePartialData();
        integer n = M.Getrow(), m = M.Getcol();
        
        if(strcmp(side, "L") == 0)
        {
            assert(length == M.Getsize()[0]);
            if(!iscomplex)
            {
                for(integer i = 0; i < n; i++)
                    scal_(&m, Space + i, resultptr + i, &n);
                return result;
            }
            realdp re, im;
            for(integer i = 0; i < n; i++)
            {
                for(integer j = 0; j < m; j++)
                {
                    re = resultptr[2 * i + j * 2 * n] * Space[2 * i] - resultptr[2 * i + 1 + j * 2 * n] * Space[2 * i + 1];
                    im = resultptr[2 * i + j * 2 * n] * Space[2 * i + 1] + resultptr[2 * i + 1 + j * 2 * n] * Space[2 * i];
                    resultptr[2 * i + j * 2 * n] = re;
                    resultptr[2 * i + 1 + j * 2 * n] = im;
                }
            }
            return result;
        }
        
        if(!iscomplex)
        {
            assert(length == M.Getsize()[1]);
            for(integer i = 0; i < m; i++)
                scal_(&n, Space + i, resultptr + i * n, &GLOBAL::IONE);
            return result;
        }
        assert(length == 2 * M.Getsize()[1]);
        realdp re, im;
        for(integer j = 0; j < m; j++)
        {
            for(integer i = 0; i < n; i++)
            {
                re = resultptr[2 * i + j * 2 * n] * Space[2 * j] - resultptr[2 * i + 1 + j * 2 * n] * Space[2 * j + 1];
                im = resultptr[2 * i + j * 2 * n] * Space[2 * j + 1] + resultptr[2 * i + 1 + j * 2 * n] * Space[2 * j];
                resultptr[2 * i + j * 2 * n] = re;
                resultptr[2 * i + 1 + j * 2 * n] = im;
            }
        }
        return result;
    };

    Element &Element::DiagTimesM(Element &M, char *side) const
    {
        assert(Space != nullptr);
        assert(M.Getiscomplex() == iscomplex);
        realdp *resultptr = M.ObtainWritePartialData();
        integer n = M.Getrow(), m = M.Getcol();
        
        if(strcmp(side, "L") == 0)
        {
            assert(length == M.Getsize()[0]);
            if(!iscomplex)
            {
                for(integer i = 0; i < n; i++)
                    scal_(&m, Space + i, resultptr + i, &n);
                return M;
            }
            realdp re, im;
            for(integer i = 0; i < n; i++)
            {
                for(integer j = 0; j < m; j++)
                {
                    re = resultptr[2 * i + j * 2 * n] * Space[2 * i] - resultptr[2 * i + 1 + j * 2 * n] * Space[2 * i + 1];
                    im = resultptr[2 * i + j * 2 * n] * Space[2 * i + 1] + resultptr[2 * i + 1 + j * 2 * n] * Space[2 * i];
                    resultptr[2 * i + j * 2 * n] = re;
                    resultptr[2 * i + 1 + j * 2 * n] = im;
                }
            }
            return M;
        }
        
        if(!iscomplex)
        {
            assert(length == M.Getsize()[1]);
            for(integer i = 0; i < m; i++)
                scal_(&n, Space + i, resultptr + i * n, &GLOBAL::IONE);
            return M;
        }
        assert(length == 2 * M.Getsize()[1]);
        realdp re, im;
        for(integer j = 0; j < m; j++)
        {
            for(integer i = 0; i < n; i++)
            {
                re = resultptr[2 * i + j * 2 * n] * Space[2 * j] - resultptr[2 * i + 1 + j * 2 * n] * Space[2 * j + 1];
                im = resultptr[2 * i + j * 2 * n] * Space[2 * j + 1] + resultptr[2 * i + 1 + j * 2 * n] * Space[2 * j];
                resultptr[2 * i + j * 2 * n] = re;
                resultptr[2 * i + 1 + j * 2 * n] = im;
            }
        }
        return M;
    };

    Element &Element::HaddRankone(realdp scalar, const Element &u, const Element &v)
    {
        assert(u.Getcol() == 1 && u.Getnum() == 1 && u.Getrow() == Getrow());
        assert(v.Getcol() == 1 && v.Getnum() == 1 && v.Getrow() == Getcol());
        assert(Space != nullptr);
        integer m = Getrow(), n = Getcol();
        const realdp *uptr = u.ObtainReadData();
        const realdp *vptr = v.ObtainReadData();

        ger_(&m, &n, &scalar, const_cast<realdp *> (uptr), &GLOBAL::IONE, const_cast<realdp *> (vptr), &GLOBAL::IONE, Space, &m);
        
        return *this;
    };

    Element &Element::HaddRankone(realdpcomplex scalar, const Element &u, const Element &v)
    {
        assert(u.Getcol() == 1 && u.Getnum() == 1 && u.Getrow() == Getrow());
        assert(v.Getcol() == 1 && v.Getnum() == 1 && v.Getrow() == Getcol());
        assert(Space != nullptr);
        assert(iscomplex == true && u.Getiscomplex() == true && v.Getiscomplex() == true);
        integer m = Getrow(), n = Getcol();
        const realdp *uptr = u.ObtainReadData();
        const realdp *vptr = v.ObtainReadData();

        ger_(&m, &n, &scalar, const_cast<realdpcomplex *> ((realdpcomplex *) uptr), &GLOBAL::IONE, const_cast<realdpcomplex *> ((realdpcomplex *) vptr), &GLOBAL::IONE,
                 (realdpcomplex *) Space, &m);
        
        return *this;
    };

    Element &Element::AlphaXaddThis(realdp alpha, const Element &X)
    {
        assert(length == X.Getlength());
        const realdp *Xptr = X.ObtainReadData();
        realdp *thisptr = this->ObtainWritePartialData();
        axpy_(&length, &alpha, const_cast<realdp *> (Xptr), &GLOBAL::IONE, thisptr, &GLOBAL::IONE);
        return (*this);
    };

    Element &Element::AlphaXaddThis(realdpcomplex alpha, const Element &X)
    {
        assert(length == X.Getlength() && X.Getiscomplex());
        const realdp *Xptr = X.ObtainReadData();
        realdp *thisptr = this->ObtainWritePartialData();
        integer L = length / 2;
        axpy_(&L, &alpha, const_cast<realdpcomplex *> ((realdpcomplex *) Xptr), &GLOBAL::IONE, (realdpcomplex *) thisptr, &GLOBAL::IONE);
        return (*this);
    };

    Element &Element::ScalarTimesThis(realdp alpha)
    {
        realdp *thisptr = this->ObtainWritePartialData();
        scal_(&length, &alpha, thisptr, &GLOBAL::IONE);
        return (*this);
    };

    Element &Element::ScalarTimesThis(realdpcomplex alpha)
    {
        realdp *thisptr = this->ObtainWritePartialData();
        integer L = length / 2;
        scal_(&L, &alpha, (realdpcomplex *) thisptr, &GLOBAL::IONE);
        return (*this);
    };

    Element &Element::AlphaABaddBetaThis(realdp alpha, const Element &A, char *transA, const Element &B, char*transB, realdp beta)
    {
        if(!A.Getiscomplex() && !B.Getiscomplex())
        { /* if both A and B are real*/
            integer rowA = A.Getrow(), colA = A.Getcol(), rowB = B.Getrow(), colB = B.Getcol(), rowOutput = Getrow(), colOutput = Getcol(), row = 0, col = 0, k = 0;
            const realdp *Aptr = A.ObtainReadData();
            const realdp *Bptr = B.ObtainReadData();
            realdp *thisptr = this->ObtainWritePartialData();
            if(strcmp(transA, "N") == 0) /* transA == "N" */
            {
                k = colA;
                row = rowA;
            }
            else
            {
                row = colA;
                k = rowA;
            }
            if(strcmp(transB, "N") == 0) /* transA == "N" */
            {
                assert(k == rowB);
                col = colB;
            }
            else
            {
                assert(k == colB);
                col = rowB;
            }
            assert(rowOutput * colOutput == row * col);

            gemm_(transA, transB, &row, &col, &k, &alpha, const_cast<realdp *> (Aptr), &rowA, const_cast<realdp *> (Bptr), &rowB, &beta, thisptr, &row);
            return (*this);
        }
        
        if (A.Getiscomplex() && A.Getiscomplex())
        { /* if both A and B are complex*/
            realdpcomplex calpha = {alpha, 0}, cbeta = {beta, 0};
            AlphaABaddBetaThis(calpha, A, transA, B, transB, cbeta);
            return *this;
        }
        assert(false); /* types of A and B are not supported. */
        return *this;
    };

    Element &Element::AlphaABaddBetaThis(realdpcomplex alpha, const Element &A, char *transA, const Element &B, char*transB, realdpcomplex beta)
    {
        assert(A.Getiscomplex() && B.Getiscomplex());
        integer rowA = A.Getrow(), colA = A.Getcol(), rowB = B.Getrow(), colB = B.Getcol(), rowOutput = Getrow(), colOutput = Getcol(), row = 0, col = 0, k = 0;
        const realdpcomplex *Aptr = (realdpcomplex *) A.ObtainReadData();
        const realdpcomplex *Bptr = (realdpcomplex *) B.ObtainReadData();
        realdp *thisptr = this->ObtainWritePartialData();
        if(strcmp(transA, "N") == 0) /* transA == "N" */
        {
            k = colA;
            row = rowA;
        }
        else
        {
            row = colA;
            k = rowA;
        }
        if(strcmp(transB, "N") == 0) /* transA == "N" */
        {
            assert(k == rowB);
            col = colB;
        }
        else
        {
            assert(k == colB);
            col = rowB;
        }
        assert(rowOutput * colOutput == row * col);
        gemm_(transA, transB, &row, &col, &k, &alpha, const_cast<realdpcomplex *> (Aptr), &rowA, const_cast<realdpcomplex *> (Bptr), &rowB, &beta, (realdpcomplex *) thisptr, &row);
        return (*this);
    };

    realdp Element::DotProduct(const Element &M) const
    {
        assert(length == M.Getlength());
        if(length == 0)
        {
            return 0;
        }
        assert(Space != nullptr);
        const realdp *Mptr = M.ObtainReadData();
        integer llength = length;
        return dot_(&llength, const_cast<realdp *> (Space), &GLOBAL::IONE, const_cast<realdp *> (Mptr), &GLOBAL::IONE);
    };

    Element Element::GetHadamardProduct(Element M) const
    {
        assert(length == M.Getlength() && iscomplex == M.Getiscomplex() && Space != nullptr);
        Element result(M);
        realdp *resultptr = result.ObtainWritePartialData();
        if(!iscomplex)
        {
            for(integer i = 0; i < length; i++)
                resultptr[i] *= Space[i];
            return result;
        }
        realdp re, im;
        for(integer i = 0; i < length; i++, i++)
        {
            re = resultptr[i] * Space[i] - resultptr[i + 1] * Space[i + 1];
            im = resultptr[i + 1] * Space[i] + resultptr[i] * Space[i + 1];
            resultptr[i] = re;
            resultptr[i + 1] = im;
        }
        return result;
    };

    Element Element::GetHadamardDivision(Element M) const
    {
        assert(length == M.Getlength() && iscomplex == M.Getiscomplex());
        Element result(*this);
        realdp *resultptr = result.ObtainWritePartialData();
        const realdp *Mptr = M.ObtainReadData();
        if(!iscomplex)
        {
            for(integer i = 0; i < length; i++)
                resultptr[i] /= Mptr[i];
            return result;
        }
        realdp re, im, de;
        for(integer i = 0; i < length; i++, i++)
        {
            de = Mptr[i] * Mptr[i] + Mptr[i + 1] * Mptr[i + 1];
            re = (resultptr[i] * Mptr[i] + resultptr[i + 1] * Mptr[i + 1]) / de;
            im = (resultptr[i + 1] * Mptr[i] - resultptr[i] * Mptr[i + 1]) / de;
            resultptr[i] = re;
            resultptr[i + 1] = im;
        }
        return result;
    };

    Element &Element::Reshape(integer r, integer c, integer n)
    {
        if(!iscomplex)
        {
            assert(r * c * n == length);
            size[0] = r;
            size[1] = c;
            size[2] = n;
            return (*this);
        }
        assert(2 * r * c * n == length);
        size[0] = r * 2;
        size[1] = c;
        size[2] = n;
        return (*this);
    };

    Element Element::GetReshape(integer r, integer c, integer n) const
    {
        Element result(*this);
        return result.Reshape(r, c, n);
    };

    Element Element::GetHaddRankone(realdp scalar, Element u, Element v) const
    {
        assert(u.Getcol() == 1 && u.Getnum() == 1 && u.Getrow() == Getrow());
        assert(v.Getcol() == 1 && v.Getnum() == 1 && v.Getrow() == Getcol());
        integer m = Getrow(), n = Getcol();
        Element result(*this);
        realdp *resultptr = result.ObtainWritePartialData();
        const realdp *uptr = u.ObtainReadData();
        const realdp *vptr = v.ObtainReadData();

        ger_(&m, &n, &scalar, const_cast<realdp *> (uptr), &GLOBAL::IONE, const_cast<realdp *> (vptr), &GLOBAL::IONE, resultptr, &m);
        
        return result;
    };

    Element Element::GetHaddRankone(realdpcomplex scalar, Element u, Element v) const
    {
        assert(u.Getcol() == 1 && u.Getnum() == 1 && u.Getrow() == Getrow());
        assert(v.Getcol() == 1 && v.Getnum() == 1 && v.Getrow() == Getcol());
        assert(iscomplex == true && u.Getiscomplex() == true && v.Getiscomplex() == true);
        integer m = Getrow(), n = Getcol();
        Element result(*this);
        realdp *resultptr = result.ObtainWritePartialData();
        const realdp *uptr = u.ObtainReadData();
        const realdp *vptr = v.ObtainReadData();

        ger_(&m, &n, &scalar, const_cast<realdpcomplex *> ((realdpcomplex *) uptr), &GLOBAL::IONE, const_cast<realdpcomplex *> ((realdpcomplex *) vptr), &GLOBAL::IONE,
                 (realdpcomplex *) resultptr, &m);
        
        return result;
    };

    Element Element::GetMax(Element u) const
    {
        assert(u.Getlength() == this->Getlength() && !iscomplex && ! u.Getiscomplex() && Space != nullptr);
        Element result(*this);
        realdp *resultptr = result.ObtainWriteEntireData();
        const realdp *uptr = u.ObtainReadData();
        for(integer i = 0; i < length; i++)
            resultptr[i] = (uptr[i] < Space[i]) ? Space[i] : uptr[i];
        return result;
    };

    Element Element::GetMin(Element u) const
    {
        assert(u.Getlength() == this->Getlength() && !iscomplex && ! u.Getiscomplex() && Space != nullptr);
        Element result(*this);
        realdp *resultptr = result.ObtainWriteEntireData();
        const realdp *uptr = u.ObtainReadData();
        for(integer i = 0; i < length; i++)
            resultptr[i] = (uptr[i] > Space[i]) ? Space[i] : uptr[i];
        return result;
    };

    Element Element::GetMax(realdp value) const
    {
        assert(Space != nullptr);
        Element result(*this);
        realdp *resultptr = result.ObtainWriteEntireData();
        for(integer i = 0; i < length; i++)
            resultptr[i] = (value < Space[i]) ? Space[i] : value;
        return result;
    };

    Element Element::GetMin(realdp value) const
    {
        assert(Space != nullptr);
        Element result(*this);
        realdp *resultptr = result.ObtainWriteEntireData();
        for(integer i = 0; i < length; i++)
            resultptr[i] = (value > Space[i]) ? Space[i] : value;
        return result;
    };

    Element Element::GetAbs(void) const
    {
        assert(Space != nullptr);
        Element result(*this);
        realdp *resultptr = result.ObtainWriteEntireData();
        for(integer i = 0; i < length; i++)
            resultptr[i] = std::fabs(Space[i]);
        return result;
    };

    Element Element::GetSqrt(void) const
    {
        assert(Space != nullptr);
        Element result(*this);
        realdp *resultptr = result.ObtainWriteEntireData();
        for(integer i = 0; i < length; i++)
            resultptr[i] = std::sqrt(Space[i]);
        return result;
    };

    Element Element::GetOrth(void) const
    {
        assert(Space != nullptr);
        Vector tmp(*this);
        tmp.SVDDecom();
        Element result = tmp.Field("_U") * tmp.Field("_Vt");
        return result;
    };

    Element Element::GetConj(void) const
    {
        assert(iscomplex && Space != nullptr);
        integer num = length / 2;
        Vector result(*this);
        realdp *resultptr = result.ObtainWritePartialData();
        scal_(&num, &GLOBAL::DNONE, resultptr + 1, &GLOBAL::ITWO);
        return result;
    };

    Element Element::Conj(void)
    {
        assert(Space != nullptr);
        assert(iscomplex);
        integer num = length / 2;
        scal_(&num, &GLOBAL::DNONE, Space + 1, &GLOBAL::ITWO);
        return (*this);
    };

#ifdef ROPTLIB_WITH_FFTW
    Element Element::GetFFT2D(int direction) const
    {
        assert(iscomplex && Space != nullptr);
        integer n1 = Getrow(), n2 = Getcol();
        unsigned int flags = FFTW_ESTIMATE;
        
        Vector result(*this);
        realdp *resultptr = result.ObtainWriteEntireData();
        /*Note that ROPTLIB uses column-major and FFTW uses row-major, therefore the
        row and column need be swapped.*/
        
#ifdef SINGLE_PRECISION
        fftwf_plan p = fftwf_plan_dft_2d(n2, n1, (fftwf_complex *) Space, (fftwf_complex *) resultptr, direction, flags);
        fftwf_execute(p);
        fftwf_destroy_plan(p);
#else
        fftw_plan p = fftw_plan_dft_2d(n2, n1, (fftw_complex *) Space, (fftw_complex *) resultptr, direction, flags);
        fftw_execute(p);
        fftw_destroy_plan(p);
#endif
        
        return result;
    };

    Element Element::FFT2D(int direction, Element *result) const
    {
        assert(iscomplex && Space != nullptr);
        integer n1 = Getrow(), n2 = Getcol();
        unsigned int flags = FFTW_ESTIMATE;
                
        realdp *resultptr = result->ObtainWriteEntireData();
        /*Note that ROPTLIB uses column-major and FFTW uses row-major, therefore the
        row and column need be swapped.*/
                
#ifdef SINGLE_PRECISION
        fftwf_plan p = fftwf_plan_dft_2d(n2, n1, (fftwf_complex *) Space, (fftwf_complex *) resultptr, direction, flags);
        fftwf_execute(p);
        fftwf_destroy_plan(p);
#else
        fftw_plan p = fftw_plan_dft_2d(n2, n1, (fftw_complex *) Space, (fftw_complex *) resultptr, direction, flags);
        fftw_execute(p);
        fftw_destroy_plan(p);
#endif
                
        return *result;
    };
#endif

    Element Element::GetHaarFWT(void) const
    {
        assert(iscomplex);
        integer n1 = Getrow(), n2 = Getcol();
        Vector result(*this);
        
        realdpcomplex *vv = (realdpcomplex *) result.ObtainWritePartialData();

        realdp r2 = static_cast<realdp> (sqrt(2.0));

        realdpcomplex *tmp = new realdpcomplex[n1 * n2];

        for (integer i = 0; i < n1 * n2; i++)
        {
            tmp[i].r = vv[i].r;
            tmp[i].i = vv[i].i;
        }

        integer k = 1;
        while (k * 2 <= n1)
        {
            k = k * 2;
        }
        while (1 < k)
        {
            k = k / 2;

            for (integer j = 0; j < n2; j++)
            {
                for (integer i = 0; i < k; i++)
                {
                    tmp[i + j * n1].r = (vv[2 * i + j * n1].r + vv[2 * i + 1 + j * n1].r) / r2;
                    tmp[i + j * n1].i = (vv[2 * i + j * n1].i + vv[2 * i + 1 + j * n1].i) / r2;
                    tmp[k + i + j * n1].r = (vv[2 * i + j*n1].r - vv[2 * i + 1 + j * n1].r) / r2;
                    tmp[k + i + j * n1].i = (vv[2 * i + j*n1].i - vv[2 * i + 1 + j * n1].i) / r2;
                }
            }
            for (integer j = 0; j < n2; j++)
            {
                for (integer i = 0; i < 2 * k; i++)
                {
                    vv[i + j * n1].r = tmp[i + j * n1].r;
                    vv[i + j * n1].i = tmp[i + j * n1].i;
                }
            }
        }
        k = 1;
        while (k * 2 <= n2)
        {
            k = k * 2;
        }
        while (1 < k)
        {
            k = k / 2;

            for (integer j = 0; j < k; j++)
            {
                for (integer i = 0; i < n1; i++)
                {
                    tmp[i + j * n1].r = (vv[i + 2 * j * n1].r + vv[i + (2 * j + 1) * n1].r) / r2;
                    tmp[i + j * n1].i = (vv[i + 2 * j * n1].i + vv[i + (2 * j + 1) * n1].i) / r2;
                    tmp[i + (k + j) * n1].r = (vv[i + 2 * j*n1].r - vv[i + (2 * j + 1) * n1].r) / r2;
                    tmp[i + (k + j) * n1].i = (vv[i + 2 * j*n1].i - vv[i + (2 * j + 1) * n1].i) / r2;
                }
            }

            for (integer j = 0; j < 2 * k; j++)
            {
                for (integer i = 0; i < n1; i++)
                {
                    vv[i + j*n1].r = tmp[i + j*n1].r;
                    vv[i + j*n1].i = tmp[i + j*n1].i;
                }
            }
        }
        delete[] tmp;

        return result;
    };

    Element Element::GetInvHaarFWT(void)
    {
        assert(iscomplex);
        integer n1 = Getrow(), n2 = Getcol();
        Vector result(*this);
        
        realdpcomplex *vv = (realdpcomplex *) result.ObtainWritePartialData();

        realdp r2 = static_cast<realdp> (sqrt(2.0));

        realdpcomplex *tmp = new realdpcomplex[n1 * n2];

        for (integer j = 0; j < n2; j++)
        {
            for (integer i = 0; i < n1; i++)
            {
                tmp[i + j * n1] = vv[i + j * n1];
            }
        }
        integer k = 1;

        while (k * 2 <= n2)
        {
            for (integer j = 0; j < k; j++)
            {
                for (integer i = 0; i < n1; i++)
                {
                    tmp[i + (2 * j) * n1].r = (vv[i + j * n1].r + vv[i + (k + j)*n1].r) / r2;
                    tmp[i + (2 * j) * n1].i = (vv[i + j * n1].i + vv[i + (k + j)*n1].i) / r2;
                    tmp[i + (2 * j + 1) * n1].r = (vv[i + j * n1].r - vv[i + (k + j)*n1].r) / r2;
                    tmp[i + (2 * j + 1) * n1].i = (vv[i + j * n1].i - vv[i + (k + j)*n1].i) / r2;
                }
            }

            for (integer j = 0; j < 2 * k; j++)
            {
                for (integer i = 0; i < n1; i++)
                {
                    vv[i + j * n1].r = tmp[i + j * n1].r;
                    vv[i + j * n1].i = tmp[i + j * n1].i;
                }
            }
            k = k * 2;
        }

        k = 1;
        while (k * 2 <= n1)
        {
            for (integer j = 0; j < n2; j++)
            {
                for (integer i = 0; i < k; i++)
                {
                    tmp[2 * i + j * n1].r = (vv[i + j * n1].r + vv[k + i + j * n1].r) / r2;
                    tmp[2 * i + j * n1].i = (vv[i + j * n1].i + vv[k + i + j * n1].i) / r2;
                    tmp[2 * i + 1 + j * n1].r = (vv[i + j * n1].r - vv[k + i + j * n1].r) / r2;
                    tmp[2 * i + 1 + j * n1].i = (vv[i + j * n1].i - vv[k + i + j * n1].i) / r2;
                }
            }

            for (integer j = 0; j < n2; j++)
            {
                for (integer i = 0; i < 2 * k; i++)
                {
                    vv[i + j * n1].r = tmp[i + j * n1].r;
                    vv[i + j * n1].i = tmp[i + j * n1].i;
                }
            }
            k = k * 2;
        }
        delete[] tmp;
        
        return result;
    };

    Element Element::GetHtimesRankone(realdp scalar, Element v, char *side) const
    {
        assert(!iscomplex && !(v.Getiscomplex()));
        Element result(*this);
        realdp *resultptr = result.ObtainWritePartialData();
        const realdp *vptr = v.ObtainReadData();
        integer m = Getrow(), n = Getcol();
        if(strcmp(side, "L") == 0) /*side == "L"*/
        {
            assert(Getrow() == v.Getrow());
            realdp *work = new realdp[n];
            larfx_(side, &m, &n, const_cast<realdp *> (vptr), &scalar, resultptr, &m, work);
            delete[] work;
            return result;
        }
        assert(Getcol() == v.Getrow());
        realdp *work = new realdp[m];
        larfx_(side, &m, &n, const_cast<realdp *> (vptr), &scalar, resultptr, &m, work);
        delete[] work;
        return result;
    };

    Element Element::GetHtimesRankone(realdpcomplex scalar, Element v, char *side) const
    {
        assert(iscomplex && v.Getiscomplex());
        Element result(*this);
        realdpcomplex *resultptr = (realdpcomplex *) result.ObtainWritePartialData();
        const realdpcomplex *vptr = (realdpcomplex *) v.ObtainReadData();
        integer m = Getrow(), n = Getcol();
        if(strcmp(side, "L") == 0) /*side == "L"*/
        {
            assert(Getrow() == v.Getrow());
            realdpcomplex *work = new realdpcomplex[n];
            larfx_(side, &m, &n, const_cast<realdpcomplex *> (vptr), &scalar, resultptr, &m, work);
            delete[] work;
            return result;
        }
        assert(Getcol() == v.Getrow());
        realdpcomplex *work = new realdpcomplex[m];
        larfx_(side, &m, &n, const_cast<realdpcomplex *> (vptr), &scalar, resultptr, &m, work);
        delete[] work;
        return result;
    };

    Element Element::GetSubmatrix(integer rstart, integer rend, integer cstart, integer cend) const
    {
        assert(rstart >= 0 && rend <= Getrow() - 1 && cstart >= 0 && cend <= Getcol() - 1 && Space != nullptr);
        if(!iscomplex)
        {
            Element result(rend - rstart + 1, cend - cstart + 1);
            realdp *resultptr = result.ObtainWriteEntireData();
            for(integer j = cstart; j <= cend; j++)
            {
                for(integer i = rstart; i <= rend; i++)
                {
                    resultptr[(i - rstart) + (j - cstart) * (rend - rstart + 1)] = Space[i + j * size[0]];
                }
            }
            return result;
        }
        
        Element result(rend - rstart + 1, cend - cstart + 1, "complex");
        realdp *resultptr = result.ObtainWriteEntireData();
        for(integer j = cstart; j <= cend; j++)
        {
            for(integer i = rstart; i <= rend; i++)
            {
                resultptr[2 * (i - rstart) + (j - cstart) * (rend - rstart + 1) * 2] = Space[2 * i + j * size[0]];
                resultptr[2 * (i - rstart) + 1 + (j - cstart) * (rend - rstart + 1) * 2] = Space[2 * i + 1 + j * size[0]];
            }
        }
        return result;
    };

    Element Element::SubmatrixAssignment(integer rstart, integer rend, integer cstart, integer cend, const Element &mat)
    {
        assert(rstart >= 0 && rend <= Getrow() - 1 && cstart >= 0 && cend <= Getcol() - 1);
        assert(rend - rstart + 1 == mat.Getrow() && cend - cstart + 1 == mat.Getcol() && mat.Getiscomplex() == iscomplex);
        const realdp *matptr = mat.ObtainReadData();
        integer matrow = mat.Getrow();
        this->NewMemoryOnWrite();
        if(!iscomplex)
        {
            for(integer j = cstart; j <= cend; j++)
            {
                for(integer i = rstart; i <= rend; i++)
                {
                    Space[i + j * size[0]] = matptr[i - rstart + (j - cstart) * matrow];
                }
            }
            return *this;
        }
        
        for(integer j = cstart; j <= cend; j++)
        {
            for(integer i = rstart; i <= rend; i++)
            {
                Space[2 * i + j * size[0]] = matptr[2 * (i - rstart) + (j - cstart) * matrow * 2];
                Space[2 * i + 1 + j * size[0]] = matptr[2 * (i - rstart) + 1 + (j - cstart) * matrow * 2];
            }
        }
        return *this;
    };

    Element Element::GetTranspose() const
    {
        assert(Space != nullptr);
        if(!iscomplex)
        {
            Element result(size[1], size[0]);
            realdp *resultptr = result.ObtainWriteEntireData();
            for(integer i = 0; i < size[0]; i++)
                for(integer j = 0; j < size[1]; j++)
                    resultptr[j + i * size[1]] = Space[i + j * size[0]];
            return result;
        }
        
        integer m = size[0] / 2, n = size[1];
        Element result(n, m, 1, "complex");
        realdp *resultptr = result.ObtainWriteEntireData();
        
        for(integer i = 0; i < n; i++)
        {
            for(integer j = 0; j < m; j++)
            {
                resultptr[2 * i + 2 * j * n] = Space[2 * j + i * 2 * m];
                resultptr[2 * i + 1 + 2 * j * n] = - Space[2 * j + 1 + i * 2 * m];
            }
        }
        
        return result;
    };

    Element Element::Transpose()
    {
        GetTranspose().CopyTo(*this);
        return (*this);
    };

    Element::Element(void)
    {
        ls = 3;
        size = new integer [3];
        size[0] = 0; size[1] = 0; size[2] = 0;
        length = 0;
        Space = nullptr;
        sharedtimes = nullptr;
        iscomplex = false;
        
        numoftypes = 0;
        powsinterval = nullptr;
        numofelements = 0;
        elements = nullptr;
    };

    Element::Element(integer r, integer l, integer n, const char *type)
    {
        if(strcmp(type, "real") == 0) /*type == 'real'*/
            iscomplex = false;
        else
            iscomplex = true;
        
        if(iscomplex)
            Initialization(3, r * 2, l, n);
        else
            Initialization(3, r, l, n);
        
        numoftypes = 0;
        powsinterval = nullptr;
        numofelements = 0;
        elements = nullptr;
    };

    Element::Element(integer r, integer l, integer n, bool iniscomplex)
    {
        iscomplex = iniscomplex;
        
        if(iscomplex)
            Initialization(3, r * 2, l, n);
        else
            Initialization(3, r, l, n);
        
        numoftypes = 0;
        powsinterval = nullptr;
        numofelements = 0;
        elements = nullptr;
    };

    Element::Element(integer r, integer l, const char *type)
    {
        if(strcmp(type, "real") == 0) /*type == 'real'*/
            iscomplex = false;
        else
            iscomplex = true;
        
        if(iscomplex)
            Initialization(3, r * 2, l, 1);
        else
            Initialization(3, r, l, 1);
        
        numoftypes = 0;
        powsinterval = nullptr;
        numofelements = 0;
        elements = nullptr;
    };

    Element::Element(integer r, const char *type)
    {
        if(strcmp(type, "real") == 0) /*type == 'real'*/
            iscomplex = false;
        else
            iscomplex = true;
        
        if(iscomplex)
            Initialization(3, r * 2, 1, 1);
        else
            Initialization(3, r, 1, 1);
        
        numoftypes = 0;
        powsinterval = nullptr;
        numofelements = 0;
        elements = nullptr;
    };

    Element::Element(integer numberoftypes, Element *ElementsTypes, integer *inpowsinterval)
    {
        iscomplex = false;
        numoftypes = numberoftypes;
        powsinterval = new integer[numoftypes + 1];
        for(integer i = 0; i <= numoftypes; i++)
            powsinterval[i] = inpowsinterval[i];
        
        Element **types = new Element *[numoftypes];
        length = 0;
        for(integer i = 0; i < numoftypes; i++)
        {
            types[i] = &ElementsTypes[i];
            length += types[i]->Getlength() * (powsinterval[i + 1] - powsinterval[i]);
        }

        numofelements = powsinterval[numoftypes];
        elements = new Element [numofelements];

        ls = 3;
        size = new integer[ls];
        size[0] = length; size[1] = 1; size[2] = 1;
        Space = nullptr;
        sharedtimes = nullptr;
        integer *isize = nullptr;

        for (integer i = 0; i < numoftypes; i++)
        {
            isize = new integer[types[i]->Getls()];
            for (integer j = 0; j < types[i]->Getls(); j++)
            {
                isize[j] = types[i]->Getsize()[j];
            }

            for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
            {
                if(types[i]->GetSpace() != nullptr)
                {
                    printf("Warning: the types of element for creating a collection of multiple elements are not empty. This will cause memory leakage!\n");
                }
                elements[j] = *(types[i]);
                delete[] elements[j].Getsize();
                elements[j].SetByParams(isize, types[i]->Getls(), types[i]->Getlength(), nullptr, nullptr);
                
            }
        }

        delete[] types;
    };

    Element::Element(integer numberoftypes, Element *FirstVar, integer Firstnum, ...)
    {
        iscomplex = false;
        numoftypes = numberoftypes;
        powsinterval = new integer[numoftypes + 1];

        Element **types = new Element *[numoftypes];

        va_list argptr;
        va_start(argptr, Firstnum);
        powsinterval[0] = 0;
        types[0] = FirstVar;
        powsinterval[1] = Firstnum;
        length = types[0]->Getlength() *Firstnum;
        for (integer i = 1; i < numoftypes; i++)
        {
            types[i] = va_arg(argptr, Element *);
            powsinterval[i + 1] = powsinterval[i] + va_arg(argptr, integer);
            length += types[i]->Getlength() * (powsinterval[i + 1] - powsinterval[i]);
        }
        va_end(argptr);

        numofelements = powsinterval[numoftypes];
        elements = new Element [numofelements];

        ls = 3;
        size = new integer[ls];
        size[0] = length; size[1] = 1; size[2] = 1;
        Space = nullptr;
        sharedtimes = nullptr;
        integer *isize = nullptr;

        for (integer i = 0; i < numoftypes; i++)
        {
            isize = new integer[types[i]->Getls()];
            for (integer j = 0; j < types[i]->Getls(); j++)
            {
                isize[j] = types[i]->Getsize()[j];
            }

            for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
            {
                if(types[i]->GetSpace() != nullptr)
                {
                    printf("Warning: the types of element for creating a collection of multiple elements are not empty. This will cause memory leakage!\n");
                }
                elements[j] = *(types[i]);
                delete[] elements[j].Getsize();
                elements[j].SetByParams(isize, types[i]->Getls(), types[i]->Getlength(), nullptr, nullptr);
                
            }
        }

        delete[] types;
    };

    realdp Element::Fnorm(void) const
    {
        assert(Space != nullptr);
        integer llength = length;
        realdp result = dot_(&llength, Space, &GLOBAL::IONE, Space, &GLOBAL::IONE);
        return std::sqrt(result);
    };

    Element Element::GetColNormsSquare(void) const
    {
        assert(Space != nullptr);
        Element result(1, Getcol());
        integer llength = size[0];
        
        realdp *resultptr = result.ObtainWriteEntireData();
        
        for(integer i = 0; i < Getcol(); i++)
        {
            resultptr[i] = dot_(&llength, Space + i * llength, &GLOBAL::IONE, Space + i * llength, &GLOBAL::IONE);
        }
        
        return result;
    };

    Element Element::GetColDotProducts(Element M) const
    {
        assert(Getrow() == M.Getrow() && Getcol() == M.Getcol() && M.Getiscomplex() == iscomplex && Space != nullptr);
        Element result(1, Getcol());
        integer llength = size[0];
        
        realdp *resultptr = result.ObtainWriteEntireData();
        const realdp *Mptr = M.ObtainReadData();
        
        for(integer i = 0; i < Getcol(); i++)
        {
            resultptr[i] = dot_(&llength, Space + i * llength, &GLOBAL::IONE, const_cast<realdp *> (Mptr + i * llength), &GLOBAL::IONE);
        }
        
        return result;
    };

    Element Element::GetRealToComplex(void) const
    {
        assert(!iscomplex && Space != nullptr);
        Element result(Getrow(), Getcol(), Getnum(), "complex");
        result.SetToZeros();
        realdp *resultptr = result.ObtainWritePartialData();
        
        for(integer i = 0; i < length; i++)
        {
            resultptr[2 * i] = Space[i];
        }
        return result;
    };

    Element Element::GetRealInComplex(void) const
    {
        assert(iscomplex && Space != nullptr);
        Element result(Getrow(), Getcol(), Getnum(), "complex");
        realdp *resultptr = result.ObtainWritePartialData();
        for(integer i = 0; i < length / 2; i++)
        {
            resultptr[i] = Space[i * 2];
        }
        return result;
    };

    Element Element::GetImagInComplex(void) const
    {
        assert(iscomplex && Space != nullptr);
        Element result(Getrow(), Getcol(), Getnum(), "complex");
        realdp *resultptr = result.ObtainWritePartialData();
        for(integer i = 0; i < length / 2; i++)
        {
            resultptr[i] = Space[i * 2 + 1];
        }
        return result;
    };

    void Element::SetToIdentity(void)
    {
        SetToZeros();
        integer minmn = (Getrow() < Getcol()) ? Getrow() : Getcol();
        if(!iscomplex)
        {
            for(integer i = 0; i < minmn; i++)
                    Space[i + i * size[0]] = 1;
            return;
        }
        for(integer i = 0; i < minmn; i++)
                Space[2 * i + i * size[0]] = 1;
    };

    Element::~Element(void)
    {
        RemoveAllFromFields();
        DeleteMultiElements();
    };

    void Element::CopyTo(Element &eta) const
    {
        SmartSpace::CopyTo(eta);
        eta.Setiscomplex(iscomplex);
        eta.ResetMultiElementsParams(powsinterval, numoftypes, elements, numofelements);
        for(integer i = 0; i < numoftypes; i++)
        {
            for(integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
            {
                elements[j].CopyFieldsTo(eta.GetElement(j));
            }
        }
        
        MAP::const_iterator thisiter = Fields.begin();
        MAP::const_iterator etaiter, etaiterpre;
        for (thisiter = Fields.begin(); thisiter != Fields.end(); thisiter++)
        {
            etaiter = eta.Fields.find(thisiter->first);
            if (etaiter != eta.Fields.end())
            {
                thisiter->second->CopyTo(*(etaiter->second));
            }
            else
            {
                Element Temp(*(thisiter->second));
                eta.AddToFields(thisiter->first, Temp);
            }
        }
        if (Fields.size() < eta.Fields.size())
        {
            for (etaiter = eta.Fields.begin(); etaiter != eta.Fields.end();)
            {
                thisiter = Fields.find(etaiter->first);
                if (thisiter == Fields.end())
                {
                    etaiterpre = etaiter;
                    etaiter++;
                    eta.RemoveFromFields(etaiterpre->first);
                }
                else
                    etaiter++;
            }
        }
    };

    void Element::CopyFieldsTo(Element &eta) const
    {
        MAP::const_iterator thisiter = Fields.begin();
        MAP::const_iterator etaiter, etaiterpre;
        for (thisiter = Fields.begin(); thisiter != Fields.end(); thisiter++)
        {
            etaiter = eta.Fields.find(thisiter->first);
            if (etaiter != eta.Fields.end())
            {
                thisiter->second->CopyTo(*(etaiter->second));
            }
            else
            {
                Element Temp(*(thisiter->second));
                eta.AddToFields(thisiter->first, Temp);
            }
        }
    };

    void Element::RandUnform(realdp start, realdp end)
    {
        RemoveAllFromFields();
        SmartSpace::RandUnform(start, end);
        ResetMultiElementsParams(powsinterval, numoftypes, elements, numofelements);
    };

    void Element::RandGaussian(realdp mean, realdp variance)
    {
        RemoveAllFromFields();
        SmartSpace::RandGaussian(mean, variance);
        ResetMultiElementsParams(powsinterval, numoftypes, elements, numofelements);
    };

    void Element::SetToZeros(void)
    {
        RemoveAllFromFields();
        SmartSpace::SetToZeros();
        ResetMultiElementsParams(powsinterval, numoftypes, elements, numofelements);
    };

    realdp *Element::ObtainWriteEntireData(void)
    {
        RemoveAllFromFields();
        return SmartSpace::ObtainWriteEntireData();
    };

    realdp *Element::ObtainWritePartialData(void)
    {
        RemoveAllFromFields();
        return SmartSpace::ObtainWritePartialData();
    };

    void Element::NewMemoryOnWrite(void)
    {
        RemoveAllFromFields();
        SmartSpace::NewMemoryOnWrite();
        ResetMultiElementsParams(powsinterval, numoftypes, elements, numofelements);
    };

    void Element::CopyOnWrite(void)
    {
        SmartSpace::CopyOnWrite();
        ResetMultiElementsParams(powsinterval, numoftypes, elements, numofelements);
    };

    void Element::Print(const char *name, bool isonlymain) const
    {
        if(numoftypes != 0)
            printf("=================Product data: %s=========================\n", name);
        
        if (Fields.size() > 0 && !isonlymain)
            printf("=================Main data: %s=========================\n", name);
        
        printf("%s:", name);
        std::cout << "(" << this->Getrow() << ", " << this->Getcol() << ", " << this->Getnum() << ")" << std::endl;
        if(numoftypes != 0)
        {
            for(integer i = 0; i < numofelements; i++)
            {
                std::stringstream strStream;
                strStream << "number " << i << " manifold";
                std::string str = strStream.str();
                elements[i].Print(str.c_str(), isonlymain);
            }
        }
        else
            std::cout << (*this) << std::endl;

        if (Fields.size() > 0 && !isonlymain)
        {
            MAP::const_iterator thisiter;
            for (thisiter = Fields.begin(); thisiter != Fields.end(); thisiter++)
            {
                printf("=================Temp data in %s ================\n", name);
                printf("%s:", thisiter->first.c_str());
                std::cout << "(" << (*thisiter->second).Getrow() << ", " << (*thisiter->second).Getcol() << ", " << (*thisiter->second).Getnum() << ")" << std::endl;
                

                if((*thisiter->second).Getnumoftypes() != 0)
                {
                    for(integer i = 0; i < (*thisiter->second).Getnumofelements(); i++)
                    {
                        std::stringstream strStream;
                        strStream << "number " << i << " manifold";
                        std::string str = strStream.str();
                        (*thisiter->second).Getelements()[i].Print(str.c_str(), isonlymain);
                    }
                }
                else
                    std::cout << (*thisiter->second) << std::endl;
                
            }
            printf("=================end of output: %s=========================\n", name);
        }
    };

    void Element::PrintSize(const char *name, bool isonlymain) const
    {
        if (Fields.size() > 0 && !isonlymain)
            printf("=================Main data: %s=========================\n", name);
        printf("%s:", name);

        std::cout << "(" << this->Getrow() << ", " << this->Getcol() << ", " << this->Getnum() << ")" << std::endl;

        if (Fields.size() > 0 && !isonlymain)
        {
            MAP::const_iterator thisiter;
            for (thisiter = Fields.begin(); thisiter != Fields.end(); thisiter++)
            {
                printf("=================Temp data in %s ================\n", name);
                printf("%s:", thisiter->first.c_str());

                std::cout << "(" << (*thisiter->second).Getrow() << ", " << (*thisiter->second).Getcol() << ", " << (*thisiter->second).Getnum() << ")" << std::endl;
            }
            printf("=================end of output: %s=========================\n", name);
        }
    };

    void Element::AddToFields(std::string name, const Element &Temp) const
    {
        MAP::iterator thisiter;
        thisiter = Fields.find(name);
        if (thisiter == Fields.end())
        {
            Element *Tmp = new Element(Temp);
            Fields.insert(std::pair<std::string, Element *>(name, Tmp));
        }
        else
        {
            Temp.CopyTo(*(thisiter->second));
        }
    };

    void Element::AddToFields(std::string name, Element &Temp) const
    {
        MAP::iterator thisiter;
        thisiter = Fields.find(name);
        if (thisiter == Fields.end())
        {
            Element *Tmp = new Element(Temp);
            Fields.insert(std::pair<std::string, Element *>(name, Tmp));
        }
        else
        {
            Temp.CopyTo(*(thisiter->second));
        }
    };

    Element &Element::Field(std::string name) const
    {
        MAP::iterator thisiter;
        Element *thisptr = const_cast<Element *> (this);
        thisiter = thisptr->Fields.find(name);
        
        if (thisiter == Fields.end())
        {
            printf("Error: Fields %s does not exist!\n", name.c_str());
            assert(false);
            return (*thisptr);
        }

        return (*(thisiter->second));
    };

    void Element::RemoveFromFields(std::string name) const
    {
        MAP::iterator thisiter;
        thisiter = Fields.find(name);
        if (thisiter != Fields.end())
        {
            delete thisiter->second;
            Fields.erase(thisiter);
        }
    };

    void Element::RemoveAllFromFields(void) const
    {
        MAP::iterator thisiter;
        for (thisiter = Fields.begin(); thisiter != Fields.end(); thisiter++)
        {
            delete thisiter->second;
        }
        Fields.clear();
        
        if(elements != nullptr)
        {
            for(integer i = 0; i < numofelements; i++)
                elements[i].RemoveAllFromFields();
        }
    };

    bool Element::FieldsExist(std::string name) const
    {
        MAP::const_iterator thisiter;
        thisiter = Fields.find(name);
        if (thisiter != Fields.end())
        {
            return true;
        }
        return false;
    };

    void Element::ScaledIdOPE(realdp scalar)
    {
        if(ls < 2 || size[0] != size[1])
        {
            printf("Warning: This is not a square matrix. It can not be assigned to be an identity matrix!\n");
            return;
        }
        SmartSpace::NewMemoryOnWrite();
        integer ell = size[0];
        for (integer i = 0; i < ell; i++)
        {
            Space[i + i * ell] = scalar;
            for (integer j = i + 1; j < ell; j++)
            {
                Space[i + j * ell] = 0;
                Space[j + i * ell] = 0;
            }
        }
    };

    void Element::Delete(void)
    {
        RemoveAllFromFields();
        SmartSpace::Delete();
    };

    void Element::ObtainTempNames(std::string *names) const
    {
        MAP::const_iterator thisiter;
        integer idx = 0;
        for (thisiter = Fields.begin(); thisiter != Fields.end(); thisiter++, idx++)
        {
            names[idx].assign(thisiter->first);
        }
    };

    void Element::ResetMultiElementsParams(const integer *inpowsinterval, integer innumoftypes, const Element *inelements, integer innumofelements)
    {
        if(innumoftypes == 0)
        {
            DeleteMultiElements();
            numoftypes = 0; powsinterval = nullptr;
            numofelements = 0; elements = nullptr;
            return;
        }
        
        if(powsinterval == inpowsinterval && elements == inelements) /*if the input parameters belong to "this element", then only reset the space in elements*/
        {
            integer shift = 0;
            integer *isharedtimes = nullptr;
            for(integer i = 0; i < numoftypes; i++)
            {
                for(integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
                {
                    isharedtimes = const_cast<integer *>(elements[j].GetSharedTimes());
                    if(isharedtimes == nullptr)
                    {
                        isharedtimes = new integer;
                        *isharedtimes = 1;
                    }
                    
                    elements[j].SetByParams(isharedtimes, Space + shift);
                    elements[j].ResetMultiElementsParams(inelements[j].Getpowsinterval(), inelements[j].Getnumoftypes(), inelements[j].Getelements(), inelements[j].Getnumofelements());
                    shift += ((Space == nullptr) ? 0 : elements[j].Getlength());
                }
            }
            return;
        }
        
        bool issamesize = true;
        if(numoftypes != innumoftypes || numofelements != innumofelements)
            issamesize = false;
        
        if(issamesize)
        {
            for(integer i = 0; i < innumofelements; i++)
            {
                if(inelements[i].Getlength() != elements[i].Getlength())
                {
                    issamesize = false;
                    break;
                }
            }
        }
        
        /*if input parameters do not belong to "this" element, but have the same size, then only update sharedtimes and elements*/
        if(issamesize)
        {
            integer shift = 0;
            integer *isharedtimes = nullptr;
            for(integer i = 0; i < numoftypes; i++)
            {
                for(integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
                {
                    isharedtimes = const_cast<integer *>(elements[j].GetSharedTimes());
                    if(isharedtimes == nullptr && Space != nullptr)
                    {
                        isharedtimes = new integer;
                        *isharedtimes = 1;
                    }
                    if(Space == nullptr && isharedtimes != nullptr)
                    {
                        delete isharedtimes;
                        isharedtimes = nullptr;
                    }
                    
                    elements[j].SetByParams(isharedtimes, Space + shift);
                    elements[j].ResetMultiElementsParams(elements[j].Getpowsinterval(), elements[j].Getnumoftypes(), elements[j].Getelements(), elements[j].Getnumofelements());
                    shift += ((Space == nullptr) ? 0 : elements[j].Getlength());
                }
            }
            return;
        }

        if(!issamesize)
        {
            DeleteMultiElements();
            
            numoftypes = innumoftypes;
            powsinterval = new integer[numoftypes + 1];
            for(integer i = 0; i < numoftypes + 1; i++)
                powsinterval[i] = inpowsinterval[i];
            
            numofelements = innumofelements;
            elements = new Element [numofelements];
            integer shift = 0;
            integer *isharedtimes = nullptr;
            integer *isize = nullptr;
            integer ils;
            for(integer i = 0; i < numoftypes; i++)
            {
                ils = inelements[i].Getls();
                isize = new integer[ils];
                
                for(integer j = 0; j < ils; j++)
                {
                    isize[j] = inelements[powsinterval[i]].Getsize()[j];
                }
                for(integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
                {
                    if(Space == nullptr)
                    {
                        isharedtimes = nullptr;
                    }
                    else
                    {
                        isharedtimes = new integer;
                        *isharedtimes = 1;
                    }
					delete[] elements[j].Getsize();
                    
                    elements[j].SetByParams(isize, ils, inelements[j].Getlength(), isharedtimes, Space + shift);
                    elements[j].Setiscomplex(inelements[j].Getiscomplex());
                    elements[j].ResetMultiElementsParams(inelements[j].Getpowsinterval(), inelements[j].Getnumoftypes(), inelements[j].Getelements(), inelements[j].Getnumofelements());
                    shift += ((Space == nullptr) ? 0 : elements[j].Getlength());
                }
            }
        }
    };

    void Element::DeleteMultiElements(void)
    {
        if(numoftypes != 0) /*if it is multiple elements, then delete related variables*/
        {
            for (integer i = 0; i < numoftypes; i++)
            {
                if (elements[powsinterval[i]].Getsize() != nullptr)
                {
                    delete[] elements[powsinterval[i]].Getsize();
                }
                
                for(integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
                {
                    if (elements[j].GetSharedTimes() != nullptr)
                    {
                        (*(elements[j].GetSharedTimes()))--;
                        if((*(elements[j].GetSharedTimes())) == 0)
                            delete elements[j].GetSharedTimes();
                    }
                }
            }
            for(integer i = 0; i < numofelements; i++)
            {
                elements[i].DeleteBySettingNull();
            }
            delete[] elements;
            delete[] powsinterval;
            numoftypes = 0;
            numofelements = 0;
        }
    };
}; /*end of ROPTLIB namespace*/

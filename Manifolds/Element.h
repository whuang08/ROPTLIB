/*
This file defines the class for storage. An object of Element type can be
used as a point on a manifold or a tangent vector. The element supports linear
algebra operations.

SmartSpace --> Element

---- WH
*/

#ifndef ELEMENT_H
#define ELEMENT_H

/*Variable, Vector, and LinearOPE are just Element*/
#define Variable Element
#define Vector Element
#define LinearOPE Element

#include "Others/randgen.h"
#include <cstdarg>
#include <map>
#include <string>
#include "Manifolds/SmartSpace.h"
#include "Others/SparseBLAS/blas_sparse.h"
#include "Others/SparseMatrix.h"
#include "Others/def.h"
#include "Others/fftw/fftw3.h"
#include <sstream>

#undef max

/*Define the namespace*/
namespace ROPTLIB{

    namespace GLOBAL{
        extern integer IZERO, IONE, ITWO;
        extern realdp DZERO, DONE, DTWO, DNONE, DNTWO;
        extern realdpcomplex ZZERO, ZONE, ZTWO, ZNONE;
        extern char *N, *T, *L, *R, *V, *C, *U, *A, *S, *O;
    };

    class Element;

    /*Defined the type MAP*/
    typedef std::map<std::string, Element *> MAP;

    class Element : public SmartSpace{
    public:
        /*Output element*/
        friend std::ostream &operator<<(std::ostream &output, const Element &x);

        /*add two elements entrywise, they have to be both real or complex*/
        friend Element operator+(Element left, Element right);

        /*add a real element to a real number*/
        friend Element operator+(Element left, realdp right);
        
        /*add a complex element to a complex number*/
        friend Element operator+(Element left, realdpcomplex right);
        
        /*add a real number to a real element*/
        friend Element operator+(realdp left, Element right);
        
        /*add a complex number to a complex element*/
        friend Element operator+(realdpcomplex left, Element right);
        
        /*left = left + right*/
        friend Element operator+=(Element left, Element right);
        
        /*substract an element from another. They have to be both real or complex*/
        friend Element operator-(Element left, Element right);
        
        /*a real element minus a real number*/
        friend Element operator-(Element left, realdp right);
        
        /*a complex element minus a complex number*/
        friend Element operator-(Element left, realdpcomplex right);
        
        /*a real number minus a real element*/
        friend Element operator-(realdp left, Element right);
        
        /*a complex number minus a complex element*/
        friend Element operator-(realdpcomplex left, Element right);
        
        /*left = left - right*/
        friend Element operator-=(Element left, Element right);
        
        /*a element times another, they have to be both real or complex
         This is matrix multiplication*/
        friend Element operator*(Element left, Element right);

        /*a real number times a real element, this is an entrywise multiplication*/
        friend Element operator*(realdp value, Element mat);
        
        /*a complex number times a complex element, this is an entrywise multiplication*/
        friend Element operator*(realdpcomplex value, Element mat);

        /*a real element times a real number, this is an entrywise multiplication*/
        friend Element operator*(Element mat, realdp value);
        
        /*a complex element times a complex number, this is an entrywise multiplication*/
        friend Element operator*(Element mat, realdpcomplex value);
        
        /*The left element times an inverse of the right element, they have to be both real or complex*/
        friend Element operator/(Element left, Element right);
        
        /*The inverse of the left element times the right element, they have to be both real or complex*/
        friend Element operator%(Element left, Element right);

        /*entrywise divide an real element by a real number*/
        friend Element operator/(Element mat, realdp value);
        
        /*entrywise divide an complex element by a complex number*/
        friend Element operator/(Element mat, realdpcomplex value);
        
        /*entrywise divide an real number by a real element*/
        friend Element operator/(realdp value, Element mat);
        
        /*entrywise divide an complex number by a complex element*/
        friend Element operator/(realdpcomplex value, Element mat);
        
        /*entrywise comparison*/
        friend Element operator==(realdp value, Element mat);
        
        /*entrywise comparison*/
        friend Element operator==(Element mat, realdp value);
        
        /*entrywise comparison*/
        friend Element operator==(Element left, Element right);
        
        /*entrywise comparison*/
        friend Element operator>(realdp value, Element mat);
        
        /*entrywise comparison*/
        friend Element operator>(Element mat, realdp value);
        
        /*entrywise comparison*/
        friend Element operator>(Element left, Element right);
        
        /*entrywise comparison*/
        friend Element operator>=(realdp value, Element mat);
        
        /*entrywise comparison*/
        friend Element operator>=(Element mat, realdp value);
        
        /*entrywise comparison*/
        friend Element operator>=(Element left, Element right);
        
        /*entrywise comparison*/
        friend Element operator<(realdp value, Element mat);
        
        /*entrywise comparison*/
        friend Element operator<(Element mat, realdp value);
        
        /*entrywise comparison*/
        friend Element operator<(Element left, Element right);
        
        /*entrywise comparison*/
        friend Element operator<=(realdp value, Element mat);
        
        /*entrywise comparison*/
        friend Element operator<=(Element mat, realdp value);
        
        /*entrywise comparison*/
        friend Element operator<=(Element left, Element right);

        /*a sparse matrix times a dense matrix*/
        friend Element operator*(const SparseMatrix &left, Element right);
        
        /*a dense matrix times a sparse matrix*/
        friend Element operator*(Element left, const SparseMatrix &right);

        /*====================================Constructors====================================*/
        
        /*Construct an empty Element. */
        Element(void);
        
        /*Copy Constructor */
        Element(const Element &eta)
        {
            Initialization(3, eta.Getsize()[0], eta.Getsize()[1], eta.Getsize()[2]);
            numoftypes = 0; powsinterval = nullptr;
            numofelements = 0; elements = nullptr;
            eta.CopyTo(*this);
        };
        
        const Element &operator=(const Element &right)
        {
            if(&right == this)
                return *this;
//            this->Print("this", false);//----
//            right.Print("right:", false);//--
            right.CopyTo(*this);
            return *this;
        };

        /*Construct an empty Element with only size information. */
        Element(integer r, integer l = 1, integer n = 1, const char * = "real");
        
        /*Construct an empty Element with only size information. */
        Element(integer r, integer l, integer n, bool iniscomplex);
        
        /*Construct an empty Element with only size information. */
        Element(integer r, integer l, const char * type);
        
        /*Construct an empty Element with only size information. */
        Element(integer r, const char * type);
        
        /*Construct a product of elements
        An example of using this constructor to generate an empty point on St(2, 3)^2 Euc(2) is:
        Variable StieX(3, 2), EucX(2);
        Variable ProdX(2, &StieX, 2, &EucX, 1);
        * The first argument indicates that there are two kinds of manifolds St(2, 3) and Euc(2).
        * The second argement indicates that first kind of point is a point on St(2, 3).
        * The third argument indicates the number for the previous element. The number of St(2, 3) is numofmani1 = 2
        * The fourth argement indicates that second kind of point is a point on Euc(2).
        * The fifth argument indicates the number for the previous element. The number of Euc(2) is numofmani2 = 1
        IMPORTANT: the types of elements must be empty, otherwise, it causes memory leakage. In above example,
        StieX and EucX need be empty. In other words, when using StieX.Print(), the output states the element is empty.
        */
        Element(integer numberoftypes, Element *FirstVar, integer Firstnum, ...);
        
        /* Construct a product of elements
        Create a product of elements that is the same as the above one
        Variable *ETypes = new Variable[2];
        ETypes[0] = StieX; ETypes[1] = EucX;
        integer *inpowsinterval = new integer [3];
        inpowsinterval[0] = 0; inpowsinterval[1] = 2; inpowsinterval[2] = 3;
        Variable ProdX(2, ETypes, inpowsinterval);
        delete[] ETypes; delete[] inpowsinterval;
        */
        Element(integer numberoftypes, Element *ElementsTypes, integer *inpowsinterval);
        
        /*====================================end of Constructors====================================*/
        
        
        /*====================================Functionos that do not update "this" object, and that allocate new memory====================================*/
        /*LU decomposition of the Element if it is a matrix. The LU information is stored in
         the field "_LU" and "_P" */
        void LUdecom(void) const;
        
        /*Cholesky decomposition of the Element if it is a matrix.
        this element = L * L^T (real) or element = L * L^H (complex).
        The low-triangle matrix L is stored in the filed "_L"*/
        void CholDecom(void) const;
        
        /* Solver a linear system L^* X = B, where L is a lower triangle matrix.
        trans = "N": solve L X = B
        trans = "T": solve L^T X = B
        trans = "C": solve L^H X = B */
        Element TriangleLinSol(Element L, char *trans = GLOBAL::N) const;
        
        /*Use *syev_ to compute the eigenvalue decomposition for a real symmetric matrix.
         Use *heev_ to compute the eigenvalue decomposition for a complex Hermitian matrix.
         By default, both eigenvalue and eigenvectors are computed, and the lower triangular part of the matrix is used.
         The eigenvalues are stored in the field "_EigVal" and the eigenvectors are stored in the field "_EigVec"
         JobZ = "N" implies only computing eigenvalues
         JobZ = "V" implies computing both eigenvalues and eigenvectors
         UorL = "L" implies only lower triangular part is used
         UorL = "U" implies only upper triangular part is used*/
        void EigenDecomSym(char *JobZ = GLOBAL::V, char *UorL = GLOBAL::L);
        
        /* Use eigenvalue decomposition to compute the exponential of a symmetric/Hermitian matrix
         UorL = "L" implies only lower triangular part is used
         UorL = "U" implies only upper triangular part is used*/
        Element ExpSym(char *UorL = GLOBAL::L);
        
        /* Use eigenvalue decomposition to compute the logarithm of a symmetric/Hermitian matrix
         UorL = "L" implies only lower triangular part is used
         UorL = "U" implies only upper triangular part is used*/
        Element LogSym(char *UorL = GLOBAL::L);
        
        /* Compute Householder reflections to the element, i.e., Q_n Q_{n-1} ... Q_1 M = R.
        The householder reflectors and R are stored in "_HHR" and tau in householder is stored in "_tau" */
        void HHRDecom(void) const;
        
        /* Return Q * Element or Q^T * Element or Element * Q or Element * Q^T, where T denote conjugate transpose
         Q is H_1 H_2 .. H_k, and H_i is the householder reflector, such as that computed by HHRDecom.
         lapack functions *unmqr_ or *ormqr_ are used.
         Trans: N, Side: L : Q * Element
         Trans: T, Side: L : Q^T * Element
         Trans: C, Side: L : Q^* * Element
         Trans: N, Side: R : Element * Q
         Trans: T, Side: R : Element * Q^T
         Trans: C, Side: R : Element * Q^*
         */
        Element HHRMtp(Element HHR, Element tau, char *Trans = GLOBAL::N, char *Side = GLOBAL::L) const;
        
        /* Compute an ecomonic QR decomposition such that M = QR. M is not a fat matrix.
        The Q and R factors are stored in the field "_Q" and "_R" */
        void QRDecom(void);

        /*SVD decomposition for the matrix. The element = U * S * Vt
         The U, S, V factors are stored in the fields "_U", "_S", "_Vt"*/
        void SVDDecom(void);

        /*Schur form of the element, the matrix of Schur vectors and the Schur form are stored
         in "_SchurVec" and "_SchurForm", respectively. */
        void SchurForm(char *jobvs = GLOBAL::V);
        
        /* solve the Sylevster equation A X + X B = C
        A and B can be a same variable*/
        Element SYL(Element A, Element B);
        
        /*Compute Diag(this) * M or M * Diag(this)
         side "L": Diag(this) * M
         side "R": M * Diag(this)
         Default: side = "L"*/
        Element GetDiagTimesM(Element M, char *side = GLOBAL::L) const;
        
        /*Compute this .* M */
        Element GetHadamardProduct(Element M) const;
        
        /*Compute this ./ M */
        Element GetHadamardDivision(Element M) const;
        
        /*Change the shape of the element*/
        Element GetReshape(integer r, integer c = 1, integer n = 1) const;
        
        /*Compute this + scalar * u * v' */
        Element GetHaddRankone(realdp scalar, Element u, Element v) const;
        
        /*Compute this + scalar * u * v' */
        Element GetHaddRankone(realdpcomplex scalar, Element u, Element v) const;
        
        /*Compute entrywise max(this, u)*/
        Element GetMax(Element u) const;
        
        /*Compute entrywise min(this, u)*/
        Element GetMin(Element u) const;
        
        /*Compute entrywise max(this, value)*/
        Element GetMax(realdp value) const;
        
        /*Compute entrywise min(this, value)*/
        Element GetMin(realdp value) const;
        
        /*Compute entrywise min(this, u)*/
        Element GetAbs(void) const;
        
        /*Compute entrywise sqrt*/
        Element GetSqrt(void) const;
        
        /*return U V^T, this = U S V^T */
        Element GetOrth(void) const;

        /*return the conjugate of this element but preserve this element*/
        Element GetConj(void) const;
        
        /*Compute
         side = "L": (I - scalar * v * v') * this
         side = "R": this * (I - scalar * v * v') */
        Element GetHtimesRankone(realdp scalar, Element v, char * side = GLOBAL::L) const;
        
        /*Compute
         side = "L": (I - scalar * v * v') * this
         side = "R": this * (I - scalar * v * v') */
        Element GetHtimesRankone(realdpcomplex scalar, Element v, char * side = GLOBAL::L) const;
        
        /*Submatrix*/
        Element GetSubmatrix(integer rstart, integer rend, integer cstart, integer cend) const;

#ifdef ROPTLIB_WITH_FFTW
        /*2D FFT:
        direction = FFTW_FORWARD: F * this, where F is the DFT matrix (kronecker from)
        direction = FFTW_BACKWARD: \bar{F} * this, where \bar denotes the conjugate operator.*/
        Element GetFFT2D(int direction) const;
        
        /*2D FFT:
        direction = FFTW_FORWARD: F * this, where F is the DFT matrix (kronecker from)
        direction = FFTW_BACKWARD: \bar{F} * this, where \bar denotes the conjugate operator.
        The output uses the memory in the input "result". */
        Element FFT2D(int direction, Element *result) const;
#endif

        /*Haar wavelet transform*/
        Element GetHaarFWT(void) const;
        
        /*Inverse Haar wavelet transform*/
        Element GetInvHaarFWT(void);
        
        /*Return the transpose of the Element if it is a matrix, without changing the original element
         If the matrix is complex, the it is the conjugate transpose of the element*/
        Element GetTranspose(void) const;
        
        Element GetColNormsSquare(void) const;
        
        Element GetColDotProducts(Element M) const;
        
        Element GetRealToComplex(void) const;
        
        Element GetRealInComplex(void) const;
        
        Element GetImagInComplex(void) const;
        
        /*====================================End of Functionos that do not update "this" object, and that allocate new memory====================================*/
        
        
        /*====================================Functionos that update "this" object or input object, and that do not allocate new memory. (efficient ones)!========================*/
        
        /*Compute Diag(this) * M or M * Diag(this)
         side "L": Diag(this) * M
         side "R": M * Diag(this)
         The matrix M is updated and the result is stored in M
         Default: side = "L"*/
        Element &DiagTimesM(Element &M, char *side = GLOBAL::L) const;
        
        /*Compute this + scalar * u * v' */
        Element &HaddRankone(realdp scalar, const Element &u, const Element &v);
        
        /*Compute this + scalar * u * v' */
        Element &HaddRankone(realdpcomplex scalar, const Element &u, const Element &v);
        
        /* this = this + alpha * X*/
        Element &AlphaXaddThis(realdp alpha, const Element &X);
        
        /* this = this + alpha * X*/
        Element &AlphaXaddThis(realdpcomplex alpha, const Element &X);
        
        /* this = alpha * this */
        Element &ScalarTimesThis(realdp alpha);
        
        /* this = alpha * this */
        Element &ScalarTimesThis(realdpcomplex alpha);
        
        /* this = alpha * A^? * B^? + beta * this*/
        Element &AlphaABaddBetaThis(realdp alpha, const Element &A, char *transA, const Element &B, char*transB, realdp beta);
        
        /* this = alpha * A^? * B^? + beta * this*/
        Element &AlphaABaddBetaThis(realdpcomplex alpha, const Element &A, char *transA, const Element &B, char*transB, realdpcomplex beta);
        
        /*Compute real(this^T M), T donotes conjuage transpose */
        realdp DotProduct(const Element &M) const;
        
        /*Change the shape of the element*/
        Element &Reshape(integer r, integer c = 1, integer n = 1);
        
        /*return the conjugate of this element and update this element*/
        Element Conj(void);
        
        /*Submatrix*/
        Element SubmatrixAssignment(integer rstart, integer rend, integer cstart, integer cend, const Element &mat);
        
        /*Return the transpose of the Element if it is a matrix, and transpose the original element
         If the matrix is complex, the it is the conjugate transport of the element*/
        Element Transpose(void);
        
        /*Compute the F-norm*/
        realdp Fnorm(void) const;
        
        /*====================================End of Functionos that update "this" object, and that do not allocate new memory. (efficient ones)!===================*/
        
        /*====================================Functionos about manipulations====================================*/
        /*Set the matrix in element to be the identity*/
        void SetToIdentity(void);
        
        /*Beside calling delete function in based class, it also removed all the tempory data.*/
        virtual ~Element(void);

        /*Copy this Element to "eta" Element. After calling this function,
        this Element and "eta" Element will use same space to store data. */
        virtual void CopyTo(Element &eta) const;

        /*Copy all the temp in this Element to "eta" element. After calling this function,
        this Element and "eta" Element will use same space to store shared temp data.*/
        virtual void CopyFieldsTo(Element &eta) const;

        /*Randomly create this Element. In other words, the space will be allocated based
        on the size. Then each entry in the space will be generated by the uniform distribution in [start, end].
        Note that all the temporary data are also removed.*/
        virtual void RandUnform(realdp start = 0, realdp end = 1);

        /*Randomly create this Element. In other words, the space will be allocated based
        on the size. Then each entry in the space will be generated by the normal distribution with mean and variance.
        Note that all the temporary data are also removed*/
        virtual void RandGaussian(realdp mean = 0, realdp variance = 1);
        
        /*Create this Element. All entries are zero*/
        virtual void SetToZeros(void);

        /*Print the data. The string "name" is to mark the output such that user can find the output easily.
        If isonlymain is true, then only output the data without outputing temporary data. Otherwise,
        all the temporary data are also output.*/
        virtual void Print(const char *name = "", bool isonlymain = true) const;
        
        /*Print size of the data. The string "name" is to mark the output such that user can find the output easily.
        If isonlymain is true, then only output the size of the data without outputing the sizes of temporary data. Otherwise,
        all sizes of the temporary data are also output.*/
        virtual void PrintSize(const char *name = "", bool isonlymain = true) const;

        /*Obtain this Element's pointer which points to the data;
        Users are encourage to call this function if they want to overwrite the data without caring about its original data.
        If the data is shared with other Element, then new memory are allocated without copying the data to the new memory.
        Note that all the temporary data are also removed. */
        virtual realdp *ObtainWriteEntireData(void);

        /*Obtain this Element's pointer which points to the data;
        If the data is shared with other Element, then new memory are allocated and the data are copied to the new memory.
        Note that all the temporary data are also removed. */
        virtual realdp *ObtainWritePartialData(void);

        /*If the data is shared with other Element, then new memory are allocated without copying the data to the new memory.
        Note that all the temporary data are also removed. */
        virtual void NewMemoryOnWrite(void);

        /*If the data is shared with other Element, then new memory are allocated and the data are copied to the new memory.
        Note that all the temporary data are also removed.*/
        virtual void CopyOnWrite(void);
        
        /*Add an object of SharedSpace type to this Element with the label name.*/
        virtual void AddToFields(std::string name, const Element &Temp) const;
        
        /*Add an object of SharedSpace type to this Element with the label name.*/
        virtual void AddToFields(std::string name, Element &Temp) const;

        /*Obtain the temp data with name. It is allowed to modify the temp data.*/
        virtual Element &Field(std::string name) const;
        
        /*Remove a temp data with name from this Element. */
        virtual void RemoveFromFields(std::string name) const;

        /*Remove all the temp data from this Element*/
        virtual void RemoveAllFromFields() const;

        /*Check whether a temp data with name exists or not*/
        virtual bool FieldsExist(std::string name) const;

        /*Obtain the names of all the temp data*/
        void ObtainTempNames(std::string *names) const;
        
        /*When the element is a matrix, assign this element to be scalar * identity */
        virtual void ScaledIdOPE(realdp scalar = 1);
        
        virtual void Delete(void);
        
        /*Obtain the number of temp data*/
        inline integer GetSizeofFields(void) const { return static_cast<integer> (Fields.size()); };
        
        inline bool Getiscomplex(void) const { return iscomplex; };
        
        inline integer Getrow(void) const { return (iscomplex ? (size[0] / 2) : size[0]); };
        
        inline integer Getcol(void) const { return size[1]; };
        
        inline integer Getnum(void) const { return size[2]; };
        
        inline Element *Getelements(void) const { return elements; };
        
        inline integer Getnumofelements(void) const { return numofelements; };
        
        inline integer *Getpowsinterval(void) const { return powsinterval; };
        
        inline integer Getnumoftypes(void) const { return numoftypes; };
        
        inline Element &GetElement(integer idx) const { return elements[idx]; };
        
        inline void Setiscomplex(bool iniscomplex) const { if (iniscomplex && ((std::floor(size[0]/2)*2 - size[0]) !=0)) {printf("Warning: cannot set to be complex!\n"); return;}; iscomplex = iniscomplex; };
        
    protected:
        
        /*The mapping which store the information of temp data*/
        mutable MAP Fields;
        
        /*the entries in element are real or complex*/
        mutable bool iscomplex;
        
        /*Reset the memory of all elements of MultiElement such that their pointers to data are consistant.*/
        virtual void ResetMultiElementsParams(const integer *powsinterval, integer numoftypes, const Element *elements, integer numofelements);
        
        integer numoftypes; /*the number of types of elements*/
        integer *powsinterval; /*The length of powsinterval is numoftypes + 1; the number of each type of elements*/
        
        integer numofelements; /*the total number of elements*/
        Element *elements; /*The length of elements is numofelements; the pointers to all the element*/
        
        void DeleteMultiElements(void); /*delete elements and powsinterval for a collection of multiple elements*/
    };
}; /*end of ROPTLIB namespace*/

#endif

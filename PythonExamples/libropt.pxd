
cdef extern from "randgen.h":
    void genrandseed(unsigned int s)
    double genrandnormal()

cdef extern from "SmartSpace.h":
    cdef cppclass SmartSpace:
        const double *ObtainReadData()

cdef extern from "Element.h":
    cdef cppclass Element(SmartSpace):
        pass

cdef extern from "Manifolds/Stiefel/StieVariable.h":
    cdef cppclass StieVariable(Element):
        StieVariable(int n, int p = 1, int num = 1) except +
        void RandInManifold()

cdef extern from "EucVariable.h":
    cdef cppclass EucVariable(Element):
        EucVariable(int r, int l = 1, int n = 1) except +
        void RandInManifold()

cdef extern from "ProductElement.h":
    cdef cppclass ProductElement:
        ProductElement(Element **elements, int numofelements, int *powsinterval, int numoftypes) except +
        void RandInManifold()

cdef extern from "Manifold.h":
    cdef cppclass Variable:
        const double *ObtainReadData()

cdef extern from "Manifold.h":
    cdef cppclass Manifold:
        pass

cdef extern from "Stiefel.h":
    cdef cppclass Stiefel(Manifold):
        Stiefel(int n, int p) except +
        void CheckParams()

cdef extern from "Euclidean.h":
    cdef cppclass Euclidean:
        Euclidean(int r, int c = 1, int n = 1) except +

cdef extern from "ProductManifold.h":
    cdef cppclass ProductManifold:
        ProductManifold(int numberofmanifolds, ...) except +

cdef extern from "Problem.h":
    cdef cppclass Problem:
        pass

cdef extern from "StieBrockett.h":
    cdef cppclass StieBrockett(Problem):
        StieBrockett(double *inB, double *inD, int inn, int inp)  except +
        void SetDomain(Manifold *inDomain) except +


cdef extern from "Solvers.h":
    cdef enum DEBUGINFO:
        NOOUTPUT, FINALRESULT, ITERRESULT, DETAILED, DEBUGLENGTH


cdef extern from "RTRNewton.h":
    cdef cppclass RTRNewton:
        # RTRNewton(const Problem *prob, const Variable *initialx, const Variable *insoln = NULL) except +
        RTRNewton(const Problem *prob, const Element *initialx) except +
        void CheckParams()
        const Variable *GetXopt()
        void Run()
        int Debug
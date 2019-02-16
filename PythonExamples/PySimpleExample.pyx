# distutils: language = c++


from libropt cimport StieVariable
from libropt cimport  Stiefel, StieBrockett, RTRNewton, DEBUGINFO, FINALRESULT
from libropt cimport genrandseed, genrandnormal
import numpy as np

cdef int py_simple_example():
    # // choose a random seed based on current clock
    genrandseed(0)

    # // size of the Stiefel manifold
    cdef int n = 12
    cdef p = 8
    cdef int i, j
    # // Generate the matrices in the Brockett problem.
    # cdef double *B = new double[n * n + p]
    B = np.zeros(n * n + p)
    cdef double[::1] arr_B = B
    # B = new double[n * n + p]
    # cdef double *D = B + n * n
    # D = &B[n**2]
    D = B[n**2:]
    cdef double[::1] arr_D = D
    # for i in range(n):
    #     for j in range(i, n):
    #         B[i + j * n] = genrandnormal()
    #         B[j + i * n] = B[i + j * n]
    #
    # for i in range(p):
    #     D[i] = static_cast<double> (i + 1)
    for i in range(n):
        for j in range(i, n):
            arr_B[i + j * n] = genrandnormal()
            arr_B[j + i * n] = arr_B[i + j * n]

    for i in range(p):
        arr_D[i] = float(i + 1)

    # // Obtain an initial iterate
    # cdef StieVariable StieX = new StieVariable(n, p)
    StieX = new StieVariable(n, p)
    StieX.RandInManifold()

    # // Define the Stiefel manifold
    # cdef Stiefel Domain = Stiefel(n, p)
    Domain = new Stiefel(n, p)

    # // Define the Brockett problem
    # cdef StieBrockett Prob = StieBrockett(B, D, n, p)
    Prob = new StieBrockett(&arr_B[0], &arr_D[0], n, p)

    # // Set the domain of the problem to be the Stiefel manifold
    Prob.SetDomain(Domain)

    # // output the parameters of the manifold of domain
    Domain.CheckParams()
    # // test RTRNewton
    print("********************************Check RTRNewton*************************************\n")
	#RTRNewton RTRNewtonsolver = RTRNewton(&Prob, &StieX)
    RTRNewtonsolver = new RTRNewton(Prob, StieX) #local variables do not need to explicitly typed (as it would be above)
    # cdef RTRNewton rtrnewton_solver = RTRNewton(&Prob, &StieX)
	# RTRNewtonsolver.Debug = FINALRESULT
    RTRNewtonsolver.CheckParams()
    RTRNewtonsolver.Run()

    # // Check gradient and Hessian
    # Prob.CheckGradHessian(&StieX)
    xopt = RTRNewtonsolver.GetXopt()


    # return result as numpy array
    out = np.zeros(n*p)
    cdef double[::1] arr_out = out
    xoptptr = xopt.ObtainReadData()
    for i in range(n*p):
        arr_out[i] = xoptptr[i]

    del B
    # very import to rememeber that data is stored in Fortran order
    result = out.reshape((n,p), order='F')
    return result
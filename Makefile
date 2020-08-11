# Makefile for ROPTLIB. Test on ubuntu 18.04 LTS

# set compiler
CC = g++
CXXFLAGS:=-O3 -ffastmath -march=native -ggdb3

# default test problem is the Brockett cost function on the Stiefel manifold
TP?=DriverCpp

#the path of ROPTLIB
ROOTPATH:=$(dir $(realpath $(firstword $(MAKEFILE_LIST))))

# set the path of Julia
# JULIA_DIR:=/home/ubuntu/Documents/julia-1.3.1
JULIA_DIR:=/Applications/Julia-1.3.app/Contents/Resources/julia/

# directories of ROPTLIB header files
INCDIRS = -I$(ROOTPATH)/
INCDIRS += -I$(ROOTPATH)/BinaryFiles/
INCDIRS += -I$(ROOTPATH)/Julia/
INCDIRS += -I$(ROOTPATH)/Manifolds/
INCDIRS += -I$(ROOTPATH)/Matlab/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/
INCDIRS += -I$(ROOTPATH)/Matlab/ForMatlab/
INCDIRS += -I$(ROOTPATH)/Others/
INCDIRS += -I$(ROOTPATH)/Others/SparseBLAS/
INCDIRS += -I$(ROOTPATH)/Others/fftw/
INCDIRS += -I$(ROOTPATH)/Others/wavelet/
INCDIRS += -I$(ROOTPATH)/Problems/
INCDIRS += -I$(ROOTPATH)/Solvers/
INCDIRS += -I$(ROOTPATH)/cwrapper/
INCDIRS += -I$(ROOTPATH)/cwrapper/blas/
INCDIRS += -I$(ROOTPATH)/cwrapper/lapack/
INCDIRS += -I$(ROOTPATH)/test/
# ROPTLIB C++ files
CPPS += $(ROOTPATH)/Manifolds/CFixedRankQ2F.cpp $(ROOTPATH)/Manifolds/CStiefel.cpp $(ROOTPATH)/Manifolds/CSymFixedRankQ.cpp $(ROOTPATH)/Manifolds/Element.cpp $(ROOTPATH)/Manifolds/Euclidean.cpp $(ROOTPATH)/Manifolds/FixedRankE.cpp $(ROOTPATH)/Manifolds/FixedRankQ2F.cpp $(ROOTPATH)/Manifolds/Grassmann.cpp $(ROOTPATH)/Manifolds/Manifold.cpp $(ROOTPATH)/Manifolds/MultiManifolds.cpp $(ROOTPATH)/Manifolds/SPDManifold.cpp $(ROOTPATH)/Manifolds/SmartSpace.cpp $(ROOTPATH)/Manifolds/Sphere.cpp $(ROOTPATH)/Manifolds/SphereTx.cpp $(ROOTPATH)/Manifolds/Stiefel.cpp $(ROOTPATH)/Manifolds/SymFixedRankQ.cpp 
CPPS += $(ROOTPATH)/Others/BlasLapackCppWrapper.cpp $(ROOTPATH)/Others/ForDebug.cpp $(ROOTPATH)/Others/MinPNormConHull.cpp $(ROOTPATH)/Others/SparseMatrix.cpp $(ROOTPATH)/Others/Spline.cpp $(ROOTPATH)/Others/Timer.cpp $(ROOTPATH)/Others/randgen.cpp 
CPPS += $(ROOTPATH)/Others/SparseBLAS/nist_spblas.cpp 
CPPS += $(ROOTPATH)/Others/wavelet/wavelet.cpp 
CPPS += $(ROOTPATH)/Problems/CFRankQ2FBlindDecon2D.cpp $(ROOTPATH)/Problems/CSFRQPhaseRetrieval.cpp $(ROOTPATH)/Problems/CStieBrockett.cpp $(ROOTPATH)/Problems/EucQuadratic.cpp $(ROOTPATH)/Problems/FRankESparseApprox.cpp $(ROOTPATH)/Problems/FRankEWeightApprox.cpp $(ROOTPATH)/Problems/FRankQ2FMatCompletion.cpp $(ROOTPATH)/Problems/GrassRQ.cpp $(ROOTPATH)/Problems/Problem.cpp $(ROOTPATH)/Problems/ProdStieSumBrockett.cpp $(ROOTPATH)/Problems/SFRQLyapunov.cpp $(ROOTPATH)/Problems/SPDKarcherMean.cpp $(ROOTPATH)/Problems/SphereConvexHull.cpp $(ROOTPATH)/Problems/SphereSparsestVector.cpp $(ROOTPATH)/Problems/SphereTxRQ.cpp $(ROOTPATH)/Problems/StieBrockett.cpp $(ROOTPATH)/Problems/StieSPCA.cpp $(ROOTPATH)/Problems/mexProblem.cpp $(ROOTPATH)/Problems/juliaProblem.cpp 
CPPS += $(ROOTPATH)/Solvers/AManPG.cpp $(ROOTPATH)/Solvers/LRBFGS.cpp $(ROOTPATH)/Solvers/LRBFGSSub.cpp $(ROOTPATH)/Solvers/LRTRSR1.cpp $(ROOTPATH)/Solvers/ManPG.cpp $(ROOTPATH)/Solvers/RBFGS.cpp $(ROOTPATH)/Solvers/RBFGSSub.cpp $(ROOTPATH)/Solvers/RBroydenFamily.cpp $(ROOTPATH)/Solvers/RCG.cpp $(ROOTPATH)/Solvers/RGS.cpp $(ROOTPATH)/Solvers/RNewton.cpp $(ROOTPATH)/Solvers/RSD.cpp $(ROOTPATH)/Solvers/RTRNewton.cpp $(ROOTPATH)/Solvers/RTRSD.cpp $(ROOTPATH)/Solvers/RTRSR1.cpp $(ROOTPATH)/Solvers/RWRBFGS.cpp $(ROOTPATH)/Solvers/Solvers.cpp $(ROOTPATH)/Solvers/SolversNSM.cpp $(ROOTPATH)/Solvers/SolversNSMPGLS.cpp $(ROOTPATH)/Solvers/SolversNSMSub.cpp $(ROOTPATH)/Solvers/SolversNSMSubLS.cpp $(ROOTPATH)/Solvers/SolversSMLS.cpp $(ROOTPATH)/Solvers/SolversSMTR.cpp $(ROOTPATH)/Solvers/SolversSM.cpp 
CPPS += $(ROOTPATH)/test/DriverMexProb.cpp $(ROOTPATH)/test/TestCFRankQ2FBlindDecon2D.cpp $(ROOTPATH)/test/TestCSFRQPhaseRetrieval.cpp $(ROOTPATH)/test/TestCStieBrockett.cpp $(ROOTPATH)/test/TestElement.cpp $(ROOTPATH)/test/TestEucQuadratic.cpp $(ROOTPATH)/test/TestFRankESparseApprox.cpp $(ROOTPATH)/test/TestFRankEWeightApprox.cpp $(ROOTPATH)/test/TestFRankQ2FMatCompletion.cpp $(ROOTPATH)/test/TestGrassRQ.cpp $(ROOTPATH)/test/TestProdStieSumBrockett.cpp $(ROOTPATH)/test/TestSFRQLyapunov.cpp $(ROOTPATH)/test/TestSPDKarcherMean.cpp $(ROOTPATH)/test/TestSphereSparsestVector.cpp $(ROOTPATH)/test/TestStieBrockett.cpp $(ROOTPATH)/test/TestStieSPCA.cpp 
# convert a string to upper case.
UPPER_TP  = $(shell echo $(TP) | tr a-z A-Z)

# make a binary file, which is called in command line
ROPTLIB:
	$(CC) -O3 -w -std=c++0x $(ROOTPATH)/test/$(TP).cpp $(CPPS) $(INCDIRS) -D$(UPPER_TP) -DROPTLIB_WITH_FFTW -llapack -lblas -lfftw3 -lm -o $(TP)

#make a library
libropt.so:
	$(CC) -w -std=c++0x -shared -fPIC -O3 $(CPPS) $(INCDIRS) -DROPTLIB_WITH_FFTW -llapack -lblas -lfftw3 -lm -o $@

JULIA_LIB:=$(JULIA_DIR)/lib
JULIA_SRC:=$(JULIA_DIR)/src
JULIA_INC:=$(JULIA_DIR)/include/julia
CPPFLAGS:=-I$(JULIA_INC) -I$(JULIA_SRC) -I$(JULIA_SRC)/support
LDFLAGS:=-L$(JULIA_LIB)
LDLIBS=-ljulia
export LD_LIBRARY_PATH:=$(JULIA_LIB):$(JULIA_LIB)/julia

# make a shared library, which is used by Julia
JuliaROPTLIB:
	$(CC) -O3 -shared -fPIC -std=c++0x $(ROOTPATH)/test/$(TP).cpp $(CPPS) $(INCDIRS) -D$(UPPER_TP) $(CPPFLAGS) $(LDFLAGS) -Wl,-rpath,$(JULIA_LIB) -lm $(LDLIBS) -DJULIA_LIB_DIR=\"$(JULIA_DIR)/lib/julia\" -DROPTLIB_WITH_FFTW -llapack -lblas -lfftw3 -o $(TP).so

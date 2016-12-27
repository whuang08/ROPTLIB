# Makefile for ROPTLIB. Test on ubuntu 16.04 LTS

# set compiler
CC = g++

# default test problem is the Brockett cost function on the Stiefel manifold
TP?=TestSimpleExample

#the path of ROPTLIB
ROOTPATH = /home/whuang/Documents/ROPTLIB

# set the path of Julia
JULIA_DIR:=/home/whuang/Documents/julia

# ROPTLIB C++ files
MANIFOLDS = $(ROOTPATH)/Manifolds/Element.cpp $(ROOTPATH)/Manifolds/LinearOPE.cpp $(ROOTPATH)/Manifolds/Manifold.cpp $(ROOTPATH)/Manifolds/ProductElement.cpp $(ROOTPATH)/Manifolds/ProductManifold.cpp $(ROOTPATH)/Manifolds/SharedSpace.cpp $(ROOTPATH)/Manifolds/SmartSpace.cpp $(ROOTPATH)/Manifolds/CpxNStQOrth/CSOVariable.cpp $(ROOTPATH)/Manifolds/CpxNStQOrth/CSOVector.cpp $(ROOTPATH)/Manifolds/CpxNStQOrth/CpxNStQOrth.cpp $(ROOTPATH)/Manifolds/ElasticShape/ElasticShape.cpp $(ROOTPATH)/Manifolds/ElasticShape/ShapeVariable.cpp $(ROOTPATH)/Manifolds/ElasticShape/ShapeVector.cpp $(ROOTPATH)/Manifolds/EucPositive/EucPosVariable.cpp $(ROOTPATH)/Manifolds/EucPositive/EucPosVector.cpp $(ROOTPATH)/Manifolds/EucPositive/EucPositive.cpp $(ROOTPATH)/Manifolds/Euclidean/EucVariable.cpp $(ROOTPATH)/Manifolds/Euclidean/EucVector.cpp $(ROOTPATH)/Manifolds/Euclidean/Euclidean.cpp $(ROOTPATH)/Manifolds/Grassmann/GrassVariable.cpp $(ROOTPATH)/Manifolds/Grassmann/GrassVector.cpp $(ROOTPATH)/Manifolds/Grassmann/Grassmann.cpp $(ROOTPATH)/Manifolds/L2Sphere/L2Sphere.cpp $(ROOTPATH)/Manifolds/L2Sphere/L2SphereVariable.cpp $(ROOTPATH)/Manifolds/L2Sphere/L2SphereVector.cpp $(ROOTPATH)/Manifolds/LowRank/LowRank.cpp $(ROOTPATH)/Manifolds/LowRank/LowRankVariable.cpp $(ROOTPATH)/Manifolds/LowRank/LowRankVector.cpp $(ROOTPATH)/Manifolds/Oblique/Oblique.cpp $(ROOTPATH)/Manifolds/Oblique/ObliqueVariable.cpp $(ROOTPATH)/Manifolds/Oblique/ObliqueVector.cpp $(ROOTPATH)/Manifolds/OrthGroup/OrthGroup.cpp $(ROOTPATH)/Manifolds/OrthGroup/OrthGroupVariable.cpp $(ROOTPATH)/Manifolds/OrthGroup/OrthGroupVector.cpp $(ROOTPATH)/Manifolds/PreShapeCurves/PSCVariable.cpp $(ROOTPATH)/Manifolds/PreShapeCurves/PSCVector.cpp $(ROOTPATH)/Manifolds/PreShapeCurves/PreShapeCurves.cpp $(ROOTPATH)/Manifolds/SPDManifold/SPDManifold.cpp $(ROOTPATH)/Manifolds/SPDManifold/SPDVariable.cpp $(ROOTPATH)/Manifolds/SPDManifold/SPDVector.cpp $(ROOTPATH)/Manifolds/SPDTensor/SPDTVariable.cpp $(ROOTPATH)/Manifolds/SPDTensor/SPDTVector.cpp $(ROOTPATH)/Manifolds/SPDTensor/SPDTensor.cpp $(ROOTPATH)/Manifolds/Sphere/Sphere.cpp $(ROOTPATH)/Manifolds/Sphere/SphereVariable.cpp $(ROOTPATH)/Manifolds/Sphere/SphereVector.cpp $(ROOTPATH)/Manifolds/Stiefel/StieVariable.cpp $(ROOTPATH)/Manifolds/Stiefel/StieVector.cpp $(ROOTPATH)/Manifolds/Stiefel/Stiefel.cpp

# ROPTLIB C++ files
OTHERS = $(ROOTPATH)/Others/ForDebug.cpp $(ROOTPATH)/Others/MinPNormConHull.cpp $(ROOTPATH)/Others/MyMatrix.cpp $(ROOTPATH)/Others/Spline.cpp $(ROOTPATH)/Others/Timer.cpp $(ROOTPATH)/Others/randgen.cpp $(ROOTPATH)/Others/nist_spblas.cpp 

# ROPTLIB C++ files
PROBLEMS = $(ROOTPATH)/Problems/Problem.cpp $(ROOTPATH)/Problems/mexProblem.cpp $(ROOTPATH)/Problems/juliaProblem.cpp $(ROOTPATH)/Problems/ElasticCurvesRO/DriverElasticCurvesRO.cpp $(ROOTPATH)/Problems/ElasticCurvesRO/ElasticCurvesRO.cpp $(ROOTPATH)/Problems/EucFrechetMean/EucFrechetMean.cpp $(ROOTPATH)/Problems/EucPosSpCd/EucPosSpCd.cpp $(ROOTPATH)/Problems/EucQuadratic/EucQuadratic.cpp $(ROOTPATH)/Problems/GrassRQ/GrassRQ.cpp $(ROOTPATH)/Problems/KarcherMean/KarcherMean.cpp $(ROOTPATH)/Problems/ObliqueTestSparsePCA/ObliqueTestSparsePCA.cpp $(ROOTPATH)/Problems/OrthBoundingBox/OrthBoundingBox.cpp $(ROOTPATH)/Problems/PreShapePathStraighten/PreShapePathStraighten.cpp $(ROOTPATH)/Problems/SPDMean/SPDMean.cpp $(ROOTPATH)/Problems/SPDTensorDL/SPDTensorDL.cpp $(ROOTPATH)/Problems/ShapePathStraighten/ShapePathStraighten.cpp $(ROOTPATH)/Problems/SphereConvexHull/SphereConvexHull.cpp $(ROOTPATH)/Problems/StieBrockett/StieBrockett.cpp $(ROOTPATH)/Problems/StieSoftICA/StieSoftICA.cpp $(ROOTPATH)/Problems/StieSparseBrockett/StieSparseBrockett.cpp $(ROOTPATH)/Problems/StieSumBrockett/StieSumBrockett.cpp $(ROOTPATH)/Problems/WeightedLowrank/WeightedLowRank.cpp 

# ROPTLIB C++ files
SOLVERS = $(ROOTPATH)/Solvers/LRBFGS.cpp $(ROOTPATH)/Solvers/LRBFGSLPSub.cpp $(ROOTPATH)/Solvers/LRTRSR1.cpp $(ROOTPATH)/Solvers/MRankAdaptive.cpp $(ROOTPATH)/Solvers/QuasiNewton.cpp $(ROOTPATH)/Solvers/RBFGS.cpp $(ROOTPATH)/Solvers/RBFGSLPSub.cpp $(ROOTPATH)/Solvers/RBroydenFamily.cpp $(ROOTPATH)/Solvers/RCG.cpp $(ROOTPATH)/Solvers/RGS.cpp $(ROOTPATH)/Solvers/RNewton.cpp $(ROOTPATH)/Solvers/RSD.cpp $(ROOTPATH)/Solvers/RTRNewton.cpp $(ROOTPATH)/Solvers/RTRSD.cpp $(ROOTPATH)/Solvers/RTRSR1.cpp $(ROOTPATH)/Solvers/RWRBFGS.cpp $(ROOTPATH)/Solvers/Solvers.cpp $(ROOTPATH)/Solvers/SolversLS.cpp $(ROOTPATH)/Solvers/SolversLSLPSub.cpp $(ROOTPATH)/Solvers/SolversTR.cpp 

# directories of ROPTLIB header files
INCDIRS = -I$(ROOTPATH)/ -I$(ROOTPATH)/BinaryFiles/ -I$(ROOTPATH)/Manifolds/ -I$(ROOTPATH)/Manifolds/CpxNStQOrth/ -I$(ROOTPATH)/Manifolds/ElasticShape/ -I$(ROOTPATH)/Manifolds/EucPositive/ -I$(ROOTPATH)/Manifolds/Euclidean/ -I$(ROOTPATH)/Manifolds/Grassmann/ -I$(ROOTPATH)/Manifolds/L2Sphere/ -I$(ROOTPATH)/Manifolds/LowRank/ -I$(ROOTPATH)/Manifolds/Oblique/ -I$(ROOTPATH)/Manifolds/OrthGroup/ -I$(ROOTPATH)/Manifolds/PreShapeCurves/ -I$(ROOTPATH)/Manifolds/SPDManifold/ -I$(ROOTPATH)/Manifolds/SPDTensor/ -I$(ROOTPATH)/Manifolds/Sphere/ -I$(ROOTPATH)/Manifolds/Stiefel/ -I$(ROOTPATH)/Matlab/ -I$(ROOTPATH)/Matlab/ForCpp/ -I$(ROOTPATH)/Matlab/ForCpp/Boundingbox/ -I$(ROOTPATH)/Matlab/ForCpp/BrockettLRBFGSVTpaper/ -I$(ROOTPATH)/Matlab/ForCpp/BrockettNonconvex/ -I$(ROOTPATH)/Matlab/ForCpp/RepaRotCurves/ -I$(ROOTPATH)/Matlab/ForCpp/SPDTensorDLandSC/ -I$(ROOTPATH)/Matlab/ForCpp/SPDTensorDLandSC/DLandSC/ -I$(ROOTPATH)/Matlab/ForCpp/SPDTensorDLandSC/EucPosSC/ -I$(ROOTPATH)/Matlab/ForCpp/SPDTensorDLandSC/EucPosSC/algos/ -I$(ROOTPATH)/Matlab/ForCpp/SPDTensorDLandSC/EucPosSC/tools/ -I$(ROOTPATH)/Matlab/ForCpp/SPDTensorDLandSC/EucPosSC/tools/others/ -I$(ROOTPATH)/Matlab/ForCpp/SPDTensorDLandSC/SPDtensorDL/ -I$(ROOTPATH)/Matlab/ForCpp/SoftICA/ -I$(ROOTPATH)/Matlab/ForCpp/SoftICA/New_folder/ -I$(ROOTPATH)/Matlab/ForCpp/SoftICA/RBFGSNonconvexResults/ -I$(ROOTPATH)/Matlab/ForCpp/SparsePCA/ -I$(ROOTPATH)/Matlab/ForMatlab/ -I$(ROOTPATH)/Matlab/ForMatlab/FromMelissa/ -I$(ROOTPATH)/Others/ -I$(ROOTPATH)/Problems/ -I$(ROOTPATH)/Problems/ElasticCurvesRO/ -I$(ROOTPATH)/Problems/EucFrechetMean/ -I$(ROOTPATH)/Problems/EucPosSpCd/ -I$(ROOTPATH)/Problems/EucQuadratic/ -I$(ROOTPATH)/Problems/GrassRQ/ -I$(ROOTPATH)/Problems/KarcherMean/ -I$(ROOTPATH)/Problems/ObliqueTestSparsePCA/ -I$(ROOTPATH)/Problems/OrthBoundingBox/ -I$(ROOTPATH)/Problems/PreShapePathStraighten/ -I$(ROOTPATH)/Problems/SPDMean/ -I$(ROOTPATH)/Problems/SPDTensorDL/ -I$(ROOTPATH)/Problems/ShapePathStraighten/ -I$(ROOTPATH)/Problems/SphereConvexHull/ -I$(ROOTPATH)/Problems/StieBrockett/ -I$(ROOTPATH)/Problems/StieSoftICA/ -I$(ROOTPATH)/Problems/StieSparseBrockett/ -I$(ROOTPATH)/Problems/StieSumBrockett/ -I$(ROOTPATH)/Problems/WeightedLowrank/ -I$(ROOTPATH)/Solvers/ -I$(ROOTPATH)/test/ -I$(ROOTPATH)/cwrapper/blas/ -I$(ROOTPATH)/cwrapper/lapack/

# convert a string to upper case.
UPPER_TP  = $(shell echo $(TP) | tr a-z A-Z)

# make a binary file, which is called in command line
ROPTLIB:
	$(CC) -O3 -w -std=c++0x $(ROOTPATH)/test/$(TP).cpp $(MANIFOLDS) $(OTHERS) $(PROBLEMS) $(SOLVERS) $(INCDIRS) -D$(UPPER_TP) -llapack -lblas -lm -o $(TP)


JULIA_LIB:=$(JULIA_DIR)/usr/lib
JULIA_SRC:=$(JULIA_DIR)/src
JULIA_INC:=$(JULIA_DIR)/usr/include
CPPFLAGS:=-I$(JULIA_INC) -I$(JULIA_SRC) -I$(JULIA_SRC)/support
LDFLAGS:=-L$(JULIA_LIB)
LDLIBS=-ljulia
export LD_LIBRARY_PATH:=$(JULIA_LIB):$(JULIA_LIB)/julia

# make a shared library, which is used by Julia
JuliaROPTLIB:
	$(CC) -O3 -shared -fPIC -std=c++0x $(ROOTPATH)/test/$(TP).cpp $(MANIFOLDS) $(OTHERS) $(PROBLEMS) $(SOLVERS) $(INCDIRS) -D$(UPPER_TP) $(CPPFLAGS) $(LDFLAGS) -Wl,-rpath,$(JULIA_LIB) -lm $(LDLIBS) -DJULIA_LIB_DIR=\"$(JULIA_DIR)/lib/julia\" -llapack -lblas -o $(TP).so


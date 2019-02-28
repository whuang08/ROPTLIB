# Makefile for ROPTLIB. Test on ubuntu 16.04 LTS

# set compiler
CC = g++
CXXFLAGS:=-O3 -ffastmath -march=native -ggdb3

# default test problem is the Brockett cost function on the Stiefel manifold
TP?=DriverCpp

#the path of ROPTLIB
ROOTPATH:=$(dir $(realpath $(firstword $(MAKEFILE_LIST))))

# set the path of Julia
JULIA_DIR:=/home/whuang/Documents/julia

# directories of ROPTLIB header files
INCDIRS = -I$(ROOTPATH)/
INCDIRS += -I$(ROOTPATH)/BinaryFiles/
INCDIRS += -I$(ROOTPATH)/Julia/
INCDIRS += -I$(ROOTPATH)/Julia/useless/
INCDIRS += -I$(ROOTPATH)/Manifolds/
INCDIRS += -I$(ROOTPATH)/Manifolds/CFixedRank2Factors/
INCDIRS += -I$(ROOTPATH)/Manifolds/CpxNStQOrth/
INCDIRS += -I$(ROOTPATH)/Manifolds/ElasticShape/
INCDIRS += -I$(ROOTPATH)/Manifolds/EucPositive/
INCDIRS += -I$(ROOTPATH)/Manifolds/Euclidean/
INCDIRS += -I$(ROOTPATH)/Manifolds/Grassmann/
INCDIRS += -I$(ROOTPATH)/Manifolds/L2Sphere/
INCDIRS += -I$(ROOTPATH)/Manifolds/LowRank/
INCDIRS += -I$(ROOTPATH)/Manifolds/Oblique/
INCDIRS += -I$(ROOTPATH)/Manifolds/OrthGroup/
INCDIRS += -I$(ROOTPATH)/Manifolds/PreShapeCurves/
INCDIRS += -I$(ROOTPATH)/Manifolds/SPDManifold/
INCDIRS += -I$(ROOTPATH)/Manifolds/SPDTensor/
INCDIRS += -I$(ROOTPATH)/Manifolds/Sphere/
INCDIRS += -I$(ROOTPATH)/Manifolds/SphereTx/
INCDIRS += -I$(ROOTPATH)/Manifolds/Stiefel/
INCDIRS += -I$(ROOTPATH)/Matlab/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/Blinddeconvolution/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/Blinddeconvolution/BlindDeconvFromKeWei/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/Blinddeconvolution/BlindDeconvFromKeWei/2D/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/Blinddeconvolution/BlindDeconvFromKeWei/Aux1/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/Blinddeconvolution/BlindDeconvFromKeWei/Main/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/Blinddeconvolution/BlindDeconv_Convex/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/Blinddeconvolution/BlindDeconv_Convex/Romberg_noiselet/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/Blinddeconvolution/BlindDeconv_Convex/Romberg_noiselet/Measurements/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/Blinddeconvolution/BlindDeconv_Convex/Romberg_noiselet/Optimization/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/Blinddeconvolution/BlindDeconv_Convex/Romberg_noiselet/Utils/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/Blinddeconvolution/BlindDeconv_Convex/minFunc/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/Blinddeconvolution/BlindDeconv_Convex/minFunc/compiled/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/Blinddeconvolution/BlindDeconv_Convex/minFunc/mex/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/Blinddeconvolution/BlindDeconv_Convex/minFunc_2012/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/Blinddeconvolution/BlindDeconv_Convex/minFunc_2012/autoDif/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/Blinddeconvolution/BlindDeconv_Convex/minFunc_2012/logisticExample/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/Blinddeconvolution/Results/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/Blinddeconvolution_basin/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/Boundingbox/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/BrockettLRBFGSVTpaper/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/BrockettNonconvex/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/PhaseRetrieval/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/RepaRotCurves/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/SPDTensorDLandSC/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/SPDTensorDLandSC/DLandSC/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/SPDTensorDLandSC/EucPosSC/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/SPDTensorDLandSC/EucPosSC/algos/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/SPDTensorDLandSC/EucPosSC/tools/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/SPDTensorDLandSC/EucPosSC/tools/others/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/SPDTensorDLandSC/SPDtensorDL/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/SoftICA/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/SoftICA/New_folder/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/SoftICA/RBFGSNonconvexResults/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/SoftICAnonconvexpaper/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/SoftICAnonconvexpaper/nonconvexResults/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/SparsePCA/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/SparsestVector/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/TOMS_ROPTLIB/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/LyapunovLR/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/MatrixCompletion/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/MatrixCompletion/Auxiliary/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/TensorCompletion/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/TensorCompletion/tensors_toolbox/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/TensorCompletion/tensors_toolbox/tensor_toolbox/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/TensorCompletion/tensors_toolbox/tensor_toolbox/doc/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/TensorCompletion/tensors_toolbox/tensor_toolbox/doc/html/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/TensorCompletion/tensors_toolbox/tensor_toolbox/doc/images/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/TensorCompletion/tensors_toolbox/tensor_toolbox/met/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/manopt_2/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/manopt_2/checkinstall/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/manopt_2/examples/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/manopt_2/manopt/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/manopt_2/manopt/core/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/manopt_2/manopt/manifolds/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/manopt_2/manopt/manifolds/complexcircle/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/manopt_2/manopt/manifolds/essential/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/manopt_2/manopt/manifolds/essential/privateessential/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/manopt_2/manopt/manifolds/euclidean/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/manopt_2/manopt/manifolds/fixedrank/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/manopt_2/manopt/manifolds/fixedranktensors/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/manopt_2/manopt/manifolds/grassmann/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/manopt_2/manopt/manifolds/multinomial/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/manopt_2/manopt/manifolds/oblique/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/manopt_2/manopt/manifolds/rotations/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/manopt_2/manopt/manifolds/specialeuclidean/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/manopt_2/manopt/manifolds/sphere/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/manopt_2/manopt/manifolds/stiefel/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/manopt_2/manopt/manifolds/symfixedrank/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/manopt_2/manopt/solvers/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/manopt_2/manopt/solvers/barzilaiborwein/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/manopt_2/manopt/solvers/bfgs/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/manopt_2/manopt/solvers/conjugategradient/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/manopt_2/manopt/solvers/gradientapproximations/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/manopt_2/manopt/solvers/hessianapproximations/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/manopt_2/manopt/solvers/linesearch/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/manopt_2/manopt/solvers/neldermead/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/manopt_2/manopt/solvers/preconditioners/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/manopt_2/manopt/solvers/pso/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/manopt_2/manopt/solvers/steepestdescent/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/manopt_2/manopt/solvers/stochasticgradient/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/manopt_2/manopt/solvers/trustregions/
INCDIRS += -I$(ROOTPATH)/Matlab/ForCpp/WithBartGena/manopt_2/manopt/tools/
INCDIRS += -I$(ROOTPATH)/Matlab/ForMatlab/
INCDIRS += -I$(ROOTPATH)/Matlab/ForMatlab/FromMelissa/
INCDIRS += -I$(ROOTPATH)/Others/
INCDIRS += -I$(ROOTPATH)/Others/SparseBLAS/
INCDIRS += -I$(ROOTPATH)/Others/fftw/
INCDIRS += -I$(ROOTPATH)/Others/wavelet/
INCDIRS += -I$(ROOTPATH)/Problems/
INCDIRS += -I$(ROOTPATH)/Problems/CFR2BlindDecon2D/
INCDIRS += -I$(ROOTPATH)/Problems/CFR2BlindDeconvolution/
INCDIRS += -I$(ROOTPATH)/Problems/CSOPhaseRetrieval/
INCDIRS += -I$(ROOTPATH)/Problems/ElasticCurvesRO/
INCDIRS += -I$(ROOTPATH)/Problems/EucBlindDeconvolution/
INCDIRS += -I$(ROOTPATH)/Problems/EucFrechetMean/
INCDIRS += -I$(ROOTPATH)/Problems/EucPosSpCd/
INCDIRS += -I$(ROOTPATH)/Problems/EucQuadratic/
INCDIRS += -I$(ROOTPATH)/Problems/GrassRQ/
INCDIRS += -I$(ROOTPATH)/Problems/KarcherMean/
INCDIRS += -I$(ROOTPATH)/Problems/LRBlindDeconvolution/
INCDIRS += -I$(ROOTPATH)/Problems/LRMatrixCompletion/
INCDIRS += -I$(ROOTPATH)/Problems/ObliqueSparsePCA/
INCDIRS += -I$(ROOTPATH)/Problems/ObliqueTestSparsePCA/
INCDIRS += -I$(ROOTPATH)/Problems/OrthBoundingBox/
INCDIRS += -I$(ROOTPATH)/Problems/PreShapePathStraighten/
INCDIRS += -I$(ROOTPATH)/Problems/SPDMean/
INCDIRS += -I$(ROOTPATH)/Problems/SPDTensorDL/
INCDIRS += -I$(ROOTPATH)/Problems/ShapePathStraighten/
INCDIRS += -I$(ROOTPATH)/Problems/SphereConvexHull/
INCDIRS += -I$(ROOTPATH)/Problems/SphereSparsestVector/
INCDIRS += -I$(ROOTPATH)/Problems/SphereTxRQ/
INCDIRS += -I$(ROOTPATH)/Problems/StieBrockett/
INCDIRS += -I$(ROOTPATH)/Problems/StieSoftICA/
INCDIRS += -I$(ROOTPATH)/Problems/StieSparseBrockett/
INCDIRS += -I$(ROOTPATH)/Problems/StieSumBrockett/
INCDIRS += -I$(ROOTPATH)/Problems/WeightedLowrank/
INCDIRS += -I$(ROOTPATH)/Solvers/
INCDIRS += -I$(ROOTPATH)/cwrapper/
INCDIRS += -I$(ROOTPATH)/cwrapper/blas/
INCDIRS += -I$(ROOTPATH)/cwrapper/lapack/
INCDIRS += -I$(ROOTPATH)/test/
# ROPTLIB C++ files
CPPS += $(ROOTPATH)/Manifolds/Element.cpp $(ROOTPATH)/Manifolds/LinearOPE.cpp $(ROOTPATH)/Manifolds/Manifold.cpp $(ROOTPATH)/Manifolds/ProductElement.cpp $(ROOTPATH)/Manifolds/ProductManifold.cpp $(ROOTPATH)/Manifolds/SharedSpace.cpp $(ROOTPATH)/Manifolds/SmartSpace.cpp 
CPPS += $(ROOTPATH)/Manifolds/CFixedRank2Factors/CFR2Variable.cpp $(ROOTPATH)/Manifolds/CFixedRank2Factors/CFR2Vector.cpp $(ROOTPATH)/Manifolds/CFixedRank2Factors/CFixedRank2Factors.cpp 
CPPS += $(ROOTPATH)/Manifolds/CpxNStQOrth/CSOVariable.cpp $(ROOTPATH)/Manifolds/CpxNStQOrth/CSOVector.cpp $(ROOTPATH)/Manifolds/CpxNStQOrth/CpxNStQOrth.cpp 
CPPS += $(ROOTPATH)/Manifolds/ElasticShape/ElasticShape.cpp $(ROOTPATH)/Manifolds/ElasticShape/ShapeVariable.cpp $(ROOTPATH)/Manifolds/ElasticShape/ShapeVector.cpp 
CPPS += $(ROOTPATH)/Manifolds/EucPositive/EucPosVariable.cpp $(ROOTPATH)/Manifolds/EucPositive/EucPosVector.cpp $(ROOTPATH)/Manifolds/EucPositive/EucPositive.cpp 
CPPS += $(ROOTPATH)/Manifolds/Euclidean/EucVariable.cpp $(ROOTPATH)/Manifolds/Euclidean/EucVector.cpp $(ROOTPATH)/Manifolds/Euclidean/Euclidean.cpp 
CPPS += $(ROOTPATH)/Manifolds/Grassmann/GrassVariable.cpp $(ROOTPATH)/Manifolds/Grassmann/GrassVector.cpp $(ROOTPATH)/Manifolds/Grassmann/Grassmann.cpp 
CPPS += $(ROOTPATH)/Manifolds/L2Sphere/L2Sphere.cpp $(ROOTPATH)/Manifolds/L2Sphere/L2SphereVariable.cpp $(ROOTPATH)/Manifolds/L2Sphere/L2SphereVector.cpp 
CPPS += $(ROOTPATH)/Manifolds/LowRank/LowRank.cpp $(ROOTPATH)/Manifolds/LowRank/LowRankVariable.cpp $(ROOTPATH)/Manifolds/LowRank/LowRankVector.cpp 
CPPS += $(ROOTPATH)/Manifolds/Oblique/Oblique.cpp $(ROOTPATH)/Manifolds/Oblique/ObliqueVariable.cpp $(ROOTPATH)/Manifolds/Oblique/ObliqueVector.cpp 
CPPS += $(ROOTPATH)/Manifolds/OrthGroup/OrthGroup.cpp $(ROOTPATH)/Manifolds/OrthGroup/OrthGroupVariable.cpp $(ROOTPATH)/Manifolds/OrthGroup/OrthGroupVector.cpp 
CPPS += $(ROOTPATH)/Manifolds/PreShapeCurves/PSCVariable.cpp $(ROOTPATH)/Manifolds/PreShapeCurves/PSCVector.cpp $(ROOTPATH)/Manifolds/PreShapeCurves/PreShapeCurves.cpp 
CPPS += $(ROOTPATH)/Manifolds/SPDManifold/SPDManifold.cpp $(ROOTPATH)/Manifolds/SPDManifold/SPDVariable.cpp $(ROOTPATH)/Manifolds/SPDManifold/SPDVector.cpp 
CPPS += $(ROOTPATH)/Manifolds/SPDTensor/SPDTVariable.cpp $(ROOTPATH)/Manifolds/SPDTensor/SPDTVector.cpp $(ROOTPATH)/Manifolds/SPDTensor/SPDTensor.cpp 
CPPS += $(ROOTPATH)/Manifolds/Sphere/Sphere.cpp $(ROOTPATH)/Manifolds/Sphere/SphereVariable.cpp $(ROOTPATH)/Manifolds/Sphere/SphereVector.cpp 
CPPS += $(ROOTPATH)/Manifolds/SphereTx/SphereTx.cpp 
CPPS += $(ROOTPATH)/Manifolds/Stiefel/StieVariable.cpp $(ROOTPATH)/Manifolds/Stiefel/StieVector.cpp $(ROOTPATH)/Manifolds/Stiefel/Stiefel.cpp 
CPPS += $(ROOTPATH)/Others/ForDebug.cpp $(ROOTPATH)/Others/MinPNormConHull.cpp $(ROOTPATH)/Others/MyMatrix.cpp $(ROOTPATH)/Others/Spline.cpp $(ROOTPATH)/Others/Timer.cpp $(ROOTPATH)/Others/randgen.cpp 
CPPS += $(ROOTPATH)/Others/SparseBLAS/nist_spblas.cpp 
CPPS += $(ROOTPATH)/Others/wavelet/wavelet.cpp 
CPPS += $(ROOTPATH)/Problems/Problem.cpp $(ROOTPATH)/Problems/juliaProblem.cpp $(ROOTPATH)/Problems/mexProblem.cpp 
CPPS += $(ROOTPATH)/Problems/CFR2BlindDecon2D/CFR2BlindDecon2D.cpp 
CPPS += $(ROOTPATH)/Problems/CFR2BlindDeconvolution/CFR2BlindDeconvolution.cpp 
CPPS += $(ROOTPATH)/Problems/CSOPhaseRetrieval/CSOPhaseRetrieval.cpp 
CPPS += $(ROOTPATH)/Problems/ElasticCurvesRO/DriverElasticCurvesRO.cpp $(ROOTPATH)/Problems/ElasticCurvesRO/ElasticCurvesRO.cpp 
CPPS += $(ROOTPATH)/Problems/EucBlindDeconvolution/EucBlindDeconvolution.cpp 
CPPS += $(ROOTPATH)/Problems/EucFrechetMean/EucFrechetMean.cpp 
CPPS += $(ROOTPATH)/Problems/EucPosSpCd/EucPosSpCd.cpp 
CPPS += $(ROOTPATH)/Problems/EucQuadratic/EucQuadratic.cpp 
CPPS += $(ROOTPATH)/Problems/GrassRQ/GrassRQ.cpp 
CPPS += $(ROOTPATH)/Problems/KarcherMean/KarcherMean.cpp 
CPPS += $(ROOTPATH)/Problems/LRBlindDeconvolution/LRBlindDeconvolution.cpp 
CPPS += $(ROOTPATH)/Problems/LRMatrixCompletion/LRMatrixCompletion.cpp 
CPPS += $(ROOTPATH)/Problems/ObliqueSparsePCA/ObliqueSparsePCA.cpp 
CPPS += $(ROOTPATH)/Problems/ObliqueTestSparsePCA/ObliqueTestSparsePCA.cpp 
CPPS += $(ROOTPATH)/Problems/OrthBoundingBox/OrthBoundingBox.cpp 
CPPS += $(ROOTPATH)/Problems/PreShapePathStraighten/PreShapePathStraighten.cpp 
CPPS += $(ROOTPATH)/Problems/SPDMean/SPDMean.cpp 
CPPS += $(ROOTPATH)/Problems/SPDTensorDL/SPDTensorDL.cpp 
CPPS += $(ROOTPATH)/Problems/ShapePathStraighten/ShapePathStraighten.cpp 
CPPS += $(ROOTPATH)/Problems/SphereConvexHull/SphereConvexHull.cpp 
CPPS += $(ROOTPATH)/Problems/SphereSparsestVector/SphereSparsestVector.cpp 
CPPS += $(ROOTPATH)/Problems/SphereTxRQ/SphereTxRQ.cpp 
CPPS += $(ROOTPATH)/Problems/StieBrockett/StieBrockett.cpp 
CPPS += $(ROOTPATH)/Problems/StieSoftICA/StieSoftICA.cpp 
CPPS += $(ROOTPATH)/Problems/StieSparseBrockett/StieSparseBrockett.cpp 
CPPS += $(ROOTPATH)/Problems/StieSumBrockett/StieSumBrockett.cpp 
CPPS += $(ROOTPATH)/Problems/WeightedLowrank/WeightedLowRank.cpp 
CPPS += $(ROOTPATH)/Solvers/LRBFGS.cpp $(ROOTPATH)/Solvers/LRBFGSLPSub.cpp $(ROOTPATH)/Solvers/LRTRSR1.cpp $(ROOTPATH)/Solvers/MRankAdaptive.cpp $(ROOTPATH)/Solvers/QuasiNewton.cpp $(ROOTPATH)/Solvers/RBFGS.cpp $(ROOTPATH)/Solvers/RBFGSLPSub.cpp $(ROOTPATH)/Solvers/RBroydenFamily.cpp $(ROOTPATH)/Solvers/RCG.cpp $(ROOTPATH)/Solvers/RGS.cpp $(ROOTPATH)/Solvers/RNewton.cpp $(ROOTPATH)/Solvers/RSD.cpp $(ROOTPATH)/Solvers/RTRNewton.cpp $(ROOTPATH)/Solvers/RTRSD.cpp $(ROOTPATH)/Solvers/RTRSR1.cpp $(ROOTPATH)/Solvers/RWRBFGS.cpp $(ROOTPATH)/Solvers/Solvers.cpp $(ROOTPATH)/Solvers/SolversLS.cpp $(ROOTPATH)/Solvers/SolversLSLPSub.cpp $(ROOTPATH)/Solvers/SolversTR.cpp 
CPPS += $(ROOTPATH)/test/ComplexToReal.cpp $(ROOTPATH)/test/RealToComplex.cpp $(ROOTPATH)/test/TestCFR2BlindDecon2D.cpp $(ROOTPATH)/test/TestCFR2BlindDeconvolution.cpp $(ROOTPATH)/test/TestCSO.cpp $(ROOTPATH)/test/TestCSOPhaseRetrieval.cpp $(ROOTPATH)/test/TestElasticCurvesRO.cpp $(ROOTPATH)/test/TestEucBlindDeconvolution.cpp $(ROOTPATH)/test/TestEucFrechetMean.cpp $(ROOTPATH)/test/TestEucPosSpCd.cpp $(ROOTPATH)/test/TestEucQuadratic.cpp $(ROOTPATH)/test/TestGrassRQ.cpp $(ROOTPATH)/test/TestKarcherMean.cpp $(ROOTPATH)/test/TestLRBlindDeconvolution.cpp $(ROOTPATH)/test/TestLRMatrixCompletion.cpp $(ROOTPATH)/test/TestMyMatrix.cpp $(ROOTPATH)/test/TestOrthBoundingBox.cpp $(ROOTPATH)/test/TestPreShapePathStraighten.cpp $(ROOTPATH)/test/TestProduct.cpp $(ROOTPATH)/test/TestProductExample.cpp $(ROOTPATH)/test/TestSPDMean.cpp $(ROOTPATH)/test/TestSPDTensorDL.cpp $(ROOTPATH)/test/TestShapePathStraighten.cpp $(ROOTPATH)/test/TestSimpleExample.cpp $(ROOTPATH)/test/TestSparsePCA.cpp $(ROOTPATH)/test/TestSphereRayQuo.cpp $(ROOTPATH)/test/TestSphereSparsestVector.cpp $(ROOTPATH)/test/TestStieBrockett.cpp $(ROOTPATH)/test/TestStieSoftICA.cpp $(ROOTPATH)/test/TestStieSparseBrockett.cpp $(ROOTPATH)/test/TestTestSparsePCA.cpp $(ROOTPATH)/test/TestWeightedLowRank.cpp 
# convert a string to upper case.
UPPER_TP  = $(shell echo $(TP) | tr a-z A-Z)

# make a binary file, which is called in command line
ROPTLIB:
	$(CC) -O3 -w -std=c++0x $(ROOTPATH)/test/$(TP).cpp $(CPPS) $(INCDIRS) -D$(UPPER_TP) -DROPTLIB_WITH_FFTW -llapack -lblas -lfftw3 -lm -o $(TP)

#make a library
libropt.so:
	$(CC) -w -std=c++0x -shared -fPIC -O3 $(CPPS) $(INCDIRS) -DROPTLIB_WITH_FFTW -llapack -lblas -lfftw3 -lm -o $@

JULIA_LIB:=$(JULIA_DIR)/usr/lib
JULIA_SRC:=$(JULIA_DIR)/src
JULIA_INC:=$(JULIA_DIR)/usr/include
CPPFLAGS:=-I$(JULIA_INC) -I$(JULIA_SRC) -I$(JULIA_SRC)/support
LDFLAGS:=-L$(JULIA_LIB)
LDLIBS=-ljulia
export LD_LIBRARY_PATH:=$(JULIA_LIB):$(JULIA_LIB)/julia

# make a shared library, which is used by Julia
JuliaROPTLIB:
	$(CC) -O3 -shared -fPIC -std=c++0x $(ROOTPATH)/test/$(TP).cpp $(CPPS) $(INCDIRS) -D$(UPPER_TP) $(CPPFLAGS) $(LDFLAGS) -Wl,-rpath,$(JULIA_LIB) -lm $(LDLIBS) -DJULIA_LIB_DIR=\"$(JULIA_DIR)/lib/julia\" -DROPTLIB_WITH_FFTW -llapack -lblas -lfftw3 -o $(TP).so

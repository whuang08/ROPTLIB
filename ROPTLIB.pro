QT += core
QT -= gui

CONFIG += c++11

TARGET = ROPTLIB
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

SOURCES += ./test/*.cpp \
    ./Solvers/*.cpp \
    ./Problems/*.cpp \
    ./Problems/ElasticCurvesRO/*.cpp \
    ./Problems/EucFrechetMean/*.cpp \
    ./Problems/EucPosSpCd/*.cpp \
    ./Problems/EucQuadratic/*.cpp \
    ./Problems/GrassRQ/*.cpp \
    ./Problems/KarcherMean/*.cpp \
    ./Problems/ObliqueTestSparsePCA/*.cpp \
    ./Problems/OrthBoundingBox/*.cpp \
    ./Problems/PreShapePathStraighten/*.cpp \
    ./Problems/ShapePathStraighten/*.cpp \
    ./Problems/SPDMean/*.cpp \
    ./Problems/SPDTensorDL/*.cpp \
    ./Problems/SphereConvexHull/*.cpp \
    ./Problems/StieBrockett/*.cpp \
    ./Problems/StieSoftICA/*.cpp \
    ./Problems/StieSparseBrockett/*.cpp \
    ./Problems/StieSumBrockett/*.cpp \
    ./Problems/WeightedLowrank/*.cpp \
    ./Others/*.cpp \
    ./Manifolds/*.cpp \
    ./Manifolds/CpxNStQOrth/*.cpp \
    ./Manifolds/ElasticShape/*.cpp \
    ./Manifolds/Euclidean/*.cpp \
    ./Manifolds/EucPositive/*.cpp \
    ./Manifolds/Grassmann/*.cpp \
    ./Manifolds/L2Sphere/*.cpp \
    ./Manifolds/LowRank/*.cpp \
    ./Manifolds/Oblique/*.cpp \
    ./Manifolds/OrthGroup/*.cpp \
    ./Manifolds/PreShapeCurves/*.cpp \
    ./Manifolds/SPDManifold/*.cpp \
    ./Manifolds/SPDTensor/*.cpp \
    ./Manifolds/Sphere/*.cpp \
    ./Manifolds/Stiefel/*.cpp \


HEADERS += ./test/*.h \
    ./Solvers/*.h \
    ./Problems/*.h \
    ./Problems/ElasticCurvesRO/*.h \
    ./Problems/EucFrechetMean/*.h \
    ./Problems/EucPosSpCd/*.h \
    ./Problems/EucQuadratic/*.h \
    ./Problems/GrassRQ/*.h \
    ./Problems/KarcherMean/*.h \
    ./Problems/ObliqueTestSparsePCA/*.h \
    ./Problems/OrthBoundingBox/*.h \
    ./Problems/PreShapePathStraighten/*.h \
    ./Problems/ShapePathStraighten/*.h \
    ./Problems/SPDMean/*.h \
    ./Problems/SPDTensorDL/*.h \
    ./Problems/SphereConvexHull/*.h \
    ./Problems/StieBrockett/*.h \
    ./Problems/StieSoftICA/*.h \
    ./Problems/StieSparseBrockett/*.h \
    ./Problems/StieSumBrockett/*.h \
    ./Problems/WeightedLowrank/*.h \
    ./Others/*.h \
    ./Manifolds/*.h \
    ./Manifolds/CpxNStQOrth/*.h \
    ./Manifolds/ElasticShape/*.h \
    ./Manifolds/Euclidean/*.h \
    ./Manifolds/EucPositive/*.h \
    ./Manifolds/Grassmann/*.h \
    ./Manifolds/L2Sphere/*.h \
    ./Manifolds/LowRank/*.h \
    ./Manifolds/Oblique/*.h \
    ./Manifolds/OrthGroup/*.h \
    ./Manifolds/PreShapeCurves/*.h \
    ./Manifolds/SPDManifold/*.h \
    ./Manifolds/SPDTensor/*.h \
    ./Manifolds/Sphere/*.h \
    ./Manifolds/Stiefel/*.h \

INCLUDEPATH += /home/whuang/julia-0.6.0_2016-08-30/include/julia/

INCLUDEPATH += ./
INCLUDEPATH += ./BinaryFiles/
INCLUDEPATH += ./Manifolds/
INCLUDEPATH += ./Manifolds/CpxNStQOrth/
INCLUDEPATH += ./Manifolds/ElasticShape/
INCLUDEPATH += ./Manifolds/EucPositive/
INCLUDEPATH += ./Manifolds/Euclidean/
INCLUDEPATH += ./Manifolds/Grassmann/
INCLUDEPATH += ./Manifolds/L2Sphere/
INCLUDEPATH += ./Manifolds/LowRank/
INCLUDEPATH += ./Manifolds/Oblique/
INCLUDEPATH += ./Manifolds/OrthGroup/
INCLUDEPATH += ./Manifolds/PreShapeCurves/
INCLUDEPATH += ./Manifolds/SPDManifold/
INCLUDEPATH += ./Manifolds/SPDTensor/
INCLUDEPATH += ./Manifolds/Sphere/
INCLUDEPATH += ./Manifolds/Stiefel/
INCLUDEPATH += ./Matlab/
INCLUDEPATH += ./Matlab/ForCpp/
INCLUDEPATH += ./Matlab/ForCpp/Boundingbox/
INCLUDEPATH += ./Matlab/ForCpp/BrockettLRBFGSVTpaper/
INCLUDEPATH += ./Matlab/ForCpp/BrockettNonconvex/
INCLUDEPATH += ./Matlab/ForCpp/RepaRotCurves/
INCLUDEPATH += ./Matlab/ForCpp/SPDTensorDLandSC/
INCLUDEPATH += ./Matlab/ForCpp/SPDTensorDLandSC/DLandSC/
INCLUDEPATH += ./Matlab/ForCpp/SPDTensorDLandSC/EucPosSC/
INCLUDEPATH += ./Matlab/ForCpp/SPDTensorDLandSC/EucPosSC/algos/
INCLUDEPATH += ./Matlab/ForCpp/SPDTensorDLandSC/EucPosSC/tools/
INCLUDEPATH += ./Matlab/ForCpp/SPDTensorDLandSC/EucPosSC/tools/others/
INCLUDEPATH += ./Matlab/ForCpp/SPDTensorDLandSC/SPDtensorDL/
INCLUDEPATH += ./Matlab/ForCpp/SoftICA/
INCLUDEPATH += ./Matlab/ForCpp/SoftICA/New_folder/
INCLUDEPATH += ./Matlab/ForCpp/SoftICA/RBFGSNonconvexResults/
INCLUDEPATH += ./Matlab/ForCpp/SparsePCA/
INCLUDEPATH += ./Matlab/ForMatlab/
INCLUDEPATH += ./Matlab/ForMatlab/FromMelissa/
INCLUDEPATH += ./Others/
INCLUDEPATH += ./Problems/
INCLUDEPATH += ./Problems/ElasticCurvesRO/
INCLUDEPATH += ./Problems/EucFrechetMean/
INCLUDEPATH += ./Problems/EucPosSpCd/
INCLUDEPATH += ./Problems/EucQuadratic/
INCLUDEPATH += ./Problems/GrassRQ/
INCLUDEPATH += ./Problems/KarcherMean/
INCLUDEPATH += ./Problems/ObliqueTestSparsePCA/
INCLUDEPATH += ./Problems/OrthBoundingBox/
INCLUDEPATH += ./Problems/PreShapePathStraighten/
INCLUDEPATH += ./Problems/SPDMean/
INCLUDEPATH += ./Problems/SPDTensorDL/
INCLUDEPATH += ./Problems/ShapePathStraighten/
INCLUDEPATH += ./Problems/SphereConvexHull/
INCLUDEPATH += ./Problems/StieBrockett/
INCLUDEPATH += ./Problems/StieSoftICA/
INCLUDEPATH += ./Problems/StieSparseBrockett/
INCLUDEPATH += ./Problems/StieSumBrockett/
INCLUDEPATH += ./Problems/WeightedLowrank/
INCLUDEPATH += ./Solvers/
INCLUDEPATH += ./test/
INCLUDEPATH += ./cwrapper/blas/
INCLUDEPATH += ./cwrapper/lapack/


% Compile all test problems in Matlab

fprintf('Compiling ComplexToReal...\n');
MyMex ComplexToReal
fprintf('Compiling DriverMexProb...\n');
MyMex DriverMexProb
fprintf('Compiling RealToComplex...\n');
MyMex RealToComplex
fprintf('Compiling TestCSO...\n');
MyMex TestCSO
fprintf('Compiling TestElasticCurvesRO...\n');
MyMex TestElasticCurvesRO
fprintf('Compiling TestEucFrechetMean...\n');
MyMex TestEucFrechetMean
fprintf('Compiling TestEucPosSpCd...\n');
MyMex TestEucPosSpCd
fprintf('Compiling TestEucQuadratic...\n');
MyMex TestEucQuadratic
fprintf('Compiling TestGrassRQ...\n');
MyMex TestGrassRQ
fprintf('Compiling TestKarcherMean...\n');
MyMex TestKarcherMean
fprintf('Compiling TestLRMatrixCompletion...\n');
MyMex TestLRMatrixCompletion
% MyMex TestMyMatrix; This can not be used in Matlab
fprintf('Compiling TestOrthBoundingBox...\n');
MyMex TestOrthBoundingBox
% MyMex TestPreShapePathStraighten; This is not done yet. Yaqing TODO
% MyMex TestShapePathStraighten; This is not done yet. Yaqing TODO
fprintf('Compiling TestProduct...\n');
MyMex TestProduct
fprintf('Compiling TestProductExample...\n');
MyMex TestProductExample
fprintf('Compiling TestSimpleExample...\n');
MyMex TestSimpleExample
fprintf('Compiling TestSparsePCA...\n');
MyMex TestSparsePCA
fprintf('Compiling TestSPDMean...\n');
MyMex TestSPDMean
fprintf('Compiling TestSPDTensorDL...\n');
MyMex TestSPDTensorDL
fprintf('Compiling TestSphereRayQuo...\n');
MyMex TestSphereRayQuo
fprintf('Compiling TestStieBrockett...\n');
MyMex TestStieBrockett
fprintf('Compiling TestStieSoftICA...\n');
MyMex TestStieSoftICA
fprintf('Compiling TestStieSparseBrockett...\n');
MyMex TestStieSparseBrockett
fprintf('Compiling TestTestSparsePCA...\n');
MyMex TestTestSparsePCA
fprintf('Compiling TestWeightedLowRank...\n');
MyMex TestWeightedLowRank


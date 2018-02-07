function MTestEucQuadratic()
    seed = floor(rand() * 100000);
    seed = 3;
    fprintf('MTestEucQuadratic seed:%d\n', seed);
    rand('state', seed);
    randn('state', seed);
    dim = 5;
    M = randn(dim, dim);
    M = M * M';
    Xinitial = randn(dim, 1);
    TestEucQuadratic(M, Xinitial);
end

function MTestEucQuadratic()
    dim = 10;
    M = randn(dim, dim);
    M = M * M';
    Xinitial = randn(dim, 1);
    TestEucQuadratic(M, Xinitial);
end

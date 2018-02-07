function MTestSphereRayQuo()
    seed = floor(rand() * 100000);
    seed = 2;
    fprintf('MTestSphereRayQuo seed:%d\n', seed);
    rand('state', seed);
    randn('state', seed);
    n = 12;
    B = randn(n, n);
    B = B + B';
    Xinitial = randn(n, 1);
    Xinitial = Xinitial / norm(Xinitial);
    TestSphereRayQuo(B, Xinitial);
end

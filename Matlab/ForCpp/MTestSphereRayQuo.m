function MTestSphereRayQuo()
    n = 12;
    B = randn(n, n);
    B = B + B';
    Xinitial = randn(n, 1);
    Xinitial = Xinitial / norm(Xinitial);
    TestSphereRayQuo(B, Xinitial);
end

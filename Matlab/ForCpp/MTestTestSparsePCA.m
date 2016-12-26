function MTestTestSparsePCA()
    p = 500; r = 2; n = 100;
    epsilons = 1e-4;
    B = randn(p, n);
    B = B - repmat((mean(B, 2)), 1, n);
    [U, S, V] = svd(B, 0);
    D = diag(S(1:r, 1:r));
    Dsq = D.^2;
    Xinitial = U(:, 1 : r);
    mus = 1e-6;
    Xopt = TestTestSparsePCA(B, Dsq, Xinitial, epsilons, mus);
end
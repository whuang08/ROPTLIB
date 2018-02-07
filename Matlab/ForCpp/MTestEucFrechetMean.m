function MTestEucFrechetMean()
    seed = floor(rand() * 100000);
    seed = 2;
    fprintf('MTestEucFrechetMean seed:%d\n', seed);
    rand('state', seed);
    randn('state', seed);
    dim = 10;
    num = 20;
    Data = randn(dim, num);
    Weight = (num : -1 : 1)';
    Xinitial = randn(dim, 1);
    TestEucFrechetMean(Data, Weight, Xinitial);
end

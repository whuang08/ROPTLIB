function MTestEucFrechetMean()
    dim = 10;
    num = 20;
    Data = randn(dim, num);
    Weight = (num : -1 : 1)';
    Xinitial = randn(dim, 1);
    TestEucFrechetMean(Data, Weight, Xinitial);
end

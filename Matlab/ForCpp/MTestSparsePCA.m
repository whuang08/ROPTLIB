function MTestSparsePCA()
    p = 500; r = 2; n = 100;
    B = randn(p, n);
    B = B - repmat((mean(B, 2)), 1, n);
    [U, S, V] = svd(B, 0);
    D = diag(S(1:r, 1:r));
    Dsq = D.^2;
    Xinitial = U(:, 1 : r);
    mus = 1e-2;
    
    SolverParams.method = 'RSD';
    SolverParams.IsCheckParams = 1;
    SolverParams.Max_Iteration = 2000;
    SolverParams.DEBUG = 2;
    SolverParams.OutputGap = 10;
    SolverParams.Accuracy = 1e50;
    SolverParams.Tolerance = 1e-6;% * norm(B);
    SolverParams.Finalstepsize = -1;
    SolverParams.Initstepsize = 1 / p / n / 10;
	SolverParams.InitSteptype = 4;
    SolverParams.IsStopped = @(x, gf, f, ngf, ngf0)IsStopped(x, gf, f, ngf, ngf0, p, r);
    
    
    Xopt = TestSparsePCA(B, Dsq, Xinitial, mus, 0, 5, SolverParams); %, soln);
end

function output = IsStopped(x, gf, f, ngf, ngf0, p, r)
    x = reshape(x.main, p, r);
    Rgf = reshape(gf.main, p, r);
    mgf = zeros(size(Rgf));
    for i = 1 : size(Rgf, 1)
        for j = 1 : size(Rgf, 2)
            if(x(i, j) == 0)
                mgf(i, j) = sign(Rgf(i, j)) * max(abs(Rgf(i, j)) - 1, 0);
            elseif(x(i, j) > 0)
                mgf(i, j) = Rgf(i, j) + 1;
            elseif(x(i, j) < 0)
                mgf(i, j) = Rgf(i, j) - 1;
            end
        end
    end

    Rmingf = mgf - x * diag(diag(x' * mgf));
%     fprintf('mingf:%.2e\n', norm(Rmingf, 'fro'));
    output = norm(Rmingf, 'fro') < 1e-9;
end

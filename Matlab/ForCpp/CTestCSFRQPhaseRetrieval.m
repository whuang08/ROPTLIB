% function [T, G] = MTestCSOPhaseRetrieval%(mask, b_exact, xtrue, xinitial)
% test file for test_PhaseRetrieval.
% 
% By Wen Huang
    r = rand();
    r = floor(r * 10000);
    r = 5028;
    fprintf('MTestCSOPhaseRetrieval seed:%d\n', r);
    rand('state', r);
    randn('state', r);
    n1 = 64;
    n2 = 64;
    n = n1 * n2;
    l = 6;
    m = l * n;
    k = 1;
    kappa = 1 / sqrt(n);
    
    %% generate masks and a true solution.
    masks = complex(randn(n1, n2, l), randn(n1, n2, l));
    phase = complex(randn(n1, n2), randn(n1, n2));
    phase = phase / norm(phase);
    b = zeros(m, 1);
    for i = 1 : l
        zi = reshape(fft2(phase .* masks(:, :, i)) / sqrt(n1 * n2), [], 1);
        b((i - 1) * n + 1 : i * n) = conj(zi) .* zi;
    end
%     bnorm = norm(b)
%     phasenorm = norm(phase, 'fro')
    
    %%
    tic
    tmp = conj(masks).* masks;
    lambda = sqrt(sum(b) / real(sum(tmp(:))) * n);
    x1.U = randn(n, k);
    [x1.U, ~] = qr(x1.U, 0);
    for i = 1 : 10
        x1.U = Ax(x1.U, masks, b, kappa);
        [x1.U, ~] = qr(x1.U, 0);
    end
    x1.U = lambda * x1.U;
    toc
    
%     x1.U = reshape(phase + complex(randn(size(phase)), randn(size(phase))) * 0.01, [], 1);
    
    initX = x1.U;
    masksn = masks;
    
    
%     SolverParams.method = 'RBFGS';
%     SolverParams.method = 'RBroydenFamily';
    SolverParams.method = 'RTRNewton';
%     SolverParams.method = 'RCG';
%     SolverParams.method = 'LRBFGS';
    SolverParams.IsCheckParams = 1;
    SolverParams.Max_Iteration = 200;
    SolverParams.LengthSY = 4;
    SolverParams.Verbose = 2;
    SolverParams.LMrestart = 0;
    SolverParams.Accuracy = 1e-6;
    SolverParams.Tolerance = 1e-9;
    SolverParams.OutputGap = 1;
    SolverParams.Finalstepsize = 1;
    SolverParams.Num_pre_funs = 1;
    SolverParams.PreFunsAccuracy = 1e6;
    HasHHR = 0;
    ParamSet = 1;
    [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times] = TestCSFRQPhaseRetrieval(b, masksn, initX, kappa, ParamSet, HasHHR, SolverParams);

% end

function output = Ax(U, masks, b, kappa)
    x.U = U;
    [n1, n2, l] = size(masks);
    n = n1 * n2;
    m = l * n;
    k = size(x.U, 2);
    ZY = zeros(m, k);
    for i = 1 : l
        for j = 1 : k
            temp = fft2(reshape(x.U(:, j), n1, n2) .* masks(:, :, i)) / sqrt(n1 * n2);
            ZY((i - 1) * n + 1 : i * n, j) = reshape(temp, [], 1);
        end
    end
    DZY = zeros(m, k);
    for i = 1 : k
        DZY(:, i) = b .* ZY(:, i);
    end
    %% Z^H *DZY
    Y = zeros(n, k);
    for i = 1 : l
        temp = zeros(n, k);
        for j = 1 : k
            temp(:, j) = reshape(ifft2(reshape(DZY((i - 1) * n + 1 : i * n, j) * sqrt(n1 * n2), n1, n2)), [], 1);
            temp(:, j) = conj(reshape(masks(:, :, i), [], 1)) .* temp(:, j);
        end
        Y = Y + temp;
    end
    output = Y;
end

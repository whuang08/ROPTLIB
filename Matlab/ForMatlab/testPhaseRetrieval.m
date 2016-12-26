function test()
    testPhaseRetrieval1(4);
%     testPhaseRetrieval1(2);
end

function [FinalX, fv, gfv, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times] = testPhaseRetrieval1(k)

    r = rand();
    r = floor(r * 10000)
%     r = 1
% r = 3665
    rand('state', r);
    randn('state', r);
    n1 = 256;
    n2 = 256;
    global count
    count = 0;
    n = n1 * n2;
    l = 6;
    m = l * n;
%     k = 2;
    p = k;
    global iscomplex
    iscomplex = 1;
    masks = complex(randn(n1, n2, l), iscomplex * randn(n1, n2, l));
    phase = complex(randn(n1, n2), iscomplex * randn(n1, n2));
    phase = phase / norm(phase);
    b = zeros(m, 1);
    for i = 1 : l
        zi = reshape(fft2(phase .* masks(:, :, i)) / sqrt(n1 * n2), [], 1);
        b((i - 1) * n + 1 : i * n) = conj(zi) .* zi;
    end
    
    isplot = 1;
    delta = 0.95;
    kappa = 1/n;
    
    inittype = 2;
    %% initial x type 1
    if(inittype == 1)
        [x1.U, ~] = qr(complex(randn(n, k), iscomplex * randn(n, k)), 0);
        x1.D = spdiags(ones(k, 1), 0, k, k);
        initialX.main = full(x1.U * sqrt(x1.D));
        initialX.main = initialX.main * GetScalar(initialX, masks, b, 0);
    end

    %% initial x type 2
    if(inittype == 2)
        Axhandle = @(U)Ax(U, masks, b, kappa);
        opts.isreal = false;
        opts.maxit = 10;
        [V, D] = eigs(Axhandle, n, k, 'lm', opts);
        initialX.main = V;
    end
    %% end initial x

    fhandle = @(x)f(x, masks, b, kappa);
    gfhandle = @(x)gf(x, masks, b, kappa);
    Hesshandle = @(x, eta)Hess(x, eta, masks, b, kappa);

%     initialf = fhandle(initialX)%%-----------
    
    SolverParams.method = 'LRBFGS';
%     SolverParams.method = 'RCG';
    SolverParams.IsCheckParams = 1;
    SolverParams.DEBUG = 2;
    SolverParams.Min_Iteration = 10;
    SolverParams.Max_Iteration = 500;
    SolverParams.OutputGap = 1;
    SolverParams.IsCheckGradHess = 0;
    SolverParams.IsStopped = @(x, gf, f, ngf, ngf0)IsStopped(x, gf, f, ngf, ngf0, delta);
%     SolverParams.LineSearch_LS = 4;
    SolverParams.LengthSY = 2;
%     SolverParams.LinesearchInput = @(x, eta, t0, s0)LinesearchInput(x, eta, t0, s0, masks, b, kappa);

    ManiParams.IsCheckParams = 1;
    ManiParams.name = 'CpxNStQOrth';
    ManiParams.n = n;
    ManiParams.p = p;
    HasHHR = 0;

    X = {}; F = []; G = []; T = []; timecost = 0; nf = 0; ng = 0; nR = 0; nH = 0; nV = 0; nVp = 0;    
    r = size(initialX.main, 2);
    while(r >= 1)
        fprintf('rank is %d\n', r);
        ManiParams.p = r;
        if(r == 1)
            SolverParams.Max_Iteration = 500;
        end
        [FinalX, fv, gfv, gfgf0, iter, nfi, ngi, nRi, nVi, nVpi, nHi, timecosti, Fi, Gi, Ti] = DriverOPT(fhandle, gfhandle, Hesshandle, SolverParams, ManiParams, HasHHR, initialX);
        % combine output information
        F = [F; Fi]; G = [G; Gi]; timecost = timecost + timecosti;
        if(length(T) > 0)
            T = [T; Ti + T(end)];
        else
            T = Ti;
        end
        nf = nf + nfi; ng = ng + ngi; nR = nR + nRi; nH = nH + nHi; nV = nV + nVi; nVp = nVp + nVpi;
        xf = FinalX.main;
        r = size(xf, 2);
        if(r == 1) % obtain solution for phaselift
            break;
        end
        
        [U, S, ~] = svd(xf, 0);
        normSdsqrk = norm(diag(S)) / sqrt(size(S, 1));
        x0 = U(:, 1 : size(U, 2) - 1) * S(1 : size(U, 2) - 1, 1 : size(U, 2) - 1);
        for i = 1 : size(U, 2)
            if(S(i, i) / normSdsqrk < delta)
                x0 = U(:, 1 : i - 1) * S(1 : i - 1, 1 : i - 1);
                break;
            end
        end

        % get next initial point
        initialX.main = x0;
%         x0
        r = size(x0, 2);
    end
    fprintf('Num. of total Iter.: %d, time cost: %f seconds, Num. of Fun, Gra, Hv, VT and R: %d, %d, %d, %d, %d \n', length(F), timecost, nf, ng, nH, nV, nR)
    if(isplot == 1)
        figure;clf
        semilogy(1 : length(F), F, '.');
        ylabel('f(x_i)');
        xlabel('iteration number');
        title(num2str(timecost));
%         figure;clf
%         semilogy(1 : length(G), G, '.');
%         ylabel('|grad(x_i)|');
%         xlabel('iteration number');
    end
end

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

function output = LinesearchInput(x, eta, t0, s0, masks, b, kappa)
% x
% eta
% t0
% s0
% tic
    [n1, n2, l] = size(masks);
    n = n1 * n2;
    k = size(x.main, 2);
    m = l * n;
    D = zeros(m, 1);
    Zeta = zeros(m, k);
    ZY = x.ZY;
    for i = 1 : l
        for j = 1 : k
            temp = reshape(fft2(reshape(eta.main(:, j), n1, n2) / sqrt(n1 * n2) .* masks(:, :, i)), [], 1);
            Zeta((i - 1) * n + 1 : i * n, j) = reshape(temp, [], 1);
        end
    end
    
    vc = -2 * (sum(real(ZY) .* real(Zeta), 2) + sum(imag(ZY) .* imag(Zeta), 2));
%     vc = zeros(m, 1);
%     for i = 1 : m
%         vc(i) = -2 * real(ZY(i, :) * Zeta(i, :)');
%     end
    output = vc' * x.cD / (vc' * vc);
end

function output = IsStopped(x, gf, f, ngf, ngf0, delta)
    M = x.main;
    r = size(M, 2);
    if(r > 1)
        [U, S, ~] = svd(M, 0);
        normSdsqrk = norm(diag(S)) / sqrt(size(S, 1));
%         [diag(S)', delta * normSdsqrk]%%----
        for i = 1 : size(U, 2)
            if(S(i, i) / normSdsqrk < delta)
                output = 1;
                return;
            end
        end

        output = 0;
    else
        output = ngf < 1e-6;
%         output = f < 1e-10;
    end
end



function output = GetScalar(x, masks, b, kappa)
    [n1, n2, l] = size(masks);
    n = n1 * n2;
    k = size(x.main, 2);
    m = l * n;
    D = zeros(m, 1);
    ZY = zeros(m, k);
    for i = 1 : l
        sZqi = zeros(n, 1);
        for j = 1 : k
            temp = reshape(fft2(reshape(x.main(:, j), n1, n2) / sqrt(n1 * n2) .* masks(:, :, i)), [], 1);
            ZY((i - 1) * n + 1 : i * n, j) = reshape(temp, [], 1);
            sZqi = sZqi + conj(temp) .* temp;
        end
        D((i - 1) * n + 1 : i * n) = reshape(sZqi, [], 1);
    end
    output = sqrt(norm(b, 'fro') / norm(D));
end


function [output, x] = f(x, masks, b, kappa)
% tic
    [n1, n2, l] = size(masks);
    n = n1 * n2;
    sqrtn1n2 = sqrt(n1 * n2);
    k = size(x.main, 2);
    m = l * n;
    D = zeros(m, 1);
    ZY = zeros(m, k);
    for i = 1 : l
        sZqi = zeros(n, 1);
        for j = 1 : k
            temp = reshape(fft2(reshape(x.main(:, j), n1, n2) / sqrtn1n2 .* masks(:, :, i)), [], 1);
            ZY((i - 1) * n + 1 : i * n, j) = reshape(temp, [], 1);
            sZqi = sZqi + conj(temp) .* temp;
        end
        D((i - 1) * n + 1 : i * n) = reshape(sZqi, [], 1) - b((i - 1) * n + 1 : i * n);
    end
    output = norm(D)^2 / norm(b, 'fro')^2;
    x.ZY = ZY;
    x.cD = D;
    if(k ~= 1)
        output = output + kappa * norm(x.main, 'fro')^2;%% add penelty term
    end
% ft = toc
%     xf = x.main%%----
%     output%%-----
end

function [output, x] = gf(x, masks, b, kappa)
% 'h1'
% xgf = x.main%%----
% xmain = x.main%%----
% tic
    [n1, n2, l] = size(masks);
    n = n1 * n2;
    m = l * n;
    k = size(x.main, 2);
    if(~isfield(x, 'ZY') || ~isfield(x, 'cD'))
        D = zeros(m, 1);
        ZY = zeros(m, k);
        for i = 1 : l
            sZqi = 0;
            for j = 1 : k
                temp = fft2(reshape(x.main(:, j), n1, n2) .* masks(:, :, i)) / sqrt(n1 * n2);
                ZY((i - 1) * n + 1 : i * n, j) = reshape(temp, [], 1);
                sZqi = sZqi + conj(temp) .* temp;
            end
            D((i - 1) * n + 1 : i * n) = reshape(sZqi, [], 1) - b((i - 1) * n + 1 : i * n);
        end
        x.cD = D;
        x.ZY = ZY;
    end
    D = x.cD;
    ZY = x.ZY;
    DZY = zeros(m, k);
    for i = 1 : k
        DZY(:, i) = D.*ZY(:, i);
    end
    
%     DZY = spdiags(D, 0, m, m) * ZY;
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
    YY.TV = 4 * Y / norm(b, 'fro')^2;
    if(k ~= 1)
        YY.TV = YY.TV + kappa * 2 * x.main;%% for penelty term
    end
    output.main = YY.TV / (x.main' * x.main);
% gft = toc
    
%         xx = x.main' * x.main;
%         xe = x.main' * YY.TV;
%         xe = xe - xe'
    
%     omain = output.main%%---
%     if(~isfield(x, 'Omega'))
%         xx = x.main' * x.main;
%         xe = x.main' * YY.TV;
%         xe = xe - xe';
% %         xx%%----
% %         xe%%---
%         x.Omega = lyap(xx, xx, -xe);
%     end
%     size(YY.TV)%%---
%     size(x.main)%%---
%     size(x.Omega)%%---
%     output.main = (YY.TV - x.main * x.Omega);
%     output.main%%---
end

function [output, x] = Hess(x, eta, masks, b, kappa)
    [n1, n2, l] = size(masks);
    k = size(x.main, 2);
    n = n1 * n2;
    m = n * l;
    Zv = zeros(m, k);
    for i = 1 : l
        for j = 1 : k
            Zv((i - 1) * n + 1 : i * n, j) = reshape(fft2(reshape(eta.main(:, j), n1, n2) .* masks(:, :, i)) / sqrt(n1 * n2), [], 1);
        end
    end
    D = x.cD;
    ZY = x.ZY;
    DD = ZY .* conj(Zv);
    if(k ~= 1)
        DD = sum(DD')';
    end
    DD = 2 * real(DD);
    DZv = full(spdiags(D, 0, m, m) * Zv);
    DDZY = full(spdiags(DD, 0, m, m) * ZY);
    DZALL = DZv + DDZY;
    
    %% Z^H *DZALL
    dotzeta = zeros(n, k);
    for i = 1 : l
        temp = zeros(n, k);
        for j = 1 : k
            temp(:, j) = reshape(ifft2(reshape(DZALL((i - 1) * n + 1 : i * n, j) * sqrt(n1 * n2), n1, n2)), [], 1);
        end
        temp = spdiags(reshape(masks(:, :, i), [], 1), 0, n, n)' * temp;
        dotzeta = dotzeta + temp;
    end
    dotzeta = dotzeta / norm(b, 'fro')^2;
    if(k ~= 1)
        dotzeta = 4 * dotzeta + kappa * 2 * eta.main;
    end
    output.main = dotzeta;
%     Vdotzeta.TV = (dotzeta - eta.main * x.Omega);
%     output = proj_q(x, Vdotzeta);
end

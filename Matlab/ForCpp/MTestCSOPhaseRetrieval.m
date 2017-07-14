function [T, G] = test_PhaseRetrieval%(mask, b_exact, xtrue, xinitial)
% test file for test_PhaseRetrieval.
% 
% By Wen Huang
    fprintf('test for Phase Retrieval problem. \n');
    r = rand();
    r = floor(r * 10000)
    r = 5028
    rand('state', r);
    randn('state', r);
    n1 = 128;
    n2 = 128;
    global count
    count = 0;
    n = n1 * n2;
    l = 20;
    m = l * n;
    k = 2;
    kappa = 1 / sqrt(n);
    
    fprintf('In this case, n1 = %d, n2 = %d, n = %d, l = %d. b, and x0 are chosen randomly. Random seed is %d \n', n1, n2, n, l, r);
    
    
    %% generate masks and a true solution.
    masks = complex(randn(n1, n2, l), randn(n1, n2, l));
    phase = complex(randn(n1, n2), randn(n1, n2));
    phase = phase / norm(phase);
    b = zeros(m, 1);
    for i = 1 : l
        zi = reshape(fft2(phase .* masks(:, :, i)) / sqrt(n1 * n2), [], 1);
        b((i - 1) * n + 1 : i * n) = conj(zi) .* zi;
    end
    bnorm = norm(b)
    phasenorm = norm(phase, 'fro')
    
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
    
    initX = zeros(2 * n, k);
    initX(1:2:end, :) = real(x1.U);
    initX(2:2:end, :) = imag(x1.U);
    masksn = zeros(2 * n1, n2, l);
    masksn(1:2:end, :, :) = real(masks);
    masksn(2:2:end, :, :) = imag(masks);
    
    delta = 0.95;
    times = 0;
    total_timecost = 0;
    while(total_timecost < 60 && times < 1)
        fprintf('%d: ', times);
        [Xopt, timecost, nfft, nn] = TestCSOPhaseRetrieval(b, masksn, initX, kappa, delta, 1e-10, 1, 1000);
        [nfft, nn]
        total_timecost = total_timecost + timecost;
        times = times + 1;
    end
    Xsoln = complex(Xopt(1:2:end, :), Xopt(2:2:end, :));
    vphase = phase(:);
    tmp = Xsoln * (Xsoln' * vphase) / norm(Xsoln' * vphase);
    norm(vphase - Xsoln * (Xsoln' * vphase) / norm(Xsoln' * vphase)) / norm(phase)
    
    fprintf('average time cost of %d times is %f. \n', times, total_timecost / times);
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
% 
% function Ax(x::Array{Complex{Float64}}, k::Int64, masks::Array{Complex{Float64}}, n::Int64, l::Int64, m::Int64)::Array{Complex{Float64}}
% 	sqrtn::Float64 = sqrt(n)
% 	ZY::Array{Complex{Float64}, 2} = zeros(m, k)
% 	temp::Array{Complex{Float64}} = zeros(n, 1)
% 	for i in 1 : l
% 		for j in 1 : k
% 			temp[:] = reshape(fft(reshape(x[:, j], n1, n2) / sqrtn .* masks[:, :, i]), n, 1)
% 			ZY[(i - 1) * n + 1 : i * n, j] = temp
% 		end
% 	end
% 	DZY::Array{Complex{Float64}, 2} = zeros(m, k)
% 	for i in 1 : k
% 		DZY[:, i] = b .* ZY[:, i]
% 	end
% 
% 	Y::Array{Complex{Float64}, 2} = zeros(n, k)
% 	tmp::Array{Complex{Float64}} = zeros(n, k)
% 	for i in 1 : l
% 		temp *= 0
% 		for j in 1 : k
% 			tmp[:, j] = reshape(ifft(reshape(DZY[(i - 1) * n + 1 : i * n, j] * sqrtn, n1, n2)), n, 1)
% 			tmp[:, j] = conj(reshape(masks[:, :, i], n, 1)) .* tmp[:, j]
% 		end
% 		Y += tmp
% 	end
% 	return Y
% end
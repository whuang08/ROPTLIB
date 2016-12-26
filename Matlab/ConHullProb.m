function [w,f] = ConHullProb(H)
% Solver for the quadractic convex program
%
%     min          w' H w
%     w 
%
%     subject to   e'w = 1          Hw is a convex combination of gradients
%                   w >= 0          in the columns of S
%  Return w and f = w'Hw.

H = H + H';
N = size(H, 1);
e = ones(1,N);
O = zeros(1,N);
w0 = ones(N,1)/N;
qp_options = optimset('Display','off','TolX', 1e-12, 'TolFun', 1e-12);%, 'Algorithm', 'active-set');
% qp_options = optimset('Display','off','TolX', 1e-12, 'TolFun', 1e-12, 'LargeScale', 'off');
% tic
[w, hf] = quadprog(H, O, [], [], e, 1, O, [], w0, qp_options);
% toc
% tic
% w = NewtonKKTqp(H,0,e,1,w0); hf = 0.5 * w' * H * w;
% toc
% error('sdf');

% full call: [w, fval, qpflag, qpoutput, lambda] = quadprog(H, O, [], [], e, 1, O, [], w0, qp_options);
if isempty(w)
    error('w is empty: MOSEK license problem?')
end
f = hf;

prtlevel = 1;
checktol = 1e-7*max(1, norm(H, inf));
if prtlevel > 0 & (min(w) < -checktol | abs(sum(w) - 1) > checktol)
    fprintf('computed w does not satisfy requirements\n')
end

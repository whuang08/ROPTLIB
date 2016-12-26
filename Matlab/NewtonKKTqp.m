function [x_final,logs] = NewtonKKTqp(H,c,A,b,x0,opts)

% NewtonKKTqp  Find a local minimizer of an indefinite quadratic
%          programming problem.
%
%   X = NewtonKKTqp(H,C,A,B,X0) returns a local minimizer X of the (indefinite)
%   quadratic function f(X) = .5*X'*H*X+C'*X subject to the linear
%   inequality constraints A*X <= B.
%   The iteration is initialized with the feasible point X0.
%
%   [X,LOGS] = NewtonKKTqp(H,C,A,B,X0) returns a structure LOGS
%   containing various information collected during runtime.
%
%   NewtonKKTqp(H,C,A,B,X0,OPTS) specifies options:
%   opts.tol: tolerance on the KKT residual [scalar | {Inf}].
%   opts.tol_cost_fn: tolerance on the decrease of the cost function
%        [scalar | {1e-6}]. Set to Inf to deactivate.
%   opts.tol_zeta: tolerance on final multipliers [scalar |
%        {Inf}]. Set to -Inf to deactivate. Suggestion: -1e-9.
%   opts.barrier: 
%        0: Algorithm A1 (affine scaling), 
%        not 0: Algorithm A2 (barrier)
%        [scalar | {1}]
%   opts.step1: 
%        0: simplified Step1, 
%        not 0: recommended Step1
%        [scalar | {1}]
%   opts.solve_condensed: 
%        0: solve unreduced Newton-KKT system
%        not 0: solve condensed system
%        [scalar | {1}]
%   opts.beta, opts.z_under, opts.z_u, opts.sigma, opts.gamma,
%   opts.z0, opts.theta, opts.phibar, ops.nu, ops.psi are the
%        parameters defined in the paper. If not specified by 
%        the user, they are assigned a default value (recommended).
%   opts.eps_constr: allowance for receding constraints 
%        [scalar | {1e-14}]
%   opts.maxit: maximum number of iterations [scalar | {1e4}]
%   opts.verbose: 
%        not 1: minimal output 
%        1: verbose output
%
%
% This file is a "cleaned" but equivalent version of indefQP34.m.
% The latter script was used to generate the numerical results reported in the
% related paper http://www.montefiore.ulg.ac.be/~absil/Publi/indefQP.htm
%
% Reference: http://www.montefiore.ulg.ac.be/~absil/Publi/indefQP.htm
% 
% Copyright P.-A. Absil and A. L. Tits, 2004-2005.
% 
% Note: This script was built for testing purposes. It is not written
% to be time efficient.
%
% Comments and suggestions are welcome.
% Authors' URLs: 
% http://www.csit.fsu.edu/~absil/ 
% http://www.ece.umd.edu/~andre/

% indQP34_fn04.m - PA - Fri 17 Jun 2005, CSIT.
%    Allow for stopping criterion based on the distance between
%    iterates.
% indQP34_fn03.m - PA - Sun 12 Jun 2005, CSIT.
%    Cleanup. Allows to return logs.
% indQP34_fn02.m - PA - Sat 11 Jun 2005, CSIT.
%    In the update for z, use zeta_minus to bound below instead of dz.

CPU_tic = cputime;

n = length(c);
m = length(b);

%  *********************************************************
%  *********************************************************
%  Process inputs and do error checking
%  *********************************************************
%  *********************************************************

if nargin < 6
  opts = [];
end

tol_stop_crit = Inf;   % Default: 1e-7, as in JHU's slides.
if isfield(opts,'tol')
  tol_stop_crit = opts.tol;
end

tol_cost_fn = 1e-6;
if isfield(opts,'tol_cost_fn')
  tol_cost_fn = opts.tol_cost_fn;
end

tol_zeta = -Inf;  % Sugg: -1e-9
if isfield(opts,'tol_zeta')
  tol_zeta = opts.tol_zeta;
end

sw_barrier = 1;
if isfield(opts,'barrier')
  sw_barrier = opts.barrier;
end

sw_step1 = 1;
if isfield(opts,'step1')
  if opts.step1==0
    sw_step1 = -2;
  end
end
% -2: Simpler and safer: compute Hmod in each iteration.
% -1: Simplest. Avoid. Take -2.
% 0: classical case
% 1: vectorial alphabar

sw_Newton_Schur = 1;
if isfield(opts,'solve_condensed')
  sw_Newton_Schur = opts.solve_condensed;
end

beta = .9;  % req: 0<beta<1  
if isfield(opts,'beta')
  beta = opts.beta;
end

zunder = 1e-4;
if isfield(opts,'z_under')
  zunder = opts.z_under;
end

zbar = 1e15;   % req: zbar>0.  Default: zbar=Inf
if isfield(opts,'zbar')
  zbar = opts.z_u;
end

sigma = 1e-5;   % req: sigma>0. Default: 1e-5 
if isfield(opts,'sigma')
  sigma = opts.sigma;
end

gamma = 1e3;   % default: 1e3
if isfield(opts,'gamma')
  gamma = opts.gamma;
end

% x0 alredy defined

% z0:
gradf0 = H*x0+c;
if isfield(opts,'z0')
  z0 = opts.z0;
else
  z0 = max(.1*ones(m,1), -pinv(A')*gradf0);
end
% Possible choice: z0 = ones(length(b),1);

theta = .8;
if isfield(opts,'theta')
  theta = opts.theta;
end

phibar = 1e6;
if isfield(opts,'phibar')
  phibar = opts.phibar;
end

nu = 3;
if isfield(opts,'nu')
  nu = opts.nu;
end

psi = 1.5;   % 1 minimizes along line; must be <2. Psi in the paper.
if isfield(opts,'psi')
  psi = opts.psi;
end

eps_constr = 1e-14;  % used in receding constraints heuristic.
if isfield(opts,'eps_constr')
  eps_constr = opts.eps_constr;
end

k_upper = 1e4;   % maximal number of iterates. Default = 1e4.
if isfield(opts,'maxit')
  k_upper = opts.maxit;
end

do_display = 0;
if isfield(opts,'verbose')
  if opts.verbose == 1;
    do_display = 1;
  end
end

sw_z_clip_below = 1;  
% Default: 1
% 0: version script before 27 Nov 2004
% 1: version script after 27 Nov 2004 (zeta_minus)
% Version "fn" of the code updated on 11 Jun 2005, in indQP34_fn02


iter_total = 0; update_E_total = 0; eig_this = 0; solve_this = 0;

%  *********************************************************
%  *********************************************************
%  Main part
%  *********************************************************
%  *********************************************************


% ******* Step 0 - Initialization *******

k = 0;
x = x0;
z = z0;
f = NaN; % xxx
% f = .5*x'*H*x + c'*x;
g0 = A*x0-b;
g = g0;
Ibar = [];
Inonbar = [1:m];
alphabar = 0;
[VH,DH] = eig(H);  eig_this=eig_this+1;
if all(diag(DH)>sigma)  
  Ebar = zeros(n,n);
  pd_case = 1;
  E = Ebar;
  W = H + E;
  mineigHmod = min(diag(DH));
  min_rat_bar = NaN;
  rat = NaN;
else
  Ebar = eye(n,n);
  pd_case = 0;
end
f_cost = .5*x'*H*x + c'*x;
test_conv = 0;   % assume stopping criterion not satisfied

while ~test_conv & k<=k_upper
  
  if do_display
    k
  end
  
  % **** Compute f, g, Schur ****
  g = A*x-b;
  g = min(g,-eps_constr);  

  S = H-A'*diag(z./g)*A;   % !!! Screws up detailed CPU counts.
  
  % ******* Step 1 - Obtain W *******
  CPU_E_tic = cputime;
  
  if ~pd_case  % If pd_case, don't do step 1.
    switch sw_step1
     case -2  % Basic and "safe".
      Hmod = S; % H-A'*diag(z./g)*A;   % !!!
      Hmodsym = (Hmod + Hmod')/2;
      [V,D] = eig(Hmodsym); eig_this = eig_this + 1;
      [dmin,indmin] = min(diag(D));
      vmin = V(:,indmin);
      mineigHmod = vmin'*Hmodsym*vmin; 
      %mineigHmod = min(real(eig(Hmod))); eig_this = eig_this + 1;
      enter_mod = 1;   % for compatibility
      rat = z./abs(g);  % idem
      min_rat_bar = min(rat); % idem
     case 1  % New case: vectorial alphabar
      warning off  % in the next 2 lines we may divide by smth very small
      rat = z./abs(g);
      rat_bar = z(Ibar)./abs(g(Ibar));
      if isempty(rat_bar) 
	rat_bar = Inf;
      end  % convention for empty case
      min_rat_bar = min(rat_bar);
      %min_rat_nonbar = min(z(Inonbar)./abs(g(Inonbar)));
      warning on
      %if isempty(min_rat_nonbar)
      %  min_rat_nonbar = -Inf;
      %end
      enter_mod = (any(rat_bar<=alphabar) | (any(rat_bar>=gamma^2*alphabar) ...
					     & norm(Ebar) ~= 0));
      if enter_mod
	warning off
	Ibar = find(z./abs(g) >= 1); Ibar=Ibar';
	%Inonbar = find(z./abs(g) < 1); Inonbar = Inonbar';
	warning on
	if isempty(Ibar)
	  alphabar = Inf;
	  Hmod = H;    % this is because H+[] is not valid
	else
	  alphabar = gamma^(-1) * (z(Ibar)./abs(g(Ibar)));
	    Hmod = H + transpose(A(Ibar,:))*diag(alphabar)*A(Ibar,:);
	end
	mineigHmod = min(real(eig(Hmod)));  eig_this=eig_this+1;
	update_E_total = update_E_total + 1;
      end
    end  % switch sw_step1
    
    if enter_mod  
	if mineigHmod > sigma
	  Ebar = zeros(n,n);
	elseif abs(mineigHmod) <= sigma
	  Ebar = (sigma-mineigHmod)*eye(n,n);
	else
	  Ebar = 2*abs(mineigHmod)*eye(n,n);
	end
    end
    if ~exist('mineigHmod')
      mineigHmod = NaN;  % for compatibility
    end
    
    E = Ebar;
    W = H + E;
  end  % if pd_case
  

  % ******* Step 2 - Newton system *******
  
  grad = H*x+c;
  if sw_Newton_Schur
    warning off
      %S = W-A'*diag(z./g)*A;  % !!!
      S = S+E;
    d = -S\(H*x+c);  solve_this=solve_this+1;
    warning on
      zeta = -diag(z./g)*A*d;
    Delta_z = zeta-z;
  else  
    M = [W,A';diag(z)*A,diag(g)];
    solu = M\[-H*x-c;zeros(m,1)];   solve_this=solve_this+1;
    d = solu(1:n,:);
    zeta = solu(n+1:n+m,:);
    Delta_z = zeta-z;
    S = 0;  % for compatibility
  end
  
  % Save variables:
  d0 = d;
  zeta0 = zeta;
  Delta_z0 = Delta_z;
  
  if d==0
    error('d = 0')
  end
  
  if sw_barrier
    delta = grad'*d;
    % new simple mu (29 Nov 2003) % DEFAULT
             % mu still scalar but multiplied by zmin
	     % It only matters when phi saturates, since when phi
             % does not saturate then mu is the "best" one.
      zmin = min(z);
      test_mu = phibar*sum((norm(d))^nu*zeta./z*zmin)<=(1-theta)*abs(delta);
      if test_mu
        phi = phibar;
	  %disp('phibar reached')
      else
        phi = (1-theta)*abs(delta)/sum((norm(d))^nu*zeta./z*zmin);
      end
      mu = phi*(norm(d))^nu*zmin*ones(m,1);
    
    if sw_Newton_Schur
      warning off
	d_mu = S\(-H*x-c+A'*(mu./g)); solve_this= ...
	       solve_this+1;
      warning on
	zeta_mu = diag(g.^(-1)) * (-mu-diag(z)*A*d_mu);
      Delta_z_mu = zeta_mu-z;
    else   
      solu_mu = M\[-H*x-c;-mu];  solve_this=solve_this+1;
      d_mu = solu_mu(1:n,:);
      zeta_mu = solu_mu(n+1:n+m,:);
      Delta_z_mu = zeta_mu-z;
    end
    d = d_mu;  % activate new values
    
    if d'*(H*x+c)>0
      disp('Not a descent direction!'); d'*(H*x+c)
    end
    zeta = zeta_mu;
    Delta_z = Delta_z_mu;
  end
  
  
  % ******* Step 3 - Updates *******

  % ** (i) **
  Ad = A*d;
  indi = find(Ad>0);
  
  if isempty(indi)
    tbar = Inf;
  else
    tbar = min(-g(indi)./Ad(indi));
  end

  dtHd = d'*H*d;
  if dtHd>0;
    t_opt = abs(d'*(H*x+c))/(d'*H*d);
  else
    t_opt = inf;
  end
  t_clip_above = min(psi*t_opt, 1);

  t = min( max(beta*tbar,tbar-norm(d)) , t_clip_above );

  xnew = x + t*d;

  % ** (ii) ** 
    switch sw_z_clip_below
     case 0
      z_clip_below = d'*d + Delta_z'*Delta_z;
     case 1
      zeta_minus = zeta(zeta<0); 
      z_clip_below = d'*d + zeta_minus'*zeta_minus;
    end
    z_clip_below = min(z_clip_below,zunder);
    znew = min(max(z_clip_below*ones(m,1), zeta ),zbar*ones(m,1));
      
    gnew = A*xnew-b;
 
  
  % **** Stopping criterion and other tests ****
  
    gnew = A*xnew-b;
    f_cost_new = .5*xnew'*H*xnew + c'*xnew;
    gap = (f_cost-f_cost_new)/(1+abs(f_cost_new));
    err_KKT_new_v = [H*xnew+c+A'*zeta; zeta.*gnew];
    resid_vec = err_KKT_new_v;
    % ! zeta instead of znew above.
  resid = norm(resid_vec,inf);
  %test_conv = (resid<=tol_stop_crit | norm(t_stepback*d,inf) <= 1e-14) & all(zeta>-1e-9);
  %test_conv = (resid<=tol_stop_crit | gap<=tol_cost_fn) & all(zeta>-1e-9); 
  %test_conv = (resid<=tol_stop_crit | gap<=tol_cost_fn);
  test_conv = (resid<=tol_stop_crit) & all(zeta>tol_zeta) & (gap<=tol_cost_fn);  
  
  if ~isfinite(resid) & flast>-1e20
    error('Bad NaN')
  end
  if test_conv & any(zeta<-1e-6)
    disp(['Computed a non-KKT stationary point: min(zeta) =', num2str(min(zeta))])
  end
  
  if ~mod(k+1,500)
    disp(k)
    %pause
  end
  
  if any(g==0) 
    disp('A comp of g is zero');
  end
  if any(isnan(g))
    disp('There is a NaN')
  end

  
  % **** Activate the updates:
  x = xnew;
  z = znew;
  g = gnew;
  f_cost_last = f_cost;
  f_cost = f_cost_new;
  k = k + 1;
  iter_total = iter_total + 1;

  if k>2e3
    disp('large k!!!')
    pause
  end
  
end  % end of loop on k
% ************************

  if k>k_upper
    disp('Too many iterates (k>k_upper).')
  end

k_final = k;
if k>=k_upper
  CPU_this = inf;
  eig_this = inf;
  solve_this = inf;
  k_final = inf;
  iter_total = inf;
end

infeas_final = max(0,max(gnew));

CPU_total = cputime - CPU_tic;

%  *********************************************************
%  *********************************************************
%  Create output
%  *********************************************************
%  *********************************************************

x_final = x;

if nargout >= 2
  logs.nb_iter = iter_total;
  logs.nb_solves = solve_this;
  logs.nb_eigsolves = eig_this;
  logs.cpu_total = CPU_total;
end

function [Sigma,err]=computeSigma_RelaxRegOT_qp2(X,Y,Gx,kappa,lambda,T,options)

% Sigma = computeSigma_RelaxRegOT_qp(X,Y,Gx,lambda,koptions);
%
%Compute the barycenter Y of Xs sets with parameters vrho
% dims: size(X) = [N,d], size(Gs) = [P,N], 
% N: number of points,
% d: number of dimensions of each point
%
% Solve:
%       min_Sigma <C,Sigma> + lambda |Gx*(diag(Sigma*1)X-Sigma*Y)|^2
% subject to
%    0 <= Sigma'*1<=kappa 
%    0 <= Sigma*1=1
%    0 <= (Sigma)_{i,j} <= 1
%
%   Copyright (c) 2013 Sira Ferradans

N = size(X,1);
d = size(X,2);
P = size(Gx,1);

options.null = 0;
Cost = getoptions(options, 'Cost', []);
if isempty(Cost)
    % use L^2 Cost  C_{i,j}=|X_i-Y_j|^2
    Cost = (repmat( sum(X'.^2)', [1 N] ) + ...
            repmat( sum(Y'.^2) , [N 1] ) - 2*X*Y'); 
end

D = repmat([1 sparse(1,N*N+N)],[1 N]);
D = reshape(D(1:N*N*N)',N*N,N);
%R = repmat(Gx,[1 d])*(ComputeR(X) * D * ComputeR(ones(N,1))-ComputeR(Y, N)); 
R =ComputeL(Gx, d)*(ComputeR(X) * D * ComputeR(ones(N,1))-ComputeR(Y, N)); %D1X= G_x(diag(Sigma 1) X - Sigma Y), with size (Nd,N N)

%R= G_x(diag(Sigma 1) X - Sigma Y ), with size (Nd,N N)

%%%%%%%%%%%
%Sigma(:) V(:)
A = [ComputeL(ones(1,N)),    sparse(N,P*d);     % 1' Sigma 
     ComputeR(ones(N,1)),    sparse(N,P*d);     %    Sigma 1
     speye(N*N,N*N),         sparse(N*N,P*d);   %    Sigma_i,j
     ones(1,N*N),            sparse(1,P*d);     % 1' Sigma 1 = T
     R,                      -speye(P*d,P*d) ]; % R Sigma - V = 0

kappa_min = getoptions(options,'kappa_min',0);
 
Amin = [    ...
    ones(N,1); ...            %    1 <= Sigma*1  
    kappa_min*ones(N,1); ...  %    0 <= 1'*Sigma
    sparse(N*N,1); ...        %    0 <= Sigma_i,j
    T;...                     %    1' Sigma 1 = T
    sparse(P*d,1)];           %    0 <= R Sigma -V 
  
Amax = [    ...
    ones(N,1); ...          %  Sigma*1  <= 1
    kappa*ones(N,1); ...    %  1'*Sigma <= k
    ones(N*N,1);  ...       %  Sigma_i,j<= 1
    T;...
    sparse(P*d,1)];         %  R Sigma -V <= 0

%%%%%%%%%%%%
%QP: <Cost,Sigma>+<Id,V> + ||R Sigma -V||
%Sigma(:) V
Xmin = [sparse(N*N,1) ; -inf*ones(P*d,1)];
Xmax = [inf*ones(N*N,1); inf*ones(P*d,1)];
C = [sparse(Cost(:)); sparse(P*d,1)];

%%
% Setup Mosek variables.
prob.qosubi = N*N+1:N*N+P*d;
prob.qosubj = N*N+1:N*N+P*d;
prob.qoval  = ones(P*d,1);

prob.c = C;
prob.a = sparse(A);
prob.blc = sparse(Amin);
prob.buc = sparse(Amax);
prob.blx = sparse(Xmin);
prob.bux = sparse(Xmax);

% Set parameters.
param = [];
% max number of iterations
param.MSK_IPAR_INTPNT_MAX_ITERATIONS = getoptions(options, 'qpprog_niter', 1000);
% tolerance, primal
param.MSK_DPAR_INTPNT_TOL_PFEAS = getoptions(options, 'qpprog_tol', 1e-14);
param.MSK_DPAR_INTPNT_TOL_REL_GAP = getoptions(options, 'qpprog_tol', 1e-14);

% verbosity level, 0=nothing is echoed, 3=all is echoed
verb = getoptions(options, 'verbose', 0);
[r,res] = mosekopt(['echo(' num2str(verb) ') minimize info'], prob,param);
if r~=0
    warning(['Mosek problem: ' res.rcodestr]);
end

err.niter = res.info.MSK_IINF_INTPNT_ITER;

w = res.sol.itr.xx;
% w=Sigma(:);
Sigma = reshape( w(1:N*N), [N N] );









  
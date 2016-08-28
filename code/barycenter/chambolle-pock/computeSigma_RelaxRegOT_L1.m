function [Sigma,Z,err] = computeSigma_RelaxRegOT_L1(X,Y,Gx,KX,lambda1,options)
% perform_ot_regul_linprog - solve the regularized OT
%
%   Sigma = perform_ot_regul_linprog(X,Y,Gx,k,lambda1,lambda2,options);
% Solve:
%       min <C,Sigma> + |Gx*(diag(Sigma*1)X-Sigma*Y)|_1 
% subject to
%       kx<=Sigma*1<=KX 
%       Sigma'*1 = 1
%       1*Sigma*1 = Total 
%  (0 <= Sigma_i,j <= 1 )? 
%
% restrictions on k: 0<=kx<=KX 
%
%   Copyright (c) 2012 Gabriel Peyre and Sira Ferradans

%%
% Recast the minimization as a LP on w=[Sigma(:); Z(:); W(:)]
%   min <C,Sigma> + lambda1*<Z,1>
%   -Z <= Gx*(Sigma*Y - diag(Sigma 1)X) <= Z, i.e. Z=|Gx*(Sigma*Y - diag(Sigma 1)X)|

N = size(X,1);
d = size(X,2);
P = size(Gx,1);

options.null = 0;
kx = getoptions(options,'ksum_min',0);
Total = getoptions(options,'Total',N);

Cost = getoptions(options, 'Cost', []);
if isempty(Cost)
   % use L^2 Cost  C_{i,j}=|X_i-Y_j|^2
    Cost = (repmat( sum(X'.^2)', [1 N] ) + ...
            repmat( sum(Y'.^2) , [N 1] ) - 2*X*Y'); 
end

D = repmat([1 zeros(1,N*N+N)],[1 N]);
D = reshape(D(1:N*N*N)',N*N,N);
D1X =ComputeL(Gx, d)*(ComputeR(Y, N) - ComputeR(X) * D * ComputeR(ones(N,1))); %D1X= G_x( Sigma Y - diag(Sigma 1) X), with size (Nd,N N)

%Sigma(:) Z
A = [ ...
     ComputeR(ones(N,1)),    sparse(N,P*d); ...    %  Sigma 1
     ComputeL(ones(1,N)),    sparse(N,P*d); ...    % 1' Sigma 
     ones(1,N*N)        ,    sparse(1,P*d); ...    % 1' Sigma 1
     speye(N*N,N*N)     ,    sparse(N*N,P*d); ...  % Sigma_i,j   
     D1X                ,    speye(P*d,P*d); ...   % Gx(diag(Sigma 1)X - Sigma Y) +Z   
     D1X                ,    -speye(P*d,P*d)];     % Gx(diag(Sigma 1)X - Sigma Y) -Z   
     
 
Amin = [    ...
    kx*ones(N,1); ...         %    kx <= Sigma*1
    ones(N,1); ...            %    1 * Sigmas = 1
    N ; ...                   %    N <= 1' Sigma 1
    sparse(N*N,1); ...        %    0 <= Sigma_i,j
    sparse(P*d,1); ...        %    diag(Sigma 1)X - Sigma Y +Z
    -inf+sparse(P*d,1)];      %    diag(Sigma 1)X - Sigma Y -Z

Amax = [    ...
    KX*ones(N,1); ...          %  Sigma*1  <= k
    ones(N,1); ...             %  1 * Sigmas = 1
    Total;                     %  1' Sigma 1< N
    ones(N*N,1); ...           %  Sigma_i,j<= 1
    inf+sparse(P*d,1); ...     %  diag(Sigma 1)X - Sigma Y +Z 
    sparse(P*d,1)];            %  diag(Sigma 1)X - Sigma Y -Z 
 
     
%Sigma(:) Z
Xmin = [sparse(N*N,1); -inf+sparse(P*d,1)];
Xmax = [];
C = [Cost(:); lambda1*ones(P*d,1)];

%%
% Setup Mosek variables.

prob.c = C;
prob.a = A;
prob.blc = Amin;
prob.buc = Amax;
prob.blx = Xmin;
prob.bux = Xmax;

%%
% Set parameters.

param = [];
% max number of iterations
param.MSK_IPAR_INTPNT_MAX_ITERATIONS = getoptions(options, 'linprog_niter', 10000);
% tolerance, primal
param.MSK_DPAR_INTPNT_TOL_PFEAS = getoptions(options, 'linprog_tol', 1e-12);
param.MSK_DPAR_INTPNT_TOL_REL_GAP = getoptions(options, 'linprog_tol', 1e-12);
% verbosity level, 0=nothing is echoed, 3=all is echoed
verb = getoptions(options, 'verbose', 0);

% Perform the optimization.
[r,res] = mosekopt(['echo(' num2str(verb) ') minimize info'], prob, param);
if r~=0
    warning(['Mosek problem: ' res.rcodestr]);
    disp(['But GAP=' num2str(res.info.MSK_DINF_MIO_OBJ_REL_GAP)]);
end
if strcmp(res.sol.bas.prosta,'PRIMAL_AND_DUAL_FEASIBLE')== 0
    warning(['Infeasible problem!:' res.sol.bas.prosta])
end 
err.niter = res.info.MSK_IINF_INTPNT_ITER;
err.rcodestr = res.rcodestr;
sol = res.sol.itr;
w   = sol.xx;
% w=[Sigma(:) Z];
Sigma = reshape( w(1:N*N), [N N] );
Z = reshape(  w(N*N+1:N*N+P*d), [P d] );


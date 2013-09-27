function [Sigma,S,Z,err] = perform_ot_regul_linprog(X,Y,G,k,lambda,options)

% perform_ot_regul_linprog - solve the regularized OT
%
%   Sigma = perform_ot_regul_linprog(X,Y,G,k,lambda,options);
%
% Solve:
%       min <C,Sigma> + |G*(X-Sigma*Y)|_1
% subject to
%       Sigma*1=1
%       Sigma'*1<=k
%
%   Copyright (c) 2012 Gabriel Peyre


%%
% Recast the minimization as a LP on w=[Sigma(:); S(:); Z(:)]
%   min <C,Sigma> + lambda*<Z,1>
%   S = Sigma*Y
%   -Z <= G*X - G*S <= Z, i.e. Z=|G*X - G*S|

N = size(X,1);
d = size(X,2);
P = size(G,1);
kmin=getoptions(options, 'kmin',0);
options.null = 0;
Cost = getoptions(options, 'Cost', []);
if isempty(Cost)
    pcost = getoptions(options,'pcost',2);
    % use L^2 Cost  C_{i,j}=|X_i-Y_j|^2
    Cost = (repmat( sum(X'.^2)', [1 N] ) + ...
            repmat( sum(Y'.^2) , [N 1] ) - 2*X*Y').^(pcost/2); 
end

flat = @(X)X(:);

A = [   ComputeR(ones(N,1)),    sparse(N,N*d),      sparse(N,P*d); ...
        ComputeL(ones(1,N)),    sparse(N,N*d),      sparse(N,P*d); ...
        sparse(P*d,N*N),        ComputeL(G, d),   -speye(P*d,P*d); ...
        sparse(P*d,N*N),        ComputeL(G, d),   +speye(P*d,P*d); ...
        ComputeR(Y, N),        -speye(N*d,N*d),     -sparse(N*d,P*d)];
        
    
Amin = [    ...
    ones(N,1); ...            % (C1)
    kmin*ones(N,1); ...           % (C2)
    zeros(P*d,1)-Inf;   ...   % (C3)
    flat(G*X);   ...          % (C4)
    zeros(N*d,1)];

Amax = [    ...
    ones(N,1); ...            % (C1)
    k*ones(N,1); ...          % (C2)
    flat(G*X);   ...          % (C3)
    zeros(P*d,1)+Inf;   ...   % (C4)
    zeros(N*d,1)];

Xmin = [ zeros(N*N,1); -Inf + zeros((N+P)*d,1) ];
Xmax = [];
C = [Cost(:); zeros(N*d,1); lambda*ones(P*d,1)];

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
param.MSK_IPAR_INTPNT_MAX_ITERATIONS = getoptions(options, 'linprog_niter', 100);
% tolerance, primal
param.MSK_DPAR_INTPNT_TOL_PFEAS = getoptions(options, 'linprog_tol', 1e-12);
param.MSK_DPAR_INTPNT_TOL_REL_GAP = getoptions(options, 'linprog_tol', 1e-12);
% verbosity level, 0=nothing is echoed, 3=all is echoed
verb = getoptions(options, 'verbose', 0);

% Perform the optimization.
[r,res] = mosekopt(['echo(' num2str(verb) ') minimize info'], prob, param);
if r~=0
    warning(['Mosek problem: ' res.rcodestr]);
end
err.niter = res.info.MSK_IINF_INTPNT_ITER;
sol = res.sol.itr;
w   = sol.xx;
% w=[Sigma(:); S(:); Z(:)];
Sigma = reshape( w(1:N*N), [N N] );
S = reshape(  w(N*N+1:N*N+N*d), [N d] );
Z = reshape(  w(N*N+N*d+1:end), [P d] );

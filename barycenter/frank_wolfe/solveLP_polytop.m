function [Sigma,err] = solveLP_polytop(C,T,k,kmin,options)

% solveLP_polytop - solve a linear program, where the variable is a
% bistochastic matrix.
%
%   Sigma = solveLP_polytop(C,T,options)
%
% Solve:
%       min <C,Sigma> 
% subject to
%       kmin<=Sigma*1<=k            (C1)
%       Sigma'*1=1                  (C2)
%       Sigma_i,j in [0,1]          (C3)
%       1'*Sigma*1 = T              (C4)
%
%  Assuming size(C)=[N,N]

%   Copyright (c) 2013 Sira Ferradans and Gabriel Peyre

N = size(C,1);

options.null = 0;

A = [   ComputeR(ones(N,1)); ...
        ComputeL(ones(1,N)); ...
        speye(N*N,N*N); ... % Sigma_i,j 
        ones(1,N*N)];
        
Amin = [    ...
    kmin*ones(N,1); ...       % (C1)
    ones(N,1); ...            % (C2)
    sparse(N*N,1);            % (C3)
    T];                       % (C4)

Amax = [    ...
    k*ones(N,1); ...          % (C1)
    ones(N,1); ...            % (C2)
    ones(N*N,1); ...          % (C3)
    T];                       % (C4)

Xmin = zeros(N*N,1);          % (C3)
Xmax = [];                    % (C3)

%%
% Setup Mosek variables.
prob.c = C(:);
prob.a = A;
prob.blc = Amin;
prob.buc = Amax;
prob.blx = Xmin;
prob.bux = Xmax;

%%
% Set parameters.
param = [];
% max number of iterations
param.MSK_IPAR_INTPNT_MAX_ITERATIONS = getoptions(options, 'linprog_niter', 1000);
% tolerance, primal
param.MSK_DPAR_INTPNT_TOL_PFEAS = getoptions(options, 'linprog_tol', eps);
param.MSK_DPAR_INTPNT_TOL_REL_GAP = getoptions(options, 'linprog_tol', eps);
% verbosity level, 0=nothing is echoed, 3=all is echoed
verb = getoptions(options, 'verbose', 0);

% Perform the optimization.
[r,res] = mosekopt(['echo(' num2str(verb) ') minimize info'], prob, param);
if r~=0
    warning(['Mosek problem: ' res.rcodestr]);
end
if (strcmp(res.sol.bas.prosta,'PRIMAL_AND_DUAL_FEASIBLE') == 0)
    warning(['Infeasible problem!:' res.sol.bas.prosta])
end
err.niter = res.info.MSK_IINF_INTPNT_ITER;
sol = res.sol.itr;
w   = sol.xx;
% w=[Sigma(:)];
Sigma = reshape( w(1:N*N), [N N] );


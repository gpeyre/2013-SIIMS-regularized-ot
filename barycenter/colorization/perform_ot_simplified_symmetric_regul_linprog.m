function [Sigma,Z,W,err] = perform_ot_simplified_symmetric_regul_linprog(X,Y,Gx,Gy,k,Total,lambda1,lambda2,options)

% perform_ot_regul_linprog - solve the regularized OT
%
%   Sigma = perform_ot_regul_linprog(X,Y,Gx,Gy,k,lambda1,lambda2,options);
%
% Solve:
%       min <C,Sigma> + |Gx*(diag(Sigma*1)X-Sigma*Y)|_1 + |Gy*(diag(1*Sigma)Y-Sigma'*X) |_1
% subject to
%       Sigma*1<=k 
%       Sigma'*1<=k
%       1*Sigma*1 = N 
%  (0 <= Sigma_i,j <= 1 )? 
%
%   Copyright (c) 2012 Gabriel Peyre and Sira Ferradans


%%
% Recast the minimization as a LP on w=[Sigma(:); Z(:); W(:)]
%   min <C,Sigma> + lambda1*<Z,1> + lambda2*<W,1>
%   -Z <= Gx*(Sigma*Y - diag(Sigma 1)X) <= Z, i.e. Z=|Gx*(Sigma*Y - diag(Sigma 1)X)|
%   -W <= Gy*( diag(1 Sigma)Y - Sigma'*X) <= W, i.e. W=|Gy*(diag(1 Sigma)Y - Sigma'*X)|


N = size(X,1);
d = size(X,2);
P = size(Gx,1);

options.null = 0;
Cost = getoptions(options, 'Cost', []);
if isempty(Cost)
    pcost = getoptions(options,'pcost',2);
    % use L^2 Cost  C_{i,j}=|X_i-Y_j|^2
    Cost = (repmat( sum(X'.^2)', [1 N] ) + ...
            repmat( sum(Y'.^2) , [N 1] ) - 2*X*Y').^(pcost/2); 
end

D = repmat([1 zeros(1,N*N+N)],[1 N]);
D = reshape(D(1:N*N*N)',N*N,N);
D1X =ComputeL(Gx, d)*(ComputeR(Y, N) - ComputeR(X) * D * ComputeR(ones(N,1))); %D1X= G_x( Sigma Y - diag(Sigma 1) X), with size (Nd,N N)

Trans=perm_mat(d,N);
D1Y =ComputeL(Gy, d)*(Trans*ComputeR(X, N) - ComputeR(Y) * D * ComputeL(ones(1,N)));%D1Y= G_y(Sigma^* X - diag(1^* Sigma)Y), with size (Nd,N N)


%Sigma(:) Z W
A = [ ...
     ComputeR(ones(N,1)),    sparse(N,P*d), sparse(N,P*d); ...    %  Sigma 1
     ComputeL(ones(1,N)),    sparse(N,P*d), sparse(N,P*d); ...    % 1' Sigma 
     ones(1,N*N)        ,    sparse(1,P*d), sparse(1,P*d);            % 1' Sigma 1
     speye(N*N,N*N)     ,    sparse(N*N,P*d), sparse(N*N,P*d); ... % Sigma_i,j 
     
     D1X                ,    speye(P*d,P*d),  sparse(P*d,P*d); ... % diag(Sigma 1)X - Sigma Y +Z   
     D1X                ,    -speye(P*d,P*d), sparse(P*d,P*d);     % diag(Sigma 1)X - Sigma Y -Z   
     D1Y                ,    sparse(P*d,P*d), speye(P*d,P*d); ...  % diag(1 Sigma)Y - Sigma X +W   
     D1Y                ,    sparse(P*d,P*d), -speye(P*d,P*d)];    % diag(1 Sigma)Y - Sigma X -W   
 
epsilon = 0;
Amin = [    ...
    epsilon*ones(N,1); ...    %    0 <= Sigma*1  
    epsilon*ones(N,1); ...    %    0 <= 1'*Sigma
    Total; ...                %    N <= 1' Sigma 1 
    sparse(N*N,1); ...        %    0 <= Sigma_i,j
    
    sparse(P*d,1); ...        %   diag(Sigma 1)X - Sigma Y +Z
    -inf+sparse(P*d,1); ...      %   diag(Sigma 1)X - Sigma Y -Z
    sparse(P*d,1); ...        %   0<=diag(1 Sigma)Y - Sigma X +W 
    -inf+sparse(P*d,1)];      %   diag(1 Sigma)Y - Sigma X -W 
 
   
    

Amax = [    ...
    k*ones(N,1); ...          %  Sigma*1  <= k
    k*ones(N,1); ...          %  1'*Sigma <= k
    Total;                    %  1' Sigma 1< N
    ones(N*N,1); ...          %  Sigma_i,j<= 1

    inf+sparse(P*d,1); ...    %  diag(Sigma 1)X - Sigma Y +Z 
    sparse(P*d,1); ...          %  diag(Sigma 1)X - Sigma Y -Z 
 
    inf+sparse(P*d,1); ...    %  diag(Sigma 1)Y - Sigma X +W 
    sparse(P*d,1)];           %  diag(Sigma 1)Y - Sigma X -W 
     

%Sigma(:) S T Z
Xmin = [sparse(N*N,1); -inf+sparse(2*P*d,1)  ];
Xmax = [];
C = [Cost(:); lambda1*ones(P*d,1); lambda2*ones(P*d,1) ];

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
verb = getoptions(options, 'verbose', 3);

% Perform the optimization.
[r,res] = mosekopt(['echo(' num2str(verb) ') minimize info'], prob, param);
if r~=0
    warning(['Mosek problem: ' res.rcodestr]);
end
err.niter = res.info.MSK_IINF_INTPNT_ITER;
sol = res.sol.itr;
w   = sol.xx;
% w=[Sigma(:) Z W];
Sigma = reshape( w(1:N*N), [N N] );
Z = reshape(  w(N*N+1:N*N+P*d), [P d] );
W = reshape(  w(N*N+P*d+1:N*N+2*P*d), [P d] );

function P = perm_mat(p,n)
% permute applied to (p,n) matrix
[b,a] = meshgrid(1:n,1:p); a = a(:); b = b(:);
P = sparse(b + (a-1)*n, a + (b-1)*p, ones(n*p,1));


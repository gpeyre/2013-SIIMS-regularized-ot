function [x,R] = perform_admm(x, K,  KS, ProxFS, ProxH, Gupdate, options)

% perform_admm - preconditionned ADMM method
%
%    [x,R] = perform_admm(x, K,  KS, ProxFS, ProxH, options);
%
%   Solves 
%       min_x F(K*x) + H(x)
%   where F and H are convex proper functions with an easy to compute proximal operator, 
%   and where K is a linear operator
%
%   Uses the Preconditioned Alternating direction method of multiplier (ADMM) method described in 
%       Antonin Chambolle, Thomas Pock, 
%       A first-order primal-dual algorithm for convex problems with applications to imaging, 
%       Preprint CMAP-685	
%
%   INPUTS:
%   ProxFS(y,sigma) computes Prox_{sigma*F^*}(y)
%   ProxG(x,tau) computes Prox_{tau*G}(x)
%   K(y) is a linear operator.
%   KS(y) compute K^*(y) the dual linear operator.
%   options.sigma and options.tau are the parameters of the
%       method, they shoudl satisfy sigma*tau*norm(K)^2<1
%   options.theta=1 for the ADMM, but can be set in [0,1].
%   options.verb=0 suppress display of progression.
%   options.niter is the number of iterations.
%   options.report(x) is a function to fill in R.
%
%   OUTPUTS:
%   x is the final solution.
%   R(i) = options.report(x) at iteration i.
%
%   Copyright (c) 2010 Gabriel Peyre

options.null = 0;
report = getoptions(options, 'report', @(x)0);
niter = getoptions(options, 'niter', 500);
theta = getoptions(options, 'theta', 1);
verb = getoptions(options, 'verb', 0);
F = getoptions(options,'F',@(x)x);
H = getoptions(options,'H',@(x)x);
k = getoptions(options,'k',2);

if isnumeric(K)
    K = @(x)K*x;
end  
if isnumeric(KS)
    KS = @(x)KS*x;
end      

%%%% ADMM parameters %%%%
sigma = getoptions(options, 'sigma', -1);
tau   = getoptions(options, 'tau', -1);
if sigma<0 || tau<0
    [L,e] = compute_operator_norm(@(x)KS(K(x)),randn(size(x)));
    if sigma<0
        sigma = 500;
    end 
    tau = .9/(sigma*L);
end

y = K(x);

x1 = x;

clear R;
for i=1:niter    
    
    if verb
        progressbar(i,niter);
         
        % record energies
        R(i) = F(K(x))+H(x);
        ff(i)= sum((sum(x,2)-1.).^2);
        hh(i)= sum(max(sum(x,1)-k,0).^2);
        
     %   plot(R);title('Energy');drawnow;
       
    end

    % update
    xold = x;
%     if (kcell) 
%         y = ProxFS( Gupdate(y,sigma,K(x1)), sigma); 
%     else 
         y = ProxFS( y+sigma*K(x1), sigma);
%    end 
    
    x = ProxH(  x-tau*KS(y), tau);
    x1 = x + theta * (x-xold);
 
end
if verb
    % figure;
    % subplot(1,3,1);
 %   plot(R);title('Energy');draw now;
    % subplot(1,3,2);plot(ff,'r');title('D1')
    % subplot(1,3,3);plot(hh,'b');title('Ck')
else
    R=-1;%not saving values
end 

function [xk,err] = frankwolfmin_onpolytop(gradf,x0,ftau,T,kappa,kappa_min,lambda,options)

%This function minimizes f using the Frank-Wolf algorithm

% gradf: function that computes the gradient of f: output 
% f: function to minimize
% x0: initialization point for the algorithm. size(x0) = (N,1)
% tau: function that gives the optimal step size for each x_k

xk = x0;
numiter = getoptions(options,'numiter',1000000);
f = getoptions(options,'f',@(x)x);
threshold = getoptions(options,'threshold',eps);

options.linprog_niter = 10000;
options.linprog_tol = 1e-10;
options.verbose = 0;
err.LPtiming =[];
for i=1:numiter

    %step 1: find direction of descent min_{s in D} gradf(x_k)' s
    c = gradf(xk);
    tstart = tic;
    [sk, ~] = solveLP_polytop(c,T,kappa,kappa_min,options); %we are solving (gradf(xk)') s
    err.LPtiming(end+1)=toc(tstart);
    %step 2: get optimal step size
    if (lambda > 1e-5)
        tk = ftau(xk,sk);%,c);
        if ( (tk > 1) || (tk < 0) )
            disp(['error']);
        end 
    else 
        tk = min(max(sum((xk(:)-sk(:)).*c(:)),0),1);
    end 
    %step 3: update
    xk = xk + tk*(sk-xk);
    
    %Debug: control the evolution of the energy
     E(i) = f(xk);
%      subplot(1,2,1);
% plot(E);drawnow;
%      title('Energy')
%      vT(i) = tk;
%      subplot(1,2,2);plot(vT);drawnow;
%      title('\tau');
% %      
    % reached minimum
     if (i > 1 ) && ((E(i-1) - E(i)) < threshold)
         err.numiterations=i;
         return
     end 
end 
err.numiterations=i;
%E(i)



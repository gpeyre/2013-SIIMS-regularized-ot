function [X,E] = knnsubsampling(X0,Xinit, options)

% knnsubsampling - K-means algorithm
%
%   [X,E] = knnsubsampling(X0,Xinit);
%
%   options.niter is the number of iterations.

% clouds
[N0,~]= size(X0); 
[N,~] = size(Xinit);

options.null = 0;
niter = getoptions(options, 'niter', 20);
verb = getoptions(options, 'verb', 20);

mynorm = @(x)norm(x(:));
%%
% K-mean algorithm
X = Xinit;
for i=1:niter
    if verb
        progressbar(i,niter);
    end
    % Step 1 : NN
    I = knnsearch(X0,X,1);
    if nargout>1
        E(i) = mynorm( X0-X(I,:) );
    end
    % Step 2 : mean
    for k=1:N
        J = find(I==k);
        if isempty(J)
            % problem, use random re-sampling
            J = floor( rand*(N0-1) ) + 1;
        end
        X(k,:) = mean( X0(J,:), 1 );
    end
end


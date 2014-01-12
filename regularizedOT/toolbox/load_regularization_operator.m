function G = load_regularization_operator(x, options)

% load_regularization_operator - load a discrete gradient operator on graph
%
%   G = load_regularization_operator(x, options);
%
%   x is an (n,d) set of points in R^d.
%   If d=1, use finite differences on adjacent point, if d>1, 
%   uses a sparse options.nn nearest neighbors graph.
%
%   options.weight='distance': weight 1/d where d is pairwise distance
%   options.weight='constant': weight 1.
%
%   Copyright (c) 


n = size(x,1);
options.null = 0;
weight = getoptions(options, 'weight', 'constant');
   
if size(x,2)==1
    %%% 1D %%%
    d = x(2:end)-x(1:end-1);
    G = eye(n) - diag(ones(n-1,1),1); G(end,:) = [];
    disp('G devrait être l"opposé ?')
    if strcmp(weight, 'distance')
        G = diag(1./d)*G;
    end
    G = sparse(G);
elseif 1
    %%% Multi-dimension %%%
    k = getoptions(options, 'nn', 2);
    radius = getoptions(options, 'radius', inf);
    D = compute_distance_matrix(x);
    D(D>=radius) = inf;
    
    [D,I] = sort(D); 
    D = sqrt(D(2:k+1,:)); D=D(:);
    I = I(2:k+1,:); I = I(:);% NN
    J = repmat( 1:n, [k 1] ); J = J(:); % indices de l'identité (matrice diagonale sparse)
    p = length(I);    
    w = 1./D;
    
    epsilon = 10e-6;
    w_norm = 1/(sqrt(mean( D.^(-2) ))+epsilon); % * 1/min(D) ; % 1./mean(D); * min(D)
 
    if strcmp(weight, 'constant')
        w = 0*w + 1/w_norm;
    end
    
    % normalization
    w = (w*w_norm)/k; % /mean(1./D) /mean(D)
    
    %This is a non-symmetric operator
    G = sparse(1:p, I(:), w, k*n,n) - sparse(1:p, J(:), w, k*n,n); % matrice p(=k*n) x n
else
    %%% TEST JULIEN : other definition using signed weight %%%
    D = compute_distance_matrix(x); % D_ij = | x_i - x_j |^2
    Delta = repmat( reshape(x,[n d 1]), [1 1 n]) - repmat( reshape(x',[1 d n]), [n 1 1]); % Delta_ij = x_i - x_j 
    
    % NN computation
    [val,I] = sort(D);
    k = getoptions(options, 'nn', 5);
    
    % compute inverse and remove division by 0 and keep only NN values
    DD = D.^(-1); DD([(1:n)+n*(0:n-1)])=0; 
    DD( I(k+2:n,:) + repmat((0:n-1)*n,[n-k-1 1]) ) = 0;
    DD=DD';
    DD( I(k+2:n,:) + repmat((0:n-1)*n,[n-k-1 1]) ) = 0;
    DD=DD';
    %figure, imagesc(DD==0)
    
    'Problème : on peut se retrouver avec moins de k voisins si on symetrise'
    'Solution : définir les voisinages apres le reshape'
    
    G = - repmat( reshape(DD,[n 1 n]), [1 d 1]).*Delta; % matrix G_ij  = - (x_i - x_j) / || x_i - x_j ||^2 \in \R^d
    %G = sparse(G);
    
    if strcmp(weight, 'constant')
        G(G>0) = 1;
        G(G<0) = -1;
    end
    
    % evaluation du gradient de x à partir de G
    for i=1:n
      Gx(i,:,:) = reshape(G(i,:,:),[d n])*x;
    end
    % visualisation du gradient 2D
    if 1 && d==2
      figure, plot(x(:,1),x(:,2),'b+'), hold on, 
      u1 = squeeze(Gx(:,1,:)); u2 = squeeze(Gx(:,2,:));
      % TEST: line([x(:,1)';x(:,1)'+randn(1,n)],[x(:,2)';x(:,2)'+randn(1,n)])
      line([x(:,1)';x(:,1)'+u1(:,1)'],[x(:,2)';x(:,2)'+u1(:,2)'])
    end
    
    % evaluation de l'adjoint à partir de G
    'à finir'
    DIV = 0;
    Gxx = zeros(n,d);
    for i=1:n
      Gxx = Gxx + reshape(G(i,:,:),[d n])' * reshape(Gx(i,:,:),[d d]);
    end
    % figure, plot(x(:,1),x(:,2),'b+')
    % hold on, plot(squeeze(Gx(:,1,:)),squeeze(Gx(:,2,:)),'rd')
    
end
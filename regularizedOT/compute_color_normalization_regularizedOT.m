load_mosek();

addpath('./toolbox/');

%% Load images
name = 'parrot'; %will load the files synthetic-1.jpg and synthetic-2.jpg

% image file
Y0 = double(imread(['../images/' name '-2.jpg']))/256;
X0 = double(imread(['../images/' name '-1.jpg']))/256;
N0 = size(X0,1);
folder = 'results';
mkdir(folder)
%subsample
N = 8; %number of points in the small image

indx = round(linspace(1,N0,N));
X = reshape(X0(indx,indx,:),N*N,3);  
Y = reshape(Y0(indx,indx,:),N*N,3); 

X0 = reshape(X0,N0*N0,3);
Y0 = reshape(Y0,N0*N0,3);

disp(['Knn clustering on X0'])
X = knnsubsampling(X0,X);
disp(['Knn clustering on Y0'])
Y = knnsubsampling(Y0,Y);

%% parameters and helper functions
lw =1;ms = 15;
myplot3 = @(X,col)plot3(X(:,1),X(:,2),X(:,3), 'o', 'MarkerSize', ms, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col, 'LineWidth', lw);
h1 = @(X)X./max(X(:)); % to imshow

% H are the original pixels of the first image
% quantify(H) are the value of the new image
% the OT is X -> U
getIdx = @(X,x)knnsearch(X,x,1);
interpolate=@(Idx,X,x,u)reshape( u(Idx,:) - x(Idx,:)+ X, [N0 N0 3] ) ;
quantify =@(X,x,u)interpolate(getIdx(X,x),X,x,u);
[N,d] = size(X);

%parameters
nnx = 4;    %neighbours in the nearest neighbour graph
weightedG = true; %weighted gradient operator?
pcost = 2; %norm to compute C
lambda = 1e-3;
ksum = 1;

%% modification of x and y
optionsG.nnx = nnx;
optionsG.weightedG = weightedG;
optionsG.plot = false;
optionsG.color_plot='b';
optionsG.epsilon = 0.1;
optionsG.name = name;

%% Compute Barycenter
%compute graph
disp(['Compute graph on X'])
Gx = computeKnnGraph(X,optionsG);
disp(['Compute graph on Y'])
Gy = computeKnnGraph(Y,optionsG);

%2.-compute the mapping
options.linprog_niter = 1000;
options.linprog_tol = 1e-12;
options.verbose = 0;
options.pcost = pcost;

%        [kx, KX, ky,KY, M,  lambda1,lambda2]
params = [0.1 ksum 0.1 ksum  N lambda lambda]; %set of parameters. For OT:[1 1 1 1 N 0 0]
% compute relaxed and regularized OT
disp(['Computing relaxed and regularized OT']);
[Sigma,Z,W,err]=perform_ot_symmetric_regul_linprog(X,Y,Gx,Gy,...
        params(1), params(2),params(3),params(4),params(5), params(6),params(7),options);

options.displayU = false;
h=plotResults(X,Y,Sigma, options);xlabel('R');ylabel('G');zlabel('B');
view([0 90]);

namefile=['./' folder '/mapping' name 'XY.eps'];
print('-depsc',namefile); 

%%transport X
SigmaX = Sigma./repmat(sum(Sigma,2),[1 N]); % X = SigmaX*Y;
h=figure;imshow(h1(quantify(X0,X,SigmaX*Y)));
namefile=['./' folder '/' name 'X.eps'];
print('-depsc',namefile); 

%transport Y
SigmaY = Sigma./repmat(sum(Sigma,1),[N 1]); % Y = SigmaY*X;
h=figure;imshow(h1(quantify(Y0,Y,SigmaY'*X)));
namefile=['./' folder '/' name 'Y.eps'];
print('-depsc',namefile); 

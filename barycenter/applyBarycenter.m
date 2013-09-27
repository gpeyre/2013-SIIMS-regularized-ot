addpath('./toolbox/');
addpath('./frank_wolfe/');
addpath('./qp/');
addpath('./chambolle-pock/');

load_mosek();

vector =@(X)reshape(X,size(X,1)*size(X,2),3);

%%load image
name = 'flowers';% to be loaded: 'name'-1.jpeg and 'name'-2.jpeg
three = false;%to obtain the barycenter of 3 images, set to true

% image file
X0 = double(imread(['../images/' name '-1.jpg']))/256;
Y0 = double(imread(['../images/' name '-2.jpg']))/256;

N0 = size(X0,1);%N_0

%% Construct data set X,Y by subsampling X0,Y0
%we subsample to compute the optimal mapping
N = 8; %N:number of samples small image

indx = round(linspace(1,N0,N));
X = reshape(X0(indx,indx,:),N*N,3);  
Y = reshape(Y0(indx,indx,:),N*N,3); 

X0 = reshape(X0,N0*N0,3);
Y0 = reshape(Y0,N0*N0,3);
disp('knn-subsampling on X')
X = knnsubsampling(X0,X);
disp('knn-subsampling on Y')
Y = knnsubsampling(Y0,Y);

X = cat(3,X,Y);
X0 = cat(3,X0,Y0); 

if (three)
    Z0 = double(imread(['../images/' name '-3.jpg']))/256;
    Z = reshape(Z0(indx,indx,:),N*N,3); 
    Z0 = reshape(Z0,N0*N0,3);
    Z = knnsubsampling(Z0,Z);
    X = cat(3,X,Z);
    X0 = cat(3,X0,Z0);
end 
clear Y;

[N,d,R] = size(X);

%% Interpolation method
% H are the original pixels of the first image
% quantify(H) are the value of the new image
% the OT is X -> U
%quantify = @(X,x,u)reshape( u(knnsearch(X,x,1),:), [l l 3] );

getIdx = @(X0,x)knnsearch(X0,x,1);
interpolate=@(Idx,X0,x,u)reshape( u(Idx,:) - x(Idx,:)+ X0, [N0 N0 3] ) ;
quantify =@(X0,x,u)interpolate(getIdx(X0,x),X0,x,u);

%functions helpers and parameters of the method
h1 = @(X)X./max(X(:)); % function for imshow
k = 1.1;               %capacity
lambda = 0.101;        %regularization parameter
M = N;                 %total capacity 
nnx = 4;               %neighbours in the nearest neighbour graph
weightedG = true;      %weighted gradient operator?

%% Compute the graphs on the datasets
%params for creating the graph
optionsG.nnx = nnx;
optionsG.weightedG = weightedG;
optionsG.plot = false;
optionsG.name = name;

Gs = zeros(N*nnx,N,R);
optionsG.epsilon=0.3;
for i=1:size(X,3)
    Gs(:,:,i) = computeKnnGraph(X(:,:,i),optionsG); 
end


%% Compute the barycenter
folder = [ './results/images/' name ];
mkdir(folder)
 
Xb = X(:,:,2);
S = [];XX=[];
for irho=[0 0.3 0.6 1]      %to create a evolution from colors of X1 to X2
    rho = [irho 1-irho];    %to compute between 3 images, need to have a rho with 3 values  
    
    options.ksum = k;       % if we want OT, this k = 1, ksum_min=1
    options.ksum_min = 0.1; 
    options.lambda = lambda;
    options.niter = 300;
    options.qpprog_tol = 1e-5;
    options.Xb = X(:,:,round(1+max(rho))); %initialization
    options.qp=true;
   
    options.T = M; %Total mass

    %% Compute the cloud of points that is the barycenter
    disp(['Create barycenter point cloud for rho=' num2str(rho(1)) ]);
    [Xb,Sigmas]=computeBarycenterRegularizedDiag(X,Gs,rho,options);

    %% Colorization: compute mapping from each data set to the barycenter
    disp(['Graph on the barycenter point cloud'])
    Gxb = computeKnnGraph(Xb,optionsG); %compute Graph on the barycenter
    optionsColor.linprog_niter = 1000;
    optionsColor.linprog_tol = 1e-12;
    optionsColor.verbose = 0;
    
    % computing match between barycenter and each cloud of points
    for r=1:size(X,3)
    
        V = X(:,:,r);
        Gv= Gs(:,:,r);
        Y0 = X0(:,:,r);
        
        %2.-compute the mapping
        k_bar = 1.1;            %mass parameter
        optionsColor.kmin=0.1;  %for OT set k_bar=kmin=1
        lambda_bar= 0.009;      %regularization parameter
        disp(['Compute the regul. transport between X' num2str(r) ' to the barycenter']);
        [Sigma,S,Z0,err] = perform_ot_regul_linprog(V,Xb,Gv,k_bar,lambda_bar,optionsColor);
        u = Sigma*Xb;
        %3.- Save image
        h=figure;imshow(h1(quantify(Y0,X(:,:,r),u)));
        namefile=[folder '/result_' name num2str(r) 'rho' num2str(rho(1)) ];
    
        disp(['saving result in ' namefile '.png']);
        saveas(h,['./' strrep(namefile,'.','') '.png']);
       % print('-depsc',['./' strrep(namefile,'.','') '.eps']);
    end
   
  
   
end %end for rho 
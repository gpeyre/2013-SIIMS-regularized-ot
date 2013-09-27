ms = 8.5; lw =1.5;
myplot = @(X,col)plot(X(:,1),X(:,2), 'o', 'MarkerSize', ms, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col, 'LineWidth', lw);
  

addpath('./toolbox/');
addpath('./frank_wolfe/');
load_mosek();

%Energy plot/test
%Define database 
%good to show the several minimum
% x=[0 1.1;0 0.05;0 1];
% x2=[1 0.1;1 0; 1 1.05];
% name = '3points';

%cluster matching
%x=[rand(10,2)*[0.1 0;0 0.4]+repmat([0 0.5],[10 1]);...
%  rand(5,2)*0.1];

%x2=[rand(5,2)*0.1+repmat([0.75 0.6],[5 1]); ...
%    rand(10,2)*0.2+repmat([0.75 0],[10 1])];

  load 'x_clusters';
%  x=[x;cat(3,[0.05 0.4],[0.8 0.2])]
%  x=[x;cat(3,[0.05 0.15],[0.8 0.6])]
%  x=[x;cat(3,[0.05 0.25;0.05 0.25],[0.85 0.3;0.85 0.5])]
% 
% name='clustersMorePointsDiag';

% name = 'PointsRegInit';
% cuad =[0 0; 0 0.2; 0.2 0; 0.2 0.2];
% trian=[0 0; 0.2 0; 0.1 0.1];
% x = [cuad; trian+repmat([0 1],[3 1])];% -0.1 0.5]; 
% x2= [cuad+repmat([1 1],[4 1]); trian+repmat([1 0],[3 1])];% 1.4 0.55];

%x=cat(3,x,x2);

[N d R]=size(x);

%% Compute the barycenter
% options.T = N;
% options.ksum = 2;
% options.lambda = 0;
% options.niter = 1000;
% options.threshold = 1e-5;
% 
% U = x(:,:,2);
% 
% %Definimos las Sigmas Optimas
% Sigmas1 = [0 0 0.5;1 1 0;0 0 0.5];
% Sigmas1 = cat(3,Sigmas1,eye(3));
% figure;
% [vrho,E1]=testingSmoothness(Sigmas1,x,x(:,:,2));
%   
% Sigmas2 = eye(3);
% Sigmas2 = cat(3,Sigmas2,[0 0.5 0 ;0 0.5 0;1 0 1]);
% figure;
% [vrho,E2]=testingSmoothness(Sigmas2,x,x(:,:,1));
% figure;
% plot(vrho,E1,'r');hold on;plot(vrho,E2,'b')
% legend('R','B')
% 
% %%Testing Sigma computation
 rho = 0.1;
% options.T =N;
% options.ksum=2;
% options.Xb =x(:,:,2);
% [Xb1,S1,E1]=computeBarycenter(x,[rho 1-rho],options);
% 
% figure;
% options.Xb =x(:,:,1);
% [Xb2,S2,E2]=computeBarycenter(x,[rho 1-rho],options);
% 
% 
% figure;
% options.threshold = 0.1;
% plotResults(x(:,:,1),Xb1,S1(:,:,1),options);
% hold on;title(['Energy ' num2str(E1)])
% plotResults(x(:,:,2),Xb1,S1(:,:,2),options);
% myplot(x(:,:,1), 'b');hold on;
% myplot(x(:,:,2), 'r');hold on;
% myplot(Xb1, 'g');hold on;
% figure;
% 
% plotResults(x(:,:,1),Xb2,S2(:,:,1),options);
% hold on;plotResults(x(:,:,2),Xb2,S2(:,:,2),options);
% title(['Energy ' num2str(E2)])
% myplot(x(:,:,1), 'b');hold on;
% myplot(x(:,:,2), 'r');hold on;
% myplot(Xb2, 'g');hold on;


nnx = min(4,size(x,1));
optionsG.nnx = nnx;
optionsG.weightedG = false;
optionsG.plot=true;
optionsG.epsilon=0.21;%0.31
Gs = zeros(N*nnx,N,R);
h=figure;
for r=1:R
    Gs(:,:,r) = computeKnnGraph(x(:,:,r),optionsG); 
end 
file=['./results/' name 'Graph-nnx' num2str(nnx) '.png']; 
 saveas(h,file);

%Gs(2,:,1)=zeros(3,1);
%Gs(3,:,2)=zeros(3,1);
flat = @(x)x(:);

%% Test with lambda
vrho=0:0.1:1;
lambda = 0.5;
d = ['TesteandoSi-lambda' num2str(lambda)];
mkdir(['./results/' d])
    
options.lambda=lambda;
options.ksum = 2;
Xb1=[];S11=[];
Xb2=[];S12=[];


options.Xb = x(:,:,2);%Xb1(:,:,end);
[Xb1,S11]=computeBarycenterRegularizedDiag(x,Gs,[0.5 0.5],options);
  
S11(S11<0.1) = 0;
S11(:,:,1) = S11(:,:,1)*diag(1./sum(S11(:,:,1)));
S11(:,:,2) = S11(:,:,2)*diag(1./sum(S11(:,:,2)));

C = @(x,y)(repmat( sum(x.^2,2), [1 N] ) + repmat( sum(y.^2,2)' , [N 1] )- 2*x*y');
flat=@(x)x(:);

for rho=vrho
     %% Obtain the associated Xb
    H= @(X)Energy(X,S11,rho,x,Gs,lambda,C);
    optionsBFGS.iterations=1000;
    [Xb,~,info] = perform_bfgs(@(X)deal(H(X),gradH(X,S11,Gs,lambda,x,rho)),flat(x(:,:,2)),optionsBFGS);
    Xb1(:,:,end+1) = reshape(Xb,N,d);

    h=figure;

    plotResults(x(:,:,1),Xb1(:,:,end),S11(:,:,1,end),options);
    hold on;plotResults(x(:,:,2),Xb1(:,:,end),S11(:,:,2,end),options);
    myplot(x(:,:,1),'b');
    myplot(x(:,:,2),'g');
    title(['Init: X2']);
    file=['./results/' d '/' name 'X2-rho' num2str(rho) '-lambda' num2str(lambda) 'nnx' num2str(nnx) '.png']; 
    saveas(h,file);
    figure;plot(sum(S11(:,:,1,end)),'r');hold on;
    plot(flat(sum(S11(:,:,2,end))),'g')
    file=['./results/' d '/weights' name 'X2-rho' num2str(rho) '-lambda' num2str(lambda) 'nnx' num2str(nnx) '.png']

end




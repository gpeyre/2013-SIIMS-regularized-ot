load_mosek();

%%load image
name='cluster_matching';
load 'cluster_matching.mat'

hold off;myplot(x,'b');hold on; myplot(y,'r');axis equal

N = size(x,1);
d = size(x,2);

%parameters
nnx = 4;    %neighbours in the nearest neighbour graph
weightedG = true; %weighted gradient operator?
pcost = 2; %norm to compute C
M = N;

%% modification of x and y
optionsG.nnx = nnx;
optionsG.weightedG = weightedG;
optionsG.plot = true;
optionsG.color_plot='b';
optionsG.name = name;
optionsG.epsilon=0.4;

disp(['Computing graphs on X'])
Gx = computeKnnGraph(x,optionsG);
disp(['and Y'])
Gy = computeKnnGraph(y,optionsG);

%2.-compute the mapping
options.linprog_niter = 1000;
options.linprog_tol = 1e-12;
options.verbose = 0;
options.pcost = pcost;



       % kx KX ky KY M lambdaX lambdaY
params =[1 1 1 1 M 0 0; ... %OT
           1 1 0.1 1.5 M 0 0;...
           1 1 0 2 M 0 0;...
           0 2 1 1 M 0 0;...
           1 1 0.1 10 M 0 0;...
           0.1 10 0.1 10 M 0 0];

   
[Sigma,Z,W,err]=perform_ot_symmetric_regul_linprog(x,y,Gx,Gy,...
    params(1), params(2),params(3),params(4),params(5), params(6),params(7),options);

options.displayU = false;
h=figure;
plotResults(y,x,Sigma', options);
 namefile = ['./results_symmetric/' name  '_l' num2str(params(i,7)) 'kx' num2str(params(i,1)) 'KX' num2str(params(i,2)) ...
                'ky' num2str(params(i,3)) 'KY' num2str(params(i,4)) '_nn' num2str(nnx) '.eps'];
 set(gcf, 'PaperPositionMode', 'auto')
 print('-depsc',namefile);

 

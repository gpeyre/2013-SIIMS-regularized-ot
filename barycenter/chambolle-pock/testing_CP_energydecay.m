load 'data/inspiringpoints.mat';
x = x(:,:,1:2);
[N d R]=size(x);

nnx = min(4,size(x,1));
optionsG.nnx = nnx;
optionsG.weightedG = false;
optionsG.plot=true;
optionsG.epsilon=0.31;
Gs = zeros(N*nnx,N,R);
h=figure;
for r=1:R
    Gs(:,:,r) = computeKnnGraph(x(:,:,r),optionsG); 
end

options.lambda=1;
options.ksum = 2;
options.ksum_min = 0.1;
Xb1=[];S11=[];
Xb2=[];S12=[];


options.Xb = x(:,:,2);
options.rho = [0.6 0.4];

optionsQP.k_min = 0;
optionsQP.qpprog_tol=getoptions(options,'qpprog_tol',1e-8);
optionsQP.qpprog_niter=getoptions(options,'qpprog_niter',10000);

for r=1:R
    [Sigmas(:,:,r),~,err]=computeSigma_RelaxRegOT_L1(x(:,:,r),options.Xb,Gs(:,:,r),options.ksum,options.lambda,optionsQP);
end 

close all;
%%Testing behaviour of admm w.r.t iterations

options.verb = true;
options.niter = 1e4;
options.sigma = 1e5;
[optimalY,optimalE] =  ComputeY_L1(Gs,Sigmas,x,options.Xb,options.lambda,options.rho, options);

Estar = optimalE(end);

options.niter = 1200;
i=1;
Ys = [];
Es = [];
color=[1 0 0; 0 1 0; 0 0 1; 1 0 1];
f1=figure;

for sigma=[1 100 1e4 1e7]
    options.sigma = sigma;
    [Y,E] = ComputeY_L1(Gs,Sigmas,x,options.Xb,options.lambda,options.rho, options);
   
    E(1)
    figure(f1);
    subplot(1,2,1);
    h(i)=plot(E,'Color',color(i,:));hold on;
    
    subplot(1,2,2);
    h2(i)=plot( real(log( (E(:)-Estar)/Estar )), 'Color', color(i,:) );hold on;
    M{i}=['sigma=' num2str(sigma)];
    
    i=i+1;
    
end

legend(h,M)
subplot(1,2,1);
title('Energy');
subplot(1,2,2);title('Relative error')
legend(h2,M);


f2=figure;i=1;
for sigma=[1 100 1e4 1e7]

    options.sigma = sigma;
    [Y,E] = ComputeY_L1(Gs,Sigmas,x,options.Xb,options.lambda,options.rho, options);
   
    h3(i)=plot( real(log( (E(:)-Estar).^2/(Estar^2) )), 'Color', color(i,:) );hold on;
    
    i=i+1;
end   
    
figure(f2);
title('relative L2 error');
legend(h3,M)
    
    
    
    
    
    
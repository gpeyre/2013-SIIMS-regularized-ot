function [Xb,Sigmas,err]=computeBarycenterRegularizedDiag(Xs,Gs,rho,options)
%% This function computes the barycenter: 
% Xb = argmin sum_r rho_r E(mu_r, mu_X)
% where
% E(mu_X,mu_Y) = min_Sigma <C_{X,Y}, Sigma> + lamdba | G_x(diag(Sigma 1)X-Sigma Y) |^p

% where the parameters have the following dimensions: 
% for N is the number of pixels, d channels, R input datasets
% Xs: (N,d,R) 
% Gs: (P,P,R), for P=Nxnn, nn-nearest neighbourhood
% rho:(R,1)
% lambda: 1x1
% p can be either 2 or 1

[N,d,R] = size(Xs);
[P,~,~] = size(Gs);
lambda = getoptions(options,'lambda',1);
T = getoptions(options,'T',N);
ksum = getoptions(options,'ksum',1.5);
optionsQP.ksum_min=getoptions(options,'ksum_min',0.1);
niter = getoptions(options,'niter',50);
threshold= getoptions(options,'threshold',1e-5);
Xb = getoptions(options,'Xb',Xs(:,:,find(rho==max(rho),1)));
Sigmas = getoptions(options,'Sigmas',repmat(eye(N,N),[1 1 R]));
qp = getoptions(options,'qp',true);
verb = getoptions(options,'verbose',false);
%auxiliary functions
C = @(x,y)(repmat( sum(x.^2,2), [1 N] ) + repmat( sum(y.^2,2)' , [N 1] )- 2*x*y');

optionsQP.qpprog_tol=getoptions(options,'qpprog_tol',1e-8);
optionsQP.qpprog_niter=getoptions(options,'qpprog_niter',10000);

timeSigma = [];
timeY = [];
En=Energy(Xb,Sigmas,rho,Xs,Gs,lambda,C,qp);

%%debugging
ms = 8.5; lw =1.5;
myplot3 = @(X,col)plot3(X(:,1),X(:,2),X(:,3), 'o', 'MarkerSize', ms, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col, 'LineWidth', lw);

if (En < 1e-8) 
  %%debuging
  if verb
    clf;
    plotResults(Xs(:,:,1),Xb,Sigmas(:,:,1),options)
    myplot3(Xs(:,:,2),'g')
    view([90 0])
    title(['lambda=' num2str(lambda) ' rho=' num2str(rho(1)) ' it=' num2str(0) ' Energy=' num2str(En(end))])
    drawnow;
%     end debuging
  end
    return;
end

for i=1:niter
    %% Obtain Sigma_r such that X = Sigma_r X_r
    for r = 1:R
        optionsQP.Cost = [];
        if qp
            [Sigmas(:,:,r)]=computeSigma_RelaxRegOT_qp2(Xs(:,:,r),Xb,Gs(:,:,r),ksum,lambda,T,optionsQP);
        else
            s=tic;
            [Sigmas(:,:,r),~,err]=computeSigma_RelaxRegOT_L1(Xs(:,:,r),Xb,Gs(:,:,r),ksum,lambda,optionsQP);
            timeSigma(end+1)=toc(s);
        end 
    end

    if verb
        En(end+1) = Energy(Xb,Sigmas,rho,Xs,Gs,lambda,C,qp);
    end 
        %% Obtain the associated Xb with L2
    if qp 
       aa = gradH(Xb,Sigmas,Gs,lambda,Xs,rho);
       if sum(abs(aa))>eps
            H= @(X)Energy(X,Sigmas,rho,Xs,Gs,lambda,C,qp);
            optionsBFGS.iterations=1000;
            [Xb,~,info] = perform_bfgs(@(X)deal(H(X),gradH(X,Sigmas,Gs,lambda,Xs,rho)),Xb(:),optionsBFGS);
            Xb = reshape(Xb,N,d);
       end
    else
          
        disp(['computeBarycenterRegularizedDiag it=' num2str(i)]);
       optionsY_L1.threshold_tau = -1;%optionsY_L1.threshold_tau / i
    
       
       Xb = ComputeY_L1(Gs,Sigmas,Xs,Xb,lambda,rho,optionsY_L1);
     
    end 
       
    if (Energy(Xb,Sigmas,rho,Xs,Gs,lambda,C,qp)<1e-8)
        disp(['Reached minimum in ' num2str(i) ' iterations.'])
        return;
    end 
    if verb
        En(end+1) = Energy(Xb,Sigmas,rho,Xs,Gs,lambda,C,qp);
        [i En(end)]
        
          %%debuging
        clf;
        plotResults(Xs(:,:,1),Xb,Sigmas(:,:,1),options)
        myplot3(Xs(:,:,2),'g')
        view([90 0])
        title(['lambda=' num2str(lambda) ' rho=' num2str(rho(1)) ' it=' num2str(i) ' Energy=' num2str(En(end))])
        drawnow;
        %end debuging
    
        
        if (En(end-1)-En(end)) < 1e-8

            if (En(end) > En(end-1))
                disp(['error: ' num2str(En(end)-En(end-1))]);
            end
            return
        end 
    end 
%     
%     if (i>1) && (abs(En(i)-En(i-1)) < threshold)
%         disp(['Reached threshold after ' num2str(i) ' iterations when qp=' num2str(qp)]);
%         err = true;
%         return;
%     end 
%     if (i > 1) && (En(i)>En(i-1)) 
%         disp(['Energy augmenting after ' num2str(i) ' iterations when qp=' num2str(qp)]);
%         err = true;
%         return;
%     end 
%     %end Debug
end

if verb
    close all;
    subplot(2,1,1);
    endE=length(En);
    plot(3:2:endE,En(3:2:end),'*r');hold on;
    plot(4:2:endE,En(4:2:end),'ob');
    plot(En,'g');
    title('Energy decay');
    
    subplot(2,1,2);
    En = En(end-10:end);
    endE=length(En);
    plot(1:2:endE,En(1:2:end),'*r');hold on;
    plot(2:2:endE,En(2:2:end),'ob');
    legend('Sigma it','Y it')
    plot(En,'g');
    title('Energy detailed(last 10 it)');
    

end 
%disp(['End: ' num2str(i) ' iterations']);
 err.val = false;
 err.timeSigma = timeSigma;
 err.timeY = timeY;
% close all


function E=Energy(X,Sigmas,rho,Xs,Gs,lambda,Cost,L2)
E = 0;
prodM = @(A,B)sum(A(:).*B(:));

D =@(S)diag(sum(S,2));
X = reshape(X,size(Xs(:,:,1)));
if L2
    norm2 = @(A)sum(A(:).^2);
    for r=1:length(rho) 
        E = E + rho(r)*0.5*( ...
             prodM(Cost(Xs(:,:,r),X), Sigmas(:,:,r)) + ...
             lambda*norm2(Gs(:,:,r)* (D(Sigmas(:,:,r))*Xs(:,:,r)-Sigmas(:,:,r)*X)) );
    end 
else
    norm1 = @(A)sum(abs(A(:)));
    for r=1:length(rho) 
        E = E + rho(r)*0.5*( ...
             prodM(Cost(Xs(:,:,r),X), Sigmas(:,:,r)) + ...
             lambda*norm1(Gs(:,:,r)* (D(Sigmas(:,:,r))*Xs(:,:,r)-Sigmas(:,:,r)*X)) );
    end 
    
end 
function gH = gradH(X,Sigmas,Gs,lambda,Xs,rho)

[N,d,R]=size(Xs);
gH = zeros(N,d);
X = reshape(X,N,d);
dC = @(x,y)repmat(reshape(x,[N 1 d]), [1 N]) - ...
           repmat(reshape(y,[1 N d]), [N 1] );%dC = @(x,y,Sigma)Sigma'*(x-Sigma*y);  
resh = @(M)reshape(sum(M,1),[N d]);
gradf=@(X,Y,Sigma)-resh( dC(X,Y).* repmat(Sigma,[1 1 d]));
D =@(S)diag(sum(S,2));
A =@(G,S,X,Y)G*(D(S)*X-S*Y);
for r =1:R
    gH = gH + rho(r)*( gradf(Xs(:,:,r),X,Sigmas(:,:,r))  ...
        - lambda * ( Sigmas(:,:,r)'*Gs(:,:,r)'*A(Gs(:,:,r),Sigmas(:,:,r),Xs(:,:,r),X)));
end 

gH = gH(:);

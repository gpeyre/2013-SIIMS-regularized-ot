function [Y,E] = ComputeY_L1(Gs,Sigmas,Xs,Xb,lambda,rho, options)

    [N,d,R] = size(Xs);
    [P,~,~] = size(Gs);
    
    cell2mat3D = @(A)reshape(cell2mat(A),size(A{1},1),size(A{1},2),length(A));
    matrixMult=@(Gs,Sigmas)cell2mat3D(arrayfun(@(ind)Gs(:,:,ind)*Sigmas(:,:,ind), 1:size(Sigmas,3),'uniformOutput',false));
  
    Bs = matrixMult(Gs,Sigmas);
    
    matrixConj =@(B)cell2mat3D(arrayfun(@(ind)B(:,:,ind)',1:size(B,3),'uniformOutput',false));
    matrixDiag =@(v)cell2mat3D(arrayfun(@(ind)diag(v(:,:,ind)),1:size(v,3),'uniformOutput',false));

    As = matrixMult(matrixMult(Gs,matrixDiag(sum(Sigmas,2))),Xs);

    K = @(Y)matrixMult(Bs,repmat(Y,[1 1 size(Bs,3)]));
    Ks= @(U)sum(matrixMult(matrixConj(Bs),U),3);
    
    rS = repmat(reshape(rho,[1 1 R]),[N N 1]).*Sigmas;
    sSX = sum(matrixMult(matrixConj(rS),Xs),3); %sum_k rho_k (S_k)' X_k
    flat = @(x)x(:);
    
    proxH = @(Y,tau)(Y+tau*sSX)./repmat(1+tau*flat(sum(sum(rS,1),3)),[1 d]);
    
    Rho = repmat(reshape(rho,[1 1 R]),[P d 1]);
    proxF = @(v,tau) As+perform_soft_thresholding(v-As, tau*lambda*Rho);
    proxFS= @(u,tau) u-tau*proxF(u/tau, 1/tau); 
    
    options.F  = @(U) sum( flat(cell2mat3D(arrayfun(@(r)rho(r)*lambda*abs(As(:,:,r)-U(:,:,r)),1:length(rho),'uniformOutput',false))));
    
    Cost = @(x,y)(repmat( sum(x.^2,2), [1 N] ) + repmat( sum(y.^2,2)' , [N 1] )- 2*x*y');
    options.H = @(Y)sum(flat( cell2mat3D(arrayfun(...
                                             @(r)rho(r)*Cost(Xs(:,:,r),Y).*Sigmas(:,:,r),...
                                             1:size(Xs,3),...
                                             'uniformOutput',false))));
 
    options.tau = getoptions(options,'threshold_tau',-1);
    [Y,E] = perform_admm(Xb, K,  Ks, proxFS, proxH, @(x)x, options);
    
end 



function E=Energy(X,Sigmas,rho,Xs,Gs,lambda,Cost)
E = 0;
prodM = @(A,B)sum(A(:).*B(:));
norm1 = @(A)sum(abs(A(:)));
D =@(S)diag(sum(S,2));
X = reshape(X,size(Xs(:,:,1)));
for r=1:length(rho) 
    E = E + rho(r)*( ...
         prodM(Cost(Xs(:,:,r),X), Sigmas(:,:,r)) + ...
         lambda*norm1(Gs(:,:,r)* (D(Sigmas(:,:,r))*Xs(:,:,r)-Sigmas(:,:,r)*X)) );
end 
end 
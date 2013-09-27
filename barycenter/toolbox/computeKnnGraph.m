function Gx = computeKnnGraph(x,optionsGx)

nnx = getoptions(optionsGx,'nnx',2);
weightedG = getoptions(optionsGx,'weightedG',true);
plotting = getoptions(optionsGx,'plot',false);
color = getoptions(optionsGx,'color_plot','b');
name = getoptions(optionsGx,'name','graphofFile');
epsilonDistance = getoptions(optionsGx,'epsilon',inf);

Nx = size(x,1);
Iknn = knnsearch(x,[],nnx);
%Ix = [1:Nx; Iknn(:,1:(nnx-1))'];
Ix = Iknn(:,1:(nnx))';
Jx = repmat(1:Nx, [nnx 1]);

norm=@(x)sqrt(sum(x.^2,2));
 
if weightedG
    Wx = 1./(norm( (x(Ix,:) - x(Jx,:))) + 1e-5); 
else 
    Wx = ones(nnx*Nx,1);
end

if epsilonDistance < inf
    D=(norm( (x(Ix,:) - x(Jx,:))));
    Wx(D>=epsilonDistance)=0;
    
end 
    
Gx = sparse( 1:nnx*Nx, Ix(:), Wx, nnx*Nx, Nx  ) - ...
     sparse( 1:nnx*Nx, Jx(:), Wx, nnx*Nx, Nx  );

if (plotting)
    Wx = reshape(Wx,size(Jx));

    if (size(x,2)>2)

     %DEBUG- Draw the graph
        myplot3 = @(X,col)plot3(X(:,1),X(:,2),X(:,3), 'o', 'MarkerSize', 15, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col, 'LineWidth', 1);
        %clf; hold on;
      %  myplot3(x, color);hold on;
        for i=1:size(Ix,2)
            for j=1:nnx
                 if (abs(Wx(j,i)) > 0) 
                     plot3([x(Jx(1,i),1) x(Ix(j,i),1)],...
                           [x(Jx(1,i),2) x(Ix(j,i),2)],...
                           [x(Jx(1,i),3) x(Ix(j,i),3)],'k', 'LineWidth', 0.5 );hold on;
                 end 
            end 
        end
        for i=1:size(x,1)
            plot3(x(i,1),x(i,2),x(i,3), 'o', 'MarkerSize', 15, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [x(i,1),x(i,2),x(i,3)], 'LineWidth', 1);
        %clf; hold on;
        end 
    else 
        myplot = @(X,col)plot(X(:,1),X(:,2),'o', 'MarkerSize', 15, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col, 'LineWidth', 1);
        %clf; hold on;
        myplot(x, color);hold on;
        for i=1:size(Ix,2)
            for j=1:nnx
                if (abs(Wx(j,i)) > 0) 
                plot([x(Jx(1,i),1) x(Ix(j,i),1)],...
                     [x(Jx(1,i),2) x(Ix(j,i),2)],'k', 'LineWidth', 0.5 );
                end
            end 
        end
    end
    xlabel('R');ylabel('G');zlabel('B');
    view(0,90);axis tight;axis off
    namefile = ['./results/' name 'graphRG.eps'];
    print('-depsc',namefile); 

% end DEBUG
end 
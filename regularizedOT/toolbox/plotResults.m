function h=plotResults(x,y,Sigma,options)

Ix = getoptions(options,'Ix',[]);
Jx = getoptions(options,'Jx',[]);
colormapping = getoptions(options,'colorbar',false);
displayU = getoptions(options,'displayU',false);
displayUX = getoptions(options,'displayUX',false);
h = getoptions(options,'h',round(rand(1)*10));
threshold = getoptions(options,'threshold',1e-7);
lwl = getoptions(options,'lwl',1);
%x = 10*x;
%y = 10*y;

%ms = 8.5; lw = 1.5; data paper
ms = 8.5; lw =1.5;
d=size(x,2);
%figure(h);
if  d< 3
    myplot = @(X,col)plot(X(:,1),X(:,2), 'o', 'MarkerSize', ms, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col, 'LineWidth', lw);
    
    if colormapping
        plotline = @(x,y,color)plot( [x(1) y(1)]', [x(2) y(2)]', 'Color', color,'LineWidth', lwl);
    else
        plotline = @(x,y,pattern)plot([x(1) y(1)]',[x(2) y(2)]', pattern,'LineWidth', lwl);
    end
else 
    myplot = @(X,col)plot3(X(:,1),X(:,2),X(:,3), 'o', 'MarkerSize', ms, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col, 'LineWidth', lw);
    if colormapping
        plotline = @(x,y,color)plot3( [x(1) y(1)]', [x(2) y(2)]',[x(3) y(3)]', 'Color', color,'LineWidth', lwl);
    else
        plotline = @(x,y,pattern)plot3( [x(1) y(1)]', [x(2) y(2)]',[x(3) y(3)]', pattern,'LineWidth', lwl);
    end

end 


myplot(x, 'b');hold on;
myplot(y, 'r');view([0 90]);axis equal;
legend('X','Y','Location','WestOutside');
for i=1:size(x,1)
    for j=1:size(y,1)
        if (Sigma(i,j) > threshold)    
            if colormapping
                color = hsv2rgb([min(max(Sigma(i,j),0),1) 1 1]);
                if ((sum(color>1) > 0) || sum(color<0)>0)
                    disp(['err']);
                else 
                   plotline(x(i,:),y(j,:),color);hold on; 
                end 
            else 
                if abs(round(Sigma(i,j))-Sigma(i,j)) < eps
                    plotline(x(i,:),y(j,:),'m');hold on;% hsv2rgb([Sigma(i,j) 1 1]));
                else 
                    if Sigma(i,j) > 0.2
                        plotline(x(i,:),y(j,:),'--y');hold on;
                    else
                      %  plotline(x(i,:),y(j,:),'--k');hold on;
                    end 
                end 
               
            end 
        end 
    end 
end 


%
%comment this if not 3D and regular plotting
% figure;hold on;
%    val = max(x(:));
%    for i = 1:size(x,1)
%          plot3(x(i,1),x(i,2),x(i,3), 'o', 'MarkerSize', ms, 'MarkerFaceColor', min(max(x(i,:)/val,0),1), 'LineWidth', lw);
%          plot3(y(i,1),y(i,2),y(i,3), 'o', 'MarkerSize', ms, 'MarkerFaceColor', min(max(y(i,:)/val,0),1), 'LineWidth', lw);
%           plot3(y(i,1),y(i,2),y(i,3), 'x', 'MarkerSize', ms, 'MarkerFaceColor', min(max(y(i,:)/val,0),1), 'LineWidth', lw);
%    end 


if (displayU)
     u = Sigma*y;
    
     myplot(u,'g');
end 

axis off;


if (displayUX)
    figure;
     u = Sigma*y;
     myplot(u,'g');hold on;
    for i=1:size(u,1)
        plot3( [x(i,1) u(i,1)]', [x(i,2) u(i,2)]',[x(i,3) u(i,3)]', 'k', 'LineWidth', 1 );
    end 
    myplot(x,'b');
    
    myplot(y,'r');%not in paper
end 

if colormapping
    colormap hsv;

    cbh = colorbar('YGrid','on');

    values = linspace(0,1,11);
    maxY = get(cbh,'YLim');
    set(cbh,'ytick',linspace(maxY(1),maxY(2),11));
    set(cbh,'yticklabel',arrayfun(@(x)num2str(x,'%1.1f'),values,'uni',false));
end 

if ( ~isempty(Ix) ) 
    axis on;
    axis equal;
   for i=1:size(Ix,2)
        for j=2:size(Ix,1)
            if d>2
                 plot3([x(Jx(1,i),1) x(Ix(j,i),1)],...
                       [x(Jx(1,i),2) x(Ix(j,i),2)],...
                       [x(Jx(1,i),3) x(Ix(j,i),3)],'k', 'LineWidth', 0.5 );
            else 
                   plot([x(Jx(1,i),1) x(Ix(j,i),1)],...
                        [x(Jx(1,i),2) x(Ix(j,i),2)],'k', 'LineWidth', 0.5 );
            end 
        end 
   end
   
%    figure;hold on;
%    val = max(x(:));
%    for i = 1:size(x,1)
%          plot3(x(i,1),x(i,2),x(i,3), 'o', 'MarkerSize', ms, 'MarkerFaceColor', min(max(x(i,:)/val,0),1), 'LineWidth', lw);
%          plot3(y(i,1),y(i,2),y(i,3), 'o', 'MarkerSize', ms, 'MarkerFaceColor', min(max(y(i,:)/val,0),1), 'LineWidth', lw);
%           plot3(y(i,1),y(i,2),y(i,3), 'x', 'MarkerSize', ms, 'MarkerFaceColor', min(max(y(i,:)/val,0),1), 'LineWidth', lw);
%    end 

   
end 
 
view(2);axis tight;

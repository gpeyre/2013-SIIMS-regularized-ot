
function plot_graph(x,G)

L = G'*G; % Laplacian
[i,j,~] = find(L);
I = find(i~=j); i = i(I); j = j(I);

myplot3 = @(X,col)plot3(X(:,1),X(:,2),X(:,3), 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col, 'LineWidth', 1);
myplot2 = @(X,col)plot(X(:,1),X(:,2), 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col, 'LineWidth', 1);
        

hold on;
if size(x,2)==3
    % plot in 3-D
    
    myplot3(x, 'b');
    for s=1:length(i)
        plot3([x(i(s),1) x(j(s),1)],...
            [x(i(s),2) x(j(s),2)],...
            [x(i(s),3) x(j(s),3)],'k', 'LineWidth', 0.5 );
    end
    
else
    % plot in 2-D
    myplot2(x, 'b');
    for s=1:length(i)
        plot([x(i(s),1) x(j(s),1)],...
            [x(i(s),2) x(j(s),2)],...
            'k', 'LineWidth', 0.5 );
    end    
end
function plot(X)

[N,~] = size(X);
lw =1;ms = 15;
myplot3 = @(X,col)plot3(X(:,1),X(:,2),X(:,3), 'o', 'MarkerSize', ms, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col, 'LineWidth', lw);

for i=1:N
    myplot3(X(i,:),X(i,:));hold on;
end 

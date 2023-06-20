function plotData(data, labels, dType, ephase)
K =length(unique(labels)); 

figure;

% cmap =colormap('lines');
% cmap = cmap(1:K,:);

newmap(1,:) = [0 1 0];
newmap(2,:) = [1 0 0];
newmap(3,:) = [0 0 0];

%Plot of the data vectors
hold on;
classLabels = unique(labels);
for ii=1:K
    ind = find(classLabels(ii) == labels);
    plot3(data(ind,1), data(ind,2),data(ind,3), 'Color',newmap(ii,:),...
        'LineStyle','none', 'Marker', '.','MarkerSize',14);
end

title([ephase ' set for ' dType ' data.']);
legend('Class 1', 'Class 2', 'Class 3');

hold off;
drawnow;
view(3);
dMin = min(data,[],1);
dMax = max(data,[],1);
axis([dMin(1) dMax(1) dMin(2) dMax(2) dMin(3) dMax(3)])
disp('');


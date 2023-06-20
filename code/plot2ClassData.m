function plot2ClassData(data, labels, dType, ephase)
K =length(unique(labels)); 

cmap =colormap('lines');
% cmap = cmap(1:K,:);

newmap(1,:) = [0 1 0];
newmap(2,:) = [1 0 0];
newmap(3,:) = [0 0 0];

figure;

%Plot of the data vectors
hold on;
classLabels = unique(labels);
for ii=1:K
    ind = find(classLabels(ii) == labels);
    plot(data(ind,1), data(ind,2), 'Color',newmap(ii,:),...
        'LineStyle','none', 'Marker', '.','MarkerSize',14);
end

title([ephase ' set for ' dType ' data.']);
legend('Class 1', 'Class 2', 'Class 3');

hold off;
drawnow;
view(2);
dMin = min(data,[],1);
dMax = max(data,[],1);
axis([dMin(1) dMax(1) dMin(2) dMax(2)])
disp('');

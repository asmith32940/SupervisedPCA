function wOptimal =solvePcaObj(classData,labels)

global VERBOSE;
if(VERBOSE); disp('solvePcaObj'); end

K =length(unique(labels)); %Get the K number of classes

[coeff,~,~] =pca(classData);

wOptimal = coeff(:,1:K);


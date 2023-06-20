function WOptimal = localGraphEmbeddingObj(classTrainData,labels)


global VERBOSE;
if(VERBOSE); disp('localGraphEmbeddingObj'); end

WOptimal =[];
K =length(unique(labels));      %Get the K number of classes

options = [];
options.NeighborMode = 'KNN';
options.k = 5;
options.WeightMode = 'Binary';
options.t = 1;
options.PCARatio = 1;

W = constructW(classTrainData,options);


[eigvector, ~] = LGE(W, [], options, classTrainData);
WOptimal = eigvector(:,1:K);

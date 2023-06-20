function WOptimal = localLinearProjectObj(classTrainData, labels)

global VERBOSE;
if(VERBOSE); disp('localLinearProjectObj'); end

K =length(unique(labels));      %Get the K number of classes

options = [];
options.NeighborMode = 'KNN';  %Supervised,KNN
options.k = 3;
options.WeightMode = 'Binary'; %Binary, HeatKernel
options.t = 1;
W = constructW(classTrainData,options);

% options.WeightMode = 'HeatKernel';
% options.t = 5;

options.PCARatio = 1;
[eigvector, ~] = LPP(W, options, classTrainData);
WOptimal = eigvector(:,1:K);




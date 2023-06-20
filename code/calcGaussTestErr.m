
function [predictionMatrix, targetMatrix] =...
                            calcGaussTestErr(classModels, ...
                                             classDataSets, labels)
global VERBOSE;
if(VERBOSE); disp('calcGaussTestErr'); end

%Gaussian classify (calculate probability of membership)
classLabels = unique(labels);
nModels =numel(classModels);
yProb = zeros(size(classDataSets,1), nModels);

%Estimate the probability for membership
for tt=1:nModels
    dMu = classModels(tt).dMu;
    dCov = classModels(tt).dCov;
    
    yProb(:,tt) = mvnpdf(classDataSets, dMu, dCov);
end
    
predictionMatrix = zeros(size(classDataSets,1), nModels);
targetMatrix = zeros(size(classDataSets,1), nModels);

for ii=1:size(classDataSets,1)
    [~,idx] = max(yProb(ii,:),[],2);  %Prediction index location
    
    predictionMatrix(ii,idx) =1;
    
    ind = find(classLabels == labels(ii));
    targetMatrix(ii, ind) =1;
end    


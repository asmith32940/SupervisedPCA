function [classModels, predictionMatrix, targetMatrix] =...
                            calcGaussTrainErr(classDataSets, labels)
global VERBOSE;
if(VERBOSE); disp('calcGaussTrainErr'); end


classLabels = unique(labels);
K =length(unique(labels));            %Get the K number of classes
classModels =struct([]);

%Generate class probability models
for kk=1:K
    ind =find(labels == classLabels(kk));
    data =classDataSets(ind,:);
    
    classModels(kk).cid   =kk;
    classModels(kk).N     =length(ind);    
    classModels(kk).dMu   =mean(data);
    classModels(kk).dCov  =cov(data);   
end

%Gaussian classify (calculate probability of membership)
nModels =numel(classModels);

yProb = zeros(size(classDataSets,1), nModels);

%Estimate the probability for membership
for tt=1:nModels
    dMu = classModels(tt).dMu;
    dCov = classModels(tt).dCov;
    
%     if(is_symmetric_matrix(dCov)), disp('Is Symmetric'); end
%     if(is_positive_semi_matrix(dCov)), disp('Is positive semi'); end
    
%      yProb(:,tt)=gauss(dMu, dCov, classDataSets);
    yProb(:,tt) = mvnpdf(classDataSets, dMu, dCov);
end
    
predictionMatrix = zeros(size(classDataSets,1), nModels);
targetMatrix = zeros(size(classDataSets,1), nModels);

for ii=1:size(classDataSets,1)
    [~,idx] = max(yProb(ii,:),[],2);  %Prediction indx location
    
    predictionMatrix(ii,idx) =1;
    
    ind = find(classLabels == labels(ii));
    targetMatrix(ii, ind) =1;
end    


function [svmMdl, predictionMatrix] =classifyMatlabSVMTrain(classDataSets, labels)

%Generate SVM models
t = templateSVM('Standardize',1,'KernelFunction','linear', 'Verbose', 0);
svmMdl = fitcecoc(classDataSets,labels,'Learners',t,'FitPosterior',0,'Verbose',0);
predictionLabels= predict(svmMdl,classDataSets);

predictionMatrix = zeros(size(classDataSets,1),1);

for ii=1:size(classDataSets,1)     
    %Update the prediction matrix
    predictionMatrix(ii) =predictionLabels(ii);

end    

disp('');
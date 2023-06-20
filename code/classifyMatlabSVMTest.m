function predictionMatrix = ...
                classifyMatlabSVMTest(svmMdl,classDataSets)

predictionLabels= predict(svmMdl,classDataSets);
predictionMatrix = zeros(size(classDataSets,1),1);

for ii=1:size(classDataSets,1)
    
    predictionMatrix(ii) =predictionLabels(ii);    
end    


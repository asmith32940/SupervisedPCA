function displayTestResults(dataSetStats)
N= numel(dataSetStats);


for ii=1:N
    disp(['      Data: ' dataSetStats(ii).dType ]);
%     disp(['ID''s: ' dataSetStats(ii).classIds]);
    disp(['Dimensions: ' num2str(dataSetStats(ii).nDimensions)]);
    disp(['   Classes: ' num2str(dataSetStats(ii).K)]);
    disp(['      Mods: ' dataSetStats(ii).modType]);
    
    
    nObj = numel(dataSetStats(ii).objFun);
    obj = dataSetStats(ii).objFun;
    for jj=1:nObj
        disp(['     Objective Func: ' obj(jj).name]);
        disp(['             Epochs: ' num2str(obj(jj).epochs)]);
        disp(['       Train Accuracy: ' num2str(obj(jj).completeTrainAccuracy)])
        disp(['Train Mean Accuracy: ' num2str(mean(obj(jj).completeTrainAccuracy))]);
        disp(' ')        
        disp(['        Validation Accuracy: ' num2str(obj(jj).completeValidationAccuracy)])
        disp([' Validation Mean Accuracy: ' num2str(mean(obj(jj).completeValidationAccuracy))]);
        disp(' ')        
        disp(['        Test Accuracy: ' num2str(obj(jj).completeTestAccuracy)])
        disp([' Test Mean Accuracy: ' num2str(mean(obj(jj).completeTestAccuracy))]);
        disp(' ');
    end
    disp(' ');
end
disp('Complete');
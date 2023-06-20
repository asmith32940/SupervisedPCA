function displayKernelTestResults(dataSetStats)
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
        disp(['              Sigma: ' num2str(obj(jj).sigma)]);
        disp(['             Trails: ' num2str(obj(jj).trials)]);
        disp(['       Train Values: ' num2str(obj(jj).completeTrainAccuracy)])
        disp(['Train Mean Accuracy: ' num2str(mean(obj(jj).completeTrainAccuracy))]);
        disp(' ')        
        disp(['        Test Values: ' num2str(obj(jj).completeTestAccuracy)])
        disp([' Test Mean Accuracy: ' num2str(mean(obj(jj).completeTestAccuracy))]);
        disp(' ');
    end
    disp(' ');
end
disp('Complete');
fname  ='TestResults_28-Jan-2016-22-53-23.mat';
load(fname);
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
        disp(['Objective Func: ' obj(jj).name]);
        disp(['        Trails: ' num2str(obj(jj).trials)]);
        disp(['    Train Mean: ' num2str(mean(obj(jj).completeTrainAccuracy))]);
        disp(['     Test Mean: ' num2str(mean(obj(jj).completeTestAccuracy))]);
        disp(['    Test Sigma: ' num2str(obj(jj).sigma)]);
        disp(['     Test Type: ' num2str(obj(jj).type)]);
        disp(' ');
    end
    disp(' ');
end
disp('Complete');
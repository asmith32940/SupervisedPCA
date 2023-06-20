%*****************************************************************
% Author: Oneil Smith
% Filename: main_iter_bola_NClass_gaussian_est.m
% Descriptions:
% Limitation:  K < D  The number of classes must be lest than
%           the number of dimensions representing a class.
%*****************************************************************
clc
clear;
close all

global VERBOSE;
VERBOSE     =1;

% global NORM_MIN;
% NORM_MIN =10e-8;

global PLOT;  %Plot object function/conditions
PLOT =1;

global THREE_D;

global MAX_ITR;
MAX_ITR=5;

% global KERNEL_TYPE;
% KERNEL_TYPE ='rbf';

global MAT_SVM;
MAT_SVM =0;

% oFun ={'K-CQS','K-CAS','K-MC-FLD' 'K-SVM', 'K-SPLINE', 'K-PCA'};
% oFun ={'K-CQS','K-CAS','K-MC-FLD', 'K-PCA'};
% oFun ={'K-MC-FLD', 'K-PCA'};
oFun ={'K-CQS'};

% 0: (K-CQS)    Kernel Bola Quadratic
% 1: (K-CAS)    Kernel Bola Absolute Value
% 2: (K-MC-FLD) Kernel Multiclass Fisher linear discriminant analysis
% 3: (K-SVM)    Kernel SVM classification
% 4: (K-SPLINE) Kernel spline classification

% dataItems ={'gen_gauss', 'iris',...
%             'letter', 'libras', 'mpeg7', 'optdigits', 'penbased','poker',...
%             'red_qual', 'satimage',  'seeds', 'segment', 'shuttle',...
%             'texture', 'vehicle', 'vertebral', 'vowel', 'wavelet',...
%             'white_qual', 'wine_rec','vowel','abalone','ecoli','glass','poker2',...
%             'zoo','thyroid','wine'};

% dataItems ={'gen_gauss', 'letter', 'libras', 'shuttle', 'texture', 'iris', ...
%             'seeds', 'wine_rec', 'vowel', 'optdigits', 'penbased','satimage', ...
%             'segment', 'vertebral', 'wavelet', 'white_qual', 'red_qual',     ...
%             'vehicle'};

% dataItems ={'abalone','ecoli','thyroid','vehicle',...
%             'satimage','vertebral','iris','seeds', 'wine_rec', 'vowel',...
%             'penbased','optdigits', 'segment','texture','wine'};

%These are small data sets
% dataItems ={'abalone','ecoli','thyroid','vehicle','satimage',...
%             'vertebral','wine', 'iris','seeds', 'wine_rec', 'vowel',};

%These are large data sets
% dataItems ={'penbased','optdigits', 'segment','texture'};

%These are data sets with 3 classes
% dataItems ={'vertebral', 'thyroid', 'wine', 'iris', 'seeds'};

dataItems ={'seeds'};


%Initialize kernel parameters
%RBF
% sigma =.2:.2:1;
sigma = .48;

TRIALS=1;
for pp=1:numel(dataItems)
    
    dType =dataItems{pp};
    [classData,labels] =loadKernelData(dType);
%     [classData,labels] =loadData(dType);
    
    %     type ='NONE';
    %     type ='D_NORM';
    %     type ='Z_SCORE';
    type ='Z_SCORE';
    classData =performDataModification(classData, type);
    
    classIds      =unique(labels);
    nDimensions   =size(classData,2);          %Dimensions of the data
    K             =length(unique(labels));      %Number of classes
    
    if(VERBOSE)
        disp(['The following experiment will be on ' num2str(TRIALS) ' trials']);
        disp(['Data set: ' dType]);
        disp(['N-Dimensions: ' num2str(nDimensions)]);
        disp(['Number of (K) classes: ' num2str(K)]);
        
        for kk=1:K
            ind = find(labels == classIds(kk));
            disp(['     Class(' num2str(kk) '):  ' num2str(numel(ind))]);
            disp('');
        end
        
        disp(' ');disp(' ');disp(' ');drawnow;
    end
    
    % Algorithm Limitation:  K < D  The number of classes must be less than
    % the number of dimensions representing a class.
    if(K > nDimensions)
        error('Program configuration error');
    end
    
    maximumMean = -999;
    for bb=1:numel(sigma)
        for tt=1:TRIALS
            [trainSet, testSet] = genRandomSubSets(labels);
            
            if(VERBOSE)
                txt =['TRIAL:  ' num2str(tt) ' of ' num2str(TRIALS)];
                disp(txt);
                drawnow;
            end
            
            
            %% Separate training and testing data
            classTrainData      =[];
            classTrainLabels    =[];
            
            classTestData       =[];
            classTestLabels     =[];
            
            for kk=1:K
                aind = trainSet(kk).train;
                classTrainData     =[classTrainData; classData(aind,:)];
                classTrainLabels   =[classTrainLabels; labels(aind)];
                
                bind = testSet(kk).test;
                classTestData      =[classTestData; classData(bind,:)];
                classTestLabels    =[classTestLabels; labels(bind)];
                
            end
            
            kernelParams.sigma = sigma(bb);
            kernelParams.type = 'rbf'; %rbf, poly
            
            % Generate the K Gram matrix.
            %Project to K category sub-space to test classification
            kClassTrainData = projectToKernelSpace([classTrainData;classTestData] , kernelParams);
            kClassTestData = projectDataToKernelSpace(classTestData, classTrainData, kernelParams);
            
            %Initialization of orthogonal basis.  Generate a random random matrix,
            %then take the exponential of the matrix
            A = orth(randn(size(kClassTrainData,2)));
            w0 = A(:,1:K);
            
            for zz=1:numel(oFun)
                disp(['Objective Function:  ' oFun{zz}])
                disp(['Sigma: ' num2str(sigma(bb))]);
                drawnow;
               %% Solve optimization for optimal weight matrix
                
                if strcmp(oFun(zz), 'K-CQS')
                    alphas =solveBolaObjQuadKernel(w0, kClassTrainData, [classTrainLabels;classTestLabels]);
                    txt ='CQS Projected Space';
                    THREE_D =1;
                elseif strcmp(oFun(zz), 'K-CAS')
                    alphas =solveBolaObjAbsoluteKernel(w0, kClassTrainData, classTrainLabels);
                    THREE_D =1;
                elseif strcmp(oFun(zz), 'K-MC-FLD')
                    alphas =solveFisherObjKernel(kClassTrainData, classTrainLabels);
                    txt ='MC-FLD Projected Space';
                    THREE_D =0;
                elseif strcmp(oFun(zz), 'K-PCA')
                    alphas =solvePcaObjKernel(kClassTrainData, classTrainLabels);
                    txt ='PCA Projected Space';
                    THREE_D =1;
                end%End if strcmp(oFun(zz), 'K-CQS')
                
                %% Project to K category sub-space to test classification
                projClassTrainSet   =projectData(alphas, kClassTrainData);
               % projClassTestSet    =projectData(alphas, kClassTestData);
                
                %Check if they go to the origin.
                
                if(PLOT)
                    if(THREE_D)
                        plotData(projClassTrainSet, [classTrainLabels;classTestLabels], dType, txt);
                    else
                        plot2ClassData(projClassTrainSet, [classTrainLabels;classTestLabels], dType, txt);
                    end
                end
                
%                 if(MAT_SVM)
%                     [classModels, trnPredictionMatrix] =...
%                         classifyMatlabSVMTrain(projClassTrainSet,classTrainLabels);
%                 else
%                     trnPredictionMatrix =...
%                         axisAngleProjTrain(alphas, projClassTrainSet,classTrainLabels);
%                 end
%                 
%                 trainConf_Mat = confusionmat(classTrainLabels,trnPredictionMatrix);
%                 if(PLOT)
%                     disp('Train Set Confusion Matrix');
%                     disp(trainConf_Mat)
%                     figure;heatmap(trainConf_Mat, classTrainLabels, classTrainLabels, 1,'Colormap','red','ShowAllTicks',1,'UseLogColorMap',true,'Colorbar',true);
%                     title('Train Set Confusion Matrix');
%                     drawnow;
%                 end
%                 
%                 trainSamples    = sum(sum(trainConf_Mat));
%                 trainCorrect    = sum(diag(trainConf_Mat));
%                 trainIncorrect  = trainSamples-trainCorrect;
%                 trnAccuracy     = (1-(trainIncorrect/trainSamples))*100;
%                 
%                 if(MAT_SVM)
%                     tstPredictionMatrix = classifyMatlabSVMTest(classModels,projClassTestSet);
%                 else
%                     tstPredictionMatrix =...
%                         axisAngleProjTrain(alphas, projClassTestSet, classTestLabels);
%                 end
%                
%                 testConf_Mat = confusionmat(classTestLabels,tstPredictionMatrix);
%                 if(PLOT)
%                     disp('Test Set Confusion Matrix');
%                     disp(testConf_Mat)
%                     figure;heatmap(testConf_Mat, classTestLabels, classTestLabels, 1,'Colormap','red','ShowAllTicks',1,'UseLogColorMap',true,'Colorbar',true);
%                     title('Test Set Confusion Matrix');
%                     drawnow;
%                 end
%                 
%                 if(PLOT)
%                     if(THREE_D)
%                         plotData(projClassTestSet, classTestLabels, dType, 'Test');
%                         plotData([projClassTrainSet; projClassTestSet], [classTrainLabels; classTestLabels], dType, 'All');
%                     else
%                         plot2ClassData(projClassTestSet, classTestLabels, dType, 'Train');
%                         plot2ClassData([projClassTrainSet; projClassTestSet], [classTrainLabels; classTestLabels], dType, 'All');
%                     end  
%                 end
%                 
%                 testSamples    = sum(sum(testConf_Mat));
%                 testCorrect    = sum(diag(testConf_Mat));
%                 testIncorrect  = testSamples-testCorrect;
%                 tstAccuracy    = (1-(testIncorrect/testSamples))*100;
%                                
%                 
%                 %% Populate accuracy vectors
%                 %             trainAccuracy(tt) =trnAccuracy;
%                 %             testAccuracy(tt)  =tstAccuracy;
%                 objFunStats(zz).name                        =oFun{zz};
%                 objFunStats(zz).trials                      =TRIALS;
%                 objFunStats(zz).completeTrainAccuracy(tt)   =trnAccuracy;
%                 objFunStats(zz).completeTestAccuracy(tt)    =tstAccuracy;
%                 objFunStats(zz).sigma                       =sigma(bb);
                
            end%End zz=1:numel(oFun)
            
        end %for tt=1:TRIALS
        
%         currentMean =0;
%         for zz=1:numel(oFun)
%             currentMean = currentMean + mean(objFunStats(zz).completeTestAccuracy);
%         end
%         
%         if(currentMean > maximumMean)
%             dataSetStats(pp).objFun         =objFunStats;
%             dataSetStats(pp).dType          =dType;
%             dataSetStats(pp).classIds       =classIds;
%             dataSetStats(pp).nDimensions    =nDimensions;
%             dataSetStats(pp).K              =K;
%             dataSetStats(pp).modType        =type;
%             
%             maximumMean = currentMean;
%         end
        
        disp(' '); disp(' ');
    end %End for bb=1:numel(sigma)
    
end %for pp=1:numel(dataItems)

% displayKernelTestResults(dataSetStats);
% 
% fname='dataSetStats-05-June-ZA.mat';
% save(fname,'dataSetStats');
% % sendolmail('anthonysmith@fit.edu','Test Complete',fname);


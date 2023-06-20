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

global NORM_MIN;
NORM_MIN =10e-4;

global PLOT;  %Plot object function/conditions
PLOT =0;

global MAX_ITR;
MAX_ITR=500;

global GAUSS_ERROR;
GAUSS_ERROR =1;

global LIB_SVM;
LIB_SVM =0;

global MAT_SVM;
MAT_SVM =0;


% oFun ={'CQS', 'CAS', 'LS_LDA', 'MC-FLD', 'SPLINE', 'SVM', 'PCA'};
% oFun ={'CQS', 'LS_LDA','MC-FLD','PCA'};
oFun ={'CAS'};

% 0: (ACQS)
% 1: (CQS)      Bola Quadratic
% 2: (CAS)      Bola Absolute Value (w/mean)
% 3: (LS_LDA)   Least squares linear discriminant analysis
% 4: (MC-FLD)   Multiclass Fisher linear discriminant analysis
% 5: (SPLINE)   Spline classification
% 6: (SVM)      SVM classification
% 7: (PCA)      Principal component analysis
% 8: (CANON)    Canonical variants

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

dataItems ={'texture','iris','seeds', 'wine_rec', 'vowel', 'optdigits', ...
            'penbased','satimage', 'segment', 'vertebral','vehicle',...
            'abalone','ecoli', 'thyroid','wine'};

% dataItems ={'iris'};
for pp=1:numel(dataItems)
    TRIALS=10;
    
    dType =dataItems{pp};
    [classData,labels] =loadData(dType);
    
    type ='NONE';
%     type ='D_NORM';
%     type ='Z_SCORE';
    classData =performDataModification(classData, type);
    
    classIds      =unique(labels);
    nDimensions   =size(classData,2);           %Dimensions of the data
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
        
        disp('');disp('');disp('');
    end
    
    % Algorithm Limitation:  K < D  The number of classes must be less than
    % the number of dimensions representing a class.
    if(K > nDimensions)
        error('Program configuration error');
    end
    
    for zz=1:numel(oFun)
        disp(['Objective Function:  ' oFun{zz}])
        
        for tt=1:TRIALS
            [trainSet, testSet] = genRandomSubSets(labels);
                                   
            if(VERBOSE)
                txt =['TRIAL:  ' num2str(tt) ' of ' num2str(TRIALS)];
                disp(txt);drawnow;
            end
                        
            %Initialization of orthogonal basis.  Generate a random random matrix,
            %then take the exponential of the matrix
            A = orth(randn(nDimensions));
            w0 = A(:,1:K);
            
            if(VERBOSE)
                %Check initialization matrix of orthogonality
                disp(['Norm of w0 is ' num2str(norm(w0))]);
                if( (norm(w0)-1.0) < 10e-6)
                    disp('Intial System is orthogonal');
                end
            end
            
            
            %% Separate training and testing data          
            classTrainData      =[];
            classTrainLabels    =[];
            
            classTestData       =[];
            classTestLabels     =[];
            
             for kk=1:K
                aind = trainSet(kk).train;
                classTrainData      =[classTrainData; classData(aind,:)];                
                classTrainLabels    =[classTrainLabels; labels(aind)];
                
                bind = testSet(kk).test;
                classTestData      =[classTestData; classData(bind,:)];
                classTestLabels    =[classTestLabels; labels(bind)];
               
            end
                     
            %% Solve optimization for optimal weight matrix          
            
            if strcmp(oFun(zz), 'PCA')
                wOptimal =solvePcaObj(classTrainData, classTrainLabels);
            elseif strcmp(oFun(zz), 'ACQS')
                wOptimal =solveFminSumOfQuadsOnStiefel(w0, classTrainData, classTrainLabels);
            elseif strcmp(oFun(zz), 'CQS')
                wOptimal =solveBolaObj(w0, classTrainData, classTrainLabels);
            elseif strcmp(oFun(zz), 'CAS')
                wOptimal =solveBolaObjAbsolute(w0, classTrainData, classTrainLabels);
            elseif strcmp(oFun(zz), 'LS_LDA')
                wOptimal =leastSquaresObj(classTrainData, classTrainLabels);
            elseif strcmp(oFun(zz), 'MC-FLD')
                wOptimal =solveFisherObj(classTrainData, classTrainLabels);
            elseif strcmp(oFun(zz), 'SPLINE')
                wOptimal =solveSplineObj(classTrainData);
            elseif strcmp(oFun(zz), 'CANON')
                wOptimal =solveCanonicalObj(classTrainData, classTrainLabels);
            end%End if(oFun(zz) ==1)
            
            if(VERBOSE)
                disp(['Norm of wOptimal is ' num2str(norm(wOptimal))]);
                if(norm(wOptimal) ==1)
                    disp('Optimal System is orthogonal');
                end
            end%End if(VERBOSE)
            
            %% Project to K category sub-space to test classification
            projClassTrainSet   =projectData(wOptimal, classTrainData);
            projClassTestSet    =projectData(wOptimal, classTestData);
                                             
            
            if(GAUSS_ERROR)
                %Gaussian classification               
                [classModels, trnPredictionMatrix, trnTargetMatrix] =...
                    calcGaussTrainErr(projClassTrainSet, classTrainLabels);
                
                %Confusion Matrix
                txt=['Gaussian Training Accuracy Trial: ' num2str(tt)];
                [~, rate]=conffig(trnPredictionMatrix, trnTargetMatrix, txt);
                trnAccuracy =rate(1);
                
                %Calculate the test accuracty
                classActualLabels = unique(classTrainLabels);
                [tstPredictionMatrix, tstTargetMatrix] =...
                      calcGaussTestErr(classModels, projClassTestSet, classTestLabels);
                  
                %Confusion Matrix
                txt=['Gaussian Testing Accuracy Trial: ' num2str(tt)];
                [~, rate]=conffig(tstPredictionMatrix, tstTargetMatrix, txt);
                tstAccuracy =rate(1);

                drawnow;
            elseif(LIB_SVM)
                %SVM classification
                [trnPredictionMatrix, trnTargetMatrix, svmClassModels] = ...
                    solveSvmObj_NClass(projClassTrainSet);
                
                txt=['libSVM Training Accuracy Trial: ' num2str(tt)];
                [~, rate]=conffig(trnPredictionMatrix, trnTargetMatrix, txt);
                trnAccuracy =rate(1);
                
                [tstPredictionMatrix, tstTargetMatrix] = ...
                    testSvmObj_NClass(projClassTestSet, svmClassModels);
                
                txt=['libSVM Testing Accuracy Trial: ' num2str(tt)];
                [~, rate]=conffig(trnPredictionMatrix, trnTargetMatrix, txt);
                trnAccuracy =rate(1);
                
               drawnow;
            elseif(MAT_SVM)
                [classModels, trnPredictionMatrix, trnTargetMatrix] =...
                    classifyMatlabSVMTrain(projClassTrainSet,classTrainLabels);
                
                txt=['Svm Training Accuracy Trail: ' num2str(tt)];
                [~, rate]=conffig(trnPredictionMatrix, trnTargetMatrix, txt);
                trnAccuracy =rate(1);
                
                [tstPredictionMatrix, tstTargetMatrix] =...
                    classifyMatlabSVMTest(classModels,projClassTestSet, classTestLabels);
                
                %Confusion Matrix
                txt=['Svm Testing Accuracy Trial: ' num2str(tt)];
                [~, rate]=conffig(tstPredictionMatrix, tstTargetMatrix, txt);
                tstAccuracy =rate(1);
              
                drawnow;
           end
            
             
            %% Populate accuracy vectors
            trainAccuracy(tt) =trnAccuracy;
            testAccuracy(tt)  =tstAccuracy;
            
        end%End for tt=1:TRIALS
        objFunStats(zz).name                    =oFun{zz};
        objFunStats(zz).trials                  =TRIALS;
        objFunStats(zz).completeTrainAccuracy   =trainAccuracy;
        objFunStats(zz).completeTestAccuracy    =testAccuracy;
        
    end %for zz=1:ND
    dataSetStats(pp).objFun         =objFunStats;
    dataSetStats(pp).dType          =dType;   
    dataSetStats(pp).classIds       =classIds;
    dataSetStats(pp).nDimensions    =nDimensions;
    dataSetStats(pp).K              =K;
    dataSetStats(pp).modType        =type;

    disp('Epoch')
end %for pp=1:numel(dataItems)

timeString =datestr(clock);
ind = find(isspace(timeString));
timeString(ind)='-';
ind = find(timeString ==':');
timeString(ind)='-';
% 
fname=['TestResults_' timeString '.mat'];
% fname=['TestResults_' oFun{1} '.mat'];
save(fname, 'dataSetStats', '-mat');
% sendolmail('anthonysmith@fit.edu','Test Complete',fname);
disp('Complete Category Dimensionality Reduction Test');

% diary off

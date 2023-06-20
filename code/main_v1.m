%*****************************************************************
% Author: Oneil Smith
% Filename: 
% Descriptions:
% Limitation:  K < D  The number of classes must be lest than
%           the number of dimensions representing a class.
%*****************************************************************
clc
clear;
close all

global VERBOSE;
VERBOSE     =0;

global NORM_MIN;
NORM_MIN =10e-4;

global PLOT;  %Plot object function/conditions
PLOT =0;

global MAX_ITR;
MAX_ITR=500;

global MAT_SVM;
MAT_SVM =1;


% oFun ={'CQS', 'CAS', 'LS_LDA', 'MC-FLD', 'SPLINE', 'SVM', 'PCA'};
% oFun ={'CQS','CAS','LS_LDA', 'MC-FLD', 'PCA','LPP','LGE','OLGE','OLPP','SDR','PP','RPCA'};
%SDR, Projection Pursuit, Robust orthogonal complement PCA
%Linear Graph EmbeddingLinear Graph Embedding(LGE)
%Orthogonal LGE (OLGE)
%Orthogonal Locality Preserving Projections(OLPP)
oFun ={'LPP'};

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

% dataItems ={'abalone','ecoli','penbased','thyroid','optdigits','vehicle',...
%             'satimage','segment','vertebral','wine','texture',...
%             'iris','seeds', 'wine_rec', 'vowel'};

% Journal test sets
% dataItems ={'vehicle','wine','iris','seeds','thyroid','satimage',...
%            'segment','vertebral'};

dataItems ={'iris'};
EPOCH = 5;
for pp=1:numel(dataItems)
    
    dType =dataItems{pp};
    [classData,labels] =loadData(dType);
    
     type ='NONE';
%     type ='D_NORM';
%    type ='Z_SCORE';
    classData =performDataModification(classData, type);
    
    classIds      =unique(labels);
    nDimensions   =size(classData,2);           %Dimensions of the data
    K             =length(unique(labels));      %Number of classes
    
    if(VERBOSE)
        disp(['The following experiment will be on ' num2str(EPOCH) ' epochs']);
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
    
    dataSet = genRandomSubSets(labels);
    
    %% Separate training and testing data
    classTrainData      =[];
    classTrainLabels    =[];
    classValidationData =[];
    
    classTestData       =[];
    classTestLabels     =[];
    classValidationLabels =[];
    
    for kk=1:K
        aind = dataSet(kk).train;
        classTrainData      =[classTrainData; classData(aind,:)];
        classTrainLabels    =[classTrainLabels; labels(aind)];
        
        bind = dataSet(kk).test;
        classTestData      =[classTestData; classData(bind,:)];
        classTestLabels    =[classTestLabels; labels(bind)];
        
        vind = dataSet(kk).valid;
        classValidationData      =[classValidationData; classData(vind,:)];
        classValidationLabels    =[classValidationLabels; labels(vind)];
    end
    
    for zz=1:numel(oFun)
        
        for iter=1:EPOCH
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
                  
            disp(['Objective Function:  ' oFun{zz} '  Epoch: ' num2str(iter) ]  )
            %% Solve optimization for optimal weight matrix
            
            if strcmp(oFun(zz), 'PCA')
                wOptimal =solvePcaObj(classTrainData, classTrainLabels);
                txt ='PCA Projected Space';
            elseif strcmp(oFun(zz), 'ACQS')
                wOptimal =solveFminSumOfQuadsOnStiefel(w0, classTrainData, classTrainLabels);
                txt ='ACQS Projected Space';
            elseif strcmp(oFun(zz), 'CQS')
                wOptimal =solveBolaObjQuad(w0, classTrainData, classTrainLabels);
                txt ='CQS Projected Space';
            elseif strcmp(oFun(zz), 'CAS')
                wOptimal =solveBolaObjAbsolute(w0, classTrainData, classTrainLabels);
                txt ='CAS Projected Space';
            elseif strcmp(oFun(zz), 'LS_LDA')
                wOptimal =leastSquaresObj(classTrainData, classTrainLabels);
                txt ='LS_LDA Projected Space';
            elseif strcmp(oFun(zz), 'MC-FLD')
                wOptimal =solveFisherObj(classTrainData, classTrainLabels);
                txt ='MC-FLD Projected Space';
            elseif strcmp(oFun(zz), 'SPLINE')
                wOptimal =solveSplineObj(classTrainData);
                txt ='SPLINE Projected Space';
            elseif strcmp(oFun(zz), 'CANON')
                wOptimal =solveCanonicalObj(classTrainData, classTrainLabels);
                txt ='CANON Projected Space';
            elseif(strcmp(oFun(zz), 'BIAS'))
                wOptimal = solveBolaBiasFreeObj(w0, classTrainData, classTrainLabels);
                txt ='BIAS Projected Space';
            elseif(strcmp(oFun(zz), 'LPP'))
                wOptimal = localLinearProjectObj(classTrainData, classTrainLabels);
                txt ='LPP-Space';
            elseif(strcmp(oFun(zz), 'LGE'))
                wOptimal = localGraphEmbeddingObj(classTrainData, classTrainLabels);
                txt ='LGE-Projected Space';
            elseif(strcmp(oFun(zz), 'OLGE'))
                %[eigvector, eigvalue, bSuccess] = OLGE(W, D, options, data);
                txt ='Orthogonal LGE Projected Space';
            elseif(strcmp(oFun(zz), 'OLPP'))
                %[eigvector, eigvalue, bSuccess] = OLPP(W, options, data);
                txt ='Orthogonal Locality Preserving Projected Space';
            elseif(strcmp(oFun(zz), 'SDR'))
                txt ='SDR Projected Space';
            elseif(strcmp(oFun(zz), 'PP'))
                
                txt ='Projection Pursuit Projected Space';
            elseif(strcmp(oFun(zz), 'RPCA'))
                %[L, S] = RobustPCA(X, lambda, mu, tol, max_iter);
                txt ='Robust Orthogonal Complement PCA Projected Space';
            elseif(strcmp(oFun(zz),'NPE'))
                %[eigvector, eigvalue] = NPE(options, data);
                txt ='Neighborhood Preserving Embedding Projected Space';                
            end%End if(oFun(zz) ==1)

            if(VERBOSE)
                disp(['Norm of wOptimal is ' num2str(norm(wOptimal))]);
                if(norm(wOptimal) ==1)
                    disp('Optimal System is orthogonal');
                end
            end%End if(VERBOSE)
            
            %% Project to K category sub-space to test classification
            projClassTrainSet         =projectData(wOptimal, classTrainData);
            projClassValidationSet    =projectData(wOptimal, classValidationData);
            
            [classModels, trnPredictionMatrix] =...
                classifyMatlabSVMTrain(projClassTrainSet,classTrainLabels);
            
            trainConfMat = confusionmat(classTrainLabels,trnPredictionMatrix);
            if(VERBOSE)
                disp(trainConfMat)
                figure;heatmap(trainConfMat, classTrainLabels, classTrainLabels, 1,'Colormap','red','ShowAllTicks',1,'UseLogColorMap',true,'Colorbar',true);
                title([txt '  (Train)'])
            end
            
            trainSamples        = sum(sum(trainConfMat));
            trainCorrect        = sum(diag(trainConfMat));
            trainIncorrect      = trainSamples-trainCorrect;
            trnAccuracy(iter)   = (1-(trainIncorrect/trainSamples))*100;
            
            valPredictionMatrix = classifyMatlabSVMTest(classModels,projClassValidationSet);
            
            valConfMat = confusionmat(classValidationLabels,valPredictionMatrix);
            if(VERBOSE)
                disp(valConfMat)
                figure;heatmap(valConfMat, classValidationLabels, classValidationLabels, 1,'Colormap','red','ShowAllTicks',1,'UseLogColorMap',true,'Colorbar',true);
                title([txt '  (Validation)'])
            end
           
           if(VERBOSE)
                plotData(classValidationData, classValidationLabels, dType, [txt '  Validation']);
                drawnow;
           end
            
            valSamples    = sum(sum(valConfMat));
            valCorrect    = sum(diag(valConfMat));
            valIncorrect  = valSamples-valCorrect;
            valAccuracy(iter)   = (1-(valIncorrect/valSamples))*100;
            dataModel(iter).model   = classModels;
            dataModel(iter).wMatrix = wOptimal;
            
        end%for iter=1:EPOCH
        
        %Compute the test accuracy
        [mVal, mInd]        = max(valAccuracy);
        projClassTestSet    = projectData(dataModel(mInd).wMatrix, classTestData);
        tstPredictionMatrix = classifyMatlabSVMTest(dataModel(mInd).model, projClassTestSet);
        tstConfMat          = confusionmat(classTestLabels,tstPredictionMatrix);
        if(VERBOSE)
            disp(tstConfMat)
            figure;heatmap(tstConfMat, classTestLabels, classTestLabels, 1,'Colormap','red','ShowAllTicks',1,'UseLogColorMap',true,'Colorbar',true);
            title([txt '  (Test)'])
        end
        
        tstSamples    = sum(sum(tstConfMat));
        tstCorrect    = sum(diag(tstConfMat));
        tstIncorrect  = tstSamples-tstCorrect;
        tstAccuracy   = (1-(tstIncorrect/tstSamples))*100;
            
        %% Populate accuracy vectors
        objFunStats(zz).name                                    =oFun{zz};
        objFunStats(zz).completeTrainAccuracy       =trnAccuracy;
        objFunStats(zz).completeValidationAccuracy  =valAccuracy;
        objFunStats(zz).dataModel                             =dataModel;
        objFunStats(zz).completeTestAccuracy        =tstAccuracy;
        objFunStats(zz).epochs                                  =EPOCH;
     
    end%End for zz=1:numel(oFun)
    
    
    dataSetStats(pp).objFun         =objFunStats;
    dataSetStats(pp).dType          =dType;
    dataSetStats(pp).classIds       =classIds;
    dataSetStats(pp).nDimensions        =nDimensions;
    dataSetStats(pp).K                          =K;
    dataSetStats(pp).modType             =type;
    
    disp('Epoch')
end %for pp=1:numel(dataItems)

displayTestResults(dataSetStats);

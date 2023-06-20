%*****************************************************************
% Author: Oneil Smith
% Filename: 
% Descriptions:
%
%*****************************************************************
clc
clear;
close all

global dPRINT;
dPRINT =0;

global VERBOSE
VERBOSE =1;

% gp=pwd;
% addpath(genpath(gp(1:strfind(pwd, 'CategoryPerceptron')-1)));

%Data need to be in DxN (i.e. for 3D data 3xN)
dataItems ={'GEN_GAUSS', 'iris','wine','poker1','poker2',...
            'poker3', 'wmpeg7', 'race', 'caltech','seeds'};
dataset ='iris';
[XiTrain, XjTrain, XkTrain, subSetIndices]=loadData3Class(dataset);

% modTypes ={'BIAS', 'UNIT_CIRCLE','NORM_DATA','SUB_MINIMUM','SUB_CENTROID'};
% type = 'NORM_DATA';
% [XiTrain, XjTrain, XkTrain] = ...
%     performDataModification3Class(XiTrain, XjTrain, XkTrain, type);

TRIALS =length(subSetIndices);
trainAccuracy =zeros(1,TRIALS);
testAccuracy  =zeros(1,TRIALS);
TrialResults  =cell(1,TRIALS) ;

for tt=1:TRIALS
    txt =['TRIAL:  ' num2str(tt) ' of ' num2str(TRIALS)];
    disp(txt);drawnow;
    
    %Seperate the training and test data based on a leave
    %one-out-method.  But the one in this case is a subset of the
    %training input vectors.
    %Get X-Y coordinate for training and test
    XiTrainInd  =subSetIndices(tt).XiTrainInd;
    XiTestInd   =subSetIndices(tt).XiTestInd;
    
    XjTrainInd  =subSetIndices(tt).XjTrainInd;
    XjTestInd   =subSetIndices(tt).XjTestInd;
    
    XkTrainInd  =subSetIndices(tt).XkTrainInd;
    XkTestInd   =subSetIndices(tt).XkTestInd;
    
    Xi =XiTrain(:,XiTrainInd);
    Xj =XjTrain(:,XjTrainInd);
    Xk =XkTrain(:,XkTrainInd);
    
    XiTest =XiTrain(:,XiTestInd);
    XjTest =XjTrain(:,XjTestInd);
    XkTest =XkTrain(:,XkTestInd);
                 
    XiValidate =Xi;
    XjValidate =Xj;
    XkValidate =Xk;
   
    KERNEL=0;
    if(KERNEL)        
        [vW1, vW2, vW3] =solveSplineObj_kernel_3Class(Xi, Xj, Xk);
    else
        [vW1, vW2, vW3, fval] =solveSplineObj_linear_3Class(Xi, Xj, Xk);
    end
    
    P =[vW1 vW2 vW3];
       
    if(VERBOSE)
        txt1 ='W1 is a unit vector ';
        disp([txt1 num2str(norm(vW1))]);
        
        txt2 ='W2 is a unit vector ';
        disp([txt2 num2str(norm(vW2))]);
        
        txt3 ='W3 is a unit vector ';
        disp([txt3 num2str(norm(vW3))]);

        txt4 ='W1 & W2 are perpendicular  ';
        disp([txt4 num2str(dot(vW1, vW2))]);
        
        txt5 ='W2 & W3 are perpendicular  ';
        disp([txt5 num2str(dot(vW2,vW3))]);        
    end%End if(VERBOSE)
    
    type ='poly';
    if(KERNEL)        
        Xi = projectToKernelSpace(Xi,[XiValidate, XjValidate, XkValidate], type);
        Xj = projectToKernelSpace(Xj,[XiValidate, XjValidate, XkValidate], type);
        Xk = projectToKernelSpace(Xk,[XiValidate, XjValidate, XkValidate], type);
    end
    
    [Zi, Zj, Zk] =projectToSpline3Class(P, Xi, Xj, Xk);
    
    %Calculate the Gaussian error
    [muZi, covZi, ...
     muZj, covZj, ...
     muZk, covZk, ...
     trnAccuracy, trnTargetMatrix, ...
     trnPredictionMatrix] =...
                        calculateTrainErrorGaussEst_3Class(Zi, Zj, Zk, tt);
    
    if(KERNEL)        
        XiTest = projectToKernelSpace(XiTest,[XiValidate, XjValidate, XkValidate], type);
        XjTest = projectToKernelSpace(XjTest,[XiValidate, XjValidate, XkValidate], type);
        XkTest = projectToKernelSpace(XkTest,[XiValidate, XjValidate, XkValidate], type);
    end
    
    [ZiTest, ZjTest, ZkTest] =projectToSubSpace3Class(P, XiTest, XjTest, XkTest);

    [tstAccuracy, tstTargetMatrix, tstPredictionMatrix] =...
        calculateTestErrorGaussEst_3Class(ZiTest, ZjTest, ZkTest, ...
                 muZi, covZi, muZj, covZj, muZk, covZk, tt);  
      
    trainAccuracy(tt) =trnAccuracy;    
    testAccuracy(tt)  =tstAccuracy;
               
    ResultStruct.Name ='SPLINE';
    ResultStruct.Class =3;
    ResultStruct.Bias =0;
    ResultStruct.Error ='GAUSS';
    ResultStruct.trainAccuracy =trnAccuracy;
    ResultStruct.testAccuracy =tstAccuracy;
        
    ResultStruct.tt =tt;
    ResultStruct.trnTargetMatrix     =uint8(trnTargetMatrix);
    ResultStruct.trnPredictionMatrix =uint8(trnPredictionMatrix);
    ResultStruct.tstTargetMatrix     =uint8(tstTargetMatrix);
    ResultStruct.tstPredictionMatrix =uint8(tstPredictionMatrix);

    TrialResults{tt} =ResultStruct;
        
    disp('');
end%End for tt=1:TRIALS

% %Calculate Training Accuracy & Error
% trueTrainAccuracy = mean(trainAccuracy);
% txt =['Mean Training Accuracy(Gauss) ' num2str(trueTrainAccuracy) '%'];
% disp(txt);
% 
% trainError = repmat(100,1,TRIALS)-trainAccuracy;
% trueTrainError = mean(trainError);
% txt =['Mean Training Error(Gauss) ' num2str(trueTrainError) '%'];
% disp(txt);

%Calculate Testing Accuracy & Error
trueTestAccuracy = mean(testAccuracy);
txt =['Mean Test Accuracy ' num2str(trueTestAccuracy) '%'];
disp(txt);

testError = repmat(100,1,TRIALS)-testAccuracy;
trueTestError = mean(testError);
txt =['Mean Test Error ' num2str(trueTestError) '%'];
disp(txt);

txt =['Test Std  ' num2str(std(testAccuracy))];
disp(txt);

disp('Complete Spline Test');


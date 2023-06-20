%*****************************************************************
% Author: Oneil Smith
% Filename: 
% Descriptions:
%
%*****************************************************************
clc
clear;
close all

dPRINT =0;
VERBOSE=1;
SAVE_RESULTS=0;

gp=pwd;
addpath(genpath(gp(1:strfind(pwd, 'CategoryPerceptron')-1)));

%Data need to be in DxN (i.e. for 3D data 3xN)
dataItems ={'GEN_GAUSS', 'Gauss','iris','wine','poker1','poker2','wmpeg7'};
dataset ='GEN_GAUSS';

[XiTrain, XjTrain, subSetIndices]=loadData2Class(dataset);

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
    
    %Get X-Y coordinate for training and test
    XjTrainInd  =subSetIndices(tt).XjTrainInd;
    XjTestInd   =subSetIndices(tt).XjTestInd;
    
    Xi =XiTrain(:,XiTrainInd);
    Xj =XjTrain(:,XjTrainInd);

    XiTest =XiTrain(:,XiTestInd);
    XjTest =XjTrain(:,XjTestInd);
                 
    XiValidate =Xi;
    XjValidate =Xj;
    XiTestValidate =XiTest;
    XjTestValidate =XjTest;
   
    [vW1, vW2, ~] =solveSplineObj_kernel_2Class(Xi, Xj);
    
%     [vW1, vW2, ~] =solveSplineObj_linear_2Class(Xi, Xj);
    P =[vW1 vW2];
       
    if(VERBOSE)
        txt1 ='W1 is a unit vector ';
        disp([txt1 num2str(norm(vW1))]);
        
        txt2 ='W2 is a unit vector ';
        disp([txt2 num2str(norm(vW2))]);
        
        txt3 ='W1 & W2 are perpendicular  ';
        disp([txt3 num2str(dot(vW1, vW2))]);
        
    end%End if(VERBOSE)
    
    type ='poly';
    Xi = projectToKernelSpace(Xi,[XiValidate, XjValidate], type);
    Xj = projectToKernelSpace(Xj,[XiValidate, XjValidate], type);
    
    [Zi, Zj] =projectToSpline2Class(P, Xi, Xj);
    
    figure;
    plot(XiValidate(1,:), XiValidate(2,:), 'r*');drawnow; hold on;
    plot(XjValidate(1,:), XjValidate(2,:), 'g*');
    plot(Zi(1,:), Zi(2,:), 'ro');hold on;
    plot(Zj(1,:), Zj(2,:), 'go');hold off;drawnow
    
    
    
%     XiTest = projectToKernelSpace(XiTest,[XiValidate, XjValidate], type);
%     XjTest = projectToKernelSpace(XjTest,[XiValidate, XjValidate], type);
    
% %     [trnAccuracy, trnTargetMatrix, trnPredictionMatrix] =calculateErrorAngle_2Class(Zi, Zj, P, tt, dPRINT, VERBOSE);
%     [trnAccuracy, ~, ~] =calculateTrainErrorGaussEst_2Class(Zi, Zj, tt, dPRINT, VERBOSE);
% 
%     [ZiTest, ZjTest] =projectToSubSpace2Class(P, XiTest, XjTest);
% %     [tstAccuracy, tstTargetMatrix, tstPredictionMatrix] =calculateErrorAngle_2Class(ZiTest, ZjTest, P, tt, dPRINT, VERBOSE);
% 
%     muZi  =mean(Zi,2)';
%     covZi =cov(Zi');
%     
%     muZj  =mean(Zj,2)';
%     covZj =cov(Zj');
%     [tstAccuracy, ~, ~] =calculateTestErrorGaussEst_2Class(ZiTest, ZjTest, muZi, covZi,...
%                                       muZj, covZj, tt, dPRINT, VERBOSE);   
%     trainAccuracy(tt) =trnAccuracy;    
%     testAccuracy(tt)  =tstAccuracy;
%                
%     ResultStruct.Name ='Bola';
%     ResultStruct.Class =2;
%     ResultStruct.Bias =0;
%     ResultStruct.Error ='ANGLE';
%         
%     TrialResults{tt} =ResultStruct;
        
    disp('');
end%End for tt=1:TRIALS

if(SAVE_RESULTS)
    c = clock;
    txt =['bola_3class_density_' dataset '_' num2str(c(2)) '_' num2str(c(3)) '_' num2str(c(1))];
    ofname =[txt '_' num2str(c(4)) '_' num2str(c(5)) '_' num2str(round(c(6))) '.mat'];
    
    save(ofname, 'TrialResults');
end


%Calculate Training Accuracy & Error
trueTrainAccuracy = mean(trainAccuracy);
txt =['Mean Training Accuracy(ANGLE) ' num2str(trueTrainAccuracy) '%'];
disp(txt);

trainError = repmat(100,1,TRIALS)-trainAccuracy;
trueTrainError = mean(trainError);
txt =['Mean Training Error(ANGLE) ' num2str(trueTrainError) '%'];
disp(txt);

%Calculate Testing Accuracy & Error
trueTestAccuracy = mean(testAccuracy);
txt =['Mean Test Accuracy ' num2str(trueTestAccuracy) '%'];
disp(txt);

testError = repmat(100,1,TRIALS)-testAccuracy;
trueTestError = mean(testError);
txt =['Mean Test Error ' num2str(trueTestError) '%'];
disp(txt);

disp('Complete Category Perceptron Test');

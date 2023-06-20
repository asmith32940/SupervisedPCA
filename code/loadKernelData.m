function [classData, labels] =loadKernelData(type)

global VERBOSE;
if(VERBOSE); disp('LOAD-KERNEL-DATA'); end

classData ={};

type = upper(type);
if(strcmp(type, 'IDEAL'))
    
    y=[];
    classData =[];
    P = 0.5*ones(3,1);
    N = 10;

    cmap = colormap(lines);hold on;
 
    dCov1 =[.02 0 0;0 .02 0;0 0 .02];
    dMu1  =[0 0 1];
    X= mvnrnd(dMu1,dCov1, fix(P(1)*N));
    classData = [classData; X];
    y = [y ones(1, size(X,1))*1];
    plot3(X(:,1),X(:,2),X(:,3),'LineStyle','none', 'Marker', '+','Color',cmap(1,:))
    
    dCov2 =[.02 0 0;0 .02 0;0 0 .02];
    dMu2  =[0 0 1];
    X= mvnrnd(dMu2, dCov2, fix(P(2)*N));
    classData = [classData; X];
    y = [y ones(1, size(X,1))*2];
    plot3(X(:,1),X(:,2),X(:,3),'LineStyle','none', 'Marker', '+','Color',cmap(2,:))
    
    dCov3 =[.02 0 0;0 .02 0;0 0 .02];
    dMu3  =[0 0 1];
    X= mvnrnd(dMu3, dCov3, fix(P(3)*N));
    classData = [classData; X];
    y = [y ones(1, size(X,1))*3];
    plot3(X(:,1),X(:,2),X(:,3),'LineStyle','none', 'Marker', '+','Color',cmap(3,:))
    
    y =y';
    hold off;
    view(3);
    
elseif(strcmp(type, 'GEN_GAUSS'))
    l = 4; %dimensions
    c = 3; %no. of classes
    N = 10;

    m = randn(l,c);  %generate clas means
    S = zeros(l,l,c);
    
    for ii=1:c
        ts = rand(l,l);
        ts = ts*ts';
        S(:,:,ii) = ts;
    end
    
    P = 0.5*ones(c,1);
    
    [~,c] = size(m);
   
%     figure; hold on;
    
    classData =[];
    y=[];
%     cmap = colormap(prism);
    for ii=1:c
        ts = rand(l,l);
        ts = ts*ts';
       
        %Generating the [p(j)*N)] = vectors from each distribution
        %The total number of points may be slightly less than N due to the fix
        %operator
        X= mvnrnd(m(:,ii), ts, fix(P(ii)*N));
%         plot3(X(:,1),X(:,2),X(:,3),'LineStyle','none', 'Marker', '+','Color',cmap(ii,:))
                
        classData = [classData; X];
        
        y = [y ones(1, fix(P(ii)*N))*ii];
    end
    y=y';
%     hold off;
%     view(3);
       
elseif(strcmp(type, 'SATIMAGE'))
    load('k-SatImageData.mat');
elseif(strcmp(type, 'VEHICLE'))  
    load('k-VehicleData.mat');
elseif(strcmp(type, 'VERTEBRAL'))
    load('k-VertebralData.mat');
elseif(strcmp(type, 'VOWEL'))
    load('k-VowelData.mat');
elseif(strcmp(type, 'ABALONE'))
    load('k-AbaloneData.mat');
elseif(strcmp(type, 'ECOLI'))
    load('k-EcoliData.mat');  
elseif(strcmp(type, 'IRIS'))
    load('k-Iris.mat');
elseif(strcmp(type, 'OPTDIGITS'))
    load('k-OptdigitsData.mat');
elseif(strcmp(type, 'PENBASED'))
    load('k-PenbasedData.mat');
elseif(strcmp(type, 'SEEDS'))
    load('k-SeedsData.mat');
elseif(strcmp(type, 'SEGMENT'))
    load('k-SegmentData.mat');
elseif(strcmp(type, 'TEXTURE'))
    load('k-TextureData.mat');
elseif(strcmp(type, 'WINE_REC'))
    load('k-WineRecognitionData.mat');
elseif(strcmp(type, 'THYROID'))
    load('k-ThyroidData.mat');
elseif(strcmp(type, 'WINE'))
    load('k-WineData.mat');  
end
classData = classData;
labels =y;

disp('COMPLETE DATA LOADING');


function [classData, labels] =loadData(type)

global VERBOSE;
if(VERBOSE); disp('LOAD-DATA'); end

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
    N = 200;

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
       
elseif(strcmp(type, 'CALTECH'))
    load('CalTechFaceGrayData.mat');
elseif(strcmp(type, 'COIL')) 
    load('CoilData.mat');
elseif(strcmp(type, 'IRIS'))
    load('IrisData.mat');
elseif(strcmp(type, 'LETTER'))  
    load('LetterData.mat');
elseif(strcmp(type, 'LIBRAS'))
    load('MovementLibrasData.mat');
elseif(strcmp(type, 'MPEG7'))
    load('Mpeg7DbData.mat');
elseif(strcmp(type, 'OPTDIGITS'))
    load('OptdigitsData.mat');
elseif(strcmp(type, 'PENBASED'))
    load('PenbasedData.mat');
elseif(strcmp(type, 'POKER'))
    load('PokerData.mat');
elseif(strcmp(type, 'RED_QUAL'))
    load('RedQualityData.mat');
elseif(strcmp(type, 'SATIMAGE'))
    load('SatImageData.mat');
elseif(strcmp(type, 'SEEDS'))
    load('SeedsData.mat');
elseif(strcmp(type, 'SEGMENT'))
    load('SegmentData.mat');
elseif(strcmp(type, 'SHUTTLE'))
    load('ShuttleData.mat');
elseif(strcmp(type, 'TEXTURE'))
    load('TextureData.mat');
elseif(strcmp(type, 'VEHICLE'))  
    load('VehicleData.mat');
elseif(strcmp(type, 'VERTEBRAL'))
    load('VertebralData.mat');
elseif(strcmp(type, 'VOWEL'))
    load('VowelData.mat');
elseif(strcmp(type, 'WAVELET'))
    load('WaveletCoefficientData.mat');
elseif(strcmp(type, 'WHITE_QUAL'))
    load('WhiteQualityData.mat');
elseif(strcmp(type, 'WINE_REC'))
    load('WineRecognitionData.mat');
elseif(strcmp(type, 'ABALONE'))
    load('AbaloneData.mat');
elseif(strcmp(type, 'ECOLI'))
    load('EcoliData.mat');
elseif(strcmp(type, 'GLASS'))
    load('GlassData.mat');
elseif(strcmp(type, 'POKER2'))
    load('PokerData2.mat');
elseif(strcmp(type, 'ZOO'))
    load('ZooData.mat');
elseif(strcmp(type, 'THYROID'))
    load('ThyroidData.mat');
elseif(strcmp(type, 'WINE'))
    load('WineData.mat');
    
end

classData = classData;
labels =y;

disp('COMPLETE DATA LOADING');


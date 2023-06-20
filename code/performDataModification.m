function modifiedClassData =performDataModification(classData, type)
global VERBOSE;
if(VERBOSE); disp('performDataModificationNClass'); end

if(~strcmp(type, 'NONE'))  
    modifiedClassData =[];
    type =upper(type);
    
    if(strcmp(type, 'BIAS'))
        %Append the bias to the begining of the feature vectors
        
        npts =size(classData,1);
        modifiedClassData =[ones(npts, 1) classData];
        
    elseif(strcmp(type, 'Z_SCORE'))
        
        modifiedClassData = zscore(classData);
        
        disp('');
        
    elseif(strcmp(type, 'D_NORM'))
        
        dMin = abs(min(classData));
        data = classData -repmat(dMin,size(classData,1),1)+eps;
        
        dMax = abs(max(data));
        data = data ./repmat(dMax,size(classData,1),1);
        
        modifiedClassData =data;
        
    end%End
else
    modifiedClassData =classData;
end

%Find and replace all NaNs
ind = isnan(modifiedClassData);
modifiedClassData(ind) =0;

disp('COMPLETE DATA MODIFICATIONS');

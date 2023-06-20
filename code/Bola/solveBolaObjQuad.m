function WOptimal = solveBolaObjQuad(w0, classDataSets, labels)

global VERBOSE;
global PLOT;
global MAX_ITR;

classLabels = unique(labels);
K =length(unique(labels));            %Get the K number of classes
classVarianceMatrices =cell(1,K);   %Hold the class variance matrices

for kk=1:K
    ind = find(labels == classLabels(kk));    
    R = classDataSets(ind,:)' *classDataSets(ind,:);
    classVarianceMatrices{kk} =R;
end
   
%newCost =bolaObj(w0, classVarianceMatrices); %Intial evaluation

% newCost =bolaObjNew(w0, classDataSets, labels); %Intial evaluation
newCost =bolaQuadObjLegendreNew(w0, classDataSets, labels); %Intial evaluation
f(1) = newCost;

if(VERBOSE)
    disp(['ITER: ' num2str(0)]);
    disp(['    Current cost = ' num2str(newCost)]);
end

[bolaVal,~, ~] =checkBolaConditions(w0, classVarianceMatrices);
bolaNorm(1) =bolaVal;

wOld    =w0;
nn      =1;

costBreakPoint(1)       =0;
bolaBreakPoint(1)       =0;

convgThresh     = 1e-4;
costFlag        =false;
criticalPointFlag =false;

while(1)
    if(VERBOSE)
        disp(['ITER: ' num2str(nn)])
    end
    
    if(isequal(nn,MAX_ITR))
        break;
    end
            
    %Perform polar decomposition
    for kk=1:K
        w =wOld(:,kk);
        R =classVarianceMatrices{kk};
        
        AX(:,kk) = R*w;
    end    
    [U,~,V]=svd(AX);
    U=U(:,1:K);
    
    WO =U*V'; %Find new weight vectors
    
    %Update with new orthonormal basis
    for kk=1:K
        wNew(:,kk) =WO(:,kk);
    end
         
    %Evaluation the objective function
    prevCost = newCost;
%     newCost=bolaQuadObjLegendreNew(wNew, classVarianceMatrices);
    newCost =bolaQuadObjLegendreNew(wNew, classDataSets, labels);
    f(nn+1) =newCost;    
        
     if(VERBOSE) disp(['    Current cost = ' num2str(newCost)]); end
    
    costDiff = abs(newCost - prevCost);
    if(costDiff < convgThresh)
        if(VERBOSE); disp('Minimum change in cost reached --TRUE'); end
        
        costBreakPoint(nn+1) =1;
        costFlag =true;
    else
        costBreakPoint(nn+1) =0;
    end
        
    [bolaVal,bola_S_symmetric_Flag, bolaGlobalCondition] =checkBolaConditions(wNew, classVarianceMatrices);
    bolaNorm(nn+1) =bolaVal;

    if(bolaVal < convgThresh)
         if(VERBOSE); disp('Bola A(X) = XS -- TRUE'); end
         bolaBreakPoint(nn+1) =1;
         criticalPointFlag = true;
    else
        bolaBreakPoint(nn+1) =0;
    end    

    if(bola_S_symmetric_Flag)
         if(VERBOSE); disp('Bola S=A(X)X is symmetric -- TRUE'); end
    end
    
    if(bolaGlobalCondition)
         if(VERBOSE); disp('Bola Global Condition Satisfied (SYMM & PSD) -- TRUE'); end
    end
               
    %Update
    wOld =wNew;
    
    nn=nn+1;

    %If all flags have been triggered break.
    if(costFlag)
        disp('Conditions Reached')
        break;
    end
       
    if(VERBOSE); disp(' '); end
end%End while(1)
if(VERBOSE); disp(' '); end

if(PLOT)
    figure;hold on;
    plot(1:length(f), f, 'r-');
    title('Objective Cost');
    ind = find(costBreakPoint >0);
    plot(ind, f(ind), 'b*');    
    hold off;
    
    figure;hold on;
    plot(bolaNorm, 'g-');
    title('');
    ind = find(bolaBreakPoint>0);
    plot(ind, bolaNorm(ind),'b*');
    hold off;

end

WOptimal =wNew;


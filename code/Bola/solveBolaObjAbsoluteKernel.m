function AOptimal =solveBolaObjAbsoluteKernel(A0, kGramMatrix, labels)

global VERBOSE;
global PLOT;
global MAX_ITR;

K =length(unique(labels)); 
isigma =1;

gHalfInv = computegHalfInv(kGramMatrix);

A0 = A0'*gHalfInv;
A0 = A0';
IC = A0'*kGramMatrix*A0;
norm(eye(size(IC,1)) - IC);
epsilon =1e-4;

Zkik = computeZvalues(A0, kGramMatrix, labels, isigma,'mabs', epsilon);

% if(VERBOSE)
%     classLabels = unique(labels);
%     for k=1:K
%         ind = find(labels == classLabels(k));  %Get items in a class
%         
%         disp(['Class ' num2str(k) ' sum ' num2str(sum(Zkik(ind)))]);
%     end
% end



Yker = computeYker(Zkik, kGramMatrix, labels, isigma);
[Uker, Sker, Vker, Lambda] = computeA(Yker, gHalfInv, labels);

norm(eye(size(kGramMatrix,1)) - (gHalfInv * kGramMatrix * gHalfInv));

newCost = bolaAbsoluteObjKernelLegendreNew(A0, kGramMatrix, Zkik, labels, isigma, epsilon);

f(1) = newCost;
if(VERBOSE)
    disp(['Iter ' num2str(0) ':']);
    disp(['    Current cost = ' num2str(newCost)]);
end

%Check the constraint 
% Val = -Yker +(Lambda*A0*G)  %According to EQ(30) this value should be zero

A = (Uker*Vker'*gHalfInv)';  %EQ(34)
A'*kGramMatrix*A;
% A = (Uker*Vker')';  %EQ(34)
AOld    =A;
nn      =1;


costBreakPoint(1) =0;

convgThresh     = 1e-6;
costFlag        =false;

while(1)
    if(VERBOSE)
        disp(['ITER: ' num2str(nn)])
    end
    
    if(isequal(nn,MAX_ITR))
        break;
    end
              
    %Perform polar decomposition   
    Zkik = computeZvalues(AOld, kGramMatrix, labels, isigma,'mabs', epsilon);
    Yker = computeYker(Zkik, kGramMatrix, labels, isigma);
    [Uker, Sker, Vker, Lambda] = computeA(Yker, gHalfInv, labels);  

    %Evaluation the objective function
    prevCost = newCost;
    newCost = bolaAbsoluteObjKernelLegendreNew(AOld, kGramMatrix, Zkik, labels, isigma, epsilon);
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
    
    
                  
    %Update
    ANew = (Uker*Vker'*gHalfInv)';  %EQ(34)
    ANew'*kGramMatrix*ANew; %Print

    
    if(norm(ANew - AOld) < convgThresh)
        disp('Norm Condition Reached');
%         break;
    end
    AOld =ANew;
    
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
       
end

AOptimal =ANew;



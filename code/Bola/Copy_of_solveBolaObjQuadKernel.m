function AOptimal =solveBolaObjQuadKernel(A0, kGramMatrix, classDataSets, labels)

global VERBOSE;
global PLOT;
global MAX_ITR;

K =length(unique(labels)); 
isigma =1;

[Yker, gHalfInv, Z] = computePolarMatrix(A0, kGramMatrix, classDataSets, labels, isigma);

%Expanding using its singular value decomposition
[Uker, Sker, Vker] = svd(Yker*gHalfInv);

%Grab the first K colums from Sker
Sker = Sker(:,1:K);

%EQ(33)
Lambda =Uker*Sker*Uker';  %IMPORTANT THIS needs to be checked with EQ(33)

Vker = Vker(:,1:K); 
A = Uker*Vker'*gHalfInv;  %EQ(34)

if(VERBOSE)
    %Check initialization matrix of orthogonality
    disp(['Norm of A is ' num2str(norm(A))]);
    if( (norm(A)-1.0) < 10e-6)
        disp('System is orthogonal');
    end
end

%Check the constraint 
Ik = A*kGramMatrix*A';  %EQ(28)
%Val = -Yker +(Lambda*A*G)  %According to EQ(30) this value should be zero


newCost =bolaQuadObjKernelLegendreNew(A, kGramMatrix, classDataSets, Z, Lambda, labels, isigma);

f(1) = newCost;

if(VERBOSE)
    disp(['Iter ' num2str(0) ':']);
    disp(['    Current cost = ' num2str(newCost)]);
end

AOld    =A';
nn      =1;

costBreakPoint(1) =0;

convgThresh     = 1e-4;
costFlag        =false;

while(1)
    if(VERBOSE)
        disp(['ITER: ' num2str(nn)])
    end
    
    if(isequal(nn,MAX_ITR))
        break;
    end
            
    %Perform polar decomposition
    [Yker, gHalfInv, Z] = computePolarMatrix(AOld, kGramMatrix, classDataSets, labels, isigma);

    [Uker, ~, Vker] = svd(Yker*gHalfInv);

    %Update with new orthonormal basis
    Vker = Vker(:,1:K);
    ANew = Uker*Vker'*gHalfInv;  %EQ(34)
        
    %Evaluation the objective function
    prevCost = newCost;
    newCost =bolaQuadObjKernelLegendreNew(ANew, kGramMatrix, classDataSets, Z, Lambda, labels, isigma);

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
    AOld =ANew';
    
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



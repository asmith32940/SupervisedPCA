function wOptimal =solveSplineObj(classDataSets)

global VERBOSE;

if(VERBOSE); disp('solveSplineObj'); end

K =numel(classDataSets);            %Get the K number of classes

classVarianceMatrices =cell(1,K);   %Hold the class variance matrices

for kk=1:K
    R =cov(classDataSets{kk}',1);

    classVarianceMatrices{kk} =R;
end

NK = numel(classVarianceMatrices);
for kk=1:NK
    K =classVarianceMatrices{kk};
    K=-(K);
    
    TK =K;
    for ii=1:NK
        if(~isequal(ii,kk))
            TK = TK+classVarianceMatrices{ii};
        end
    end
    
    [w, ~]=eigs(TK,1,'lm');
    wOptimal(:,kk) =w;
end


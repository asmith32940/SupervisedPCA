function  [EW, classVarianceMatrices]=bolaAbsObjLegendreNew(W, classDataSets, labels, epsilon)

classLabels = unique(labels);
K =length(unique(labels));          %Get the K number of classes
classVarianceMatrices =cell(1,K);

%Compute the overall data mean
EW=0;
for kk=1:K
   
    ind = find(labels == classLabels(kk));  %Get items in a class
    w = W(:,kk);                            %Get class weight vector  
    Ck = numel(ind);
    
    T = 0;
    for ii=1:numel(ind)
        x = classDataSets(ind(ii),:)';
        T =T+ (w'*x)/sqrt((w'*x)^2+epsilon);
    end
    
    ZkikCons = zeros(Ck,1);
    LW =0;
    for ik=1:numel(ind)

        xik = classDataSets(ind(ik),:)';
        
        %EQ(26)
        zkik = ((w'*xik)/sqrt((w'*xik)^2+epsilon)) -((1/Ck)*T);
        
        %EQ(7),EQ(25)  Absolute value of the inner product
        LW = LW + (-zkik*(w'*xik) - (epsilon*sqrt(1-zkik^2)));
        if(~isreal(LW))
            disp('Imaginary');
        end    
        ZkikCons(ik) = zkik;
    end
    
    classVarianceMatrices{kk} = classDataSets(ind,:)' *classDataSets(ind,:);
   
%     disp(['Constraint:  '  num2str(sum(ZkikCons))]);
    
    EW = EW+LW;
end

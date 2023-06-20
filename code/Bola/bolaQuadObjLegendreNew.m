function EW=bolaQuadObjLegendreNew(W, classDataSets, labels)

classLabels = unique(labels);
K =length(unique(labels));          %Get the K number of classes

%Compute the overall data mean
EW=0;
for kk=1:K
   
    ind = find(labels == classLabels(kk));  %Get items in a class
    w = W(:,kk);                            %Get class weight vector
    
    Ck = numel(ind);
    SumXik = sum(classDataSets(ind,:))';
    
    ZkikCons = zeros(Ck,1);
    LW =0;
    for ik=1:Ck
        xik = classDataSets(ind(ik),:)';

        %EQ(24)
        zkik = w'* (xik-((1/Ck)*SumXik));
        
        %Quadratic
        LW = LW + (-zkik*(w'*xik) + (.5*zkik^2));  %EQ(4),EQ(23)
        
        if(~isreal(LW))
            disp('Imaginary')
        end
            
        ZkikCons(ik) = zkik;
    end
    
%     disp(['Constraint:  '  num2str(sum(ZkikCons))]);
    
    EW = EW+LW;
end




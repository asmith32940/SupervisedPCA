function Zkik = computeZvalues(A, classDataSets, labels, isigma,obj, epsilon)


classLabels = unique(labels);
K =length(unique(labels));          %Get the K number of classes
N = size(classDataSets,1);

Zkik = zeros(1,N);

if(strcmp(obj,'quad')==1)
 
    for k=1:K
        ind = find(labels == classLabels(k));  %Get items in a class
        Ck = numel(ind);
        CkMean =0;
        
        for ik=1:Ck
            xik =classDataSets(ind(ik),:)';
            
            xikSum =0;
            for j=1:N
                xj = classDataSets(j,:)';
                xikSum = xikSum+ (A(j,k)*gaussianKernel(xj,xik, isigma));
            end
            Zkik(ind(ik)) = xikSum;
        end
        
        CkMean = mean(Zkik(ind(:)));
        Zkik(ind(:)) = Zkik(ind(:)) - CkMean;
        
        disp('');
        
        
    end
elseif(strcmp(obj, 'mabs')==1)
    
    for k=1:K
        ind = find(labels == classLabels(k));  %Get items in a class
        Ck = numel(ind);
        CkMean =0;
        
        for ik=1:Ck
            xik =classDataSets(ind(ik),:)';
            
            xikSum =0;
            for j=1:N
                xj = classDataSets(j,:)';
                xikSum = xikSum+ (A(j,k)*gaussianKernel(xj,xik, isigma));
            end
            
            Zkik(ind(ik)) = xikSum./sqrt(xikSum*xikSum+epsilon);
                        
        end
        
        CkMean = mean(Zkik(ind(:)));
        Zkik(ind(:)) = Zkik(ind(:)) - CkMean;
        
        disp('');
        
        
    end
end

function Yker = computeYker(Zkik, classDataSets, labels, isigma)

classLabels = unique(labels);
K =length(unique(labels));          %Get the K number of classes
N = size(classDataSets,1);
Yker = zeros(K,N);

for k=1:K
    ind = find(labels == classLabels(k));  %Get items in a class
    Ck = numel(ind);
    
    for j=1:N
        xj = classDataSets(j,:)';
        
        xikSum =0;
        for ik=1:Ck
            xik =classDataSets(ind(ik),:)';         
            xikSum = xikSum + (Zkik(ind(ik))*gaussianKernel(xj,xik,isigma));
        end
        Yker(k,j) = xikSum;
    end
    
    disp('');
end


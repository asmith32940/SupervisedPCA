function Lker=bolaAbsoluteObjKernelLegendreNew(A, GramMatrix, Z, labels, isigma, epsilon)

classLabels = unique(labels);
K =length(unique(labels));          %Get the K number of classes
N = size(GramMatrix,1);

Lker =0;

for k=1:K
    ind = find(labels == classLabels(k));  %Get items in a class
    Ck = numel(ind);
    
    classSum =0;
    for ik=1:Ck
        xik =GramMatrix(ind(ik),:)';
        dataSum =0;
        for j=1:N
            xj = GramMatrix(j,:)';
            
            %EQ(29) & EQ(11) Extension to the kernel
            %The constraints can be expressed using a Lagrange parameter
            %matrix to obtain the following Lagrangian
%             dataSum = dataSum + ( Z(ind(ik))* A(j,k)*gaussianKernel(xj,xik, isigma);
            dataSum = dataSum + ( Z(ind(ik))* A(j,k)*gaussianKernel(xj,xik, isigma) - epsilon*(1-Z(ind(ik))*Z(ind(ik)) ) );
            
        end%End for j=1:N
        
        classSum = classSum + dataSum;
    end%End for ik=1:Ck
    
    Lker = Lker + classSum;
end%End for k=1:K

Lker = -Lker;

disp('');
% Lker = -Lker + trace(Lambda * (A'*GramMatrix*A - Ik));


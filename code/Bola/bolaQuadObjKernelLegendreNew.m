function Lker=bolaQuadObjKernelLegendreNew(A, GramMatrix, classDataSets, Z, Lambda, labels, isigma)

classLabels = unique(labels);
K =length(unique(labels));          %Get the K number of classes
N = size(classDataSets,1);

Lker =0;

for k=1:K
    ind = find(labels == classLabels(k));  %Get items in a class
    Ck = numel(ind);
    
    for ik=1:Ck
        xik =classDataSets(ind(ik),:)';
        for j=1:N
            xj = classDataSets(j,:)';
            
            %EQ(29) & EQ(11) Extension to the kernel
            %The constraints can be expressed using a Lagrange parameter
            %matrix to obtain the following Lagrangian
            Lker = Lker + ( Z(k,ind(ik))* ...
                            A(k,j)*gaussianKernel(xj,xik, isigma)+...
                            trace(Lambda*(A*GramMatrix*A'-ik)))+...
                            ((.5)*Z(k,ind(ik))^2 );
            
        end%End for j=1:N
    end%End for ik=1:Ck
    
end%End for k=1:K

Lker = -Lker;



function [Yker, gHalfInv, Z] = computePolarMatrix(A0, kGramMatrix, classDataSets, labels, isigma)

classLabels = unique(labels);
K =length(unique(labels));            %Get the K number of classes
N = size(classDataSets,1);

Yker = zeros(K,N);
Z = zeros(K,N);

for k=1:K
    ak = A0(:,k);
    ind = find(labels == classLabels(k));  %Get items in a class
    Ck = numel(ind);
                   
        for j=1:N
            
            Yker(k,j) =0;
            xj =classDataSets(j,:)';
 
            for ik=1:Ck
                xik =classDataSets(ind(ik),:)';
                
                %Calculate zkik
                zkik =0;
                for jj=1:N
                    akj = ak(jj);                    
                    xjj =classDataSets(jj,:)';
                    
                    termTwo =0;  %Second term summation of EQ(35)
                    for ikk=1:Ck
                        xkikik = classDataSets(ind(ikk),:)';
                        termTwo = termTwo +(gaussianKernel(xjj, xkikik, isigma));
                    end
                    
                    %EQ(35)
                    zkik = zkik + (akj * (gaussianKernel(xjj, xik, isigma)-((1/Ck)*termTwo)));
                    
                end%End for jj=1:N
                Z(k,ind(ik)) = zkik;
                
                %EQ(31)  Yker of size KxN
                Yker(k,j) = Yker(k,j) + (zkik*gaussianKernel(xj, xik, isigma));
            end%End for ik=1:Ck
        end%End for j=1:N
end%for k=1:K

%Take the square root of the matrix
gHalf =sqrtm(kGramMatrix);

if ~isreal(gHalf)
    disp('Not a good matrix (~imaginary)');
    gHalf =real(gHalf);
end

if ~is_positive_semi_matrix(gHalf)
    disp('Not a good matrix (~psd)');
end

gHalf =gHalf ./norm(gHalf,'fro');
gHalfInv =pinv(gHalf);


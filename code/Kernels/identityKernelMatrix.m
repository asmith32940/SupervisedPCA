function ikm = identityKernelMatrix(points1, points2)

n1 =size(points1,2);
n2 =size(points2,2);

ikm = zeros(n2, n1);
for ii =1:n1
   
    x1= points1(:,ii);
    for jj=1:n2
        x2= points2(:,jj);
        ikm(jj, ii) = dot(x1,x2);
    
    end
    disp('');

end

disp('');

function mkm= mlpKernelMatrix(points1, points2, t, s)

n1 =size(points1,2);
n2 =size(points2,2);

mkm = zeros(n2, n1);
for ii =1:n1
    
    x1= points1(:,ii);
    for jj=1:n2
        x2= points2(:,jj);
        mkm(jj, ii) = mlpKernel(x1, x2, t, s);
        
    end
    disp('');
    
end

disp('');
end

function val = mlpKernel(x1, x2, t, s)

    val =tanh(s*dot(x1,x2) + t^2);

end
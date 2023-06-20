function pkm = polynomialKernelMatrix(points1, points2, t, d)

n1 =size(points1,1);
n2 =size(points2,1);

pkm = zeros(n2, n1);
for ii =1:n1
   
    x1= points1(ii,:);
    for jj=1:n2
        x2= points2(jj,:);
        pkm(jj, ii) = polynomialKernel(x1, x2, t, d);
    
    end
    disp('');

end

disp('');
end

function val = polynomialKernel(x1, x2, t, d)

    val =(dot(x1,x2) +t)^d;

end
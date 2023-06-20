function pkm = tangentKernelMatrix(points1, points2, t, d)

n1 =size(points1,2);
n2 =size(points2,2);

pkm = zeros(n2, n1);
for ii =1:n1
   
    x1= points1(:,ii);
    for jj=1:n2
        x2= points2(:,jj);
        pkm(jj, ii) = tanKernel(x1, x2, t, d);
    
    end
    disp('');

end

disp('');
end

function val = tanKernel(x1, x2, a, c)

    val =tanh((a*dot(x1,x2)) +c);

end
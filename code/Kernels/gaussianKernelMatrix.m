function gkm = gaussianKernelMatrix(classData, isigma)

n =size(classData,1);
gkm = zeros(n, n);

for ii =1:n    
    x1= classData(ii,:);
    
    for jj=1:n
        x2= classData(jj,:);
%         gkm(ii, jj) = exp(-(dot(x1,x1))/(2*isigma^2)) *  ...
%             exp((dot(x1,x2))/(isigma^2)) * ...
%             exp((-dot(x2,x2))/(2*isigma^2));
        
        gkm(ii,jj) = exp(-0.5*sqrt(sum((x1-x2).^2))/(isigma^2));
        
        disp('');
    end   
end

disp('')
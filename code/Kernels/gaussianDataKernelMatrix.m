function gkm = gaussianDataKernelMatrix(classData, projectionData, isigma)

N =size(classData,1);
NP =size(projectionData,1);
gkm = zeros(NP, NP);

for ii =1:N    
    x1= classData(ii,:);
    
    for jj=1:NP
        x2= projectionData(jj,:);
%         gkm(ii, jj) = exp(-(dot(x1,x1))/(2*isigma^2)) *  ...
%             exp((dot(x1,x2))/(isigma^2)) * ...
%             exp((-dot(x2,x2))/(2*isigma^2));
        
        gkm(ii,jj) = exp(-0.5*sqrt(sum((x1-x2).^2))/(isigma^2));
    end   
end

disp('')
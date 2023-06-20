function km=projectToKernelSpace(classData, kernelParams)

type = kernelParams.type;
if(strcmp(type, 'rbf'))    
    %RBF Kernel
    km = gaussianKernelMatrix(classData, kernelParams.sigma);
    
%     if is_positive_semi_matrix(km)
%         disp('Good matrix');
%     else
%         disp('');
%     end
    
    
    [V,D] = eig(km);
    
    ev=diag(D);
    thresh = 1.0e-4;
    ind = find(ev< thresh);
    
    if (numel(ind) > 0)
        ev(ind) = thresh;
        Dp = diag(ev);
        
        km = V*Dp*V';
    end
    
    km = real((km+km')/2);
    
    eig(km);
%     disp('');

    
elseif(strcmp(type, 'poly'))
    km = polynomialKernelMatrix(classData, classData, 0, 2);
%     if is_positive_semi_matrix(km)
%         disp('Good matrix');
%     else
%         disp('');
%     end
    
    
    [V,D] = eig(km);
    
    ev=diag(D);
    thresh = 1.0e-4;
    ind = find(ev< thresh);
    
    if (numel(ind) > 0)
        ev(ind) = thresh;
        Dp = diag(ev);
        
        km = V*Dp*V';
    end
    
    km = real((km+km')/2);
    
    eig(km);
%     disp('');

    

end

    

function km=projectDataToKernelSpace(classData, projectionData, kernelParams)

type = kernelParams.type;
if(strcmp(type, 'rbf'))    
    %RBF Kernel
    km = gaussianDataKernelMatrix(classData, projectionData, kernelParams.sigma);      

elseif(strcmp(type, 'poly'))
    km = polynomialKernelMatrix(projectionData, classData, 0, 2);

    disp('')
end

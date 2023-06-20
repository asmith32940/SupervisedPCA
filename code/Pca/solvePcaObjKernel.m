function wOptimal =solvePcaObjKernel(kClassDataSets, labels)

K =length(unique(labels));            %Get the K number of classes

% May not be needed, if so this must be moved to right after the data has
% been projected to the kernel space.
if ~is_positive_semi_matrix(kClassDataSets)
%     disp('Gram matrix is: (~psd)');
    
     
    [V,D] = eigs(kClassDataSets);
    evals=diag(D);

    eigAdj = 1-evals(end);  %Get the smallest eigenvalue.
    D = D + (eye(size(D)).*(eigAdj + eps));
    
    kClassDataSets =V*D*V';
end

[coeff,~,~,~] =pca(kClassDataSets, 'NumComponents',K);

wOptimal =coeff;



function I = axisAngleProjTrain(alphas, classDataSets,labels)


K =length(unique(labels));          %Get the K number of classes
N = size(classDataSets,1);

% acos(dot()/dot(norm(), norm()))
catSpaceAxis = eye(K);
dang = zeros(N,K);
% tdang = zeros(N,K);

for k=1:K
    
    daxis = catSpaceAxis(:,k);
%     catAxis  = abs(alphas' * alphas(:,k));
    
    for ii=1:N
        xi = abs(classDataSets(ii,:));
%         xi = abs(xi/norm(xi));
        
        dang(ii,k) = acos(dot(daxis, xi)/dot(norm(daxis), norm(xi)))/pi*180;
%         tdang(ii,k) = acos(dot(catAxis, xi)/dot(norm(catAxis), norm(xi)))/pi*180;
    end
    
end

[M, I]= min(dang,[],2);

disp('');
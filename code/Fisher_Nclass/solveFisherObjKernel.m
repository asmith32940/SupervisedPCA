function WOptimal=solveFisherObjKernel(kClassDataSets, labels)
global VERBOSE;

classLabels = unique(labels);
K =length(unique(labels));            %Get the K number of classes

D = size(kClassDataSets,2);


%Get the class mean for all data sets
M =zeros(D,K);
for kk=1:K
    ind = find(labels == classLabels(kk));    

    M(:,kk) =mean(kClassDataSets(ind,:),1)';
end
GM =mean(kClassDataSets,2);   %Global mean

SW =zeros(D,D);
for kk=1:K
    ind = find(labels == classLabels(kk));        
    S =cov(kClassDataSets(ind,:),1);
        
%     if(VERBOSE)
%         disp(['Rank S(' num2str(kk) '): ' num2str(rank(S))]);
%     end
    
    %The generalization of the within-class covariance matrix to the case
    %of K classes
    SW = SW+S;
end

% if(VERBOSE)
%     disp(['Rank SW: ' num2str(rank(SW))]);
% end

%Generate the generalization of the total covariance matrix
SB =zeros(D,D);
for kk=1:K
    ind = find(labels == classLabels(kk));        
    N =numel(ind);

    T =N*((M(:,kk)-GM)*(M(:,kk)-GM)');
    
    SB =SB+T;
end

% if(VERBOSE)
%     disp(['Rank SB: ' num2str(rank(SB))]);
% end


%Calculate the weight matrix.
%SB has rank at most equal to (K-1) and so there are at most (K-1)
%nonzero eigenvalues.
[W, ~] = eigs(pinv(SW)*SB, K-1, 'LM');

WOptimal =W;




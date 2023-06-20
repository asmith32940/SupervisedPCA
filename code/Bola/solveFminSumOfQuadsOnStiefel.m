function optBasis = solveFminSumOfQuadsOnStiefel(initBasis, classData, labels)

maxIters = 100;

convgThresh = 1e-4;

% Form xx' quadratic sums for each the terms in the heterogenous mix of
% quadratics.
K = length(unique(labels));
classLabels = unique(labels);
quadTerms = cell(1,K);
for kk=1:K

    %ADD
    ind = find(labels == classLabels(kk));
    data = classData(ind,:);
    quadTerms{kk} = (data'*data);
end

prevCost = 99999999999;
currCost = computeCost(initBasis, quadTerms, K);
f(1) = currCost;
%     disp(['Iter ' num2str(0) ':']);
%     disp(['    Current cost = ' num2str(currCost)]);

iter = 0;
% Run the Bola iteration.
optBasis = initBasis;
costDiff = abs(currCost - prevCost);
while((costDiff > convgThresh) && (iter < maxIters))
    iter = iter + 1;
    f(iter+1) = currCost;
    % Form A(X), should be numDims x numClasses
    AX = [];
    for kk = 1 : K
        AX = [AX quadTerms{kk}*optBasis(:,kk)];
    end
    % Get the polar decomp of AX via SVD (economy SVD).
    [U, ~, V] = svd(AX,0);
    
    %S = V*S*V'; % No need to actully compute S from Bola paper.
    
    % Should be numDims x numClasses
    optBasis = U*V';
    prevCost = currCost;
    currCost = computeCost(optBasis, quadTerms, K);
    costDiff = abs(currCost - prevCost);
    disp(['Iter ' num2str(iter) ':']);
    disp(['    Current cost = ' num2str(currCost) ', Cost diff = ' num2str(costDiff)]);
end

function cost = computeCost(basis,quadTerms, numClasses)

cost = 0;
for kk = 1 : numClasses
    cost = cost - basis(:,kk)'*quadTerms{kk}*basis(:,kk);
end


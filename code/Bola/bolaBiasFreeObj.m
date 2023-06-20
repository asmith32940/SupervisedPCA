function [fval, classVarianceMatrices]=bolaBiasFreeObj(W, classDataSets, labels)

classLabels = unique(labels);
K =length(unique(labels));            %Get the K number of classes
classVarianceMatrices =cell(1,K);

fval =0;
%
for kk=1:K
    ind = find(labels == classLabels(kk));
    wk =W(:,kk);
    nk = numel(ind);
    
    data = classDataSets(ind,:);
    
    %Get class projected mean
    xo=0;
    for jj=1:nk
        xi = data(jj,:);
        xo=xo+(dot(wk,xi)/nk);
    end
    
    %Calculate total class 
    fcval = 0;
    for jj=1:nk
        xi = data(jj,:);
        fcval = fcval +((dot(wk,xi)-xo))^2;
    end

    classVarianceMatrices{kk} = classDataSets(ind,:)' *classDataSets(ind,:);

    fval = fval+fcval;
end
fval =-(0.5)*fval;

disp('');



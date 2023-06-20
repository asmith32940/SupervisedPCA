function dataSet =genRandomSubSets(labels)

dataSet =[];

%Training 60%
%Validation 20%
%Test 20%
trainPercent = .60;
testPercent = .20;

K = unique(labels);

for kk=1:length(K)
    ind = find(labels == K(kk));
    pind = ind(randperm(length(ind))); %Randomly shuffle elements in a vector
    
    nElements = numel(pind);
    
    nTrainSamps = round(nElements*trainPercent);    
    dataSet(kk).train = pind(1:nTrainSamps);
    pind(1:nTrainSamps) = [];
    
    nTestSamps = round(nElements*testPercent);
    dataSet(kk).test = pind(1:nTestSamps);
    pind(1:nTestSamps) = [];
    
    dataSet(kk).valid = pind;
   
end


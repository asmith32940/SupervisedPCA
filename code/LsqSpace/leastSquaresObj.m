function WOptimal =leastSquaresObj(classDataSets, labels)

global VERBOSE;

if(VERBOSE); disp('leastSquaresObj'); end

classLabels = unique(labels);
K =length(unique(labels));            %Get the K number of classes

TN =size(classDataSets,1);

D =size(classDataSets,2);  %Get the data dimensions 
X =zeros(TN, D);
tgt =zeros(TN,K);

idx =1;
for kk=1:K
    ind = find(labels == classLabels(kk));    
    N =length(ind);
    
    tgtVector =zeros(1,K);
    tgtVector(kk) =1;
    
    X(ind,:)    =classDataSets(ind,:);
    tgt(ind,:)  =repmat(tgtVector,N,1);
end

%Calculate the weight vector
lambda = 1;

PT=(X'*X);
I = eye(size(PT));
I(size(PT,2),size(PT,2)) = 0;      %Regularizes the lambda

%This is the pseudo inverse solution for w
WOptimal = inv(PT + (lambda * I))*X'*tgt;


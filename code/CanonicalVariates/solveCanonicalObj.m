function WOptimal =solveCanonicalObj(classDataSets, labels)

global VERBOSE;

if(VERBOSE); disp('solveCanonicalObj'); end

classLabels = unique(labels);
K =length(unique(labels));      %Get the K number of classes

D =size(classDataSets,2);       %Get the data dimensions 
N =size(classDataSets,1);       %Total number os samples

M  =zeros(D,1);
B = zeros(D,D);
for kk =1:K
    ind = find(labels == classLabels(kk));
    data =classDataSets(ind,:);
    Nk =length(ind);
    
    
    %Get the class mean for all data sets
    mk =mean(data,1)';    
    B = B + Nk*cov(data);
    
    %Generalization of the between-class covariance matrix
    M = M+(Nk*mk);
end

%Normalize the between-class mean
M =(1/N)*M;

for kk = 1:K
    ind = find(labels == classLabels(kk));
    data =classDataSets(ind,:);
    Nk =length(ind);
    
    mk =mean(data,1)';
    A(:,kk) =sqrt(Nk)*(mk-M);
end
Q       =orth(classDataSets');  % basis to represent the solution
AA      =Q'*A*A'*Q;             % modified between scatter
B       =Q'*B*Q;                % modified within scatter
cB      =chol(B); 
invcB   =inv(cB);

% [Wt evals]=eigs(invcB'*AA*invcB,K);

options.tol = 1e-50;
% [Wt, ~, ~]  =svds(invcB'*AA*invcB,K, 'L', options);
[Wt, ~, ~]  =svd(full(invcB'*AA*invcB));
WOptimal    =Q*(Wt(:,1:K));
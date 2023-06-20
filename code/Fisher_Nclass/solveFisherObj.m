function WOptimal=solveFisherObj(classDataSets, labels)

global VERBOSE;

if(VERBOSE); disp('solveFisherObj'); end

classLabels = unique(labels);
K =length(unique(labels));          %Get the K number of classes
D =size(classDataSets,2);           %Get the data dimensions 
N =size(classDataSets,1);           %Total number os samples

SW =zeros(D,D);
m  =zeros(D,1);

for kk=1:K
    ind = find(labels == classLabels(kk));
    data =classDataSets(ind,:)';
    Nk =length(ind);
    
    %Get the class mean for all data sets
    mk =mean(data,2);
    
    %Class covariance matrix
    mData = data - repmat(mk,1,Nk);
    Sk = mData*mData';
    
    if(VERBOSE)
        disp(['Rank S(' num2str(kk) '): ' num2str(rank(Sk))]);
    end
    
    %The generalization of the within-class covariance matrix to the case
    %of K classes
    SW = SW+Sk;
    
    
    %Generalization of the between-class covariance matrix
    m = m+(Nk*mk);
end

if(VERBOSE)
    disp(['Rank SW: ' num2str(rank(SW))]);
end

%Normalize the between-class mean
m =(1/N)*m;

SB =zeros(D,D);
for kk=1:K
    ind = find(labels == classLabels(kk));
    data =classDataSets(ind,:);
    Nk =length(ind);
    
    %Get the class mean for all data sets
    mk =mean(data,1)';
    SB =SB +(Nk*((mk-m)*(mk-m)'));
end

if(VERBOSE)
    disp(['Rank SB: ' num2str(rank(SB))]);
end


%Calculate the weight matrix.
%SB has rank at most equal to (K-1) and so there are at most (K-1)
%nonzero eigenvalues.
WOptimal =[];
[WOptimal, ~] = eigs(pinv(SW)*SB, K-1, 'LM');





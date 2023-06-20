function [rapSymmFlag, rapPosSemiFlag, rapVal] = ...
                            checkRapcsakCondition(W, classVarianceMatrices)
                        
global NORM_MIN;
global VERBOSE;

%Create the KD X KD block diagonal matrix, where K=number of classes and
%D=dimensions of the data
K =numel(classVarianceMatrices);
D =size(classVarianceMatrices{1},1);  %Check the size of a variance 
                                      %matrix to get the data dimension
A = eye(K*D);

stidx=1;
edidx=D;
for kk=1:K
     A(stidx:edidx,stidx:edidx) = classVarianceMatrices{kk};
     stidx =edidx+1;
     edidx =(kk+1)*D;     
end

%RAPCSAK CONDITIONS
X0 =W(:);

n=numel(X0);
S1 =zeros(n,n);
S1 =A*(X0*X0');

rstidx=1;
redidx=D;
SA=zeros(size(S1));

%Construct Rapcsak symmetric S1 matrix
for kk=1:K
    cstidx=1;
    cedidx=D;
    for jj=1:K
        ST=S1(rstidx:redidx,cstidx:cedidx);
        TST=trace(ST)*eye(D);
        SA(rstidx:redidx,cstidx:cedidx)=TST;
        
        cstidx =cedidx+1;
        cedidx =(jj+1)*D;
    end
    rstidx =redidx+1;
    redidx =(kk+1)*D;
end
S1=(SA+SA')./2;

rapSymmFlag =0;
rapPosSemiFlag =0;
if(is_symmetric_matrix(S1))
    rapSymmFlag =1;
    
    rapVal =norm((A*X0)-(S1*X0),'fro');
%     if( (rapVal < NORM_MIN) && (VERBOSE) )
%         disp('A*X0 = S1*X0 -- TRUE');
%     end
    
    C =(A-S1);
    if(is_positive_semi_matrix(C))
        rapPosSemiFlag =1;
    end
    
end


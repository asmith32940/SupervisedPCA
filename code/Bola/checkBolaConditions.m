function [bolaVal,bolaFlag, bolaGlobalCondition] =checkBolaConditions(W, classVarianceMatrices)

K =numel(classVarianceMatrices);

%BOLA CONDITIONS
for kk=1:K
    w =W(:,kk);
    R =classVarianceMatrices{kk};
    
    AX(:,kk) =R*w;
    X(:,kk)  =w;
end

%Construct Bola symmetric S matrix
S  =X'*AX;  %kxk symmetric matrix

%Following the frobenius norm to zero display how the optimization is
%converging
bolaVal =norm(AX-(X*S),'fro');  %Bolla Eq(2.2)

bolaFlag =is_symmetric_matrix(S);

%To characterize the critical points of the functional, let us denote by
%A(X) =[Ax, ....,Ax] and X=[x, ....x]  the nxk matrices containing the
%enumerated vectors as their columns..  I will be shown that for an optimal
%orthonormal system [x,...x]

%Lemma 3.1 If X is a global maximum of the functional then the
%corresponding matrix S=XA(X) is positive semidefinite
bolaGlobalCondition = false;
if( (bolaFlag) && (is_positive_semi_matrix(S)) )
    bolaGlobalCondition = true;
end

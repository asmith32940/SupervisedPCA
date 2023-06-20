function [C1, C2, f] =solveSplineObj_linear_2Class(Xi, Xj)

Ni =size(Xi,2);
Nj =size(Xj,2);

K1 =zeros(size(Xi,1),size(Xi,1));
for ii=1:Ni
    ix =Xi(:,ii);
    K1 = K1 + (ix * ix');
end
K1 = K1 .* (1/Ni);

K2 =zeros(size(Xj,1),size(Xj,1));
for jj=1:Nj
    ij =Xj(:,jj);
    K2 = K2 + (ij * ij');
end
K2 = K2 .* (1/Nj);


[C1, ~]=eigs((K2-K1),1,'sm');
[C2, ~]=eigs((K1-K2),1,'lm');

f = -(C1'*K1*C1) + (C2'*K1*C2) -(C2'*K2*C2) +(C1'*K2*C1)

disp('Complete')
function [C1, C2, C3, f] =solveSplineObj_linear_3Class(Xi, Xj, Xk)

Ni =size(Xi,2);
Nj =size(Xj,2);
Nk =size(Xk,2);

Kj =zeros(size(Xi,1),size(Xi,1));
for ii=1:Ni
    ix =Xi(:,ii);
    Kj = Kj + (ix * ix');
end
Kj = Kj .* (1/Ni);

Kk =zeros(size(Xj,1),size(Xj,1));
for jj=1:Nj
    ij =Xj(:,jj);
    Kk = Kk + (ij * ij');
end
Kk = Kk .* (1/Nj);

Kl =zeros(size(Xk,1),size(Xk,1));
for kk=1:Nk
    ik =Xk(:,kk);
    Kl = Kl + (ik * ik');
end
Kl = Kl .* (1/Nk);

[C1, ~]=eigs((-Kj+Kk+Kl),1,'lm');
[C2, ~]=eigs((Kj-Kk+Kl),1,'lm');
[C3, ~]=eigs((Kj+Kk-Kl),1,'lm');

f = -(C1'*Kj*C1) + (C2'*Kj*C2) + (C3'*Kj*C3) ...
    -(C2'*Kk*C2) + (C1'*Kk*C1) + (C3'*Kk*C3) ...
    -(C3'*Kl*C3) + (C1'*Kl*C1) + (C2'*Kl*C2);

disp('Complete')

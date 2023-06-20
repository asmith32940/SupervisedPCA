function [C1, C2, C3, f] =solveSplineObj_kernel_NClass(Xi, Xj, Xk)

Ni =size(Xi,2);
Nj =size(Xj,2);
Nk =size(Xk,2);

type ='poly';
XiPrj = projectToKernelSpace(Xi,[Xi, Xj, Xk], type);
Kj =zeros(size(XiPrj,1),size(XiPrj,1));
for ii=1:Ni
    ix =XiPrj(:,ii);
    Kj = Kj + (ix * ix');
end
Kj = Kj .* (1/Ni);

XjPrj = projectToKernelSpace(Xj,[Xi, Xj, Xk], type);
Kk =zeros(size(XjPrj,1),size(XjPrj,1));
for jj=1:Nj
    ij =XjPrj(:,jj);
    Kk = Kk + (ij * ij');
end
Kk = Kk .* (1/Nj);

XkPrj = projectToKernelSpace(Xk,[Xi, Xj, Xk], type);
Kl =zeros(size(XkPrj,1),size(XkPrj,1));
for kk=1:Nk
    ik =XkPrj(:,kk);
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

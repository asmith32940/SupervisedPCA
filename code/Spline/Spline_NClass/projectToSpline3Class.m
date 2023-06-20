function [Zi, Zj, Zk] =projectToSpline3Class(P, Xi, Xj, Xk)

dim =size(P,2);

%Project Xi
Ni = size(Xi,2);
Zi = zeros(dim,Ni);
% for ii=1:Ni
%     b = Xi(:,ii);
%     Zi(:,ii) =(P'*b);
% end

for ii=1:Ni
    ix = Xi(:,ii);
    Zi(1,ii) =dot(ix,P(:,1));
    Zi(2,ii) =dot(ix,P(:,2));
end




%Project Xj
Nj = size(Xj,2);
Zj = zeros(dim,Nj);
% for jj=1:Nj
%     b = Xj(:,jj);
%     Zj(:,jj) =(P'*b);
% end

for jj=1:Nj
    ij = Xj(:,jj);
    Zj(1,jj) =dot(ij,P(:,1));
    Zj(2,jj) =dot(ij,P(:,2));
end


%Project Xk
Nk = size(Xk,2);
Zk = zeros(dim,Nk);
% for jj=1:Nj
%     b = Xj(:,jj);
%     Zj(:,jj) =(P'*b);
% end

for kk=1:Nk
    ik = Xk(:,kk);
    Zk(1,kk) =dot(ik,P(:,1));
    Zk(2,kk) =dot(ik,P(:,2));
end

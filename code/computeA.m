function [Uker, Sker, Vker, Lambda, gHalfInv] =computeA(Yker, gHalfInv, labels)

K =length(unique(labels)); 

% %Expanding using its singular value decomposition
[Uker, Sker, Vker] = svd(Yker*gHalfInv);

Sker = Sker(:,1:K);
Vker = Vker(:,1:K); 
Lambda =Uker*Sker*Uker';  %IMPORTANT THIS needs to be checked with EQ(33)



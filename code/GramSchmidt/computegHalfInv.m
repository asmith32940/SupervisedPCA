function gHalfInv = computegHalfInv(kGramMatrix)


[V,D] = eig(kGramMatrix);

ev = 1./sqrt(diag(D));

gHalfInv = V*diag(ev)*V';


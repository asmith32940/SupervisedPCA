function result = is_positive_semi_matrix(A)


[~,D] = eigs(A);

D=diag(D);

%Is positive-semi definite
% if D(length(D)) >= 0
%     result = 1;
% else
%     result = 0;
% end


%Is positive definite
if D(length(D)) > 0
    result = 1;
else
    result = 0;
end

function result = is_positive_matrix(A)


[~,p] = chol(A);

if p == 0
    result = 1;
else
    result = 0;
end

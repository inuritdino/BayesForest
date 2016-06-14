function B = mat_vec_subtraction(A,v)

% Subtracts from each row of the matrix A the vector v.
% If A is (n x m)-matrix, then v needs to be m-vector.

s = size(A);
B = zeros(s(1),s(2));
for i = 1:s(2)
    B(:,i) = A(:,i)-v(i);
end
function a = is_structured(obj)
% A tree is structured if every branch except one has a parent.

n = obj.number_of_branches;

if nnz(obj.parent(1:n)) >= n - 1
    a = true;
else
    a = false;
end
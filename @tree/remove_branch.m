function obj = remove_branch(obj, k)
% Remove the last <k> branches of the tree. The coresponding rows of the
% member variables are set to empty and the number of objects is decreased
% by <k>.

n = obj.number_of_branches;

% If one tries to remove more than the number of branches, remove all.
if k > n
    k = n;
end

obj.radius(n-k+1:n) = 0;
obj.length(n-k+1:n) = 0;
obj.start_point(n-k+1:n,:) = 0;
obj.end_point(n-k+1:n,:) = 0;
obj.axis(n-k+1:n,:) = 0;
obj.is_trunk(n-k+1:n) = 0;

for i = 0:k-1
    obj.children{n-i} = [];
    
    if obj.parent(n-i) ~= 0
        
        if obj.extension(obj.parent(n-i)) >= n-k+1
            obj.extension(obj.parent(n-i)) = 0;
            
        end
        
        obj.children{obj.parent(n-i)} = setdiff(obj.children{obj.parent(n-i)},n-i);
    end   
end

obj.parent(n-k+1:n) = 0;
obj.extension(n-k+1:n) = 0;

obj.number_of_branches = n - k;
function obj = append_branch(obj, new_branch)
% Append a new tree to the tree. The new branch has to be class BRANCH. 
%
% See also BRANCH.

    obj.number_of_branches = obj.number_of_branches + 1;
    n = obj.number_of_branches;
    
    obj.radius(n) = new_branch.radius;
    obj.length(n) = new_branch.length;
    obj.axis(n,:) = new_branch.axis;
    obj.start_point(n,:) = new_branch.start_point;
    obj.is_trunk(n) = new_branch.is_trunk;
    
    p = new_branch.start_point(:);
    q = p + new_branch.axis(:)*new_branch.length;
    obj.end_point(n,:) = q;
end
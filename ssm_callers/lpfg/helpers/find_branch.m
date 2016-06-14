function branch_ind = find_branch(branches,order,parent_ind,cyl)
% Find branches
% USE:
% 	branch_ind = find_branch(branches,order,parent_ind,cyl)
% Find the branches of the given order in the struct <branches>. The struct
% must have <.order> element in it. For details, look into gen_tree()
% function. If <parent_ind> is specified, the function gives the indices of
% the branches having the parent branches with index = <parent_ind>.
% <.parent_ind> must exist in the struct <branches>. If <cyl> is specified,
% the function gives the branch indices of the branches having the
% specified cyl indices.
%
% Found branches correspond to ALL specified criteria simultaneously 
% (maximum of them is 3: order, parent_ind and cyl). To ignore some 
% criteria put them empty, i.e. [].
%
% All branch indices are within <branches> structure.
%
% SEE ALSO: intersect

if(~isstruct(branches))
    error('Error: not a struct!');
end
if(nargin < 4)
    cyl = [];
    if(nargin < 3)
        parent_ind = [];
        if(nargin < 2)
            order = [];
        end
    end
end
if(isempty(cyl) && isempty(parent_ind) && isempty(order))
    disp('Error: searching criteria are not set.')
    return;
end


branch_ind = [];
x1 = [];% searching results on order
x2 = [];% searching results on parent
x3 = [];% searching results on cyl
for ii=1:length(branches)
    if((~isempty(order)) && (branches(ii).order(1) == order))
        x1 = cat(2,x1,ii);
    end
    if((~isempty(parent_ind)) && (~isempty(branches(ii).parent_num)) &&...
            (branches(ii).parent_num == parent_ind))
        x2 = cat(2,x2,ii);
    end
    if((~isempty(cyl)) && (any(cyl == branches(ii).cyl_ind)))
        x3 = cat(2,x3,ii);
    end
end
if(~isempty(order))
    branch_ind = x1;
end
if(~isempty(parent_ind))
    if(~isempty(branch_ind))
        branch_ind = intersect(branch_ind,x2);
    else
        branch_ind = x2;
    end
end
if(~isempty(cyl))
    if(~isempty(branch_ind))
        branch_ind = intersect(branch_ind,x3);
    else
        branch_ind = x3;
    end
end

end

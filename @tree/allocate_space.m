function obj = allocate_space(obj, min_space)
% Function that allocates more space for the member variables of a tree
% object. Zeros are added to the end of all member variables.

n = max(100, min_space);

obj.radius = [obj.radius; zeros(n,1)];
obj.length = [obj.length; zeros(n,1)];
obj.start_point = [obj.start_point; zeros(n,3)];
obj.end_point = [obj.end_point; zeros(n,3)];
obj.axis = [obj.axis; zeros(n,3)];
obj.is_trunk = [obj.is_trunk; sparse(n,1)];
obj.parent = [obj.parent; sparse(n,1)];
obj.extension = [obj.extension; sparse(n,1)];
obj.children = [obj.children; cell(n,1)];
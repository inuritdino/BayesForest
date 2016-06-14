function obj = squeeze(obj)
% Function deletes the empty rows of the member variables to free memory.

n = obj.number_of_branches;

obj.radius = obj.radius(1:n);
obj.length = obj.length(1:n);
obj.start_point = obj.start_point(1:n,:);
obj.end_point = obj.end_point(1:n,:);
obj.axis = obj.axis(1:n,:);
obj.is_trunk = obj.is_trunk(1:n);
obj.parent = obj.parent(1:n);
obj.extension = obj.extension(1:n);
obj.children = obj.children(1:n);
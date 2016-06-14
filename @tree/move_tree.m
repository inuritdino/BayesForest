function obj = move_tree(obj,pos)
% Move the tree-object to the new position POS=(X,Y,Z), namely, the first
% cylinder of the tree (seed) will be moved to the specified position.

N = obj.number_of_branches;
vec = pos - obj.start_point(1,:);

obj.start_point(1:N,:) = obj.start_point(1:N,:) + repmat(vec,N,1);
obj.end_point(1:N,:) = obj.start_point(1:N,:) + bsxfun(@times,obj.length(1:N),obj.axis(1:N,:));

end
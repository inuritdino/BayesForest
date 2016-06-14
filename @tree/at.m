function br = at(obj,i)

if length(i) > 1
   error('Only one index allowed.');
end

br = branch(obj.radius(i),obj.length(i),obj.start_point(i,:),obj.axis(i,:));
br.is_trunk = obj.is_trunk(i);
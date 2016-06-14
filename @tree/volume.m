function v = volume(obj,varargin)
% Compute the volume of the tree. The volume is defined as the sum of the
% volumes of the cylinders that form the tree. The function accepts
% additional parameters.
% 'Trunk'       Only compute volume of the trunk.
% 'Branches'    Volume of the branches.

n = obj.number_of_branches;
I = 1:n;

i = 1;
while i <= size(varargin,2)

    if ischar(varargin{i})

        switch lower(varargin{i})
            case 'trunk'
                I = logical(obj.is_trunk(1:n));
            case 'branches'
                I = not(logical(obj.is_trunk(1:n)));
            otherwise
                disp(['Ignoring unknown parameter: ' varargin{i}]);
        end
    end
    i = i + 1;
end

if nnz(I) == 0
    v = 0;
else
    v = pi*sum(obj.radius(I).^2.*obj.length(I));
end
function obj = set_parent(obj, whose, to_what, is_child)
% Function to set the parent of a branch. The function also updates the
% child and extension info automatically. 
% Inputs:
%   <whose>     N x 1 double    Indexes of the branches whose parent
%                               information should be changed. 
%   <to_what>   1 x 1 double    Index of the new parent branch.
%   <is_child>  N x 1 boolean   Truth value indicating whether the branch
%                               is a child of the <to_what> branch or an
%                               extension.

% If the <is_child> parameter is not given assume that all <whose> are
% children of the <to_what>.
if nargin == 2
    is_child = true(size(whose));
end

if max(size(to_what)) == 1
    to_what = to_what(ones(max(size(whose))));
end    

% Get possible old parents of the branches <whose>.
old_parents = obj.parent(whose);

% Set the new parent value.
obj.parent(whose) = to_what;

for i = 1:length(whose)
    
    % If the new parent value is not zero, add the branch <whose> as a new
    % child  or an extension of the <to_what>.
    if to_what(i)
        
        % If <whose> is a child then append new child.
        if is_child(i)

            obj.children{to_what(i)} = [obj.children{to_what(i)}; whose(i)];
        
        % Otherwise check if <to_what> already had an extension and replace
        % it. Set the parent of the old extension to zero to ensure that no
        % branch has two extensions.
        else
            if obj.extension(to_what(i))
                obj.parent(obj.extension(to_what(i))) = 0;
            end

            % Set new extension value.
            obj.extension(to_what(i)) =  whose(i);
        end
    end
    
    % If <whose> had an old parent. Remove <whose> as a child or an
    % extension of the old parent.
    if old_parents(i) ~= 0
        
        % If <whose> was an extension, set extension of the old parent to
        % zero.
        if obj.extension(old_parents(i)) == whose(i)
            obj.extension(old_parents(i)) = 0;
        % If <whose> was a child, remove it but keep the rest.
        else
            obj.children{old_parents(i)} = setdiff(obj.children{old_parents(i)},whose(i));
        end
        
    end    
end
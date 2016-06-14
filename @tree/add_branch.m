function obj = add_branch(obj,varargin)
% Add a branch or branches to a tree. If there is only one
% additional argument it is assumed to be a class branch
% object or a cell of such objects. Otherwise there should
% be 5 parameters in a specific order: radius, length,
% start_point, axis and is_trunk.

    switch nargin
        
        case 1
            return;
            
        case 2
            new_branch = varargin{1};
            
            if iscell(new_branch)
                for i = 1:max(size(new_branch))

                    obj = obj.append_branch(new_branch{i});
                end
            else
                obj = obj.append_branch(new_branch);
            end
            
        otherwise
            r = varargin{1};
            l = varargin{2};
            p = varargin{3};
            a = varargin{4};
            t = varargin{5};
            
            if size(p,2) ~= 3
                p = p';
            end
            
            if size(a,2) ~= 3
                a = a';
            end
            
            a_n = sqrt(sum(a.^2,2));
            
            a = a./a_n(:,ones(1,3));
            
            nn = max(size(r));
            no = obj.number_of_branches;
            
            if nn+no > size(obj.radius,1)
                obj = obj.allocate_space(nn);
            end
            
            obj.radius(no+1:no+nn) = r(:);
            obj.length(no+1:no+nn) = l(:);
            obj.start_point(no+1:no+nn,:) = p;
            obj.end_point(no+1:no+nn,:) = p + bsxfun(@times,a,l(:));
            obj.axis(no+1:no+nn,:) = a;
            obj.is_trunk(no+1:no+nn) = t(:);
            
            obj.number_of_branches = obj.number_of_branches + nn;
    end
end
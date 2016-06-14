classdef tree
% Class that describes a single tree with multiple branches. A tree can
% have a type and numerous branches which are presented by cylinders. The
% class stores the radii, lengths, start points and axis directions of each
% cylinder.
%
% See also 
% TREE.TREE                     TREE.ADD_BRANCH,        TREE.AT, 
% TREE.BRANCH_DISTRIBUTION,     TREE.DRAW,              TREE.DRAW_BRANCHES,
% TREE.FILL_GAPS,               TREE.HEIGHT,            TREE.HIGH_POINT, 
% TREE.IS_STRUCTURED,           TREE.LOW_POINT,         TREE.ORGANIZE,
% TREE.REMOVE_BRANCH,           TREE.ROOT,              TREE.SET_PARENT, 
% TREE.SET_TYPE,                TREE.SQUEEZE,           TREE.VOLUME
  
    properties
        type = '';
        number_of_branches = 0;

        radius;
        length;
        start_point;
        end_point;
        axis;
        is_trunk;
        
        parent;
        extension;
        children;
        
        trunk_direction = [0 0 1];
    end
    
    methods
        function obj = tree(varargin)
        % Constructor for class TREE. Optional arguments can be a string,
        % which is used for the type of the tree and a three-element vector
        % which is used for the direction of the trunk.
        %
        % Examples:
        % p = tree();
        % p = tree('Pine');
        % p = tree('Spruce',[0 0.1 0.9]);
            
            n = 2000;
            
            obj.radius = zeros(n,1);
            obj.length = zeros(n,1);
            obj.start_point = zeros(n,3);
            obj.end_point = zeros(n,3);
            obj.axis = zeros(n,3);
            obj.is_trunk =sparse(n,1);
            
            obj.parent = sparse(n,1);
            obj.extension = sparse(n,1);
            obj.children = cell(n,1);
            
            for i = 1:nargin

                if ischar(varargin{i})
                    obj.type = varargin{i};
                elseif isnumeric(varargin{i}) && length(varargin{i}) == 3
                    obj.trunk_direction = varargin{i};
                end
            end
        end
        
        function [h,i] = low_point(obj)
            % Compute the lowest point of the tree. Only the start and end
            % points of the branches are used in the computations.
            % h = r.low_point()      returns the height of the lowest point
            %                        in the trunk direction of the tree r.
            % [h,i] = r.low_point()  returns also the index of the lowest
            %                        branch.
            
            n = obj.number_of_branches;
            
            h = [obj.start_point(1:n,:); obj.end_point(1:n,:)]*obj.trunk_direction(:);
            [h,i] = min(h);
            
            if i > n
                i = i - n;
            end
        end
        
        function [h,i] = high_point(obj)
            % Compute the highest point of the tree. Only the start and end
            % points of the branches are used in the computations.
            % h = r.low_point()      returns the height of the highest point
            %                        in the trunk direction of the tree r.
            % [h,i] = r.low_point()  returns also the index of the highest
            %                        branch.
            
            n = obj.number_of_branches;
            
            h = [obj.start_point(1:n,:); obj.end_point(1:n,:)]*obj.trunk_direction(:);
            [h,i] = max(h);
            
            if i > n
                i = i - n;
            end
        end
        
        function I = root(obj)
            % Returns the index of the root branch(es) of the tree. The
            % branches that do not have a set parent value are classified
            % as root. If there are multiple roots, the indexes are sorted
            % according to the height of the corresponding branches. The
            % first element of the return value always corresponds with the
            % lowest root branch.
            
            n = obj.number_of_branches;
            I = find(obj.parent(1:n) == 0);
            
            if length(I) > 1
                h = [obj.start_point(I,:); obj.end_point(I,:)]*obj.trunk_direction(:);
                h = min(h(1:n),h(n+1:end));
                [~,I] = sort(h);
            end
        end
    end
    
    methods(Access=protected)
        
        obj = allocate_space(obj, min_space)
    
        obj = append_branch(obj, new_branch)
%        draw_branches(obj,I, varargin)
    end
end
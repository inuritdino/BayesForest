function varargout = sorted_trunk(obj, varargin)
% Return tree stem cylinder properties sorted according to height. By
% default the sorting order is ascending. There can be one to four output
% arguments:
%
% SP = SOMETREE.SORTED_TRUNK(...)
% [SP,R] = SOMETREE.SORTED_TRUNK(...)
% [SP,R,H] = SOMETREE.SORTED_TRUNK(...)
% [SP,R,H,AX] = SOMETREE.SORTED_TRUNK(...)
%
% where SP is the start point, R the radius, H the lenght, AX the axis
% direction. The sorting order can be defined as either 'ascend' or
% 'descend' as follows.
%
% [...] = SOMETREE.SORTED_TRUNK('descend')
%

n = nargout;

% Select cylinders that are in the trunk.
I = find(obj.is_trunk(1:obj.number_of_branches));

if not(isempty(I)) && n > 0
    
    if nargin > 1 && ischar(varargin{1})
        order = varargin{1};
    else
        order = 'ascend';
    end
    
    % Starting points of the cylinders.
    sp = obj.start_point(I,:);
    [~, J] = sort(sp(:,3),1,order);    
    varargout{1} = sp(J,:);
    
    if n > 1
        % Radius of the selected cylinders.
        r = obj.radius(I);
        varargout{2} = r(J);
    end
    
    if n > 2
        % Radius of the selected cylinders.
        h = obj.length(I);
        varargout{3} = h(J);
    end
    
    if n > 3
        % Radius of the selected cylinders.
        ax = obj.axis(I,:);
        varargout{4} = ax(J,:);
    end
    
end
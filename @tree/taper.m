function varargout = taper(obj,varargin)
% Taper of the tree stem. Function either draws the taper as a function of
% height or returns the corresponding values. Plotting is done in the
% current set of axis if no output arguments are defined. All input
% arguments are passed to the plot command. If one output argument is 
% defined the diameter and corresponding height values are returned as 
% columns of a matrix. In the case of two output arguments the two vectors
% are returned separately. Examples:
%
%   SOMETREE.TAPER()
%   SOMETREE.TAPER('-','Color',[0.4 0.1 0.9])
%   T = SOMETREE.TAPER()
%   [H, D] = SOMETREE.TAPER()
%
% SEE ALSO TREE.LEAN, TREE.SWEEP

% Get sorted trunk position and radius values.
[sp,r] = obj.sorted_trunk;

if size(varargin,2) > 0 && isnumeric(varargin{1})
    zimax = find(sp(:,3)-sp(1,3) > varargin{1},1,'first');
    if ~isempty(zimax)
        sp = sp(1:zimax-1,:);
        r = r(1:zimax-1);
    end
    if size(varargin,2) > 1
        varargin = varargin(2:end);
    end
end

if not(isempty(sp))
    
    % Compute cylinder diameters.
    d = 2*r;

    % Move the height values to start from zero.
    h = sp(:,3) - sp(1,3);

    % Select output according to number of output arguments.
    switch nargout

        % If one output, return matrix.
        case 1
            varargout{1} = [h, d];

        % Return height and diameter values separately.
        case 2
            varargout{1} = h;
            varargout{2} = d;
            
        % By default plot taper as a function of the height.
        otherwise
            if size(varargin,1) > 0
                plot(h,d,varargin{:});
            else
                plot(h,d);
            end

    end
else
    % If there are no cylinders in the tree classified as part of the stem,
    % raise an error.
    error('No stem in tree. Unable to form taper.')
end
function varargout = sweep(obj,varargin)
% Analyze the sweep of the stem of a tree. Sweep is denined as the maximum
% pitch per unit, which by default is 1.0 (m). There are three different
% ways to call the function depending on the number of output arguments:
%
%   SOMETREE.SWEEP(...)             (1)
%   SW = SOMETREE.SWEEP(...)        (2)
%   [SW, PH] = SOMETREE.SWEEP(...)  (3)
%
% With no output arguments (1), the pitch is plotted as a function of 
% height. The sweep location is highlighted and the value is shown as the 
% title of the figure. In forms (2) and (3) SW is the sweep and PH is a
% N-by-2-matrix where the first column contains the pitch values and the
% second column the corresponding height values.
%
% The function recognizes the following input arguments as pairs:
% 'Unit', <UNIT>        Set the unit against which the pitch is computed.
%                       By default <UNIT> is 1.0.
%
% 'Intervals', <INTS>   Each unit is divided into intervals whose number is
%                       defined by <INTS>. This means that pitch is
%                       computed at frequency <UNIT>/<INTS>. The number of
%                       intervals must be a positive even interger. By
%                       default <INTS> is 10.
%
% 'MaxH', <MAXH>        Maximum height to which the pitch analysis can
%                       advance. By default the whole stem is analyzed.
%
% 'Plot', ...           This input argument can be used to pass arguments
%                       to the plot command when no output arguments are
%                       defined. All the input arguments after 'Plot' are
%                       passed to the plot command, so naturally all other
%                       settings should be defined before this one.
%
% SEE ALSO TREE.TAPER, TREE.LEAN

% Get sorted trunk position values.
sp = obj.sorted_trunk('descend');

if not(isempty(sp))

    % Default values for settings.
    intervals = 10;
    unit = 1;
    limith = max(sp(:,3));
    params{1} = '-';

    % Process additional arguments.
    i = 1;
    while i <= size(varargin,2)

        if ischar(varargin{i})

            switch lower(varargin{i})
                case 'intervals'
                    intervals = varargin{i+1};
                    % Ensure even number.
                    intervals = intervals + mod(intervals,2);
                    i = i + 1;
                case 'unit'
                    unit = varargin{i+1};
                    i = i + 1;
                case 'maxh'
                    limith = varargin{i+1};
                    i = i + 1;
                case 'plot'
                    params = varargin(i+1:end);
                    break;
                otherwise
                    disp(['Ignoring unknown parameter: ' varargin{i}]);
            end
        end
        i = i + 1;
    end
    
    % Half of the number of intervals.
    half = intervals/2;

    % Height values for the stem cylinders.
    h = sp(:,3);
    
    % Extreme values for height.
    minh = min(h);
    maxh = max(h);
    
    % Height positions at which the pitch is computed at.
    H = minh:(unit/intervals):min(maxh,limith);
    N = length(H);
    
    % Estimate stem location at height values <H>.
    X = obj.properties_at(H);
    
    H = H - minh;
    
    % Pitch values.
    ph = zeros(N-intervals,1);
    
    % Start at half a unit above bottom and advance to half a unit below
    % height limit.
    for i = half+1:N-half
        
        % Direction vector for the line from the cylinder half a unit below
        % to the cylinder half a unit above the current point.
        dir = X(i+half,:)-X(i-half,:);
        dir = dir/norm(dir);
        
        % Compute distance to line which equals the pitch.
        ph(i-half) = distances_to_line(X(i,:),dir,X(i-half,:));
    end
    
    % Find maximum which is the sweep.
    [sw, j] = max(ph);
    % Normalize by unit.
    sw = sw/unit;
    
    % Height values corresponding to the computed pitch values.
    x = H(half+1:N-half);
    
    % Select output according to number of output arguments.
    switch nargout

        % If one output, return single value: the sweep.
        case 1
            varargout{1} = sw;

        % Return sweep and pitch values.
        case 2
            varargout{1} = sw;
            varargout{2} = [ph, x'];
            
        % By default plot sweep as a function of the height.
        otherwise
            set_hold_off = ~ishold;
            
            % Plot pitch.
            plot(x, ph,params{:});
            hold on;
            % Plot sweep.
            plot(x(j),ph(j),'ro');
            title(['Sweep is ' num2str(ph(j)) ' at the height ' num2str(x(j))]);
            
            if set_hold_off
                hold off;
            end
    end
    
else
    % If there are no cylinders in the tree classified as part of the stem,
    % raise an error.
    error('No stem in tree. Unable to analyse sweep.')
end
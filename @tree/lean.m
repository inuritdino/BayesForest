function varargout = lean(obj,varargin)
% Analyze the lean of a the stem of a tree. The function computes the
% center points of each cylinder. These points are then transformed into
% a coordinate system where the z-axis is defined by the direction of the
% first cylinder, and projected on a plane perpendicular to that direction.
% Polar coodinates of the center points are recorded. The following three
% forms are available depending on number of output arguments:
%
%   SOMETREE.LEAN(...)
%   L = SOMETREE.LEAN(...)
%   [R, A] = SOMETREE.LEAN(...)
%
% Where the Nx2 matrix L has two columns R and A, and these vectors contain
% the radii and angles of the projected center points, respectively. With
% no output arguments a visualization of the lean property is plotted.
%
% The function recognizes the following input arguments as pairs:
% 
% 'North', <NORTH>      Direction of North. The projected coordinates are
%                       projected so that the y-axis is aligned with the
%                       North-direction. <NORTH> should be a 3x1 or 1x3
%                       vector. The default value is [0 1 0].
%
% 'Circles', <CIRC>     This setting only affects the plot created if no
%                       output arguments are defined. By default four
%                       evenly distributed circles are drawn no highlight
%                       the scale of the lean. The number and radii of
%                       these circles can be controlled with this setting.
%                       <CIRC> can be either a vector or radii values or
%                       the key string 'off' in which case no circles are
%                       drawn.
%
% SEE ALSO TREE.TAPER, TREE.SWEEP

% Get sorted trunk position and radius values.
[sp,r,h,ax] = obj.sorted_trunk('descend');

if not(isempty(sp))
    
    % Process additional arguments.
    i = 1;
    
    % Default settings.
    north = [0 1 0];
    draw_circles = true;
    circles = [];
    interval = 1;

    while i <= size(varargin,2)

        if ischar(varargin{i})

            switch lower(varargin{i})
                case 'north'
                    north = varargin{i+1};
                    north = (north(:)/norm(north))';
                    i = i + 1;
                case 'circles'
                    if ischar(varargin{i+1}) && strcmp(varargin{i+1},'off')
                        draw_circles = false;
                    else
                        circles = varargin{i+1};
                        if isempty(circles)
                            draw_circles = false;
                        end
                    end
                    i = i + 1;
                case 'maxheight'
                    zimax = find(sp(:,3)-sp(end,3) < varargin{i+1},1,'first');
                    if ~isempty(zimax)
                        h = h(zimax:end);
                        sp = sp(zimax:end,:);
                        ax = ax(zimax:end,:);
                        r = r(zimax:end);
                    end
                    i = i + 1;
                case 'interval'
                    interval = floor(varargin{i+1});
                    i = i + 1;
                otherwise
                    disp(['Ignoring unknown parameter: ' varargin{i}]);
            end
        end
        i = i + 1;
    end
    
    % Compute cylinder center points.
    cp = sp + bsxfun(@times,ax,h)/2;
    % Move lowest to Origin.
    cp = bsxfun(@minus,cp,cp(end,:));
    % Axis direction of the bottom cylinder.
    a = ax(end,:);
    
    % Generate a coordinate system for the lowest cylinder.
    if a(1) == 0 && a(2) == 0
        M = eye(3);
    else
        a1 = [a(2) -a(1) 0];
        a1 = a1/norm(a1);
        a2 = cross(a,a1);
        a2 = a2/norm(a2);

        M = [a1(:) a2(:) a(:)];
    end
    
    % Transform center points to the new coordinate system.
    cp = cp*M;
    
    % Transform north vector into the new system.
    north = north*M;
    north = north(1:2);
    north = north/norm(north);
    east = [north(2) -north(1)];
    
    % Rotate points so that north is up.
    cp2 = cp(:,1:2)*(-[east(:), north(:)]');
    
    % Compute distance to Origin.
    d = sqrt(cp2(:,1).^2 + cp2(:,2).^2);
    
    % Compute angle of data points around the Origin.
    alfa = asin(cp2(:,2)./d);

    % Correct values for second and third quaters.
    q2 = (cp2(:,1) < 0 & cp2(:,2) >= 0);
    q3 = (cp2(:,1) < 0 & cp2(:,2) <  0);
    alfa(q2) =  pi - alfa(q2);
    alfa(q3) = -pi - alfa(q3);

    % Select output according to number of output arguments.
    switch nargout

        % If one output, return matrix.
        case 1
            varargout{1} = [d, alfa];

        % Return height and diameter values separately.
        case 2
            varargout{1} = d;
            varargout{2} = alfa;
            
        % By default plot taper as a function of the height.
        otherwise
            
            % Variables for drawing circles.
            N = 31;
            ang = linspace(0,2*pi,N);
            cosang = cos(ang);
            sinang = sin(ang);
            
            % Total length of the tree. Used to select color.
            lenh = abs(cp(end,3)-cp(1,3));
            
            set_hold_off = ~ishold;
            
            % Plot segmented line of center points.
            plot(cp2(1:interval:end,1)',cp2(1:interval:end,2)','k-');
            hold on;
            
            % Draw circles at each center point presenting the width of the
            % stem at that point.
            for i = 1:interval:size(cp,1)
                
                % Percent of total height.
                hval = cp(i,3)/lenh;
                % Select color.
                color = ones(1,3)*(0.2+hval*0.6);
                
                % Circle circumference.
                xp = r(i)*cosang;
                yp = r(i)*sinang;
                
                % Plot circle at center point.
                plot(cp2(i,1) + xp,cp2(i,2) + yp,'-','color',color);%,'linewidth',2);
            end
            
            % Draw scale highlights.
            if draw_circles
                
                % Default circle radii are determined from maximum
                % distance.
                if isempty(circles)
                    n = 2;
                    
                    maxd = max(d+r);
                    
                    circles = linspace(0,maxd,n+1);
                    circles = circles(2:end);
                end
                
                for i = 1:length(circles);
                    xp = circles(i)*cosang;
                    yp = circles(i)*sinang;
                    plot(xp,yp,'-','color',[1 0 0]);
                    text(cos(pi/4)*circles(i)+0.05*circles(end),-cos(pi/4)*circles(i),sprintf('%.2f',circles(i)),...
                          'VerticalAlignment','top','color',[1 0 0]);
                end
            end
            
            if set_hold_off
                hold off;
            end
            
            axis equal;
    end
else
    % If there are no cylinders in the tree classified as part of the stem,
    % raise an error.
    error('No stem in tree. Unable to form taper.')
end
function draw(obj, varargin)

    param = cell(7,1);
    n = 0;
    trunk = false;
    both = true;
    min_radius = 0;

    i = 1;
    while i <= size(varargin,2)

        if ischar(varargin{i})

            switch lower(varargin{i})

                case 'trunk'
                    trunk = true;
                    both = false;
                case 'branches'
                    both = false;
                case {'minradius', 'radius'}
                    min_radius = varargin{i+1};
                case {'points','curve','color'}
                    n = n + 1;
                    param{n} = varargin{i};
                    n = n + 1;
                    param{n} = varargin{i+1};
                    i = i + 1;
                case 'vector'
                    n = n + 1;
                    param{n} = varargin{i};
            end
        elseif isnumeric(varargin{i}) && ishandle(varargin{i})
            %figure(h);
            figure(varargin{i});
        end
        i = i + 1;
    end

    param = param(1:n);
    
    if min_radius > 0 
        I = obj.radius(1:obj.number_of_branches) >= min_radius;
    else
        I = true(obj.number_of_branches,1);
    end

    if both
        J = true(obj.number_of_branches,1);
    elseif trunk
        J = obj.is_trunk(1:obj.number_of_branches);
    else
        J = not(obj.is_trunk(1:obj.number_of_branches));
    end

    obj.draw_branches(I & J, param{:});
    axis equal;
end
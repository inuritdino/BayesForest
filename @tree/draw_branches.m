function h = draw_branches(obj,I, varargin)
% Draw a branch (cylinder). The additional arguments in
% varargin are passed to the cylinder-command which computes
% points for a unit cylinder. The first argument is assumed
% to be the number of points to be used in the circles of
% the cylinder. The second is a generator curve which is
% used for the profile of the cylinder. 
% 
% See also CYLINDER.

cdata = [];
color = [];
curve = 1;
points = 20;
vector = false;

% Process additional arguments.
i = 1;
while i <= size(varargin,2)

    if ischar(varargin{i})

        switch lower(varargin{i})
            case 'points'
                points = varargin{i+1};
                i = i + 1;
            case 'curve'
                curve = varargin{i+1};
                i = i + 1;
            case 'color'
                cdata = varargin{i+1};
                i = i + 1;
            case 'vector'
                vector = true;
                if i < size(varargin,2) && isnumeric(varargin{i+1})
                    color = varargin{i+1};
                    i = i + 1;
                end
            otherwise
                disp(['Ignoring unknown parameter: ' varargin{i}]);
        end
    end
    i = i + 1;
end

[x y z] = cylinder(curve,points);

s = size(x);

if isscalar(cdata)
    map = colormap;
    if isempty(color)
        color = map(cdata);
    end
    cdata = repmat(cdata,s(1),s(2));
elseif min(size(cdata)) == 1 && max(size(cdata)) == 3
    if isempty(color)
        color = cdata(:);
    end
    temp(:,:,1) = repmat(cdata(1),s(1),s(2));
    temp(:,:,2) = repmat(cdata(2),s(1),s(2));
    temp(:,:,3) = repmat(cdata(3),s(1),s(2));
    cdata = temp;
end

if isempty(color) || isscalar(color)
    color = [0 0 1];
end

n = nnz(I);

x = repmat(x(:)',n,1);
y = repmat(y(:)',n,1);
z = repmat(z(:)',n,1);

x = bsxfun(@times,x,obj.radius(I));
y = bsxfun(@times,y,obj.radius(I));
z = bsxfun(@times,z,obj.length(I));

u = [0; 0; 1];

start_points = obj.start_point(I,:);
ax = obj.axis(I,:);

raxis = (bsxfun(@cross,u,ax'))';
angle = acos(ax(:,3)./sqrt(sum(ax.^2,2)));

set_hold_off = ~ishold;

h = zeros(n,1);

if vector
    end_points = obj.end_point(I,:);
end

for i = 1:n
    
    if i == 2
        hold on
    end
   
    if any(raxis(i,:))
        X = [x(i,:)' y(i,:)' z(i,:)']*rotation_matrix(raxis(i,:)/norm(raxis(i,:)),angle(i))';
    elseif ax(i,3) == -1
        X = [x(i,:)' y(i,:)' -z(i,:)'];
    else
        X = [x(i,:)' y(i,:)' z(i,:)'];
    end

    % Reshape to size which is required by the surf-command.
    xi = reshape(X(:,1),s);
    yi = reshape(X(:,2),s);
    zi = reshape(X(:,3),s);

    % Transfer computed cylinder to the starting point of the
    % branch.
    xi = xi + start_points(i,1);
    yi = yi + start_points(i,2);
    zi = zi + start_points(i,3);

    % Draw the surface of the branch.
    if ~isempty(cdata)
        h(i) = surf(xi,yi,zi,cdata,'EdgeColor','none');
    else
        h(i) = surf(xi,yi,zi,'EdgeColor','none','FaceColor',[0.3 0.3 0.3]);
    end
    
    if vector
        hold on;        
        line = [start_points(i,1), start_points(i,2), start_points(i,3); ...
                end_points(i,1),end_points(i,2),end_points(i,3)];
        plot3(line(:,1),line(:,2),line(:,3),'-','color',color);
        plot3(line(end,1),line(end,2),line(end,3),'*','color',color);
    end
end

axis equal;

if set_hold_off
    hold off;
end

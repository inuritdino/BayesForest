function varargout = stem_branches(obj,varargin)


% Select cylinders that are in the trunk.
P = find(obj.is_trunk(1:obj.number_of_branches));

if not(isempty(P))
       
    % Process additional arguments.
    i = 1;
    
    % Default settings.
    north = [0 1 0];
    relative = true;
    radiusrelative = false;
    
    % Estimate tree axis direction as the mean of all its cylinders
    % axes.
    tree_axis = mean(obj.axis(P,:),1);

    while i <= size(varargin,2)

        if ischar(varargin{i})

            switch lower(varargin{i})
                case 'north'
                    north = varargin{i+1};
                    north = (north(:)/norm(north))';
                    i = i + 1;
                case 'relative'
                    relative = varargin{i+1};
                    i = i + 1;
                case 'axis'
                    tree_axis = varargin{i+1};
                    i = i + 1;
                case 'radiusrelative'
                    radiusrelative = true;
                otherwise
                    disp(['Ignoring unknown parameter: ' varargin{i}]);
            end
        end
        i = i + 1;
    end
    
    tree_axis = tree_axis/norm(tree_axis);

    % Stem children.
    try
        C = cat(1,obj.children{P});
    catch
        C = cat(2,obj.children{P});
        C = C(:);
    end
    
    % Properties of the children.
    sp = obj.start_point(C,:);
    ax = obj.axis(C,:);
    
    % Compute elevation and rotation relative to the parent.
    if relative
     
        j = 1;
        
        % Go through every stem cylinder with children.
        for i = 1:length(P)
            
            chi = obj.children{P(i)};
            n = length(chi);
            
            if n
                a = obj.axis(P(i),:);
                
                % Generate a coordinate system for the tree axis.
                if a(1) == 0 && a(2) == 0
                    M = eye(3);
                else
                    a1 = [a(2) -a(1) 0];
                    a1 = a1/norm(a1);
                    a2 = cross(a,a1);
                    a2 = a2/norm(a2);

                    M = [a1(:) a2(:) a(:)];
                end

                % Correct north direction.
                north = north*M;
                north = [north(1:2), 0]/norm(north(1:2));
                east = cross(a,north);

                M = [east(:) north(:) a(:)];
                
                ax(j:j+n,:) = ax(j:j+n,:)*M;
            end
            
        end
        
    % Use same parent axis direction for all child cylinders.
    else

        % Generate a coordinate system for the tree axis.
        if tree_axis(1) == 0 && tree_axis(2) == 0
            M = eye(3);
        else
            a1 = [tree_axis(2) -tree_axis(1) 0];
            a1 = a1/norm(a1);
            a2 = cross(tree_axis,a1);
            a2 = a2/norm(a2);

            M = [a1(:) a2(:) tree_axis(:)];
        end
        
        % Correct north direction.
        north = north*M;
        north = [north(1:2), 0]/norm(north(1:2));
        east = cross(tree_axis,north);
        
        M = [east(:) north(:) tree_axis(:)];
        
        % Project child axis into new coordinate system.
        ax = ax*M;
        
    end
    
    % Rotation beta.
    
    % Compute distance to Origin.
    d = sqrt(ax(:,1).^2 + ax(:,2).^2);

    % Compute angle of data points around the Origin.
    beta = asin(ax(:,2)./d);

    % Correct values for second and third quaters.
    q2 = (ax(:,1) < 0 & ax(:,2) >= 0);
    q3 = (ax(:,1) < 0 & ax(:,2) <  0);
    beta(q2) =  pi - beta(q2);
    beta(q3) = -pi - beta(q3);
    
    % Elevation alfa.
    alfa = atan(ax(:,3)./d);
    
    % Tree base height.
    [~,i] = min(obj.start_point(P,3));
    hb = obj.start_point(P(i),:)*tree_axis';
    
    % Height h.
    h = sp*tree_axis' - hb;
    
    r = obj.radius(C);
    
    % Scale with parent radius.
    if radiusrelative
        r = r./obj.radius(obj.parent(C));
    end
    
    % Sort child cylinders according to accending height.
    [h, J] = sort(h);
    beta = beta(J);
    alfa = alfa(J);
    r = r(J);
    
    % Select output according to number of output arguments.
    switch nargout

        % If one output, return matrix.
        case 1
            varargout{1} = [h, beta, alfa, r];

        % Return height and rotation values separately.
        case 2
            varargout{1} = h;
            varargout{2} = beta;
            
        % Return height, rotation and elevation values separately.
        case 3
            varargout{1} = h;
            varargout{2} = beta;
            varargout{3} = alfa;
            
        % Return height, rotation, elevation and radius values separately.
        case 4
            varargout{1} = h;
            varargout{2} = beta;
            varargout{3} = alfa;
            varargout{4} = r;
            
        % By default plot properties as functions of height.
        otherwise
            
            % Rotation.
            subplot(1,3,1);
            plot(beta,h,'o');
            xlabel('Rotation')
            ylabel('Height')
            
            % Elevation.
            subplot(1,3,2);
            plot(alfa,h,'o');
            xlabel('Elevation')
            ylabel('Height')
            
            % Radius.
            subplot(1,3,3);
            plot(r,h,'o');
            xlabel('Radius')
            ylabel('Height')
            
    end
else
    % If there are no cylinders in the tree classified as part of the stem,
    % raise an error.
    error('No stem in tree. Unable to analyze stem properties.')
end
function A = lin_1d_transforms(dim, nr, data, dir_gen_method, pr, pl) 
% Computes linear, 1d transforms a^T : R^dim -> R where ||a||=1 using NTM
% (or random projections or pca).
% - MJ

% input:
% dim, dimension of the transformation vector
% nr, number of these 1d transforms requested
% data, datapoints that must be given if one wants to use pca to determine 
% the directions
% dir_gen_method, -1 if compute directions using principal components 
% estimated from provided data, 0 if random directions, 1...5 quasi-MC methods
% pr, 1 if print info of the generated directions
% pl, 1 if plot the directions for visual inspection (only in 2d or 3d)
%
% output:
% A, the requested transforms, each row == one transformation vector


% ------ some debug settings -------
check_normalisation = 0;
check_parallel = 0;
% ------ end of settings -----------

if nargin < 6
    pl = 0;
end
if nargin < 5
    pr = 0;
end
if nargin < 4;
    dir_gen_method = 2; % default method: gp-2 set of points
end
if nargin < 3
    data = [];
end

A = [];
if nr <= 0 || dim <= 0
    disp('Error in 1d transformation function!');
    return;
end

%% generate random directions
if dir_gen_method == 0
    A = rand_sphere(dim,nr,pl);
    A = A';

%% compute principal components from data and use them as directions
elseif dir_gen_method == -1
    if ~isempty(data)
        A = pca(data'); % columns contain the directions by pca
        A = A';
    else
        disp('Error! Dataset was not given.');
        return;
    end

%% use the old code: uniform angles and compute rotation matrices
% elseif dir_gen_method == -2 && dim <= 3
%     thetas = give_thetas(nr, dim, 0);
%     nr = size(thetas,2);
%     A = NaN(nr,dim);
%     for i = 1:nr
%         R = rot_matrix(thetas(:,i), dim);
%         A(i,:) = R(1,:);
%     end

%% use the new code that is based on number-theoretic methods
% (pp. 167-170 in Number-theoretic methods in statistics book)
else    
    if dim == 1
        A = 1;
        
    elseif dim == 2 % as a special case
        angles = pi/2*(2*(1:nr)'-1)/nr;
        A = [cos(angles) sin(angles)];
        
    % different algorithms for odd and even dim
    elseif rem(dim,2) == 0 % even dimension
        s = dim;
        A = NaN(nr,s);
        ck = generate_NT_net(nr,s-1,dir_gen_method);
        %ck
        m = s/2;
        g = zeros(nr,m);
        g(:,m) = 1;
        for j = m-1:-1:1
            g(:,j) = g(:,j+1).*ck(:,j).^(1/j);
        end
        d = NaN(nr,m);
        for l = m:-1:2
            d(:,l) = sqrt(g(:,l) - g(:,l-1));
        end
        d(:,1) = sqrt(g(:,1));
        for l = 1:m
            A(:,2*l-1) = d(:,l).*cos(2*pi*ck(:,m+l-1));
            A(:,2*l) = d(:,l).*sin(2*pi*ck(:,m+l-1));
        end
    else % odd dimension 
        s = dim;
        A = NaN(nr,s);
        ck = generate_NT_net(nr,s-1,dir_gen_method);
        %ck
        m = (s-1)/2;
        g = zeros(nr,m);
        g(:,m) = 1;
        for j = m-1:-1:1
            g(:,j) = g(:,j+1).*ck(:,j).^(2/(2*j+1));
        end
        d = NaN(nr,m);
        for j = m:-1:2
            d(:,j) = sqrt(g(:,j) - g(:,j-1));
        end
        d(:,1) = sqrt(g(:,1));
%         A(:,1) = d(:,1).*(1-2*ck(:,m));
%         A(:,2) = d(:,1).*sqrt(ck(:,m).*(1-ck(:,m))).*cos(2*pi*ck(:,m+1));
%         A(:,3) = d(:,1).*sqrt(ck(:,m).*(1-ck(:,m))).*sin(2*pi*ck(:,m+1));
        % there is an error in the book (?), the algorithm should be as below 
        % since otherwise the norms of the direction vectors do not equal 1
        A(:,1) = d(:,1).*(1-2*ck(:,m));
        A(:,2) = 2*d(:,1).*sqrt(ck(:,m).*(1-ck(:,m))).*cos(2*pi*ck(:,m+1));
        A(:,3) = 2*d(:,1).*sqrt(ck(:,m).*(1-ck(:,m))).*sin(2*pi*ck(:,m+1));
        for l = 2:m
            A(:,2*l) = d(:,l).*cos(2*pi*ck(:,2*l));
            A(:,2*l+1) = d(:,l).*sin(2*pi*ck(:,2*l));
        end
    end
end

%% check normalisation - just in case
if check_normalisation
    A
    for i = 1:size(A,1)
        % print length of direction vectors!
        norm(A(i,:))
    end
end

%% leave out directions that are the same (shouldn't happen) or
%% parallel (which might be possible)
if check_parallel
    B = A; A = [];
    tol = 10*eps;
    for i = 1:nr
        is_parallel = 0;
        for j = i+1:nr
            if norm(B(i,:) + B(j,:)) < tol || norm(B(i,:) - B(j,:)) < tol
                is_parallel = 1;
                break;
            end
        end
        if ~is_parallel
            A = [A; B(i,:)];
        end
    end
    nr = size(A,1);
end

%% draw the samples for visual inspection in 2d or 3d case
if pl
    if dim == 2
        figure;
        plot(A(:,1),A(:,2),'*'); hold on; % points on circle
        t = linspace(0,2*pi,1000);
        plot(sin(t),cos(t),'r--'); % circle
        for i = 1:size(A,1)
            plot([-A(i,1), A(i,1)],[-A(i,2), A(i,2)],'b--'); % lines
        end
        axis square; hold off;
    elseif dim == 3
        figure;
        plot3(A(:,1),A(:,2),A(:,3),'*'); hold on;
        for i = 1:size(A,1)
            plot3([-A(i,1), A(i,1)],[-A(i,2), A(i,2)],[-A(i,3),A(i,3)],'b--'); % lines
        end
        k = 5;
        n = 2^k-1;
        [x,y,z] = sphere(n);
        surf(x,y,z,'FaceColor','red','EdgeColor','none');
        camlight left; lighting phong;
        alpha(0.5);
        axis square; hold off;
        grid on;
    elseif dim > 3
        disp('Dimension > 3, cannot plot directions!');
    end
end

if pr
    disp(['Total # of "directions": ', num2str(size(A,1))]);
end

end






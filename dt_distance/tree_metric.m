function d = tree_metric(data1, data2, settings)
% Compute 'distance' between two datasets, this wrapper is for testing the
% similarity/dissimilarity between forest trees
% - MJ

% INPUT:
% data1: first dataset, size must be dim*n, each column == one datapoint. 
% data1 can also be cell array where each cell contains dim*n matrix. 
% In this case data2 must be also cell array with equal length and each
% cell must contain datamatrix with correct size
%
% data2: second dataset, size must be dim*m, each column == one datapoint
%
% settings.stat1d: 1d "discrepancy" i.e. distance between the two datasets 
% in 1d:
% 1 == Kolmogorov-Smirnov i.e. sup-distance (default)
% 2 == Cramer von Mises i.e. L_2 distance
% 3 == Anderson-Darling, weighted L_2 distance
% 4 == Watson, modified Cramer von Mises distance
% 5 == L^1, L_1 distance
% 6 == L^2, basically the same as Cramer von Mises
% 7 == weighted L^2 basically the same as Anderson-Darling, 
% 8 == ecdf area^2, the area between empirical cdf's squared
% 9 == ecdf area, the area between empirical cdf's
% 10 == 1d energy distance
% -1 == NN distance
% -2 == Distance based on empirical moments of data
% -3 == Energy distance
% -4 == chi^2 type NN distance
%
% settings.directions: amount of line projections
% If data1 and data2 are cell arrays then it can be also cell array that 
% has the same cell array length as data1 and data2 and contains the 
% corresponding amounts of projections lines
%
% settings.norm_method: 1 if compute the average of all statistics
% corresponding the different directions, 2 if compute the maximum of those
%
% OUTPUT:
% d, the requested distance


% --- some settings/fixed values:
smooth_ecdf = 0; % 1 or 2 if smooth the empirical cdf for gradient based optimization
two_sample = 1; % 1 if compare two discrete sets of data, 0 if model cdf is known
d1_scaling = 1; % additional scaling
d2_scaling = 1; % additional scaling
max_directions = 1000; % upper bound for the number of the directions
% ---

if nargin < 3
    settings.stat1d = 1;
    settings.directions = 100;
    settings.norm_method = 1;
end

if iscell(data1) % multiple dataset pairs were given
    n = numel(data1); % # datasets
    d = 0;
    
    % set up the correct amounts of projection lines
    if ~iscell(settings.lines)
        settings.directions = num2cell(settings.directions*ones(n,1));
    end
    directions = NaN(n,1);
    for i = 1:n
        directions(i) = settings.directions{i};
        if directions(i) > max_directions
            disp(['Too many directions! Reduced the directions to ',...
                num2str(max_directions),'.']);
            directions(i) = max_directions;
        end
    end
    tot_directions = sum(directions);
    
    % go though the different, independent datasets
    for i = 1:n
        [dn_sum, dn_max] = dt_gof_statistic(data1{i}, data2{i}, [], [], ...
            directions(i), settings.stat1d, two_sample, smooth_ecdf);
        
        % take into account the amount of projection lines in the weights
        % of different datasets
        if settings.norm_method == 1
            d = d + d1_scaling * dn_sum * directions(i)/tot_directions;
        else
            d = max(d, d2_scaling * dn_max);
        end
    end
    
else % only one dataset pair was given
    directions = settings.directions;
    if directions > max_directions
        disp(['Too many directions! Reduced the directions to ',...
            num2str(max_directions),'.']);
        directions = max_directions;
    end
    
    [dn_sum, dn_max] = dt_gof_statistic(data1, data2, [], [], ...
        directions, settings.stat1d, two_sample, smooth_ecdf);
    
    if settings.norm_method == 1
        d = d1_scaling * dn_sum;
    else
        d = d2_scaling * dn_max;
    end
end

end






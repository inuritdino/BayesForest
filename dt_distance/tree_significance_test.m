function [H, p_value, p_dist] = tree_significance_test(data1, data2, settings)
% Compute the p-value for the hypothesis test that two datasets (trees) are
% the same 
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
% (This is the same setting as the 'sum_max'.)
%
% OUTPUT:
% H, 1 if H_0 hypothesis is true i.e. data1 and data2 come from the same
% distribution, 0 otherwise (if several datasets are provided, then the
% decision for each dataset is returned)
% p_value, the p-value for the hypothesis test (if several datasets are 
% provided, then all p-values are provided)
% p_dist, 1 - p_value
%
% NOTE:
% This function is not finished and tested thoroughfully but the 
% implemented algorithms etc. should work.


% --- some settings/fixed values:
alpha = 0.05; % confidence level of the test
max_directions = 1000; % upper bound for the number of the directions
% ---


if nargin < 3
    settings.stat1d = 1;
    settings.directions = 100;
    settings.norm_method = 1;
end

if iscell(data1) % multiple dataset pairs were given
    
    n = numel(data1); % # datasets
    p_value = NaN(n,1);
    H = NaN(n,1);
    
    % set up the correct amounts of projection lines
    if ~iscell(settings.directions)
        settings.directions = num2cell(settings.directions*ones(n,1));
    end
    
    % go though the different, independent datasets
    for i = 1:n
        directions = settings.directions{i};
        if directions > max_directions
            disp(['Too many directions! Reduced the directions to ',...
                num2str(max_directions),'.']);
            directions = max_directions;
        end
        p_value(i) = mc_permutation_test(data1{i}, data2{i}, settings.stat1d, ...
            directions, settings.norm_method);
    end

    % check the hypothesis rejection for each dataset separately for now...
    for i = 1:n
        if p_value(i) < alpha
            H(i) = 0;
        else
            H(i) = 1;
        end
    end
    
%     % use bonferroni correction to take multiple testing into account
%     for i = 1:n
%         if p_value(i) < alpha/n
%             H(i) = 0;
%         else
%             H(i) = 1;
%         end
%     end
    
    % todo: could sidak correction be used??? Or just count the distances
    % together as in one paper??? But the datasets are not of the same size
    % and not inpedendent -> ??? Should this approach be used with tree
    % data at all in the first place??? 
    
else % only one dataset pair was given
    
    directions = settings.directions;
    if directions > max_directions
        disp(['Too many directions! Reduced the directions to ',...
            num2str(max_directions),'.']);
        directions = max_directions;
    end
    p_value = mc_permutation_test(data1, data2, settings.stat1d, ...
        directions, settings.norm_method);
        
    if p_value < alpha
        H = 0;
    else
        H = 1;
    end
end

p_dist = 1 - p_value;

end





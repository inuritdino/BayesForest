function [D_sum, D_max, D_all] = dt_distance(samples_model, samples_exp, ...
    directions, stat1d, smoothing, scaling) 
% Wrapper for computing the "distance" between two (structural) datasets.
% The distance is calculated using the notion of distribution tomography.
% See elsewhere. 
% Author: Marko J?rvenp??.
%
% INPUT:
% samples_model, size dim*n, samples from the model, each column == one
% datapoint
% samples_exp, size dim*m, experimental data corresponding the model, each 
% column == one datapoint
% directions, number of directions (projections) for which the 1d statistics 
% are computed (default 100)
% stat1d, 1d "discrepancy" i.e. the distance between the two datasets in 1d:
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
% smoothing, 0 no smoothing, 1 if smooth the empirical cdfs for gradient 
% based optimization using kernel density estimation, 2 if smooth using 
% interpolation (default 0)
%
% OUTPUT:
% D_sum, average of 1d statistics computed to different directions (when stat1d > 0)
% D_max, max value of 1d statistics computed to different directions (when stat1d > 0)
% D_all, all the distances from all directions
%
% NOTE:
% Recommended settings: 
% 1) set directions to atleast 100; the higher the number the better D_sum
% and D_max approximate the theoretical distance functions but, unfortunately,
% computational cost increases linearly on the number of directions -> some
% compromise, like 100 directions, must be made 
% 2) Kolmogorov-Smirnov or Cramer von Mises are probably the best choices but 
% all distances should work; they just compare the datasets somewhat 
% differently; from the non-projection type distances the energy distance 
% is perhaps the most useful and it has nice theory behind it, unfortunately 
% it gets very slow if the dataset is large 
% 3) smoothing 2 seems to be  the best option for approximating the exact 
% value of the distance (smoothing 0 == exact value) with continuous function 
% when minimizing the distance with gradient-based methods. If one simply 
% wants to compute the distance then there is no reason to set smoothing to
% other value than 0. 
% 4) D_sum is the recommended option for the distance measure, d_max can be
% somewhat inaccurate unless the amount of directions is high, like >1000. 


% --- some settings/fixed values:
two_sample = 1; % 1 if compare two finite sets of data, 0 if model cdf is known
d1_scaling = 1; % additional scaling
d2_scaling = 1; % additional scaling
max_directions = 100000; % upper bound for the number of the directions
if(nargin < 6)
    scaling = 0;
end
if(scaling)
    use_normalised_stat = 2; % 1 no scaling, 2 scaling ~ sqrt(n) or ~ sqrt(n1*n2/(n1+n2)) 
    % so that the value of the distance does not approximately depend on sample size
else
    use_normalised_stat = 1;
end
% ---

% check input
dim = size(samples_model,1);
if nargin < 5
    smoothing = 0;
end
if nargin < 4
    % default 1d-statistic:
    stat1d = 1;
end
if nargin < 3
    % default amount of projection lines:
    directions = 100;
    if dim == 1
        directions = 1;
    end
end

% just in case, check the amount of directions
if directions > max_directions
    disp(['Too many directions! Reduced the directions to ',...
        num2str(max_directions),'.']);
    directions = max_directions;
end

% compute the statistic!
[dn_sum, dn_max, dn_all] = dt_gof_statistic(samples_exp, samples_model, [], [], ...
    directions, stat1d, two_sample, smoothing, use_normalised_stat);

% the average of the discrepancies to different directions:
D_sum = d1_scaling * dn_sum; 
% the maximum discrepancy of the different directions:
D_max = d2_scaling * dn_max; 
% Distance/statistic values from all directions
D_all = dn_all;
end





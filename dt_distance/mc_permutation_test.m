function [p_value, d, p_err, cv_err] = mc_permutation_test(data1, data2, stat1d, ...
    directions, sum_max)
% Test the hypothesis that data1 and data2 come from the same unknown
% distribution using some distance based on computing 1d test statistics to
% different 1d projections. Other distances can also be used. 
% Monte Carlo permutation test (or actually bootstrap test) is used to 
% estimate the p-value because the test depends on the unknown
% distribution. These tests are only distribution-free in 1d case with e.g. 
% K-S or CvM distance.
% - MJ

% input:
% data1, the first dataset dim*n matrix
% data2, the second dataset dim*m matrix
% stat1d, 1d statistic
% directions, amount of projection lines
% sum_max, 1 if average the 1d statistics, 2 if maximum 1d discrepancy
%
% output:
% p_value, (approximative) p-value of the hypothesis test
% d, distances generated from pooled data with Monte Carlo
% p_err, (approximative) standard error in p-value
% cv_err, (approximative) coefficient of variation in p-value


% ****** settings ***************

if nargin < 3
    % stat1d: 1 == ks, 2 == cvm, 3 == anderson-darling, 4 == Watson, 5 == L^1,
    % 6 == L^2, 7 == weighted L^2, 8 == ecdf area^2, 9 == ecdf area, 
    % -1 == NN, -2 == SMM, -3 == e-distance, -4 == chi^2 NN
    stat1d = 1;
    directions = 100;
    sum_max = 1;
end

pl = 1; % 2 if plot/print everything, 1 if print only the results, 
% 0 print nothing

% Monte Carlo error control:
adapt_err_control = 1;
std_tol = 0.025;
cv_tol = 0.5; % see p. 211 in the 1993 Bootstrap book by Efron & Tibshirani
max_draws = 500; % maximum amount of MC draws in any case

% ****** settings end here ******


two_sample = 1;
use_normalised_stat = 1;
n = size(data1,2); 
m = size(data2,2); 
N = n + m;
[d0_sum, d0_max] = dt_gof_statistic(data1, data2, [], [], directions, ...
    stat1d, two_sample, 0, use_normalised_stat);
if sum_max == 1
    d0 = d0_sum;
else
    d0 = d0_max;
end
% should this value be included in the MC sample?


if adapt_err_control
    %mc_draws = min(ceil(0.25/std_tol^2), max_draws);
    mc_draws = max_draws;
else
    mc_draws = max_draws;
end

d = NaN(mc_draws,1);
pooled_data = [data1 data2];
j = 0;
for i = 1:mc_draws
    if pl >= 2
        disp(['iteration: ',num2str(i)]);
    end
    j = j + 1;
    
    % split the pooled data into two sets
    p = randperm(N);
    pooled_data_i = pooled_data(:,p);
    data1 = pooled_data_i(:,1:n);
    data2 = pooled_data_i(:,n+1:N);
    
    % compute distance and add it to the result
    [d_sum,d_max] = dt_gof_statistic(data1, data2, [], [], ...
        directions, stat1d, two_sample, 0, use_normalised_stat);
    
    if sum_max == 1
        d(i) = d_sum;
    else
        d(i) = d_max;
    end
    
    % compute p-value based on all MC draws this far
    p_value = sum(d(1:i) > d0)/length(d(1:i));
    
    % approximate the current standard error and coefficient of variation
    p_err = sqrt(p_value*(1 - p_value)/i);
    cv_err = sqrt(((1 - p_value)/p_value)/i);
    
    % stop computing if the approximated MC error is low enough
    if i >= 100 && adapt_err_control && p_err <= std_tol && cv_err <= cv_tol
        break;
    end
end
d = d(~isnan(d));

if adapt_err_control && pl >= 1
    disp(['MC draws performed: ', num2str(j)]);
    disp(['Estimated p-value: ', num2str(p_value)]);
    disp(['Estimated standard error in p-value: ', num2str(p_err)]);
    disp(['Estimated coefficient of variation in p-value: ', num2str(cv_err)]);
end

if pl >= 2
    figure
    %hist(d)
    [f,xi] = ksdensity(d);
    plot(xi,f,'-b'); hold on;
    plot(d0,0,'*b'); hold off;
end

end





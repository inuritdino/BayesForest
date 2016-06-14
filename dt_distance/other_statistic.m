function [Z, D] = other_statistic(meas_cdf, real_cdf, two_sample, v)
% Compute some obvious 1d statistics between sample and model cdf or 
% between two sets of samples. In some cases these correspond some 
% statistically motivated distances (L^2 ~ Cramer von Mises, 
% weighted L^2 ~ Anderson-Darlingm "squared area" statistic ~ energy distance")
% - MJ

% input:
% meas_cdf, the empirical cdf of the first dataset (if two_sample == 0), 
% or the first set of datapoints in increasing order (if two_sample == 1)
% real_cdf, the values of the model cdf (if two_sample == 0),
% the second set of datapoints also in increasing order (if two_sample == 1)
% two_sample, 1 if two sample version of the test, 0 if one sample version 
% i.e. the model cdf is known and is given as input
%
% output:
% Z, D, both are the same values of the statistic
%
% Note:
% 1) this function is for testing only, no reason to compute these instead
% of CvM etc.
% 2) if one sample test then meas_cdf and real_cdf are of same size, the 
% points correspond to each other and cdf's are normalised to one
% 3) if two sample test then meas_cdf and real_cdf are the points and not the 
% cdf values (!) and are not necessarily of the same size
% 4) no special scaling is used for the values of the statistics i.e. Z == D
% 5) weighted L^2 is not symmetric since it used cdf (or the second
% dataset) for weighting


if nargin < 4
    v = 2;
end
if nargin < 3
    two_sample = 0;
end

Z = Inf;
D = Inf;
if ~two_sample
    if v == 1 % L^1 norm
        D = sum(abs(meas_cdf - real_cdf));
        Z = D;
    elseif v == 2 % L^2 norm 
        D = sum((meas_cdf - real_cdf).^2);
        Z = D;
    elseif v == 3 % weighted L^2 norm 
        tol = 1e-12;
        % tol is used to make sure that no division with zero happens
        D = sum((meas_cdf - real_cdf).^2./(max(real_cdf.*(1-real_cdf), tol)));
        Z = D;
    elseif v == 4
        % computing this statistic would require (numerical) integration!
        disp('Requested 1d statistic can be computed only in 2-sample case!');
    elseif v == 5
        % computing this statistic would require (numerical) integration!
        disp('Requested 1d statistic can be computed only in 2-sample case!');
    else
        disp('Requested invalid 1d statistic!');
    end
else
    % two-sample case
    n1 = length(meas_cdf);
    n2 = length(real_cdf);
    N = n1 + n2;
    c = sqrt(n1*n2/N);
    
    % pool the samples 
    all = [meas_cdf, real_cdf];
    all = sort(all);
    ecdf1 = emp_cdf(meas_cdf, all);
    ecdf2 = emp_cdf(real_cdf, all);
    
    if v == 1 % L^1 norm
        Z = sum(abs(ecdf1 - ecdf2));
        D = Z;
    elseif v == 2 % L^2 norm 
        % Note! This is the same as (unscaled) CvM 1d statistic!
        Z = sum((ecdf1 - ecdf2).^2);
        D = Z;
    elseif v == 3 % weighted L^2 norm 
        tol = 1e-12;
        % tol is used to make sure that no division with zero happens
        Z = sum((ecdf1 - ecdf2).^2./(max(ecdf2.*(1-ecdf2), tol)));
        D = Z;
    elseif v == 4
        % This statistics is the "squared area" between empirical cdfs
        % Note! This is the same as the energy distance, only scaling is
        % different!
        Z = sum((all(2:end) - all(1:(end-1))) .* ...
            (ecdf1(1:(end-1)) - ecdf2(1:(end-1))).^2);
        D = Z;
    elseif v == 5
        % this statistics is the "area" between empirical cdfs
        Z = sum((all(2:end) - all(1:(end-1))) .* ...
            abs(ecdf1(1:(end-1)) - ecdf2(1:(end-1))));
        D = Z;
    else
        disp('Requested invalid 1d statistic!');
    end
end

end





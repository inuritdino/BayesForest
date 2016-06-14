function [Z, D] = ks_statistic(meas_cdf, real_cdf, two_sample) % ok
% computes K-S statistic in one dimensional case between one dataset and the
% known model cdf or between two datasets
% - MJ

% input:
% meas_cdf, the empirical cdf of the first dataset (if two_sample == 0), 
% or the first set of datapoints (in increasing order!) (if two_sample == 1)
% real_cdf, the values of the model cdf (if two_sample == 0),
% the second set of datapoints (also in increasing order!) (if two_sample == 1)
% two_sample, 1 if two sample version of the test, 0 if one sample version 
% i.e. the model cdf is known and is given as input
% 
% output:
% Z, value of ks statistic so that it belongs to interval [0,1]
% D, as above but scaled with sqrt(n) or similarly in 2 sample case
%
% Note:
% 1) if one sample test then meas_cdf and real_cdf are of same size, the 
% points correspond to each other and cdf's are normalised to one
% 2) if two sample test then meas_cdf and real_cdf are the points and not the 
% cdf values (!) and are not necessarily of the same size 

if nargin < 3
    two_sample = 0;
end

if ~two_sample
    % one sample case:
    n = length(meas_cdf);
    meas_cdf2 = [0, meas_cdf(1:end-1)];
    d1 = max(abs(meas_cdf - real_cdf));
    d2 = max(abs(meas_cdf2 - real_cdf));
    Z = max([d1,d2]);
    D = sqrt(n) * Z;
else
    % two sample case:
    n1 = length(meas_cdf);
    n2 = length(real_cdf);
    N = n1 + n2;
    all = [meas_cdf, real_cdf];
    all = sort(all); 
    m1 = emp_cdf(meas_cdf, all);
    m2 = emp_cdf(real_cdf, all);
    % unlike in 1 sample test, it is enough to 
    % check only either one of the differences
    d1 = max(abs(m1 - m2));
    Z = d1;
    D = sqrt(n1*n2/N) * Z;
end

end





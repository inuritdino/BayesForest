function [v, stdev] = total_var(sa)
% Compute the "total variance" i.e. the sum of variances (diagonal elements)
% in a covariance matrix estimated using samples sa
% - MJ

% input:
% sa, samples that are used for estimating the cov matrix, size dim x n
%
% output:
% v, requested total variance
% stdev, total standard deviation

v = sum(var(sa,[],2));
stdev = sqrt(v);

end





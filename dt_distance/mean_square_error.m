function [mse, avg_bias] = mean_square_error(true_param, samples)
%Compute the mean square error (and average bias) between true value 
%true_param and some given estimates for it
% - MJ

% input:
% true_param, true parameter (as scalar or vector)
% samples, some computed estimates that are to be compared to true_param,
% each column must be one estimate
%
% output:
% mse, mean square error computed using samples
% avg_bias, average bias to true value computed using samples

% check input etc.
true_param = true_param(:);
m = length(true_param);
n = size(samples,1);
N = size(samples,2);

if m ~= n
    disp('Error! Sizes do not match.');
    mse = NaN;
    return;
end

% compute the 'distance' (which should be close to 0 if the estimate is 
% unbiased and sample size is large enough)
M = samples - repmat(true_param, 1, N);
d = sqrt(sum(abs(M).^2,1));
mse = mean(d);

% compute (average) bias
avg_bias = mean(M,2);

end





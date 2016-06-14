function resampled_data = bootstrap_resample(data,m)
% Draw a bootstrap resample (i.e. sample from the original data with 
% replacement) with m points. Some points of the original data can thus 
% appear several times in the resulting data.

% input:
% data, datapoints as dim*n matrix i.e. one column == one data point
% m, size of the resampled data, often the same as the size of the given 
% data
%
% output:
% resampled_data, resampled data points

if nargin < 2
    m = size(data,2);
end
m = max(0,m);
m = floor(m);

n = size(data,2);
inds = rnd_categorical(m, ones(1,n));
resampled_data = data(:,inds);

end





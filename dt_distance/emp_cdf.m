function z = emp_cdf(data, points) % ok
% Computes 1d empirical cdf using vector 'data' and returns the values as 
% given in vector 'points'
% - MJ

% input:
% data: row vector containing the data points (or matrix, each row computed
% separately)
% points: evaluates the empcdf at these points and return these values. If
% points vector is not given, then computes the values at given data points
%
% output:
% z: results returned as a row vector (or matrix)
%
% Note:
% Input data and points are assumed to be ordered


if nargin < 2 || isempty(points)
    points = data;
end

dim = size(data,1);
npoints = size(points,2);
ndata = size(data,2);
z = NaN(dim,npoints);

% this code also handles values that are possibly requested at different 
% points than the data
for i = 1:dim
    edges = [-inf, points(i,:).*(1+sign(points(i,:))*eps), inf]; 
    % 0 does not work quite right 
    % -> replace this anyway with better (faster) code later
    % Check for monotonicity
%     if(any(diff(edges) < 0.0))
%         if(any(isnan(points(i,:))))
%             disp(points(i,:));
%             disp('NaN points!');
%         end
%         error('Error:EMP_CDF:Non-monotonical values!');
%     end
    bincounts = histc(data(i,:), edges);
    sumcounts = cumsum(bincounts)/ndata;
    z(i,:) = sumcounts(1:end-2);
end

end






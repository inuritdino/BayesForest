function d = nn_distance(data1,data2,k,dist1d,chi2variant)
%Computes the 'nearest neigbours distance' (Shilling 1986) or chi^2-type 
%nearest neighbours distance between datasets data1 and data2 using k 
%nearest neighbours
% - MJ

% input:
% data1, the first dataset, each column is data point
% data2, the second dataset, each column is data point
% k, amount of nearest neighbours to compare, default == 3
% dist1d, which distance to use between the vectors, e.g. Euclidean == 1
% variant, 1 if compute the chi^2 type variant of NN-distance in which the
% proportions of the points in each NN-ball is compared
%
% output:
% d, the distance to be computed
%
% note:
% Ties in data are handled by taking the first value found (and not by choosing
% it randomly fro all distances with the same value)


if nargin < 5
    chi2variant = 0;
end
if nargin < 4
    dist1d = 1;
end
if nargin < 3
    k = 3;
end

% sizes of data
n1 = size(data1,2);
n2 = size(data2,2);
N = n1 + n2;

% distances between data points
all_data = [data1 data2]'; % pooled data

% find k nearest neighbours for each point in the pooled sample
if dist1d == 1
    nn_ind = knnsearch(all_data,all_data,'K',k + 1);
elseif dist1d == 2
    e = 1e-6;
    distfun = @(x,y) -log(e + col_2norm((repmat(x,size(y,1),1) - y)'))'; 
    nn_ind = knnsearch(all_data,all_data,'K',k + 1,'Distance',distfun);
else
    disp('Incorrect distance!');
    d = NaN;
    return;
end

% ignore the same point to be compared since this distance is 0, of course!
nn_ind = nn_ind(:,2:end); 

% check whether they come from sample 1 or 2 and compute the final distance
% measure
if ~chi2variant
    d = sum(sum(nn_ind(1:n1,:) <= n1));
    d = d + sum(sum(nn_ind(n1+1:N,:) >= n1 + 1));
    d = 1/(N*k) * d;
else
    d = sum((sum(nn_ind <= n1,2)/k - sum(nn_ind > n1,2)/k).^2);
    d = 1/N * d;
end

end






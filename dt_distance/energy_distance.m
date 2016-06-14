function d = energy_distance(data1,data2,dist1d)
%Computes the 'energy-distance' (Szekely et al. 2006) between datasets
%data1 and data2
% - MJ

% input:
% data1, the first dataset, each column is data point
% data2, the second dataset, each column is data point
% dist1d, which distance to use between the vectors, e.g. Euclidean == 1
%
% output:
% d, the distance to be computed
%
% note:
% Ties in data are handled by taking the first value found (and not by choosing
% it randomly fro all distances with the same value).


if nargin < 3
    dist1d = 1;
end

%data1,data2

% sizes of data
n1 = size(data1,2);
n2 = size(data2,2);
N = n1 + n2;

% distances between data points
all_data = [data1 data2]; % pooled data
if dist1d == 1
    d = pdist(all_data');
elseif dist1d == 2
    e = 1e-6;
    distfun = @(x,y) -log(e + col_2norm((repmat(x,size(y,1),1) - y)'))'; 
    d = pdist(all_data',distfun);
    %...
else
    disp('Incorrect distance!');
    d = NaN;
    return;
end
D = squareform(d);

% compute the distance!
d12 = sum(sum(D(1:n1,n1+1:N)));
d1 = sum(sum(D(1:n1,1:n1)));
d2 = sum(sum(D(n1+1:N,n1+1:N)));
d = n1*n2/N * (2/(n1*n2)*d12 - 1/n1^2*d1 - 1/n2^2*d2);

end







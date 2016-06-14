function [] = test_different_distances()
% Test the distances and their time consumption
% - MJ


% ---- settings -------

use_profiler = 0;

dim = 10; 
n = 1000000; 
nr_directions = 100;
smooth_ecdf = 0;

%stat1ds = [-4,-3,-2,-1,1,2,3,8,9,10]; % all methods
stat1ds = [1,2,3]; % useful pp-methods only

if n >= 10000 % ignore methods that are slow for >10000 points!
    stat1ds = stat1ds(stat1ds ~= -3 & stat1ds ~= 10); 
end
if smooth_ecdf % ignore methods that can't be computed in 2-sample case
    stat1ds = stat1ds(stat1ds ~= 8 & stat1ds ~= 9); 
end
% ---------------------

% methods
dists_nr = [-4,-3,-2,-1,1,2,3,8,9,10];
dists = {'chi^2-NN','e-distance','MM','NN','PP-KS','PP-CvM','PP-AD',...
    'PP-area^2','PP-area','PP-energy'};

% generate the data
my_seed = 1;
rng(my_seed);

data1 = randn(dim,n);
data2 = randn(dim,n);

if use_profiler
    profile on
end

% compute distances
for i = stat1ds
    disp(['distance ', num2str(i), ': ', dists{find(dists_nr==i)}]);
    tic;
    dt_gof_statistic(data1, data2, [], [], nr_directions, i, 1, smooth_ecdf)
    toc
    disp('---------------------------------');
end


% test also pdist that is used for energy-distance
if n <= 10000
    X = [data1,data2]';
    tic;
    d = pdist(X);
    toc, size(d), n*(n-1)/2
end


if use_profiler
    profile viewer
end





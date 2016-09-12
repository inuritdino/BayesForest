function [Davg, Ds, n, D_all] = optim_avg_distance(set1,set2,varargin)
% Calculate average (distribution tomography) distance (dt_distance) for 
% any number of pairs of sub-data sets within two sets SET1 and SET2.
% USAGE:
%       [Davg,Ds,n,Dall] = OPTIM_AVG_DISTANCE(SET1,SET2,...)
%
% SET1/SET2: if number of data sets > 1, the cell array of these data sets,
% otherwise, can be a matrix
% Number of sub-sets in SET1 and SET2 must be the same. The dt_distance()
% compares pairs of sub-sets giving a distance value for each pair.
%
% Davg: average distance between many sub-set pairs
% Ds: array of individual distances between each sub-set
% n: 2xN matrix with number of samples in each sub-set of N sub-sets in
%   SET1/SET2
% Dall: all NxNdir matrix with all individual 1D metrics obtained in
%   dt_distance by its characteristic transformation (see dt_distance).
%   Ndir is the number of directions to perform the distribution
%   tomography.
% VARARGIN:
%   'w': next element determines the weghts for the sub-sets, the weighted
%       average is then computed.
%   'dirs': next element is number of directions for the distribution
%       tomography (dt_distance).
%   'stat': 1D statistic to measure, a numerical value (see dt_distance for
%       values).
%   'smooth': next element is boolean, whether to perform smoothing on the
%       1D cumulative distribution.
%   'scaling': next element determines whether to do scaling of the
%       distance by the number of samples.
%   'subsets': the array of indices between 1 and N indicating which
%       sub-sets to choose for the analysis.
%   'max': calculate maximum of the 1D statistics instead of the default
%       average.
%
% See also dt_distance


%% Initials
% Look up this at dt_distance() function
w = [];
ndirs = 100;
stat1d = 1;
smooth = 0;
scaling = 0;
subsets = [];
max_dist = false;
multiobj = false;

%% Process the input
tf = strcmp('w',varargin);
if(find(tf))
    w = varargin{find(tf)+1};
end
tf = strcmp('dirs',varargin);
if(find(tf))
    ndirs = varargin{find(tf)+1};
end
tf = strcmp('stat',varargin);
if(find(tf))
    stat1d = varargin{find(tf)+1};
end
tf = strcmp('smooth',varargin);
if(find(tf))
    smooth = varargin{find(tf)+1};
end
tf = strcmp('scaling',varargin);
if(find(tf))
    scaling = varargin{find(tf)+1};
end
tf = strcmp('subsets',varargin);
if(find(tf))
    subsets = varargin{find(tf)+1};
end
tf = strcmp('max',varargin);
if(find(tf))
    max_dist = true;
end
tf = strcmp('multiobj',varargin);
if(find(tf))
    multiobj = true;
end

%% Make cell arrays from data sets
% Normally it should happen when there is only one data set in set1 or set2
if(~isa(set1,'cell'))
    set1 = {set1};
end
if(~isa(set2,'cell'))
    set2 = {set2};
end

%% Dimension check
empty1 = all(cellfun(@isempty,set1));
empty2 = all(cellfun(@isempty,set2));
N1 = length(set1);
N2 = length(set2);
% Form the array of weights
if(isempty(w))
    w1 = repmat(1/N1,1,N1);% equal weights by default
    w2 = repmat(1/N2,1,N2);% equal weights by default
else% rescale the weights, which were supplied by the user
    w = w(:)./sum(w(:));% normalize the weights
    w1 = w;
    w2 = w;
end

% Both data sets are empty --> ERROR
if(empty1 && empty2)
    error('Error: Both data sets are EMPTY.');
end

% Is one data set empty?
% Empty data set indicates some problem ==> maximum distance
Ds = [];
n = [];
if(empty1 || empty2)
    fprintf('Warning: empty set. D = 1.0\n');
    % Set max distance
    if(~scaling)
        Davg = 1.0;
    else
        % when scaling is present the prefactor is the "max distance",
        % which is sqrt(n^2/(2*n)), n - number of samples in the non-empty
        % data set. This is equivalent to make a complementary data set
        % with n samples of Inf's instead of the empty data set.
        % The actual distance is 1.0, but we need to calculate the
        % prefactor.
        if(empty1)
            N = N2;
            w = w2;
            Davg = 0.0;
            if(numel(w) ~= N)
                error('Error: weights dimension does not match data.');
            end
            for ii = 1:N
                n = size(set2{ii},2);
                Davg = Davg + w(ii)*sqrt(n^2/(2*n));
            end
        elseif(empty2)
            N = N1;
            w = w1;
            if(numel(w) ~= N)
                error('Error: weights dimension does not match data.');
            end
            Davg = 0.0;
            for ii = 1:N
                n = size(set1{ii},2);
                Davg = Davg + w(ii)*sqrt(n^2/(2*n));
            end
        end
    end
    % must return to interrupt the calculation
    return;
else% There is none empty set, check the dimensions
    if(N1 ~= N2)
        error('Error: dimensions of SET1 and SET2 do not match.');
    end
    % Dimensions are Ok, proceed with calculation of Davg.
    Davg = 0.0;
    N = N1;% N1 == N2
    w = w1;% w1 == w2
    if(numel(w) ~= N)
        error('Error: weights dimension does not match data.');
    end
end

%% Collect distances
Ds = zeros(1,N);
D_all = zeros(N,ndirs);
n = zeros(2,N);
for ii = 1:N
    if(~isempty(subsets) && isempty(find(subsets == ii,1)))
        continue;
    end
    if(~isempty(set1{ii}) && ~isempty(set2{ii}))
        n(:,ii) = [numel(set1{ii}); numel(set2{ii})];
        % compare the data with modelled data
        if(max_dist)
            [~,Ds(ii),D_all(ii,:)] = dt_distance(set1{ii},set2{ii},ndirs,stat1d,smooth,scaling);
        else
            [Ds(ii),~,D_all(ii,:)] = dt_distance(set1{ii},set2{ii},ndirs,stat1d,smooth,scaling);
        end
%         if(scaling)
%             % set2 is DATA and FIXED during optimization.
%             inf_set = Inf(1,n(2,ii));
%             Dmax = dt_distance(inf_set,set2{ii},ndirs,stat1d,smooth,scaling);
%             Ds(ii) = Ds(ii) + abs(n(1,ii)-n(2,ii))/n(2,ii).*Dmax.*sqrt(n(2,ii)/2);
%         end
        Ds(ii) = w(ii)*Ds(ii);% weight
        Davg = Davg + Ds(ii);% add up
    else
        Ds(ii) = w(ii);
        Davg = Davg + Ds(ii);% 1.0 max score
        n(:,ii) = [numel(set1{ii}); numel(set2{ii})];
    end
end

if(multiobj)% A trick to make the multiobjective outcome being the each U-distance
    Dswap = Davg;
    Davg = Ds;
    Ds = Dswap;
end
end
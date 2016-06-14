function [dn_sum, dn_max, dn_all] = dt_gof_statistic(samples, samples2, dens, params, ...
    directions, stat1d, two_sample, smooth_ecdf, use_normalised_stat) 
% Compute a distance between two datasets. It is assumed that 
% 1) the model cdf 'dens' with parameters 'params' is given (currently 
% this works only in the case of Gaussian i.e. dens == 1 or GM distribution 
% i.e. dens == 5) or
% 2) the model is characterized only by second set of samples 'samples2'
% which must be provided and the cdf information is unnecessary
% - MJ

% INPUT:
% samples: experimental i.e. measured data
%
% samples2: another set of data simulated from the model (in the case model 
% cdf is not known i.e. two_sample == 1)
%
% dens: if two_sample == 1 then 1 gaussian, 2 GM
% params: current set of parameters for the known model above
%
% directions: #directions (can also be computed beforehand and in this case 
% this is a matrix of the directions)
%
% stat1d: 1d statistic (1d 'distance') as below:
% 1 == Kolmogorov-Smirnov i.e. sup-distance (default)
% 2 == Cramer von Mises i.e. L_2 distance
% 3 == Anderson-Darling, weighted L_2 distance
% 4 == Watson, modified Cramer von Mises distance
% 5 == L^1, L_1 distance
% 6 == L^2, basically the same as Cramer von Mises
% 7 == weighted L^2 basically the same as Anderson-Darling, 
% 8 == ecdf area^2, the area between empirical cdf's squared
% 9 == ecdf area, the area between empirical cdf's
% 10 == 1d energy distance
% -1 == NN distance
% -2 == Distance based on empirical moments of data
% -3 == Energy distance
% -4 == chi^2 type NN distance
%
% two_sample: 0 if the model cdf is known (and in this case dens must be given)
% or 1 if second set of samples is given and model cdf is not to be used
%
% smooth_ecdf: 0, no smoothing i.e. compute the usual statistic, 1 if 
% estimate the model cdf from the samples2 using kernel density estimation 
% and proceed as if model cdf would be known, 2 if use interpolation for 
% estimating the cdf (This is useful when using the distance in 
% gradient-based optimization.) 
%
% use_normalised_stat: 1 no scaling, 2 scaling ~ sqrt(n) or ~ sqrt(n1*n2/(n1+n2))
%
%
% OUTPUT:
% dn_sum, average of the 1d statistics computed into different directions
%
% dn_max, maximum of the 1d statistics computed into different directions
%
% NOTE: 
% 1) This is not random but deterministic function!!!
% 2) If NN, SMM or energy distance is computed then dn_sum == dn_max and the 
% input values for directions, stat1d etc. are not used.
% 3) Recommended choice for smoothing is setting 2 i.e. interpolation method, 
% setting 1 i.e. kernel method can be inaccurate.


% extra settings
kernel_type = 1; % 1 gaussian, 2 laplace, 3 epanechnikov 
interp_type = 2; % 1 linear interpolation, 2 spline/cubic interpolation
use_loc_kernel = 0; % 1 if use local bandwidth in kernel case (currently slow)
moments2match = 2; % which moments to match (for SMM only)
dist1d = 1; % distance for NN and Energy distances, 1 == 'euclidean'
k_nn = 3; % amount of nearest neighbours for NN
k_nn_old = 10; % amount of nearest neighbours for chi^2-NN
% --------------

if nargin < 9
    use_normalised_stat = 1;
end
if nargin < 8
    smooth_ecdf = 0;
end
if nargin < 7
    two_sample = 0;
end

dn_sum = NaN; 
dn_max = NaN;
dim = size(samples,1); % dimension of the density


%% extract parameters if exact model cdf is used
if ~two_sample
    if dens == 1 % gaussian
        mu = params{1};
        sigma = params{2};
    elseif dens == 2 % GM
        w = params{1};
        nr_mixt = length(w);
        mus = params{2};
        sigmas = params{3};
    else
        disp('Error! Formulas for the cdf are not coded for given distribution.');
        return;
    end
end


%%-------------------------------------------------------------------------
%% SOME STATISTICALLY MOTIVATED DISTANCES BETWEEN TWO DATASETS:
%%-------------------------------------------------------------------------

%% Nearest-neighbour method as in Schilling 1986 (does only work for two-sample case)
if stat1d == -1 && two_sample
    dn_sum = nn_distance(samples,samples2,k_nn,dist1d);
    dn_max = dn_sum;
    return;
end


%% simulated/generalized method of moments (SMM/GMM)
% idea: match the (empirical) moments between two datasets or between one 
% data set and known model cdf! 
% NOTE: Those moments that should be matched depends on the particular problem
% NOTE2: Moments could also be weighted but no weighting is implemented here
% NOTE3: Comparing the covariance matrices could be done differently
if stat1d == -2
    
    % exact formulas for the moments of the model are known (currently only
    % gaussian or GM case)
    if ~two_sample && (dens == 1 || dens == 5) 
        % this simplifies to MLE:
        dn_sum = norm(mean(samples,2) - mu,2)^2;
        if moments2match == 2
            dn_sum = dn_sum + norm(cov(samples') - sigma,2)^2;
        end
    end
    
    % moments are estimated from data for both datasets
    if two_sample
        dn_sum = norm(mean(samples,2) - mean(samples2,2),2)^2;
        if moments2match == 2
            dn_sum = dn_sum + norm(cov(samples') - cov(samples2'),2)^2;
        end
    end
    
    dn_max = dn_sum;
    return;
end


%% Energy distance as in Szekely & Rizzo 2004 (does only work for two-sample case)
if stat1d == -3 && two_sample
    dn_sum = energy_distance(samples,samples2,dist1d);
    dn_max = dn_sum;
    return;
end


%% (Another variant of) nearest-neighbour method (does only work for two-sample case)
% idea: for any measurement point, take t nearest points and compute the
% fractions of the points in these cells and compute the sum of these
% differences
if stat1d == -4 && two_sample
    t_nn = min([size(samples,2),size(samples2,2),k_nn_old]); 
    dn_sum = nn_distance(samples,samples2,t_nn,dist1d,1);
    dn_max = dn_sum;
    return;
end


%%-------------------------------------------------------------------------
%% THE METHOD OF 1D PROJECTIONS I.E. OUR DISTANCE METRIC:
%% a.k.a. projection pursuit distance
%%-------------------------------------------------------------------------

%% compute rotations if they are not given as input
if dim > 1 && isscalar(directions)
    Rs = lin_1d_transforms(dim, directions, samples);
else
    Rs = directions;
end
directions = size(Rs,1); % make sure # directions is ok
d = NaN(directions,1);

if smooth_ecdf
    two_sample = 0; % model cdf is estimated from the data and is thus 'known'
end
    
%% loop for computing the 1d statistics into different directions
for j = 1:directions 
    
    %% rotate first set of samples to the current direction
    R = Rs(j,:);
    x = R * samples;
    x = sort(x);
    
    %% the case where data is compared to other data set
    if two_sample 
        
        % rotate second set of samples to the current direction
        x2 = R * samples2;
        x2 = sort(x2);
        exact_cdf = x2;
        e_ecdf = x;
    end
    
    %% the case where data in compared to known model cdf
    if ~two_sample
        
        % compute empirical cdf of data (it is needed for some 1d stats)
        if stat1d == 1 || (stat1d >= 5 && stat1d <= 9)
            e_ecdf = emp_cdf(x);
        else
            e_ecdf = [];
        end
        
        % 1) estimate cdf from the second data set (using either kernel 
        % density estimation or interpolation)
        if smooth_ecdf 
            
            x2 = R * samples2;
            x2 = sort(x2);
            if smooth_ecdf == 1
                exact_cdf = ksdensity_cdf(x, x2, kernel_type, use_loc_kernel);
            else
                exact_cdf = interp_cdf(x, x2, interp_type);
            end
        
        % 2) use the analytic cdf formulas (in some rare cases)
        elseif dens == 1 % Gaussian
            
            RSRT = R*sigma*R';
            dRSRT = sqrt(2*diag(RSRT));
            exact_cdf = 0.5*(1 + erf((x - R*mu)./dRSRT));
            
        elseif dens == 2 % GM
            
            exact_cdf = 0;
            for i = 1:nr_mixt
                RSRT = R*sigmas(:,dim*i-dim+1:dim*i)*R';
                dRSRT = sqrt(2*diag(RSRT));
                exact_cdf = exact_cdf + ...
                    w(i)*0.5*(1 + erf((x - R*mus(:,i))./dRSRT));
            end  
        else
            disp('Error! Unknown cdf, cannot compute the distance.');
            return;
        end 
    end
    
    %% compute the required 1d statistic to the current direction j
    if stat1d == 1
        [dn(1),dn(2)] = ks_statistic(e_ecdf, exact_cdf, two_sample);
    elseif stat1d == 2
        [dn(1),dn(2)] = cvm_statistic(e_ecdf, exact_cdf, two_sample);
    elseif stat1d == 3
        [dn(1),dn(2)] = anderson_statistic(e_ecdf, exact_cdf, two_sample);
    elseif stat1d == 4
        [dn(1),dn(2)] = cvm_statistic(e_ecdf, exact_cdf, two_sample, 1);
    elseif stat1d >= 5 && stat1d <= 9
        [dn(1),dn(2)] = other_statistic(e_ecdf, exact_cdf, two_sample, stat1d - 4);
    elseif stat1d == 10
        dn(1) = energy_distance(e_ecdf, exact_cdf, 1); % todo: pre-sorting the data is unnecessary!
        dn(2) = dn(1);
    else
        disp('Error! Invalid 1d statistic, cannot compute the distance.');
        return;
    end
    
    %% select the value with desired scaling
    if use_normalised_stat == 1
        d(j) = dn(1);
    else
        d(j) = dn(2);
    end
end

%% compute means and maximums of the discrepancies to different directions 
dn_sum = sum(d)/length(d);
dn_max = max(d);
dn_all = d;
end


% -------------------------------------------------------------------------

function eicdf = interp_cdf(points, samples, interp_type)
% Use interpolation to approximate the empirical cdf at given points. The
% interpolation estimate for empirical cdf is continuous unlike the
% original empirical cdf that has jumps at the data points.

% input:
% points, datapoints where the smoothed ecdf is to be computed
% samples, points that are used to generate ecdf
% interp_type, which interpolation method is used: 1 if linear interp., 
% 2 if spline or cubic interpolation
%
% output:
% eicdf, values of the "smoothed" empirical cdf at given points

tol = 1e-6;

if interp_type == 1
    e = emp_cdf(samples) - 0.5/length(samples); 
    
    % deal with possible duplicates:
    e = unique(e);
    samples1 = unique(samples);
    
    e1 = [0; e(:); 1];
    mins = min(min(samples1),min(points));
    maxs = max(max(samples1),max(points));
    samples1 = [mins - abs(mins)*tol; samples1(:); maxs + abs(maxs)*tol];
    % no input checking, not recommended in matlab's help:
    eicdf = interp1q(samples1, e1, points')'; 

elseif interp_type == 2
    %method = 'spline';
    method = 'pchip';
    
    e = emp_cdf(samples) - 0.5/length(samples);
    
    % deal with possible duplicates:
    e = unique(e);
    samples1 = unique(samples);
    
    e1 = [0, e, 1];
    mins = min(min(samples1),min(points));
    maxs = max(max(samples1),max(points));
    samples1 = [mins - abs(mins)*tol, samples1, maxs + abs(maxs)*tol];
    eicdf = interp1(samples1, e1, points, method); 
    
else
    disp('Incorrect interpolation algorithm!');
    eicdf = NaN(1,length(points));
    return;
    
end

% check that empirical cdf values are indeed valid
eicdf(eicdf > 1) = 1;
eicdf(eicdf < 0) = 0;

end


function eccdf = ksdensity_cdf(points, samples, kernel_type, loc_kernel)
% Computes smoothed empirical cdf using some cdf kernel and bandwith 
% (i.e. stdev) h estimated from the data

% input:
% points, datapoints where the smoothed ecdf is to be computed
% samples, points that are used to generate ecdf
% kernel_type, which cdf kernel is used for smoothing
% loc_kernel, 1 if estimate bandwidth locally (using a very simple and slow 
% method for that currently)
%
% output:
% eccdf, values of the "smoothed" empirical cdf at given points


tol = 1e-36;
ns = 20; % how many points to use for local kernel method

n = length(samples);
m = length(points);

% estimate bandwidth from data: if the underlying model is Gaussian, h is 
% optimal ('Silverman's rule of thumb')
% here we use this choice of h for all problems anyway...

if loc_kernel % this is slow currently! -> TODO: WRITE WITHOUT LOOP!
    
    h = NaN(1,m);
    ns = max(1,ns);
    ns = min(ns,m);
    nsf = ceil(ns/2);
    
    for i = 1:m
        ind = max(1,i-nsf):min(n,i+nsf);
        h(i) = sqrt(var(samples(1,ind))) * (4/(3*n))^(1/5);
        h(i) = min(h(i),quantile(samples(1,ind),0.75) - ...
            quantile(samples(1,ind),0.25));
        h(i) = max(h(i),tol);
    end
    h = repmat(h,m,1);
else
    h = sqrt(var(samples)) * (4/(3*n))^(1/5);
    h = min(h,quantile(samples,0.75) - quantile(samples,0.25));
    h = max(h,tol);
end

eccdf = 1/n * sum(cdf_kernel(repmat(points(:),1,n) - repmat(samples,m,1),...
    h, kernel_type), 2)';
end


function y = cdf_kernel(x,h,kernel_type)
% compute the value of the cdf kernel

% input:
% x, datapoints
% h, bandwidth (standard deviation) for the kernel
% type, which kernel to use
%
% output:
% y, values of the smoothing kernel cdf at given points

if kernel_type == 1 % Gaussian kernel
    y = 0.5 * (1 + erf(x./(sqrt(2)*h)));
    
elseif kernel_type == 2 % Laplace kernel
    s = sqrt(2./h);
    y = (x < 0) .* 0.5 .* exp(s.*x) + (x >= 0) .* (1 - 0.5*exp(-s.*x));
    
elseif kernel_type == 3 % Epanechnikov kernel
    y = (x >= -h & x <= h) .* x./(4*h) .* (3 - (x./h).^2) + (x > h) * 1;
    
end
end







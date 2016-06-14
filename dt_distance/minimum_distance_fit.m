function [] = minimum_distance_fit()
% Estimate parameters of some model by fitting it to data using some
% distance between this data and data generated from the model.

% FIT PARAMETERS TO DATA USING THE PROJECTION DISTANCE!

try
    cd /Users/jarvenpm/work/distribution' tomography'/dt/
end
close all; clear all; format compact;

%%-------------------------------------------------------------------------
%% initial settings related to the profiler and random generator
%%-------------------------------------------------------------------------

use_profiler = 0;
if use_profiler
    profile on
end
tic

% initialize random numbers, keep at global scope (be careful with this!)
global myseed;
global fixed_seed;
fixed_seed = 0;
myseed = clock; myseed = myseed(end); %myseed = 111;
rng(myseed);

%%-------------------------------------------------------------------------
%% initial settings related to the testing and algorithms
%%-------------------------------------------------------------------------

nr_simul = 1;
n = 100; % # "real", observed data samples (sampled only once)
m = 100; % # samples to be generated from the known model
% test_problem: 1 gaussian w/ sigma known, 2 gaussian w/ sigma unknown, 3 gm, 4
% uniform, 5 t w/ sigma known, 6 t w/ sigma unknown, 10 tree model
test_problem = 1; 
% stat1d: 1 == ks, 2 == cvm, 3 == anderson-darling, 4 == Watson, 5 == L^1,
% 6 == L^2, 7 == weighted L^2, 8 == ecdf area^2, 9 == ecdf area, 
% -1 == NN, -2 == SMM, -3 == e-distance, -4 == chi^2 NN
stat1d = 1; 
nr_directions = 100; 
dir_gen_method = 2; % how to generate the directions
% !!!:
smooth_ecdf = 2; % 2 interp smoothing, 1 kernel smoothing, 0 no smoothing
model_sampling = 1; % 0 if use known cdf, 1 if only sampling from cdf's is possible
sum_max = 1; % 1 sum, 2 max

use_mcmc = 0; % EXPERIMENTAL method
use_abc = 0;
use_opt_abc = 0; % EXPERIMENTAL method
algs = [3]; % 1 == genetic, 2 == simulated annealing, 3 == quasi-Newton, 4 == pattern search

% !!!:
add_outliers = 1; % add a few outliers to test robustness
load_real_data = 0; % load real data instead of generating data
plotd = 1;

%% optimization settings
fix_rng_for_opt = 1; % 1 if fix rng seed for optimization
reps_opt = 0; % # restarts from the previous optimum
N_opt = 5; % # initial values for optimization, drawed from the prior

%% opt_abc settings
N_optabc = 25; % how may "samples"
reps_optabc = 0; % # restarts from the previous optimum
alg_optabc = [3]; % which algorithm to use for optimization

%% abc settings
abc_tol = [.75 .5 .25 .15 .12 .09];

%% mcmc settings
trh = 0.13; % "threshold" for defining the posterior, depends on the distance
N = 20000;
burnin = 10000;

%%-------------------------------------------------------------------------
%% initial settings related to the model
%%-------------------------------------------------------------------------

A = []; b = [];
if test_problem == 1 % gaussian (sigma known)
    % set the real parameters
    dim = 1; % dimension of the density
    dens = 1;
    if dim == 1
        mu = 1;
        sigma = 1;
    elseif dim == 2
        mu = [-3; 2];
        sigma = 1/2*[1 1/2; 1/2 2];
    elseif dim == 3
        mu = [-3; 2; 3];
        sigma = [1 .2 .2;.2 2 .2; .2 .2 1.5];
    else
        mu = ones(dim,1);
        sigma = eye(dim,dim);
    end
    params{1} = mu;
    params{2} = sigma;
    
    % dealing with fixed parameters
    theta2p = @(th) theta2params_gaussian1(dim, params, th);
    p2theta = @(params) params2theta_gaussian1(dim, params);
    
    % set prior
    l = -5*ones(dim,1);
    u = 5*ones(dim,1);
    prior_gen = @(s) rand_ab(s,dim,l,u);
    prior_value = @(th) prior_value_uniform(th,l,u);
    
    % generative model
    model = @(s,th) samples_from_gaussian_model(theta2p(th),dim,s);
    
elseif test_problem == 2 % gaussian (sigma unknown)
    dim = 2; % dimension of the density (not the same as # parameters)
    dens = 1;
    if dim == 1
        mu = 1;
        sigma = 1;
    elseif dim == 2
        mu = [-3; 2];
        sigma = 1/2*[1 1/2; 1/2 2];
    elseif dim == 3
        mu = [-3; 2; 3];
        sigma = [1 .2 .2;.2 2 .2; .2 .2 1.5];
    else
        mu = ones(dim,1);
        sigma = eye(dim,dim);
    end
    params{1} = mu;
    params{2} = sigma;
    
    % dealing with fixed parameters
    theta2p = @(th) theta2params_gaussian2(dim, params, th);
    p2theta = @(params) params2theta_gaussian2(dim, params);
    
    % set prior, note: constraints on covariances are not fully taken into
    % account!!!
    l = -5*ones(dim,1);
    u = 5*ones(dim,1);
    prior_gen = @(s) prior_gen_gaussian2(s,dim,l,u);
    prior_value = @(th) prior_value_gaussian2(th,dim,l,u);
    
    % generative model
    model = @(s,th) samples_from_gaussian_model(theta2p(th),dim,s);
    
elseif test_problem == 3 % gaussian mixture (2 components)
    dim = 1; % dimension of the density (not the same as # parameters)
    dens = 2;
    p = 0.7;
    mu1 = zeros(dim,1);
    mu2 = 3.1*ones(dim,1);
    params{1} = [p, 1-p]; % components fixed to 2
    params{2} = [mu1(:), mu2(:)];
    params{3}(1:dim,1:dim) = eye(dim,dim);
    params{3}(1:dim,dim+(1:dim)) = eye(dim,dim);
    
    % dealing with fixed parameters
    theta2p = @(th) theta2params_gm(dim, params, th);
    p2theta = @(params) params2theta_gm(dim, params);
    
    % set prior
    l = -5*ones(2*dim,1);
    u = 5*ones(2*dim,1);
    prior_gen = @(s) rand_ab(s,2*dim,l,u);
    prior_value = @(th) prior_value_uniform(th,l,u);
    
    % generative model
    model = @(s,th) samples_from_gm_model(theta2p(th),dim,s);
    
elseif test_problem == 4 % uniform distribution with both bounds unknown
    dim = 1; % dimension of the density (not the same as # parameters)
    dens = NaN;
    params = repmat([0;1],dim,1);
    
    % dealing with fixed parameters
    theta2p = @(th) th;
    p2theta = @(params) params;
    
    % set prior
    l = -5*ones(2*dim,1);
    u = 5*ones(2*dim,1);
    prior_gen = @(s) prior_gen_2bd_uni(s,dim,l,u);
    prior_value = @(th) prior_value_2bd_uni(th,l,u);
    b = 1e-6*ones(dim,1);
    A = zeros(dim,2*dim);
    for i = 1:dim
        A(i,2*i-1:2*i) = [1 -1];
    end
    
    % generative model
    model = @(s,th) samples_from_2bd_uni(theta2p(th),dim,s);
    
elseif test_problem == 5 % t-distribution (sigma known)
    % set the real parameters
    dim = 1; % dimension of the density
    dens = 4;
    df = 3;
    if dim == 1
        mu = 1;
        sigma = 1;
    elseif dim == 2
        mu = [-3; 2];
        sigma = 1/2*[1 1/2; 1/2 2];
    elseif dim == 3
        mu = [-3; 2; 3];
        sigma = [1 .2 .2;.2 2 .2; .2 .2 1.5];
    else
        mu = ones(dim,1);
        sigma = eye(dim,dim);
    end
    params{1} = mu;
    params{2} = sigma;
    params{3} = df; 
    
    % dealing with fixed parameters
    theta2p = @(th) theta2params_gaussian1(dim, params, th);
    p2theta = @(params) params2theta_gaussian1(dim, params);
    
    % set prior
    l = -5*ones(dim,1);
    u = 5*ones(dim,1);
    prior_gen = @(s) rand_ab(s,dim,l,u);
    prior_value = @(th) prior_value_uniform(th,l,u);
    
    % generative model
    model = @(s,th) samples_from_t_model(theta2p(th),dim,s);
    
elseif test_problem == 6 % t-distribution (sigma unknown)
    dim = 1; % dimension of the density (not the same as # parameters)
    dens = 4;
    df = 3;
    if dim == 1
        mu = 1;
        sigma = 1;
    elseif dim == 2
        mu = [-3; 2];
        sigma = 1/2*[1 1/2; 1/2 2];
    elseif dim == 3
        mu = [-3; 2; 3];
        sigma = [1 .2 .2;.2 2 .2; .2 .2 1.5];
    else
        mu = ones(dim,1);
        sigma = eye(dim,dim);
    end
    params{1} = mu;
    params{2} = sigma;
    params{3} = df;
    
    % dealing with fixed parameters
    theta2p = @(th) theta2params_gaussian2(dim, params, th);
    p2theta = @(params) params2theta_gaussian2(dim, params);
    
    % set prior, note: constraints on covariances are not fully taken into
    % account!!!
    l = -5*ones(dim,1);
    u = 5*ones(dim,1);
    prior_gen = @(s) prior_gen_gaussian2(s,dim,l,u);
    prior_value = @(th) prior_value_gaussian2(th,dim,l,u);
    
    % generative model
    model = @(s,th) samples_from_t_model(theta2p(th),dim,s);
    
    
elseif test_problem == 10 % tree model
    dim = 2; % keep at 2
    dens = NaN;
    np = 2; % how many params to estimate (1,2,3,4)
    params = [40 2 40 0.5]';
    params = params(1:np);
    
    % dealing with fixed parameters
    theta2p = @(th) th;
    p2theta = @(pars) pars;
    
    % set prior
    l = [-90 0 0 0]'+1e-4;  l = l(1:np);
    u = [90 5 90 5]'-1e-4;  u = u(1:np);
    prior_gen = @(s) rand_ab(s,np,l,u);
    prior_value = @(th) prior_value_uniform(th,l,u);
    
    % generative model
    model = @(s,th) samples_from_tree_model(theta2p(th),s);
    
elseif test_problem == 11
    % new tests here....
end

%%-------------------------------------------------------------------------
%% end of settings
%%-------------------------------------------------------------------------

% some checks
if ~model_sampling && test_problem > 3
    disp('Error! Tried to use model cdf when it is not known!');
    return;
end
if (use_abc || use_opt_abc) && m ~= n
    disp('Error! m must equal n => OTHERWISE ISSUES WITH "OVERPRECISE ABC"');
    return;
end
if load_real_data
    disp('Using "real" data...');
    nr_simul = 1;
end

% setting up some values
theta_real = p2theta(params);
dim_th = length(theta_real); % dimension of the parameter space (not dim of the model density)
algss = {'genetic','simulated annealing','quasi-Newton','pattern search'};

% results are saved here
theta_mcmc = NaN(dim_th,nr_simul);
theta_abc = NaN(dim_th,1);
theta_opt_abc = NaN(dim_th,1);
theta_alg = cell(length(algs));
theta_ml = NaN(dim_th,nr_simul);

% print some info about the settings to be used
disp('----------------------------');
test_problem,stat1d,model_sampling
if dim > 1
    nr_directions
end
dim,theta_real,dim_th
disp('----------------------------'); disp(' '); 

%%-------------------------------------------------------------------------
%% generate "real" data, directions and initial points for optimization
%%-------------------------------------------------------------------------

real_data = cell(nr_simul);
directions = cell(nr_simul);
for i = 1:nr_simul
    if load_real_data
        % load real data from file etc.
        % real_data{i} = ...
        real_data{i} = [-2*ones(dim,n/2),2*rand(dim,n/2)];
%         real_data{i} = 4*rand(dim,n); % rand data now here...!
%         real_data{i} = [-1-1e-6*rand(dim,n/2),2*rand(dim,n/2)]; 
%         real_data{i} = [-1*ones(dim,1000),2*rand(dim,n-1000)]; 
        
    else
        % generate data that is used as "real" observations
        real_data{i} = model(n,theta_real);
    end
    
    if add_outliers % contaminate about 10% of data
        am = floor(n/10);
        ma = max(real_data{i}(:));
        real_data{i} = [real_data{i}(:,1:n-am), 10*ma + 1.1*rand(dim,am)];
        %real_data{i} = [real_data{i}(:,1:n-1), 4.5*ones(dim,1)];
    end
    
    % compute transformations:
    pr = (i == 1);
    directions{i} = lin_1d_transforms(dim, nr_directions, real_data{i},...
        dir_gen_method, pr, pr);
end
disp([num2str(nr_simul*n), ' samples generated from the "real" model!']);

if ~isempty(algs)
    inputs = prior_gen(N_opt);
end

%%-------------------------------------------------------------------------
%% actual estimation loop starts here
%%-------------------------------------------------------------------------

for i = 1:nr_simul
    disp(['test run ', num2str(i), '/', num2str(nr_simul)]);
    
    % define functions that depend on the data
    distance = @(s1,s2) dt_gof_statistic(s1, s2, [], [], directions{i}, stat1d, 1, smooth_ecdf); % for abc only
    
    posterior = @(th) posterior_dt(th,theta2p,model,prior_value,dens,...
        directions{i},stat1d,real_data{i},model_sampling,m,sum_max,smooth_ecdf,trh); % for mcmc
    
    f = @(th) optimisation_dt(th,theta2p,model,prior_value,dens,...
        directions{i},stat1d,real_data{i},model_sampling,m,sum_max,smooth_ecdf); % for optimization
    
    f_not_smooth = @(th) optimisation_dt(th,theta2p,model,prior_value,dens,...
        directions{i},stat1d,real_data{i},model_sampling,m,sum_max,0); % for gradient-based optimization
    
    %% SOLUTION ALGORITHMS: ***********************************************
    
    %% 1) --- solve using MCMC (compute at most once) ---
    if use_mcmc && i == 1
        P = 0.1 * eye(dim_th);
        theta_0 = theta_real;
        [mcmc_sa, post_values] = metropolis_adaptive(posterior, theta_0, P, N, burnin, 1, 1);
        theta_mcmc(:,i) = mean(mcmc_sa,2);
        show_mcmc_results(mcmc_sa, post_values, 2);
    end
    
    %% 2) --- solve using ABC approach (compute at most once) ---
    if use_abc && i == 1
        disp('Running ABC...');
        [theta_abc,~,h1] = approximate_bayesian_computation(prior_gen,...
            prior_value, model, distance, real_data{i}, dim_th, abc_tol, m);
        disp('ABC finished!');
    end
    
    %% 3) --- solve using (so-called) Opt-ABC (compute at most once) ---
    if use_opt_abc && i == 1
        disp('Running a variant of ABC...');
        theta_opt_abc = optimized_abc(f, prior_gen, dim_th, reps_optabc,...
            N_optabc, alg_optabc(1)); % add constraints
        disp('Variant of ABC finished!');
    end
    
    %% 4) --- solve using deterministic or stochastic optimization --- 
    
    %% !!! FIXED RNG SEED FOR OPTIMISATION PROBLEMS -> f IS DETERMINISTIC!!!
    if fix_rng_for_opt
        fixed_seed = 1;
        myseed = myseed + 456;
    end
    
    for a = algs
        [theta_opt, f_opt] = optimize_with_prior_inputs(f, prior_gen,...
            dim_th, reps_opt, N_opt, a, inputs, l, u, A, b, (nr_simul==1)); 
        theta_alg{a}(:,i) = theta_opt;
    end
    
    %% END OF SOLUTION ALGORITHMS: ****************************************
    
    %% COMPARISON TO 'EXACT' METHODS: ML, BAYES
    % ml
    theta_ml(:,i) = ml_estimate_for_dt_test(real_data{i}, p2theta, test_problem);

    % bayes
    if use_abc
        bayesian_estimate_for_dt_test(h1, real_data{i}, theta_ml(:,i), theta_real,...
            params, l, u, test_problem);
    elseif use_mcmc
        bayesian_estimate_for_dt_test(50, real_data{i}, theta_ml(:,i), theta_real,...
            params, l, u, test_problem);
    end
end % simul
toc
disp(' ');


%% show objective function values f
if nr_simul == 1
    % objective function values at different points
    disp('***************************************');
    disp('Objective function values at optimum:');
    f_at_real = f(theta_real)
    if use_mcmc
        disp('MCMC posterior mean:');
        f(theta_mcmc)
    end
    for a = algs
        disp([algss{a}, ' algorithm:']);
        f(theta_alg{a})
    end
    try
        disp('ML:');
        f(theta_ml)
    end
    disp('***************************************'); 
end

%% show final results!
disp('#######################################');
disp('Real parameter:');
theta_real'
if use_mcmc
    disp('MCMC:');
    m_theta_mcmc = mean(theta_mcmc,2)'
end
if use_opt_abc
    disp('Opt-ABC:');
    theta_opt_abc'
end
for a = algs
    disp([algss{a}, ' algorithm:']);
    mean(theta_alg{a},2)'
end
if use_abc
    disp('ABC:');
    theta_abc(:)'
end
disp('ML:');
m_theta_ml = mean(theta_ml,2)'
disp('#######################################');

%% compare bias & variance of different estimators
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('--- Bias of estimates:');
for a = algs
    disp([algss{a}, ' algorithm:']);
    [mse,bias] = mean_square_error(theta_real, theta_alg{a})
end
disp('ML:');
[mse,bias] = mean_square_error(theta_real, theta_ml)
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');


%% plot objective function
if plotd
    titl = {'Objective function f','','error %'};
    plot_the_obj_function(f,f_not_smooth,smooth_ecdf,theta_real,[l';u'],...
        theta_real,titl); 
end

if use_profiler
    profile viewer
end
end


%%-------------------------------------------------------------------------
%% FUNCTIONS:

%%-------------------------------------------------------------------------
%% computing the objective function/posterior values:

function d = optimisation_dt(theta,th2p,model,prior_value,dens,...
    directions,stat1d,realdata,model_sampling,n,sum_max,smooth_ecdf)
% optimisation criteria

% ----- settings -.-------
rep = 1; % # datasets to generate
reg_pen = 0; % 1 if add regularisation penalty to the distance
lambda = 1; % regularisation parameter
% ------------------------

params = th2p(theta);

% check parameter bounds
if prior_value(theta) <= 0
    d = Inf;
    return;
end

if model_sampling
    % use several datasets to comparison (not tested yet!!!)
    d1 = NaN(rep,1);
    d2 = NaN(rep,1);
    for i = 1:rep
        % generate the data corresponding the stoch model
        samples2 = model(n,theta);
        % compute the distance
        [d1(i),d2(i)] = dt_gof_statistic(realdata, samples2, [], [],...
            directions, stat1d, 1, smooth_ecdf);
    end
else % model known and cdf can be computed exactly -> no need to sample from the model
    [d1,d2] = dt_gof_statistic(realdata, [], dens, params,...
        directions, stat1d, 0);
end

if sum_max == 1
    d = mean(d1);
else
    d = mean(d2);
end

% test: add a regularisation penalty term to distance
if reg_pen
    dim_th = length(theta);
    L = eye(dim_th,dim_th);
    b = ones(dim_th,1);
    d = d + lambda*norm(L*theta - b);
end
end


function p = posterior_dt(theta,th2p,model,prior_value,dens,directions,stat1d,...
    realdata,model_sampling,n,sum_max,smooth_ecdf,trh)
% compute the dt posterior value -> use for mcmc sampling

res = optimisation_dt(theta,th2p,model,prior_value,dens,directions,stat1d,...
    realdata,model_sampling,n,sum_max,smooth_ecdf);

% return log posterior value
if res > trh || isinf(res) || isnan(res)
    p = -inf;
else
    p = 0;
end

end


%%-------------------------------------------------------------------------
%% sampling from the models:

function [] = set_global_prng_seed()

global fixed_seed;
global myseed;
if fixed_seed == 1
    try
        rng(myseed);
        %disp('Using fixed seed when generating data form the model!');
    catch 
        disp('Error with rng when generating samples from the model!');
    end
end

end

function sa = samples_from_gaussian_model(params,dim,n)
% generate n samples from the gaussian model

set_global_prng_seed();
sa = generate_samples(n, dim, 1, params);
end


function sa = samples_from_gm_model(params,dim,n)
% generate n samples from the gaussian mixture model

set_global_prng_seed();
sa = generate_samples(n, dim, 2, params);
end


function sa = samples_from_t_model(params,dim,n)
% generate n samples from the gaussian model

set_global_prng_seed();
sa = generate_samples(n, dim, 4, params);
end


function sa = samples_from_tree_model(params,n)
% generate n samples from the "tree" model 

set_global_prng_seed();
sa = simple_tree_model(params, min(4,n), 1, 0); % note n is not amount of data here!!!
end


function sa = samples_from_2bd_uni(params,dim,n)
% generate n samples from dim-dimensional uniform with bounds in params

set_global_prng_seed();
sa = rand_ab(n,dim,params(1:2:end),params(2:2:end));
end


%%-------------------------------------------------------------------------
%% priors values and generating values from the prior: 
% (uniform prior for all, exception is the gaussian w/ unknown variance)

function pr = prior_value_uniform(th, l, u)
% general function, if some limits are missing, don't check them; if too
% many limits, ignore the extra ones!

lim1 = 1:min(length(th),length(l));
lim2 = 1:min(length(th),length(u));
pr = ~(sum(l(lim1) > th(lim1)) || sum(u(lim2) < th(lim2))); 
end


function pr = prior_value_gaussian2(th, dim, l, u)
% check the bounds, return 1 if ok, 0 if out of bounds

% first check bounds
pr = prior_value_uniform(th, l, u);
if pr == 0
    return;
end
% then, check if cov matrix is valid
S = theta2params_gaussian2(dim, [], th);
S = S{2};
% test for NaNs and Infs
p1 = sum(isnan(S(:))) + sum(isinf(S(:)));
% test if the matrix is symmetric
p2 = ~isequal(S,S');
% test if the matrix is pos.def.
[~,p3] = chol(S);
if p1 + p2 + p3
    pr = 0;
end
end


function pr = prior_value_2bd_uni(th,l,u)

% first check bounds
pr = prior_value_uniform(th, l, u);
if pr == 0
    return;
end
% check that for any index i it hold that a_i < b_i
if(any(th(1:2:end) >= th(2:2:end)))
    pr = 0;
    return;
end
end


function sa = prior_gen_gaussian2(s, dim, l, u)
%s, amount of samples

sa = NaN(dim + dim*(dim+1)/2,s);
sa(1:dim,:) = rand_ab(s, dim, l, u); % generate mean

% NOTE!!!!
% Covariance must be spd matrix! Prior for cov could be e.g. Wishart 
% distribution but some random spd matrix is generated for now!

for i = 1:s % generate cov matrices
    sigma = randn(dim,dim);
    sigma = sigma'*sigma;
    sa(dim+1:end,i) = collect_upperdiag(sigma);
end
end


function sa = prior_gen_2bd_uni(s,dim,l,u)
% Generate s samples from the prior with bounds l and u and it must hold
% a_i < b_i for any i

sa = NaN(2*dim,s);
for i = 1:dim
    lo = 2*i-1;
    up = 2*i;
    sa(lo:up,:) = rand_uni_simplex(2,s);
    sa(lo,:) = l(i) * sa(lo,:);
    sa(up,:) = u(i) * (1 - sa(up,:)); 
end
end


%%-------------------------------------------------------------------------
%% conversion functions between estimated and the full set of parameters:

function theta = params2theta_gaussian1(dim, params)
% Extract the parameters to be estimated from the full parameter set

theta = params{1}(:);
end


function theta = params2theta_gaussian2(dim, params)
% Extract the parameters to be estimated from the full parameter set

mu = params{1};
sigma = params{2};
s = collect_upperdiag(sigma);
theta = [mu; s];
end


function theta = params2theta_gm(dim, params)
% Extract the parameters to be estimated from the full parameter set

theta = params{2}(:);
end



function params = theta2params_gaussian1(dim, params, theta)
% Add the current parameters to be estimated (theta) to the full parameter
% set (params)

params{1} = theta;
% params{2} contains still the cov matrix
end


function params = theta2params_gaussian2(dim, params, theta)
% Add the current parameters to be estimated (theta) to the full parameter
% set (params)

params{1} = theta(1:dim);
sigma = make_cov_matrix(theta(dim+1:end));
params{2} = sigma;
end


function params = theta2params_gm(dim, params, theta)
% Add the current parameters to be estimated (theta) to the full parameter
% set (params)

params{2} = reshape(theta(:),dim,2);
% params{1} and params{3} contain still weights and cov matrices
end






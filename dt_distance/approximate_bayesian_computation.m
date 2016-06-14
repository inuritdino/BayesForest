function [theta_mean, theta_cov, h1] = approximate_bayesian_computation(prior,...
    prior_value, model, distance, real_data, dim_th, tol, n, method) % ok
% Use ABC for simulating from the posterior
% - MJ

% input:
% prior, handle to function that generates samples from the prior prior =
% prior(N)
% prior_value, handle to function that computes the value of the prior 
% (prior value can be unnormalized) prior_value = prior_value(theta)
% model, handle to function (i.e. generative model) that generates data, 
% model = model(n,theta)
% distance, handle to function that computes the "distance" between 2
% generated datasets, distance = distance(data1,data2)
% realdata, observed experimental data
% dim_th, dimension of the parameter space
% tol, same as epsilon i.e. tolerance level(s) for the distance if they are
% not computed adaptively
% n, #samples from the generative model for each theta value 
% method, which ABC method to use, % 1 rejection algorithm, 2 PMC (with or
% without adaptive tolerance iteration)
%
% output:
% theta_mean, (approximative) posterior mean
% theta_cov, (approximative) posterior covariance matrix
% h1, handle to the plot of the results (use for comparing ABC to the exact solution)
%
% Note:
% THIS FUNCTION IS NOT FINISHED ALTHOUGH IT SEEMS TO WORK!!!
%
% TODO: 
% check and finish dealing with the weights
% check and finish my idea of adaptive tol control

if nargin < 9
    method = 2; 
end
if nargin < 8
    n = size(real_data,2); 
end
if nargin < 7
    tol = 0.1; 
end

%% settings ---------------------------------------------------------------
N = 100; % #samples from prior (at most)
plot_w = 0; % plot weights at each step
quant_tol = 1; % simple,adaptive tol based on quantiles and acceptance prob
p_accmin = 0.001; % stopping criteria
popul_rel_err_tol = 0.0; % stop computing new population if the error is smaller than this value
alpha = 0.9; % cut-off quantile of distance for the next PMC iteration
%% ------------------------------------------------------------------------

theta_mean = []; theta_cov = [];
dim = size(real_data,1);

%% Basic ABC rejection algorithm
if method == 1
    tol = tol(1);
    i = 0; % # iterations
    s = 0; % current sample
    max_tries = 100000;
    thetas = NaN(dim_th,N); % values from prior
    d = NaN(1,N); % distances corresponding thetas
    pr_all = NaN(dim_th,max_tries); % save all prior values
    while s < N && i < max_tries
        i = i + 1;
        % generate samples from prior
        pr = prior(1);
        pr_all(:,i) = pr;
        % generate data from the model
        model_data = model(n, pr);
        % compute the distance
        d1 = distance(model_data, real_data);
        % check acceptance
        if d1 <= tol
            s = s + 1;
            thetas(:,s) = pr;
            d(s) = d1;
        end
    end
    if i >= max_tries
        disp('Computing was stopped, tol probably too small!');
    end
    
    % check that no NaN is returned
    thetas = reshape(thetas(~isnan(thetas)),dim_th,s);
    d = d(~isnan(d));
    pr_all = reshape(pr_all(~isnan(pr_all)),dim_th,i);
    
    
    % illustrate the results of the rejection algorithm!
    tol
    tot_iter_computed = i
    accepted = length(d)
    acc_ratio = accepted/tot_iter_computed
    minimum_distance = min(d)
    [theta_mean, theta_cov] = show_mcmc_results(thetas, -d, 1);
    
    plot_abc_samples(pr_all, 1, 0);
    plot_abc_samples(thetas, 1, 1);
    return;
end
    
%--------------------------------------------------------------------------

%% ABC PMC (population Monte Carlo) algorithm
if method == 2
    %% my adaptive tol
    if quant_tol 
        tol = Inf(max(length(tol),100),1);
        popul_rel_err = Inf;
    end
    
    %% initialize
    max_tries = 1000000; % exit if this is not satisfied
    i = 0; % iteration in population
    j = 0; 
    tot = 0;
    t = 1; % current population
    T = length(tol); % # populations to compute
    % elements: (1) vec element (1...dim), (2) iter in population (1...N), 
    % (3) current population (1...T)
    thetas = NaN(dim_th,N,T); 
    d = NaN(1,N,T); 
    disp([num2str(T), ' populations to simulate!']);
    disp('Computing population 1...');
    disp(['eps=', num2str(tol(1))]);
    
    %% iteration t = 1
    while i < N && j < max_tries
        j = j + 1; tot = tot + 1;
        % generate samples from prior
        pr = prior(1);
        % generate data from the model
        model_data = model(n, pr);
        % compute the distance
        d1 = distance(model_data, real_data);
        % check acceptance
        if d1 < tol(1)
            i = i + 1;
            thetas(:,i,t) = pr;
            d(1,i,t) = d1;
        end
    end
    if j >= max_tries
        disp('Computing was stopped already in first population!');
        disp('Initial tol probably too small!');
        disp('Quitting...');
        return;
    end
    
    w = NaN(N,T);
    w(:,t) = 1/N; % set initial weights
    w(:,t) = check_scaling(w(:,t));
    sigma2 = 2*cov(thetas(:,:,t)'); % compute initial cov, weights are the same here
    
    % plot the prior
    acc_ratio = N/j;
    disp(['Tries: ', num2str(j)]);
    disp(['Acceptance ratio: ', num2str(acc_ratio)]);
    disp(' ');
    
    %plot_abc_samples(thetas(:,:,t), 2, t);
    h1 = plot_abc_samples_kernel(thetas(:,:,t), 2, t);
    
    %% iteration t > 1
    for t = 2:T
        if quant_tol
            tol(t) = quantile(sort(d(1,i,t-1)),1-alpha);
        end
        
        disp(['Computing population ', num2str(t), '...']);
        disp(['eps=', num2str(tol(t))]);
        h = waitbar(0,['Computing population ', num2str(t), '...'],...
            'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
        setappdata(h,'canceling',0);

        i = 0; % iteration in population
        j = 0;
        sigma2
        cancelled = 0;
        chol_sigma2 = chol(sigma2,'lower');
        while i < N && j < max_tries
            j = j + 1; tot = tot + 1;
            % generate samples from previous iteration
            ipr = rnd_categorical(1,w(:,t-1));
            pr = thetas(:,ipr,t-1);
            % perturb current sample using adaptive Gaussian kernel
            pr = pr + chol_sigma2 * randn(dim_th,1);
            % generate data from the model
            model_data = model(n, pr);
            % compute the distance
            d1 = distance(model_data, real_data);
            % check acceptance
            if d1 < tol(t)
                i = i + 1;
                thetas(:,i,t) = pr;
                d(1,i,t) = d1;
                % update weight
                w(i,t) = prior_value(thetas(:,i,t)) / sum(w(:,t-1) .* ...
                    mvnpdf(thetas(:,:,t-1)', thetas(:,i,t)', sigma2));
                if rem(i,floor(N/100)) == 0
                    waitbar(i/N,h);
                end
            end
            if getappdata(h,'canceling')
                cancelled = 1;
                break;
            end
        end
        delete(h); % delete the waitbar
        if j >= max_tries || cancelled
            if j >= max_tries
                disp('Computing was stopped as maximum number iterations were performed!'); 
            else
                disp('Computation was cancelled!');
            end
            t = t-1;
            h1 = plot_abc_samples_kernel(thetas(:,:,t), 2, t, h1, 1, w(:,t));
            % plot also non-weighted version:
            %h1 = plot_abc_samples_kernel(thetas(:,:,t), 2, t, h1, 2);
            plot_weights(w(:,t),t,thetas(:,:,t));
            break;
        end
        w(:,t) = check_scaling(w(:,t)); % check weights
        sigma2 = 2*weighted_covariance(thetas(:,:,t)',w(:,t)); % compute weighted emp cov
        
        % current acceptance ratio
        acc_ratio = N/j;
        
        % plot current population and some info
        disp(['Tries: ', num2str(j)]);
        disp(['Acceptance ratio: ', num2str(acc_ratio)]);
        disp(' ');
        if plot_w
            plot_weights(w(:,t),t,thetas(:,:,t));
        end
        
        % compute the difference between populations in two last population
        if quant_tol && t >= 3
            old_dist = 1;%distance(thetas(:,:,t-1),thetas(:,:,t-2));
            new_dist = 1;%distance(thetas(:,:,t),thetas(:,:,t-1));
            %popul_rel_err = abs(new_dist - old_dist)/old_dist;
            popul_rel_err = new_dist/old_dist;
            disp(['Difference of last two populations: ', num2str(new_dist)]);
            disp(['Relative error in populations: ', num2str(popul_rel_err)]);
        end
        
        %plot_abc_samples(thetas(:,:,t), 2, t);
        h1 = plot_abc_samples_kernel(thetas(:,:,t), 2, t, h1, 0, w(:,t));
        
        if quant_tol && (acc_ratio <= p_accmin || popul_rel_err <= popul_rel_err_tol)
            disp('Computing was stopped, potential final posterior was found!');
            h1 = plot_abc_samples_kernel(thetas(:,:,t), 2, t, h1, 1, w(:,t));
            % plot also non-weighted version:
            %h1 = plot_abc_samples_kernel(thetas(:,:,t), 2, t, h1, 2);
            plot_weights(w(:,t),t,thetas(:,:,t));
            break;
        end
    end
    
    % prepare the results i.e. take the final population computed fully
    thetas = thetas(:,:,t);
    weights = w(:,t);
    d = d(:,:,t);
    
    % show results & plotting
    disp(' ');
    tot_iter_computed = tot
    minimum_distance = min(d)
    disp('Mean (weighted):');
    theta_mean = weighted_mean(thetas', weights)
    disp('Covariance (weighted)');
    theta_cov = weighted_covariance(thetas', weights)
    disp('Correlation matrix:');
    [~,corr_matrix] = cov2corr(theta_cov);
    corr_matrix
    disp('Min,max:');
    theta_min = min(thetas,[],2)'
    theta_max = max(thetas,[],2)'
    disp('Normal approximation for quantiles:');
    theta_mean - 1.96*sqrt(diag(theta_cov))'
    theta_mean + 1.96*sqrt(diag(theta_cov))'
    
    %[theta_mean, theta_cov] = show_mcmc_results(thetas, -d, 2); % this ignores weights
    plot_cov_triangle_abc(thetas,weights);
    
    % plot tolerances
    figure;
    tol = tol(isfinite(tol));
    n = length(tol);
    maxy = 1.1*tol(min(2,n));
    miny = 0;
    bar(1:n,tol);
    ylim([miny,maxy]);
    title('Tolerances');
    xlabel('population');
end

end


%%-------------------------------------------------------------------------

function p = check_scaling(p)
% check that weights are normalised i.e. sum to 1 and if not normalise them
% and print a notification.

pr = 0;
tol = 1e-12;
n = sum(p);
if pr && abs(n-1) > tol
    disp(['Norm constant is ', num2str(n)]);
    disp('Weights are now normalised!');
end
p = p/n;

end


%%-------------------------------------------------------------------------

function [] = plot_abc_samples(s, method, popul)
% plot abc samples as a histogram (ignoring possible weighting)

dim = size(s,1);
n = size(s,2);
methods = {'Rejection ABC','ABC PMC','ABC PMC adapt tol','ABC APMC'};
bins = min(20, ceil(n/20));
figure;
for i = 1:dim
    subplot(1,dim,i);
    hist(s(i,:), bins);
    xlabel(['\theta_', num2str(i)]);
    if nargin < 3 || popul <= 1
        title(['Prior (', methods{method}, ')']);
    else
        title(['Posterior (', methods{method}, '), population ', ...
            num2str(popul)]);
    end
end
end

function [] = plot_cov_triangle_abc(s,weights)

samples = s(:,weights > 1e-6); % ignore samples with very small weight
save_res = 0;
comps_cor = 1:size(s,1);
params = [];
C = [];
plot_cov_triangle(C, samples, params, comps_cor, save_res);

end


function h = plot_abc_samples_kernel(s, method, popul, h, final, weights)
% plot abc samples in the same figure using kernel estimation

if nargin < 6 || isempty(weights)
    weights = [];
end
if nargin < 5 || isempty(final)
    final = 0;
end
if nargin < 4 || isempty(h)
    h = figure;
else
    figure(h);
end

if final == 1 % black color to plot if final population
    c = 'k';
elseif final == 2 % blue plot: debug option...
    c = 'b';
else % otherwise many colors indicating current population
    if popul < 10 
        c = '--m';
    elseif popul < 20
        c = '--g';
    elseif popul < 30
        c = '--c';
    else
        c = '--b';
    end
end

dim = size(s,1);
for i = 1:dim
    subplot(dim,1,i); 
    hold on;
    kern = 'normal'; % 'epanechnikov'
    if isempty(weights)
        [y,x] = ksdensity(s(i,:),'kernel',kern);
    else
        [y,x] = ksdensity(s(i,:),'weights',weights,'kernel',kern);
    end
    plot(x,y,c);
    xlabel(['\theta_', num2str(i)]);
    %suptitle('Populations');
    if popul == 1
        minx = min(s(i,:));
        maxx = max(s(i,:));
        xlim([minx - 0.1*abs(minx), maxx + 0.1*abs(maxx)]);
    end
    hold off;
end
 
end


function [] = plot_weights(w, popul, theta)
% plot weights, w must be a (row or col) vector that contains them

if nargin >= 3
    [~,i] = sort(theta(1,:)); % sort weights according to the 1st param of samples
    w = w(i);
end
n = length(w);

figure;
bar(1:n,w); % todo: put thetas on x-axis!!!
xlim([1,n]);
title(['Weights, population ', num2str(popul)]);

end





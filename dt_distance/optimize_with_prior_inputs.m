function [theta_opt, f_opt] = optimize_with_prior_inputs(f, prior,...
    dim_th, reps, N, alg, inputs, lb,ub,A,b, pr)
% Solve an optimization problem by using multistart by drawing values from
% the "prior"
% - MJ

% input:
% f: handle to the function being minimized [f = f(theta)]
% prior: handle to the function that generates parameter values for optimization 
% i.e. prior returns dim*N matrix of parameter values that are used
% separately as initial points for multistart optimization [prior = prior(N)]
% reps: amount of times the optimization is repeated using the previous result
% as a new initial point (useful for stochastic optimization only, can be 
% set to 0)
% N: how many initial points for optimization
% alg: which algorithm is used for optimisation: 1 genetic algorithm, 
% 2 simulated annealing, 3 Quasi-Nnewton, 4 pattern search
% inputs: if provided, use these values as initial points instead of drawing
% values from the prior
% lb, ub: lower and upper bounds for the parameter so that lb < x < ub
% A, b: inequality constraints for the parameters so that A*x < b
% pr, 1 if print extra info and do some plottings, 0 no printing/plotting
%
% output:
% theta_opt: optimum value found
% f_opt: objective function at theta_opt


% ---- settings ----
show_waitbar = 1;
% ------------------

if nargin < 12
    pr = 1;
end
if nargin < 11
    A = [];  b = [];
end
if nargin < 9
    lb = [];  ub = [];
end

% draw initial points from the "prior" if the initial points are not given 
% as input
if nargin < 7
    inputs = prior(N);
end

algss = {'genetic','simulated annealing','quasi-Newton','pattern search','',''};

% initialize
thetas = NaN(dim_th,N);
fs = NaN(1,N);

%--------------------------------------------------------------------------

if show_waitbar
    h = waitbar(0,['Optimization with ',num2str(N*(1+reps)),' initial states...'],...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    setappdata(h,'canceling',0);
end

% optimize
for i = 1:size(inputs,2) % this loop is embarrassingly parallellizable 
    
    if show_waitbar && getappdata(h,'canceling')
        break;
    end
    
    theta0 = inputs(:,i);
    theta_old = theta0;
    f_at_theta_old = Inf;
    
    for j = 1:reps + 1
        % todo: in genetic alg might be sensible to use all prior values as
        % initial population
        [theta_new, f_at_theta_new] = stochastic_optimisation(f, theta_old,...
            alg, pr, lb,ub,A,b); 
        if pr
            f_at_theta_new
        end
        % check if better point after the repeat
        if f_at_theta_new <= f_at_theta_old
            theta_old = theta_new;
            f_at_theta_old = f_at_theta_new;
        elseif pr
            disp('The old point was better than the new one!');
        end
    end
    thetas(:,i) = theta_old;
    fs(1,i) = f_at_theta_old;
    
    if show_waitbar
        waitbar(i/N,h);
    end
end

if show_waitbar
    delete(h); % delete the waitbar
end

%--------------------------------------------------------------------------

% return only the best value
[f_opt, ind] = min(fs);
theta_opt = thetas(:,ind);

% plot all values of objective function
if pr && length(fs) > 1
    %fs
    figure;
    bar(1:length(fs),fs);
    xlim([1-0.5,length(fs)+0.5]);
    ylim([0,min(1,1.1*max(fs))]);
    xlabel('Prior draw');
    ylabel('Values of f');
    msg = ['Algorithm ', num2str(alg), ' (', algss{alg}, ')'];
    title(msg);
end

% print ending info
disp(['Algorithm #', num2str(alg), ' finished!']);
disp(' ');

end






function [theta, f_at_theta] = stochastic_optimisation(f, theta0, alg, pr, lb,ub,A,b)
% Wrapper for solving the minimum of given objective function f using 
% stochastic (or non-stochastic) optimization algorithms
% - MJ

% input:
% f: handle to the function to be minimized
% theta0: initial point for optimization (as column vector)
% alg: 1 genetic algorithm, 2 simulated annealing, 3 Quasi-Nnewton, 4 pattern search
% pr: 2 print all info, 1 if print extra info, 0 print no info
% lb and ub: bounds so that lb < x < ub
% A and b: constraints so that A*x < b 
% 
% output:
% theta, optimal value found
% f_at_theta, objective function at theta
%
% Note: one needs to transform the problem so that is it minimization 
% problem before calling this function


% --- general settings
max_time_min = 30;
% --------------------

if nargin < 8
    A = [];  b = [];
end
if nargin < 6
    lb = [];  ub = [];
end
if nargin < 4
    pr = 1;
end
if nargin < 3
    alg = 1;
end

dim = length(theta0);
theta = NaN(size(theta0));
f_at_theta = NaN(1,1);

if ~isempty(lb) % complete possible missing limits
    lb = [lb(:);-inf(dim-length(lb),1)];
end
if ~isempty(ub)
    ub = [ub(:);inf(dim-length(ub),1)];
end

if alg == 1 %% genetic algorithm
    % set options etc
    disp('Solving with genetic algorithm...');
    
    options = gaoptimset('Display','final',... % off, iter, diagnose, final 
        'FitnessLimit',0,...
        'InitialPopulation',theta0',...
        'PopulationSize',50,... % default was only 20
        'Timelimit',max_time_min*60,...
        'UseParallel','never',...
        'Vectorized','off');
    g = @(x) f(x');
    nvars = length(theta0);
    
    % compute
    [theta,f_at_theta,exitflag,output] = ga(g, nvars,A,b,[],[],lb,ub,[],options);
    
    % print output
    if pr >= 1
        output
    end
    theta = theta';
    
elseif alg == 2 %% simulated annealing algorithm
    % set options etc
    % note: no possibility to set A*x < b type constraint
    disp('Solving with simulated annealing algorithm...');
    
    options = gaoptimset('Display','final',... % off, iter, diagnose, final
        'Timelimit',max_time_min*60);
    
    % compute
    [theta,f_at_theta,exitflag,output] = simulannealbnd(f, theta0,lb,ub,options);
    
    % print output
    if pr >= 1
        output
    end
    
elseif alg == 3 %% deterministic, gradient-based optimization
    if isempty(ub) && isempty(b)
        % no constraints -> use medium-scale algorithm i.e. BFGS Quasi-Newton
        
        % set options etc
        disp('Solving with medium-scale Quasi-Newton method...');
        options = optimset('Diagnostics','off',...
            'Display','iter',... % off, iter, iter-detailed, notify, notify-detailed, final, final-detailed
            'FinDiffType','forward',... % 'central' more accurate
            'FunValCheck','off',... %on
            'LargeScale','off',...
            'MaxFunEvals', 500*dim,...
            'HessUpdate','bfgs');
        
        % compute
        [theta,f_at_theta,exitflag,output,grad,hessian] = fminunc(f, theta0,options);
        
        % print output
        if pr >= 1
            output
        end
        if pr >= 2
            disp('gradient:');
            grad
            disp('hessian:');
            hessian
        end
        
        % TEST: print initial value if optimization failed for some reason
        if ~isfinite(f_at_theta)
            disp('initial value for optimization was:');
            theta0
        end
    else
        % use active-set, interior-point or sqp method (default active-set)
        % set options etc
        
        algor = 'active-set';
        %algor = 'interior-point';
        %algor = 'sqp';
        disp('Solving with medium-scale Quasi-Newton method with constraints...');
        options = optimset('Algorithm',algor,...
            'Diagnostics','off',...
            'Display','iter',... % off, iter, iter-detailed, notify, notify-detailed, final, final-detailed
            'FinDiffType','forward',... % 'central' more accurate, problems with bounds?
            'FunValCheck','off',... %on
            'MaxFunEvals', 500*dim,...
            'UseParallel','never');
        
        % compute
        [theta,f_at_theta,exitflag,output,lambda,grad,hessian] = ...
            fmincon(f, theta0,A,b,[],[],lb,ub,[],options);
        
        % print output
        if pr >= 1
            output
        end
        if pr >= 2
            disp('gradient:');
            grad
            disp('hessian:');
            hessian
        end
        
        % TEST: print initial value if optimization failed for some reason
        if ~isfinite(f_at_theta)
            disp('initial value for optimization was:');
            theta0
        end
    end 
    
elseif alg == 4 %% pattern search algorithm
    % set options etc
    disp('Solving with pattern search...');
    
    options = psoptimset('Display','final',... % off, iter, diagnose, final
        'Timelimit',max_time_min*60,...
        'UseParallel','never',...
        'Vectorized','off');
    
    % compute
    [theta,f_at_theta,exitflag,output] = patternsearch(f, theta0,...
        A,b,[],[],lb,ub,[],options);
    
    % print output
    if pr >= 1
        output
    end
    
else
    disp('Error, incorrect algorithm id was given!');
end

end





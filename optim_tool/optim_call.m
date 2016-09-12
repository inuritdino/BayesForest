function [xVal, FVal, Output, problem, X0] = optim_call(data,model,varargin)
% General interface for calling an optimizational procedure: fitting the
% "data" with model.
% USAGE:
%       [xVal, FVal, Output, problem, X0] = optim_call(data,model,varargin)
%
% XVAL, FVAL, OUTPUT are the general output of the MATLAB optimization
% procedures, namely, the best point, the function value at the best point,
% and the output structure.
% X0 is a vector of the initial points, if the algorithm is gradient-based
% one.
% PROBLEM is the problem structure, containing, besides all, the options
% structure.
%
% DATA is the data distribution function (DF). size(DATA) = [M N], where M
% is the number of variates and N is the number of samples. Multiple DF's
% can be set by putting single DF's into a cell-array.
%
% MODEL is a function that generates the model distribution function to fit
% to the DATA. Characteristics and dimensions of MODEL are the same as
% those for DATA.
%
% VARARGIN inputs:
%
% Solvers: the default solver is GA, FMINCON - the Newton-based algorithm
% with constraints, FMINSEARCH - the simplex algorithm. The gradient-based
% methods might use multiple initial points evenly distributed within the
% user-specified boundaries, the number of such points determined by the
% 'intcon' value in the problem struct (since there is no integer constraints
% for the gradient-based methods), see the general MATLAB documentation for
% setting up the problem struct and description of the problem struct
% fields specification below.
%
% OPTS defines the options structure for the problem (see general MATLAB
% documentation for optimization). Syntax ...,'opts',{<your_opts_specifications>},...
%
% FITNESS defines an external fitness function (through function handle).
% The default fitness function is the dt_distance(...) with 100 directions
% and Kolmogorov-Smirnov statistic to use. The default fitness uses the
% interpolation-based smoothing for the gradient-based methods, for the GA
% no smoothing is applied.
%
% PROBLEM struct fields can be directly set up from VARARGIN. Examples:
% 'ub','lb' for the upper and lower bounds, respectively, 'intcon' for the
% integer constraints,'x0' for the initial point etc. See the main MATLAB
% documentation.
%
% OPTIM_CALL(...,'w',WGHT,...) specifies the weight coefficients WGHT for the
% weighted sum of the distance of multiple DF's. SUM(WGHT) should be 1.0,
% but no checking on this is performed in the function.
%
% OPTIM_CALL(...,'dirs',ND,...) specifies the number of directions for the
% dt_distance() function. Default is 100.
%
% OPTIM_CALL(...,'stat',S,...) specifies the statistic type. S is a numeric
% value, for possible choices see dt_distance() description. Default is 1
% (Kolmogorov-Smirnov).
%% Preliminaries
if(isempty(varargin))% Workaround to have varargin non-empty for strcmp
    varargin = {''};
end
% Set up the problem struct
problem = struct;
% Weight values for multiple DF's
wght = [];
tf = strcmp('w',varargin);
if(find(tf))
    wght = varargin{find(tf)+1};
end
% Number of directions for the dt_distance()
dirs = 100;
tf = strcmp('dirs',varargin);
if(find(tf))
    dirs = varargin{find(tf)+1};
end
% Statistic to compute for 1D comparisons within dt_distance()
stat = 1;
tf = strcmp('stat',varargin);
if(find(tf))
    stat = varargin{find(tf)+1};
end
scale = false;
tf = strcmp('scale',varargin);
if(find(tf))
    scale = varargin{find(tf)+1};
end
%% Set up the solvers.
% GA id=1, default solver
solver = 1;
% FMINCON solver(Newton's method). id=2.
tf = strcmp('fmincon',varargin);
if(find(tf))
    solver = 2;
end
% FMINSEARCH (SIMPLEX algorithm)
tf = strcmp('fminsearch',varargin);
if(find(tf))
    solver = 3;
end
% GA for multi-objective function
tf = strcmp('gamultiobj',varargin);
if(find(tf))
    solver = 4;
end

%% Set up the options
tf = strcmp('opts',varargin);
if(find(tf))
    opts = set_opts(solver,varargin{find(tf)+1}{:});
else
    opts = set_opts(solver,[]);
end

%% Set up the fitness function
tf = strcmp('fitness',varargin);
fun = [];
if(find(tf))
    fun = varargin{find(tf)+1};
end

%% Solver specific options and problem structures
if(solver == 1)
    % Set some options to the GA
    %opts = set_opts(solver,varargin{:});
    problem.solver = 'ga';
    problem.fitnessfcn = @fitness;
elseif(solver == 2)
    %opts = set_opts(solver,varargin{:});
    problem.solver = 'fmincon';
    problem.objective = @fitness;
elseif(solver == 3)% simplex
    % FMINSEARCH does not use problem, we use it to specify objfun only.
    problem.objective = @fitness;
    %opts = set_opts(solver,varargin{:});
elseif(solver == 4)% gamultiobj
    problem.solver = 'gamultiobj';
    problem.fitnessfcn = @fitness;
end

%% General options and problem fields
problem = set_problem(problem,varargin{:});
% set options from the algorithm specific opts structure (see prev. cell)
problem.options = opts;
% disp('Problem struct:');
% disp('(hit Enter after you are done)');
% disp(problem);
% pause;
%% CALL THE OPTIMIZATION ROUTINE
[xVal,FVal,exitFlag,Output,X0] = optim_routine(solver,problem,opts);

%% REPORT RESULTS
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
if(solver == 1 || solver == 4)
    fprintf('%% The problem type was: %s\n', Output.problemtype);
    fprintf('%% The best point found: %g\n',xVal);
    fprintf('%% The fitness at the best point: %g\n',FVal);
    fprintf('%% I have stopped because: %d =>\n',exitFlag);
    fprintf('%% %s\n',Output.message);
end
if(solver == 2)
    fprintf('%% Multi-start results:\n');
    fprintf('%% Algorithm: %s\n',Output(1).algorithm);
    fprintf('%% BEST SOLUTION(s)\n');
    bii = find( abs(FVal - min(FVal)) < 1e-04 );
    for ii = bii
        fprintf('#%d) Initial conditions (X0): ',ii);
        for jj = 1:problem.nvars
            fprintf('%g ',X0(jj,ii));
        end
        fprintf('\n');
        fprintf('Best values (X): ');
        for jj = 1:problem.nvars
            fprintf('%g ',xVal(jj,ii));
        end
        fprintf('\n');
        fprintf('Distance: ');
        fprintf('%g\n',FVal(ii));
        fprintf('Function count: ');
        fprintf('%d\n',Output(ii).funcCount);
    end
    [~,I] = min(FVal);
    FVal = FVal(I);
    xVal = xVal(:,I);
    Output = Output(I);
end
% Performance of the algorithm
if(solver == 1)
    fprintf('%% The number of generations was : %d\n', Output.generations);
end
%fprintf('%% The number of function evaluations was : %d\n', Output.funccount);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

%% FITNESS == OBJECTIVE FUNCTION
% NESTED FUNCTION for calculating the fitness
    function z = fitness(X)
        % Fitness is calculated for being optimized value(s) in X.
        if(isempty(fun))
            z = calc_fitness(X,data,model,solver,wght,dirs,stat,scale);
        else
            z = fun(X,data,model);
        end
    end
% END OF NESTED FUNCTION for the fitness calculation
end

%% Fitness calculation
function [z,n] = calc_fitness(X,data,model,solver,w,ndirs,stat1d,scaling)
% Calculate fitness given the data and model function.

% X is an input to the model() function.
model_data = model(X);% generate model data given the function

% Calculate the average distance:
if(solver == 2)%Gradient-based methods
    smooth = 2;
else
    smooth = 0;
end
if (solver == 4)
    [z,~,n] = optim_avg_distance(model_data,data,'w',w,'dirs',ndirs,'stat',stat1d,...
        'scaling',scaling,'smooth',smooth,'multiobj');
else
    [z,~,n] = optim_avg_distance(model_data,data,'w',w,'dirs',ndirs,'stat',stat1d,...
        'scaling',scaling,'smooth',smooth);
end
% if(~isempty(n))
%     fun = @(x,y)sqrt(x.*y./(x+y));
%     disp(['Parameters: ' num2str(X) ':']);
%     disp(['Prefactors = ' num2str(bsxfun(fun,n(1,:),n(2,:)))]);
% end
end

%% Optimization routine call
function [xVal,Fval,exitFlag,Output,X0] = optim_routine(solver,problem,opts)
% Call the specific solver for the optimization problem.
if(solver == 1)% GA
    [xVal,Fval,exitFlag,Output] = ga(problem);
    X0 = [];
elseif(solver == 2)
    if(isempty(problem.intcon) || problem.intcon == 1)
        X0 = problem.x0;
        [xVal,Fval,exitFlag,Output] = fmincon(problem);
    else
        % NOTE: problem.intcon is used as number of initial points
        xVal = zeros(problem.nvars,problem.intcon);
        Fval = zeros(1,problem.intcon);
        exitFlag = zeros(problem.nvars,problem.intcon);
        X0 = zeros(problem.nvars,problem.intcon);
        for ii=1:problem.intcon
            %Random init condition generation
            X0(:,ii) = problem.lb' + (problem.ub-problem.lb)'.*rand(problem.nvars,1);
            problem.x0 = X0(:,ii);
            [xVal(:,ii),Fval(ii),exitFlag(:,ii),out] = fmincon(problem);
            if(ii == 1)
                Output = out;
            else
                Output(ii) = out;
            end
        end
    end
elseif(solver == 3)
    X0 = problem.x0;
    [xVal,Fval,exitFlag,Output] = fminsearch(problem.objective,problem.x0,opts);
elseif(solver == 4)
    [xVal,Fval,exitFlag,Output] = gamultiobj(problem);
    X0 = [];
end
end
%% Problem structure forming
function problem = set_problem(problem,varargin)
tf = strcmp('nvars',varargin);
if(find(tf))
    problem.nvars = varargin{find(tf)+1};
else
    problem.nvars = 2;
end
tf = strcmp('x0',varargin);
if(find(tf))
    problem.x0 = varargin{find(tf)+1};
else
    problem.x0 = ones(1,problem.nvars);
end
tf = strcmp('lb',varargin);
if(find(tf))
    problem.lb = varargin{find(tf)+1};
end
tf = strcmp('ub',varargin);
if(find(tf))
    problem.ub = varargin{find(tf)+1};
end
tf = strcmp('nonlcon',varargin);
if(find(tf))
    problem.nonlcon = varargin{find(tf)+1};
end
tf = strcmp('intcon',varargin);
if(find(tf))
    problem.intcon = varargin{find(tf)+1};
end
end
%% Options setting for the solvers
function opts = set_opts(solver,varargin)
% We can use the same interface to set options structure for all solvers.
% Some solver specific conflicts may be resolved using the solver
% information as the input (see for example 'PlotFcns' option).

if(isempty(varargin))
    if(solver == 1 || solver == 4)
        opts = gaoptimset(@ga);
    elseif(solver == 2)
        opts = optimset('fminunc');
    elseif(solver == 3)
        opts = optimset('fminsearch');
    elseif(solver == 4)
        opts = gaoptimset(@gamultiobj);
    end
else
    opts = struct;
    if(solver == 1 || solver == 4)
        opts = gaoptimset(opts,varargin{:});
    else
        opts = optimset(opts,varargin{:});
    end
end
end
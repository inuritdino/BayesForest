function optim_multi_start(data,model)

data = 10 + 10*randn(1,1000);
data = [data 60 + 10*randn(1,1000)];

xstart = [0 1 6 1];
opts = optimset();
problem = createOptimProblem('fminunc','x0',xstart,'objective',@objective,...
    'options',opts,'lb',[0 0 0 0],'ub',[100 100 100 100]);
%[x fval eflag out] = fminunc(problem)


%gs = GlobalSearch;
ms = MultiStart;


[xmin,fmin,flag,out,allmins] = run(ms,problem,100)

    function d = objective(X)
        d = calc_dist(X,data,model);
    end

end

function d = calc_dist(X,data,model)

% Generate some data
x = X(1) + X(2)*randn(1,1000);
y = X(3) + X(4)*randn(1,1000);
z = [x y];

d = dt_distance(z,data);

end
function conf = bf_process_input(input_file)
% Process the input file.
% Usage:
%	conf = bf_process_input(input_file)
% conf: the output configuration struct processed further in BayesForest.m
% input_file: the input file with general syntaxing scheme 
%	OPTION = VALUE1 [, VALUE2, VALUE3, ... ] (comma can be substituted with
% a space.
% Available OPTIONS:
% target_dir = the target directory where to perform all calculations,
%		the algorithm returns back to the current directory after all
%		calculations are done. Default is './', [string in Unix format].
% scatter = 'segment' or 'branch', or both; indicates which feature table is
%		to be used. Default is 'segment'.
% order = list of topological orders to be used with feature types indicated
%		in 'scatter' option. Default is 1. NOTE: order count is from 0 (trunk).
% qsm_mat_file = the .mat file for the QSM features. The format understandable by
%		gen_scatter2() function (see qsm_data/ folder). Default: empty.
% qsm_cyl_table = a variable name containing the matrix with segment-related 
%		features. Dimension: number of segments x 11. Note the variable
% 		must exist in the workspace. The value may be the name of a file 
% 		containing data for such matrix (`importdata` is used). 
%		Default: []. For the format of the matrix see `help import_qsm_data`
%		at the Matlab prompt.
% qsm_br_table = a variable name containing the matrix with branch-related
%		features. Dimension: number of branches x 5. Note the variable
%		must exist in the workspace. The value may be the name of a file 
%		containing data for such matrix (`importdata` is used). 
%		Default: []. For the format of the matrix see `help import_qsm_data`
%		at the Matlab prompt. 
% qsm_merge = whether to use merging of scatters for the QSM data. Note: SSM
%		scatters must be also produced with the same merging. Default: false.
% segment = list of topological orders to be used with segment-related tables.
%		The option is useful when different orders for segment- and 
% 		branch-related features are to be used. Do not use along with 'scatter' 
% 		and 'order' options. This option has a preference over 'scatter'.
% branch = list of topological orders to be used with branch-related tables.
%		The option is useful when different orders for segment- and 
% 		branch-related features are to be used. Do not use along with 'scatter'
%		and 'order' options. This option has a preference over 'scatter.
% ssm_fun = a function handle to run SSM and produce the scatters. For example:
% 		@(x)ssm_lpfg2(x,...) See LPFG callers from ssm_callers/.
% ssm_fun_best = same as 'ssm_fun', but called upon the final best-fit SSM 
%		simulation. This could be useful when generating output MTG file
%		(time consuming) for visualizations outside of Matlab.
% ga_init_lb = lower boundaries of the parameters for the optimization search at
%		the first run. >= the global lower boundary and usually contains
%		the best solution. Corresponds to the lower value of 'PopInitRange'
%		in Genetic Algorithm(GA) of Matlab Optimization Toolbox.
% ga_init_ub = 'PopInitRange' upper boundary; <= the global upper boundary.
% ga_lb = global lower boundary, corresponds to 'lb' in GA.
% ga_ub = global upper boudnary.
% ga_int_con = indices of the integer variables of the SSM.
% ga_pop_size = Population size for GA.
% ga_gens = number of generations in GA.
% ga_stall = stall limit in GA.
% ga_elite = GA Elite size.
% dt_stat1d = number indicating which 1D statistic to use in the distance code.
%		See dt_distance() code (stat1d input).Default: 1 (KS-test).
% dt_dirs = number of directions to used when calculating the structural distance.
%		Default: 100.
% dt_scale = whether to use scaling in distance calculation. Default: false.
% dt_w = vector of weights if weighted average over the several distance values
%		is to be calculated. Default: [].
% ssm_tree_fun = a function that generates a tree object after best-fit SSM
%		is run (see ssm_fun_best option).
% ga_out_dat = output file name for the GA progress report. Default: 'gaOut.dat'.
% ga_use_par = whether to use parallel optimization. Default: false. Any non-zero
%		value is treated as true.
% ga_rng = random number generator seed value. To reproduce the optimization run.
%		Default: [].
% movie = whether to create a move of evolving structure over optimization runs.
%		Default: true.
% plot = whether to plot the tree structures and, perhaps, to make movie (controlled
%		by 'movie' option). Default: true. Note that scatter plots are created
%		anyways.
% ga_out_fun = a function to run after each GA iteration. Could be a printing
% 		of a progress or similar.
% ga_tol_fun = distance tolerance threshold determines which changes in distance
% 		are considered significant.
% ga_multi = whether to use Multi-Objective GA (Pareto front). Default: false.
%		Note: the Multi-Objective GA is not fully implemented yet.
% ================================================================================
% Note that the new options can be easily implemented here, but their processing
% needs to be coded into BayesForest.m
%
% See also: ga

conf = struct('target_dir','.','scatter','segment','order',1,'qsm_mat_file',[],...
    'qsm_cyl_table',[],'qsm_br_table',[],'segment',[],'branch',[],...
    'ssm_fun',[],'ga_init_lb',[],'ga_init_ub',[],'ga_lb',[],'ga_ub',[],'ga_int_con',[],...
    'ga_pop_size',[],'ga_gens',[],'ga_stall',[],'ga_elite',[],...
    'dt_stat1d',1,'dt_dirs',100,'dt_scale',false,'dt_w',[],...
    'ssm_tree_fun','@()read_mtg(''out.mtg'')','ga_out_dat','gaOut.dat',...
    'ga_use_par',0,'ga_rng',[],'ssm_fun_best',[],'movie',true,'qsm_merge',false,'ga_out_fun',[],...
    'ga_tol_fun',1e-6,'ga_multi',0,'plot',true);

fprintf('Start processing ''%s''...',input_file);

fid = fopen(input_file,'r');

line = fgetl(fid);

% Get continued lines
line = get_continued(fid,line);

while(ischar(line))
    if(~isempty(line) && line(1) ~= '#')
        C = textscan(line,'%s','delimiter','=');% process line
        left = C{1}{1}; left = strtrim(left);% left with no spaces
        right = C{1}{2}; right = strtrim(right);% right, no spaces
        %fprintf('Read: `%s'' = `%s''\n',left,right);
        conf = write_conf(left,right,conf);% save to the conf
    end
    line = fgetl(fid);
    line = get_continued(fid,line);
end

fclose(fid);
fprintf('Done.\n');
end

function line = get_continued(fid,line)
% Check and get a continued line.
k = strfind(line,'...');
while( ischar(line) && ~isempty(k) )
    line = line(1:k(1)-1);
    line = cat(2,line,strtrim(fgetl(fid)));
    k = strfind(line,'...');
end
end

function conf = write_conf(left,right,conf)
% Write the conf structure
if(strcmpi(left,'target_dir'))
    conf.target_dir = right;
elseif(strcmpi(left,'scatter'))
    C = textscan(right,'%s','delimiter',',');
    conf.scatter = C{1};
elseif(strcmpi(left,'order'))
    C = textscan(right,'%s','delimiter',',');
    conf.order = [];% to erase the default value
    for ii = 1:length(C{1})
        conf.order = cat(2,conf.order,str2double(strtrim(C{1}{ii})));
    end
elseif(strcmpi(left,'qsm_mat_file'))
    conf.qsm_mat_file = right;
elseif(strcmpi(left,'qsm_cyl_table'))
    conf.qsm_cyl_table = right;
elseif(strcmpi(left,'qsm_br_table'))
    conf.qsm_br_table = right;
elseif(strcmpi(left,'qsm_merge'))
    conf.qsm_merge = str2double(right) > 0;
elseif(strcmpi(left,'segment'))
    C = textscan(right,'%s','delimiter',',');
    for ii = 1:length(C{1})
        conf.segment = cat(2,conf.segment,str2double(strtrim(C{1}{ii})));
    end
elseif(strcmpi(left,'branch'))
    C = textscan(right,'%s','delimiter',',');
    for ii = 1:length(C{1})
        conf.branch = cat(2,conf.branch,str2double(strtrim(C{1}{ii})));
    end
elseif(strcmpi(left,'ssm_fun'))
    conf.ssm_fun = right;
elseif(strcmpi(left,'ssm_fun_best'))
    conf.ssm_fun_best = right;
elseif(strcmpi(left,'ga_init_lb'))
    C = textscan(right,'%s','delimiter',', ','MultipleDelimsAsOne',true);
    for ii = 1:length(C{1})
        conf.ga_init_lb = cat(2,conf.ga_init_lb,str2double(strtrim(C{1}{ii})));
    end
elseif(strcmpi(left,'ga_init_ub'))
    C = textscan(right,'%s','delimiter',', ','MultipleDelimsAsOne',true);
    for ii = 1:length(C{1})
        conf.ga_init_ub = cat(2,conf.ga_init_ub,str2double(strtrim(C{1}{ii})));
    end
elseif(strcmpi(left,'ga_out_dat'))
    conf.ga_out_dat = right;
elseif(strcmpi(left,'ga_lb'))
    C = textscan(right,'%s','delimiter',', ','MultipleDelimsAsOne',true);
    for ii = 1:length(C{1})
        conf.ga_lb = cat(2,conf.ga_lb,str2double(strtrim(C{1}{ii})));
    end
elseif(strcmpi(left,'ga_ub'))
    C = textscan(right,'%s','delimiter',', ','MultipleDelimsAsOne',true);
    for ii = 1:length(C{1})
        conf.ga_ub = cat(2,conf.ga_ub,str2double(strtrim(C{1}{ii})));
    end
elseif(strcmpi(left,'ga_int_con'))
    C = textscan(right,'%s','delimiter',',');
    for ii = 1:length(C{1})
        conf.ga_int_con = cat(2,conf.ga_int_con,str2double(strtrim(C{1}{ii})));
    end
elseif(strcmpi(left,'ga_out_fun'))
    C = textscan(right,'%s','delimiter',',');
    for ii = 1:length(C{1})
        conf.ga_out_fun = [conf.ga_out_fun {str2func(strtrim(C{1}{ii}))}];
    end
elseif(strcmpi(left,'ga_pop_size'))
    conf.ga_pop_size = str2double(right);
elseif(strcmpi(left,'ga_gens'))
    conf.ga_gens = str2double(right);
elseif(strcmpi(left,'ga_stall'))
    conf.ga_stall = str2double(right);
elseif(strcmpi(left,'ga_elite'))
    conf.ga_elite = str2double(right);
elseif(strcmpi(left,'ga_tol_fun'))
    conf.ga_tol_fun = str2double(right);
elseif(strcmpi(left,'ga_use_par'))
    conf.ga_use_par = str2double(right) > 0;
elseif(strcmpi(left,'ga_rng'))
    conf.ga_rng = uint64(str2double(right));
elseif(strcmpi(left,'ga_multi'))
    conf.ga_multi = str2double(right) > 0;
elseif(strcmpi(left,'dt_stat1d'))
    conf.dt_stat1d = str2double(right);
elseif(strcmpi(left,'dt_dirs'))
    conf.dt_dirs = str2double(right);
elseif(strcmpi(left,'dt_scale'))
    conf.dt_scale = str2double(right) > 0;
elseif(strcmpi(left,'dt_w'))
    C = textscan(right,'%s','delimiter',', ','MultipleDelimsAsOne',true);
    for ii = 1:length(C{1})
        conf.dt_w = cat(2,conf.dt_w,str2double(strtrim(C{1}{ii})));
    end
elseif(strcmpi(left,'movie'))
    conf.movie = str2double(right) > 0;
elseif(strcmpi(left,'plot'))
    conf.plot = str2double(right) > 0;
else
    fprintf('Unknown line: %s = %s\n',left,right);
end

end

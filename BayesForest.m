function BayesForest(input_file)
% Optimize Stochastic Structure Model (SSM) to the Quantitative Structure
% Model (QSM). The former is a simulated model, the latter is (usually) the
% model of a real tree (data).
% USAGE:
%       BayesForest(input_file)
%
% The input_file describes the configuration of the morphological data,
% optimization configuration, SSM and additional details of the procedure.
% Refer to the help pages of `bf_process_input`.
%
% BayesForest creates several figures during optimization (e.g. optimization
% progress) and some figures after the program has finished. The output
% folder (coded by date and time) is created with all figures and data files
% at the end.
%
% SEE ALSO: bf_process_input, import_qsm_data, gen_scatter2

%% Read the configuration file
config = bf_process_input(input_file);
% Copy configuration file
rnd_arr = char(['A':'Z' 'a':'z' '0':'9']);
tmp_input_file = rnd_arr(randi(length(rnd_arr),1,4));
tmp_input_file = ['do_not_change_input_' tmp_input_file '.temporary'];
% copyfile(input_file,tmp_input_file);% copy to temporary file
% movefile(tmp_input_file,config.target_dir);
[~,basename,ext] = fileparts(input_file);
input_file_name = [basename ext];

%% Some preparations
% Enter target directory
currDir = pwd;% remember current directory to return at the end
disp(['Entering target directory:' pwd '/' config.target_dir])
cd(config.target_dir);% move to the target directory
copyfile([currDir '/' input_file],tmp_input_file);
%% Define the input: scatter types, order etc.
if(isempty(config.scatter))
    scat = {'branch','segment'};
else
    scat = config.scatter;
end
if(isempty(config.order))
    order = 1;
else
    order = config.order;
end
%% Load the data
if(isempty(config.qsm_table))
    if(~isempty(config.qsm_mat_file))
        [qsm_bra,qsm_seg,qsm_tree] = gen_scatter2(config.qsm_mat_file);
    elseif(~isempty(config.qsm_cyl_table) && ~isempty(config.qsm_br_table))
        if(exist(config.qsm_cyl_table,'file') == 2 && exist(config.qsm_br_table,'file') == 2)
            fprintf('Reading QSM data from the data files...\n');
            import_qsm_data(config.qsm_br_table,config.qsm_cyl_table,'temp.mat');
        else
            fprintf('Getting QSM data from the workspace...\n');
            import_qsm_data(evalin('base',config.qsm_br_table),evalin('base',config.qsm_cyl_table),'temp.mat');
        end
        [qsm_bra,qsm_seg,qsm_tree] = gen_scatter2('temp.mat');
        delete('temp.mat');
    else
        delete(tmp_input_file);
        error('Error: data (QSM) file is not specified.');
    end
    if(isempty(config.segment) && isempty(config.branch))
        if(config.qsm_merge)
            qsm_scatter = arrange_scatter(qsm_bra,qsm_seg,'scat',scat,...
                'order',order,'merge');
        else
            qsm_scatter = arrange_scatter(qsm_bra,qsm_seg,'scat',scat,...
                'order',order);
        end
    else
        if(config.qsm_merge)
            qsm_scatter = arrange_scatter(qsm_bra,qsm_seg,'branch',config.branch,...
                'segment',config.segment,'merge');
        else
            qsm_scatter = arrange_scatter(qsm_bra,qsm_seg,'branch',config.branch,...
                'segment',config.segment);
        end
    end
end
%% Arbitrary data sets
if(~isempty(config.qsm_table))% arbitrary format scatter
    if(exist(config.qsm_table,'file') == 2)
        fprintf('Reading user data from file...');
        qsm_scatter = importdata(config.qsm_table);
        fprintf('Done.\n');
    else
        fprintf('Reading user data from workspace...');
        qsm_scatter = evalin('base',config.qsm_table);
        fprintf('Done.\n');
    end
end
%% Define the trial model
ssm_fun = str2func(config.ssm_fun);
if(~isempty(config.ssm_fun_best))
	ssm_fun_best = str2func(config.ssm_fun_best);
else
	ssm_fun_best = ssm_fun;
end
%% Technical parameters for the GA
InitRange = [config.ga_init_lb; config.ga_init_ub];
LB = config.ga_lb;
UB = config.ga_ub;

DIM = length(LB);
INTCON = config.ga_int_con;
POPSIZE = config.ga_pop_size;
ELITE = config.ga_elite;
GENS = config.ga_gens;
STALL = config.ga_stall;
TOLFUN = config.ga_tol_fun;
USEPAR = config.ga_use_par;
VECT = 'off';
OUTFUN = config.ga_out_fun;

if(~isempty(config.ga_rng))
    disp(['Fixing RNG seed to ' num2str(config.ga_rng)]);
    rng(config.ga_rng);
end
%% Technical parameters for the Distance
STAT1D = config.dt_stat1d;
DIRS = config.dt_dirs;
SCALE = config.dt_scale;
WGHT = config.dt_w;
%% Optimize
if( config.ga_multi )
    [X,~,Output,Problem,~] = optim_call(qsm_scatter,ssm_fun,'nvars',DIM,'gamultiobj',...
        'opts',{'PlotFcns',{@gaplotpareto,@gaplotparetodistance},'PopulationSize',POPSIZE,...
        'EliteCount',ELITE,'Generations',GENS,'StallGenLimit',STALL,...
        'PopInitRange',InitRange,'UseParallel',USEPAR,'Vectorized',VECT,...
        'TolFun',TOLFUN},...
        'lb',LB,'ub',UB,'intcon',INTCON,'stat',STAT1D,'dirs',DIRS,'scale',SCALE,...
        'w',WGHT);
else
    [X,~,Output,Problem,~] = optim_call(qsm_scatter,ssm_fun,'nvars',DIM,...
        'opts',{'PlotFcns',{@gaplotbestf,@gaplotdistance},'PopulationSize',POPSIZE,...
        'EliteCount',ELITE,'Generations',GENS,'StallGenLimit',STALL,...
        'PopInitRange',InitRange,'UseParallel',USEPAR,'Vectorized',VECT,...
        'OutputFcns',OUTFUN,'TolFun',TOLFUN},...
        'lb',LB,'ub',UB,'intcon',INTCON,'stat',STAT1D,'dirs',DIRS,'scale',SCALE,...
        'w',WGHT);
end
%% Run the best solution
% if(isempty(config.ssm_fun_best))
%     ssm_best_scatter = ssm_fun(X);
% else
%     ssm_best_scatter = ssm_fun_best(X);
% end
%% Create folder for the results
t = clock;
out_name = [num2str(t(3)) '.' num2str(t(2)) '.' num2str(t(1)) '_' ...
    num2str(t(4)) '.' num2str(t(5))];
mkdir(out_name);

%% Plot the trees
if(config.plot)
    %%% Make the final movie
    if(config.movie)
        mov = optim_plot_generations(config.ga_out_dat,ssm_fun_best,qsm_tree);
        %cd(out_name);% from now on we are in the final results directory
        obj = VideoWriter('gaMov.avi');
        obj.FrameRate = 0.5;% 2 sec between frames, nice looking
        open(obj);
        writeVideo(obj,mov);
        close(obj);
        movefile('gaMov.avi',out_name);
    else
        %%%
        qsm_tree = qsm_tree.move_tree([0 0 0]);
        [xmin0,xmax0,ymin0,ymax0,zmax0] = span(qsm_tree);
        ssm_scatter = ssm_fun_best(X);
        ssm_tree = read_mtg('out.mtg');% must be encoded in config
        ssm_tree = ssm_tree.move_tree([0 0 0]);
        [xmin1,xmax1,ymin1,ymax1,zmax1] = span(ssm_tree);
        %%%
        figure('position',[100 100 1000 500]);
        % XZ-view
        subplot(211);
        qsm_tree.draw;
        view(0,0);% SIDE view of the data
        xlim([min(xmin0,xmin1) max(xmax0,xmax1)]);
        zlim([0 max(zmax0,zmax1)]);
        subplot(212);
        new_ssm_tree = ssm_tree;
        max_o = 4;
        while( new_ssm_tree.number_of_branches > 100000 )
            max_o = max_o - 1;
            new_ssm_tree = branches_by_order(new_ssm_tree,max_o);
        end
        new_ssm_tree.draw;title(['w <= ' num2str(max_o)]);
        view(0,0);% SIDE view of the model
        xlim([min(xmin0,xmin1) max(xmax0,xmax1)]);
        zlim([0 max(zmax0,zmax1)]);
        saveas(gcf,'trees-xz.png');
        % XY-view
        subplot(211);
        view(0,90);% TOP view of the data
        xlim([min(xmin0,xmin1) max(xmax0,xmax1)]);
        ylim([min(ymin0,ymin1) max(ymax0,ymax1)]);
        subplot(212);
        view(0,90);% TOP view of the model
        xlim([min(xmin0,xmin1) max(xmax0,xmax1)]);
        ylim([min(ymin0,ymin1) max(ymax0,ymax1)]);
        saveas(gcf,'trees-xy.png');
        % Scatter plot
        plot_scatter(qsm_scatter,ssm_scatter);
        saveas(gcf,'scatters.png');
        %%%
        movefile('trees-xz.png',out_name);
        movefile('trees-xy.png',out_name);
        movefile('scatters.png',out_name);
        %cd(out_name);% we are in the final results folder
    end
else% No plot for trees, but scatters
    ssm_scatter = ssm_fun_best(X);
    plot_scatter(qsm_scatter,ssm_scatter);
    saveas(gcf,'scatters.png');
    movefile('scatters.png',out_name);
end
%% Finalize
% Moving gaOut.dat
movefile(config.ga_out_dat,out_name);
% Save the results
save('bf_out.mat','X','ssm_fun','ssm_fun_best','Problem','Output');
movefile('bf_out.mat',out_name);
% Move the input configuration file renaming it to the original name
movefile(tmp_input_file,[out_name '/' input_file_name]);
% Save genetic algorithm figure
hall = findall(0,'type','figure');
saveas(hall(strcmp({hall(:).Name},'Genetic Algorithm')),'ga');
movefile('ga.fig',out_name);
fprintf('The output directory:\n%s\n',[pwd '/' out_name]);% report the results folder
% Return to the original directory
cd(currDir);
end

function [xmin,xmax,ymin,ymax,zmax] = span(tr)
% Find the span in directions
ADDUP = 0.2;
if(min(tr.end_point(:,1)) < min(tr.start_point(:,1)))
    [~,xmin] = min(tr.end_point(:,1));
    xmin = tr.end_point(xmin,1) - ADDUP;
else
    [~,xmin] = min(tr.start_point(:,1));
    xmin = tr.start_point(xmin,1) - ADDUP;
end
if(max(tr.end_point(:,1)) > max(tr.start_point(:,1)))
    [~,xmax] = max(tr.end_point(:,1));
    xmax = tr.end_point(xmax,1) + ADDUP;
else
    [~,xmax] = max(tr.start_point(:,1));
    xmax = tr.start_point(xmax,1) + ADDUP;
end

if(min(tr.end_point(:,2)) < min(tr.start_point(:,2)))
    [~,ymin] = min(tr.end_point(:,2));
    ymin = tr.end_point(ymin,2) - ADDUP;
else
    [~,ymin] = min(tr.start_point(:,2));
    ymin = tr.start_point(ymin,2) - ADDUP;
end
if(max(tr.end_point(:,2)) > max(tr.start_point(:,2)))
    [~,ymax] = max(tr.end_point(:,2));
    ymax = tr.end_point(ymax,2) + ADDUP;
else
    [~,ymax] = max(tr.start_point(:,2));
    ymax = tr.start_point(ymax,2) + ADDUP;
end

if(max(tr.end_point(:,3)) > max(tr.start_point(:,3)))
    [~,zmax] = max(tr.end_point(:,3));
    zmax = tr.end_point(zmax,3) + ADDUP;
else
    [~,zmax] = max(tr.start_point(:,3));
    zmax = tr.start_point(zmax,3) + ADDUP;
end
end

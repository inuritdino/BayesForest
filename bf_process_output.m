function out = bf_process_output(folder,varargin)
% Process the output of BayesForest.
% USAGE:
%       OUT = BF_PROCESS_OUTPUT(FOLDER,...)
% FOLDER: the output folder created after BayesForest execution.
% OUT: output struct with fields:
%   ssm: the SSM function handle calling the growth simulation.
%   x0: the best-fit parameter vector of the SSM.
%   Um: the model/SSm data sets Um requested in the BayesForest input
%       configuration
%   Ud: the data/QSM Ud requested in BayesForest input configuration
%   D: the best distance reached during optimization
%   Dall: vector of individual distances from scatter vs. scatter
%   comparisons (usually, D is mean(Dall))
%   td: the data/QSM tree structure stored and manipulated via the tree-class
%   tm: the model/SSM tree structure stored and manipulated via the tree-class
%
% VARARGIN:
%   'inputFile': the input file name to use, default 'input.txt'
%   'gaOutFile': the output file produced by the optimization with Genetic
%       Algorithm by means of the function OPTIM_BEST_OUTPUT(). Default is
%       'gaOut.dat'.
%   'gen': the generation to extract from gaOutFile, default is the last,
%       i.e. the best score producing generation.
%
% See also optim_best_output

%% Initials
input_f = 'input.txt';
gaOut_f = 'gaOut.dat';
gen = 0;% the best/last generation by default

%% Input analysis
tf = strcmpi('inputFile',varargin);
if(find(tf))
    input_f = varargin{find(tf)+1};
end
tf = strcmpi('gaOutFile',varargin);
if(find(tf))
    gaOut_f = varargin{find(tf)+1};
end
tf = strcmpi('gen',varargin);
if(find(tf))
    gen = varargin{find(tf)+1};
end

%% Get the best function handle
c = bf_process_input([folder '/' input_f]);
ssm = str2func(c.ssm_fun_best);

%% Get the best/intermediate score and associated parameter vector
fid = fopen([folder '/' gaOut_f]);
last = [];
while 1
    line = fgetl(fid);
    if ~ischar(line), break, end
    last = line;
    if(gen)
        scan = textscat(last,'%f');
        if(scan{1}(1) == gen)
            x0 = scan{1}(4:end);
        end
    end
end
if(~gen)
    scan = textscan(last,'%f');
    x0 = scan{1}(4:end)';
end

disp(['Gen: ' num2str(scan{1}(1)) ' extracted.']);
fclose(fid);

%% Get the data sets Um and the tree structure
Um = ssm(x0);%NOTE: this command MUST produce out.mtg file!
tm = read_mtg('out.mtg');

%% Get the Ud 
if(~isempty(c.qsm_mat_file))
    [qsm_bra,qsm_seg,qsm_tree] = gen_scatter2(c.qsm_mat_file);
elseif(~isempty(c.qsm_cyl_table) && ~isempty(c.qsm_br_table))
    import_qsm_data(evalin('base',c.qsm_br_table),evalin('base',c.qsm_cyl_table),'temp.mat');
    [qsm_bra,qsm_seg,qsm_tree] = gen_scatter2('temp.mat');
    delete('temp.mat');
else
    error('Error: data (QSM) file is not specified.');
end
if(isempty(c.segment) && isempty(c.branch))
    if(c.qsm_merge)
        Ud = arrange_scatter(qsm_bra,qsm_seg,'scat',scat,...
            'order',order,'merge');
    else
        Ud = arrange_scatter(qsm_bra,qsm_seg,'scat',scat,...
            'order',order);
    end
else
    if(c.qsm_merge)
        Ud = arrange_scatter(qsm_bra,qsm_seg,'branch',c.branch,...
            'segment',c.segment,'merge');
    else
        Ud = arrange_scatter(qsm_bra,qsm_seg,'branch',c.branch,...
            'segment',c.segment);
    end
end
qsm_tree = qsm_tree.move_tree([0 0 0]);
%% Get the distances
[D,Dall] = optim_avg_distance(Um,Ud,'dirs',c.dt_dirs,'stat',c.dt_stat1d,...
    'smooth',0,'scaling',c.dt_scale,'w',c.dt_w);
N = length(Dall);
if(isempty(c.dt_w))
    w = repmat(1/N,1,N);% equal weights by default
else% rescale the weights, which were supplied by the user
    w = c.dt_w;
    w = w(:)./sum(w(:));% normalize the weights
end
Dall = Dall./w;

%% Output struct
out = struct('ssm',[],'x0',[],'Um',[],'Ud',[],'D',[],'Dall',[],'td',[],'tm',[]);
out.ssm = ssm;
out.x0 = x0;
out.Um = Um;
out.Ud = Ud;
out.D = D;
out.Dall = Dall;
out.td = qsm_tree;
out.tm = tm;


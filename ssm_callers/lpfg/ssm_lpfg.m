function [Scatter , Overall, Outfname] = ssm_lpfg(X,varargin)
% This is the MODEL caller for the Stochastic Structural Model (SSM)
% implemented in the LPFG simulator. The function calls the model and
% collects the STRUCTURAL data output from it. The data must be organized
% in a certain way (see read_scatter_dat2() function and scatter/data
% generation functions in some SSM's LPFG implementations, like LIGNUM-LPFG
% and SelfOrganizingTrees).
%
% USAGE:
%      [SCATTER,OVERALL,OUTFNAME] = SSM_LPFG(X,...)
% SCATTER is the output scatter(s). All scatters are reduced in two large
% data sets called BRANCH and SEGMENT data sets (see read_scatter_dat2(), 
% for example)
%
% OVERALL is the overall information about the generated tree: height, 
% diameters etc.
%
% OUTFNAME is the unique output folder name one wants to save the results to.
% Usually called at the end of the optimization procedure to fixate the
% best solution for the subsequent analysis.
%
% X is the vector of input parameter values. The names of the parameters to
% vary are specified in 'ARGS' argument to varargin (see below). To know
% the parameter names one needs to look into specific LPFG implementation
% of the model. Description as to how to write the LPFG models susceptible
% to the parameter variations through this interface can be found
% elsewhere.
% 
% VARARGIN inputs:
%
% 'outf' - create the unique output filename.
%
% 'visual' - display the tree from within LPFG simulator. When iterating in
% the optimization, one does not need that.
%
% 'argsConst' - meta-description of the parameter values that are constant
% throughout the simulation. See the simple description of the
% meta-language below in the head of the function ARGS_TO_SSMARGS(...).
%
% 'C' - array of the constant values described by the argsConst.
%
% 'args' - the meta-description of the values in X.
%
% 'scat' - name of the target scatter to obtain from the simulation. The
% two types are supported: 'Branch' and/or 'Segment'. 
% See READ_SCATTER_DAT2() for the details.
%
% 'order' - vector of branch orders for the scatters.
%
% 'rndseed' - the additional possibility to specify seed to the SSM
% (i.e. RNDSEED parameter). It will overwrite the normal usage of RNDSEED
% through the argsConst. RNDSEED parameter is used to control the
% randomness of the model.

%% Initializations
outf_flag = 0;
visual = 0;
args = {};
argsConst = {};
ConstVals = [];
Overall = [];
rndseed = [];
% Default output is empty.
Scatter = [];
% Max order preallocation for read_scatter_dat2.
max_order = 10;
% Orders to consider
scat_order = 1;
% Type of scatter
scat_type = {'branch','segment'};
% Fine tuning of the scatters
branch_type_orders = [];
segment_type_orders = [];
%% Reading input
ii = 1;
while ( ii <= length(varargin) )
    if(strcmp('outf',varargin{ii}))
        outf_flag = 1;
    elseif(strcmp('visual',varargin{ii}))
        visual = 1;
    elseif(strcmp('rndseed',varargin{ii}))
        rndseed = varargin{ii+1};
        ii = ii + 1;
    elseif(strcmp('argsConst',varargin{ii}))
        argsConst = varargin{ii+1};
        if(~isa(argsConst,'cell'))
            argsConst = {argsConst};
        end
        ii = ii + 1;
    elseif(strcmp('C',varargin{ii}))
        ConstVals = varargin{ii+1};
        ii = ii + 1;
    elseif(strcmpi('MaxOrder',varargin{ii}))
        max_order = varargin{ii+1};
        ii = ii + 1;
    elseif(strcmpi('order',varargin{ii}))
        scat_order = varargin{ii+1};
        ii = ii + 1;
    elseif(strcmp('scat',varargin{ii}))
        scat_type = varargin{ii+1};
        ii = ii + 1;
    elseif(strcmp('branch',varargin{ii}))
        branch_type_orders = varargin{ii+1};
        ii = ii + 1;
    elseif(strcmp('segment',varargin{ii}))
        segment_type_orders = varargin{ii+1};
        ii = ii + 1;
    elseif(strcmp('args',varargin{ii}))
        args = varargin{ii+1};
        if(~isa(args,'cell'))
            args = {args};
        end
        ii = ii + 1;
    end
    ii = ii + 1;
end

%% Print info of the input parameter values
fprintf('Simulating with: \n');
for ii=1:length(X)
    fprintf('X(%d) = %G\n',ii,X(ii));
end

%% Process the meta-language
ssmmArgsConst = args_to_ssmArgs(ConstVals,argsConst);
if(~isempty(ssmmArgsConst))
    ssmArgs = [args_to_ssmArgs(X,args) ssmmArgsConst];
else
    ssmArgs = args_to_ssmArgs(X,args);
end
% Process rndseed separately. Note this will overwrite RNDSEED in argsConst
if(~isempty(rndseed))
    ssmArgs = [ssmArgs {'RNDSEED'} rndseed];
end
%disp(lgmArgs);

%% If we need an output to save the optim results, this will get a unique name
if(outf_flag)
    t = clock;
    Outfname = [num2str(t(3)) '.' num2str(t(2)) '.' num2str(t(1)) '_' ...
        num2str(t(4)) '.' num2str(t(5))];
    fprintf('Output: %s\n',Outfname);
    mkdir(Outfname);
    % Save the arguments to the specified folder
    save([Outfname '/args.mat'],'args','argsConst','ConstVals');
    % 
end

%% RUN SSM in LPFG
lpfg_run(visual,0,ssmArgs{:});

% Determine whether everything was smooth
fid = fopen('out.log','r');
ok = true;
while 1
    line = fgetl(fid);
    if(str2double(line) > 0)
        ok = false;
    end
    if(~ischar(line)), break; end
end
fclose(fid);

% Return from the function if simulation was not OK
if(~ok)
    return;
end

%% Read the overall info
Overall = read_info_dat('info.dat');

%% Generate scatter from the output file
[Branch, Segment] = read_scatter_dat2('scatter.dat','MaxOrder',max_order);

if(~isempty(branch_type_orders) || ~isempty(segment_type_orders))
    Scatter = arrange_scatter(Branch,Segment,'branch',branch_type_orders,...
        'segment',segment_type_orders);
else
    Scatter = arrange_scatter(Branch,Segment,'order',scat_order,'scat',scat_type);
end

end

%% Meta-language processing function
function ssmArgs = args_to_ssmArgs(vals,args)
% Simple transfer function from parameter description meta language to the
% LPFG SSM notation. Meta description: PARNAME_DISTR_PARNUM, i.e. a string
% with PARNAME - par name, DISTR - distribution it follows (empty in case
% of a constant), and PARNUM is the par number (e.g., for gauss 1 is mean, 2 is
% std). Two distributions are supported at the moment: GAUSS and UNIFORM.
% Examples: 
% LR_GAUSS_2_3: par name is LR, following gauss with mean and std as 2nd
% and 3rd values in VALS.

jj = 1;
for ii = 1:length(args)
    C = textscan(args{ii},'%s','delimiter','_');
    ssmArgs{jj} = C{1}{1};% PARNAME
    jj = jj + 1;
    if(isempty(C{1}{2}))% No distribution, just const
        ssmArgs{jj} = vals(str2double(C{1}{3}));
    else% some distribution
        if(strcmp('GAUSS',C{1}{2}))% Gauss
            ssmArgs{jj} = ['gauss(' num2str(vals(str2double(C{1}{3}))) ...
                ',' num2str(vals(str2double(C{1}{4}))) ')'];
        elseif(strcmp('UNIFORM',C{1}{2}))
            ssmArgs{jj} = ['uniform(' num2str(vals(str2double(C{1}{3}))) ...
                ',' num2str(vals(str2double(C{1}{4}))) ')'];
        end
    end
    jj = jj + 1;
end
if(isempty(args))
    ssmArgs = [];
end


end
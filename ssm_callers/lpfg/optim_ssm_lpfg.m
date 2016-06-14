function [scatOut,overall,outfname] = optim_ssm_lpfg(X,varargin)
% This is the MODEL caller for the a Stochastic Structural Model (SSM)
% implemented in the LPFG simulator.
% USAGE:
%      [s,ov,outf] = optim_ssm_lpfg(X,...)
% S is the target scatter (default: tapering function).
%
% OV is the overal information about the tree: height, diameters etc.
%
% OUTF is the output folder name that is unique one to save the results to.
% Usually called at the end of the optimization procedure to fixate the
% best solution.
%
% X is the vector of input parameter values.
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
% 'scat' - name of the target scatter to obtain from the simulation. See
% READ_SCATTER_DAT(...) for the supported scatters. The default is tapering
% distribution function. Additionally, the multiple scatters can be
% specified by putting the names into a cell-array. For example,
% OPTIM_SSM_LPFG(...,'scat',{'taper','bra'},...) would produce two
% scatters: tapering and branching angle distribution functions. Similarly,
% larger number of scatters can be used. New scatters can be introduced by
% modifying the SSM and READ_SCATTER_DAT().
%
% 'order' - order of the scatter determined by 'scat'. See
% READ_SCATTER_DAT(...).
%
% 'rndseed' - the additional possibility to specify seed to the SSM
% (i.e. RNDSEED parameter). It will overwrite the normal usage of RNDSEED
% through the argsConst.

%% Initializations
outf_flag = 0;
visual = 0;
args = {};
argsConst = {};
ConstVals = [];
scat_type = {'taper'};% default scatter type is tapering.
scat_order = 1;
scatOut = cell(1,length(scat_type));
overall = [];
rndseed = [];
trunk_scatter_flag = 0;

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
            %error('Error: argsConst must be cell array.');
        end
        ii = ii + 1;
    elseif(strcmp('C',varargin{ii}))
        ConstVals = varargin{ii+1};
        ii = ii + 1;
    elseif(strcmp('args',varargin{ii}))
        args = varargin{ii+1};
        if(~isa(args,'cell'))
            args = {args};
            %error('Error: args must be cell array.');
        end
        ii = ii + 1;
    elseif(strcmp('order',varargin{ii}))
        scat_order = varargin{ii+1};
        ii = ii + 1;
    elseif(strcmp('scat',varargin{ii}))
        scat_type = varargin{ii+1};
        if(isa(scat_type,'char'))% Make the cell array if not
            scat_type = {scat_type};
        end
        scatOut = cell(1,length(scat_type));
        ii = ii + 1;
    elseif(strcmp('trunk',varargin{ii}))
        trunk_scatter_flag = 1;
    end
    ii = ii + 1;
end

%% Print info of the input parameter values
fprintf('Simulating with: \n');
for ii=1:length(X)
    fprintf('X(%d) = %G\n',ii,X(ii));
end

%% Process the meta-language
lgmArgsConst = args_to_ssmArgs(ConstVals,argsConst);
if(~isempty(lgmArgsConst))
    lgmArgs = [args_to_ssmArgs(X,args) lgmArgsConst];
else
    lgmArgs = args_to_ssmArgs(X,args);
end
% Process rndseed separately. Note this will overwrite RNDSEED in argsConst
if(~isempty(rndseed))
    lgmArgs = [lgmArgs {'RNDSEED'} rndseed];
end
%disp(lgmArgs);

%% If we need an output to save the optim results, this will get a unique name
if(outf_flag)
    t = clock;
    outfname = [num2str(t(3)) '.' num2str(t(2)) '.' num2str(t(1)) '_' ...
        num2str(t(4)) '.' num2str(t(5))];
    fprintf('Output: %s\n',outfname);
    mkdir(outfname);
    % Save the arguments to the specified folder
    save([outfname '/args.mat'],'args','argsConst','ConstVals');
    % 
end

%% RUN LPFG-LGM
lpfg_run(visual,0,lgmArgs{:});

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
    scatOut = []; % make the output empty to track the error
    return;
end

%% Read the overall info
overall = read_info_dat('info.dat');

%% Generate scatter from the output file
[scatter,scat_names] = read_scatter_dat('scatter.dat','MaxOrder',max(scat_order));

% if(isa(scat_type,'char'))% Make the cell array if not
%     scat_type = {scat_type};
% end
% 
% scatOut = cell(1,length(scat_type));
for jj = 1:length(scat_order)
    for ii = 1:length(scat_type)
        idx = strcmp(scat_type{ii},scat_names);
        if(all(~idx))
            %error('Error: the scatter name ** %s ** does not exist.',scat_type{ii});
            scatOut{ii} = [];
        else
            scatOut{(jj-1)*length(scat_type)+ii} = (scatter{idx}{scat_order(jj)+1})';
        end
    end
end
if(trunk_scatter_flag)
    m = length(scatOut);
    idx_tap = strcmp('taper',scat_names);
    idx_curv = strcmp('curv',scat_names);
    scatOut{m+1} = (scatter{idx_tap}{1})';
    scatOut{m+2} = (scatter{idx_curv}{1})';
end

end

%% Meta-language processing function
function lgmArgs = args_to_ssmArgs(vals,args)
% Simple transfer function from parameter description meta language to the
% LPFG LGM notation. Meta description: PARNAME_DISTR_PARNUM, i.e. a string
% with PARNAME - par name, DISTR - distribution it follows (empty in case
% of a constant), and PARNUM is the par number (e.g., for gauss 1 is mean, 2 is
% std).
% Examples: 
% LR_GAUSS_2_3: par name is LR, following gauss with mean and std as 2nd
% and 3rd values in VALS.

jj = 1;
for ii = 1:length(args)
    C = textscan(args{ii},'%s','delimiter','_');
    lgmArgs{jj} = C{1}{1};% PARNAME
    jj = jj + 1;
    if(isempty(C{1}{2}))% No distribution, just const
        lgmArgs{jj} = vals(str2double(C{1}{3}));
    else% some distribution
        if(strcmp('GAUSS',C{1}{2}))% Gauss
            lgmArgs{jj} = ['gauss(' num2str(vals(str2double(C{1}{3}))) ...
                ',' num2str(vals(str2double(C{1}{4}))) ')'];
        elseif(strcmp('UNIFORM',C{1}{2}))
            lgmArgs{jj} = ['uniform(' num2str(vals(str2double(C{1}{3}))) ...
                ',' num2str(vals(str2double(C{1}{4}))) ')'];
        end
    end
    jj = jj + 1;
end
if(isempty(args))
    lgmArgs = [];
end


end
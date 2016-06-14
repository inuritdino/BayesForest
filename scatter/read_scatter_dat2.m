function [Branch,Segment] = read_scatter_dat2(fn,varargin)
% Read the structured scatter data from LPFG-LIGNUM or other Stochastic
% Structure Model (SSM).
% USAGE:
%       [BRANCH, SEGMENT] = READ_SCATTER_DAT2(FN)
%
% CONVENTION: The new convention implies only two types of data sets:
% Branch and Segment based.
% *** Branch data sets describe the whole-branch characteristics such as 
% branching and azimuthal angles, base radius of the branch, its total 
% length etc.
% *** Segment data sets describe the segment based characteristics of the
% branches, where segments are the smaller elements of the branch. The
% Segment data set consists, for example, of tapering of the branches, 3D
% curvature etc.
%
% BRANCH - Branch data set: a struct with fields .data and .info. Info
% field describes the columns in the data. Data field contains the actual
% data points. Data field is a cell array, where each cell contains the
% data points for the corresponding order of the tree. Namely,
% Branch.data{1} is for order = 0(trunk), Branch.data{2} - for order = 1
% etc.
%
% SEGMENT - Segment data set: a struct with fields .data and .info. Info
% field describes the columns in the data. Data field contains the actual
% data points. Data field is a cell array, where each cell contains the
% data points for the corresponding order of the tree. Namely,
% Segment.data{1} is for order = 0(trunk), Segment.data{2} - for order = 1
% etc. 
%
% FN is the file name for the scatter data.
%
% ... = READ_SCATTER_DAT2(...,'VARARGIN','MaxOrder') specifies the maximum
% order of the data. It is used only for the preallocation purposes, i.e.
% to avoid memory constraints. Default: 10.

% Find out the order
MAXORDER = 10;% now MAXORDER is served for the preallocation purposes only
tf = strcmpi('MaxOrder',varargin);
if(find(tf))
    MAXORDER = varargin{find(tf)+1};
end
fid = fopen(fn);

% Initiate the output
Branch = struct('info','data');
Branch.data = cell(1,MAXORDER+1);
Branch.info = {};

Segment = struct('info','data');
Segment.data = cell(1,MAXORDER+1);
Segment.info = {};

order = [];
type = 0; % 0 - segment data set, 1 - branch data set.
nl_count = 0;
while 1
    line = fgetl(fid);% get line
    if(~ischar(line)), break; end;% EOF break
    if(isempty(line)) 
        nl_count = nl_count + 1;
        if(nl_count == 2)
            type = 0;% change type after 2 consecutive newlines.
            nl_count = 0;
        end
        continue;
    end
    if(line(1) == '#')% comment or instruction
        C = textscan(line,'%s','delimiter',' ');
        if(strcmp(C{1}{2},'Branch:'))
            Branch.info = C{1}(3:end)';
        elseif(strcmp(C{1}{2},'Segment:'))
            Segment.info = C{1}(3:end)';
        elseif(strcmp(C{1}{2},'order'))
            order = str2double(C{1}(3));
            if(order > MAXORDER)
                Branch.data{order+1} = [];
                Segment.data{order+1} = [];
            end
            % First data set type after order has been selected
            if(order > 0)
                % order=0 has only Segment data
                % order>0: Branch data sets go first
                type = 1;
            end
        end
    else% data string
        C = textscan(line,'%f','delimiter',' ');
        if(isempty(order))
            error('Error: order must be specified.\n');
        elseif(order == 0)
            Branch.data{order+1} = [];
            Segment.data{order+1} = cat(1,Segment.data{order+1},C{1}(:)');
        else
            if(type == 1)
                Branch.data{order+1} = cat(1,Branch.data{order+1},C{1}(:)');
            end
            if(type == 0)
                Segment.data{order+1} = cat(1,Segment.data{order+1},C{1}(:)');
            end
        end
        
    end
end
% Max order detected
%disp(['Max order detected: ' num2str(order) ]);

% Clean the data cell-arrays
Branch.data = Branch.data(1:order+1);
Segment.data = Segment.data(1:order+1);

fclose(fid);

end

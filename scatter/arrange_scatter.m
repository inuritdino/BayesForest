function Scatter = arrange_scatter(Branch,Segment,varargin)
% Arranges the scatters according to the order and scatter type. Two types
% of scatters are supported: Branch and Segment based data sets. Use with
% function gen_scatter2().
% USAGE:
% 
%       SCATTER = ARRANGE_SCATTER(BRANCH,SEGMENT,VARARGIN)
%
% SCATTER -- arranged scatter.
% 
% BRANCH -- Branch data set.
%
% SEGMENT -- Segment data set.
%
% VARARGIN inputs:
%
% 'scat' - name of the target scatter to obtain from the simulation. The
% two types are supported: 'Branch' and/or 'Segment'. 
% See READ_SCATTER_DAT2() for the details.
%
% 'order' - vector of branch orders for the scatters.
%
% The output order is Order1-ScatType1-ScatType2,
% Order2-ScatType1-ScatType2, ... if two scatter types are specified or
% similarly for the one-type arrangement.
%
% Another alternative to specify the scatters with orders is to use
% 'branch' and/or 'segment' varargin's. This take precedence over the
% tradition 'scat'/'order' way. This way also gives more flexibility.
%
% 'branch' - array of order values for the Branch data sets.
%
% 'segment' - array of order values for the Segment data sets.
%
% NOTE: NaN's are removed from the output SCATTER.


% Default output is empty.
Scatter = [];
% Orders to consider
scat_order = 1;
% Type of scatter
scat_type = {'branch','segment'};
% Order and types: more customizable option
branch_type_orders = [];
segment_type_orders = [];
% Merge orders: merging all-order segment and branch scatters: max two
% huge tables to the output containing scatters of all specified orders
merge_flag = false;

ii = 1;
while ( ii <= length(varargin) )
    if(strcmpi('order',varargin{ii}))
        scat_order = varargin{ii+1};
        ii = ii + 1;
    elseif(strcmp('scat',varargin{ii}))
        scat_type = varargin{ii+1};
        if(isa(scat_type,'char'))% Make the cell array if not cell array
            scat_type = {scat_type};
        end
        idx1 = strcmpi('branch',scat_type);
        idx2 = strcmpi('segment',scat_type);
        if(length(scat_type) > 2 || (all(~idx1) && all(~idx2)) || ...
                ~all(idx1 | idx2))
            error('Error: scatter types are two: branch and/or segment.');
        elseif(length(scat_type) == 2 && strcmp(scat_type{1},scat_type{2}))
            scat_type = {scat_type(1)};
        end
        ii = ii + 1;
    elseif(strcmp('branch',varargin{ii}))
        branch_type_orders = varargin{ii+1};
        ii = ii + 1;
    elseif(strcmp('segment',varargin{ii}))
        segment_type_orders = varargin{ii+1};
        ii = ii + 1;
    elseif(strcmp('merge',varargin{ii}))
        merge_flag = true;
    end
    ii = ii + 1;
end

% Check for availability of the orders
if(((max(scat_order)+1) > length(Branch.data)) || ...
        ((max(scat_order)+1) > length(Segment.data)) || ...
        (~isempty(branch_type_orders) && (max(branch_type_orders)+1) > length(Branch.data)) || ...
        (~isempty(segment_type_orders) && (max(segment_type_orders)+1) > length(Segment.data)))
    fprintf('Warning: requested order(s) are not available.\n');
    return;
end

% Allocate the scatter
if(~isempty(branch_type_orders) || ~isempty(segment_type_orders))
    scat_arrange = 1;
    % Specific orders for specific scatters
    BOS = length(branch_type_orders);
    SOS = length(segment_type_orders);
    Scatter = cell(1,BOS+SOS);
    for ii = 1:BOS
        Scatter{ii} = Branch.data{branch_type_orders(ii)+1}';
    end
    for ii = 1:SOS
        Scatter{BOS+ii} = Segment.data{segment_type_orders(ii)+1}';
    end
else% specified orders for all specified scatters
    scat_arrange = 2;
    NS = length(scat_type);
    Scatter = cell(1,NS*length(scat_order));
    for ii = 1:length(scat_order)
        for jj = 1:NS
            if(strcmpi(scat_type{jj},'branch'))
                Scatter{(ii-1)*NS+jj} = Branch.data{scat_order(ii)+1}';
            elseif(strcmpi(scat_type{jj},'segment'))
                Scatter{(ii-1)*NS+jj} = Segment.data{scat_order(ii)+1}';
            end
        end
    end
end

% Remove NaN's and Inf's which obstruct the analysis. NOTE: Scatter{ii} is 
% a matrix with fixed number of rows (variables) and varying number of 
% columns (samples). So, identify the columns with NaN's or Inf's and 
% remove them.
for ii = 1:length(Scatter)
    [r,~] = size(Scatter{ii});
    nans_infs = find(isnan(Scatter{ii}(:)) | isinf(Scatter{ii}(:)));
    %target_rows = mod(nans_infs,r);
    %target_rows(target_rows == 0) = r;
    target_cols = ceil(nans_infs./r);
    Scatter{ii}(:,target_cols) = [];% remove target columns
end

% Merging the scatters.
if(merge_flag)
    Scatter_merged = cell(1,2);
    if(scat_arrange == 1)
        Scatter_merged{1} = cat(2,Scatter_merged{1},Scatter{1:BOS});
        Scatter_merged{2} = cat(2,Scatter_merged{2},Scatter{BOS+1:BOS+SOS});
    elseif(scat_arrange == 2)
        if(NS == 2)
            idx = 1:2:NS*length(scat_order);
            Scatter_merged{1} = cat(2,Scatter_merged{1},Scatter{idx});
            idx = 2:2:NS*length(scat_order);
            Scatter_merged{2} = cat(2,Scatter_merged{2},Scatter{idx});
        elseif(NS == 1)
            Scatter_merged{1} = cat(2,Scatter_merged{1},Scatter{:});
            Scatter_merged = Scatter_merged(1);
        end
    end
    if(isempty(Scatter_merged{1}))
        Scatter_merged = Scatter_merged(2);
    elseif(isempty(Scatter_merged{2}))
        Scatter_merged = Scatter_merged(1);
    end
    Scatter = Scatter_merged;
end
end
function varargout = branch_distribution(obj,varargin)
% Compute and return branch distribution data of a tree.
% The function returns the branch length and/or volume 
% distributions.
%
% dists return value is a matrix whos first row contains the
% radius distribution. The following rows contain
% corresponding distributions such as the total length of
% branches with the given radius.
%
% The function understands the following additional
% parameters which are given as name-value -pairs. The names
% are not case-sensitive.
%
% 'Radius' or  defines the minimum radius. Branches with 
% 'MinRadius'  radius under this value are grouped together. 
%              By default the minimum radius is set to 5cm.
%              To disable minimum radius it must be set to
%              zero.
% 
% 'Precision' defines the accuracy of the analysis. By
%             default the value is 2 which means that the
%             analysis is done with a centimeters accuracy.
%             Similarly the value 3 would mean millimeter.
%
% 'Interval' defines either the number of subintervals or
%            the upper bounds of the intervals to be used.
%            If the value is scalar the region between the
%            smallest and larges radius is divided into that
%            many subintervals using the linspace-command.
%            If the value is a vector then the elements of
%            that vector are used as the interval limits.
%            The elements are assumed to be sorted
%            ascendingly.
%
% 'Relative' return relative distribution instead of
%            absolute.
%
% Examples of usage:
%
% d = obj.branch_distribution('Radius',0.01,'Interval',15);
% d = obj.branch_distribution('Interval',[5 10 15 20]);
% d = obj.branch_distribution('Radius',0);
% d = obj.branch_distribution();
%

if obj.number_of_branches == 0
    error('No cylinders in the tree model.');
%elseif obj.number_of_branches == 1
%    error('Only one cylinder in the tree model.');
end

% Minimum radius. Branches with radius lower than this are
% all summed up.
r_min = 0.005;

% Precision of rounding. By default the unit of length is
% assumed to be a meter, so the default nd = 2 would mean
% the precision of centimeters. Note that in any case the
% possible minimum radius should still be given in meters.
nd = 2;

% Number of intervals. The variable determines to how many
% subintervals the radii are divided into. By default the
% value is zero which means that the precision unit is used
% as the interval width.
ni = 0;
intervals = [];

% Truth value which decides whether the distribution is
% absolute (false) or relative (true). By default the
% distribution is absolute.
relative = false;

vol = false;
len = true;
color = 'b';

% Number of properties to return.
k = 2;

% Inspection of optional arguments.
i = 1;
while i <= size(varargin,2)

    if ischar(varargin{i})

        switch lower(varargin{i})

            case {'minradius', 'radius'}
                r_min = varargin{i+1};
                i = i + 1;
            case 'precision'
                nd = varargin{i+1};
                i = i + 1;
            case 'intervals'
                if isscalar(varargin{i+1})
                    ni = varargin{i+1};
                else
                    intervals = varargin{i+1};
                    ni = length(intervals);
                end
                i = i + 1;
            case 'relative'
                relative = true;
            case 'volume'
                vol = true;
                len = false;
            case 'length'
                len = true;
                vol = false;
            case 'both'
                len = true;
                vol = true;
                k = 3;
            case 'color'
                color = varargin{i+1}(1);
                i = i + 1;
        end
    end
    i = i + 1;
end

% Matrix to which tree data is collected for further
% processing.
A = zeros(k,obj.number_of_branches);

% Collect data from the tree to the matrix.
A(1,:) = obj.radius(1:obj.number_of_branches);
A(2,:) = obj.length(1:obj.number_of_branches);

% To short hand further computations.
b = A(1,:);

% Round to given or default precision. N.B. nb+1
b = floor(10^(nd+1)*b);

% Radii smaller that r_min are set to zero.
b(b < r_min*10^(nd+1)) = 0;

% Divide by 10.
b = floor(b./10);

% Sort according to accending radius.
[b, I] = sort(b);

A = A(:,I);

sp = find(b,1,'first');
if isempty(sp)
    error('Minimum radius too high. No branches left.');
end

% Limits of the distribution subintervals. If intervals
% argument is not set, move from minimum to maximum with
% increments of one.
if ni == 0
    intervals = b(sp):b(end);
    ni = length(intervals);
elseif isempty(intervals)
    intervals = linspace(b(sp),b(end),ni);
else
    % If the interval maximum is smaller than the largest
    % actual value, add new interval.
    if intervals(end) < b(end)
        intervals = [intervals(:); b(end)]';
        ni = ni + 1;
    end
    
end

if length(intervals) > 1
    dint = (intervals(2)-intervals(1))/2;
else
    dint = 0.5;
end

% If the minumum radius is set.
if r_min ~= 0 && intervals(1) > r_min*10^(nd)
    intervals = [r_min*10^(nd); intervals(:)]';
    ni = ni + 1;
end

%First row of dists are the upper bounds of the radius
%subintervals. Rest of the matrix is initialized as zeros.
dists = [intervals; zeros(k-1,ni)];

i = 1; % Column of the property matrix.
j = 1; % Number of interval.

% Go through all the radii.
while i <= length(b)

    % If the radius is smaller than the current interval
    % bound do the proper calculations and move to the next
    % radius.
    if b(i) <= intervals(j) + (j~=1)*dint
        if len
            dists(2,j) = dists(2,j) + A(2,i);
        end
        if vol
            dists(end,j) = dists(end,j) + A(2,i)*A(1,i)^2*pi;
        end
        i = i + 1;
    
    % Otherwise the radius is larger than the bound so move
    % on to the next interval.
    else
        j = j + 1;
    end
end

% If relative distribution is required divide each element
% of every row except the first with the sum of the elements
% in that row. 
if relative
    
    dists(2:k,:) = bsxfun(@rdivide,dists(2:k,:),sum(dists(2:k,:),2));
end

% Select output according to number of output arguments.
switch nargout
    case 2
        varargout{1} = dists(1,:);
        varargout{2} = dists(2:end,:);
    case 1
        varargout{1} = dists;
    otherwise
        if r_min
            bar([0 dists(1,2:end)],dists(2:end,:),color);
        else
            bar(dists(1,:),dists(2:end,:),color);
        end
end
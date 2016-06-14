function [scatter, trunk, tree_out] = gen_scatter(input_data_file)
% USAGE: scatter = gen_scatter(input_data_file)
%
% The function computes the all neccessary to the moment scatters needed for
% the subsequent use in the bayes tree generation as the source for the
% distribution functions (DF). The prior knowledge of the certain
% characteristics (INPUT DATA) is needed. The characteristics come from the
% analysis of the measurements/simulations.
% -------- INPUT DATA ---------
% The following are cylinder related data. Each cylinder is indexed and the 
% properties of the cylinder with index "I", can be found from the Ith row
% of these inputs:
% Rad       Radii of the cylinders, (n_cyls x 1)-vector
% Len       Lengths of the cylinders, (n_cyls x 1)-vector
% Sta       Starting points of the cylinders (the axis point at the bottom
%               of the cylinder), (n_cyls x 3)-matrix
% Axe       Axis direction vectors (unit vectors) of the cylinders, 
%               (n_cyls x 3)-matrix
% CPar      Parent cylinders (the index of the parent cylinder), (n_cyls x 1)-vector
% Added     Logigal vector indicating if the cylinder is added to fill gaps
%               (n_cyls x 1)-vector
%
% The following are branch related data. Each branch is indexed and the 
% properties of the branch with index "I", can be found from the Ith row
% of these inputs:
% BOrd      Order of the branch (0 = trunk, 1 = branches originating from 
%               the trunk, etc.), (n_branch x 1)-vector
% BChi      Child branches, (n_branch x 1)-cell array
% CiB       Cylinders in the branches (gives the indexes of the cylinders 
%               in the branch), (n_branch x 1)-cell array
% BLen      Branch lengths (sum of the cylinder lengths), (n_branch x 1)-vector
% BAng      Angle between branch and its parent branch, (n_branch x 1)-vector
% FCB       Index of the first cylinder in the branch (not the "Added"
%               cylinder), (n_branch x 1)-vector
% ------------------------------
% The output is the structure scatter with the corresponding scatters in its
% fields. For example, scatter.idist gives the all data points for the
% 'inter-branch distances' scatter. The dimension of the scatter struct is
% the highest order of the tree analyzed, e.g. scatter(2).idist gives the
% inter-branch distances scatter for the 2-nd order branches.
% 
% Additionally, the function calculates all the characteristics needed for
% the trunk generation. Here, the trunk is the branch of order 0.
% 


% ----- INITIALIZATION -----

% Max order of the interest. Non-automatic.
MAX_W = 50;

% Initialization of the scatter structure with all its fields.
scatter(1) = struct('idist',[],'az',[],'azd_idist',[],'idist_thght',[],...
    'rchipar_rpar',[],'mo_thght',[],'curv',[],'tap',[],'ltot_rini',[],...
    'bra',[],'lchi_lapar',[],'lchi_bra_lapar',[]);
% just for initialization of the whole structure array
scatter(MAX_W).idist = [];

% Load the input data from the specified file
in = load(input_data_file);

% Check the all needed fields in the input
if(~isfield(in,'Rad'))
    error('Error: RAD is not found');
end
if(~isfield(in,'Len'))
    error('Error: LEN is not found');
end
if(~isfield(in,'Sta'))
    error('Error: STA is not found');
end
if(~isfield(in,'Axe'))
    error('Error: Axe is not found');
end
if(~isfield(in,'CPar'))
    error('Error: CPAR is not found');
end
if(~isfield(in,'Added'))
    error('Error: ADDED is not found');
end
if(~isfield(in,'BOrd'))
    error('Error: BORD is not found');
end
if(~isfield(in,'BChi'))
    if(isfield(in,'BPar'))% if we have BPar info, we use it to retrieve BChi
        for ii=1:length(in.BPar)% iterate over all branches
            ind = find(in.BPar == ii);% find branches whose parent is ii
            if(isempty(ind))% if no one is found, then an empty matrix
                in.BChi{ii} = [];
            else% put all found branches as chilren to ii.
                in.BChi{ii} = ind;
            end
        end
    else
        error('Error: BCHI is not found and cannot be resolved');
    end
end
if(~isfield(in,'CiB'))
    if(isfield(in,'BoC'))% if we have BoC, we can retrieve CiB
        for ii=1:length(in.BOrd)% iterate over all branches
            ind = find(in.BoC == ii);% find cylinders of the ii-th branch
            if(isempty(ind))% if not found, put an empty matrix
                in.CiB{ii} = [];
            else% put all found cylinders as belonging to the ii-th branch
                in.CiB{ii} = ind;
            end
        end
    else
        error('Error: CIB is not found and cannot be resolved');
    end
end
if(~isfield(in,'BLen'))
    error('Error: BLEN is not found');
end
if(~isfield(in,'BAng'))
    error('Error: BANG is not found');
end
if(~isfield(in,'FCB'))
    if(isfield(in,'BoC'))
        for ii=1:length(in.BOrd)% iterate over all branches
            if(~isempty(in.CiB{ii}))
                ind = in.CiB{ii}(1);% take the 1st in CiB
            end
            % Find the parent cyl. until we reach the start of the branch
            % The condition: must be the same branch and not zero value
            while((in.CPar(ind) ~= 0) && (in.BoC(in.CPar(ind)) == ii))
                ind = in.CPar(ind);% change ind to be the parent
            end
            in.FCB(ii) = ind;
        end
    else
        error('Error: FCB is not found and cannot be resolved');
    end
end

% convert the names for convenience
Rad = in.Rad;
Len = in.Len;
Sta = in.Sta;
Axe = in.Axe;
CPar = in.CPar;
Added = in.Added;
BOrd = in.BOrd;
BChi = in.BChi;
CiB = in.CiB;
BLen = in.BLen;
BAng = in.BAng;
FCB = in.FCB;

% Some constants
N_CYL = length(Rad);% number of cylinders
N_BR = length(BOrd);% number of branches
U = [1 0 0]';% reference branch direction
V = [0 1 0]';% reference vector for gamma (angle in xy-plane)
W = [0 0 1]';% reference vector for zeta (angle in xz-plane)

% ----- END of INITIALIZATION -----

% ----- TRUNK INFORMATION EXTRACT -----

trunk_ind = find(BOrd == 0);% branch index of the trunk
if(isempty(trunk_ind))
    trunk = [];
else
    trunk = struct('tree',[],'scatter',[]);
    trunk.tree = tree();
    trunk.tree.number_of_branches = CiB{trunk_ind}(end) - FCB(trunk_ind) + 1;
    trunk.tree.radius = Rad(FCB(trunk_ind):CiB{trunk_ind}(end));
    trunk.tree.length = Len(FCB(trunk_ind):CiB{trunk_ind}(end));
    trunk.tree.start_point = Sta(FCB(trunk_ind):CiB{trunk_ind}(end),:);
    trunk.tree.axis = Axe(FCB(trunk_ind):CiB{trunk_ind}(end),:);
    trunk.tree.end_point = trunk.tree.start_point + bsxfun(@times,trunk.tree.length,trunk.tree.axis);
    % The only statistics (if any) one can get from the trunk is for
    % 'taper' and 'curv' scatters
    trunk.scatter = struct('tap',[],'curv',[]);
    % TAP
    C = CiB{trunk_ind};% find all cylinders in the branch
    L = cumsum(Len(C));% cumulative sum of the cyl. lengths
    for ii = 1:length(C)
        trunk.scatter.tap = cat(1,trunk.scatter.tap,[L(ii) Rad(C(ii))]);
    end
    % CURV
    I = Added(C(1));% check if the 1st cylinder is Added one
    if I% if so...
        C = C(2:end);% start from the 2nd one
    end
    trunk.scatter.curv = curv(C,Axe,Len,W);
end

% ----- END of TRUNK INFORMATION EXTRACT -----

% ----- TREE CONVERT TO TREE-OBJECT -----
% Simple code to covert to the fields of tree-object. Do not care about the
% topology information, just for visualization purposes.
tree_out = tree();
tree_out.number_of_branches = length(Rad);
tree_out.radius = Rad;
tree_out.length = Len;
tree_out.start_point = Sta;
tree_out.axis = Axe;
tree_out.end_point = Sta + bsxfun(@times,Len,Axe);

% ----- END of TREE CONVERT TO TREE-OBJECT -----


% ----- SCATTERS CALCULATION ------

% ***1*** order-1 -> order-n branch map
% Needed for so called 'Trunk-height'(thght) DF's

br = find(BOrd == 1);% order-1 branches
map = [];% the instance of the map

% Find the immediate (order-2) children of the order-1 branches
for ii = 1:length(br)% iterate over 1-order branches
    child = BChi{br(ii)};% all the children of the iteration
    map = [map; repmat(br(ii),length(child),1) child];% change the map
end

% Associate the other descendants (order > 2) with the order-1 branches
for j = 2:MAX_W% iterate over the orders >= 2
    br = find(BOrd == j);% find j-th order branches
    for ii = 1:length(br)% iterate over the parent branches br
        child = BChi{br(ii)};% find children of the current parent br(ii)
        if(isempty(child))% no children were found, go on to the next
            continue;
        end
        % find the parent in the map
        br_in_map = find(br(ii) == map(:,2));
        %  br_in_map must be 1x1
        if(isempty(br_in_map))% not found in the map
            continue;
        elseif(length(br_in_map) ~= 1)
            disp('Warning: BR_IN_MAP is not scalar');
        end
        % assign 1-order branch to the children as for their parents
        map = [map; repmat(map(br_in_map,1),length(child),1) child];
    end
end

% ***2*** Inter-branch distance: DF of the relative distances between the
% child branch bases along the parent branch. CODE NAME: IDIST (1D)
% ***3*** Azimuth angle: DF of azimuth angles around the parent branch, in
% radians. CODE NAME: AZ (1D)
% ***4*** Radii ratio: the ratio between radius of the child at the base of
% the branch and radius of the parent (in place where the child emerges) 
% vs. the latter (parent radius). CODE NAME: RCHIPAR_RPAR (2D)
% ***5*** Azimuth vs. inter-branch distance: difference between azimuth 
% angles of the consecutive branches vs. the inter-branch distances. 
% CODE NAME: AZD_IDIST (2D)
% ***6*** Inter-branch vs. trunk height/thght: inter-branch distances vs. 
% the height along the trunk. CODE NAME: IDIST_THGHT (2D)
% ***7*** Max order vs. thght: max order of the tree branches vs. the trunk
% height. CODE NAME: MO_THGHT (2D)

% Collect azimuth and distance data

for k = 1:MAX_W
    m = nnz(BOrd == k);
    % some data collection matrices
    data = zeros(m,3);
    %dataR = zeros(m,2);
    dataRD = zeros(m,2);
    dataD = [];
    dataH = [];
    dataAD = [];
    % some helping counters
    t = 0;
    t0 = 1;
    n_par = 0;% number of the current order parents, which have a child
    for j = 1:N_BR% over all branches
        if BOrd(j) == k-1% A possible parent branch j for the current order k
            C = CiB{j};% Cylinders in the branch j, cylinder index
            Ind = 1:1:length(C);% Cylinder index list
            BC = BChi{j};% All child branches of the branch j
            n = length(BC);% Number of the child branches of the branch j
            for i = 1:n% iterate over the children of order k
                f = FCB(BC(i));% First cylinder in the child i
                % Find emergence point (cylinder) of child i on parent j
                I = CPar(f) == C;
                CP = C(I);% parent cylinder index (CP)
                I = Ind(I);% row number of the parent cylinder
                if isempty(CP)% empty CP: the parent is Added cylinder
                    if(~Added(CPar(f)))% if it is NOT Added, it is STRANGE!
                        disp('Strange!');
                    end
                    CP = CPar(CPar(f));% find parent of parent j of child i
                    I = CP == C;% parent cylinder index
                    I = Ind(I);% row number of the parent cylinder
                end
                % rotation vector from X-axis to parent axis (new x-axis)
                D = cross(Axe(CP,:),U);
                if norm(D) > 0
                    % angle between x-axis and new x-axis, y- and z-axis 
                    % are rotated for this angle as well
                    a = acos(Axe(CP,:)*U);
                    R = rotation_matrix(D,-a);% rotation matrix
                    v = R*V;% rotate V
                    w = R*W;% rotate W
                else% no rotation is needed
                    v = V;
                    w = W;
                end
                % Look into distances_to_line code. h is the segment of
                % Axe(CP,:) line between the Sta(CP,:) point and the point 
                % where the perpendicular from Sta(f,:) falls.
                % Normalization of the Axe(CP,:) takes place in the
                % function distances_to_line.
                [~,~,h] = distances_to_line(Sta(f,:),Axe(CP,:),Sta(CP,:));
                % y and z components, giving the direction in yz-plane
                a = Axe(f,:)*[v w];
                t = t+1;% increase counter
                % length along the parent, starting from the 1st cylinder
                data(t,1) = sum(Len(C(1:I-1)))+h;
                % normalize the length along the parent
                data(t,2) = data(t,1)/sum(Len(C));
                % find the azimuth, reference is v (new y-axis)
                data(t,3) = 180/pi*atan2(a(2),a(1));
                % initial radius of child i vs. distance along parent j
                %dataR(t,:) = [data(t,2) Rad(f)];
                % ratio between initial radius of child i and its parent
                % cylinder radius vs. radius of the parent cylinder
                dataRD(t,:) = [Rad(CP) Rad(f)/Rad(CP)];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Max order distribution.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Need to check only the order-1 child branches
                if(k == 1)
                    p = 1;% Max order p is 1, initially.
                    bchild = BC(i);% get the current child branch
                    while(p)% Iterate
                        p_old = p;% Save the old value of p
                        % Iterate over the target child branches: iterative
                        % process starting from the order-1 bchild and then
                        % progressively moving towards its descendants.
                        for q = 1:length(bchild)
                            % if there is a child of the higher order
                            % then increase p and break the for loop moving
                            % to the next order.
                            if(~isempty(BChi{bchild(q)}))
                                p = p+1;
                                break;
                            end
                        end
                        % if the p was not increased, break while loop
                        if(p == p_old), break; end
                        % create the new list of the target child branches.
                        bchild = BChi{bchild};
                    end
                    % Finally, write down data to the scatter.
                    dataO(t,:) = [data(t,2) p];
                    tf = t;% mark the final t, for the extraction later.
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Inter-branch difference DF's
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % If there was no child branch, t was not increased and t0 was.
            % So we need to subtract 1 from t0
            if(t0 == t+1)
                t0 = t0 - 1;
            elseif(t0 == t)
                % If ONE child branch was added then t0 == t, if more than 
                % one, then t0 < t. If t0 == t, then no inter-branch 
                % distance can be calculated and we skip this and increase 
                % number of the parents having a child.
                n_par = n_par + 1;
            end
            if(t0 < t)% More than 1 child branches were added.
                [data(t0:t,1),I] = sort(data(t0:t,1));% Sort the distances
                % If we assume that the parent branch is more or less
                % straight then we can take the subsequent values of
                % azimuth. We need this assumption since the reference
                % coord. system is always changing depending on the parent
                % orientation. Moreover, it is plausible if we use this
                % distribution only for child branches of the same year(in
                % case of coniferous trees), i.e. located very close to 
                % each other on the parent. The indicator of this is the 
                % bimodal character of the distance distribution (see 
                % scatter.idist).
                
                % account for the sorted information and the shift t0
                data(t0:t,3) = data(I+t0-1,3);
                % Azimuth difference
                dataAD = [dataAD; diff(data(t0:t,3))];
                % Distance difference
                dataD = [dataD; diff(data(t0:t,1))/BLen(j)];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% IDist-Thght distribution.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if(k == 1)% Similar to idist DF.
                    % Relative height along the trunk. data is sorted, 
                    % hence we can use index of thght to get the branch 
                    % number, taking into account that trunk's branch 
                    % number is 1.
                    thght = data(t0:t,1)/BLen(j);
                    dataH = [thght(1:end-1) dataD];% DF thght_idist.
                else
                    n_par = n_par + 1;% increase the raw number of parent
                    % BC is all children of the j-branch. The children are 
                    % of order k.
                    % length(BC) is the number of children for a specific
                    % j-branch. dataD is the length(BC)-1 size of the
                    % differences for those children.
                    
                    % iterate over the children, skip the last since dataD is of size length(BC)-1.
                    for i = 1:length(BC)-1
                        % Find a child in the map.
                        Ind = find(BC(i) == map(:,2));
                        % if not found, then skip the iteration
                        if(isempty(Ind))
                            %disp(['Warning: Ind is empty for parent ' num2str(j) ' and child ' num2str(BC(i)) '.']);
                            break;
                        end
                        % Ind must be unique!
                        if(length(Ind) ~= 1)
                            disp(Ind);
                            disp('Warning: length(Ind) is not 1.');
                        end
                        %disp(['t=' num2str(t) ':' num2str(length(dataD)) ':' num2str(t0) '+' num2str(i) '-' num2str(n_par)]);
                        % Finally, write down the data for the scatter.
                        % thght(... - 1) since order-1 branches start 
                        % with 2 (trunk = 1).
                        dataH = [dataH; thght(map(Ind,1)-1) dataD(t0+i-n_par)];
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            t0 = t + 1;% Increasing t0 in ANY case.
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Write down the OUTPUT scatter structure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Remove not plausible >1 values in the dataD
    I = dataD > 1;
    dataD(I) = 1;
    scatter(k).idist = dataD;
    scatter(k).idist_thght = dataH;
    if(k==1)
        scatter(k).mo_thght = dataO(1:tf,:);
    end
    scatter(k).az = (pi/180).*data(1:t,3);
    scatter(k).azd_idist = [dataD (pi/180).*dataAD];
    % Remove not plausible >1 values in the dataRD
    I = dataRD(:,2) > 1;
    dataRD(I,2) = 1;
    scatter(k).rchipar_rpar = dataRD(1:t,:);
end

% ***8*** Branch curvature: how a branch curves in the space is described
% by the angles gamma (in xy-plane) and zeta (in xz-plane) and the relative
% distance along the branch. The x-,y- and z-directions are determined by
% the 2 consequtive rotations of the global coordinate system. First, the
% XY-projection of the j-1 cyl form new x-axis, the Y-axis is rotated
% correspondingly to give y-axis => gamma is calculated via taking the
% projection of j cyl on xy-plane. Second, the x-axis is rotated to
% coincide with j-1 cyl and the Z-axis (not touched so far) is rotated
% correspondingly to give z-axis => zeta is calculated via taking the 
% projection onto the new xz-plane.
% CODE NAME: CURV (3D)

for k = 1:MAX_W% iterate until the max order is reached
    data = zeros(1,3);% initialize the helping variable data
    t = 0;% some counter
    for i = 1:N_BR% iterate over all branches
        if BOrd(i) == k% find the branch of order k
            C = CiB{i};% cylinders in the branch
            I = Added(C(1));% check if the 1st cylinder is Added one
            if I% if so...
                C = C(2:end);% start from the 2nd one
            end
            m = length(C);% number of cylinders in the branch
            L = sum(Len(C));% length of the branch
            j = 1;
            while (j <= m)
                j = j+1;% increase the counter
                if(j <= m)% to make sure we are still OK after increasing j
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % OLD CODE
%                     % Determine the reference
%                     D = cross(Axe(C(j-1),:),U);% find rotation axis
%                     if (norm(D) > 0)% if it exists...
%                         a = acos(Axe(C(j-1),:)*U);% rotation angle
%                         R = rotation_matrix(D,-a);% rotation matrix
%                         v = R*V;% rotate y
%                         w = R*W;% rotate z
%                     else% otherwise, leave them untouched
%                         v = V;
%                         w = W;
%                     end
                    % Compute the angles and distance
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % First rotation
                    newX = Axe(C(j-1),:);
                    newX(3) = 0.0;% projection on XY is new x-axis
                    newX = newX ./ sqrt(sum(newX.^2));% normalize
                    newY = cross(W,newX);% rotate (around Z) Y-axis => y-axis
                    % Z is not rotated, so xy=XY
                    % Orthogonal check
                    if((abs(newX*newY') > 1e-05) || (abs(newX*[0 0 1]') > 1e-05)...
                            || (abs(newY*[0 0 1]') > 1e-05))
                        disp('ERROR: Basis 1 is not orthogonal');
                    end
                    if(sum(newX.^2) > (1.0+1e-05) || sum(newY.^2) > (1.0+1e-05))
                        disp('ERROR: Basis 1 length is wrong');
                    end
                    targetAxe = Axe(C(j),:);
                    targetAxe(3) = 0.0;%Z-coord = 0 => XY=xy-projection
                    targetAxe = targetAxe ./ sqrt(sum(targetAxe.^2));%normalize
                    if(sum(targetAxe.^2) < 1e-05)
                        gamma = 0.0;
                    else
                        gamma = atan2(targetAxe * newY', targetAxe * newX');
                    end
                    
                    % Second rotation
                    newX = Axe(C(j-1),:);%new x-axis, not projection
                    newZ = cross(newX,newY);% rotate (around y) Z-axis => z-axis
                    % y is not rotated, although not original Y (newY
                    % modified above)
                    % Orthogonal check
                    if((abs(newX*newY') > 1e-05) || (abs(newX*newZ') > 1e-05)...
                            || (abs(newY*newZ') > 1e-05))
                        disp('ERROR: Basis 2 is not orthogonal');
                    end
                    if(sum(newX.^2) > (1.0+1e-05) || sum(newY.^2) > (1.0+1e-05) ...
                            || sum(newZ.^2) > (1.0+1e-05))
                        disp('ERROR: Basis 2 length is wrong');
                    end
                    targetAxe = Axe(C(j),:);
                    % The next two cross products give the projection of
                    % the target onto xz-plane.
                    targetAxe = cross(targetAxe,newY);
                    targetAxe = cross(newY,targetAxe);
                    targetAxe = targetAxe ./ sqrt(sum(targetAxe.^2));%normalize
                    if(sum(targetAxe.^2) < 1e-05)
                        zeta = 0.0;
                    else
                        zeta = atan2(targetAxe * newZ', targetAxe * newX');
                    end
                    t = t+1;% increase the counter
%                     gamma = Axe(C(j),:)*[Axe(C(j-1),:)' v];% xy-components
%                     zeta = Axe(C(j),:)*[Axe(C(j-1),:)' w];% xz-components
%                     data(t,1) = atan2(gamma(2),gamma(1));% gamma angle
%                     data(t,2) = atan2(zeta(2),zeta(1));% zeta angle
                    data(t,1) = gamma;
                    data(t,2) = zeta;
                    data(t,3) = sum(Len(C(1:j-1)))/L;% relative br. length
                end
            end
        end
    end
    data = data(1:t,:);% take only first t values
    scatter(k).curv = data(1:t,:);% copy to the output scatter struct
end

% ***9*** Tapering of the branches: radius of the cylinders vs. the 
% relative length along the branch. CODE NAME: TAP (2D)
% ***10*** Total length of the branch vs. its initial radius. CODE NAME:
% LTOT_RINI (2D)

% Constants
BR_FR = 0.1;% threshold fraction of the branch, see the code below

for k = 1:MAX_W
    p = 0;% some counter
    t = 0;% some counter
    data = zeros(1,2);% helping matrix
    data1 = zeros(1,2);% helping matrix
    for j = 1:N_BR% iterate over all branches
        if BOrd(j) == k% select only the particular order branch
            C = CiB{j};% find all cylinders in the branch
            L = cumsum(Len(C));% cumulative sum of the cyl. lengths
            %%% LTOT_RINI DF
            p = p+1;% increase counter
            data1(p,1:2) = [Rad(C(1)) L(end)];
            %%% TAP DF
            for i = 1:length(C)% iterate over all cylinders of the branch
                t = t+1;% increase the counter
                % Finally, get the tapering: radius of the cyl. vs. the
                % relative length of the branch at the place where the cyl
                % is situated (at the cyl's end).
%                 edge = L(i);% set the edge
%                 if(i > 1)% treat i == 1 case separately
%                     while((edge - L(i-1)) > BR_FR*L(end))
%                         edge = edge - (edge - L(i-1))/2;
%                     end
%                     incr = edge - L(i-1);% increment segment
%                     m = 1;
%                     while((L(i-1) + m*incr) < L(i))
%                         data(t,1:2) = [(L(i-1)+m*incr)./L(end) Rad(C(i))];
%                         if(data(t,1) < 0.1 && data(t,2) < 0.005)
%                             disp(data(t,:))
%                         end
%                         m = m+1;
%                         t = t+1;
%                     end
%                 elseif(i == 1)
%                     while(edge > BR_FR*L(end))
%                         edge = edge - edge/2;
%                     end
%                     incr = edge;
%                     m = 1;
%                     while((m*incr) < L(i))
%                         data(t,1:2) = [m*incr./L(end) Rad(C(i))];
%                         if(data(t,1) < 0.1 && data(t,2) < 0.005)
%                             disp(data(t,:))
%                         end
%                         m = m+1;
%                         t = t+1;
%                     end
%                 end
                data(t,1:2) = [L(i)./L(end) Rad(C(i))];
                data2(t,1:2) = [L(i) Rad(C(i))];
            end
        end
    end
    scatter(k).tap_norm = data(1:t,:);% copy to the output struct
    scatter(k).tap = data2(1:t,:);% non-normalized tapering
    scatter(k).ltot_rini = data1(1:p,:);% copy to the output struct
end

% ***11*** Branch angles: inclination angles at which branches emerge from
% the parent. CODE NAME: BRA (1D)

for k = 1:MAX_W
    % We consider only branches with parents
    P = find(BOrd == (k-1));% list of parent branches
    for ii = 1:length(P)
        C = BChi{P(ii)};% children branches of parent P(ii)
        % Somehow, during the process, children's order does not correspond
        % to the parent's order, that is not k+1 (if k is order of the par)
        C = C(BOrd(C) == k);% take only the right order branches
        scatter(k).bra = [scatter(k).bra; (pi/180).*BAng(C)];% pi/180 makes angles in radians
    end
end

% ***12*** Length of the CHIld vs. absolute Length Along the PARent branch.
% CODE NAME: LCHI_LAPAR (2D).

for k = 1:MAX_W
    P = find(BOrd == (k-1));% list of parent branches
    data = zeros(100,2);
    t = 0;
    for ii = 1:length(P)
        C = BChi{P(ii)};% children of the parent P(ii);
        % Somehow, during the process, children's order does not correspond
        % to the parent's order, that is not k+1 (if k is order of the par)
        C = C(BOrd(C) == k);% take only the right order branches
        CP = CiB{P(ii)};
        for jj = 1:length(C)
            t = t + 1;
            f = FCB(C(jj));
            pc = CPar(f);
            ind = find(CP == pc,1);
            [~,~,h] = distances_to_line(Sta(f,:),Axe(pc,:),Sta(pc,:));
            I = CP(1:ind-1);
            data(t,:) = [sum(Len(I))+h, sum(Len(CiB{C(jj)}))];
        end
    end
    scatter(k).lchi_lapar = data(1:t,:);
end

% ***13*** The combination scatter of LCHI_LAPAR and BRA to mimick a single
% scatter: to hide 2 in 1. LCHI and BRA as a function of LAPAR.
% CODE NAME: LCHI_BRA_LAPAR (3D)

for k = 1:MAX_W
    scatter(k).lchi_bra_lapar = [scatter(k).lchi_lapar scatter(k).bra];
end

% ----- END of SCATTERS CALCULATION ------

end


% ----- AUXILLARY FUNCTIONS -----

function s = curv(C,Axe,Len,W)
% Caclulate the curvature scatters
s = [];
m = length(C);% number of cylinders in the branch
L = sum(Len(C));% length of the branch
j = 1;
while (j < m)
    j = j+1;% increase the counter
    % First rotation
    newX = Axe(C(j-1),:);
    newX(3) = 0.0;% projection on XY is new x-axis
    newX = newX ./ sqrt(sum(newX.^2));% normalize
    newY = cross(W,newX);% new y-axis
    % Z is not rotated, so xy=XY
    % Orthogonal check
    if((abs(newX*newY') > 1e-05) || (abs(newX*W) > 1e-05)...
            || (abs(newY*W) > 1e-05))
        disp('ERROR: Basis 1 is not orthogonal');
    end
    if(sum(newX.^2) > (1.0+1e-05) || sum(newY.^2) > (1.0+1e-05))
        disp('ERROR: Basis 1 length is wrong');
    end
    targetAxe = Axe(C(j),:);
    targetAxe(3) = 0.0;%Z-coord = 0 => XY=xy-projection
    targetAxe = targetAxe ./ sqrt(sum(targetAxe.^2));%normalize
    if(sum(targetAxe.^2) < 1e-05)
        gamma = 0.0;
    else
        gamma = atan2(targetAxe * newY', targetAxe * newX');
    end
    
    % Second rotation
    newX = Axe(C(j-1),:);%new x-axis, not projection
    newZ = cross(newX,newY);% new z-axis
    % y is not rotated, although not original Y (newY inherited from above)
    % Orthogonal check
    if((abs(newX*newY') > 1e-05) || (abs(newX*newZ') > 1e-05)...
            || (abs(newY*newZ') > 1e-05))
        disp('ERROR: Basis 2 is not orthogonal');
    end
    if(sum(newX.^2) > (1.0+1e-05) || sum(newY.^2) > (1.0+1e-05) ...
            || sum(newZ.^2) > (1.0+1e-05))
        disp('ERROR: Basis 2 length is wrong');
    end
    targetAxe = Axe(C(j),:);
    % The next two cross products give the projection of
    % the target onto xz-plane.
    targetAxe = cross(targetAxe,newY);
    targetAxe = cross(newY,targetAxe);
    targetAxe = targetAxe ./ sqrt(sum(targetAxe.^2));%normalize
    if(sum(targetAxe.^2) < 1e-05)
        zeta = 0.0;
    else
        zeta = atan2(targetAxe * newZ', targetAxe * newX');
    end
    s = cat(1,s,[gamma zeta sum(Len(C(1:j-1)))/L]);
end
end

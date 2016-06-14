function [Branch, Segment, tree_out] = gen_scatter2(input_data_file)
% USAGE: [Branch, Segment, tree_out] = gen_scatter2(input_data_file)
%
% The function (version 2) computes the all neccessary scatters required by
% the subsequent analysis in the Bayes Forest tree generation as the source
% for the distributions. The prior knowledge of the certain characteristics
% (INPUT DATA) is needed. The characteristics come from the analysis of the
% measurements/simulations.
% IMPORTANT: as opposed to the 1st version, here the scatters are roughly
% divided into two groups: BRANCH and SEGMENT oriented data sets. The
% Branch data set contains information of the whole branches like total
% length, radis at the base, branching and azimuthal angles etc. The
% Segment data set contains information of the segments in any particular
% branch such as tapering of the branches and curvature in space. Both data
% sets are sorted by the order. Trunk (order 0) has only the Segment data 
% set.
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


%% ----- LOAD INPUT -----

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

% ----- END of LOAD INPUT -----

%% ----- INITIALIZATION OF THE DATA SETS -----

Branch = struct('info','data');
Branch.data = cell(1,max(BOrd)+1);
Branch.info = {'bra','az','ltot','rini','lapar'};

Segment = struct('info','data');
Segment.data = cell(1,max(BOrd)+1);
Segment.info = {'rad','len','gamma','zeta'};

% Some constants
%N_CYL = length(Rad);% number of cylinders
N_BR = length(BOrd);% number of branches
U = [1 0 0]';% reference branch direction
V = [0 1 0]';% reference vector for gamma (angle in xy-plane)
W = [0 0 1]';% reference vector for zeta (angle in xz-plane)

% ----- END of INITIALIZATION OF THE DATA SETS -----

%% ----- TREE-OBJECT -----
% Simple code to covert to the fields of tree-object. Do not care about the
% topology information, just for visualization purposes.
tree_out = tree();
tree_out.number_of_branches = length(Rad);
tree_out.radius = Rad;
tree_out.length = Len;
tree_out.start_point = Sta;
tree_out.axis = Axe;
tree_out.end_point = Sta + bsxfun(@times,Len,Axe);

% ----- END of TREE-OBJECT -----

%% ----- SCATTERS CALCULATION ------
n_unprocessed = 0;
unprocessed = [];
for order = 0:max(BOrd)
    Br = find(BOrd == order);% all branches of the order = order.
    for ii = 1:length(Br)% iterate over branches
        %% BRANCH data sets
        if(order > 0)% Exclude trunk from the Branch data set
            % Branch must have parent determined, if this was not
            % determined we omit the branch from the analysis
            if(CPar(FCB(Br(ii))) ~= 0)% branch must have a parent
                bra = BAng(Br(ii));
                az = compute_az(Br(ii));
                ltot = BLen(Br(ii));
                rini = Rad(FCB(Br(ii)));
                lapar = compute_lapar(Br(ii),order);
                Branch.data{order+1} = cat(1,Branch.data{order+1},[bra az ltot rini lapar]);
            else
                n_unprocessed = n_unprocessed + 1;
                unprocessed = cat(2,unprocessed,Br(ii));
                %disp(['Branch ' num2str(Br(ii)) ' has no parent']);
                %disp(['BAng = ' num2str(BAng(Br(ii))) '.']);
            end
        end
        %% SEGMENT data sets
        if(order > 0 && CPar(FCB(Br(ii))) == 0)
            % REMOVE this if you want to account for the branches, not
            % processed for the Branch data-sets, in the Segment data-sets
            continue;
        end
        len = 0.0;
        for jj = 1:length(CiB{Br(ii)})% segments in the branch
            seg = CiB{Br(ii)}(jj);
            if(jj == 1) 
                gamma = 0; 
                zeta = 0; 
            else
                k = CiB{Br(ii)}(jj-1); 
                [gamma, zeta] = curvature(seg,k);
            end
            Segment.data{order+1} = cat(1,Segment.data{order+1},[Rad(seg) len gamma zeta]);
            len = len + Len(seg);
        end
    end
end
if(n_unprocessed)
    disp(['No unprocessed branches: ' num2str(n_unprocessed) ' (' ...
        num2str(n_unprocessed/length(BOrd)) '%).']);
    disp(['Total no branches: ' num2str(length(BOrd)) '.']);
    order_unprocessed = unique(BOrd(unprocessed));
    disp('Order of the unprocessed (order/no branches): ');
    for ii = 1:length(order_unprocessed)
        fprintf('%d/%d ',order_unprocessed(ii),length(find(BOrd(unprocessed) == order_unprocessed(ii))));
    end
    fprintf('\n');
end


%% ----- INLINE FUNCTIONS -----
    function az = compute_az(br)
        % Calculate azimuth angles. BR - branch of the given order.
        rot_axis = cross(Axe(CPar(FCB(br)),:),U);
        first = FCB(br);
        ang = acos(Axe(CPar(first),:)*U);
        R = rotation_matrix(rot_axis,-ang);% rotation matrix
        if (norm(R) > 0)
            v = R*V;% rotate V
            w = R*W;% rotate W
        else
            v = V;
            w = W;
        end
        az = atan2d(Axe(first,:)*w,Axe(first,:)*v);% in degrees!
    end
    function lapar = compute_lapar(br,order)
        % Calculate Distance from the parent's base at which the branch BR
        % emanates
        % 1. Find the parent
        parSeg = CPar(FCB(br));% parent segment of the parent branch
        potParBrs = find(BOrd == order-1);% all potential parent branches
        for m = 1:length(potParBrs)
            if(~isempty(find(CiB{potParBrs(m)} == parSeg,1)))
                parBr = potParBrs(m);
                break;
            end
        end
        % 2. Find the length along the parent
        [~,~,h] = distances_to_line(Sta(FCB(br),:),Axe(parSeg,:),Sta(parSeg,:));
        I = find(CiB{parBr} == parSeg);
        lapar = sum(Len(CiB{parBr}(1:I-1)))+h;
    end
    function [gamma, zeta] = curvature(next,prev)
        % Calculate curvature angles gamma and zeta based on the next and
        % previous segments. Two rotation method.
        % First rotation
        newX = Axe(prev,:);
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
        targetAxe = Axe(next,:);
        targetAxe(3) = 0.0;%Z-coord = 0 => XY=xy-projection
        targetAxe = targetAxe ./ sqrt(sum(targetAxe.^2));%normalize
        if(sum(targetAxe.^2) < 1e-05)
            gamma = 0.0;
        else
            gamma = atan2d(targetAxe * newY', targetAxe * newX');% in deg!
        end
        % Second rotation
        newX = Axe(prev,:);%new x-axis, not projection
        newZ = cross(newX,newY);% rotate (around y) Z-axis => z-axis
        % y is not rotated, although not original Y (newY modified above)
        % Orthogonal check
        if((abs(newX*newY') > 1e-05) || (abs(newX*newZ') > 1e-05)...
                || (abs(newY*newZ') > 1e-05))
            disp('ERROR: Basis 2 is not orthogonal');
        end
        if(sum(newX.^2) > (1.0+1e-05) || sum(newY.^2) > (1.0+1e-05) ...
                || sum(newZ.^2) > (1.0+1e-05))
            disp('ERROR: Basis 2 length is wrong');
        end
        targetAxe = Axe(next,:);
        % The next two cross products give the projection of
        % the target onto xz-plane.
        targetAxe = cross(targetAxe,newY);
        targetAxe = cross(newY,targetAxe);
        targetAxe = targetAxe ./ sqrt(sum(targetAxe.^2));%normalize
        if(sum(targetAxe.^2) < 1e-05)
            zeta = 0.0;
        else
            zeta = atan2d(targetAxe * newZ', targetAxe * newX');% in deg!
        end
    end
end


% for k = 1:MAX_W
%     m = nnz(BOrd == k);
%     % some data collection matrices
%     data = zeros(m,3);
%     %dataR = zeros(m,2);
%     dataRD = zeros(m,2);
%     dataD = [];
%     dataH = [];
%     dataAD = [];
%     % some helping counters
%     t = 0;
%     t0 = 1;
%     n_par = 0;% number of the current order parents, which have a child
%     for j = 1:N_BR% over all branches
%         if BOrd(j) == k-1% A possible parent branch j for the current order k
%             C = CiB{j};% Cylinders in the branch j, cylinder index
%             Ind = 1:1:length(C);% Cylinder index list
%             BC = BChi{j};% All child branches of the branch j
%             n = length(BC);% Number of the child branches of the branch j
%             for i = 1:n% iterate over the children of order k
%                 f = FCB(BC(i));% First cylinder in the child i
%                 % Find emergence point (cylinder) of child i on parent j
%                 I = CPar(f) == C;
%                 CP = C(I);% parent cylinder index (CP)
%                 I = Ind(I);% row number of the parent cylinder
%                 if isempty(CP)% empty CP: the parent is Added cylinder
%                     if(~Added(CPar(f)))% if it is NOT Added, it is STRANGE!
%                         disp('Strange!');
%                     end
%                     CP = CPar(CPar(f));% find parent of parent j of child i
%                     I = CP == C;% parent cylinder index
%                     I = Ind(I);% row number of the parent cylinder
%                 end
%                 % rotation vector from X-axis to parent axis (new x-axis)
%                 D = cross(Axe(CP,:),U);
%                 if norm(D) > 0
%                     % angle between x-axis and new x-axis, y- and z-axis 
%                     % are rotated for this angle as well
%                     a = acos(Axe(CP,:)*U);
%                     R = rotation_matrix(D,-a);% rotation matrix
%                     v = R*V;% rotate V
%                     w = R*W;% rotate W
%                 else% no rotation is needed
%                     v = V;
%                     w = W;
%                 end
%                 % Look into distances_to_line code. h is the segment of
%                 % Axe(CP,:) line between the Sta(CP,:) point and the point 
%                 % where the perpendicular from Sta(f,:) falls.
%                 % Normalization of the Axe(CP,:) takes place in the
%                 % function distances_to_line.
%                 [~,~,h] = distances_to_line(Sta(f,:),Axe(CP,:),Sta(CP,:));
%                 % y and z components, giving the direction in yz-plane
%                 a = Axe(f,:)*[v w];
%                 t = t+1;% increase counter
%                 % length along the parent, starting from the 1st cylinder
%                 data(t,1) = sum(Len(C(1:I-1)))+h;
%                 % normalize the length along the parent
%                 data(t,2) = data(t,1)/sum(Len(C));
%                 % find the azimuth, reference is v (new y-axis)
%                 data(t,3) = 180/pi*atan2(a(2),a(1));
%                 % initial radius of child i vs. distance along parent j
%                 %dataR(t,:) = [data(t,2) Rad(f)];
%                 % ratio between initial radius of child i and its parent
%                 % cylinder radius vs. radius of the parent cylinder
%                 dataRD(t,:) = [Rad(CP) Rad(f)/Rad(CP)];
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 %%% Max order distribution.
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 % Need to check only the order-1 child branches
%                 if(k == 1)
%                     p = 1;% Max order p is 1, initially.
%                     bchild = BC(i);% get the current child branch
%                     while(p)% Iterate
%                         p_old = p;% Save the old value of p
%                         % Iterate over the target child branches: iterative
%                         % process starting from the order-1 bchild and then
%                         % progressively moving towards its descendants.
%                         for q = 1:length(bchild)
%                             % if there is a child of the higher order
%                             % then increase p and break the for loop moving
%                             % to the next order.
%                             if(~isempty(BChi{bchild(q)}))
%                                 p = p+1;
%                                 break;
%                             end
%                         end
%                         % if the p was not increased, break while loop
%                         if(p == p_old), break; end
%                         % create the new list of the target child branches.
%                         bchild = BChi{bchild};
%                     end
%                     % Finally, write down data to the scatter.
%                     dataO(t,:) = [data(t,2) p];
%                     tf = t;% mark the final t, for the extraction later.
%                 end
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%
%             end
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             %%% Inter-branch difference DF's
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             % If there was no child branch, t was not increased and t0 was.
%             % So we need to subtract 1 from t0
%             if(t0 == t+1)
%                 t0 = t0 - 1;
%             elseif(t0 == t)
%                 % If ONE child branch was added then t0 == t, if more than 
%                 % one, then t0 < t. If t0 == t, then no inter-branch 
%                 % distance can be calculated and we skip this and increase 
%                 % number of the parents having a child.
%                 n_par = n_par + 1;
%             end
%             if(t0 < t)% More than 1 child branches were added.
%                 [data(t0:t,1),I] = sort(data(t0:t,1));% Sort the distances
%                 % If we assume that the parent branch is more or less
%                 % straight then we can take the subsequent values of
%                 % azimuth. We need this assumption since the reference
%                 % coord. system is always changing depending on the parent
%                 % orientation. Moreover, it is plausible if we use this
%                 % distribution only for child branches of the same year(in
%                 % case of coniferous trees), i.e. located very close to 
%                 % each other on the parent. The indicator of this is the 
%                 % bimodal character of the distance distribution (see 
%                 % scatter.idist).
%                 
%                 % account for the sorted information and the shift t0
%                 data(t0:t,3) = data(I+t0-1,3);
%                 % Azimuth difference
%                 dataAD = [dataAD; diff(data(t0:t,3))];
%                 % Distance difference
%                 dataD = [dataD; diff(data(t0:t,1))/BLen(j)];
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 %%% IDist-Thght distribution.
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 if(k == 1)% Similar to idist DF.
%                     % Relative height along the trunk. data is sorted, 
%                     % hence we can use index of thght to get the branch 
%                     % number, taking into account that trunk's branch 
%                     % number is 1.
%                     thght = data(t0:t,1)/BLen(j);
%                     dataH = [thght(1:end-1) dataD];% DF thght_idist.
%                 else
%                     n_par = n_par + 1;% increase the raw number of parent
%                     % BC is all children of the j-branch. The children are 
%                     % of order k.
%                     % length(BC) is the number of children for a specific
%                     % j-branch. dataD is the length(BC)-1 size of the
%                     % differences for those children.
%                     
%                     % iterate over the children, skip the last since dataD is of size length(BC)-1.
%                     for i = 1:length(BC)-1
%                         % Find a child in the map.
%                         Ind = find(BC(i) == map(:,2));
%                         % if not found, then skip the iteration
%                         if(isempty(Ind))
%                             %disp(['Warning: Ind is empty for parent ' num2str(j) ' and child ' num2str(BC(i)) '.']);
%                             break;
%                         end
%                         % Ind must be unique!
%                         if(length(Ind) ~= 1)
%                             disp(Ind);
%                             disp('Warning: length(Ind) is not 1.');
%                         end
%                         %disp(['t=' num2str(t) ':' num2str(length(dataD)) ':' num2str(t0) '+' num2str(i) '-' num2str(n_par)]);
%                         % Finally, write down the data for the scatter.
%                         % thght(... - 1) since order-1 branches start 
%                         % with 2 (trunk = 1).
%                         dataH = [dataH; thght(map(Ind,1)-1) dataD(t0+i-n_par)];
%                     end
%                 end
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             end
%             t0 = t + 1;% Increasing t0 in ANY case.
%         end
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%% Write down the OUTPUT scatter structure
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Remove not plausible >1 values in the dataD
%     I = dataD > 1;
%     dataD(I) = 1;
%     scatter(k).idist = dataD;
%     scatter(k).idist_thght = dataH;
%     if(k==1)
%         scatter(k).mo_thght = dataO(1:tf,:);
%     end
%     scatter(k).az = (pi/180).*data(1:t,3);
%     scatter(k).azd_idist = [dataD (pi/180).*dataAD];
%     % Remove not plausible >1 values in the dataRD
%     I = dataRD(:,2) > 1;
%     dataRD(I,2) = 1;
%     scatter(k).rchipar_rpar = dataRD(1:t,:);
% end
% 
% % ***8*** Branch curvature: how a branch curves in the space is described
% % by the angles gamma (in xy-plane) and zeta (in xz-plane) and the relative
% % distance along the branch. The x-,y- and z-directions are determined by
% % the 2 consequtive rotations of the global coordinate system. First, the
% % XY-projection of the j-1 cyl form new x-axis, the Y-axis is rotated
% % correspondingly to give y-axis => gamma is calculated via taking the
% % projection of j cyl on xy-plane. Second, the x-axis is rotated to
% % coincide with j-1 cyl and the Z-axis (not touched so far) is rotated
% % correspondingly to give z-axis => zeta is calculated via taking the 
% % projection onto the new xz-plane.
% % CODE NAME: CURV (3D)
% 
% for k = 1:MAX_W% iterate until the max order is reached
%     data = zeros(1,3);% initialize the helping variable data
%     t = 0;% some counter
%     for i = 1:N_BR% iterate over all branches
%         if BOrd(i) == k% find the branch of order k
%             C = CiB{i};% cylinders in the branch
%             I = Added(C(1));% check if the 1st cylinder is Added one
%             if I% if so...
%                 C = C(2:end);% start from the 2nd one
%             end
%             m = length(C);% number of cylinders in the branch
%             L = sum(Len(C));% length of the branch
%             j = 1;
%             while (j <= m)
%                 j = j+1;% increase the counter
%                 if(j <= m)% to make sure we are still OK after increasing j
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     % OLD CODE
% %                     % Determine the reference
% %                     D = cross(Axe(C(j-1),:),U);% find rotation axis
% %                     if (norm(D) > 0)% if it exists...
% %                         a = acos(Axe(C(j-1),:)*U);% rotation angle
% %                         R = rotation_matrix(D,-a);% rotation matrix
% %                         v = R*V;% rotate y
% %                         w = R*W;% rotate z
% %                     else% otherwise, leave them untouched
% %                         v = V;
% %                         w = W;
% %                     end
%                     % Compute the angles and distance
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     % First rotation
%                     newX = Axe(C(j-1),:);
%                     newX(3) = 0.0;% projection on XY is new x-axis
%                     newX = newX ./ sqrt(sum(newX.^2));% normalize
%                     newY = cross(W,newX);% rotate (around Z) Y-axis => y-axis
%                     % Z is not rotated, so xy=XY
%                     % Orthogonal check
%                     if((abs(newX*newY') > 1e-05) || (abs(newX*[0 0 1]') > 1e-05)...
%                             || (abs(newY*[0 0 1]') > 1e-05))
%                         disp('ERROR: Basis 1 is not orthogonal');
%                     end
%                     if(sum(newX.^2) > (1.0+1e-05) || sum(newY.^2) > (1.0+1e-05))
%                         disp('ERROR: Basis 1 length is wrong');
%                     end
%                     targetAxe = Axe(C(j),:);
%                     targetAxe(3) = 0.0;%Z-coord = 0 => XY=xy-projection
%                     targetAxe = targetAxe ./ sqrt(sum(targetAxe.^2));%normalize
%                     if(sum(targetAxe.^2) < 1e-05)
%                         gamma = 0.0;
%                     else
%                         gamma = atan2(targetAxe * newY', targetAxe * newX');
%                     end
%                     
%                     % Second rotation
%                     newX = Axe(C(j-1),:);%new x-axis, not projection
%                     newZ = cross(newX,newY);% rotate (around y) Z-axis => z-axis
%                     % y is not rotated, although not original Y (newY
%                     % modified above)
%                     % Orthogonal check
%                     if((abs(newX*newY') > 1e-05) || (abs(newX*newZ') > 1e-05)...
%                             || (abs(newY*newZ') > 1e-05))
%                         disp('ERROR: Basis 2 is not orthogonal');
%                     end
%                     if(sum(newX.^2) > (1.0+1e-05) || sum(newY.^2) > (1.0+1e-05) ...
%                             || sum(newZ.^2) > (1.0+1e-05))
%                         disp('ERROR: Basis 2 length is wrong');
%                     end
%                     targetAxe = Axe(C(j),:);
%                     % The next two cross products give the projection of
%                     % the target onto xz-plane.
%                     targetAxe = cross(targetAxe,newY);
%                     targetAxe = cross(newY,targetAxe);
%                     targetAxe = targetAxe ./ sqrt(sum(targetAxe.^2));%normalize
%                     if(sum(targetAxe.^2) < 1e-05)
%                         zeta = 0.0;
%                     else
%                         zeta = atan2(targetAxe * newZ', targetAxe * newX');
%                     end
%                     t = t+1;% increase the counter
% %                     gamma = Axe(C(j),:)*[Axe(C(j-1),:)' v];% xy-components
% %                     zeta = Axe(C(j),:)*[Axe(C(j-1),:)' w];% xz-components
% %                     data(t,1) = atan2(gamma(2),gamma(1));% gamma angle
% %                     data(t,2) = atan2(zeta(2),zeta(1));% zeta angle
%                     data(t,1) = gamma;
%                     data(t,2) = zeta;
%                     data(t,3) = sum(Len(C(1:j-1)))/L;% relative br. length
%                 end
%             end
%         end
%     end
%     data = data(1:t,:);% take only first t values
%     scatter(k).curv = data(1:t,:);% copy to the output scatter struct
% end
% 

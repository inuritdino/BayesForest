function obj = delete_branch(obj,cyls)
% USAGE: OBJ = DELETE_BRANCH(OBJ,CYLS)
% Delete ANY subset of cyl's CYLS from the tree object. The function is
% divided into two distinct pieces: 1. to delete corresponding properties
% of the cyl's to remove and 2. adjust the row position to remove
% zero-value gaps in the object.
% Parent-extension-children information is preserved within the cyl's that
% are left after removal.

%% Deleting properties
N = obj.number_of_branches;

if(any(cyls > N))
    fprintf('Warning: some cyl''s out of the range.');
    cyls = cyls(cyls <= N);
end
if(isempty(cyls))
    fprintf('Warning: no cyl''s to be removed. Exit.');
    return;
end
% Make cyls unique, in case use supplies double values
cyls = unique(cyls);

% Set properties to zero.
remove_prop(cyls);

%% Shifting the rows of the object
% Note N was not modified and contains the old number of branches=cyl's

gap = [];% the newly created gaps during the shifting process
ii = 1;% the indicator of the empty positions=gaps
while (ii <= N)
    if(ismember(ii,cyls) || ismember(ii,gap))% We are at the empty position
        % Iterate from current gap to the first meaningful position
        jj = ii + 1;% JJ will hold the meaningful position
        while(ismember(jj,cyls) || ismember(jj,gap))% Iterate over all removed position=gaps
            jj = jj + 1;
        end
        % JJ contains first meaningful position OR outside of the range
        if(jj <= N)
            % Move the properties
            move_prop(ii,jj);
            % We moved from JJ to II, JJ creates a gap
            gap = cat(2,gap,jj);
        end
    end
    ii = ii + 1;
end

%% Make the inline (embedded) functions for faster operation when OBJ is BIG
% Look at the end of the code for the isolated functions.
    function remove_prop(cyls)
        % Set properties of the object to zero
        
        obj.radius(cyls) = 0;
        obj.length(cyls) = 0;
        obj.start_point(cyls,:) = 0;
        obj.end_point(cyls,:) = 0;
        obj.axis(cyls,:) = 0;
        obj.is_trunk(cyls) = 0;
        
        for kk = cyls
            % II has parent
            if(obj.parent(kk))
                % Remove extension information of the parent
                if(obj.extension(obj.parent(kk)) == kk)% II is the extension
                    obj.extension(obj.parent(kk)) = 0;
                else% II is the child, remove it from the children list of the parent
                    obj.children{obj.parent(kk)} = setdiff(obj.children{obj.parent(kk)},kk);
                end
            end
            % Delete parent
            obj.parent(kk) = 0;
            % II is parent
            if(obj.extension(kk))% parent of extension
                obj.parent(obj.extension(kk)) = 0;% no parent anymore for the extension
                obj.extension(kk) = 0;% remove the extension itself
            end
            if(~isempty(obj.children{kk}))% parent of children
                obj.parent(obj.children{kk}) = 0;%no parent for the children
                obj.children{kk} = [];% remove the children themselves
            end
        end
        
        % Adjust number of cyl's information. The very last operation.
        obj.number_of_branches = obj.number_of_branches - length(cyls);
        
    end
    function move_prop(dst,src)
        % Move properties from DST row in OBJ to its SRC row.
        
        if(numel(dst) ~= 1 || numel(src) ~= 1)
            error('Error:MoveProp: ambiguity in the command.');
        end
        % Copy SRC position to the DST gap and remove the corresponding position in
        % SRC
        obj.radius(dst) = obj.radius(src);
        obj.radius(src) = 0;
        obj.length(dst) = obj.length(src);
        obj.length(src) = 0;
        obj.start_point(dst,:) = obj.start_point(src,:);
        obj.start_point(src,:) = 0;
        obj.end_point(dst,:) = obj.end_point(src,:);
        obj.end_point(src,:) = 0;
        obj.axis(dst,:) = obj.axis(src,:);
        obj.axis(src,:) = 0;
        obj.is_trunk(dst) = obj.is_trunk(src);
        obj.is_trunk(src) = 0;
        obj.parent(dst) = obj.parent(src);
        obj.parent(src) = 0;
        obj.extension(dst) = obj.extension(src);
        obj.extension(src) = 0;
        obj.children{dst} = obj.children{src};
        obj.children{src} = [];
        % SRC is child or extension, change its parent content
        if(obj.parent(dst))% There is parent for the cyl
            if(obj.extension(obj.parent(dst)) == src)% dst is the extension of the parent
                obj.extension(obj.parent(dst)) = dst;% change to dst (prev value was src)
            else% dst is the child of the parent
                % First identify the position of SRC in parent's child array and change it to DST
                obj.children{obj.parent(dst)}(obj.children{obj.parent(dst)} == src) = dst;
                % Sort the children (may be not necessary)
                obj.children{obj.parent(dst)} = sort(obj.children{obj.parent(dst)});
            end
        end
        % SRC is parent, change its offspring content
        if(obj.extension(dst))
            obj.parent(obj.extension(dst)) = dst;% it was SRC before, change to DST
        end
        if(~isempty(obj.children{dst}))
            obj.parent(obj.children{dst}) = dst;% it was SRC before, change to DST
        end
        
    end
end
%% Remove the properties of one or several rows in OBJ
% function obj = remove_prop(obj,cyls)
% % Set properties of the object to zero
% 
% obj.radius(cyls) = 0;
% obj.length(cyls) = 0;
% obj.start_point(cyls,:) = 0;
% obj.end_point(cyls,:) = 0;
% obj.axis(cyls,:) = 0;
% obj.is_trunk(cyls) = 0;
% 
% for kk = cyls
%     % II has parent
%     if(obj.parent(kk))
%         % Remove extension information of the parent
%         if(obj.extension(obj.parent(kk)) == kk)% II is the extension
%             obj.extension(obj.parent(kk)) = 0;
%         else% II is the child, remove it from the children list of the parent
%             obj.children{obj.parent(kk)} = setdiff(obj.children{obj.parent(kk)},kk);
%         end
%     end
%     % Delete parent
%     obj.parent(kk) = 0;
%     % II is parent
%     if(obj.extension(kk))% parent of extension
%         obj.parent(obj.extension(kk)) = 0;% no parent anymore for the extension
%         obj.extension(kk) = 0;% remove the extension itself
%     end
%     if(~isempty(obj.children{kk}))% parent of children
%         obj.parent(obj.children{kk}) = 0;%no parent for the children
%         obj.children{kk} = [];% remove the children themselves
%     end
% end
% 
% % Adjust number of cyl's information. The very last operation.
% obj.number_of_branches = obj.number_of_branches - length(cyls);
% 
% end
%% Move one row in OBJ to other place
% function obj = move_prop(obj,dst,src)
% % Move properties from DST row in OBJ to its SRC row.
% 
% if(numel(dst) ~= 1 || numel(src) ~= 1)
%     error('Error:MoveProp: ambiguity in the command.');
% end
% % Copy SRC position to the DST gap and remove the corresponding position in
% % SRC
% obj.radius(dst) = obj.radius(src);
% obj.radius(src) = 0;
% obj.length(dst) = obj.length(src);
% obj.length(src) = 0;
% obj.start_point(dst,:) = obj.start_point(src,:);
% obj.start_point(src,:) = 0;
% obj.end_point(dst,:) = obj.end_point(src,:);
% obj.end_point(src,:) = 0;
% obj.axis(dst,:) = obj.axis(src,:);
% obj.axis(src,:) = 0;
% obj.is_trunk(dst) = obj.is_trunk(src);
% obj.is_trunk(src) = 0;
% obj.parent(dst) = obj.parent(src);
% obj.parent(src) = 0;
% obj.extension(dst) = obj.extension(src);
% obj.extension(src) = 0;
% obj.children{dst} = obj.children{src};
% obj.children{src} = [];
% % SRC is child or extension, change its parent content
% if(obj.parent(dst))% There is parent for the cyl
%     if(obj.extension(obj.parent(dst)) == src)% dst is the extension of the parent
%         obj.extension(obj.parent(dst)) = dst;% change to dst (prev value was src)
%     else% dst is the child of the parent
%         % First identify the position of SRC in parent's child array and change it to DST
%         obj.children{obj.parent(dst)}(obj.children{obj.parent(dst)} == src) = dst;
%         % Sort the children (may be not necessary)
%         obj.children{obj.parent(dst)} = sort(obj.children{obj.parent(dst)});
%     end
% end
% % SRC is parent, change its offspring content
% if(obj.extension(dst))
%     obj.parent(obj.extension(dst)) = dst;% it was SRC before, change to DST
% end
% if(~isempty(obj.children{dst}))
%     obj.parent(obj.children{dst}) = dst;% it was SRC before, change to DST
% end
% 
% end

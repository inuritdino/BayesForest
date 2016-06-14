function br = get_branches(tree)
% Extract the branching topology from a tree object TREE
% USE:
%     BR = GET_BRANCHES(TREE)
%
% TREE is a tree class (see elsewhere)
% BR is an array of branching structs with fields .order, .parent_num, and .cyl_ind.
% 	BR(I).ORDER is the branch I order (starting with 0 for trunk)
%	BR(I).PARENT_NUM is the I's parent branch number in the same struct BR
% 	BR(I).CYL_IND is the index of all cylinders belonging to the branch I.

% Br init
br = struct('order',[],'parent_num',[],'cyl_ind',[]);
% Start with the first cyl, which is not supposed to have parent
br_offset = 1;
% Br number
brN = 0;
% Order
order = 0;
while(br_offset)
    chi = [];
    for ii=1:length(br_offset)
        brN = brN + 1;
        br(brN).order = order;% Order
        tmp = find_branch(br,[],[],tree.parent(br_offset(ii)));% Parent
        if(isempty(tmp))
            br(brN).parent_num = 0;
        else
            br(brN).parent_num = tmp;
        end
        % Iterate to get cyl's of the branch
        tmp = br_offset(ii);
        while(tmp)
            br(brN).cyl_ind = cat(2,br(brN).cyl_ind,tmp);
            tmp = tree.extension(tmp);
        end
        %br(brN).start = tree.start_point(br_offset(ii),:);% Start position
        % Get aggregatedly children => new branch offsets
        for jj=1:length(br(brN).cyl_ind)
            chi = cat(1,chi,tree.children{br(brN).cyl_ind(jj)});
        end
    end
    br_offset = chi;
    order = order + 1;% increase order for children.
end


end

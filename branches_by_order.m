function [newtr, num_bra_by_order, num_cyl_by_order] = branches_by_order(tr, Order)
% Very rude (Matlab-ically) algorithm to extract branches up to the given
% topological order, calculated empirically based on the thickness of the
% branches (thickest is extension, rest are children).
% USAGE:
%     NEWTR = BRANCHES_BY_ORDER(TR, ORDER)
% TR: tree object
% ORDER: max order of the extracted branches, i.e. <=ORDER (default 1)
% NEWTR: new tree with orders up to the given order ORDER
%
% NOTE:
% Some discrepancy was found between similar algorithms in different
% languages/implementations: sometimes the algorithm chooses different
% branch pathways due to the equal possibilities between options (some
% internal double precision numbers representation affect this process).
% Obviously, this happens for higher order (>1) branches where the
% probability of similar width segments is high. This affects the number of
% branches for the given (high) order, which may different for different
% implementations.
%
% See also : tree

if(nargin < 2)
    Order = 1;
end

newtr = tree();

bases = 1;% first cyl is initiator of a branch
orders = 0;
num_bra_by_order = zeros(1,Order+1);
num_cyl_by_order = zeros(1,Order+1);
while (~isempty(bases))
    cc = bases(1);
    curr_order = orders(1);
    if( curr_order > Order )
        break;
    end
    tmp_br = cc;
    rad = -1;
    k = -1;
    while( tr.extension(cc) || ~isempty(tr.children{cc}) )
       % k = 0 => extension, 1,2,3... =>  children
       % Find the thickest branch
       if ( tr.extension(cc) )
           rad = tr.radius(tr.extension(cc));
           k = 0;
       end
       for ii = 1:length(tr.children{cc})
           if( rad < tr.radius(tr.children{cc}(ii)) )
               rad = tr.radius(tr.children{cc}(ii));
               k = ii;
           end
       end
       % Form new topology
       if( k == 0 )
           % nominal children are real children
           bases = cat(2,bases,tr.children{cc});
           orders = cat(2,orders,repmat(curr_order+1,1,length(tr.children{cc})));
           % nominal extension is real extension
           cc = tr.extension(cc);
       elseif ( k > 0 )
           if(tr.extension(cc))
               bases = cat(2,bases,tr.extension(cc));
               orders = cat(2,orders,curr_order+1);
           end
           bases = cat(2,bases,[tr.children{cc}(1:k-1) tr.children{cc}(k+1:end)]);
           orders = cat(2,orders,repmat(curr_order+1,1,length(tr.children{cc})-1));
           cc = tr.children{cc}(k);
       else
           break;
       end
       k = -1;
       rad = -1;
       tmp_br = cat(2,tmp_br,cc);
    end
    bases = bases(2:end);
    orders = orders(2:end);
    num_bra_by_order(curr_order+1) = num_bra_by_order(curr_order+1) + 1;
    num_cyl_by_order(curr_order+1) = num_cyl_by_order(curr_order+1) + length(tmp_br);
    % Copy the extracted branch
    nbr0 = newtr.number_of_branches;
    newtr.number_of_branches = newtr.number_of_branches + length(tmp_br);
    nbr = newtr.number_of_branches;
    newtr.radius(nbr0+1:nbr) = tr.radius(tmp_br);
    newtr.length(nbr0+1:nbr) = tr.length(tmp_br);
    newtr.start_point(nbr0+1:nbr,:) = tr.start_point(tmp_br,:);
    newtr.end_point(nbr0+1:nbr,:) = tr.end_point(tmp_br,:);
    newtr.axis(nbr0+1:nbr,:) = tr.axis(tmp_br,:);
end

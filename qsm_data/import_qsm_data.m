function tr = import_qsm_data(branch_data_table,cyl_data_table,out_fn)
% Import Quantitative Structural Model (QSM) data.
% USE:
%     TREE = IMPORT_QSM_DATA(BRANCH_DATA_TABLE,CYL_DATA_TABLE,OUT_FN)
%
%
% INPUT:
%
% BRANCH_DATA_FN - branches table-data file name or data matrix
% CYL_DATA_FN - cylinder table-data file name or data matrix
% OUT_FN - (optional) output .mat file name (default is 'tree_out.mat')
%
%
% OUTPUT:
%
% TREE - the tree object tree, which can be, for example, drawn with
% tree.draw command.

if(nargin < 3)
    out_fn = 'qsm_out.mat';
end

if(isa(branch_data_table,'char'))
    A = importdata(branch_data_table);
else
    A = branch_data_table;
end
BOrd = A(:,1);
BPar = A(:,2);
BVol = A(:,3);
BLen = A(:,4);
BAng = A(:,5);
clear A;
if(isa(cyl_data_table,'char'))
    A = importdata(cyl_data_table);
else
    A = cyl_data_table;
end
Rad = A(:,1);
Len = A(:,2);
Sta = [A(:,3) A(:,4) A(:,5)];
Axe = [A(:,6) A(:,7) A(:,8)];
CPar = A(:,9);
CExt = A(:,10);
BoC = A(:,11); %A(:,12) A(:,13)];
n = size(A,2);
if(n == 14)
    Added = A(:,14);
elseif(n == 12)
    Added = A(:,12);
else
    fprintf('Number of columns in data odd\n');
end
clear A branch_data_table cyl_data_table n;

save(out_fn);

% Make tree object out of it
tr = tree();
tr.radius = Rad;
tr.length = Len;
tr.start_point = Sta;
tr.axis = Axe;
tr.end_point = tr.start_point + bsxfun(@times,tr.length,tr.axis);
tr.parent = CPar;
tr.extension = CExt;
tr.number_of_branches = size(Rad,1);

end
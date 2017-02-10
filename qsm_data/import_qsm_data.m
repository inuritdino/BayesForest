function tr = import_qsm_data(branch_data_table,cyl_data_table,out_fn)
% Import Quantitative Structural Model (QSM) data.
% USE:
%     TREE = IMPORT_QSM_DATA(BRANCH_DATA_TABLE,CYL_DATA_TABLE,OUT_FN)
%
% INPUT:
%
% BRANCH_DATA_FN - branches table-data file name or data matrix
% CYL_DATA_FN - cylinder table-data file name or data matrix
% OUT_FN - (optional) output .mat file name (default is 'tree_out.mat')
% 
% The following format for the tables is dictate by the QSM algorithm
% generating the tables for purposes other than use in BayesForest Toolbox.
% This format, hence, can be changed in the future to accomodate the Toolbox
% needs specifically.
%
% BRANCH TABLE(num branches x at least 5):
% Br Order | Br Parent | Br Volume (not needed) | Br Length | Br Angle |
%   ...    |   ...     |         ...            |   ...     |    ...   |
%
% CYLINDER/SEGMENT TABLE(num cylinders x at least 11):
% Seg Radius | Seg Len | Seg Start (X,Y,Z) | Seg Axis (X,Y,Z) | Seg Parent | Seg Extension | Seg Branch |
%   ...      |  ...    |  ... | ... | ...  |  ... | ... | ... |  ...       |      ...      |     ...    |
%
% OUTPUT:
%
% TREE - the tree object tree, which can be, for example, drawn with
% tree.draw command.
% 
% Additionally, this function generates the output .mat file with the information
% read from the tables (OUT_FN, default 'qsm_out.mat').


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
BVol = A(:,3);%l
BLen = A(:,4);%cm or m
BAng = A(:,5);%deg
clear A;
if(isa(cyl_data_table,'char'))
    A = importdata(cyl_data_table);
else
    A = cyl_data_table;
end
Rad = A(:,1);%m
Len = A(:,2);%cm or m
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

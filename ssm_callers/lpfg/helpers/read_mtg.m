function [tr,TS,br] = read_mtg(mtgfn)
% Read the Multi-scale Tree Graph format file and convert it to the TREE
% object, TS and BR structs (see LIGMAT package).
% Cyl's numbering does not change by deleted cyl's that have empty values.
%
% USAGE:
%      [TR,TS,BR] = READ_MTG(MTGFN)
%
% MTGFN - MTG file name.
%
% TR is the output tree.
% TS is the output tree segment struct.
% BR is the output branching struct.

% Create the pipe
fid = fopen(mtgfn);
% Create the fileformat
%M = 15;% number of fields in the format
format_line = '%s %d %d %f %f %s %s %s %f %f %f %d %d %s %d';
%format_line = '%s %d %d %f %f %s %s %s %d %d %s %d';
% Read the input file name
T = textscan(fid,format_line,'delimiter',' ');
% Number of lines read
N = length(T{1});
% Create output struct's
tr = tree();
TS = struct('age',[],'Wf0',[],'Wf',[],'hwR',[]);
ncyl = 0;% max cyl index to cover all cyl's
disp('Start reading MTG... (this may take some time)');
for ii=1:N
    % Cyl number and topology
    tok = regexp(T{1}{ii},'(\/|<|(C\d+\+))C(?<cyl>\d+)','names');
    %tok = regexp(T{1}{ii},'(\/|<|C(\d+)\+)C(\d+)','tokens');
%     if(T{1}{ii}(1) == '/') %First cyl
%         disp(['first cyl is ' tok.cyl]);
%         last_cyl = int32(str2double(tok.cyl)) + 1;
%     elseif(T{1}{ii}(1) == '<') % extension
%         disp(['Cyl is ' tok.cyl ' parent is ' num2str(last_cyl)]);
%     elseif(T{1}{ii}(1) == 'C')% child
%         p = find(T{1}{ii}(1:length('C'):end) == '+');
%         disp(['Parent is ' T{1}{ii}(1+length('C'):p-1) ' child is ' tok.cyl]);
%         last_cyl = int32(str2double(tok.cyl)) + 1;
%     end
    cyl = str2double(tok.cyl) + 1;
    if(cyl > ncyl)
        ncyl = cyl;
    end
    % Age
    TS(cyl).age = T{2}(ii);
    % Order
    tr.is_trunk(cyl) = ~T{3}(ii);
    % Length
    tr.length(cyl) = T{4}(ii);
    % Radius
    tr.radius(cyl) = T{5}(ii);
    % Start
    tr.start_point(cyl,:) = get_vector(T{6}{ii},0);
    % End
    tr.end_point(cyl,:) = get_vector(T{7}{ii},0);
    % Axis
    tr.axis(cyl,:) = get_vector(T{8}{ii},0);
    % Initial foliage mass
    TS(cyl).Wf0 = T{9}(ii);
    % Current foliage mass
    TS(cyl).Wf = T{10}(ii);
    % Heartwood radius
    TS(cyl).hwR = T{11}(ii);
    % Parent
    tr.parent(cyl) = T{12}(ii) + 1;% convert from C-style
    % Extension
    tr.extension(cyl) = T{13}(ii) + 1;
    % Children
    tr.children{cyl} = (get_vector(T{14}{ii},1))' + 1;
    % Printing for debugging
%     fprintf('Cyl %d:\n',cyl);
%     fprintf('\tR = %G, ',tr.radius(cyl));
%     fprintf('L = %G, ',tr.length(cyl));
%     fprintf('Sta = (%G,%G,%G), ',tr.start_point(cyl,1),tr.start_point(cyl,2),tr.start_point(cyl,3));
%     fprintf('End = (%G,%G,%G), ',tr.end_point(cyl,1),tr.end_point(cyl,2),tr.end_point(cyl,3));
%     fprintf('Axi = (%G,%G,%G), ',tr.axis(cyl,1),tr.axis(cyl,2),tr.axis(cyl,3));
%     disp(['Par = ' num2str(tr.parent(cyl))]);
%     disp(['Ext = ' num2str(tr.extension(cyl))]);
%     fprintf('Chi = \n');
%     for jj=1:length(tr.children{cyl})
%         disp(num2str(tr.children{cyl}(jj)));
%     end
%     fprintf('\n');
end
tr.number_of_branches = ncyl;
%disp(['I''ve read ' num2str(N) ' lines and ' num2str(ncyl) ' cylinders.']);
% Get branch information
br = get_branches(tr);

disp('Done');
fclose(fid);
end

function vec = get_vector(str,chi)
% Reads the vector in the format (X,Y,Z)

str = str(2:end-1);% strip the parenthesis
if(isempty(str))
    vec = [];
    return;
end
%if(chi)% This is for children array
%    c = textscan(str,'%d','delimiter',',');
%else% this is for vectors
    c = textscan(str,'%f','delimiter',',');
%end
vec = (cell2mat(c))';
if(~chi && length(vec) == 3)% Rotation of the coordinates.
    vec = [vec(1) -vec(3) vec(2)];
end
end

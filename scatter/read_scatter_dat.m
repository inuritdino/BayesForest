function [s,type_names] = read_scatter_dat(fn,varargin)
% Read the structured scatter data from LPFG-LIGNUM.
% USAGE:
%       [S,TYPES] = READ_SCATTER_DAT(FN)
%
% S is the cell array of scatters/distribution functions.
%
% TYPES is the cell array of the names of the corresponding scatters.
%
% FN is the file name for the scatter data.
%
% SCATTER TYPES SUPPORTED:
% 'taper' - tapering function
% 'bra' - branching angles
% 'curv' - 3D curvature in space
% 'lchi_lapar' - length of children branches as a function of length along
%   the parent.
% 'lchi_bra_lapar' - the same as previous but augmented with the branching
%   angle the children emanate at from the parent.
% 'az' - azimuthal angles of the offsets around their parent segments.
% 'ltot_rini' - totla length of a branch vs. radius at its base.

MAXORDER = 10;
tf = strcmpi('MaxOrder',varargin);
if(find(tf))
    MAXORDER = varargin{find(tf)+1};
end
fid = fopen(fn);

order = 0;
s = cell(1,10);% up to 10 types
type_names = cell(1,10);
type_id = 0;
n_types = 0;
for ii=1:10
    s{ii} = cell(1,MAXORDER+1);% up to MAXORDER orders (+ 0 order)
end

while 1
    line = fgetl(fid);% get line
    if(~ischar(line)), break; end;% EOF break
    if(isempty(line)), continue; end;% empty line skip
    if(line(1) == '#')% Comment-instruction
        line = strrep(line,' ','');%remove blanks
        line = strrep(line,'#','');%remove #
        if(isempty(strfind(line,'order')))% if NOT order-instruction
            type = line;%get type
            type_id = find(strcmp(type,type_names),1);
            if(isempty(type_id))
                type_names{n_types+1} = type;
                n_types = n_types + 1;
                type_id = n_types;
            end
            
            %fprintf('Scatter type: %s: %d\n',type,type_id);
        else% order-instruction
            order = str2double(line((strfind(line,'order')+5):end));
            if(order > MAXORDER)%do not load order higher than MAXORDER
                fprintf('WARNING: max order reached (%d). Use ''MaxOrder'' option.\n',MAXORDER);
                break;
            end
            %fprintf('Order: %G\n',order);
        end
    else
        if(strcmp(type,'bra') || strcmp(type,'az'))
            C = textscan(line,'%f');
            s{type_id}{order+1} = cat(1,s{type_id}{order+1},C{1});
        elseif(strcmp(type,'taper') || strcmp(type,'lchi_lapar') || ...
                strcmp(type,'ltot_rini'))
            C = textscan(line,'%f %f');
            s{type_id}{order+1} = cat(1,s{type_id}{order+1},[C{1} C{2}]);
        elseif(strcmp(type,'curv') || strcmp(type,'lchi_bra_lapar'))
            C = textscan(line,'%f %f %f');
            s{type_id}{order+1} = cat(1,s{type_id}{order+1},[C{1} C{2} C{3}]);
        else
            fprintf('I do not know the type `%s''.\n',type);
        end
    end
    
end

fclose(fid);

end

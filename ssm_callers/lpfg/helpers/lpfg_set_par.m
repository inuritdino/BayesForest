function lpfg_set_par(fname,varargin)
% Set user defined parameters to a file (e.g. 'user.h') for an LPFG model
% USAGE:
%      LPFG_SET_PAR(FNAME,...)
%
% FNAME is the filename to write to.
% 
% ...,'PARNAME',VALUE,... pairs determine the "#define PARNAME VALUE"
% entries to the FNAME

% Open the file
fid = fopen(fname,'w');

% INPUT must go in pairs (...,'PARNAME',VALUE,...)
ii = 1;
while(ii <= length(varargin))
    if(isa(varargin{ii},'char'))
        if(all(isstrprop(varargin{ii},'lower')))%Convert lower case input
            fprintf(fid,'#define %s ',upper(varargin{ii}));
        else
            fprintf(fid,'#define %s ',varargin{ii});
        end
        if(isa(varargin{ii+1},'double'))
            fprintf(fid,'%G\n',varargin{ii+1});
        elseif(isa(varargin{ii+1},'char'))% SPECIAL CONSTRUCTS
            if(strncmp(varargin{ii+1},'gauss',5))% Gaussian variate
                par = get_par(varargin{ii+1}(6:end));
                fprintf(fid,'(ran_gauss(%s,%s))\n',num2str(par(1)),num2str(par(2)));
            elseif(strncmp(varargin{ii+1},'uniform',7))% Uniform variate
                par = get_par(varargin{ii+1}(8:end));
                fprintf(fid,'(ran_uniform(%s,%s))\n',num2str(par(1)),num2str(par(2)));
            else
                fprintf('ERROR: unrecognized option %s\n',varargin{ii+1});
            end
        end
    end
    ii = ii + 2;
end

fclose(fid);
end

function vec = get_par(str)
% Reads the vector of parameters in the format (X,Y,Z,...)

str = str(2:end-1);% strip the parenthesis
if(isempty(str))
    vec = [];
    return;
end
c = textscan(str,'%f','delimiter',',');
vec = (cell2mat(c))';

end